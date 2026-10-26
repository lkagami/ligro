#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
'''
O raio não é uma constante fixa. Ele vem de self.radius na 
SticksRepresentation (que carrega de sticks_radius no config), 
enviado ao shader pelo uniform vert_rad. O que é constante é 
só o multiplicador de afinamento das múltiplas (o 0.6). 

Então: para mudar a espessura geral, ajuste sticks_radius/self.radius; 
para mudar só as múltiplas em relação às simples, ajuste o 0.6 no shader.
[2026-10-14] Superado: a tripla agora e' eixo + dois cilindros lado a lado
(como em lines) e as duplas ficam no plano do anel/substituintes
(u_bond_side_tbo). O texto abaixo descreve a versao anterior.
Tripla triangular: aplicado. Em vez de três tubos numa linha, 
agora os três ficam nos vértices de um triângulo equilátero, 
no plano perpendicular à ligação. Para isso criei um segundo 
vetor perpendicular, g_up_dir = cross(eixo, g_side_dir), formando 
o plano da seção transversal, e cada passada ocupa um vértice a 120° 
(90° + 120°·u_pass). A dupla continua linear (dois tubos lado a lado),
 só a tripla mudou.
 
 
Um detalhe para calibrar: o raio do triângulo usa o mesmo u_separation
da dupla, o que pode espalhar os três tubos um pouco mais do que você
quer. Se a tripla ficar "aberta" demais, há duas saídas — reduzir
o 1.8 da separação (no representations.py), que afeta dupla e 
tripla juntas, ou, se quiser controlar a tripla independentemente,
multiplicar o u_separation só no ramo order == 3 (por exemplo * 0.7).
Me diga como ficou visualmente e eu ajusto o fator certo. Também 
dá para orientar o triângulo pela câmera (uma face sempre voltada 
para você) se preferir, mas o triângulo fixo no espaço do objeto 
costuma ler melhor para tripla.
'''



vertex_shader_ribbon_stick = """
#version 330

layout(std140) uniform CameraMatrices {
    mat4 view_mat;
    mat4 proj_mat;
};
precision highp float; 
precision highp int;
uniform mat4 model_mat;

in vec3 vert_coord;
in vec3 vert_color;
const float vert_rad = 0.4;

out vec3 geom_color;
out vec4 geom_coord;
out float geom_rad;

void main(){
    geom_color = vert_color;
    geom_rad = vert_rad;
    geom_coord = view_mat * model_mat * vec4(vert_coord, 1.0);
}
"""

vertex_shader_sticks = """
#version 330

layout(std140) uniform CameraMatrices {
    mat4 view_mat;
    mat4 proj_mat;
};
precision highp float; 
precision highp int;
uniform mat4 model_mat;

in vec3 vert_coord;
in vec3 vert_color;
//const float vert_rad = 0.10;
uniform float vert_rad; // float uniform
// [ATUALIZACAO] A ordem da ligacao NAO vem mais como atributo por-vertice
// (vert_bond_order foi removido). O problema do atributo por-vertice: ele e
// indexado por ATOMO (via GL_ELEMENT_ARRAY_BUFFER = index_bonds), entao um
// atomo com duas ligacoes de ordens diferentes so podia guardar UM valor
// (por isso a solucao antiga usava a MAIOR ordem entre as ligacoes do atomo
// -- efetivamente "valencia maxima", nao a ordem real de cada ligacao).
//
// Agora a ordem e lida diretamente no GEOMETRY SHADER via uma Buffer Texture
// (u_bond_order_tbo) + gl_PrimitiveIDIn, que da o indice da ligacao (0-based,
// mesma ordem de self.index_bonds) sendo desenhada em cada invocacao do GS.
// Isso pareia bond_order_list[k] com o k-esimo PAR de self.index_bonds de
// forma exata, sem nenhuma ambiguidade por atomo compartilhado. Ver
// geometry_shader_sticks abaixo e representations.py (SticksRepresentation).

out vec3 geom_color;
out vec4 geom_coord;
out float geom_rad;

void main(){
    geom_color = vert_color;
    geom_rad = vert_rad;
    geom_coord = view_mat * model_mat * vec4(vert_coord, 1.0);
}
"""

geometry_shader_sticks = """
#version 330

layout(std140) uniform CameraMatrices {
    mat4 view_mat;
    mat4 proj_mat;
};
precision highp float; 
precision highp int;
layout (lines) in;
// UM cilindro por invocacao do GS (40 vertices) -> sempre cabe no limite
// GL_MAX_GEOMETRY_TOTAL_OUTPUT_COMPONENTS (1024 e comum). Ligacoes multiplas
// sao desenhadas em PASSADAS separadas no lado Python (u_pass = 0,1,2), cada
// passada deslocando o cilindro lateralmente; passadas que nao se aplicam a
// ordem da ligacao nao emitem nada. Assim duplas/triplas funcionam sem
// estourar o limite de saida do geometry shader.
layout (triangle_strip, max_vertices = 40) out;

// Passada atual (0 = central/primeiro, 1 e 2 = laterais) e separacao base
// entre cilindros de uma mesma ligacao. Setados por draw_representation.
uniform int u_pass;
uniform float u_separation;
// [EN] 2026-10-15: radius of each cylinder and spacing of double and triple
// bonds, in Angstrom (Preferences, shared by sticks and dynamic bonds; the
// single-bond radius is vert_rad). Spacing: double = distance between the two
// cylinder axes, triple = distance from the bond axis to each cylinder (radius
// of the triangle). <= 0: fall back to the old rules (0.4 x vert_rad,
// u_separation).
uniform float u_rad_double;
uniform float u_rad_triple;
uniform float u_sep_double;
uniform float u_sep_triple;

// [ATUALIZACAO] Ordem da ligacao lida por PRIMITIVA (nao mais por vertice).
// u_bond_order_tbo e uma Buffer Texture (GL_TEXTURE_BUFFER, formato R32UI)
// com exatamente uma entrada por ligacao, na MESMA ordem em que os pares
// aparecem em self.index_bonds -- ou seja, texel k == bond_order_list[k].
// gl_PrimitiveIDIn e o indice (0-based) da primitiva GL_LINES sendo
// processada nesta invocacao do geometry shader dentro do glDrawElements
// atual; como cada par de index_bonds vira exatamente uma primitiva
// GL_LINES, gl_PrimitiveIDIn == k. Isso pareia a ordem 1:1 com cada par de
// index_bonds, sem qualquer ambiguidade por atomo compartilhado (o problema
// do esquema antigo, que usava um atributo por-vertice/por-atomo).
// [EN] 2026-09-28 (macOS report, "unsigned-integer texture warning during
// picking"): the picking program never binds this texture, and in a core
// profile an unbound sampler does NOT read as 0 -- it samples whatever
// texture sits on unit 0 (type mismatch -> driver warning/undefined
// result). The picking variant (sel_geometry_shader_sticks, below) is
// compiled with PICKING_PASS defined, which removes the sampler entirely.
#ifndef PICKING_PASS
uniform usamplerBuffer u_bond_order_tbo;
// [EN] 2026-10-14: plane of the double bonds (same style as "lines"). One
// RGBA32F texel per bond: xyz = unit vector perpendicular to the bond in the
// ring / substituent plane (object space), w = 0 ring double (axis cylinder +
// shorter inner cylinder), 1 double in the substituent plane, 2 offset on the
// screen (triple bonds, no substituent). u_has_side = 0: no texture bound,
// every bond uses the screen.
uniform samplerBuffer u_bond_side_tbo;
uniform int u_has_side;
uniform mat4 model_mat;
#endif

in vec3 geom_color[];
in vec4 geom_coord[];
in float geom_rad[];

out vec3 frag_coord;
out vec3 frag_color;
out vec3 frag_norm;

// This data is used for the cylinder vertices, i.e. the points that form the
// circle in a horizontal cut of the cylinder. The first half are the points
// at the begining of the cylinder, and the last half are the points at the end.
// The quantity of points should be changed to get a smoother cylinder, but that
// will result in more resources used.
vec3 bs_0 = vec3( 1.000000000000, 0.000000000000, 0.000000000000); // base
vec3 bs_1 = vec3( 0.766044443119, 0.000000000000, 0.642787609687); // base
vec3 bs_2 = vec3( 0.173648177667, 0.000000000000, 0.984807753012); // base
vec3 bs_3 = vec3(-0.500000000000, 0.000000000000, 0.866025403784); // base
vec3 bs_4 = vec3(-0.939692620786, 0.000000000000, 0.342020143326); // base
vec3 bs_5 = vec3(-0.939692620786, 0.000000000000,-0.342020143326); // base
vec3 bs_6 = vec3(-0.500000000000, 0.000000000000,-0.866025403784); // base
vec3 bs_7 = vec3( 0.173648177667, 0.000000000000,-0.984807753012); // base
vec3 bs_8 = vec3( 0.766044443119, 0.000000000000,-0.642787609687); // base
vec3 up_0 = vec3( 0.939692620786, 0.000000000000, 0.342020143326); // up
vec3 up_1 = vec3( 0.500000000000, 0.000000000000, 0.866025403784); // up
vec3 up_2 = vec3(-0.173648177667, 0.000000000000, 0.984807753012); // up
vec3 up_3 = vec3(-0.766044443119, 0.000000000000, 0.642787609687); // up
vec3 up_4 = vec3(-1.000000000000, 0.000000000000, 0.000000000000); // up
vec3 up_5 = vec3(-0.766044443119, 0.000000000000,-0.642787609687); // up
vec3 up_6 = vec3(-0.173648177667, 0.000000000000,-0.984807753012); // up
vec3 up_7 = vec3( 0.500000000000, 0.000000000000,-0.866025403784); // up
vec3 up_8 = vec3( 0.939692620786, 0.000000000000,-0.342020143326); // up

// The rotation matrix used for translating the cylinder points to their correct
// places is rot_mat. This matrix is created using the my_glRotatef function,
// see the function documentation to get more information.
mat3 rot_mat;

// mid_point is the middle point in the line.
vec3 mid_point;

// These points are the vertices calculated for the cylinder.
vec3 p_00, p_01, p_02, p_03, p_04, p_05, p_06, p_07, p_08;
vec3 p_09, p_10, p_11, p_12, p_13, p_14, p_15, p_16, p_17;
vec3 p_18, p_19, p_20, p_21, p_22, p_23, p_24, p_25, p_26;

// Centros efetivos do cilindro atual. Para ligacao simples sao iguais a
// geom_coord[0]/[1] e ao ponto medio. Para dupla/tripla, cada cilindro e
// deslocado lateralmente por 'bond_offset' (calculado no main a partir de
// um vetor perpendicular a ligacao, orientado pela camera). calculate_points
// e a emissao usam center0/center1/mid_c em vez de geom_coord/mid_point.
vec3 center0;
vec3 center1;
vec3 mid_c;
// Raio efetivo do cilindro atual. Igual a geom_rad[0] para ligacoes simples;
// reduzido para ligacoes multiplas (definido no main antes de emit_cylinder).
float eff_rad;

float get_angle(vec3 vec_A, vec3 vec_B){
    // Returns the angle in radians formed between the vectors A and B. The
    // vectors are initially normalized to avoid errors.
    // The initial result is clamped to the [-1,+1] range to avoid errors in
    // the arc cosine function.
    vec3 vecA_u = normalize(vec_A);
    vec3 vecB_u = normalize(vec_B);
    return acos(clamp(dot(vecA_u, vecB_u), -1.0, 1.0));
}

mat3 get_rot_mat(float my_angle, vec3 dir_vec){
    // The get_rot_mat creates a rotation matrix using an angle and a direction
    // vector. This matrix is used for move points acording to a defined angle
    // in a defined position. We'll use it to obtain the cylinder vertices at
    // the correct orientations.
    mat3 my_mat = mat3(0.0);
    // ndv stands for normalized direction vector
    vec3 ndv = normalize(dir_vec);
    float cosa = cos(my_angle);
    float sina = sin(my_angle);
    my_mat[0][0] = ndv.x*ndv.x*(1-cosa)+cosa;
    my_mat[1][0] = ndv.x*ndv.y*(1-cosa)+ndv.z*sina;
    my_mat[2][0] = ndv.x*ndv.z*(1-cosa)-ndv.y*sina;
    my_mat[0][1] = ndv.x*ndv.y*(1-cosa)-ndv.z*sina;
    my_mat[1][1] = ndv.y*ndv.y*(1-cosa)+cosa;
    my_mat[2][1] = ndv.y*ndv.z*(1-cosa)+ndv.x*sina;
    my_mat[0][2] = ndv.x*ndv.z*(1-cosa)+ndv.y*sina;
    my_mat[1][2] = ndv.y*ndv.z*(1-cosa)-ndv.x*sina;
    my_mat[2][2] = ndv.z*ndv.z*(1-cosa)+cosa;
    return my_mat;
}

void calculate_points(){
    // This void function fills the vertices of the sticks with data.
    p_00.x = (bs_0.x*rot_mat[0][0] + bs_0.y*rot_mat[0][1] + bs_0.z*rot_mat[0][2])*eff_rad + center0.x;
    p_00.y = (bs_0.x*rot_mat[1][0] + bs_0.y*rot_mat[1][1] + bs_0.z*rot_mat[1][2])*eff_rad + center0.y;
    p_00.z = (bs_0.x*rot_mat[2][0] + bs_0.y*rot_mat[2][1] + bs_0.z*rot_mat[2][2])*eff_rad + center0.z;
    p_01.x = (bs_1.x*rot_mat[0][0] + bs_1.y*rot_mat[0][1] + bs_1.z*rot_mat[0][2])*eff_rad + center0.x;
    p_01.y = (bs_1.x*rot_mat[1][0] + bs_1.y*rot_mat[1][1] + bs_1.z*rot_mat[1][2])*eff_rad + center0.y;
    p_01.z = (bs_1.x*rot_mat[2][0] + bs_1.y*rot_mat[2][1] + bs_1.z*rot_mat[2][2])*eff_rad + center0.z;
    p_02.x = (bs_2.x*rot_mat[0][0] + bs_2.y*rot_mat[0][1] + bs_2.z*rot_mat[0][2])*eff_rad + center0.x;
    p_02.y = (bs_2.x*rot_mat[1][0] + bs_2.y*rot_mat[1][1] + bs_2.z*rot_mat[1][2])*eff_rad + center0.y;
    p_02.z = (bs_2.x*rot_mat[2][0] + bs_2.y*rot_mat[2][1] + bs_2.z*rot_mat[2][2])*eff_rad + center0.z;
    p_03.x = (bs_3.x*rot_mat[0][0] + bs_3.y*rot_mat[0][1] + bs_3.z*rot_mat[0][2])*eff_rad + center0.x;
    p_03.y = (bs_3.x*rot_mat[1][0] + bs_3.y*rot_mat[1][1] + bs_3.z*rot_mat[1][2])*eff_rad + center0.y;
    p_03.z = (bs_3.x*rot_mat[2][0] + bs_3.y*rot_mat[2][1] + bs_3.z*rot_mat[2][2])*eff_rad + center0.z;
    p_04.x = (bs_4.x*rot_mat[0][0] + bs_4.y*rot_mat[0][1] + bs_4.z*rot_mat[0][2])*eff_rad + center0.x;
    p_04.y = (bs_4.x*rot_mat[1][0] + bs_4.y*rot_mat[1][1] + bs_4.z*rot_mat[1][2])*eff_rad + center0.y;
    p_04.z = (bs_4.x*rot_mat[2][0] + bs_4.y*rot_mat[2][1] + bs_4.z*rot_mat[2][2])*eff_rad + center0.z;
    p_05.x = (bs_5.x*rot_mat[0][0] + bs_5.y*rot_mat[0][1] + bs_5.z*rot_mat[0][2])*eff_rad + center0.x;
    p_05.y = (bs_5.x*rot_mat[1][0] + bs_5.y*rot_mat[1][1] + bs_5.z*rot_mat[1][2])*eff_rad + center0.y;
    p_05.z = (bs_5.x*rot_mat[2][0] + bs_5.y*rot_mat[2][1] + bs_5.z*rot_mat[2][2])*eff_rad + center0.z;
    p_06.x = (bs_6.x*rot_mat[0][0] + bs_6.y*rot_mat[0][1] + bs_6.z*rot_mat[0][2])*eff_rad + center0.x;
    p_06.y = (bs_6.x*rot_mat[1][0] + bs_6.y*rot_mat[1][1] + bs_6.z*rot_mat[1][2])*eff_rad + center0.y;
    p_06.z = (bs_6.x*rot_mat[2][0] + bs_6.y*rot_mat[2][1] + bs_6.z*rot_mat[2][2])*eff_rad + center0.z;
    p_07.x = (bs_7.x*rot_mat[0][0] + bs_7.y*rot_mat[0][1] + bs_7.z*rot_mat[0][2])*eff_rad + center0.x;
    p_07.y = (bs_7.x*rot_mat[1][0] + bs_7.y*rot_mat[1][1] + bs_7.z*rot_mat[1][2])*eff_rad + center0.y;
    p_07.z = (bs_7.x*rot_mat[2][0] + bs_7.y*rot_mat[2][1] + bs_7.z*rot_mat[2][2])*eff_rad + center0.z;
    p_08.x = (bs_8.x*rot_mat[0][0] + bs_8.y*rot_mat[0][1] + bs_8.z*rot_mat[0][2])*eff_rad + center0.x;
    p_08.y = (bs_8.x*rot_mat[1][0] + bs_8.y*rot_mat[1][1] + bs_8.z*rot_mat[1][2])*eff_rad + center0.y;
    p_08.z = (bs_8.x*rot_mat[2][0] + bs_8.y*rot_mat[2][1] + bs_8.z*rot_mat[2][2])*eff_rad + center0.z;
    
    p_09.x = (up_0.x*rot_mat[0][0] + up_0.y*rot_mat[0][1] + up_0.z*rot_mat[0][2])*eff_rad + mid_c.x;
    p_09.y = (up_0.x*rot_mat[1][0] + up_0.y*rot_mat[1][1] + up_0.z*rot_mat[1][2])*eff_rad + mid_c.y;
    p_09.z = (up_0.x*rot_mat[2][0] + up_0.y*rot_mat[2][1] + up_0.z*rot_mat[2][2])*eff_rad + mid_c.z;
    p_10.x = (up_1.x*rot_mat[0][0] + up_1.y*rot_mat[0][1] + up_1.z*rot_mat[0][2])*eff_rad + mid_c.x;
    p_10.y = (up_1.x*rot_mat[1][0] + up_1.y*rot_mat[1][1] + up_1.z*rot_mat[1][2])*eff_rad + mid_c.y;
    p_10.z = (up_1.x*rot_mat[2][0] + up_1.y*rot_mat[2][1] + up_1.z*rot_mat[2][2])*eff_rad + mid_c.z;
    p_11.x = (up_2.x*rot_mat[0][0] + up_2.y*rot_mat[0][1] + up_2.z*rot_mat[0][2])*eff_rad + mid_c.x;
    p_11.y = (up_2.x*rot_mat[1][0] + up_2.y*rot_mat[1][1] + up_2.z*rot_mat[1][2])*eff_rad + mid_c.y;
    p_11.z = (up_2.x*rot_mat[2][0] + up_2.y*rot_mat[2][1] + up_2.z*rot_mat[2][2])*eff_rad + mid_c.z;
    p_12.x = (up_3.x*rot_mat[0][0] + up_3.y*rot_mat[0][1] + up_3.z*rot_mat[0][2])*eff_rad + mid_c.x;
    p_12.y = (up_3.x*rot_mat[1][0] + up_3.y*rot_mat[1][1] + up_3.z*rot_mat[1][2])*eff_rad + mid_c.y;
    p_12.z = (up_3.x*rot_mat[2][0] + up_3.y*rot_mat[2][1] + up_3.z*rot_mat[2][2])*eff_rad + mid_c.z;
    p_13.x = (up_4.x*rot_mat[0][0] + up_4.y*rot_mat[0][1] + up_4.z*rot_mat[0][2])*eff_rad + mid_c.x;
    p_13.y = (up_4.x*rot_mat[1][0] + up_4.y*rot_mat[1][1] + up_4.z*rot_mat[1][2])*eff_rad + mid_c.y;
    p_13.z = (up_4.x*rot_mat[2][0] + up_4.y*rot_mat[2][1] + up_4.z*rot_mat[2][2])*eff_rad + mid_c.z;
    p_14.x = (up_5.x*rot_mat[0][0] + up_5.y*rot_mat[0][1] + up_5.z*rot_mat[0][2])*eff_rad + mid_c.x;
    p_14.y = (up_5.x*rot_mat[1][0] + up_5.y*rot_mat[1][1] + up_5.z*rot_mat[1][2])*eff_rad + mid_c.y;
    p_14.z = (up_5.x*rot_mat[2][0] + up_5.y*rot_mat[2][1] + up_5.z*rot_mat[2][2])*eff_rad + mid_c.z;
    p_15.x = (up_6.x*rot_mat[0][0] + up_6.y*rot_mat[0][1] + up_6.z*rot_mat[0][2])*eff_rad + mid_c.x;
    p_15.y = (up_6.x*rot_mat[1][0] + up_6.y*rot_mat[1][1] + up_6.z*rot_mat[1][2])*eff_rad + mid_c.y;
    p_15.z = (up_6.x*rot_mat[2][0] + up_6.y*rot_mat[2][1] + up_6.z*rot_mat[2][2])*eff_rad + mid_c.z;
    p_16.x = (up_7.x*rot_mat[0][0] + up_7.y*rot_mat[0][1] + up_7.z*rot_mat[0][2])*eff_rad + mid_c.x;
    p_16.y = (up_7.x*rot_mat[1][0] + up_7.y*rot_mat[1][1] + up_7.z*rot_mat[1][2])*eff_rad + mid_c.y;
    p_16.z = (up_7.x*rot_mat[2][0] + up_7.y*rot_mat[2][1] + up_7.z*rot_mat[2][2])*eff_rad + mid_c.z;
    p_17.x = (up_8.x*rot_mat[0][0] + up_8.y*rot_mat[0][1] + up_8.z*rot_mat[0][2])*eff_rad + mid_c.x;
    p_17.y = (up_8.x*rot_mat[1][0] + up_8.y*rot_mat[1][1] + up_8.z*rot_mat[1][2])*eff_rad + mid_c.y;
    p_17.z = (up_8.x*rot_mat[2][0] + up_8.y*rot_mat[2][1] + up_8.z*rot_mat[2][2])*eff_rad + mid_c.z;
    
    p_18.x = (bs_0.x*rot_mat[0][0] + bs_0.y*rot_mat[0][1] + bs_0.z*rot_mat[0][2])*eff_rad + center1.x;
    p_18.y = (bs_0.x*rot_mat[1][0] + bs_0.y*rot_mat[1][1] + bs_0.z*rot_mat[1][2])*eff_rad + center1.y;
    p_18.z = (bs_0.x*rot_mat[2][0] + bs_0.y*rot_mat[2][1] + bs_0.z*rot_mat[2][2])*eff_rad + center1.z;
    p_19.x = (bs_1.x*rot_mat[0][0] + bs_1.y*rot_mat[0][1] + bs_1.z*rot_mat[0][2])*eff_rad + center1.x;
    p_19.y = (bs_1.x*rot_mat[1][0] + bs_1.y*rot_mat[1][1] + bs_1.z*rot_mat[1][2])*eff_rad + center1.y;
    p_19.z = (bs_1.x*rot_mat[2][0] + bs_1.y*rot_mat[2][1] + bs_1.z*rot_mat[2][2])*eff_rad + center1.z;
    p_20.x = (bs_2.x*rot_mat[0][0] + bs_2.y*rot_mat[0][1] + bs_2.z*rot_mat[0][2])*eff_rad + center1.x;
    p_20.y = (bs_2.x*rot_mat[1][0] + bs_2.y*rot_mat[1][1] + bs_2.z*rot_mat[1][2])*eff_rad + center1.y;
    p_20.z = (bs_2.x*rot_mat[2][0] + bs_2.y*rot_mat[2][1] + bs_2.z*rot_mat[2][2])*eff_rad + center1.z;
    p_21.x = (bs_3.x*rot_mat[0][0] + bs_3.y*rot_mat[0][1] + bs_3.z*rot_mat[0][2])*eff_rad + center1.x;
    p_21.y = (bs_3.x*rot_mat[1][0] + bs_3.y*rot_mat[1][1] + bs_3.z*rot_mat[1][2])*eff_rad + center1.y;
    p_21.z = (bs_3.x*rot_mat[2][0] + bs_3.y*rot_mat[2][1] + bs_3.z*rot_mat[2][2])*eff_rad + center1.z;
    p_22.x = (bs_4.x*rot_mat[0][0] + bs_4.y*rot_mat[0][1] + bs_4.z*rot_mat[0][2])*eff_rad + center1.x;
    p_22.y = (bs_4.x*rot_mat[1][0] + bs_4.y*rot_mat[1][1] + bs_4.z*rot_mat[1][2])*eff_rad + center1.y;
    p_22.z = (bs_4.x*rot_mat[2][0] + bs_4.y*rot_mat[2][1] + bs_4.z*rot_mat[2][2])*eff_rad + center1.z;
    p_23.x = (bs_5.x*rot_mat[0][0] + bs_5.y*rot_mat[0][1] + bs_5.z*rot_mat[0][2])*eff_rad + center1.x;
    p_23.y = (bs_5.x*rot_mat[1][0] + bs_5.y*rot_mat[1][1] + bs_5.z*rot_mat[1][2])*eff_rad + center1.y;
    p_23.z = (bs_5.x*rot_mat[2][0] + bs_5.y*rot_mat[2][1] + bs_5.z*rot_mat[2][2])*eff_rad + center1.z;
    p_24.x = (bs_6.x*rot_mat[0][0] + bs_6.y*rot_mat[0][1] + bs_6.z*rot_mat[0][2])*eff_rad + center1.x;
    p_24.y = (bs_6.x*rot_mat[1][0] + bs_6.y*rot_mat[1][1] + bs_6.z*rot_mat[1][2])*eff_rad + center1.y;
    p_24.z = (bs_6.x*rot_mat[2][0] + bs_6.y*rot_mat[2][1] + bs_6.z*rot_mat[2][2])*eff_rad + center1.z;
    p_25.x = (bs_7.x*rot_mat[0][0] + bs_7.y*rot_mat[0][1] + bs_7.z*rot_mat[0][2])*eff_rad + center1.x;
    p_25.y = (bs_7.x*rot_mat[1][0] + bs_7.y*rot_mat[1][1] + bs_7.z*rot_mat[1][2])*eff_rad + center1.y;
    p_25.z = (bs_7.x*rot_mat[2][0] + bs_7.y*rot_mat[2][1] + bs_7.z*rot_mat[2][2])*eff_rad + center1.z;
    p_26.x = (bs_8.x*rot_mat[0][0] + bs_8.y*rot_mat[0][1] + bs_8.z*rot_mat[0][2])*eff_rad + center1.x;
    p_26.y = (bs_8.x*rot_mat[1][0] + bs_8.y*rot_mat[1][1] + bs_8.z*rot_mat[1][2])*eff_rad + center1.y;
    p_26.z = (bs_8.x*rot_mat[2][0] + bs_8.y*rot_mat[2][1] + bs_8.z*rot_mat[2][2])*eff_rad + center1.z;
}

// Vetor perpendicular a ligacao, orientado pela camera, usado para deslocar
// os cilindros de ligacoes multiplas para os lados (do ponto de vista da
// camera). Mesma ideia do shader de lines: cross(eixo_da_ligacao, dir_camera).
// Calculado uma vez no main e reutilizado por cilindro.
vec3 g_side_dir;

// Eixo/rotacao do cilindro: identicos para todos os cilindros de uma mesma
// ligacao, calculados uma vez no main.
// (rot_mat ja e global; aqui so documentamos.)

void emit_cylinder(){
    calculate_points();
    // Now we send the vertices to the fragment shader in a defined order
    // base-> 0, 9, 1, 10, 2, 11, 3, 12, 4, 13, 5, 14, 6, 15, 7, 16, 8, 17, 0, 9
    gl_Position = proj_mat * vec4(p_00, 1);
    frag_coord = p_00;
    frag_color = geom_color[0];
    frag_norm = p_00 - center0.xyz;
    EmitVertex();
    gl_Position = proj_mat * vec4(p_09, 1);
    frag_coord = p_09;
    frag_color = geom_color[0];
    frag_norm = p_09 - mid_c;
    EmitVertex();
    gl_Position = proj_mat * vec4(p_01, 1);
    frag_coord = p_01;
    frag_color = geom_color[0];
    frag_norm = p_01 - center0.xyz;
    EmitVertex();
    gl_Position = proj_mat * vec4(p_10, 1);
    frag_coord = p_10;
    frag_color = geom_color[0];
    frag_norm = p_10 - mid_c;
    EmitVertex();
    gl_Position = proj_mat * vec4(p_02, 1);
    frag_coord = p_02;
    frag_color = geom_color[0];
    frag_norm = p_02 - center0.xyz;
    EmitVertex();
    gl_Position = proj_mat * vec4(p_11, 1);
    frag_coord = p_11;
    frag_color = geom_color[0];
    frag_norm = p_11 - mid_c;
    EmitVertex();
    gl_Position = proj_mat * vec4(p_03, 1);
    frag_coord = p_03;
    frag_color = geom_color[0];
    frag_norm = p_03 - center0.xyz;
    EmitVertex();
    gl_Position = proj_mat * vec4(p_12, 1);
    frag_coord = p_12;
    frag_color = geom_color[0];
    frag_norm = p_12 - mid_c;
    EmitVertex();
    gl_Position = proj_mat * vec4(p_04, 1);
    frag_coord = p_04;
    frag_color = geom_color[0];
    frag_norm = p_04 - center0.xyz;
    EmitVertex();
    gl_Position = proj_mat * vec4(p_13, 1);
    frag_coord = p_13;
    frag_color = geom_color[0];
    frag_norm = p_13 - mid_c;
    EmitVertex();
    gl_Position = proj_mat * vec4(p_05, 1);
    frag_coord = p_05;
    frag_color = geom_color[0];
    frag_norm = p_05 - center0.xyz;
    EmitVertex();
    gl_Position = proj_mat * vec4(p_14, 1);
    frag_coord = p_14;
    frag_color = geom_color[0];
    frag_norm = p_14 - mid_c;
    EmitVertex();
    gl_Position = proj_mat * vec4(p_06, 1);
    frag_coord = p_06;
    frag_color = geom_color[0];
    frag_norm = p_06 - center0.xyz;
    EmitVertex();
    gl_Position = proj_mat * vec4(p_15, 1);
    frag_coord = p_15;
    frag_color = geom_color[0];
    frag_norm = p_15 - mid_c;
    EmitVertex();
    gl_Position = proj_mat * vec4(p_07, 1);
    frag_coord = p_07;
    frag_color = geom_color[0];
    frag_norm = p_07 - center0.xyz;
    EmitVertex();
    gl_Position = proj_mat * vec4(p_16, 1);
    frag_coord = p_16;
    frag_color = geom_color[0];
    frag_norm = p_16 - mid_c;
    EmitVertex();
    gl_Position = proj_mat * vec4(p_08, 1);
    frag_coord = p_08;
    frag_color = geom_color[0];
    frag_norm = p_08 - center0.xyz;
    EmitVertex();
    gl_Position = proj_mat * vec4(p_17, 1);
    frag_coord = p_17;
    frag_color = geom_color[0];
    frag_norm = p_17 - mid_c;
    EmitVertex();
    gl_Position = proj_mat * vec4(p_00, 1);
    frag_coord = p_00;
    frag_color = geom_color[0];
    frag_norm = p_00 - center0.xyz;
    EmitVertex();
    gl_Position = proj_mat * vec4(p_09, 1);
    frag_coord = p_09;
    frag_color = geom_color[0];
    frag_norm = p_09 - mid_c;
    EmitVertex();
    EndPrimitive();
    
    // up-> 9, 18, 10, 19, 11, 20, 12, 21, 13, 22, 14, 23, 15, 24, 16, 25, 17, 26, 9, 18
    gl_Position = proj_mat * vec4(p_09, 1);
    frag_coord = p_09;
    frag_color = geom_color[1];
    frag_norm = p_09 - mid_c;
    EmitVertex();
    gl_Position = proj_mat * vec4(p_18, 1);
    frag_coord = p_18;
    frag_color = geom_color[1];
    frag_norm = p_18 - center1.xyz;
    EmitVertex();
    gl_Position = proj_mat * vec4(p_10, 1);
    frag_coord = p_10;
    frag_color = geom_color[1];
    frag_norm = p_10 - mid_c;
    EmitVertex();
    gl_Position = proj_mat * vec4(p_19, 1);
    frag_coord = p_19;
    frag_color = geom_color[1];
    frag_norm = p_19 - center1.xyz;
    EmitVertex();
    gl_Position = proj_mat * vec4(p_11, 1);
    frag_coord = p_11;
    frag_color = geom_color[1];
    frag_norm = p_11 - mid_c;
    EmitVertex();
    gl_Position = proj_mat * vec4(p_20, 1);
    frag_coord = p_20;
    frag_color = geom_color[1];
    frag_norm = p_20 - center1.xyz;
    EmitVertex();
    gl_Position = proj_mat * vec4(p_12, 1);
    frag_coord = p_12;
    frag_color = geom_color[1];
    frag_norm = p_12 - mid_c;
    EmitVertex();
    gl_Position = proj_mat * vec4(p_21, 1);
    frag_coord = p_21;
    frag_color = geom_color[1];
    frag_norm = p_21 - center1.xyz;
    EmitVertex();
    gl_Position = proj_mat * vec4(p_13, 1);
    frag_coord = p_13;
    frag_color = geom_color[1];
    frag_norm = p_13 - mid_c;
    EmitVertex();
    gl_Position = proj_mat * vec4(p_22, 1);
    frag_coord = p_22;
    frag_color = geom_color[1];
    frag_norm = p_22 - center1.xyz;
    EmitVertex();
    gl_Position = proj_mat * vec4(p_14, 1);
    frag_coord = p_14;
    frag_color = geom_color[1];
    frag_norm = p_14 - mid_c;
    EmitVertex();
    gl_Position = proj_mat * vec4(p_23, 1);
    frag_coord = p_23;
    frag_color = geom_color[1];
    frag_norm = p_23 - center1.xyz;
    EmitVertex();
    gl_Position = proj_mat * vec4(p_15, 1);
    frag_coord = p_15;
    frag_color = geom_color[1];
    frag_norm = p_15 - mid_c;
    EmitVertex();
    gl_Position = proj_mat * vec4(p_24, 1);
    frag_coord = p_24;
    frag_color = geom_color[1];
    frag_norm = p_24 - center1.xyz;
    EmitVertex();
    gl_Position = proj_mat * vec4(p_16, 1);
    frag_coord = p_16;
    frag_color = geom_color[1];
    frag_norm = p_16 - mid_c;
    EmitVertex();
    gl_Position = proj_mat * vec4(p_25, 1);
    frag_coord = p_25;
    frag_color = geom_color[1];
    frag_norm = p_25 - center1.xyz;
    EmitVertex();
    gl_Position = proj_mat * vec4(p_17, 1);
    frag_coord = p_17;
    frag_color = geom_color[1];
    frag_norm = p_17 - mid_c;
    EmitVertex();
    gl_Position = proj_mat * vec4(p_26, 1);
    frag_coord = p_26;
    frag_color = geom_color[1];
    frag_norm = p_26 - center1.xyz;
    EmitVertex();
    gl_Position = proj_mat * vec4(p_09, 1);
    frag_coord = p_09;
    frag_color = geom_color[1];
    frag_norm = p_09 - mid_c;
    EmitVertex();
    gl_Position = proj_mat * vec4(p_18, 1);
    frag_coord = p_18;
    frag_color = geom_color[1];
    frag_norm = p_18 - center1.xyz;
    EmitVertex();
    EndPrimitive();
}

void main(){
    // Ordem da ligacao: um texel por ligacao, indexado pela primitiva atual
    // (gl_PrimitiveIDIn == k-esimo par de self.index_bonds == bond_order_list[k]).
    // No picking (PICKING_PASS) a ordem nao importa -- desenha so a passada 0
    // com u_separation=0 -- entao e fixada em 1 (ligacao simples), sem ler
    // nenhuma textura (uma textura NAO vinculada nao retorna 0 no core
    // profile; le o que estiver na unidade 0).
#ifdef PICKING_PASS
    int order = 1;
#else
    int order = int(texelFetch(u_bond_order_tbo, gl_PrimitiveIDIn).r);
#endif
    if (order < 1) order = 1;
    if (order > 3) order = 3;

    // Esta passada produz geometria para esta ligacao?
    //   u_pass 0: sempre (cilindro central/primeiro)
    //   u_pass 1: ordem >= 2
    //   u_pass 2: ordem >= 3
    if (u_pass > 0 && order < u_pass + 1) {
        return; // nada a emitir nesta passada
    }

    // Eixo da ligacao e matriz de rotacao do cilindro.
    vec3 vec_p0_p1 = geom_coord[1].xyz - geom_coord[0].xyz;
    vec3 ort_vec = normalize(cross(vec3(0,1,0), vec_p0_p1));
    float g_angle = get_angle(vec3(0,1,0), vec_p0_p1);
    rot_mat = get_rot_mat(g_angle, ort_vec);

    // Direcao lateral orientada pela camera (separa os cilindros lado a lado
    // na tela em qualquer orientacao). Fallback para ort_vec se degenerar.
    vec3 mid_view = (geom_coord[0].xyz + geom_coord[1].xyz) * 0.5;
    vec3 to_cam = normalize(-mid_view);
    g_side_dir = cross(normalize(vec_p0_p1), to_cam);
    if (length(g_side_dir) < 1e-4) {
        g_side_dir = ort_vec;
    }
    g_side_dir = normalize(g_side_dir);
    // [EN] 2026-10-14: double bonds in the plane of the ring / substituents
    // (u_bond_side_tbo, see above); otherwise g_side_dir (screen).
    float side_mode = 2.0;
#ifndef PICKING_PASS
    if (order > 1 && u_has_side == 1) {
        vec4 s = texelFetch(u_bond_side_tbo, gl_PrimitiveIDIn);
        if (s.w < 1.5 || s.w > 2.5) {
            vec3 axis_v = normalize(vec_p0_p1);
            vec3 plane = mat3(view_mat * model_mat) * s.xyz;
            plane = plane - dot(plane, axis_v) * axis_v;
            if (length(plane) > 1e-4) {
                g_side_dir = normalize(plane);
                side_mode = s.w;
            }
        }
    }
#endif

    // Offset no plano perpendicular, conforme a ordem e a passada.
    //   double (ring or not): two thin cylinders of full length at +-sep/2, in
    //     the ring / substituent plane (2026-10-15, user: the axis cylinder +
    //     shorter inner cylinder used for ring doubles looked worse in sticks
    //     than the old symmetric pair -- "lines" keeps the inner line);
    //   triple: three thin cylinders on a triangle fixed to the molecule
    //     (u_bond_side_tbo w = 3); without that data, axis + two at +-sep.
    vec3 off = vec3(0.0);
    float shrink = 0.0;
    eff_rad = geom_rad[0];
    if (order == 2) {
        float sep2 = (u_sep_double > 0.0) ? u_sep_double : u_separation;
        float d = (u_pass == 0) ? (-sep2 * 0.5) : (sep2 * 0.5);
        off = g_side_dir * d;
        eff_rad = (u_rad_double > 0.0) ? u_rad_double : geom_rad[0] * 0.4;
    } else if (order == 3) {
        float sep3 = (u_sep_triple > 0.0) ? u_sep_triple : u_separation;
        if (side_mode > 2.5) {
            // [EN] 2026-10-15: triangle fixed to the molecule (w = 3): the three
            // cylinders at 90/210/330 degrees in the plane perpendicular to the
            // bond, starting from the object-space vector of u_bond_side_tbo,
            // so the triangle turns with the molecule like the double bonds.
            vec3 up_dir = normalize(cross(normalize(vec_p0_p1), g_side_dir));
            float ang = radians(90.0 + 120.0 * float(u_pass));
            off = (g_side_dir * cos(ang) + up_dir * sin(ang)) * sep3;
        } else {
            if (u_pass == 1) off = g_side_dir * sep3;
            if (u_pass == 2) off = -g_side_dir * sep3;
        }
        eff_rad = (u_rad_triple > 0.0) ? u_rad_triple : geom_rad[0] * 0.4;
    }

    center0 = geom_coord[0].xyz + off + vec_p0_p1 * shrink;
    center1 = geom_coord[1].xyz + off - vec_p0_p1 * shrink;
    mid_c = (center0 + center1) * 0.5;
    emit_cylinder();
}
"""

fragment_shader_sticks = """
#version 330
precision highp float; 
precision highp int;

#include <vm_lighting>

uniform vec4 fog_color;
uniform float fog_start;
uniform float fog_end;
uniform float fog_intensity;

in vec3 frag_coord;
in vec3 frag_color;
in vec3 frag_norm;

out vec4 final_color;

vec4 calculate_color(vec3 fnrm, vec3 fcrd, vec3 fcol){
    return vec4(vm_lighting(normalize(fnrm), normalize(fcrd), fcol, 0.0), 1.0);
}

void main(){
    vec4 out_color = calculate_color(frag_norm, frag_coord, frag_color);
    
    float dist = abs(frag_coord.z);
    if(dist>=fog_start){
        float fog_factor = 1.0 - fog_intensity * (1.0 - clamp((fog_end-dist)/(fog_end-fog_start), 0.0, 1.0));
        final_color = mix(fog_color, out_color, fog_factor);
    }
    else{
       final_color = out_color;
    }
}
"""

sel_fragment_shader_sticks = """
#version 330
precision highp float; 
precision highp int;

in vec3 frag_coord;
in vec3 frag_color;
in vec3 frag_norm;

out vec4 final_color;

void main(){
    float alpha = length(normalize(frag_coord + frag_norm));
    final_color = vec4(frag_color, alpha);
    //final_color =  frag_color
}
"""


vertex_shader_impostor_sticks = """
#version 330

precision highp float;

layout(std140) uniform CameraMatrices {
    mat4 view_mat;
    mat4 proj_mat;
};

uniform mat4 model_mat;

in vec3 vert_coord;
in vec3 vert_color;
//in float vert_rad;
const float vert_rad = 0.1;

uniform vec3 u_campos;

out vec3 geom_color;
out vec3 geom_coord;
out vec3 geom_center;
out vec3 geom_cam;
out float geom_radius;

void main() {
    geom_color = vert_color;
    geom_coord = vert_coord;
    geom_center = vert_coord;
    geom_cam = u_campos;
    geom_radius = vert_rad;
}
"""

geometry_shader_impostor_sticks = """
#version 330

layout(std140) uniform CameraMatrices {
    mat4 view_mat;
    mat4 proj_mat;
};

layout (lines) in;
layout (triangle_strip, max_vertices = 18) out;

uniform mat4 model_mat;

in vec3 geom_color[];
in vec3 geom_coord[];
in vec3 geom_center[];
in vec3 geom_cam[];
in float geom_radius[];

out vec3 frag_color;
out vec3 frag_coord;
out vec3 frag_pA;
out vec3 frag_pB;
out vec3 frag_cam;
out float frag_radius;

void main(){
    vec3 AB_dir = normalize(geom_coord[1] - geom_coord[0]);
    vec3 side_dir = normalize(geom_cam[0] - geom_coord[0]);
    vec3 base_dir = normalize(cross(AB_dir, side_dir));
    side_dir = normalize(cross(AB_dir, base_dir));
    
    float rad_fac = geom_radius[0] * 1.415;
    vec3 base_1 = geom_coord[0] + base_dir * rad_fac - AB_dir * rad_fac;
    vec3 base_2 = geom_coord[0] + side_dir * rad_fac - AB_dir * rad_fac;
    vec3 base_3 = geom_coord[0] - base_dir * rad_fac - AB_dir * rad_fac;
    vec3 base_4 = geom_coord[0] - side_dir * rad_fac - AB_dir * rad_fac;
    vec3 top_1 = geom_coord[1] + base_dir * rad_fac + AB_dir * rad_fac;
    vec3 top_2 = geom_coord[1] + side_dir * rad_fac + AB_dir * rad_fac;
    vec3 top_3 = geom_coord[1] - base_dir * rad_fac + AB_dir * rad_fac;
    vec3 top_4 = geom_coord[1] - side_dir * rad_fac + AB_dir * rad_fac;
    
    mat4 pvm_mat = proj_mat * view_mat * model_mat;
    mat4 vm_mat = model_mat;
    
    gl_Position = pvm_mat * vec4(base_1, 1.0);
    frag_color = geom_color[0];
    frag_coord = vec3(vm_mat * vec4(base_1, 1.0));
    frag_pA = vec3(vm_mat * vec4(geom_center[0], 1.0));
    frag_pB = vec3(vm_mat * vec4(geom_center[1], 1.0));
    frag_cam = vec3(vm_mat * vec4(geom_cam[0], 1.0));
    frag_radius = geom_radius[0];
    EmitVertex();
    gl_Position = pvm_mat * vec4(top_1, 1.0);
    frag_color = geom_color[1];
    frag_coord = vec3(vm_mat * vec4(top_1, 1.0));
    frag_pA = vec3(vm_mat * vec4(geom_center[0], 1.0));
    frag_pB = vec3(vm_mat * vec4(geom_center[1], 1.0));
    frag_cam = vec3(vm_mat * vec4(geom_cam[0], 1.0));
    frag_radius = geom_radius[0];
    EmitVertex();
    gl_Position = pvm_mat * vec4(base_2, 1.0);
    frag_color = geom_color[0];
    frag_coord = vec3(vm_mat * vec4(base_2, 1.0));
    frag_pA = vec3(vm_mat * vec4(geom_center[0], 1.0));
    frag_pB = vec3(vm_mat * vec4(geom_center[1], 1.0));
    frag_cam = vec3(vm_mat * vec4(geom_cam[0], 1.0));
    frag_radius = geom_radius[0];
    EmitVertex();
    gl_Position = pvm_mat * vec4(top_2, 1.0);
    frag_color = geom_color[1];
    frag_coord = vec3(vm_mat * vec4(top_2, 1.0));
    frag_pA = vec3(vm_mat * vec4(geom_center[0], 1.0));
    frag_pB = vec3(vm_mat * vec4(geom_center[1], 1.0));
    frag_cam = vec3(vm_mat * vec4(geom_cam[0], 1.0));
    frag_radius = geom_radius[0];
    EmitVertex();
    gl_Position = pvm_mat * vec4(base_3, 1.0);
    frag_color = geom_color[0];
    frag_coord = vec3(vm_mat * vec4(base_3, 1.0));
    frag_pA = vec3(vm_mat * vec4(geom_center[0], 1.0));
    frag_pB = vec3(vm_mat * vec4(geom_center[1], 1.0));
    frag_cam = vec3(vm_mat * vec4(geom_cam[0], 1.0));
    frag_radius = geom_radius[0];
    EmitVertex();
    gl_Position = pvm_mat * vec4(top_3, 1.0);
    frag_color = geom_color[1];
    frag_coord = vec3(vm_mat * vec4(top_3, 1.0));
    frag_pA = vec3(vm_mat * vec4(geom_center[0], 1.0));
    frag_pB = vec3(vm_mat * vec4(geom_center[1], 1.0));
    frag_cam = vec3(vm_mat * vec4(geom_cam[0], 1.0));
    frag_radius = geom_radius[0];
    EmitVertex();
    gl_Position = pvm_mat * vec4(base_4, 1.0);
    frag_color = geom_color[0];
    frag_coord = vec3(vm_mat * vec4(base_4, 1.0));
    frag_pA = vec3(vm_mat * vec4(geom_center[0], 1.0));
    frag_pB = vec3(vm_mat * vec4(geom_center[1], 1.0));
    frag_cam = vec3(vm_mat * vec4(geom_cam[0], 1.0));
    frag_radius = geom_radius[0];
    EmitVertex();
    gl_Position = pvm_mat * vec4(top_4, 1.0);
    frag_color = geom_color[1];
    frag_coord = vec3(vm_mat * vec4(top_4, 1.0));
    frag_pA = vec3(vm_mat * vec4(geom_center[0], 1.0));
    frag_pB = vec3(vm_mat * vec4(geom_center[1], 1.0));
    frag_cam = vec3(vm_mat * vec4(geom_cam[0], 1.0));
    frag_radius = geom_radius[0];
    EmitVertex();
    gl_Position = pvm_mat * vec4(base_1, 1.0);
    frag_color = geom_color[0];
    frag_coord = vec3(vm_mat * vec4(base_1, 1.0));
    frag_pA = vec3(vm_mat * vec4(geom_center[0], 1.0));
    frag_pB = vec3(vm_mat * vec4(geom_center[1], 1.0));
    frag_cam = vec3(vm_mat * vec4(geom_cam[0], 1.0));
    frag_radius = geom_radius[0];
    EmitVertex();
    gl_Position = pvm_mat * vec4(top_1, 1.0);
    frag_color = geom_color[1];
    frag_coord = vec3(vm_mat * vec4(top_1, 1.0));
    frag_pA = vec3(vm_mat * vec4(geom_center[0], 1.0));
    frag_pB = vec3(vm_mat * vec4(geom_center[1], 1.0));
    frag_cam = vec3(vm_mat * vec4(geom_cam[0], 1.0));
    frag_radius = geom_radius[0];
    EmitVertex();
    
    EndPrimitive();
    gl_Position = pvm_mat * vec4(base_1, 1.0);
    frag_color = geom_color[0];
    frag_coord = vec3(vm_mat * vec4(base_1, 1.0));
    frag_pA = vec3(vm_mat * vec4(geom_center[0], 1.0));
    frag_pB = vec3(vm_mat * vec4(geom_center[1], 1.0));
    frag_cam = vec3(vm_mat * vec4(geom_cam[0], 1.0));
    frag_radius = geom_radius[0];
    EmitVertex();
    gl_Position = pvm_mat * vec4(base_2, 1.0);
    frag_color = geom_color[0];
    frag_coord = vec3(vm_mat * vec4(base_2, 1.0));
    frag_pA = vec3(vm_mat * vec4(geom_center[0], 1.0));
    frag_pB = vec3(vm_mat * vec4(geom_center[1], 1.0));
    frag_cam = vec3(vm_mat * vec4(geom_cam[0], 1.0));
    frag_radius = geom_radius[0];
    EmitVertex();
    gl_Position = pvm_mat * vec4(base_4, 1.0);
    frag_color = geom_color[0];
    frag_coord = vec3(vm_mat * vec4(base_4, 1.0));
    frag_pA = vec3(vm_mat * vec4(geom_center[0], 1.0));
    frag_pB = vec3(vm_mat * vec4(geom_center[1], 1.0));
    frag_cam = vec3(vm_mat * vec4(geom_cam[0], 1.0));
    frag_radius = geom_radius[0];
    EmitVertex();
    gl_Position = pvm_mat * vec4(base_3, 1.0);
    frag_color = geom_color[0];
    frag_coord = vec3(vm_mat * vec4(base_3, 1.0));
    frag_pA = vec3(vm_mat * vec4(geom_center[0], 1.0));
    frag_pB = vec3(vm_mat * vec4(geom_center[1], 1.0));
    frag_cam = vec3(vm_mat * vec4(geom_cam[0], 1.0));
    frag_radius = geom_radius[0];
    EmitVertex();
    EndPrimitive();

    gl_Position = pvm_mat * vec4(top_1, 1.0);
    frag_color = geom_color[1];
    frag_coord = vec3(vm_mat * vec4(top_1, 1.0));
    frag_pA = vec3(vm_mat * vec4(geom_center[1], 1.0));
    frag_pB = vec3(vm_mat * vec4(geom_center[1], 1.0));
    frag_cam = vec3(vm_mat * vec4(geom_cam[0], 1.0));
    frag_radius = geom_radius[0];
    EmitVertex();
    gl_Position = pvm_mat * vec4(top_2, 1.0);
    frag_color = geom_color[1];
    frag_coord = vec3(vm_mat * vec4(top_2, 1.0));
    frag_pA = vec3(vm_mat * vec4(geom_center[1], 1.0));
    frag_pB = vec3(vm_mat * vec4(geom_center[1], 1.0));
    frag_cam = vec3(vm_mat * vec4(geom_cam[0], 1.0));
    frag_radius = geom_radius[0];
    EmitVertex();
    gl_Position = pvm_mat * vec4(top_4, 1.0);
    frag_color = geom_color[1];
    frag_coord = vec3(vm_mat * vec4(top_4, 1.0));
    frag_pA = vec3(vm_mat * vec4(geom_center[1], 1.0));
    frag_pB = vec3(vm_mat * vec4(geom_center[1], 1.0));
    frag_cam = vec3(vm_mat * vec4(geom_cam[0], 1.0));
    frag_radius = geom_radius[0];
    EmitVertex();
    gl_Position = pvm_mat * vec4(top_3, 1.0);
    frag_color = geom_color[1];
    frag_coord = vec3(vm_mat * vec4(top_3, 1.0));
    frag_pA = vec3(vm_mat * vec4(geom_center[1], 1.0));
    frag_pB = vec3(vm_mat * vec4(geom_center[1], 1.0));
    frag_cam = vec3(vm_mat * vec4(geom_cam[0], 1.0));
    frag_radius = geom_radius[0];
    EmitVertex();
    EndPrimitive();
}
"""

fragment_shader_impostor_sticks = """
#version 330

layout(std140) uniform CameraMatrices {
    mat4 view_mat;
    mat4 proj_mat;
};
#extension GL_EXT_frag_depth: enable
precision highp float;

#include <vm_lighting>

//uniform float u_depth;

uniform vec4 fog_color;
uniform float fog_start;
uniform float fog_end;
uniform float fog_intensity;

in vec3 frag_color;
in vec3 frag_coord;
in vec3 frag_pA;
in vec3 frag_pB;
in vec3 frag_cam;
in float frag_radius;

out vec4 final_color;

vec3 caps_normal( vec3 p, vec3 a, vec3 b, float r ){
  vec3 pa = p - a;
  vec3 ba = b - a;
  float h = clamp( dot(pa,ba)/dot(ba,ba), 0.0, 1.0 );
  return ( pa - ba*h ) / r;
}

float caps_intersect( in vec3 ro, in vec3 rd, in vec3 pa, in vec3 pb, in float ra ){
    vec3  ba = pb - pa;
    vec3  oa = ro - pa;
    float baba = dot(ba,ba);
    float bard = dot(ba,rd);
    float baoa = dot(ba,oa);
    float rdoa = dot(rd,oa);
    float oaoa = dot(oa,oa);
    float a = baba      - bard*bard;
    float b = baba*rdoa - baoa*bard;
    float c = baba*oaoa - baoa*baoa - ra*ra*baba;
    float h = b*b - a*c;
    if( h >= 0.0 )
    {
        float t = (-b-sqrt(h))/a;
        float y = baoa + t*bard;
        // body
        if( y>0.0 && y<baba ) return t;
        // caps
        vec3 oc = (y <= 0.0) ? oa : ro - pb;
        b = dot(rd,oc);
        c = dot(oc,oc) - ra*ra;
        h = b*b - c;
        if( h>0.0 ) return -b - sqrt(h);
    }
    return -1.0;
}

vec4 calculate_color(vec3 fnrm, vec3 fcrd, vec3 fcol){
    return vec4(vm_lighting(normalize(fnrm), normalize(fcrd), fcol, 0.0), 1.0);
}

void main() {
    vec3 ray_orig = frag_cam;
    vec3 ray_dir = normalize(frag_coord - frag_cam);
    float dist_to_caps = caps_intersect(ray_orig, ray_dir, frag_pA, frag_pB, frag_radius);
    if (dist_to_caps < 0.0) discard;
    
    vec3 coord_on_caps = frag_cam + ray_dir * dist_to_caps;
    vec3 normal = normalize(caps_normal(frag_coord, frag_pA, frag_pB, frag_radius));
    
    final_color = calculate_color(normal, coord_on_caps, frag_color);
    //gl_FragDepth = depth.z/depth.w;
}
"""

sel_fragment_shader_impostor_sticks = """
#version 330

layout(std140) uniform CameraMatrices {
    mat4 view_mat;
    mat4 proj_mat;
};
#extension GL_EXT_frag_depth: enable
precision highp float;


in vec3 frag_color;
in vec3 frag_coord;
in vec3 frag_pA;
in vec3 frag_pB;
in vec3 frag_cam;
in float frag_radius;

out vec4 final_color;

float caps_intersect( in vec3 ro, in vec3 rd, in vec3 pa, in vec3 pb, in float ra ){
    vec3  ba = pb - pa;
    vec3  oa = ro - pa;
    float baba = dot(ba,ba);
    float bard = dot(ba,rd);
    float baoa = dot(ba,oa);
    float rdoa = dot(rd,oa);
    float oaoa = dot(oa,oa);
    float a = baba      - bard*bard;
    float b = baba*rdoa - baoa*bard;
    float c = baba*oaoa - baoa*baoa - ra*ra*baba;
    float h = b*b - a*c;
    if( h >= 0.0 )
    {
        float t = (-b-sqrt(h))/a;
        float y = baoa + t*bard;
        // body
        if( y>0.0 && y<baba ) return t;
        // caps
        vec3 oc = (y <= 0.0) ? oa : ro - pb;
        b = dot(rd,oc);
        c = dot(oc,oc) - ra*ra;
        h = b*b - c;
        if( h>0.0 ) return -b - sqrt(h);
    }
    return -1.0;
}

void main() {
    vec3 ray_orig = frag_cam;
    vec3 ray_dir = normalize(frag_coord - frag_cam);
    float dist_to_caps = caps_intersect(ray_orig, ray_dir, frag_pA, frag_pB, frag_radius);
    if (dist_to_caps < 0.0) discard;
    final_color = vec4(frag_color, 1.0);
    
    //vec4 depth = proj_mat * vec4(coord_on_caps, 1.0);
    //gl_FragDepth = depth.z/depth.w;
}
"""

# Picking variant of geometry_shader_sticks: identical except PICKING_PASS is
# defined right after the #version line (GLSL requires #version first), which
# drops the bond-order sampler -- see the comment at u_bond_order_tbo above.
sel_geometry_shader_sticks = geometry_shader_sticks.replace(
    "#version 330", "#version 330\n#define PICKING_PASS", 1)

shader_type ={0: {"vertex_shader"      : vertex_shader_sticks,
                  "geometry_shader"    : geometry_shader_sticks,
                  "fragment_shader"    : fragment_shader_sticks,
                  "sel_vertex_shader"  : vertex_shader_sticks,
                  "sel_geometry_shader": sel_geometry_shader_sticks,
                  "sel_fragment_shader": sel_fragment_shader_sticks
                 },
              1: {"vertex_shader"      : vertex_shader_impostor_sticks,
                  "geometry_shader"    : geometry_shader_impostor_sticks,
                  "fragment_shader"    : fragment_shader_impostor_sticks,
                  "sel_vertex_shader"  : vertex_shader_impostor_sticks  ,
                  "sel_geometry_shader": geometry_shader_impostor_sticks,
                  "sel_fragment_shader": sel_fragment_shader_impostor_sticks
                  },
              2: {"vertex_shader"      : vertex_shader_ribbon_stick, 
                  "geometry_shader"    : geometry_shader_sticks,     
                  "fragment_shader"    : fragment_shader_sticks,     
                  "sel_vertex_shader"  : vertex_shader_ribbon_stick, 
                  "sel_geometry_shader": sel_geometry_shader_sticks,     
                  "sel_fragment_shader": sel_fragment_shader_sticks  
                    },
               3: {"vertex_shader"     : vertex_shader_ribbon_stick,
                  "geometry_shader"    : geometry_shader_sticks,
                  "fragment_shader"    : fragment_shader_sticks,
                  "sel_vertex_shader"  : vertex_shader_sticks,
                  "sel_geometry_shader": sel_geometry_shader_sticks,
                  "sel_fragment_shader": sel_fragment_shader_sticks
                 },  
 
             
  
        
             }
