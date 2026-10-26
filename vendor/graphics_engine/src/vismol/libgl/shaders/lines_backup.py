#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#

geometry_shader_lines_bad = """
#version 330 core

layout (lines) in;
layout (triangle_strip, max_vertices = 4) out;

uniform mat4 proj_mat;
uniform mat4 view_mat;

in vec4 geom_coord[];
in vec3 geom_color[];
in float geom_radius[];

out vec2 frag_uv;
out vec3 frag_color;
out vec3 frag_viewdir;

void main() {
    vec3 p0 = geom_coord[0].xyz;
    vec3 p1 = geom_coord[1].xyz;

    // Direção da ligação
    vec3 bond_dir = normalize(p1 - p0);
    // Vetor ortogonal para expandir o quad (tela)
    vec3 view_dir = normalize(vec3(inverse(view_mat)[3]) - p0);
    vec3 ortho = normalize(cross(bond_dir, view_dir));
    
    float radius = 0.5 * (geom_radius[0] + geom_radius[1]) * 0.5;

    vec3 offset = ortho * radius;

    vec4 pos0 = vec4(p0 + offset, 1.0);
    vec4 pos1 = vec4(p0 - offset, 1.0);
    vec4 pos2 = vec4(p1 + offset, 1.0);
    vec4 pos3 = vec4(p1 - offset, 1.0);

    frag_color = geom_color[0];
    frag_uv = vec2(0.0, 0.0); frag_viewdir = view_dir;
    gl_Position = proj_mat * view_mat * pos0; EmitVertex();

    frag_uv = vec2(0.0, 1.0);
    gl_Position = proj_mat * view_mat * pos1; EmitVertex();

    frag_color = geom_color[1];
    frag_uv = vec2(1.0, 0.0);
    gl_Position = proj_mat * view_mat * pos2; EmitVertex();

    frag_uv = vec2(1.0, 1.0);
    gl_Position = proj_mat * view_mat * pos3; EmitVertex();

    EndPrimitive();
}

"""
fragment_shader_lines_bad = """
#version 330 core

in vec2 frag_uv;
in vec3 frag_color;
in vec3 frag_viewdir;

out vec4 out_color;

void main() {
    float x = frag_uv.x * 2.0 - 1.0;  // -1 to 1
    float y = frag_uv.y * 2.0 - 1.0;

    float dist_sq = x*x;
    if (dist_sq > 1.0) discard;  // fora do cilindro

    // Simula iluminação difusa em torno do cilindro
    vec3 normal = vec3(x, 0.0, sqrt(1.0 - dist_sq));
    vec3 light_dir = normalize(vec3(0.3, 0.5, 1.0));  // direcional
    float diffuse = max(dot(normal, light_dir), 0.0);

    vec3 color = frag_color * (0.3 + 0.7 * diffuse);
    out_color = vec4(color, 1.0);
}

"""




vertex_shader_lines   = """
#version 330
precision highp float;
precision highp int;

uniform mat4 model_mat;
uniform mat4 view_mat;

in vec3 vert_coord;
in vec3 vert_color;
in float vert_dot_size;  // new:   van der Waals radii

out vec3 geom_color;
out vec4 geom_coord;
out float geom_radius; // send to geometry shader

void main() {
    // Transformação para coordenadas da câmera (view space)
    vec4 world_pos = model_mat * vec4(vert_coord, 1.0);
    geom_coord = view_mat * world_pos;
    geom_color = vert_color;
    geom_radius = vert_dot_size;
}
"""

geometry_shader_lines = """
#version 330
precision highp float;
precision highp int;

layout (lines) in;
layout (line_strip, max_vertices = 8) out;  // Corrigido: máximo de 8 vértices emitidos

uniform mat4 proj_mat;
uniform mat4 view_mat;     // Adicionado
const float separation = 0.06;

in vec3 geom_color[];
in vec4 geom_coord[];
in float geom_radius[];

out vec3 frag_color;
out vec4 frag_coord;

vec4 get_vec(vec3 point_A, vec3 point_B, vec3 campos){
    vec3 tmp = (point_B - point_A);
    tmp = cross(campos, tmp);
    tmp = normalize(tmp);
    return vec4(tmp, 0.0);
}

void main() {
    float r0 = geom_radius[0];
    float r1 = geom_radius[1];
    float sum_radius = (r0 + r1)*0.87;
    
    //sum_radius = sum_radius;
    
    
    vec4 p0 = geom_coord[0];
    vec4 p1 = geom_coord[1];
    //float dist = length(p1.xyz - p0.xyz);
    float dist = sum_radius+1;
    
    
    if (dist < sum_radius) {
        vec4 mid_coord = vec4((p0.xyz + p1.xyz) * 0.5, 1.0);
        vec3 campos = vec3(inverse(view_mat)[3]);  // Corrigido: pegando posição da câmera
        vec3 vort = get_vec(p0.xyz, p1.xyz, campos).xyz;  // Corrigido: extrair vec3 do vec4

        // Segmento inferior
        gl_Position = proj_mat * (p0 - vec4(vort * separation, 0.0));
        frag_color = geom_color[0];
        frag_coord = p0 - vec4(vort * separation, 0.0);
        EmitVertex();

        gl_Position = proj_mat * (mid_coord - vec4(vort * separation, 0.0));
        frag_color = geom_color[0];
        frag_coord = mid_coord - vec4(vort * separation, 0.0);
        EmitVertex();
        EndPrimitive();

        gl_Position = proj_mat * (mid_coord - vec4(vort * separation, 0.0));
        frag_color = geom_color[1];
        frag_coord = mid_coord - vec4(vort * separation, 0.0);
        EmitVertex();

        gl_Position = proj_mat * (p1 - vec4(vort * separation, 0.0));
        frag_color = geom_color[1];
        frag_coord = p1 - vec4(vort * separation, 0.0);
        EmitVertex();
        EndPrimitive();

        // Segmento superior
        gl_Position = proj_mat * (p0 + vec4(vort * separation, 0.0));
        frag_color = geom_color[0];
        frag_coord = p0 + vec4(vort * separation, 0.0);
        EmitVertex();

        gl_Position = proj_mat * (mid_coord + vec4(vort * separation, 0.0));
        frag_color = geom_color[0];
        frag_coord = mid_coord + vec4(vort * separation, 0.0);
        EmitVertex();
        EndPrimitive();

        gl_Position = proj_mat * (mid_coord + vec4(vort * separation, 0.0));
        frag_color = geom_color[1];
        frag_coord = mid_coord + vec4(vort * separation, 0.0);
        EmitVertex();

        gl_Position = proj_mat * (p1 + vec4(vort * separation, 0.0));
        frag_color = geom_color[1];
        frag_coord = p1 + vec4(vort * separation, 0.0);
        
        EmitVertex();
        EndPrimitive();
    
    } else {
        vec4 mid = vec4((p0.xyz + p1.xyz) * 0.5, 1.0);
        vec3 c0 = geom_color[0];
        vec3 c1 = geom_color[1];

        frag_color = c0;
        frag_coord = p0;
        gl_Position = proj_mat * p0;
        EmitVertex();

        frag_coord = mid;
        gl_Position = proj_mat * mid;
        EmitVertex();
        EndPrimitive();

        frag_color = c1;
        frag_coord = mid;
        gl_Position = proj_mat * mid;
        EmitVertex();

        frag_coord = p1;
        gl_Position = proj_mat * p1;
        //frag_color = vec3(geom_radius[0], 0.0, 0.0); // apagar depois
        EmitVertex();
        EndPrimitive();
    }
}
"""
fragment_shader_lines = """
#version 330 core
precision highp float;
precision highp int;

uniform vec4 fog_color;
uniform float fog_start;
uniform float fog_end;

in vec3 frag_color;
in vec4 frag_coord;

out vec4 final_color;

void main() {
    float dist = abs(frag_coord.z); // ou length(frag_coord.xyz), se preferir

    // Cálculo do fator de névoa linear
    float fog_factor = 1.0;
    if (dist >= fog_start) {
        fog_factor = clamp((fog_end - dist) / (fog_end - fog_start), 0.0, 1.0);
    }

    // Interpolação linear entre a cor da névoa e a cor do fragmento
    final_color = mix(fog_color, vec4(frag_color, 1.0), fog_factor);
}
"""





geometry_shader_lines_new  = """
#version 330

layout (lines) in;
layout (triangle_strip, max_vertices = 4) out;

uniform mat4 proj_mat;
uniform mat4 view_mat;
uniform vec3 camera_pos;

in vec3 geom_color[];
in vec4 geom_coord[];
in float geom_radius[];

out vec3 frag_color;
out vec3 frag_pos;
out vec3 frag_axis;
out float frag_radius;

void main() {
    vec3 p0 = geom_coord[0].xyz;
    vec3 p1 = geom_coord[1].xyz;

    float radius = (geom_radius[0] + geom_radius[1]) * 0.5;

    vec3 axis = normalize(p1 - p0);
    vec3 view_dir = normalize(camera_pos - p0);

    // vetor perpendicular (base do billboard)
    vec3 right = normalize(cross(axis, view_dir)) * radius;

    vec3 v0 = p0 - right;
    vec3 v1 = p0 + right;
    vec3 v2 = p1 - right;
    vec3 v3 = p1 + right;

    frag_axis = axis;
    frag_radius = radius;

    // tri strip
    frag_color = geom_color[0];
    frag_pos = v0;
    gl_Position = proj_mat * view_mat * vec4(v0, 1.0);
    EmitVertex();

    frag_pos = v1;
    gl_Position = proj_mat * view_mat * vec4(v1, 1.0);
    EmitVertex();

    frag_color = geom_color[1];
    frag_pos = v2;
    gl_Position = proj_mat * view_mat * vec4(v2, 1.0);
    EmitVertex();

    frag_pos = v3;
    gl_Position = proj_mat * view_mat * vec4(v3, 1.0);
    EmitVertex();

    EndPrimitive();
}"""

fragment_shader_lines_new  = """
#version 330

in vec3 frag_color;
in vec3 frag_pos;
in vec3 frag_axis;
in float frag_radius;

uniform vec3 camera_pos;

out vec4 out_color;

void main() {
    // vetor da câmera até o fragmento
    vec3 view_dir = normalize(frag_pos - camera_pos);

    // projetar ponto no eixo do cilindro
    vec3 axis = normalize(frag_axis);

    vec3 to_frag = frag_pos;

    float proj = dot(to_frag, axis);
    vec3 closest = proj * axis;

    vec3 radial = to_frag - closest;
    float dist = length(radial);

    if (dist > frag_radius) discard;

    // normal do cilindro
    vec3 normal = normalize(radial);

    // iluminação simples (Phong)
    vec3 light_dir = normalize(vec3(1.0, 1.0, 1.0));

    float diff = max(dot(normal, light_dir), 0.0);
    float spec = pow(max(dot(reflect(-light_dir, normal), -view_dir), 0.0), 32.0);

    vec3 color = frag_color * (0.2 + 0.8 * diff) + vec3(1.0) * spec * 0.5;

    out_color = vec4(color, 1.0);
}"""













################################## SELECTION ###################################

sel_vertex_shader_lines = """
#version 330
precision highp float; 
precision highp int;

uniform mat4 model_mat;
uniform mat4 view_mat;

in vec3 vert_coord;
in vec3 vert_color;

out vec3 geom_color;
out vec4 geom_coord;

void main(){
    geom_color = vert_color;
    geom_coord = view_mat * model_mat * vec4(vert_coord, 1.0);
}
"""
sel_geometry_shader_lines = """
#version 330
precision highp float; 
precision highp int;

layout (lines) in;
layout (line_strip, max_vertices = 4) out;

uniform mat4 proj_mat;

in vec3 geom_color[];
in vec4 geom_coord[];

out vec3 frag_color;

void main(){
    vec4 mid_coord = vec4((geom_coord[0].xyz + geom_coord[1].xyz)/2, 1.0);
    gl_Position = proj_mat * geom_coord[0];
    frag_color = geom_color[0];
    EmitVertex();
    gl_Position = proj_mat * mid_coord;
    frag_color = geom_color[0];
    EmitVertex();
    EndPrimitive();
    gl_Position = proj_mat * mid_coord;
    frag_color = geom_color[1];
    EmitVertex();
    gl_Position = proj_mat * geom_coord[1];
    frag_color = geom_color[1];
    EmitVertex();
    EndPrimitive();
}
"""
sel_fragment_shader_lines = """
#version 330
precision highp float; 
precision highp int;

in vec3 frag_color;

out vec4 final_color;

void main(){
    final_color = vec4(frag_color, 1.0);
}
"""




vertex_shader_nlines = """
#version 330
precision highp float; 
precision highp int;

uniform mat4 model_mat;
uniform mat4 view_mat;

in vec3 vert_coord;
in vec3 vert_color;
const float vert_width = 0.01810;
const float antialias_length = 0.040;

out vec3 geom_color;
out vec4 geom_coord;
out float geom_width;
out float geom_anti_l;

void main(){
    geom_color = vert_color;
    geom_width = vert_width;
    geom_anti_l = antialias_length;
    geom_coord = view_mat * model_mat * vec4(vert_coord, 1.0);
}
"""
geometry_shader_nlines = """
#version 330
precision highp float; 
precision highp int;

layout (lines) in;
layout (triangle_strip, max_vertices = 20) out;

uniform mat4 proj_mat;

in vec3 geom_color[];
in vec4 geom_coord[];
in float geom_width[];
in float geom_anti_l[];

varying vec4 point_a[4];
varying vec4 point_b[4];
varying vec4 point_c[4];

out vec3 frag_coord;
out vec3 frag_color;

void main(){
    vec4 mid_point = (geom_coord[0] + geom_coord[1])/2;
    vec2 dir_vec = normalize((geom_coord[1] - geom_coord[0]).xy);
    vec3 orto_vec = normalize(cross(vec3(dir_vec, 0), vec3(0, 0, 1)));
    vec3 alias_color = vec3(0, 0, 0);
    point_a[0] = vec4(geom_coord[0].xyz + orto_vec*geom_width[0] + orto_vec*geom_anti_l[0], 1.0);
    point_a[1] = vec4(geom_coord[0].xyz + orto_vec*geom_width[0] + orto_vec*geom_anti_l[0]/2, 1.0);
    point_a[2] = vec4(geom_coord[0].xyz - orto_vec*geom_width[0] - orto_vec*geom_anti_l[0]/2, 1.0);
    point_a[3] = vec4(geom_coord[0].xyz - orto_vec*geom_width[0] - orto_vec*geom_anti_l[0], 1.0);
    point_b[0] = vec4(mid_point.xyz + orto_vec*geom_width[0] + orto_vec*geom_anti_l[0], 1.0);
    point_b[1] = vec4(mid_point.xyz + orto_vec*geom_width[0] + orto_vec*geom_anti_l[0]/2, 1.0);
    point_b[2] = vec4(mid_point.xyz - orto_vec*geom_width[0] - orto_vec*geom_anti_l[0]/2, 1.0);
    point_b[3] = vec4(mid_point.xyz - orto_vec*geom_width[0] - orto_vec*geom_anti_l[0], 1.0);
    point_c[0] = vec4(geom_coord[1].xyz + orto_vec*geom_width[0] + orto_vec*geom_anti_l[0], 1.0);
    point_c[1] = vec4(geom_coord[1].xyz + orto_vec*geom_width[0] + orto_vec*geom_anti_l[0]/2, 1.0);
    point_c[2] = vec4(geom_coord[1].xyz - orto_vec*geom_width[0] - orto_vec*geom_anti_l[0]/2, 1.0);
    point_c[3] = vec4(geom_coord[1].xyz - orto_vec*geom_width[0] - orto_vec*geom_anti_l[0], 1.0);
    
    gl_Position = proj_mat * point_a[0];
    frag_coord = point_a[0].xyz;
    frag_color = alias_color;
    EmitVertex();
    gl_Position = proj_mat * point_b[0];
    frag_coord = point_b[0].xyz;
    frag_color = alias_color;
    EmitVertex();
    gl_Position = proj_mat * point_a[1];
    frag_coord = point_a[1].xyz;
    frag_color = geom_color[0];
    EmitVertex();
    gl_Position = proj_mat * point_b[1];
    frag_coord = point_b[1].xyz;
    frag_color = geom_color[0];
    EmitVertex();
    gl_Position = proj_mat * geom_coord[0];
    frag_coord = geom_coord[0].xyz;
    frag_color = geom_color[0];
    EmitVertex();
    gl_Position = proj_mat * mid_point;
    frag_coord = mid_point.xyz;
    frag_color = geom_color[0];
    EmitVertex();
    gl_Position = proj_mat * point_a[2];
    frag_coord = point_a[2].xyz;
    frag_color = geom_color[0];
    EmitVertex();
    gl_Position = proj_mat * point_b[2];
    frag_coord = point_b[2].xyz;
    frag_color = geom_color[0];
    EmitVertex();
    gl_Position = proj_mat * point_a[3];
    frag_coord = point_a[3].xyz;
    frag_color = alias_color;
    EmitVertex();
    gl_Position = proj_mat * point_b[3];
    frag_coord = point_b[3].xyz;
    frag_color = alias_color;
    EmitVertex();
    EndPrimitive();
    
    gl_Position = proj_mat * point_b[0];
    frag_coord = point_b[0].xyz;
    frag_color = alias_color;
    EmitVertex();
    gl_Position = proj_mat * point_c[0];
    frag_coord = point_c[0].xyz;
    frag_color = alias_color;
    EmitVertex();
    gl_Position = proj_mat * point_b[1];
    frag_coord = point_b[1].xyz;
    frag_color = geom_color[1];
    EmitVertex();
    gl_Position = proj_mat * point_c[1];
    frag_coord = point_c[1].xyz;
    frag_color = geom_color[1];
    EmitVertex();
    gl_Position = proj_mat * mid_point;
    frag_coord = mid_point.xyz;
    frag_color = geom_color[1];
    EmitVertex();
    gl_Position = proj_mat * geom_coord[1];
    frag_coord = geom_coord[1].xyz;
    frag_color = geom_color[1];
    EmitVertex();
    gl_Position = proj_mat * point_b[2];
    frag_coord = point_b[2].xyz;
    frag_color = geom_color[1];
    EmitVertex();
    gl_Position = proj_mat * point_c[2];
    frag_coord = point_c[2].xyz;
    frag_color = geom_color[1];
    EmitVertex();
    gl_Position = proj_mat * point_b[3];
    frag_coord = point_b[3].xyz;
    frag_color = alias_color;
    EmitVertex();
    gl_Position = proj_mat * point_c[3];
    frag_coord = point_c[3].xyz;
    frag_color = alias_color;
    EmitVertex();
    EndPrimitive();
}
"""
fragment_shader_nlines = """
#version 330
precision highp float; 
precision highp int;

uniform vec4 fog_color;
uniform float fog_start;
uniform float fog_end;

in vec3 frag_coord;
in vec3 frag_color;

out vec4 final_color;

void main() {
    float dist = abs(frag_coord.z);
    float fog_factor = smoothstep(fog_end, fog_start, dist);
    final_color = mix(fog_color, vec4(frag_color, 1.0), fog_factor);
}
"""

################################## SELECTION ###################################
'''
sel_vertex_shader_nlines = """
#version 330

uniform mat4 model_mat;
uniform mat4 view_mat;

in vec3 vert_coord;
in vec3 vert_color;
const float vert_width = 0.15;

out vec3 geom_color;
out vec4 geom_coord;
out float geom_width;

void main(){
    geom_color = vert_color;
    geom_width = vert_width;
    geom_coord = view_mat * model_mat * vec4(vert_coord, 1.0);
}
"""
sel_geometry_shader_nlines = """
#version 330

layout (lines) in;
layout (triangle_strip, max_vertices = 20) out;

uniform mat4 proj_mat;
const float antialias_length = 0.1;

in vec3 geom_color[];
in vec4 geom_coord[];
in float geom_width[];

varying vec4 point_a[4];
varying vec4 point_b[4];
varying vec4 point_c[4];

out vec3 frag_coord;
out vec3 frag_color;

void main(){
    vec4 mid_point = (geom_coord[0] + geom_coord[1])/2;
    vec2 dir_vec = normalize((geom_coord[1] - geom_coord[0]).xy);
    vec3 orto_vec = normalize(cross(vec3(dir_vec, 0), vec3(0, 0, 1)));
    vec3 alias_color = vec3(0, 0, 0);
    point_a[0] = vec4(geom_coord[0].xyz + orto_vec * antialias_length, 1.0);
    point_a[1] = vec4(geom_coord[0].xyz + orto_vec * antialias_length/2, 1.0);
    point_a[2] = vec4(geom_coord[0].xyz - orto_vec * antialias_length/2, 1.0);
    point_a[3] = vec4(geom_coord[0].xyz - orto_vec * antialias_length, 1.0);
    point_b[0] = vec4(mid_point.xyz + orto_vec * antialias_length, 1.0);
    point_b[1] = vec4(mid_point.xyz + orto_vec * antialias_length/2, 1.0);
    point_b[2] = vec4(mid_point.xyz - orto_vec * antialias_length/2, 1.0);
    point_b[3] = vec4(mid_point.xyz - orto_vec * antialias_length, 1.0);
    point_c[0] = vec4(geom_coord[1].xyz + orto_vec * antialias_length, 1.0);
    point_c[1] = vec4(geom_coord[1].xyz + orto_vec * antialias_length/2, 1.0);
    point_c[2] = vec4(geom_coord[1].xyz - orto_vec * antialias_length/2, 1.0);
    point_c[3] = vec4(geom_coord[1].xyz - orto_vec * antialias_length, 1.0);
    
    gl_Position = proj_mat * point_a[0];
    frag_color = alias_color;
    EmitVertex();
    gl_Position = proj_mat * point_b[0];
    frag_color = alias_color;
    EmitVertex();
    gl_Position = proj_mat * point_a[1];
    frag_color = geom_color[0];
    EmitVertex();
    gl_Position = proj_mat * point_b[1];
    frag_color = geom_color[0];
    EmitVertex();
    gl_Position = proj_mat * geom_coord[0];
    frag_color = geom_color[0];
    EmitVertex();
    gl_Position = proj_mat * mid_point;
    frag_color = geom_color[0];
    EmitVertex();
    gl_Position = proj_mat * point_a[2];
    frag_color = geom_color[0];
    EmitVertex();
    gl_Position = proj_mat * point_b[2];
    frag_color = geom_color[0];
    EmitVertex();
    gl_Position = proj_mat * point_a[3];
    frag_color = alias_color;
    EmitVertex();
    gl_Position = proj_mat * point_b[3];
    frag_color = alias_color;
    EmitVertex();
    EndPrimitive();
    
    gl_Position = proj_mat * point_b[0];
    frag_color = alias_color;
    EmitVertex();
    gl_Position = proj_mat * point_c[0];
    frag_color = alias_color;
    EmitVertex();
    gl_Position = proj_mat * point_b[1];
    frag_color = geom_color[1];
    EmitVertex();
    gl_Position = proj_mat * point_c[1];
    frag_color = geom_color[1];
    EmitVertex();
    gl_Position = proj_mat * mid_point;
    frag_color = geom_color[1];
    EmitVertex();
    gl_Position = proj_mat * geom_coord[1];
    frag_color = geom_color[1];
    EmitVertex();
    gl_Position = proj_mat * point_b[2];
    frag_color = geom_color[1];
    EmitVertex();
    gl_Position = proj_mat * point_c[2];
    frag_color = geom_color[1];
    EmitVertex();
    gl_Position = proj_mat * point_b[3];
    frag_color = alias_color;
    EmitVertex();
    gl_Position = proj_mat * point_c[3];
    frag_color = alias_color;
    EmitVertex();
    EndPrimitive();
}
"""
sel_fragment_shader_nlines = """
#version 330

in vec3 frag_color;
out vec4 final_color;

void main(){
    final_color = vec4(frag_color, 1.0);
}
"""
'''










vertex_shader_sticks = """
#version 330
precision highp float; 
precision highp int;

uniform mat4 model_mat;
uniform mat4 view_mat;

in vec3 vert_coord;
in vec3 vert_color;
//const float vert_rad = 0.12;
const float vert_rad = 0.025;

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
precision highp float; 
precision highp int;

layout (lines) in;
layout (triangle_strip, max_vertices = 40) out;

uniform mat4 proj_mat;

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
varying mat3 rot_mat;

// mid_point is the middle point in the line.
varying vec3 mid_point;

// These points are the vertices calculated for the cylinder.
varying vec3 p_00, p_01, p_02, p_03, p_04, p_05, p_06, p_07, p_08;
varying vec3 p_09, p_10, p_11, p_12, p_13, p_14, p_15, p_16, p_17;
varying vec3 p_18, p_19, p_20, p_21, p_22, p_23, p_24, p_25, p_26;

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
    p_00.x = (bs_0.x*rot_mat[0][0] + bs_0.y*rot_mat[0][1] + bs_0.z*rot_mat[0][2])*geom_rad[0] + geom_coord[0].x;
    p_00.y = (bs_0.x*rot_mat[1][0] + bs_0.y*rot_mat[1][1] + bs_0.z*rot_mat[1][2])*geom_rad[0] + geom_coord[0].y;
    p_00.z = (bs_0.x*rot_mat[2][0] + bs_0.y*rot_mat[2][1] + bs_0.z*rot_mat[2][2])*geom_rad[0] + geom_coord[0].z;
    p_01.x = (bs_1.x*rot_mat[0][0] + bs_1.y*rot_mat[0][1] + bs_1.z*rot_mat[0][2])*geom_rad[0] + geom_coord[0].x;
    p_01.y = (bs_1.x*rot_mat[1][0] + bs_1.y*rot_mat[1][1] + bs_1.z*rot_mat[1][2])*geom_rad[0] + geom_coord[0].y;
    p_01.z = (bs_1.x*rot_mat[2][0] + bs_1.y*rot_mat[2][1] + bs_1.z*rot_mat[2][2])*geom_rad[0] + geom_coord[0].z;
    p_02.x = (bs_2.x*rot_mat[0][0] + bs_2.y*rot_mat[0][1] + bs_2.z*rot_mat[0][2])*geom_rad[0] + geom_coord[0].x;
    p_02.y = (bs_2.x*rot_mat[1][0] + bs_2.y*rot_mat[1][1] + bs_2.z*rot_mat[1][2])*geom_rad[0] + geom_coord[0].y;
    p_02.z = (bs_2.x*rot_mat[2][0] + bs_2.y*rot_mat[2][1] + bs_2.z*rot_mat[2][2])*geom_rad[0] + geom_coord[0].z;
    p_03.x = (bs_3.x*rot_mat[0][0] + bs_3.y*rot_mat[0][1] + bs_3.z*rot_mat[0][2])*geom_rad[0] + geom_coord[0].x;
    p_03.y = (bs_3.x*rot_mat[1][0] + bs_3.y*rot_mat[1][1] + bs_3.z*rot_mat[1][2])*geom_rad[0] + geom_coord[0].y;
    p_03.z = (bs_3.x*rot_mat[2][0] + bs_3.y*rot_mat[2][1] + bs_3.z*rot_mat[2][2])*geom_rad[0] + geom_coord[0].z;
    p_04.x = (bs_4.x*rot_mat[0][0] + bs_4.y*rot_mat[0][1] + bs_4.z*rot_mat[0][2])*geom_rad[0] + geom_coord[0].x;
    p_04.y = (bs_4.x*rot_mat[1][0] + bs_4.y*rot_mat[1][1] + bs_4.z*rot_mat[1][2])*geom_rad[0] + geom_coord[0].y;
    p_04.z = (bs_4.x*rot_mat[2][0] + bs_4.y*rot_mat[2][1] + bs_4.z*rot_mat[2][2])*geom_rad[0] + geom_coord[0].z;
    p_05.x = (bs_5.x*rot_mat[0][0] + bs_5.y*rot_mat[0][1] + bs_5.z*rot_mat[0][2])*geom_rad[0] + geom_coord[0].x;
    p_05.y = (bs_5.x*rot_mat[1][0] + bs_5.y*rot_mat[1][1] + bs_5.z*rot_mat[1][2])*geom_rad[0] + geom_coord[0].y;
    p_05.z = (bs_5.x*rot_mat[2][0] + bs_5.y*rot_mat[2][1] + bs_5.z*rot_mat[2][2])*geom_rad[0] + geom_coord[0].z;
    p_06.x = (bs_6.x*rot_mat[0][0] + bs_6.y*rot_mat[0][1] + bs_6.z*rot_mat[0][2])*geom_rad[0] + geom_coord[0].x;
    p_06.y = (bs_6.x*rot_mat[1][0] + bs_6.y*rot_mat[1][1] + bs_6.z*rot_mat[1][2])*geom_rad[0] + geom_coord[0].y;
    p_06.z = (bs_6.x*rot_mat[2][0] + bs_6.y*rot_mat[2][1] + bs_6.z*rot_mat[2][2])*geom_rad[0] + geom_coord[0].z;
    p_07.x = (bs_7.x*rot_mat[0][0] + bs_7.y*rot_mat[0][1] + bs_7.z*rot_mat[0][2])*geom_rad[0] + geom_coord[0].x;
    p_07.y = (bs_7.x*rot_mat[1][0] + bs_7.y*rot_mat[1][1] + bs_7.z*rot_mat[1][2])*geom_rad[0] + geom_coord[0].y;
    p_07.z = (bs_7.x*rot_mat[2][0] + bs_7.y*rot_mat[2][1] + bs_7.z*rot_mat[2][2])*geom_rad[0] + geom_coord[0].z;
    p_08.x = (bs_8.x*rot_mat[0][0] + bs_8.y*rot_mat[0][1] + bs_8.z*rot_mat[0][2])*geom_rad[0] + geom_coord[0].x;
    p_08.y = (bs_8.x*rot_mat[1][0] + bs_8.y*rot_mat[1][1] + bs_8.z*rot_mat[1][2])*geom_rad[0] + geom_coord[0].y;
    p_08.z = (bs_8.x*rot_mat[2][0] + bs_8.y*rot_mat[2][1] + bs_8.z*rot_mat[2][2])*geom_rad[0] + geom_coord[0].z;
    
    p_09.x = (up_0.x*rot_mat[0][0] + up_0.y*rot_mat[0][1] + up_0.z*rot_mat[0][2])*geom_rad[0] + mid_point.x;
    p_09.y = (up_0.x*rot_mat[1][0] + up_0.y*rot_mat[1][1] + up_0.z*rot_mat[1][2])*geom_rad[0] + mid_point.y;
    p_09.z = (up_0.x*rot_mat[2][0] + up_0.y*rot_mat[2][1] + up_0.z*rot_mat[2][2])*geom_rad[0] + mid_point.z;
    p_10.x = (up_1.x*rot_mat[0][0] + up_1.y*rot_mat[0][1] + up_1.z*rot_mat[0][2])*geom_rad[0] + mid_point.x;
    p_10.y = (up_1.x*rot_mat[1][0] + up_1.y*rot_mat[1][1] + up_1.z*rot_mat[1][2])*geom_rad[0] + mid_point.y;
    p_10.z = (up_1.x*rot_mat[2][0] + up_1.y*rot_mat[2][1] + up_1.z*rot_mat[2][2])*geom_rad[0] + mid_point.z;
    p_11.x = (up_2.x*rot_mat[0][0] + up_2.y*rot_mat[0][1] + up_2.z*rot_mat[0][2])*geom_rad[0] + mid_point.x;
    p_11.y = (up_2.x*rot_mat[1][0] + up_2.y*rot_mat[1][1] + up_2.z*rot_mat[1][2])*geom_rad[0] + mid_point.y;
    p_11.z = (up_2.x*rot_mat[2][0] + up_2.y*rot_mat[2][1] + up_2.z*rot_mat[2][2])*geom_rad[0] + mid_point.z;
    p_12.x = (up_3.x*rot_mat[0][0] + up_3.y*rot_mat[0][1] + up_3.z*rot_mat[0][2])*geom_rad[0] + mid_point.x;
    p_12.y = (up_3.x*rot_mat[1][0] + up_3.y*rot_mat[1][1] + up_3.z*rot_mat[1][2])*geom_rad[0] + mid_point.y;
    p_12.z = (up_3.x*rot_mat[2][0] + up_3.y*rot_mat[2][1] + up_3.z*rot_mat[2][2])*geom_rad[0] + mid_point.z;
    p_13.x = (up_4.x*rot_mat[0][0] + up_4.y*rot_mat[0][1] + up_4.z*rot_mat[0][2])*geom_rad[0] + mid_point.x;
    p_13.y = (up_4.x*rot_mat[1][0] + up_4.y*rot_mat[1][1] + up_4.z*rot_mat[1][2])*geom_rad[0] + mid_point.y;
    p_13.z = (up_4.x*rot_mat[2][0] + up_4.y*rot_mat[2][1] + up_4.z*rot_mat[2][2])*geom_rad[0] + mid_point.z;
    p_14.x = (up_5.x*rot_mat[0][0] + up_5.y*rot_mat[0][1] + up_5.z*rot_mat[0][2])*geom_rad[0] + mid_point.x;
    p_14.y = (up_5.x*rot_mat[1][0] + up_5.y*rot_mat[1][1] + up_5.z*rot_mat[1][2])*geom_rad[0] + mid_point.y;
    p_14.z = (up_5.x*rot_mat[2][0] + up_5.y*rot_mat[2][1] + up_5.z*rot_mat[2][2])*geom_rad[0] + mid_point.z;
    p_15.x = (up_6.x*rot_mat[0][0] + up_6.y*rot_mat[0][1] + up_6.z*rot_mat[0][2])*geom_rad[0] + mid_point.x;
    p_15.y = (up_6.x*rot_mat[1][0] + up_6.y*rot_mat[1][1] + up_6.z*rot_mat[1][2])*geom_rad[0] + mid_point.y;
    p_15.z = (up_6.x*rot_mat[2][0] + up_6.y*rot_mat[2][1] + up_6.z*rot_mat[2][2])*geom_rad[0] + mid_point.z;
    p_16.x = (up_7.x*rot_mat[0][0] + up_7.y*rot_mat[0][1] + up_7.z*rot_mat[0][2])*geom_rad[0] + mid_point.x;
    p_16.y = (up_7.x*rot_mat[1][0] + up_7.y*rot_mat[1][1] + up_7.z*rot_mat[1][2])*geom_rad[0] + mid_point.y;
    p_16.z = (up_7.x*rot_mat[2][0] + up_7.y*rot_mat[2][1] + up_7.z*rot_mat[2][2])*geom_rad[0] + mid_point.z;
    p_17.x = (up_8.x*rot_mat[0][0] + up_8.y*rot_mat[0][1] + up_8.z*rot_mat[0][2])*geom_rad[0] + mid_point.x;
    p_17.y = (up_8.x*rot_mat[1][0] + up_8.y*rot_mat[1][1] + up_8.z*rot_mat[1][2])*geom_rad[0] + mid_point.y;
    p_17.z = (up_8.x*rot_mat[2][0] + up_8.y*rot_mat[2][1] + up_8.z*rot_mat[2][2])*geom_rad[0] + mid_point.z;
    
    p_18.x = (bs_0.x*rot_mat[0][0] + bs_0.y*rot_mat[0][1] + bs_0.z*rot_mat[0][2])*geom_rad[0] + geom_coord[1].x;
    p_18.y = (bs_0.x*rot_mat[1][0] + bs_0.y*rot_mat[1][1] + bs_0.z*rot_mat[1][2])*geom_rad[0] + geom_coord[1].y;
    p_18.z = (bs_0.x*rot_mat[2][0] + bs_0.y*rot_mat[2][1] + bs_0.z*rot_mat[2][2])*geom_rad[0] + geom_coord[1].z;
    p_19.x = (bs_1.x*rot_mat[0][0] + bs_1.y*rot_mat[0][1] + bs_1.z*rot_mat[0][2])*geom_rad[0] + geom_coord[1].x;
    p_19.y = (bs_1.x*rot_mat[1][0] + bs_1.y*rot_mat[1][1] + bs_1.z*rot_mat[1][2])*geom_rad[0] + geom_coord[1].y;
    p_19.z = (bs_1.x*rot_mat[2][0] + bs_1.y*rot_mat[2][1] + bs_1.z*rot_mat[2][2])*geom_rad[0] + geom_coord[1].z;
    p_20.x = (bs_2.x*rot_mat[0][0] + bs_2.y*rot_mat[0][1] + bs_2.z*rot_mat[0][2])*geom_rad[0] + geom_coord[1].x;
    p_20.y = (bs_2.x*rot_mat[1][0] + bs_2.y*rot_mat[1][1] + bs_2.z*rot_mat[1][2])*geom_rad[0] + geom_coord[1].y;
    p_20.z = (bs_2.x*rot_mat[2][0] + bs_2.y*rot_mat[2][1] + bs_2.z*rot_mat[2][2])*geom_rad[0] + geom_coord[1].z;
    p_21.x = (bs_3.x*rot_mat[0][0] + bs_3.y*rot_mat[0][1] + bs_3.z*rot_mat[0][2])*geom_rad[0] + geom_coord[1].x;
    p_21.y = (bs_3.x*rot_mat[1][0] + bs_3.y*rot_mat[1][1] + bs_3.z*rot_mat[1][2])*geom_rad[0] + geom_coord[1].y;
    p_21.z = (bs_3.x*rot_mat[2][0] + bs_3.y*rot_mat[2][1] + bs_3.z*rot_mat[2][2])*geom_rad[0] + geom_coord[1].z;
    p_22.x = (bs_4.x*rot_mat[0][0] + bs_4.y*rot_mat[0][1] + bs_4.z*rot_mat[0][2])*geom_rad[0] + geom_coord[1].x;
    p_22.y = (bs_4.x*rot_mat[1][0] + bs_4.y*rot_mat[1][1] + bs_4.z*rot_mat[1][2])*geom_rad[0] + geom_coord[1].y;
    p_22.z = (bs_4.x*rot_mat[2][0] + bs_4.y*rot_mat[2][1] + bs_4.z*rot_mat[2][2])*geom_rad[0] + geom_coord[1].z;
    p_23.x = (bs_5.x*rot_mat[0][0] + bs_5.y*rot_mat[0][1] + bs_5.z*rot_mat[0][2])*geom_rad[0] + geom_coord[1].x;
    p_23.y = (bs_5.x*rot_mat[1][0] + bs_5.y*rot_mat[1][1] + bs_5.z*rot_mat[1][2])*geom_rad[0] + geom_coord[1].y;
    p_23.z = (bs_5.x*rot_mat[2][0] + bs_5.y*rot_mat[2][1] + bs_5.z*rot_mat[2][2])*geom_rad[0] + geom_coord[1].z;
    p_24.x = (bs_6.x*rot_mat[0][0] + bs_6.y*rot_mat[0][1] + bs_6.z*rot_mat[0][2])*geom_rad[0] + geom_coord[1].x;
    p_24.y = (bs_6.x*rot_mat[1][0] + bs_6.y*rot_mat[1][1] + bs_6.z*rot_mat[1][2])*geom_rad[0] + geom_coord[1].y;
    p_24.z = (bs_6.x*rot_mat[2][0] + bs_6.y*rot_mat[2][1] + bs_6.z*rot_mat[2][2])*geom_rad[0] + geom_coord[1].z;
    p_25.x = (bs_7.x*rot_mat[0][0] + bs_7.y*rot_mat[0][1] + bs_7.z*rot_mat[0][2])*geom_rad[0] + geom_coord[1].x;
    p_25.y = (bs_7.x*rot_mat[1][0] + bs_7.y*rot_mat[1][1] + bs_7.z*rot_mat[1][2])*geom_rad[0] + geom_coord[1].y;
    p_25.z = (bs_7.x*rot_mat[2][0] + bs_7.y*rot_mat[2][1] + bs_7.z*rot_mat[2][2])*geom_rad[0] + geom_coord[1].z;
    p_26.x = (bs_8.x*rot_mat[0][0] + bs_8.y*rot_mat[0][1] + bs_8.z*rot_mat[0][2])*geom_rad[0] + geom_coord[1].x;
    p_26.y = (bs_8.x*rot_mat[1][0] + bs_8.y*rot_mat[1][1] + bs_8.z*rot_mat[1][2])*geom_rad[0] + geom_coord[1].y;
    p_26.z = (bs_8.x*rot_mat[2][0] + bs_8.y*rot_mat[2][1] + bs_8.z*rot_mat[2][2])*geom_rad[0] + geom_coord[1].z;
}

void main(){
    mid_point = (geom_coord[0].xyz + geom_coord[1].xyz)/2;
    // vec_p0_p1 is the vector defined by the line.
    vec3 vec_p0_p1 = geom_coord[1].xyz - geom_coord[0].xyz;
    // ort_vec is the orthogonal vector between the line vector and the Y axis.
    vec3 ort_vec = normalize(cross(vec3(0,1,0), vec_p0_p1));
    // g_angle is the angle between the line vector and the Y axis.
    float g_angle = get_angle(vec3(0,1,0), vec_p0_p1);
    // g_length is the line vector length or simply the line length.
    float g_length = length(vec_p0_p1);
    rot_mat = get_rot_mat(g_angle, ort_vec);
    calculate_points();
    // Now we send the vertices to the fragment shader in a defined order
    // base-> 0, 9, 1, 10, 2, 11, 3, 12, 4, 13, 5, 14, 6, 15, 7, 16, 8, 17, 0, 9
    gl_Position = proj_mat * vec4(p_00, 1);
    frag_coord = p_00;
    frag_color = geom_color[0];
    frag_norm = p_00 - geom_coord[0].xyz;
    EmitVertex();
    gl_Position = proj_mat * vec4(p_09, 1);
    frag_coord = p_09;
    frag_color = geom_color[0];
    frag_norm = p_09 - mid_point;
    EmitVertex();
    gl_Position = proj_mat * vec4(p_01, 1);
    frag_coord = p_01;
    frag_color = geom_color[0];
    frag_norm = p_01 - geom_coord[0].xyz;
    EmitVertex();
    gl_Position = proj_mat * vec4(p_10, 1);
    frag_coord = p_10;
    frag_color = geom_color[0];
    frag_norm = p_10 - mid_point;
    EmitVertex();
    gl_Position = proj_mat * vec4(p_02, 1);
    frag_coord = p_02;
    frag_color = geom_color[0];
    frag_norm = p_02 - geom_coord[0].xyz;
    EmitVertex();
    gl_Position = proj_mat * vec4(p_11, 1);
    frag_coord = p_11;
    frag_color = geom_color[0];
    frag_norm = p_11 - mid_point;
    EmitVertex();
    gl_Position = proj_mat * vec4(p_03, 1);
    frag_coord = p_03;
    frag_color = geom_color[0];
    frag_norm = p_03 - geom_coord[0].xyz;
    EmitVertex();
    gl_Position = proj_mat * vec4(p_12, 1);
    frag_coord = p_12;
    frag_color = geom_color[0];
    frag_norm = p_12 - mid_point;
    EmitVertex();
    gl_Position = proj_mat * vec4(p_04, 1);
    frag_coord = p_04;
    frag_color = geom_color[0];
    frag_norm = p_04 - geom_coord[0].xyz;
    EmitVertex();
    gl_Position = proj_mat * vec4(p_13, 1);
    frag_coord = p_13;
    frag_color = geom_color[0];
    frag_norm = p_13 - mid_point;
    EmitVertex();
    gl_Position = proj_mat * vec4(p_05, 1);
    frag_coord = p_05;
    frag_color = geom_color[0];
    frag_norm = p_05 - geom_coord[0].xyz;
    EmitVertex();
    gl_Position = proj_mat * vec4(p_14, 1);
    frag_coord = p_14;
    frag_color = geom_color[0];
    frag_norm = p_14 - mid_point;
    EmitVertex();
    gl_Position = proj_mat * vec4(p_06, 1);
    frag_coord = p_06;
    frag_color = geom_color[0];
    frag_norm = p_06 - geom_coord[0].xyz;
    EmitVertex();
    gl_Position = proj_mat * vec4(p_15, 1);
    frag_coord = p_15;
    frag_color = geom_color[0];
    frag_norm = p_15 - mid_point;
    EmitVertex();
    gl_Position = proj_mat * vec4(p_07, 1);
    frag_coord = p_07;
    frag_color = geom_color[0];
    frag_norm = p_07 - geom_coord[0].xyz;
    EmitVertex();
    gl_Position = proj_mat * vec4(p_16, 1);
    frag_coord = p_16;
    frag_color = geom_color[0];
    frag_norm = p_16 - mid_point;
    EmitVertex();
    gl_Position = proj_mat * vec4(p_08, 1);
    frag_coord = p_08;
    frag_color = geom_color[0];
    frag_norm = p_08 - geom_coord[0].xyz;
    EmitVertex();
    gl_Position = proj_mat * vec4(p_17, 1);
    frag_coord = p_17;
    frag_color = geom_color[0];
    frag_norm = p_17 - mid_point;
    EmitVertex();
    gl_Position = proj_mat * vec4(p_00, 1);
    frag_coord = p_00;
    frag_color = geom_color[0];
    frag_norm = p_00 - geom_coord[0].xyz;
    EmitVertex();
    gl_Position = proj_mat * vec4(p_09, 1);
    frag_coord = p_09;
    frag_color = geom_color[0];
    frag_norm = p_09 - mid_point;
    EmitVertex();
    EndPrimitive();
    
    // up-> 9, 18, 10, 19, 11, 20, 12, 21, 13, 22, 14, 23, 15, 24, 16, 25, 17, 26, 9, 18
    gl_Position = proj_mat * vec4(p_09, 1);
    frag_coord = p_09;
    frag_color = geom_color[1];
    frag_norm = p_09 - mid_point;
    EmitVertex();
    gl_Position = proj_mat * vec4(p_18, 1);
    frag_coord = p_18;
    frag_color = geom_color[1];
    frag_norm = p_18 - geom_coord[1].xyz;
    EmitVertex();
    gl_Position = proj_mat * vec4(p_10, 1);
    frag_coord = p_10;
    frag_color = geom_color[1];
    frag_norm = p_10 - mid_point;
    EmitVertex();
    gl_Position = proj_mat * vec4(p_19, 1);
    frag_coord = p_19;
    frag_color = geom_color[1];
    frag_norm = p_19 - geom_coord[1].xyz;
    EmitVertex();
    gl_Position = proj_mat * vec4(p_11, 1);
    frag_coord = p_11;
    frag_color = geom_color[1];
    frag_norm = p_11 - mid_point;
    EmitVertex();
    gl_Position = proj_mat * vec4(p_20, 1);
    frag_coord = p_20;
    frag_color = geom_color[1];
    frag_norm = p_20 - geom_coord[1].xyz;
    EmitVertex();
    gl_Position = proj_mat * vec4(p_12, 1);
    frag_coord = p_12;
    frag_color = geom_color[1];
    frag_norm = p_12 - mid_point;
    EmitVertex();
    gl_Position = proj_mat * vec4(p_21, 1);
    frag_coord = p_21;
    frag_color = geom_color[1];
    frag_norm = p_21 - geom_coord[1].xyz;
    EmitVertex();
    gl_Position = proj_mat * vec4(p_13, 1);
    frag_coord = p_13;
    frag_color = geom_color[1];
    frag_norm = p_13 - mid_point;
    EmitVertex();
    gl_Position = proj_mat * vec4(p_22, 1);
    frag_coord = p_22;
    frag_color = geom_color[1];
    frag_norm = p_22 - geom_coord[1].xyz;
    EmitVertex();
    gl_Position = proj_mat * vec4(p_14, 1);
    frag_coord = p_14;
    frag_color = geom_color[1];
    frag_norm = p_14 - mid_point;
    EmitVertex();
    gl_Position = proj_mat * vec4(p_23, 1);
    frag_coord = p_23;
    frag_color = geom_color[1];
    frag_norm = p_23 - geom_coord[1].xyz;
    EmitVertex();
    gl_Position = proj_mat * vec4(p_15, 1);
    frag_coord = p_15;
    frag_color = geom_color[1];
    frag_norm = p_15 - mid_point;
    EmitVertex();
    gl_Position = proj_mat * vec4(p_24, 1);
    frag_coord = p_24;
    frag_color = geom_color[1];
    frag_norm = p_24 - geom_coord[1].xyz;
    EmitVertex();
    gl_Position = proj_mat * vec4(p_16, 1);
    frag_coord = p_16;
    frag_color = geom_color[1];
    frag_norm = p_16 - mid_point;
    EmitVertex();
    gl_Position = proj_mat * vec4(p_25, 1);
    frag_coord = p_25;
    frag_color = geom_color[1];
    frag_norm = p_25 - geom_coord[1].xyz;
    EmitVertex();
    gl_Position = proj_mat * vec4(p_17, 1);
    frag_coord = p_17;
    frag_color = geom_color[1];
    frag_norm = p_17 - mid_point;
    EmitVertex();
    gl_Position = proj_mat * vec4(p_26, 1);
    frag_coord = p_26;
    frag_color = geom_color[1];
    frag_norm = p_26 - geom_coord[1].xyz;
    EmitVertex();
    gl_Position = proj_mat * vec4(p_09, 1);
    frag_coord = p_09;
    frag_color = geom_color[1];
    frag_norm = p_09 - mid_point;
    EmitVertex();
    gl_Position = proj_mat * vec4(p_18, 1);
    frag_coord = p_18;
    frag_color = geom_color[1];
    frag_norm = p_18 - geom_coord[1].xyz;
    EmitVertex();
    EndPrimitive();
}
"""
fragment_shader_sticks = """
#version 330
precision highp float; 
precision highp int;

struct Light {
   vec3 position;
   //vec3 color;
   vec3 intensity;
   //vec3 specular_color;
   float ambient_coef;
   float shininess;
};

uniform Light my_light;

uniform vec4 fog_color;
uniform float fog_start;
uniform float fog_end;

in vec3 frag_coord;
in vec3 frag_color;
in vec3 frag_norm;

out vec4 final_color;

void main(){
    vec3 normal = normalize(frag_norm);
    vec3 vert_to_light = normalize(my_light.position);
    vec3 vert_to_cam = normalize(frag_coord);
    
    // Ambient Component
    vec3 ambient = my_light.ambient_coef * frag_color * my_light.intensity;
    
    // Diffuse component
    float diffuse_coef = max(0.0, dot(normal, vert_to_light));
    vec3 diffuse = diffuse_coef * frag_color * my_light.intensity;
    
    // Specular component
    float specular_coef = 0.0;
    if (diffuse_coef > 0.0)
        specular_coef = pow(max(0.0, dot(vert_to_cam, reflect(vert_to_light, normal))), my_light.shininess);
    vec3 specular = specular_coef * my_light.intensity;
    specular = specular * (vec3(1) - diffuse);
    vec4 my_color = vec4(ambient + diffuse + specular, 1.0);
    
    float dist = abs(frag_coord.z);
    if(dist>=fog_start){
        float fog_factor = (fog_end-dist)/(fog_end-fog_start);
        final_color = mix(fog_color, my_color, fog_factor);
    }
    else{
       final_color = my_color;
    }
}
"""











shader_type ={
            0: { "vertex_shader"      : vertex_shader_lines  ,
                 "geometry_shader"    : geometry_shader_lines,
                 "fragment_shader"    : fragment_shader_lines,
                 "sel_vertex_shader"  : sel_vertex_shader_lines  ,
                 "sel_geometry_shader": sel_geometry_shader_lines,
                 "sel_fragment_shader": sel_fragment_shader_lines

               },
            
            
            1: {"vertex_shader"      : vertex_shader_nlines  ,
                "geometry_shader"    : geometry_shader_nlines,
                "fragment_shader"    : fragment_shader_nlines,
                "sel_vertex_shader"  : sel_vertex_shader_lines  ,
                "sel_geometry_shader": sel_geometry_shader_lines,
                "sel_fragment_shader": sel_fragment_shader_lines
                },
            
            2: {"vertex_shader"      : vertex_shader_sticks  ,
                "geometry_shader"    : geometry_shader_sticks,
                "fragment_shader"    : fragment_shader_sticks,
                "sel_vertex_shader"  : sel_vertex_shader_lines  ,
                "sel_geometry_shader": sel_geometry_shader_lines,
                "sel_fragment_shader": sel_fragment_shader_lines
                },
}
