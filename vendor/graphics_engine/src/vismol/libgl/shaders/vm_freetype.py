#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#

vertex_shader_freetype =  """
#version 330

layout(std140) uniform CameraMatrices {
    mat4 view_mat;
    mat4 proj_mat;
};


mat4 identityMatrix = mat4(1.0, 0.0, 0.0, 0.0,
                           0.0, 1.0, 0.0, 0.0,
                           0.0, 0.0, 1.0, 0.0,
                           0.0, 0.0, 0.0, 1.0);


in vec3 vert_coord;
in vec4 vert_uv;
// [EN] Character slot index (0, 1, 2, ...) of this glyph inside its
// string. Every character of a given label shares the SAME vert_coord
// (the string's anchor point, in world/model space) -- the actual
// per-character horizontal advance is no longer baked into vert_coord
// on the CPU (that used to displace each glyph along the WORLD X axis,
// which only looks right when the camera happens to be looking exactly
// down -Z; any other camera orientation skewed/sheared the whole
// label). Instead we pass the anchor once and let the geometry shader
// compute the advance along the camera's own screen-space right/up
// axes, so labels always read as straight, upright, camera-facing text
// regardless of view orientation. See geometry_shader_freetype below.
in float vert_char_idx;

out vec4 geom_coord;
out vec4 geom_text_uv;
out float geom_char_idx;

void main(){
    //geom_coord = identityMatrix* vec4(vert_coord, 1.0);
    geom_coord = view_mat * vec4(vert_coord, 1.0);
    geom_text_uv = vert_uv;
    geom_char_idx = vert_char_idx;
}
"""
geometry_shader_freetype = """
#version 330

layout(std140) uniform CameraMatrices {
    mat4 view_mat;
    mat4 proj_mat;
};

layout (points) in;
layout (triangle_strip, max_vertices = 4) out;

// [EN] "offset" is the glyph quad half-size (half-width/half-height),
// "char_advance" is the horizontal step between consecutive glyphs of
// the same string, and "string_shift" is a small constant (x, y) nudge
// -- in character-size units -- applied to the whole string (used to
// center/offset a label next to the atom/dot it names). All three are
// calibrated at camera distance == depth_ref (see calculate_points()).
uniform vec2 offset;
uniform float char_advance;
uniform vec2 string_shift;
// [EN] Reference camera distance (world units) at which "offset" and
// "char_advance" are exactly the sizes given -- i.e. how big the label
// looks when the camera is depth_ref units away. Kept as a uniform
// (rather than a shader constant) so every VismolFont can share the
// same Python-side calibration value (see LABEL_DEPTH_REFERENCE in
// vismol_font.py).
uniform float depth_ref;
// [EN] Blends between two behaviors for how label size responds to
// camera distance ("zoom"):
//   0.0 -- constant size on screen, independent of camera distance
//          (the default since the billboard refactor: labels never
//          shrink when you zoom/dolly out, nor grow when you zoom in).
//   1.0 -- natural perspective size, i.e. the ORIGINAL (pre-refactor)
//          behavior: the label is a fixed size in world space, so it
//          shrinks as the camera moves away and grows as it gets
//          closer, exactly like any other piece of 3D geometry.
// Values in between smoothly blend the two. See calculate_points().
uniform float zoom_sensitivity;

in vec4 geom_coord[];
in vec4 geom_text_uv[];
in float geom_char_idx[];

out vec2 frag_text_uv;

void calculate_uv(in vec4 coord, out vec2 uvA, out vec2 uvB, out vec2 uvC, out vec2 uvD){
    // Taking the data from the upper-left and lower-right points of the quad,
    // it generates the coordinates of two triangles to form the letter.
    // The quad is created using the following pattern:
    //                                
    // uvA       uvD     uvA      uvD 
    //  |         |       |     /  |  
    //  |         |  -->  |   /    |  
    //  |         |       | /      |  
    // uvB-------uvC     uvB      uvC 
    //                                
    uvA = vec2(coord.xy);
    uvB = vec2(coord.xw);
    uvC = vec2(coord.zw);
    uvD = vec2(coord.zy);
}

void calculate_points(in vec4 coord, in float char_idx, out vec4 pA, out vec4 pB, out vec4 pC, out vec4 pD){
    // Creates the coordinates for a quad using the coord vector as center and
    // the xyz_offset as margins, defining the letter size. You can change the
    // xyz_value to get a bigger letter, but this will reduce the resolution.
    // Using # as a coordinate example, the quad is constructed as:
    //                                
    //                \|/         ┌--┐ 
    //      #   -->   -#-   -->   |  | 
    //                /|\         └--┘ 
    //                                
    //
    // [EN] "coord" is already in VIEW space (see vertex shader). In view
    // space the X/Y axes ARE the camera's screen-space right/up axes no
    // matter how the camera is rotated, so any offset we add to coord.x/
    // coord.y here (before the perspective projection below) is, by
    // construction, screen-aligned -- this is the classic view-space
    // billboard trick, and it's what makes the quad always face the
    // camera. Both the per-character advance (char_idx * char_advance)
    // and the small per-string nudge (string_shift) are computed here,
    // in the SAME screen-aligned space, instead of on the CPU in world
    // space -- that's what keeps a whole label reading as straight,
    // upright text instead of shearing when the camera isn't looking
    // straight down -Z.
    //
    // "depth" is the (positive) distance from the camera to this point,
    // in view space -Z. The projection matrix divides X/Y by (roughly)
    // this same depth during the perspective divide, which is exactly
    // why distant geometry looks smaller.
    //
    // depth_factor is what we actually multiply offset/char_advance by.
    // It's designed so that, after the perspective divide (roughly
    // "size / depth"), the resulting apparent size is:
    //   zoom_sensitivity == 0  ->  offset / depth_ref            (CONSTANT,
    //       independent of the camera's actual distance -- labels never
    //       shrink/grow with zoom.)
    //   zoom_sensitivity == 1  ->  offset / depth                (matches
    //       plain, uncompensated perspective -- exactly how these labels
    //       behaved before the billboard refactor, and how any other
    //       piece of 3D geometry scales with distance.)
    // Both cases (and everything in between) agree at depth == depth_ref,
    // which is what keeps "offset"/"char_advance" meaning the same
    // physical on-screen size regardless of zoom_sensitivity.
    float depth = -coord.z;
    float s = clamp(zoom_sensitivity, 0.0, 1.0);
    float depth_factor = pow(depth / depth_ref, 1.0 - s);
    vec2 scaled_offset = offset * depth_factor;
    float advance = (char_idx + string_shift.x) * char_advance * depth_factor;
    float center_x = coord.x + advance;
    float center_y = coord.y + string_shift.y * char_advance * depth_factor;
    pA = vec4(center_x - scaled_offset.x, center_y + scaled_offset.y, coord.z, 1.0);
    pB = vec4(center_x - scaled_offset.x, center_y - scaled_offset.y, coord.z, 1.0);
    pC = vec4(center_x + scaled_offset.x, center_y - scaled_offset.y, coord.z, 1.0);
    pD = vec4(center_x + scaled_offset.x, center_y + scaled_offset.y, coord.z, 1.0);
}

void main(){
    vec2 textA, textB, textC, textD;
    vec4 pointA, pointB, pointC, pointD;
    calculate_uv(geom_text_uv[0], textA, textB, textC, textD);
    calculate_points(geom_coord[0], geom_char_idx[0], pointA, pointB, pointC, pointD);
    gl_Position = proj_mat * pointA; frag_text_uv = textA; EmitVertex();
    gl_Position = proj_mat * pointB; frag_text_uv = textB; EmitVertex();
    gl_Position = proj_mat * pointD; frag_text_uv = textD; EmitVertex();
    gl_Position = proj_mat * pointC; frag_text_uv = textC; EmitVertex();
    
    //gl_Position =  pointA; frag_text_uv = textA; EmitVertex();
    //gl_Position =  pointB; frag_text_uv = textB; EmitVertex();
    //gl_Position =  pointD; frag_text_uv = textD; EmitVertex();
    //gl_Position =  pointC; frag_text_uv = textC; EmitVertex();
    
    EndPrimitive();
}
"""
fragment_shader_freetype = """
#version 330

uniform sampler2D textu;
uniform vec4 text_color;

//uniform float border_size
const float border_size = 0.05;

in vec2 frag_text_uv;

out vec4 final_color;

void main(){
    vec4 sampled = vec4(1.0, 1.0, 1.0, texture(textu, frag_text_uv).r);
    if (sampled.a==0.0)
        discard;
    final_color = text_color * sampled;
}



//void main() {
//    // Amostra a textura do caractere
//    vec4 sampled = texture(textu, frag_text_uv);
//    
//    // Verifica se o fragmento está na borda do caractere
//    float border_distance = fwidth(length(frag_text_uv - 0.5)); // Distância do fragmento ao centro do caractere
//    float alpha = sampled.a;
//    float border_alpha = smoothstep(0.5 - border_size, 0.5, border_distance) - smoothstep(0.5, 0.5 + border_size, border_distance);
//    
//    // Combina a cor do caractere com a cor da borda preta
//    vec4 border_color = vec4(0.0, 0.0, 0.0, 1.0); // Cor preta
//    vec4 final_alpha = mix(sampled, border_color, border_alpha);
//    
//    // Combina a cor do caractere com a cor do texto e a cor da borda
//    final_color = (text_color * final_alpha);
//}

//void main() {
//    // Amostra a textura do caractere
//    vec4 sampled = texture(textu, frag_text_uv);
//    
//    // Calcula a dilatação da textura do caractere com a borda preta
//    float border_distance = fwidth(length(frag_text_uv - 0.5)); // Distância do fragmento ao centro do caractere
//    float dilated_alpha = 1.0 - smoothstep(0.5 - border_size, 0.5 + border_size, border_distance);
//    
//    // Calcula a cor final da borda
//    vec4 border_color = vec4(0.0, 0.0, 0.0, 1.0); // Cor preta
//    vec4 border = dilated_alpha * border_color;
//    
//    // Calcula a cor final do caractere com a borda preta
//    vec4 final_alpha = mix(border, sampled, step(0.0, sampled.a));
//    
//    // Combina a cor do caractere com a cor do texto e a cor da borda
//    final_color = text_color * final_alpha;
//}

//void main(){
//    vec4 sampled = vec4(1.0, 1.0, 1.0, texture(textu, frag_text_uv).r); // Amostra a cor do texto
//    if (sampled.a == 0.0) // Se o fragmento for transparente, descartamos
//        discard;
//
//    // Calcula a distância do fragmento ao centro da textura
//    vec2 center = vec2(0.5, 0.5);
//    float distance_to_center = distance(frag_text_uv, center);
//
//    // Calcula a cor da borda (semi-transparente cinza)
//    vec4 border_color = vec4(0.5, 0.5, 0.5, 0.5); // Cinza semi-transparente
//
//    // Se o fragmento estiver dentro do círculo circunscrito, usa a cor da borda
//    if (distance_to_center > 0.5 - border_size && distance_to_center < 0.5) {
//        final_color = border_color;
//    } else {
//        final_color = text_color * sampled; // Caso contrário, usa a cor do texto original
//    }
//}

"""



static_vertex_shader_freetype =  """
#version 330

layout(std140) uniform CameraMatrices {
    mat4 view_mat;
    mat4 proj_mat;
};


mat4 identityMatrix = mat4(1.0, 0.0, 0.0, 0.0,
                           0.0, 1.0, 0.0, 0.0,
                           0.0, 0.0, 1.0, 0.0,
                           0.0, 0.0, 0.0, 1.0);


in vec3 vert_coord;
in vec4 vert_uv;

out vec4 geom_coord;
out vec4 geom_text_uv;

void main(){
    geom_coord = identityMatrix* vec4(vert_coord, 1.0);
    //geom_coord = view_mat * vec4(vert_coord, 1.0);
    geom_text_uv = vert_uv;
}
"""
static_geometry_shader_freetype = """
#version 330

layout(std140) uniform CameraMatrices {
    mat4 view_mat;
    mat4 proj_mat;
};

layout (points) in;
layout (triangle_strip, max_vertices = 4) out;

uniform vec2 offset;

in vec4 geom_coord[];
in vec4 geom_text_uv[];

out vec2 frag_text_uv;

void calculate_uv(in vec4 coord, out vec2 uvA, out vec2 uvB, out vec2 uvC, out vec2 uvD){
    // Taking the data from the upper-left and lower-right points of the quad,
    // it generates the coordinates of two triangles to form the letter.
    // The quad is created using the following pattern:
    //                                
    // uvA       uvD     uvA      uvD 
    //  |         |       |     /  |  
    //  |         |  -->  |   /    |  
    //  |         |       | /      |  
    // uvB-------uvC     uvB      uvC 
    //                                
    uvA = vec2(coord.xy);
    uvB = vec2(coord.xw);
    uvC = vec2(coord.zw);
    uvD = vec2(coord.zy);
}

void calculate_points(in vec4 coord, out vec4 pA, out vec4 pB, out vec4 pC, out vec4 pD){
    // Creates the coordinates for a quad using the coord vector as center and
    // the xyz_offset as margins, defining the letter size. You can change the
    // xyz_value to get a bigger letter, but this will reduce the resolution.
    // Using # as a coordinate example, the quad is constructed as:
    //                                
    //                \|/         ┌--┐ 
    //      #   -->   -#-   -->   |  | 
    //                /|\         └--┘ 
    //                                
    pA = vec4(coord.x - offset.x, coord.y + offset.y, coord.z, 1.0);
    pB = vec4(coord.x - offset.x, coord.y - offset.y, coord.z, 1.0);
    pC = vec4(coord.x + offset.x, coord.y - offset.y, coord.z, 1.0);
    pD = vec4(coord.x + offset.x, coord.y + offset.y, coord.z, 1.0);
}

void main(){
    vec2 textA, textB, textC, textD;
    vec4 pointA, pointB, pointC, pointD;
    calculate_uv(geom_text_uv[0], textA, textB, textC, textD);
    calculate_points(geom_coord[0], pointA, pointB, pointC, pointD);
    //gl_Position = proj_mat * pointA; frag_text_uv = textA; EmitVertex();
    //gl_Position = proj_mat * pointB; frag_text_uv = textB; EmitVertex();
    //gl_Position = proj_mat * pointD; frag_text_uv = textD; EmitVertex();
    //gl_Position = proj_mat * pointC; frag_text_uv = textC; EmitVertex();
    
    gl_Position =  pointA; frag_text_uv = textA; EmitVertex();
    gl_Position =  pointB; frag_text_uv = textB; EmitVertex();
    gl_Position =  pointD; frag_text_uv = textD; EmitVertex();
    gl_Position =  pointC; frag_text_uv = textC; EmitVertex();
    
    EndPrimitive();
}
"""
static_fragment_shader_freetype = """
#version 330

uniform sampler2D textu;
uniform vec4 text_color;

in vec2 frag_text_uv;

out vec4 final_color;

void main(){
    vec4 sampled = vec4(1.0, 1.0, 1.0, texture(textu, frag_text_uv).r);
    if (sampled.a==0.0)
        discard;
    final_color = text_color * sampled;
}
"""








v_shader_freetype =  """
#version 330

layout(std140) uniform CameraMatrices {
    mat4 view_mat;
    mat4 proj_mat;
};

uniform mat4 model_mat;

in vec3 vert_coord;
in vec4 vert_uv;

out vec4 geom_coord;
out vec4 geom_text_uv;

void main(){
    geom_coord = model_mat * vec4(vert_coord.xy, 0.0, 1.0);
    //geom_coord = view_mat * model_mat * vec4(vert_coord.xy, 0.0, 1.0);
    geom_text_uv = vert_uv;
}
"""
g_shader_freetype = """
#version 330

layout(std140) uniform CameraMatrices {
    mat4 view_mat;
    mat4 proj_mat;
};

layout (points) in;
layout (triangle_strip, max_vertices = 4) out;

uniform vec2 offset;

in vec4 geom_coord[];
in vec4 geom_text_uv[];

out vec2 frag_text_uv;

void calculate_uv(in vec4 coord, out vec2 uvA, out vec2 uvB, out vec2 uvC, out vec2 uvD){
    // Taking the data from the upper-left and lower-right points of the quad,
    // it generates the coordinates of two triangles to form the letter.
    // The quad is created using the following pattern:
    //                                
    // uvA       uvD     uvA      uvD 
    //  |         |       |     /  |  
    //  |         |  -->  |   /    |  
    //  |         |       | /      |  
    // uvB-------uvC     uvB      uvC 
    //                                
    uvA = vec2(coord.xy);
    uvB = vec2(coord.xw);
    uvC = vec2(coord.zw);
    uvD = vec2(coord.zy);
}

void calculate_points(in vec4 coord, out vec4 pA, out vec4 pB, out vec4 pC, out vec4 pD){
    // Creates the coordinates for a quad using the coord vector as center and
    // the xyz_offset as margins, defining the letter size. You can change the
    // xyz_value to get a bigger letter, but this will reduce the resolution.
    // Using # as a coordinate example, the quad is constructed as:
    //                                
    //                \|/         ┌--┐ 
    //      #   -->   -#-   -->   |  | 
    //                /|\         └--┘ 
    //                                
    pA = vec4(coord.x - offset.x, coord.y + offset.y, coord.z, 1.0);
    pB = vec4(coord.x - offset.x, coord.y - offset.y, coord.z, 1.0);
    pC = vec4(coord.x + offset.x, coord.y - offset.y, coord.z, 1.0);
    pD = vec4(coord.x + offset.x, coord.y + offset.y, coord.z, 1.0);
}

void main(){
    vec2 textA, textB, textC, textD;
    vec4 pointA, pointB, pointC, pointD;
    calculate_uv(geom_text_uv[0], textA, textB, textC, textD);
    calculate_points(geom_coord[0], pointA, pointB, pointC, pointD);
    gl_Position = proj_mat * pointA; frag_text_uv = textA; EmitVertex();
    gl_Position = proj_mat * pointB; frag_text_uv = textB; EmitVertex();
    gl_Position = proj_mat * pointD; frag_text_uv = textD; EmitVertex();
    gl_Position = proj_mat * pointC; frag_text_uv = textC; EmitVertex();
    EndPrimitive();
}
"""
f_shader_freetype = """
#version 330

uniform sampler2D textu;
uniform vec4 text_color;

in vec2 frag_text_uv;

out vec4 final_color;

void main(){
    vec4 sampled = vec4(1.0, 1.0, 1.0, texture(textu, frag_text_uv).r);
    if (sampled.a==0.0)
        discard;
    final_color = text_color * sampled;
}
"""
