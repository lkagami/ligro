#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#

vertex_shader_surface = """
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
in vec3 vert_normal;

out vec3 geom_color;
out vec3 geom_normal;
out vec4 geom_coord;

void main(){
    geom_color = vert_color;
    // transforma a normal pro mesmo espaco (view space) que geom_coord ja
    // usa -- consistente com o norm_vec plano calculado no geometry
    // shader a partir das arestas de geom_coord (tambem em view space).
    // mat3() descarta a translacao; assume-se transformacao rigida
    // (rotacao + translacao, sem escala nao-uniforme), que e o caso de
    // model_mat neste visualizador.
    geom_normal = mat3(view_mat * model_mat) * vert_normal;
    geom_coord = view_mat * model_mat * vec4(vert_coord, 1.0);
}

"""
geometry_shader_surface = """
#version 330

layout(std140) uniform CameraMatrices {
    mat4 view_mat;
    mat4 proj_mat;
};
precision highp float; 
precision highp int;
layout (triangles) in;
layout (triangle_strip, max_vertices = 3) out;


in vec3 geom_color[];
in vec3 geom_normal[];
in vec4 geom_coord[];

uniform int smooth_shading;

out vec3 frag_coord;
out vec3 frag_color;
out vec3 frag_norm;

void main(){
    vec3 vec_p0_p1 = geom_coord[1].xyz - geom_coord[0].xyz;
    vec3 vec_p0_p2 = geom_coord[2].xyz - geom_coord[0].xyz;
    vec3 norm_vec = cross(vec_p0_p1, vec_p0_p2);
    norm_vec = normalize(norm_vec);

    // smooth_shading == 0 (default): flat shading, usa norm_vec (normal
    // de face constante pro triangulo inteiro, calculada acima a partir
    // das proprias arestas) -- mostra as facetas do marching cubes.
    // smooth_shading != 0: usa geom_normal[i], a normal por vertice
    // (media das faces adjacentes, ver compute_smooth_normals em
    // surface_analysis_window.py) interpolada pelo rasterizador --
    // shading suave, sem facetas aparentes.
    vec3 n0 = (smooth_shading != 0) ? normalize(geom_normal[0]) : norm_vec;
    vec3 n1 = (smooth_shading != 0) ? normalize(geom_normal[1]) : norm_vec;
    vec3 n2 = (smooth_shading != 0) ? normalize(geom_normal[2]) : norm_vec;

    gl_Position = proj_mat * geom_coord[0];
    frag_coord  = geom_coord[0].xyz;
    frag_color  = geom_color[0];
    frag_norm   = n0;
    EmitVertex();
    
    gl_Position = proj_mat * geom_coord[1];
    frag_coord  = geom_coord[1].xyz;
    frag_color  = geom_color[1];
    frag_norm   = n1;
    EmitVertex();
    
    gl_Position = proj_mat * geom_coord[2];
    frag_coord  = geom_coord[2].xyz;
    frag_color  = geom_color[2];
    frag_norm   = n2;
    EmitVertex();
    
    EndPrimitive();
}

"""
fragment_shader_surface = """
#version 330
precision highp float; 
precision highp int;
#include <vm_lighting>
uniform float surf_alpha;

uniform vec4 fog_color;
uniform float fog_start;
uniform float fog_end;
uniform float fog_intensity;

in vec3 frag_coord;
in vec3 frag_color;
in vec3 frag_norm;

out vec4 final_color;

void main(){
    vec3 normal = normalize(frag_norm);
    vec4 my_color = vec4(vm_lighting(normal, normalize(frag_coord), frag_color, 1.0), surf_alpha);
    
    float dist = abs(frag_coord.z);
    if(dist>=fog_start){
        float fog_factor = 1.0 - fog_intensity * (1.0 - clamp((fog_end-dist)/(fog_end-fog_start), 0.0, 1.0));
        final_color = mix(fog_color, my_color, fog_factor);
    }
    else{
       final_color = my_color;
    }
    
 
}
"""








vertex_shader_lines = """
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

out vec3 geom_color;
out vec4 geom_coord;

void main(){
    geom_color = vert_color;
    geom_coord = view_mat * model_mat * vec4(vert_coord, 1.0);
}
"""
geometry_shader_lines = """
#version 330

layout(std140) uniform CameraMatrices {
    mat4 view_mat;
    mat4 proj_mat;
};
precision highp float; 
precision highp int;
layout (lines) in;
layout (line_strip, max_vertices = 4) out;


in vec3 geom_color[];
in vec4 geom_coord[];

out vec3 frag_color;
out vec4 frag_coord;

void main(){
    vec4 mid_coord = vec4((geom_coord[0].xyz + geom_coord[1].xyz)/2, 1.0);
    gl_Position = proj_mat * geom_coord[0];
    frag_color = geom_color[0];
    frag_coord = geom_coord[0];
    EmitVertex();
    gl_Position = proj_mat * mid_coord;
    frag_color = geom_color[0];
    frag_coord = mid_coord;
    EmitVertex();
    EndPrimitive();
    gl_Position = proj_mat * mid_coord;
    frag_color = geom_color[1];
    frag_coord = mid_coord;
    EmitVertex();
    gl_Position = proj_mat * geom_coord[1];
    frag_coord = geom_coord[1];
    frag_color = geom_color[1];
    EmitVertex();
    EndPrimitive();
}
"""
fragment_shader_lines = """
#version 330
precision highp float; 
precision highp int;
uniform vec4 fog_color;
uniform float fog_start;
uniform float fog_end;
uniform float fog_intensity;

in vec3 frag_color;
in vec4 frag_coord;

out vec4 final_color;

void main(){
    float dist = abs(frag_coord.z);
    if(dist>=fog_start){
        float fog_factor = 1.0 - fog_intensity * (1.0 - clamp((fog_end-dist)/(fog_end-fog_start), 0.0, 1.0));
        final_color = mix(fog_color, vec4(frag_color, 1.0), fog_factor);
    }
    else{
       final_color = vec4(frag_color, 1.0);
    }
}
"""


################################## SELECTION ###################################

sel_vertex_shader_lines = """
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

out vec3 geom_color;
out vec4 geom_coord;

void main(){
    geom_color = vert_color;
    geom_coord = view_mat * model_mat * vec4(vert_coord, 1.0);
}
"""
sel_geometry_shader_lines = """
#version 330

layout(std140) uniform CameraMatrices {
    mat4 view_mat;
    mat4 proj_mat;
};
precision highp float; 
precision highp int;
layout (lines) in;
layout (line_strip, max_vertices = 4) out;


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






vertex = """
#version 330

layout(std140) uniform CameraMatrices {
    mat4 view_mat;
    mat4 proj_mat;
};
precision highp float; 
precision highp int;
uniform mat4 model_mat;
uniform mat4 normal_mat;
attribute vec3 vert_coord;
attribute vec3 normal;
//vec3 normal = vec3(0.0 , 0.0, 1.0);
varying vec3 v_normal;
varying vec3 v_position;

void main()
{
    gl_Position = <transform>;
    vec4 P = view_mat * model_mat* vec4(vert_coord, 1.0);
    v_position = P.xyz / P.w;
    v_normal = vec3(normal_mat * vec4(normal,0.0));
}
"""

fragment = """
#version 330
precision highp float; 
precision highp int;
varying vec3 v_normal;
varying vec3 v_position;

const vec3 light_position = vec3(1.0,1.0,1.0);
const vec3 ambient_color  = vec3(0.1, 0.0, 0.0);
const vec3 diffuse_color  = vec3(0.75, 0.125, 0.125);
const vec3 specular_color = vec3(1.0, 1.0, 1.0);
const float shininess = 128.0;
const float gamma = 2.2;

void main()
{
    vec3 normal= normalize(v_normal);
    vec3 light_direction = normalize(light_position - v_position);
    float lambertian = max(dot(light_direction,normal), 0.0);
    float specular = 0.0;
    if (lambertian > 0.0)
    {
        vec3 view_direction = normalize(-v_position);
        vec3 half_direction = normalize(light_direction + view_direction);
        float specular_angle = max(dot(half_direction, normal), 0.0);
        specular = pow(specular_angle, shininess);
    }
    vec3 color_linear = ambient_color +
                        lambertian * diffuse_color +
                        specular * specular_color;
    vec3 color_gamma = pow(color_linear, vec3(1.0/gamma));
    gl_FragColor = vec4(color_gamma, 1.0);
}
"""
