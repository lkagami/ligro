#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#


vertex_shader_spheres = """
#version 330

layout(std140) uniform CameraMatrices {
    mat4 view_mat;
    mat4 proj_mat;
};

uniform mat4 model_mat;

in vec3 vert_coord;
in vec3 vert_color;
in vec3 vert_instance;
in float vert_radius;

vec3 vert_norm;

out vec3 frag_coord;
out vec3 frag_color;
out vec3 frag_norm;

void main(){
    mat4 modelview = view_mat * model_mat;
    vec3 offset_coord = vert_coord * vert_radius + vert_instance;
    gl_Position = proj_mat * modelview * vec4(offset_coord, 1.0);
    
    vert_norm = normalize(offset_coord - vert_instance);
    frag_coord = vec3(modelview * vec4(offset_coord, 1.0));
    frag_norm = mat3(transpose(inverse(model_mat))) * vert_norm;
    frag_color = vert_color;
}
"""
fragment_shader_spheres = """
#version 330

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
        final_color = mix(fog_color, out_color, fog_factor );
    }
    else{
       final_color = out_color;
    }
}
"""


sel_vertex_shader_spheres = """
#version 330

layout(std140) uniform CameraMatrices {
    mat4 view_mat;
    mat4 proj_mat;
};

uniform mat4 model_mat;

in vec3 vert_coord;
in vec3 vert_color;
in vec3 vert_instance;
in float vert_radius;

out vec3 frag_color;

void main(){
    frag_color = vert_color;
    vec3 offset_coord = vert_coord * vert_radius + vert_instance;
    gl_Position = proj_mat * view_mat * model_mat * vec4(offset_coord, 1.0);
}
"""

sel_fragment_shader_spheres = """
#version 330

in vec3 frag_color;

out vec4 final_color;

void main(){
    final_color = vec4(frag_color, 1.0);
}
"""



#picking_spheres
vertex_shader_picking_spheres = """
#version 330

uniform mat4 model_mat;
uniform mat4 view_mat;
uniform mat4 proj_mat;

in vec3 vert_coord;
in vec3 vert_color;
in vec3 vert_instance;
in float vert_radius;

out vec3 frag_color;

void main(){
    frag_color = vert_color;
    vec3 offset_coord = vert_coord * vert_radius + vert_instance;
    gl_Position = proj_mat * view_mat * model_mat * vec4(offset_coord, 1.0);
}
"""

fragment_shader_picking_spheres = """
#version 330

in vec3 frag_color;

out vec4 final_color;

void main(){
    final_color = vec4(frag_color, 0.5);
}
"""




