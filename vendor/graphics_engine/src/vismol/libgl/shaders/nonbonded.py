#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#

vertex_shader_non_bonded = """
#version 330

in vec3 vert_coord;
in vec3 vert_color;

out vec3 geom_color;
out vec4 geom_coord;

void main(){
    geom_color = vert_color;
    geom_coord = vec4(vert_coord, 1.0);
}
"""
geometry_shader_non_bonded = """
#version 330

const float xyz_offset = 0.3;

layout (points) in;
layout (line_strip, max_vertices = 6) out;

layout(std140) uniform CameraMatrices {
    mat4 view_mat;
    mat4 proj_mat;
};
uniform mat4 model_mat;


in vec3 geom_color[];
in vec4 geom_coord[];

out vec3 frag_color;
out vec4 frag_coord;

void main(){
    vec4 p1 = geom_coord[0];
    vec4 p2 = geom_coord[0];
    vec4 p3 = geom_coord[0];
    vec4 p4 = geom_coord[0];
    vec4 p5 = geom_coord[0];
    vec4 p6 = geom_coord[0];
    p1.x -= xyz_offset;
    p2.x += xyz_offset;
    p3.y -= xyz_offset;
    p4.y += xyz_offset;
    p5.z -= xyz_offset;
    p6.z += xyz_offset;
    
    frag_coord = view_mat * model_mat * p1;
    frag_color = geom_color[0];
    gl_Position = proj_mat * frag_coord;
    EmitVertex();
    frag_coord = view_mat * model_mat * p2;
    frag_color = geom_color[0];
    gl_Position = proj_mat * frag_coord;
    EmitVertex();
    EndPrimitive();
    
    frag_coord = view_mat * model_mat * p3;
    frag_color = geom_color[0];
    gl_Position = proj_mat * frag_coord;
    EmitVertex();
    frag_coord = view_mat * model_mat * p4;
    frag_color = geom_color[0];
    gl_Position = proj_mat * frag_coord;
    EmitVertex();
    EndPrimitive();
    
    frag_coord = view_mat * model_mat * p5;
    frag_color = geom_color[0];
    gl_Position = proj_mat * frag_coord;
    EmitVertex();
    frag_coord = view_mat * model_mat * p6;
    frag_color = geom_color[0];
    gl_Position = proj_mat * frag_coord;
    EmitVertex();
    EndPrimitive();
}
"""
fragment_shader_non_bonded = """
#version 330

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

sel_vertex_shader_non_bonded = """
#version 330

in vec3 vert_coord;
in vec3 vert_color;

out vec3 geom_color;
out vec4 geom_coord;

void main(){
    geom_color = vert_color;
    geom_coord = vec4(vert_coord, 1.0);
}
"""
sel_geometry_shader_non_bonded = """
#version 330

const float xyz_offset = 0.3;

layout (points) in;
layout (line_strip, max_vertices = 6) out;

layout(std140) uniform CameraMatrices {
    mat4 view_mat;
    mat4 proj_mat;
};
uniform mat4 model_mat;


in vec3 geom_color[];
in vec4 geom_coord[];

out vec3 frag_color;

void main(){
    vec4 p1 = geom_coord[0];
    vec4 p2 = geom_coord[0];
    vec4 p3 = geom_coord[0];
    vec4 p4 = geom_coord[0];
    vec4 p5 = geom_coord[0];
    vec4 p6 = geom_coord[0];
    p1.x -= xyz_offset;
    p2.x += xyz_offset;
    p3.y -= xyz_offset;
    p4.y += xyz_offset;
    p5.z -= xyz_offset;
    p6.z += xyz_offset;
    
    frag_color = geom_color[0];
    gl_Position = proj_mat * view_mat * model_mat * p1;
    EmitVertex();
    frag_color = geom_color[0];
    gl_Position = proj_mat * view_mat * model_mat * p2;
    EmitVertex();
    EndPrimitive();
    
    frag_color = geom_color[0];
    gl_Position = proj_mat * view_mat * model_mat * p3;
    EmitVertex();
    frag_color = geom_color[0];
    gl_Position = proj_mat * view_mat * model_mat * p4;
    EmitVertex();
    EndPrimitive();
    
    frag_color = geom_color[0];
    gl_Position = proj_mat * view_mat * model_mat * p5;
    EmitVertex();
    frag_color = geom_color[0];
    gl_Position = proj_mat * view_mat * model_mat * p6;
    EmitVertex();
    EndPrimitive();
}
"""
sel_fragment_shader_non_bonded = """
#version 330

in vec3 frag_color;

out vec4 final_color;

void main(){
    final_color = vec4(frag_color, 1.0);
}
"""
