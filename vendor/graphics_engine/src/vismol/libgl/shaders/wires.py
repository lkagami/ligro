#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#

vertex_shader_wires = """
#version 330
precision highp float; 
precision highp int;

layout(std140) uniform CameraMatrices {
    mat4 view_mat;
    mat4 proj_mat;
};
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
geometry_shader_wires = """
#version 330
precision highp float; 
precision highp int;
layout (triangles) in;
layout (line_strip, max_vertices = 3) out;

layout(std140) uniform CameraMatrices {
    mat4 view_mat;
    mat4 proj_mat;
};

in vec3 geom_color[];
in vec4 geom_coord[];

out vec3 frag_coord;
out vec3 frag_color;

void main(){
    gl_Position = proj_mat * geom_coord[0];
    frag_coord = geom_coord[0].xyz;
    frag_color = geom_color[0];
    EmitVertex();
    gl_Position = proj_mat * geom_coord[1];
    frag_coord = geom_coord[1].xyz;
    frag_color = geom_color[1];
    EmitVertex();
    gl_Position = proj_mat * geom_coord[2];
    frag_coord = geom_coord[2].xyz;
    frag_color = geom_color[2];
    EmitVertex();
    EndPrimitive();
}

"""
fragment_shader_wires = """
#version 330
precision highp float; 
precision highp int;
uniform vec4 fog_color;
uniform float fog_start;
uniform float fog_end;
uniform float fog_intensity;

in vec3 frag_coord;
in vec3 frag_color;

out vec4 final_color;

void main(){
    float dist = abs(frag_coord.z);
    if(dist>=fog_start){
        float fog_factor = 1.0 - fog_intensity * (1.0 - clamp((fog_end-dist)/(fog_end-fog_start), 0.0, 1.0));
        final_color = mix(fog_color, vec4(frag_color, 1.0), fog_factor);
    }
    else{
       final_color = vec4(frag_color, 1.0);
       //final_color = vec4(0.5 , 0.5, 0.5 , 1.0);
    }
}
"""

