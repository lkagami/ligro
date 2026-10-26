#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#

vertex_shader_dashed_lines = """
#version 330
precision highp float; 
precision highp int;

layout(std140) uniform CameraMatrices {
    mat4 view_mat;
    mat4 proj_mat;
};
uniform mat4 model_mat;
uniform vec3 uniform_color;
//uniform int test_int;

in vec3 vert_coord;
in vec3 vert_color;

out vec3 geom_color;
out vec4 geom_coord;

const float vert_rad = 0.15;
out float geom_rad;

//const float antialias_length = 0.058;

void main(){
    //geom_width = vert_width;
    //geom_color = vert_color;
    geom_color = mix(uniform_color, vert_color, 0.01);
    geom_coord = view_mat * model_mat * vec4(vert_coord, 1.0);
    geom_rad = vert_rad;
}
"""



#"""
##version 330
#
#uniform mat4 model_mat;
#uniform mat4 view_mat;
#precision highp float; 
#precision highp int;
#
#in vec3 vert_coord;
#in vec3 vert_color;
#
#out vec3 geom_color;
#out vec4 geom_coord;
#
#void main(){
#    
#    geom_color = vert_color;
#    //geom_color = vec3(0.5, 0.5, 0.5);
#    geom_coord = view_mat * model_mat * vec4(vert_coord, 1.0);
#}
#"""
geometry_shader_dashed_lines =  """
#version 330
precision highp float; 
precision highp int;

layout (lines) in;
layout (line_strip, max_vertices = 4) out;

layout(std140) uniform CameraMatrices {
    mat4 view_mat;
    mat4 proj_mat;
};

in vec3 geom_color[];
in vec4 geom_coord[];
in float geom_rad[];

out vec3 frag_color;
out vec4 frag_coord;
out float line_dot_value; 
out float  scalar_distance;
void main(){
    scalar_distance = distance(geom_coord[0].xyz,geom_coord[1].xyz);
    vec4 mid_coord = vec4((geom_coord[0].xyz + geom_coord[1].xyz)/2, 1.0);
    gl_Position = proj_mat * geom_coord[0];
    frag_color = geom_color[0];
    frag_coord = geom_coord[0];
    line_dot_value = 1;
    EmitVertex();
    gl_Position = proj_mat * mid_coord;
    frag_color = geom_color[0];
    frag_coord = mid_coord;
    line_dot_value = 0;
    EmitVertex();
    EndPrimitive();
    gl_Position = proj_mat * mid_coord;
    frag_color = geom_color[1];
    frag_coord = mid_coord;
    line_dot_value = 0;
    EmitVertex();
    gl_Position = proj_mat * geom_coord[1];
    frag_coord = geom_coord[1];
    frag_color = geom_color[1];
    line_dot_value = 1;
    EmitVertex();
    EndPrimitive();
}
"""

#"""
##version 330
#precision highp float; 
#precision highp int;
#
#layout (lines) in;
#layout (line_strip, max_vertices = 4) out;
#
#uniform mat4 proj_mat;
#
#in vec3 geom_color[];
#in vec4 geom_coord[];
#
#out vec3 frag_color;
#out vec4 frag_coord;
#out float line_dot_value;
#
#void main(){
#    vec4 mid_coord = vec4((geom_coord[0].xyz + geom_coord[1].xyz)/2, 1.0);
#    gl_Position = proj_mat * geom_coord[0];
#    frag_color = geom_color[0];
#    frag_coord = geom_coord[0];
#    line_dot_value = 1;
#    EmitVertex();
#    gl_Position = proj_mat * mid_coord;
#    frag_color = geom_color[0];
#    frag_coord = mid_coord;
#    line_dot_value = 0;
#    EmitVertex();
#    EndPrimitive();
#    gl_Position = proj_mat * mid_coord;
#    frag_color = geom_color[1];
#    frag_coord = mid_coord;
#    line_dot_value = 0;
#    EmitVertex();
#    gl_Position = proj_mat * geom_coord[1];
#    frag_coord = geom_coord[1];
#    frag_color = geom_color[1];
#    line_dot_value = 1;
#    EmitVertex();
#    EndPrimitive();
#}
#"""
fragment_shader_dashed_lines = """
#version 330
precision highp float; 
precision highp int;

uniform vec4 fog_color;
uniform float fog_start;
uniform float fog_end;
uniform float fog_intensity;
in float scalar_distance;

in vec3 frag_color;
in vec4 frag_coord;
in float line_dot_value; 

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

    final_color = vec4(frag_color, 1.0);

    if(  mod( round(line_dot_value * 5*scalar_distance) , 2 )  > 0 )
         discard;
    
}
"""

#"""
##version 330
#
#uniform vec4 fog_color;
#uniform float fog_start;
#uniform float fog_end;
#
#
#in vec3 frag_color;
#in vec4 frag_coord;
#in float line_dot_value;
#
#out vec4 final_color;
#
#void main(){
#    //frag_color2 = vec3(0.5,0.5,0.5 );
#    if(mod(round(line_dot_value * 20), 4) > 0)
#         discard;
#    
#    float dist = abs(frag_coord.z);
#    if(dist>=fog_start){
#        float fog_factor = (fog_end-dist)/(fog_end-fog_start);
#        final_color = mix(fog_color, vec4(frag_color, 1.0), fog_factor);
#    }
#    else{
#       final_color = vec4(frag_color, 1.0);
#    }
#}
#"""

sel_vertex_shader_dashed_lines = """
#version 330

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
sel_geometry_shader_dashed_lines = """
#version 330

layout (lines) in;
layout (line_strip, max_vertices = 4) out;

layout(std140) uniform CameraMatrices {
    mat4 view_mat;
    mat4 proj_mat;
};

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
sel_fragment_shader_dashed_lines = """
#version 330

in vec3 frag_color;

out vec4 final_color;

void main(){
    final_color = vec4(frag_color, 1.0);
}
"""
