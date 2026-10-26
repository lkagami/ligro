#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#

vertex_shader_nlines = """
#version 330
precision highp float; 
precision highp int;

uniform mat4 model_mat;
uniform mat4 view_mat;

in vec3 vert_coord;
in vec3 vert_color;
const float vert_width = 0.001;
const float antialias_length = 0.058;

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
    }
}
"""


sel_vertex_shader_nlines = """
#version 330
precision highp float; 
precision highp int;

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
precision highp float; 
precision highp int;

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
precision highp float; 
precision highp int;

in vec3 frag_color;
out vec4 final_color;

void main(){
    final_color = vec4(frag_color, 1.0);
}
"""
