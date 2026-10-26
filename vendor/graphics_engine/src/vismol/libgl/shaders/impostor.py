#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#


vertex_shader_glumpy = """
#version 330

layout(std140) uniform CameraMatrices {
    mat4 view_mat;
    mat4 proj_mat;
};

uniform mat4 model_mat;

uniform vec3 u_campos;

in vec3 vert_coord;        // attribute vec3 position;
in vec3 vert_color;        // attribute vec3 color;
in float vert_dot_size;   // attribute float radius;
float hw_ratio;

out vec3 frag_color;       // varying vec3 v_color;
out float f_radius;        // varying float v_radius;
out float f_size;          // varying float v_size;
out vec4 frag_coord;       // varying vec4 v_eye_position;

void main (void){
    hw_ratio = proj_mat[0][0] * proj_mat[1][1];
    frag_color = vert_color;
    f_radius = vert_dot_size;
    frag_coord = view_mat * model_mat * vec4(vert_coord, 1.0);
    gl_Position = proj_mat * frag_coord;
    vec4 p = proj_mat * vec4(vert_dot_size, vert_dot_size, frag_coord.z, frag_coord.w);
    f_size = 256.0 * hw_ratio * vert_dot_size / p.w;
    gl_PointSize = f_size;
}
"""

fragment_shader_glumpy = """
#version 330

layout(std140) uniform CameraMatrices {
    mat4 view_mat;
    mat4 proj_mat;
};

#include <vm_lighting>

uniform mat4 model_mat;

uniform vec4 fog_color;
uniform float fog_start;
uniform float fog_end;
uniform float fog_intensity;

vec4 outline(float distance, float linewidth, float antialias, vec4 fg_color, vec4 bg_color){
    vec4 frag_color;
    float t = linewidth/2.0 - antialias;
    float signed_distance = distance;
    float border_distance = abs(signed_distance) - t;
    float alpha = border_distance/antialias;
    alpha = exp(-alpha*alpha);

    if( border_distance < 0.0 )
        frag_color = fg_color;
    else if( signed_distance < 0.0 )
        frag_color = mix(bg_color, fg_color, sqrt(alpha));
    else {
        if( abs(signed_distance) < (linewidth/2.0 + antialias) ) {
            frag_color = vec4(fg_color.rgb, fg_color.a * alpha);
        } else {
            discard;
        }
    }
    return frag_color;
}

in vec3 frag_color;       // varying vec3 v_color;
in float f_radius;        // varying float v_radius;
in float f_size;          // varying float v_size;
in vec4 frag_coord;       // varying vec4 v_eye_position;

out vec4 final_color;

vec4 calculate_color(vec3 fnrm, vec3 fcrd, vec3 fcol){
    return vec4(vm_lighting(normalize(fnrm), normalize(fcrd), fcol, 0.0), 1.0);
}

void main(){
    vec2 P = gl_PointCoord.xy - vec2(0.5,0.5);
    float distance = length(P*f_size) - f_size/2;
    vec2 texcoord = gl_PointCoord* 2.0 - vec2(1.0);
    float x = texcoord.x;
    float y = texcoord.y;
    float d = 1.0 - x*x - y*y;
    if (d <= 0.0){
        discard;
    }
    float z = sqrt(d);
    vec4 pos = frag_coord;
    pos.z += f_radius*z;
    vec3 pos2 = pos.xyz;
    pos = proj_mat * pos;
    gl_FragDepth = 0.5*(pos.z / pos.w)+0.5;
    vec3 normal = vec3(x,-y,z);
    vec4 color = calculate_color(normal, frag_coord.xyz, frag_color);
    //vec4 temp_color = outline(distance, 1.0, 1.0, vec4(0,0,0,1), color);
    float dist = abs(frag_coord.z);
    if(dist>=fog_start){
        float fog_factor = 1.0 - fog_intensity * (1.0 - clamp((fog_end-dist)/(fog_end-fog_start), 0.0, 1.0));
        final_color = mix(fog_color, color, fog_factor);
    }
    else{
       final_color = color;
    }
}
"""

sel_fragment_shader_glumpy = """
#version 330

layout(std140) uniform CameraMatrices {
    mat4 view_mat;
    mat4 proj_mat;
};

uniform mat4 model_mat;

in vec3 frag_color;       // varying vec3 v_color;
in float f_radius;        // varying float v_radius;
in float f_size;          // varying float v_size;
in vec4 frag_coord;       // varying vec4 v_eye_position;

out vec4 final_color;

void main()
{
    vec2 P = gl_PointCoord.xy - vec2(0.5,0.5);
    float distance = length(P*f_size) - f_size/2;
    vec2 texcoord = gl_PointCoord* 2.0 - vec2(1.0);
    float x = texcoord.x;
    float y = texcoord.y;
    float d = 1.0 - x*x - y*y;
    if (d <= 0.0){
        discard;
    }
    float z = sqrt(d);
    vec4 pos = frag_coord;
    pos.z += f_radius*z;
    pos = proj_mat * pos;
    gl_FragDepth = 0.5*(pos.z / pos.w)+0.5;
    final_color = vec4(frag_color, 1.0);
}
"""


vertex_shader_impostor = """
#version 330

layout(std140) uniform CameraMatrices {
    mat4 view_mat;
    mat4 proj_mat;
};
precision highp float;

uniform mat4 model_mat;

in vec3 vert_coord;
in vec3 vert_color;

//in float vert_dot_size;
//const float vert_dot_size = 0.5;
in float vert_dot_size;   // attribute float radius;

uniform vec3 u_campos;

out vec3 geom_color;
out vec3 geom_coord;
out vec3 geom_center;
out vec3 geom_cam;
out float geom_radius;

void main() {
    geom_color = vert_color;
    geom_coord = (view_mat * model_mat * vec4(vert_coord, 1.0)).xyz;
    geom_center = (view_mat * model_mat * vec4(vert_coord, 1.0)).xyz;
    geom_cam = (view_mat * model_mat * vec4(u_campos, 1.0)).xyz;
    geom_radius = vert_dot_size;
}
"""

geometry_shader_impostor = """
#version 330

layout(std140) uniform CameraMatrices {
    mat4 view_mat;
    mat4 proj_mat;
};

layout (points) in;
layout (triangle_strip, max_vertices = 18) out;

uniform mat4 model_mat;

in vec3 geom_color[];
in vec3 geom_coord[];
in vec3 geom_center[];
in vec3 geom_cam[];
in float geom_radius[];

out vec3 frag_color;
out vec3 frag_coord;
out vec3 frag_center;
out vec3 frag_cam;
out float frag_radius;

vec3 p_1 = vec3(-1.0,-1.0,-1.0);
vec3 p_2 = vec3(-1.0,-1.0, 1.0);
vec3 p_3 = vec3( 1.0,-1.0, 1.0);
vec3 p_4 = vec3( 1.0,-1.0,-1.0);
vec3 p_5 = vec3(-1.0, 1.0,-1.0);
vec3 p_6 = vec3(-1.0, 1.0, 1.0);
vec3 p_7 = vec3( 1.0, 1.0, 1.0);
vec3 p_8 = vec3( 1.0, 1.0,-1.0);

void main(){
    vec3 point1 = geom_coord[0] + p_1 * geom_radius[0];
    vec3 point2 = geom_coord[0] + p_2 * geom_radius[0];
    vec3 point3 = geom_coord[0] + p_3 * geom_radius[0];
    vec3 point4 = geom_coord[0] + p_4 * geom_radius[0];
    vec3 point5 = geom_coord[0] + p_5 * geom_radius[0];
    vec3 point6 = geom_coord[0] + p_6 * geom_radius[0];
    vec3 point7 = geom_coord[0] + p_7 * geom_radius[0];
    vec3 point8 = geom_coord[0] + p_8 * geom_radius[0];
    
    gl_Position = proj_mat * vec4(point1, 1.0);
    frag_color = geom_color[0];
    frag_coord = point1;
    frag_center = geom_center[0];
    frag_cam = geom_cam[0];
    frag_radius = geom_radius[0];
    EmitVertex();
    gl_Position = proj_mat * vec4(point6, 1.0);
    frag_color = geom_color[0];
    frag_coord = point6;
    frag_center = geom_center[0];
    frag_cam = geom_cam[0];
    frag_radius = geom_radius[0];
    EmitVertex();
    gl_Position = proj_mat * vec4(point2, 1.0);
    frag_color = geom_color[0];
    frag_coord = point2;
    frag_center = geom_center[0];
    frag_cam = geom_cam[0];
    frag_radius = geom_radius[0];
    EmitVertex();
    gl_Position = proj_mat * vec4(point7, 1.0);
    frag_color = geom_color[0];
    frag_coord = point7;
    frag_center = geom_center[0];
    frag_cam = geom_cam[0];
    frag_radius = geom_radius[0];
    EmitVertex();
    gl_Position = proj_mat * vec4(point3, 1.0);
    frag_color = geom_color[0];
    frag_coord = point3;
    frag_center = geom_center[0];
    frag_cam = geom_cam[0];
    frag_radius = geom_radius[0];
    EmitVertex();
    gl_Position = proj_mat * vec4(point8, 1.0);
    frag_color = geom_color[0];
    frag_coord = point8;
    frag_center = geom_center[0];
    frag_cam = geom_cam[0];
    frag_radius = geom_radius[0];
    EmitVertex();
    gl_Position = proj_mat * vec4(point4, 1.0);
    frag_color = geom_color[0];
    frag_coord = point4;
    frag_center = geom_center[0];
    frag_cam = geom_cam[0];
    frag_radius = geom_radius[0];
    EmitVertex();
    gl_Position = proj_mat * vec4(point5, 1.0);
    frag_color = geom_color[0];
    frag_coord = point5;
    frag_center = geom_center[0];
    frag_cam = geom_cam[0];
    frag_radius = geom_radius[0];
    EmitVertex();
    gl_Position = proj_mat * vec4(point1, 1.0);
    frag_color = geom_color[0];
    frag_coord = point1;
    frag_center = geom_center[0];
    frag_cam = geom_cam[0];
    frag_radius = geom_radius[0];
    EmitVertex();
    gl_Position = proj_mat * vec4(point6, 1.0);
    frag_color = geom_color[0];
    frag_coord = point6;
    frag_center = geom_center[0];
    frag_cam = geom_cam[0];
    frag_radius = geom_radius[0];
    EmitVertex();
    
    EndPrimitive();
    gl_Position = proj_mat * vec4(point1, 1.0);
    frag_color = geom_color[0];
    frag_coord = point1;
    frag_center = geom_center[0];
    frag_cam = geom_cam[0];
    frag_radius = geom_radius[0];
    EmitVertex();
    gl_Position = proj_mat * vec4(point2, 1.0);
    frag_color = geom_color[0];
    frag_coord = point2;
    frag_center = geom_center[0];
    frag_cam = geom_cam[0];
    frag_radius = geom_radius[0];
    EmitVertex();
    gl_Position = proj_mat * vec4(point4, 1.0);
    frag_color = geom_color[0];
    frag_coord = point4;
    frag_center = geom_center[0];
    frag_cam = geom_cam[0];
    frag_radius = geom_radius[0];
    EmitVertex();
    gl_Position = proj_mat * vec4(point3, 1.0);
    frag_color = geom_color[0];
    frag_coord = point3;
    frag_center = geom_center[0];
    frag_cam = geom_cam[0];
    frag_radius = geom_radius[0];
    EmitVertex();
    EndPrimitive();

    gl_Position = proj_mat * vec4(point5, 1.0);
    frag_color = geom_color[0];
    frag_coord = point5;
    frag_center = geom_center[0];
    frag_cam = geom_cam[0];
    frag_radius = geom_radius[0];
    EmitVertex();
    gl_Position = proj_mat * vec4(point6, 1.0);
    frag_color = geom_color[0];
    frag_coord = point6;
    frag_center = geom_center[0];
    frag_cam = geom_cam[0];
    frag_radius = geom_radius[0];
    EmitVertex();
    gl_Position = proj_mat * vec4(point8, 1.0);
    frag_color = geom_color[0];
    frag_coord = point8;
    frag_center = geom_center[0];
    frag_cam = geom_cam[0];
    frag_radius = geom_radius[0];
    EmitVertex();
    gl_Position = proj_mat * vec4(point7, 1.0);
    frag_color = geom_color[0];
    frag_coord = point7;
    frag_center = geom_center[0];
    frag_cam = geom_cam[0];
    frag_radius = geom_radius[0];
    EmitVertex();
    EndPrimitive();
}
"""

fragment_shader_impostor = """
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
in vec3 frag_center;
in vec3 frag_cam;
in float frag_radius;

out vec4 final_color;

float sph_intersect(vec3 ro, vec3 rd, vec3 sph, float rad){
    vec3 oc = ro - sph;
    float b = dot(oc, rd);
    float c = dot(oc, oc) - rad*rad;
    float h = b*b - c;
    if( h<0.0 ) return -1.0;
    return -b - sqrt(h);
}

void main() {
    vec3 ray_orig = frag_cam;
    vec3 ray_dir = normalize(frag_coord - frag_cam);
    float dist_to_sph = sph_intersect(ray_orig, ray_dir, frag_center, frag_radius);
    if (dist_to_sph < 0.0) discard;
    vec3 coord_on_sph = ray_orig + ray_dir * dist_to_sph;
    vec3 normal = normalize(coord_on_sph - frag_center);
    vec4 lit_color = vec4(vm_lighting(normal, ray_dir, frag_color, 0.0), 1.0);
    //final_color = vec4(ambient + diffuse, 1.0);
    
    //vec4 depth = proj_mat * vec4(frag_center + normal * frag_radius, 1.0);
    //gl_FragDepthEXT = depth.z/depth.w;
    vec4 depth = proj_mat * vec4(coord_on_sph, 1.0);
    //gl_FragDepth = depth.z/depth.w;
    
    float dist = abs(depth.z);
    if(dist>=fog_start){
        float fog_factor = 1.0 - fog_intensity * (1.0 - clamp((fog_end-dist)/(fog_end-fog_start), 0.0, 1.0));
        final_color = mix(fog_color, lit_color, fog_factor);
    }
    else{
       final_color = lit_color;
    }
}
"""

sel_fragment_shader_impostor = """
#version 330
#extension GL_EXT_frag_depth: enable
precision highp float;

in vec3 frag_color;
in vec3 frag_coord;
in vec3 frag_center;
in vec3 frag_cam;
in float frag_radius;

out vec4 final_color;

float sph_intersect(vec3 ro, vec3 rd, vec3 sph, float rad){
    vec3 oc = ro - sph;
    float b = dot(oc, rd);
    float c = dot(oc, oc) - rad*rad;
    float h = b*b - c;
    if( h<0.0 ) return -1.0;
    return -b - sqrt(h);
}

void main() {
    vec3 ray_orig = frag_cam;
    vec3 ray_dir = normalize(frag_coord - frag_cam);
    float dist_to_sph = sph_intersect(ray_orig, ray_dir, frag_center, frag_radius);
    if (dist_to_sph < 0.0) discard;
    final_color = vec4(frag_color, 1.0);
}
"""




shader_type ={0: { "vertex_shader"       : vertex_shader_glumpy,
                   "fragment_shader"     : fragment_shader_glumpy,
                   "geometry_shader"     : None,
                   "sel_vertex_shader"   : vertex_shader_glumpy,
                   "sel_fragment_shader" : sel_fragment_shader_glumpy,
                   "sel_geometry_shader" : None,
                 },
              
              1: {"vertex_shader"       : vertex_shader_impostor,
                  "fragment_shader"     : fragment_shader_impostor,
                  "geometry_shader"     : geometry_shader_impostor,
                  "sel_vertex_shader"   : vertex_shader_impostor,
                  "sel_fragment_shader" : sel_fragment_shader_impostor,
                  "sel_geometry_shader" : geometry_shader_impostor,
                  }
}
