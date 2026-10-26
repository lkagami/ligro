#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
#  dots.py
#  
#  Copyright 2022 Carlos Eduardo Sequeiros Borja <casebor@gmail.com>
#  
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#  
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#  
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#  MA 02110-1301, USA.
#  
#  







vertex_shader_dot_simple = """
# version 330

layout(std140) uniform CameraMatrices {
    mat4 view_mat;
    mat4 proj_mat;
};

uniform mat4 model_mat;

in vec3 vert_coord;
in vec3 vert_color;
out vec3 frag_coord;
out vec3 frag_color;

void main(){
    gl_Position = proj_mat * view_mat * model_mat * vec4(vert_coord, 1.0);
    frag_coord = (view_mat * model_mat * vec4(vert_coord, 1.0)).xyz;
    frag_color = vert_color;
}
"""

fragment_shader_dot_simple = """
# version 330

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





sel_fragment_shader_dot_simple = """
# version 330

in vec3 frag_color;
out vec4 final_color;

void main(){
    final_color = vec4(frag_color, 1.0);
}
"""


vertex_shader_dot_circle = """
# version 330

layout(std140) uniform CameraMatrices {
    mat4 view_mat;
    mat4 proj_mat;
};

uniform mat4 model_mat;
uniform float point_size;

in vec3 vert_coord;
in vec3 vert_color;
out vec3 frag_coord;
out vec3 frag_color;

void main(){
    gl_Position = proj_mat * view_mat * model_mat * vec4(vert_coord, 1.0);
    // NOTA: com GL_VERTEX_PROGRAM_POINT_SIZE habilitado (necessario para
    // usar gl_PointCoord no fragment shader), o glPointSize() do lado
    // Python e IGNORADO pela especificacao OpenGL - o tamanho do ponto
    // so pode vir daqui. Sem essa linha, o tamanho fica indefinido pelo
    // driver: alguns (NVIDIA) mostram algo visivel por sorte, outros
    // (Mesa/Intel) renderizam tamanho efetivamente zero (invisivel, sem
    // erro nenhum). Ver representations.py: point_size e enviado via
    // glUniform1f a cada draw, com o mesmo calculo que antes ia so pro
    // glPointSize.
    gl_PointSize = point_size;
    frag_coord = (view_mat * model_mat * vec4(vert_coord, 1.0)).xyz;
    frag_color = vert_color;
}
"""

fragment_shader_dot_circle = """
# version 330

uniform vec4 fog_color;
uniform float fog_start;
uniform float fog_end;
uniform float fog_intensity;

in vec3 frag_coord;
in vec3 frag_color;
out vec4 final_color;

void main(){
    float dist = length(gl_PointCoord.xy - vec2(0.5,0.5));
    if (dist > 0.5)
        discard;
    
    dist = abs(frag_coord.z);
    if(dist>=fog_start){
        float fog_factor = 1.0 - fog_intensity * (1.0 - clamp((fog_end-dist)/(fog_end-fog_start), 0.0, 1.0));
        final_color = mix(fog_color, vec4(frag_color, 1.0), fog_factor);
    }
    else{
       final_color = vec4(frag_color, 1.0);
    }
}
"""

sel_fragment_shader_dot_circle = """
# version 330

in vec3 frag_color;
out vec4 final_color;

void main(){
    float dist = length(gl_PointCoord.xy - vec2(0.5,0.5));
    if (dist > 0.5)
        discard;
    final_color = vec4(frag_color, 1.0);
}
"""


vertex_shader_dot_disc = """
# version 330

layout(std140) uniform CameraMatrices {
    mat4 view_mat;
    mat4 proj_mat;
};

uniform mat4 model_mat;

in vec3 vert_coord;
in vec3 vert_color;
out vec3 frag_coord;
out vec3 frag_color;

void main(){
    gl_Position = proj_mat * view_mat * model_mat * vec4(vert_coord, 1.0);
    frag_coord = (view_mat * model_mat * vec4(vert_coord, 1.0)).xyz;
    frag_color = vert_color;
}
"""

fragment_shader_dot_disc = """
# version 330

uniform vec4 fog_color;
uniform float fog_start;
uniform float fog_end;
uniform float fog_intensity;

in vec3 frag_coord;
in vec3 frag_color;
out vec4 final_color;

void main(){
    vec4 circle_color;
    float dist = length(gl_PointCoord.xy - vec2(0.5,0.5));
    if (dist > 0.5)
        discard;
    float ligth_dist = length(gl_PointCoord - vec2(0.3, 0.3));
    circle_color = mix(vec4(frag_color, 1), vec4(0, 0, 0, 1), sqrt(ligth_dist)*.78);
    
    dist = abs(frag_coord.z);
    if(dist>=fog_start){
        float fog_factor = 1.0 - fog_intensity * (1.0 - clamp((fog_end-dist)/(fog_end-fog_start), 0.0, 1.0));
        final_color = mix(fog_color, circle_color, fog_factor);
    }
    else{
       final_color = circle_color;
    }
}
"""

sel_fragment_shader_dot_disc = """
# version 330

in vec3 frag_color;
out vec4 final_color;

void main(){
    float dist = length(gl_PointCoord.xy - vec2(0.5,0.5));
    if (dist > 0.5)
        discard;
    final_color = vec4(frag_color, 1.0);
}
"""


vertex_shader_dot_extra = """
# version 330

layout(std140) uniform CameraMatrices {
    mat4 view_mat;
    mat4 proj_mat;
};

uniform mat4 model_mat;

in vec3 vert_coord;
in vec3 vert_color;
out vec4 frag_coord;
out vec3 frag_color;

void main(){
    gl_Position = proj_mat * view_mat * model_mat * vec4(vert_coord, 1.0);
    frag_coord = view_mat * model_mat * vec4(vert_coord, 1.0);
    frag_color = vert_color;
}
"""

fragment_shader_dot_extra = """
# version 330

layout(std140) uniform CameraMatrices {
    mat4 view_mat;
    mat4 proj_mat;
};

in vec4 frag_coord;
in vec3 frag_color;
out vec4 final_color;

void main(){
    vec2 P = gl_PointCoord.xy - vec2(0.5,0.5);
    float point_size = 5.0;
    float distance = length(P);
    vec2 texcoord = gl_PointCoord* 2.0 - vec2(1.0);
    float x = texcoord.x;
    float y = texcoord.y;
    float d = 1.0 - x*x - y*y;
    if (d <= 0.0)
        discard;
    float z = sqrt(d);
    vec4 pos = frag_coord;
    pos.z += z;
    pos = proj_mat * pos;
    gl_FragDepth = 0.5*(pos.z / pos.w)+0.5;
    float dist = length(gl_PointCoord.xy - vec2(0.5,0.5));
    if (dist > 0.5)
        discard;
    float ligth_dist = length(gl_PointCoord - vec2(0.3, 0.3));
    final_color = mix(vec4(frag_color, 1), vec4(0, 0, 0, 1), sqrt(ligth_dist)*.78);
}
"""

sel_fragment_shader_dot_extra = """
# version 330

in vec3 frag_color;
out vec4 final_color;

void main(){
    float dist = length(gl_PointCoord.xy - vec2(0.5,0.5));
    if (dist > 0.5)
        discard;
    final_color = vec4(frag_color, 1.0);
}
"""








vertex_shader_posdot_simple = """
# version 330

layout(std140) uniform CameraMatrices {
    mat4 view_mat;
    mat4 proj_mat;
};

uniform mat4 model_mat;

in vec3 vert_coord;
in vec3 vert_color;
out vec3 frag_coord;
out vec3 frag_color;

void main(){
    gl_Position = proj_mat * view_mat * model_mat * vec4(vert_coord, 1.0);
    frag_coord = (view_mat * model_mat * vec4(vert_coord, 1.0)).xyz;
    frag_color = vert_color;
    //frag_color = vec3(1.0, 0.0, 1.0);
}
"""

fragment_shader_posdot_simple = """
#version 330

uniform vec4 fog_color;
uniform float fog_start;
uniform float fog_end;
uniform float fog_intensity;

in vec3 frag_coord;
in vec3 frag_color;
out vec4 final_color;

void main(){
    final_color = vec4(1.0, 0.0, 1.0, 1.0); // roxo/magenta
}
"""









shader_type = {0: { "vertex_shader"      : vertex_shader_dot_simple,
                    "fragment_shader"    : fragment_shader_dot_simple,
                    "sel_vertex_shader"  : vertex_shader_dot_simple,
                    "sel_fragment_shader": sel_fragment_shader_dot_simple
                  },
               1: {"vertex_shader"      : vertex_shader_dot_circle,
                   "fragment_shader"    : fragment_shader_dot_circle ,
                   "sel_vertex_shader"  : vertex_shader_dot_circle,
                   "sel_fragment_shader": sel_fragment_shader_dot_circle
                   },
               2: {"vertex_shader"      : vertex_shader_dot_disc,
                   "fragment_shader"    : fragment_shader_dot_disc,
                   "sel_vertex_shader"  : vertex_shader_dot_disc,
                   "sel_fragment_shader": sel_fragment_shader_dot_disc
                   },
               3: {"vertex_shader"      : vertex_shader_dot_extra,
                   "fragment_shader"    : fragment_shader_dot_extra,
                   "sel_vertex_shader"  : vertex_shader_dot_extra,
                   "sel_fragment_shader": sel_fragment_shader_dot_extra
                   },
               4: {"vertex_shader"      : vertex_shader_posdot_simple,
                   "fragment_shader"    : fragment_shader_posdot_simple,
                   "sel_vertex_shader"  : vertex_shader_posdot_simple,
                   "sel_fragment_shader": fragment_shader_posdot_simple
                   },
}
