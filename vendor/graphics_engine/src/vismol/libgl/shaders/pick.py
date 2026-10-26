#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
#  vismol_shaders.py
#  
#  Copyright 2016 Carlos Eduardo Sequeiros Borja <casebor@gmail.com>
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


#
#vertex_shader_picked = """
##version 330
#
#uniform mat4 model_mat;
#uniform mat4 view_mat;
#uniform mat4 proj_mat;
#
#in vec3 vert_coord;
#in vec3 vert_centr;
#in vec3 vert_color;
#
#varying vec3 vert_norm;
#
#out vec3 frag_coord;
#out vec3 frag_color;
#out vec3 frag_norm;
#
#void main(){
#    mat4 modelview = view_mat * model_mat;
#    gl_Position = proj_mat * modelview * vec4(vert_coord, 1.0);
#    vert_norm = normalize(vert_coord - vert_centr);
#    frag_coord = -vec3(modelview * vec4(vert_coord, 1.0));
#    frag_norm = mat3(transpose(inverse(model_mat))) * vert_norm;
#    frag_color = vert_color;
#}
#"""
#fragment_shader_picked = """
##version 330
#
#struct Light {
#    vec3 position;
#    //vec3 color;
#    vec3 intensity;
#    //vec3 specular_color;
#    float ambient_coef;
#    float shininess;
#};
#
#uniform Light my_light;
#
#uniform vec4 fog_color;
#uniform float fog_start;
#uniform float fog_end;
#
#const float alpha = 0.5;
#
#in vec3 frag_coord;
#in vec3 frag_color;
#in vec3 frag_norm;
#
#out vec4 final_color;
#
#void main(){
#    vec3 normal = normalize(frag_norm);
#    vec3 vert_to_light = normalize(my_light.position);
#    vec3 vert_to_cam = normalize(frag_coord);
#    
#    // Ambient Component
#    vec3 ambient = my_light.ambient_coef * frag_color * my_light.intensity;
#    
#    // Diffuse component
#    float diffuse_coef = max(0.0, dot(normal, vert_to_light));
#    vec3 diffuse = diffuse_coef * frag_color * my_light.intensity;
#    
#    // Specular component
#    float specular_coef = 0.0;
#    if (diffuse_coef > 0.0)
#        specular_coef = pow(max(0.0, dot(vert_to_cam, reflect(-vert_to_light, normal))), my_light.shininess);
#    vec3 specular = specular_coef * my_light.intensity;
#    specular = specular * (vec3(1) - diffuse);
#    
#    vec4 my_color = vec4(ambient + diffuse + specular, alpha);
#    
#    float dist = abs(frag_coord.z);
#    if(dist>=fog_start){
#        float fog_factor = (fog_end-dist)/(fog_end-fog_start);
#        final_color = mix(fog_color, my_color, fog_factor);
#    }
#    else{
#        final_color = my_color;
#    }
#}
#"""
#
#


'''
In this modified vertex shader:
    We calculate the vertices of a square around the input vert_position. The square is centered on vert_position, and we use a small offset (0.1) to define the size of the square.
    We check whether the current vert_position is inside the square using the insideSquare boolean variable.
    The color for each vertex is set based on whether it's inside the square. If it's inside, we use the vert_color input attribute; otherwise, we set the color to black (vec3(0.0)).
    The vertex position is set to vert_position.
With this modification, the vertex shader creates a square around the input position, and the color of each vertex is determined by whether it's inside the square or not. The fragment shader can remain the same as the one you originally provided because it simply uses the color set in the vertex shader.
You can use this modified vertex shader along with your existing fragment shader to draw squares based on the input positions and colors.
'''


#vertex_shader_picking_dots_safe = '''
#version 330
#precision highp float;
#precision highp int;
#
#uniform mat4 model_mat;
#uniform mat4 view_mat;
#uniform mat4 proj_mat;
#
#
#in vec3 vert_coord;
#in vec3 vert_color;
#
#out vec3 frag_color;
#
#void main() {
#    // Define the size of the square
#    float squareSize = 0.1;
#
#    // Calculate the position of the square's vertices
#    vec3 bottomLeft = vert_coord - vec3(squareSize, squareSize, 0.0);
#    vec3 bottomRight = vert_coord + vec3(squareSize, -squareSize, 0.0);
#    vec3 topLeft = vert_coord + vec3(-squareSize, squareSize, 0.0);
#    vec3 topRight = vert_coord + vec3(squareSize, squareSize, 0.0);
#
#    // Set the vertex position and color based on the square's vertices
#    if (gl_VertexID == 0) {
#        gl_Position = vec4(bottomLeft, 1.0);
#        frag_color = vert_color;
#    }
#    else if (gl_VertexID == 1) {
#        gl_Position = vec4(bottomRight, 1.0);
#        frag_color = vert_color;
#    }
#    else if (gl_VertexID == 2) {
#        gl_Position = vec4(topLeft, 1.0);
#        frag_color = vert_color;
#    }
#    else if (gl_VertexID == 3) {
#        gl_Position = vec4(topRight, 1.0);
#        frag_color = vert_color;
#    }
#}
#'''


'''These are the correct shaders.'''
vertex_shader_picking_dots_safe = """
#version 330
precision highp float; 
precision highp int;
uniform mat4 model_mat;
uniform mat4 view_mat;
uniform mat4 proj_mat;


uniform int _size; // Integer uniform

in vec3  vert_coord;
in vec3  vert_color;
out vec3 index_color;

void main(){
    gl_Position  = proj_mat * view_mat * model_mat * vec4(vert_coord, 1.0);
    gl_PointSize = _size;
    index_color = vert_color;
}
"""
#fragment_shader_picking_dots_safe ='''
##version 330 core
#precision highp float;
#precision highp int;
#
#in vec3 index_color;
#out vec4 frag_color;
#
#void main() {
#    // Coordenadas do ponto normalizadas de (0,0) a (1,1)
#    vec2 coord = gl_PointCoord.xy;
#    float dist = length(coord - vec2(0.5));
#
#    // Raio do círculo interno (cor do ponto)
#    float inner_radius = 0.5;
#    // Largura do contorno (ajustável)
#    float border_thickness = 0.08;
#    float outer_radius = inner_radius + border_thickness;
#
#    if (dist > outer_radius) {
#        discard;
#    } else if (dist > inner_radius) {
#        // Região do contorno
#        frag_color = vec4(0.0, 0.0, 0.0, 1.0); // preto
#    } else {
#        // Interior do ponto
#        frag_color = vec4(index_color, 1.0);
#    }
#}
#'''

fragment_shader_picking_dots_safe = """
#version 330
precision highp float;
precision highp int;
in vec3 index_color;
out vec4 frag_color;

void main(){
    float dist = length(gl_PointCoord.xy - vec2(0.5,0.5));
    if (dist > 0.6)
        discard;
    frag_color = vec4(index_color,1.0);
}

"""






'''These shaders are obsolete. We should avoid using the glPointSize 
function (this function is very old and should not be used in modern openGL)'''

vertex_shader_picking_dots = """
#version 330
precision highp float; 
precision highp int;
uniform mat4 model_mat;
uniform mat4 view_mat;
uniform mat4 proj_mat;
  


in vec3  vert_coord;
in vec3  vert_color;
out vec3 index_color;

void main(){
    gl_Position  = proj_mat * view_mat * model_mat * vec4(vert_coord, 1.0);
    //gl_PointSize = 15;
    index_color = vert_color;
}
"""

fragment_shader_picking_dots = """
#version 330
precision highp float;
precision highp int;

in vec3 index_color;
out vec4 frag_color;

void main(){
    float dist = length(gl_PointCoord.xy - vec2(0.5,0.5));
    if (dist > 0.6)
        discard;
    frag_color = vec4(index_color,1.0);
}

"""

############################### VisMolDrawWidget ###############################













