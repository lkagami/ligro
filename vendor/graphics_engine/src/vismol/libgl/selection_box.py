#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
#  selection_box.py
#  
#  Copyright 2017 Carlos Eduardo Sequeiros Borja <casebor@gmail.com>
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

from vismol.utils.debug import dprint
import ctypes
import numpy as np
from OpenGL import GL

class SelectionBox:
    """ This class represent the selection box when you press shift in VisMol
        and drag the mouse. This class is meant to be as independent as possible
        from the VisMol widget, just like the GLAxis class.
    """
    
    def __init__ (self):
        """ Constructor that defines the color of the box, the indices to form
            the box and the triangles for the fill box. The color can be changed
            here, but you should not change the other parameters, since that
            will cause unwanted behavior.
            Despite the box only needs 4 points to create a quad, we use 5 to
            obtain a closed quad, this is because the GL_LINE_STRIP function do
            not join the first item with the last.
            For the box lines we use the indices disposed as:
            
                       0,4----------3
                        |           |
                        1-----------2
            
            For the triangles we use the same indices but the order will be:
            
                  0,4----------3        0       3
                   |           |        | \     |
                   |           |        |  \    |
                   |           |   ==>  |   \   |
                   |           |        |    \  |
                   |           |        |     \ |
                   1-----------2        1       2
            
        """
        self.vao = None
        self.buffers = None
        self.selection_box_program = None
        self.start = None
        self.end = None
        self.points = None
        self.color = np.array([0.0, 0.5, 0.5]*5, dtype=np.float32)
        self.indices = np.array([0, 1, 2, 3, 4], dtype=np.uint32)
        self.triangles = np.array([1, 0, 2, 3], dtype=np.uint32)
        self.vertex_shader_sb = """
#version 330

in vec2 vert_coord;
in vec3 vert_color;

out vec3 frag_color;

void main(){
    gl_Position = vec4(vert_coord, 0.0, 1.0);
    frag_color = vert_color;
}
"""
        self.fragment_shader_sb = """
#version 330

in vec3 frag_color;

out vec4 final_color;

void main(){
    final_color = vec4(frag_color, 0.6);
}
"""
    
    def initialize_gl(self):
        """ Gives the order to create and compile the OpenGL program. Maybe this
            function is redundant, only putted this way to maintain some format
            with GLAxis class.
        """
        self._make_selection_box_program()
        return True
    
    def _make_selection_box_program(self):
        """ Creates, compiles and attach the shaders to the selection box
            program.
        """
        v_shader = GL.glCreateShader(GL.GL_VERTEX_SHADER)
        GL.glShaderSource(v_shader, self.vertex_shader_sb)
        GL.glCompileShader(v_shader)
        f_shader = GL.glCreateShader(GL.GL_FRAGMENT_SHADER)
        GL.glShaderSource(f_shader, self.fragment_shader_sb)
        GL.glCompileShader(f_shader)
        self.selection_box_program = GL.glCreateProgram()
        GL.glAttachShader(self.selection_box_program, v_shader)
        GL.glAttachShader(self.selection_box_program, f_shader)
        GL.glLinkProgram(self.selection_box_program)
        if GL.glGetShaderiv(v_shader, GL.GL_COMPILE_STATUS) != GL.GL_TRUE:
            dprint("Error compiling the shader: ", "GL_VERTEX_SHADER")
            raise RuntimeError(GL.glGetShaderInfoLog(v_shader))
        if GL.glGetShaderiv(f_shader, GL.GL_COMPILE_STATUS) != GL.GL_TRUE:
            dprint("Error compiling the shader: ", "GL_FRAGMENT_SHADER")
            raise RuntimeError(GL.glGetShaderInfoLog(f_shader))
        return True
    
    def update_points(self):
        """ Updates the points that form the selection box with information from
            the start and end coordinates, e.g. imagine the start point (sX,sY)
            and the end (eX,eY), we create the other points as follow:
            
            (sX,sY)                              (sX,sY)   <---   (eX,sY)
                                        ==>         |                |  
                            (eX,eY)              (sX,eY)   --->   (eX,eY)
        """
        assert(self.start is not None)
        assert(self.end is not None)
        self.points = np.array([self.start[0], self.start[1], self.start[0], self.end[1],
                                self.end[0], self.end[1], self.end[0], self.start[1],
                                self.start[0], self.start[1]], dtype=np.float32)
    
    def _make_gl_selection_box(self):
        """ Creates the Vertex Array Object, Vertex Buffer Objects and fill the
            shaders with the data of the corresponding points. Since we will use
            two kind of drawing methods, lines and triangles, the index buffer
            will set to GL_DYNAMIC_DRAW, so as the coordinates, but not the
            colors, which will be set as GL_STATIC_DRAW.
        """
        self.vao = GL.glGenVertexArrays(1)
        GL.glBindVertexArray(self.vao)
        
        ind_vbo = GL.glGenBuffers(1)
        GL.glBindBuffer(GL.GL_ELEMENT_ARRAY_BUFFER, ind_vbo)
        GL.glBufferData(GL.GL_ELEMENT_ARRAY_BUFFER, self.indices.itemsize*int(len(self.indices)), self.indices, GL.GL_DYNAMIC_DRAW)
        
        coord_vbo = GL.glGenBuffers(1)
        GL.glBindBuffer(GL.GL_ARRAY_BUFFER, coord_vbo)
        GL.glBufferData(GL.GL_ARRAY_BUFFER, self.points.itemsize*int(len(self.points)), self.points, GL.GL_DYNAMIC_DRAW)
        att_position = GL.glGetAttribLocation(self.selection_box_program, "vert_coord")
        GL.glEnableVertexAttribArray(att_position)
        GL.glVertexAttribPointer(att_position, 2, GL.GL_FLOAT, GL.GL_FALSE, 2*self.points.itemsize, ctypes.c_void_p(0))
        
        col_vbo = GL.glGenBuffers(1)
        GL.glBindBuffer(GL.GL_ARRAY_BUFFER, col_vbo)
        GL.glBufferData(GL.GL_ARRAY_BUFFER, self.color.itemsize*int(len(self.color)), self.color, GL.GL_STATIC_DRAW)
        att_colors = GL.glGetAttribLocation(self.selection_box_program, "vert_color")
        GL.glEnableVertexAttribArray(att_colors)
        GL.glVertexAttribPointer(att_colors, 3, GL.GL_FLOAT, GL.GL_FALSE, 3*self.color.itemsize, ctypes.c_void_p(0))
        
        # [EN] BUG FIX (regression from an earlier "macOS core-profile
        # fix" pass): this VAO (self.vao) is only ever re-bound for
        # drawing via "GL.glBindVertexArray(self.vao)" in _draw() below
        # -- glEnableVertexAttribArray is called ABOVE, ONCE, during
        # setup, and never again before a draw call. Vertex attribute
        # enable/disable state is part of a VAO's OWN stored state (not
        # global GL state), so calling glDisableVertexAttribArray here
        # permanently disables these attributes on THIS vao, and nothing
        # ever re-enables them: every subsequent draw silently renders
        # nothing (the selection box). There was never a good reason to
        # disable them here -- glBindVertexArray(0) below already fully
        # detaches this VAO's state from whatever gets bound/drawn next.
        GL.glBindVertexArray(0)
        GL.glBindBuffer(GL.GL_ARRAY_BUFFER, 0)
        GL.glBindBuffer(GL.GL_ELEMENT_ARRAY_BUFFER, 0)
        
        self.buffers = (ind_vbo, coord_vbo, col_vbo)
    
    def _draw(self):
        """ The drawing method for the selection box. Initially we will draw the
            box boundaries with lines, and then fill the interior with two
            triangles.
            
            IMPORTANT!!!
            THIS FUNCTION MUST BE CALLED ONLY WHEN AN OPENGL CONTEXT WINDOW HAS
            BEEN CREATED AND INITIALIZED, OTHERWISE WILL RAISE AN ERROR IN THE
            OPENGL WRAPPER!!!
            YOU HAVE BEEN WARNED
        """
        GL.glUseProgram(self.selection_box_program)
        GL.glLineWidth(1)
        GL.glBindVertexArray(self.vao)
        GL.glBindBuffer(GL.GL_ARRAY_BUFFER, self.buffers[1])
        GL.glBufferData(GL.GL_ARRAY_BUFFER, self.points.itemsize*int(len(self.points)), self.points, GL.GL_DYNAMIC_DRAW)
        GL.glBindBuffer(GL.GL_ARRAY_BUFFER, self.buffers[0])
        GL.glBufferData(GL.GL_ARRAY_BUFFER, self.indices.itemsize*int(len(self.indices)), self.indices, GL.GL_DYNAMIC_DRAW)
        GL.glDrawElements(GL.GL_LINE_STRIP, int(len(self.indices)), GL.GL_UNSIGNED_INT, None)
        GL.glEnable(GL.GL_DEPTH_TEST)
        GL.glEnable(GL.GL_BLEND)
        GL.glBlendFunc(GL.GL_SRC_ALPHA, GL.GL_SRC_ALPHA)
        GL.glBindBuffer(GL.GL_ARRAY_BUFFER, self.buffers[0])
        GL.glBufferData(GL.GL_ARRAY_BUFFER, self.triangles.itemsize*int(len(self.triangles)), self.triangles, GL.GL_DYNAMIC_DRAW)
        GL.glDrawElements(GL.GL_TRIANGLE_STRIP, int(len(self.triangles)), GL.GL_UNSIGNED_INT, None)
        GL.glDisable(GL.GL_BLEND)
        GL.glDisable(GL.GL_DEPTH_TEST)
        GL.glBindVertexArray(0)
        GL.glUseProgram(0)
        return True
    
