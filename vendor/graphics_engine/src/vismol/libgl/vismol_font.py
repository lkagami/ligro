#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
#  vismol_font.py
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

from vismol.utils.debug import dprint
import numpy as np
import freetype as ft
import ctypes
from OpenGL import GL
import os
import vismol.libgl.glaxis as glaxis
_fontdir = os.path.split(glaxis.__file__)[:-1]
FONTS_DIR = os.path.join(*_fontdir, "fonts")
#fontpath = os.path.join(*fontpath, "fonts", "VeraMono.ttf")
fontpath = os.path.join(FONTS_DIR, "Amiko-SemiBold.ttf")
DEFAULT_FONT_FILE = "Amiko-SemiBold.ttf"
DEFAULT_FONT_SIZE = 0.35
# [EN] The freetype geometry shader now scales glyph size/advance by the
# label's actual distance from the camera, so that labels keep a
# CONSTANT size on screen regardless of zoom by default (see the
# comment block in shaders/vm_freetype.py for the full reasoning, and
# VismolFont.zoom_sensitivity below to make a font's labels respond to
# zoom again, partially or fully). That distance-scaling needs to be
# calibrated against a reference distance, or every "char_width"/
# "char_height" value that used to mean "this many world units" (tuned
# by eye, including DEFAULT_FONT_SIZE above, back when label size was
# NOT distance-scaled) would suddenly render at the wrong size at the
# app's typical viewing distance. LABEL_DEPTH_REFERENCE is that
# calibration distance (sent to the shader as the "depth_ref" uniform):
# it matches GLCamera's default starting position (pos=(0, 0, 10), i.e.
# 10 world units from the origin -- see glcamera.py), so a freshly
# opened glArea looks the same size as before, no matter what
# zoom_sensitivity is set to.
LABEL_DEPTH_REFERENCE = 10.0


def list_available_fonts():
    """ Returns a sorted list of the .ttf font filenames bundled with
        VisMol (found in the libgl/fonts folder), e.g.
        ['Amiko-Bold.ttf', 'Amiko-Regular.ttf', ...].
    """
    try:
        fonts = [f for f in os.listdir(FONTS_DIR) if f.lower().endswith(".ttf")]
    except OSError:
        fonts = [DEFAULT_FONT_FILE]
    return sorted(fonts)


def resolve_font_path(font_name_or_path):
    """ Resolves a font "name" (as stored in the preferences, e.g.
        'Amiko-SemiBold.ttf') to a full path. If a full/relative path is
        already given (and exists) it is returned unchanged. Falls back to
        the default bundled font when the requested one can't be found.
    """
    if not font_name_or_path:
        return fontpath
    if os.path.isabs(font_name_or_path) and os.path.isfile(font_name_or_path):
        return font_name_or_path
    candidate = os.path.join(FONTS_DIR, os.path.basename(font_name_or_path))
    if os.path.isfile(candidate):
        return candidate
    return fontpath


class VismolFont():
    """ VismolFont stores the data created using the freetype python binding
        library, such as filename, character width, character height, character
        resolution, font color, etc.
    """
    
    def __init__(self, vismol_object=None, font_file=fontpath, char_res=64,
                 char_width=0.35, char_height=0.35, color=None):
                            
        """ Class initialiser
        """
        if color is None:
            color = [1, 1, 1, 1]
        self.vm_object = vismol_object
        self.font_file = font_file
        self.char_res = char_res
        self.char_width = char_width
        self.char_height = char_height
        self.offset = np.array([char_width/1.5, char_height/1.5], dtype=np.float32)
        self.color = np.array(color, dtype=np.float32)
        self.font_buffer = None
        self.texture_id = None
        self.text_u = None
        self.text_v = None
        self.vao = None
        self.text_vbo = None
        self.coord_vbo = None
        # [EN] Per-character slot index (0, 1, 2, ...), consumed by the
        # geometry shader to compute the horizontal advance between
        # glyphs of the same string in screen-aligned (view) space --
        # see the comment in shaders/vm_freetype.py for why this
        # replaced baking the advance into world-space coordinates on
        # the CPU.
        self.char_idx_vbo = None
        # [EN] Small constant (x, y) nudge applied to an entire string,
        # in character-size units, e.g. to shift a label so it doesn't
        # sit exactly on top of the atom/dot it names. Consumed by the
        # "string_shift" uniform. Callers set this right before drawing
        # (see VismolGLCore._draw_text_labels).
        self.string_shift = np.array([0.0, 0.0], dtype=np.float32)
        # [EN] How sensitive THIS font's labels are to camera zoom/dolly,
        # from 0.0 (constant size on screen, regardless of distance --
        # the default, matching the billboard refactor) to 1.0 (natural
        # perspective size, i.e. the label behaves like a fixed-size
        # object in world space and shrinks/grows with camera distance,
        # same as before that refactor). See depth_ref/zoom_sensitivity
        # in shaders/vm_freetype.py for the exact blend. Different
        # VismolFont instances (picking labels vs. distance labels vs.
        # atom labels) can each set their own value.
        self.zoom_sensitivity = 0.0
    
    def make_freetype_font(self):
        """ Function doc
        """
        face = ft.Face(self.font_file)
        face.set_char_size(self.char_res*64)
        # Determine largest glyph size
        width, height, ascender, descender = 0, 0, 0, 0
        for c in range(32,128):
            face.load_char(chr(c), ft.FT_LOAD_RENDER | ft.FT_LOAD_FORCE_AUTOHINT)
            bitmap = face.glyph.bitmap
            width = max(width, bitmap.width)
            ascender = max(ascender, face.glyph.bitmap_top)
            descender = max(descender, bitmap.rows-face.glyph.bitmap_top)
        height = ascender+descender
        # Generate texture data
        #self.font_buffer = np.zeros((height*6, width*16), dtype=np.ubyte)
        self.font_buffer = np.zeros((height*6, width*16), dtype=np.ubyte)
        for j in range(6):
            for i in range(16):
                face.load_char(chr(32+j*16+i), ft.FT_LOAD_RENDER | ft.FT_LOAD_FORCE_AUTOHINT )
                bitmap = face.glyph.bitmap
                x = i*width  + face.glyph.bitmap_left
                y = j*height + ascender - face.glyph.bitmap_top
                self.font_buffer[y:y+bitmap.rows,x:x+bitmap.width].flat = bitmap.buffer
        # Bound texture
        GL.glPixelStorei(GL.GL_UNPACK_ALIGNMENT, 1)
        self.texture_id = GL.glGenTextures(1)
        GL.glBindTexture(GL.GL_TEXTURE_2D, self.texture_id)
        GL.glTexImage2D(GL.GL_TEXTURE_2D, 0, GL.GL_RED, self.font_buffer.shape[1], self.font_buffer.shape[0], 0, GL.GL_RED, GL.GL_UNSIGNED_BYTE, self.font_buffer)
        GL.glTexParameterf(GL.GL_TEXTURE_2D, GL.GL_TEXTURE_MAG_FILTER, GL.GL_LINEAR)
        GL.glTexParameterf(GL.GL_TEXTURE_2D, GL.GL_TEXTURE_MIN_FILTER, GL.GL_LINEAR)
        # Fill the font variables with data
        self.text_u = width/float(self.font_buffer.shape[1])
        self.text_v = height/float(self.font_buffer.shape[0])
    
    def make_freetype_texture(self, program):
        """ Function doc
        """
        coords = np.zeros(3, dtype=np.float32)
        uv_pos = np.zeros(4, dtype=np.float32)
        
        vao = GL.glGenVertexArrays(1)
        GL.glBindVertexArray(vao)
        
        coord_vbo = GL.glGenBuffers(1)
        GL.glBindBuffer(GL.GL_ARRAY_BUFFER, coord_vbo)
        GL.glBufferData(GL.GL_ARRAY_BUFFER, coords.itemsize*len(coords), coords, GL.GL_DYNAMIC_DRAW)
        gl_coord = GL.glGetAttribLocation(program, "vert_coord")
        GL.glEnableVertexAttribArray(gl_coord)
        GL.glVertexAttribPointer(gl_coord, 3, GL.GL_FLOAT, GL.GL_FALSE, 3*coords.itemsize, ctypes.c_void_p(0))
        
        text_vbo = GL.glGenBuffers(1)
        GL.glBindBuffer(GL.GL_ARRAY_BUFFER, text_vbo)
        GL.glBufferData(GL.GL_ARRAY_BUFFER, uv_pos.itemsize*len(uv_pos), uv_pos, GL.GL_DYNAMIC_DRAW)
        gl_texture = GL.glGetAttribLocation(program, "vert_uv")
        GL.glEnableVertexAttribArray(gl_texture)
        GL.glVertexAttribPointer(gl_texture, 4, GL.GL_FLOAT, GL.GL_FALSE, 4*uv_pos.itemsize, ctypes.c_void_p(0))
        
        # [EN] Third attribute: the character's slot index inside its
        # string (float, one value per point/glyph). See the class-level
        # comment on char_idx_vbo and the geometry shader for why this
        # exists -- it lets the advance between glyphs be computed in
        # screen-aligned space on the GPU instead of world space on the
        # CPU.
        char_idx = np.zeros(1, dtype=np.float32)
        char_idx_vbo = GL.glGenBuffers(1)
        GL.glBindBuffer(GL.GL_ARRAY_BUFFER, char_idx_vbo)
        GL.glBufferData(GL.GL_ARRAY_BUFFER, char_idx.itemsize*len(char_idx), char_idx, GL.GL_DYNAMIC_DRAW)
        gl_char_idx = GL.glGetAttribLocation(program, "vert_char_idx")
        GL.glEnableVertexAttribArray(gl_char_idx)
        GL.glVertexAttribPointer(gl_char_idx, 1, GL.GL_FLOAT, GL.GL_FALSE, char_idx.itemsize, ctypes.c_void_p(0))
        
        # [EN] BUG FIX (regression from an earlier "macOS core-profile
        # fix" pass): this VAO is only ever re-bound for drawing via
        # "GL.glBindVertexArray(self.vm_object.vm_font.vao)" in
        # representations.py -- glEnableVertexAttribArray is called
        # ABOVE, ONCE, during setup, and never again before a draw call.
        # Vertex attribute enable/disable state is part of a VAO's OWN
        # stored state (not global GL state), so calling
        # glDisableVertexAttribArray here -- REGARDLESS of whether it
        # happens before or after unbinding -- permanently disables
        # these attributes on THIS vao, and nothing ever re-enables them:
        # every subsequent draw silently renders nothing (no crash, no
        # error -- the attributes are just off). This is what broke atom
        # labels/picking labels/distance labels rendering entirely.
        # There was never a good reason to disable them here in the
        # first place -- glBindVertexArray(0) below already fully
        # detaches this VAO's state from whatever gets bound/drawn next.
        GL.glBindVertexArray(0)
        GL.glBindBuffer(GL.GL_ARRAY_BUFFER, 0)

        self.vao = vao
        self.text_vbo = text_vbo
        self.coord_vbo = coord_vbo
        self.char_idx_vbo = char_idx_vbo
    
    def _get_uniform_location(self, program, name):
        """ Cached wrapper around glGetUniformLocation. Uniform locations
            don't change once a program is linked, so we look each one up
            once per (program, name) pair instead of every frame.
        """
        if not hasattr(self, "_uniform_loc_cache"):
            self._uniform_loc_cache = {}
        key = (program, name)
        loc = self._uniform_loc_cache.get(key)
        if loc is None:
            loc = GL.glGetUniformLocation(program, name)
            self._uniform_loc_cache[key] = loc
        return loc
    
    def load_matrices(self, program, view_mat, proj_mat):
        """ Function doc """
        view = self._get_uniform_location(program, "view_mat")
        GL.glUniformMatrix4fv(view, 1, GL.GL_FALSE, view_mat)
        proj = self._get_uniform_location(program, "proj_mat")
        GL.glUniformMatrix4fv(proj, 1, GL.GL_FALSE, proj_mat)
    
    def load_font_params(self, program):
        """ Loads the uniform parameters for the OpenGL program, such as the
            offset coordinates (X,Y) to calculate the quad and the color of
            the font. Also loads char_advance (the per-glyph horizontal
            step within a string), string_shift (a small constant
            per-string nudge), depth_ref (the calibration distance) and
            zoom_sensitivity (how much this font's size responds to
            camera distance) -- all consumed by the geometry shader to
            build camera-facing labels. See the comments in
            shaders/vm_freetype.py for the full reasoning.
        """
        offset = self._get_uniform_location(program, "offset")
        GL.glUniform2fv(offset, 1, self.offset)
        color = self._get_uniform_location(program, "text_color")
        GL.glUniform4fv(color, 1, self.color)
        char_advance = self._get_uniform_location(program, "char_advance")
        GL.glUniform1f(char_advance, self.char_width)
        string_shift = self._get_uniform_location(program, "string_shift")
        GL.glUniform2fv(string_shift, 1, self.string_shift)
        depth_ref = self._get_uniform_location(program, "depth_ref")
        GL.glUniform1f(depth_ref, LABEL_DEPTH_REFERENCE)
        zoom_sensitivity = self._get_uniform_location(program, "zoom_sensitivity")
        GL.glUniform1f(zoom_sensitivity, self.zoom_sensitivity)
        return True
    
    def print_all(self):
        """ Function created only with debuging purposes.
        """
        dprint("#############################################")
        dprint(self.font_file, "font_file")
        dprint(self.char_res, "char_res")
        dprint(self.char_width, "char_width")
        dprint(self.char_height, "char_height")
        dprint(self.offset, "offset")
        dprint(self.color, "color")
        dprint(self.font_buffer, "font_buffer")
        dprint(self.texture_id, "texture_id")
        dprint(self.text_u, "text_u")
        dprint(self.text_v, "text_v")
        dprint(self.vao, "vao")
        dprint(self.text_vbo, "text_vbo")
        dprint(self.coord_vbo, "coord_vbo")
    
    
    def set_dimensions (self, width, height ):
        """ Function doc """
        self.char_width  =  width
        self.char_height =  height
        self.offset = np.array([ width/2.0,  height/2.0], dtype=np.float32)
    
    def apply_settings(self, font_file=None, size=None):
        """ Updates the font file and/or the character size (width/height
            are kept equal, driven by a single "size" value) and marks the
            OpenGL texture/VAO for regeneration on the next draw call
            (vao=None). Used by the Preferences window to let the user
            change the font family and font size used for every label
            drawn in the glArea (atom labels, picking labels, distance
            labels).
        """
        if font_file is not None:
            self.font_file = resolve_font_path(font_file)
        if size is not None:
            self.set_dimensions(size, size)
        # Force the texture/VAO/font bitmap to be rebuilt on next use
        self.vao = None
        self.font_buffer = None
        self.texture_id = None
    
    def set_color (self, r = 1.0, g = 1.0, b = 1.0):      
        self.color = np.array([r,g,b], dtype=np.float32)
    
    
    def draw_labels(self):
        """ Function doc """
        pass
