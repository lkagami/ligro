#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
#  representations.py
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
#  update Bachega 2026

from vismol.utils.debug import dprint
import ctypes
import numpy as np
from OpenGL import GL
from logging import getLogger
from vismol.libgl.vismol_font import VismolFont
from vismol.libgl.vismol_font import resolve_font_path, DEFAULT_FONT_FILE, DEFAULT_FONT_SIZE

logger = getLogger(__name__)


def compute_atom_label_text(atom, content):
    """ Computes the text to display for `atom`'s label given a `content`
        selector -- one of 'name', 'symbol', 'index', 'mm_charge',
        'residue_name', 'residue_index', 'chain'. Mirrors the per-atom
        logic already used by the atom-labeling context menu (see
        gui/eSession.py, menu_show_atom_name/.../menu_show_chain), so the
        Preferences window ("Atom Labels (glArea)") and that menu always
        agree on how each option is computed. Never raises: falls back to
        the atom name (or an empty string) if anything is missing (e.g.
        no residue link, or no MM charge available for this atom yet).
    """
    try:
        if content == 'symbol':
            return atom.symbol
        elif content == 'index':
            return str(atom.index)
        elif content == 'mm_charge':
            vm_object = atom.vm_object
            p_session = vm_object.vm_session.main.p_session
            charge = float(p_session.psystem[vm_object.e_id].mmState.charges[atom.index - 1])
            return '%4.3f' % charge
        elif content == 'residue_name':
            return atom.residue.name
        elif content == 'residue_index':
            return str(atom.residue.index)
        elif content == 'chain':
            return str(atom.chain.name)
        else:
            # 'name' (default) and any unrecognized value
            return atom.name
    except Exception:
        return atom.name or ''


class Representation:
    """ Class doc """
    
    def __init__ (self, vismol_object, vismol_glcore, name, active, indexes, is_dynamic = False):
        self.vm_object = vismol_object
        self.vm_session = vismol_object.vm_session
        self.vm_glcore = vismol_glcore
        self.name = name
        self.active = active
        # indexes pode chegar como None (ex: build_core_representations() em
        # vismol_object.py passa indexes=self.index_bonds pra
        # DashedLinesRepresentation, e index_bonds e None por padrao pra
        # objetos sem atomos/ligacoes -- como os VismolObject "so-superficie"
        # criados em surface_analysis_window.py pra orbital/potential/
        # density/MEP/external cube. np.array(None, dtype=uint32) da
        # TypeError em vez de um array vazio, entao tratamos isso aqui em
        # vez de em cada chamador individualmente.
        # [EN] indexes may arrive as None here (e.g. build_core_representations()
        # in vismol_object.py passes indexes=self.index_bonds to
        # DashedLinesRepresentation, and index_bonds defaults to None for
        # objects with no atoms/bonds -- such as the surface-only VismolObject
        # instances created in surface_analysis_window.py for orbital/potential/
        # density/MEP/external-cube surfaces). np.array(None, dtype=uint32)
        # raises TypeError instead of producing an empty array, so it is
        # normalised here once, rather than at every call site.
        if indexes is None:
            indexes = []
        self.indexes = np.array(indexes, dtype=np.uint32)
        self.elements = np.uint32(self.indexes.shape[0])
        
        self.is_dynamic = is_dynamic

        # Marca se esta representacao usa uma cor UNIFORME/fixa (nao
        # segue vobject.colors, o array de cores por-atomo do objeto).
        # Codigo generico que reaplica cores por atomo em todas as
        # representacoes de um vobject (ex: set_color_by_index em
        # eSession.py) deve checar essa flag e pular representacoes
        # marcadas, ou vai sobrescrever a cor fixa delas.
        self.uses_uniform_color = False

        self.was_rep_modified = False
        self.was_sel_modified = False
        self.was_col_modified = False
        self.was_rep_coord_modified = False
        self.was_sel_coord_modified = False
        self.was_rep_ind_modified = False
        self.was_sel_ind_modified = False
        self.was_rep_col_modified = False
        # Cache do ultimo conjunto de IDs de selecao enviado ao VBO de
        # picking_dots. Evita reconstruir indices e re-subir o buffer todo
        # frame quando a selecao nao mudou (Gargalo 1). None = nunca subiu.
        self._last_uploaded_sel_ids = None
        # representation
        self.vao = None
        self.ind_vbo = None
        self.coord_vbo = None
        self.col_vbo = None
        self.size_vbo = None
        # selection
        self.sel_vao = None
        self.sel_ind_vbo = None
        self.sel_coord_vbo = None
        self.sel_col_vbo = None
        self.sel_size_vbo = None
        # shaders
        self.shader_program = None
        self.sel_shader_program = None
    
    def _check_vao_and_vbos(self):
        #print(self.name)
        self.shader_program = self.vm_glcore.shader_programs[self.name]
        self.sel_shader_program = self.vm_glcore.shader_programs[self.name + "_sel"]
        if self.vao is None:
            self._make_gl_representation_vao_and_vbos()
        if self.sel_vao is None:
            self._make_gl_sel_representation_vao_and_vbos()
    
    def _make_gl_representation_vao_and_vbos(self):
        """ Function doc """
        logger.debug("building '{}' representation VAO and VBOs".format(self.name))
        self.vao = self._make_gl_vao()
        self.ind_vbo = self._make_gl_index_buffer(self.indexes)
        self.coord_vbo = self._make_gl_coord_buffer(self.vm_object.frames[0], self.shader_program)
        self.col_vbo = self._make_gl_color_buffer(self.vm_object.colors, self.shader_program)
    
    def _make_gl_sel_representation_vao_and_vbos(self):
        """ Function doc """
        logger.debug("building '{}' background selection VAO and VBOs".format(self.name))
        self.sel_vao = self._make_gl_vao()
        self.sel_ind_vbo = self._make_gl_index_buffer(self.indexes)
        self.sel_coord_vbo = self._make_gl_coord_buffer(self.vm_object.frames[0], self.sel_shader_program)
        self.sel_col_vbo = self._make_gl_color_buffer(self.vm_object.color_indexes, self.sel_shader_program)
    
    def _make_gl_vao(self):
        """ Function doc """
        vao = GL.glGenVertexArrays(1)
        GL.glBindVertexArray(vao)
        return vao
    
    def _make_gl_index_buffer(self, indexes):
        """ Function doc """
        ind_vbo = GL.glGenBuffers(1)
        GL.glBindBuffer(GL.GL_ELEMENT_ARRAY_BUFFER, ind_vbo)
        GL.glBufferData(GL.GL_ELEMENT_ARRAY_BUFFER, indexes.nbytes, indexes, GL.GL_DYNAMIC_DRAW)
        return ind_vbo
    
    def _make_gl_coord_buffer(self, coords, program):
        """ Function doc """
        coord_vbo = GL.glGenBuffers(1)
        GL.glBindBuffer(GL.GL_ARRAY_BUFFER, coord_vbo)
        GL.glBufferData(GL.GL_ARRAY_BUFFER, coords.nbytes, coords, GL.GL_STATIC_DRAW)
        att_position = GL.glGetAttribLocation(program, "vert_coord")
        GL.glEnableVertexAttribArray(att_position)
        GL.glVertexAttribPointer(att_position, 3, GL.GL_FLOAT, GL.GL_FALSE, 3*coords.itemsize, ctypes.c_void_p(0))
        return coord_vbo
    
    def _make_gl_color_buffer(self, colors, program, instances=False):
        """ Function doc """
        col_vbo = GL.glGenBuffers(1)
        GL.glBindBuffer(GL.GL_ARRAY_BUFFER, col_vbo)
        GL.glBufferData(GL.GL_ARRAY_BUFFER, colors.nbytes, colors, GL.GL_STATIC_DRAW)
        att_colors = GL.glGetAttribLocation(program, "vert_color")
        GL.glEnableVertexAttribArray(att_colors)
        GL.glVertexAttribPointer(att_colors, 3, GL.GL_FLOAT, GL.GL_FALSE, 3*colors.itemsize, ctypes.c_void_p(0))
        if instances:
            GL.glVertexAttribDivisor(att_colors, 1)
        return col_vbo
    
    def _make_gl_bond_order_buffer(self, orders, program):
        """ Cria o VBO do atributo inteiro 'vert_bond_order' (ordem da ligacao
            por atomo). Usa glVertexAttribIPointer porque o atributo e 'in int'
            no shader -- glVertexAttribPointer (float) faria o driver entregar
            valores convertidos/zerados. Se o shader nao declara o atributo
            (att == -1), nao faz nada e retorna None (representacoes sem
            suporte a ordem de ligacao continuam funcionando). """
        att = GL.glGetAttribLocation(program, "vert_bond_order")
        if att == -1:
            return None
        orders = np.ascontiguousarray(orders, dtype=np.int32)
        bo_vbo = GL.glGenBuffers(1)
        GL.glBindBuffer(GL.GL_ARRAY_BUFFER, bo_vbo)
        GL.glBufferData(GL.GL_ARRAY_BUFFER, orders.nbytes, orders, GL.GL_STATIC_DRAW)
        GL.glEnableVertexAttribArray(att)
        GL.glVertexAttribIPointer(att, 1, GL.GL_INT, 0, ctypes.c_void_p(0))
        return bo_vbo
    
    def _make_gl_bond_order_texture_buffer(self, orders):
        """ Cria uma Buffer Texture (GL_TEXTURE_BUFFER, formato R32UI) com uma
            entrada de ordem de ligacao POR PRIMITIVA (nao por atomo/vertice).

            Diferenca-chave em relacao ao antigo _make_gl_bond_order_buffer:
            aquele era um atributo por-VERTICE, indexado pelo GL_ELEMENT_ARRAY_
            BUFFER (index_bonds) -- ou seja, por ATOMO. Um atomo com duas
            ligacoes de ordens diferentes so guardava UM valor (a solucao
            antiga usava a maior ordem entre elas, o que na pratica parecia
            "valencia maxima do atomo" em vez da ordem real de cada ligacao).

            Aqui, 'orders' e simplesmente bond_order_list (uma ordem por
            ligacao, na MESMA sequencia dos pares em index_bonds). No
            geometry shader, o texel k e lido via gl_PrimitiveIDIn, que da o
            indice da primitiva GL_LINES (== o k-esimo par de index_bonds)
            sendo desenhada. Isso pareia ordem <-> par de index_bonds 1:1,
            sem nenhuma ambiguidade por atomo compartilhado.

            Retorna (buffer_id, texture_id). Ambos devem ser mantidos vivos
            (guardados em self) e liberados/recriados quando o numero de
            ligacoes mudar (ex.: Dynamic Bonds trocando de frame).
        """
        orders = np.ascontiguousarray(orders, dtype=np.uint32)
        if orders.size == 0:
            # Buffer Texture vazia (0 bytes) e invalida em alguns drivers;
            # garante ao menos 1 elemento dummy.
            orders = np.zeros(1, dtype=np.uint32)
        tbo_buf = GL.glGenBuffers(1)
        GL.glBindBuffer(GL.GL_TEXTURE_BUFFER, tbo_buf)
        GL.glBufferData(GL.GL_TEXTURE_BUFFER, orders.nbytes, orders, GL.GL_DYNAMIC_DRAW)
        tex_id = GL.glGenTextures(1)
        GL.glBindTexture(GL.GL_TEXTURE_BUFFER, tex_id)
        GL.glTexBuffer(GL.GL_TEXTURE_BUFFER, GL.GL_R32UI, tbo_buf)
        GL.glBindTexture(GL.GL_TEXTURE_BUFFER, 0)
        GL.glBindBuffer(GL.GL_TEXTURE_BUFFER, 0)
        return tbo_buf, tex_id

    def _update_gl_bond_order_texture_buffer(self, tbo_buf, orders):
        """ Reenvia o conteudo de uma Buffer Texture ja existente (criada por
            _make_gl_bond_order_texture_buffer), sem recriar buffer/textura.
            Usado quando bond_order_list muda mas o numero de ligacoes e
            compativel com o buffer atual (ex.: recarregar a mesma molecula).
            Se o TAMANHO em bytes mudar, o chamador deve recriar do zero
            (ver _refresh_bond_order_tbo em SticksRepresentation). """
        orders = np.ascontiguousarray(orders, dtype=np.uint32)
        if orders.size == 0:
            orders = np.zeros(1, dtype=np.uint32)
        GL.glBindBuffer(GL.GL_TEXTURE_BUFFER, tbo_buf)
        GL.glBufferData(GL.GL_TEXTURE_BUFFER, orders.nbytes, orders, GL.GL_DYNAMIC_DRAW)
        GL.glBindBuffer(GL.GL_TEXTURE_BUFFER, 0)

    def _make_gl_radius_buffer(self, radii, program, instances=False):
        """ Function doc """
        rad_vbo = GL.glGenBuffers(1)
        GL.glBindBuffer(GL.GL_ARRAY_BUFFER, rad_vbo)
        GL.glBufferData(GL.GL_ARRAY_BUFFER, radii.nbytes, radii, GL.GL_STATIC_DRAW)
        att_rads = GL.glGetAttribLocation(program, "vert_radius")
        GL.glEnableVertexAttribArray(att_rads)
        GL.glVertexAttribPointer(att_rads, 1, GL.GL_FLOAT, GL.GL_FALSE, radii.itemsize, ctypes.c_void_p(0))
        if instances:
            GL.glVertexAttribDivisor(att_rads, 1)
        return rad_vbo
    
    def _make_gl_instance_buffer(self, instances, program):
        """ Function doc """
        insta_vbo = GL.glGenBuffers(1)
        GL.glBindBuffer(GL.GL_ARRAY_BUFFER, insta_vbo)
        GL.glBufferData(GL.GL_ARRAY_BUFFER, instances.nbytes, instances, GL.GL_STATIC_DRAW)
        gl_insta = GL.glGetAttribLocation(program, "vert_instance")
        GL.glEnableVertexAttribArray(gl_insta)
        GL.glVertexAttribPointer(gl_insta, 3, GL.GL_FLOAT, GL.GL_FALSE, 0, ctypes.c_void_p(0))
        GL.glVertexAttribDivisor(gl_insta, 1)
        return insta_vbo
    
    def _make_gl_impostor_buffer(self, impostors_radii, program):
        """ Function doc """
        size_vbo = GL.glGenBuffers(1)
        GL.glBindBuffer(GL.GL_ARRAY_BUFFER, size_vbo)
        GL.glBufferData(GL.GL_ARRAY_BUFFER, impostors_radii.nbytes, impostors_radii, GL.GL_STATIC_DRAW)
        att_size = GL.glGetAttribLocation(program, "vert_dot_size")
        GL.glEnableVertexAttribArray(att_size)
        GL.glVertexAttribPointer(att_size, 1, GL.GL_FLOAT, GL.GL_FALSE, impostors_radii.itemsize, ctypes.c_void_p(0))
        
        self.ratio = self.vm_glcore.width / self.vm_glcore.height
        ratio_vbo = 1
        # ratio = np.repeat(self.ratio, impostors_radii.shape[0])
        # ratio_vbo = GL.glGenBuffers(1)
        # GL.glBindBuffer(GL.GL_ARRAY_BUFFER, ratio_vbo)
        # GL.glBufferData(GL.GL_ARRAY_BUFFER, ratio.nbytes, ratio, GL.GL_STATIC_DRAW)
        # att_ratio = GL.glGetAttribLocation(program, "hw_ratio")
        # GL.glEnableVertexAttribArray(att_ratio)
        # GL.glVertexAttribPointer(att_ratio, 1, GL.GL_FLOAT, GL.GL_FALSE, ratio.itemsize, ctypes.c_void_p(0))
        return size_vbo, ratio_vbo
    
    def _load_coord_vbo(self, coord_vbo=False, sel_coord_vbo=False):
        """ This function assigns the coordinates to 
        be drawn by the function  draw_representation
        
        NOTE on usage hint: this path is re-executed every time the
        coordinates change (e.g. trajectory/MD playback), unlike
        _make_gl_coord_buffer which only runs once at VBO creation.
        GL_DYNAMIC_DRAW tells the driver this buffer's contents are
        updated frequently, which avoids the reallocation overhead
        GL_STATIC_DRAW can trigger on repeated glBufferData calls to
        the same buffer (the same byte size is reused here, so the
        driver can update in place instead of reallocating storage).
        """
        frame, f = self.vm_glcore._safe_frame_coords(self.vm_object)
        if coord_vbo:
            GL.glBindBuffer(GL.GL_ARRAY_BUFFER, self.coord_vbo)
            GL.glBufferData(GL.GL_ARRAY_BUFFER, frame.nbytes, frame, GL.GL_DYNAMIC_DRAW)
        
        if sel_coord_vbo:
            GL.glBindBuffer(GL.GL_ARRAY_BUFFER, self.sel_coord_vbo)
            GL.glBufferData(GL.GL_ARRAY_BUFFER, frame.nbytes, frame, GL.GL_DYNAMIC_DRAW)
    
    def _load_ind_vbo(self, ind_vbo=False, sel_ind_vbo=False):
        """ Function doc """
        if self.is_dynamic:
            frame, f = self.vm_glcore._safe_frame_coords(self.vm_object)
            # [EN] BUG FIX (reported by the user: "IndexError: list index
            # out of range" here while viewing a Normal Modes trajectory
            # -- root cause was vismol_object.dynamic_bonds being built
            # for the wrong (much smaller) frame count at representation-
            # setup time, now fixed at the source in import_trajectory.py
            # 's _import_normal_modes_data()). Guarding here too: dynamic_
            # bonds is a plain Python list kept in sync with frames by
            # whoever last called vm_session.define_dynamic_bonds(), so
            # nothing here can guarantee it never falls behind again --
            # falling back to the last available frame's bonds (like
            # _safe_frame_coords already does for out-of-range frame
            # indices) is a harmless, better failure mode than crashing
            # the whole render loop.
            if len(self.vm_object.dynamic_bonds) > 0:
                f = min(f, len(self.vm_object.dynamic_bonds) - 1)
                self.define_new_indexes_to_vbo(input_indexes = self.vm_object.dynamic_bonds[f])

        if ind_vbo:
            GL.glBindBuffer(GL.GL_ELEMENT_ARRAY_BUFFER, self.ind_vbo)
            GL.glBufferData(GL.GL_ELEMENT_ARRAY_BUFFER, self.indexes.nbytes, self.indexes, GL.GL_DYNAMIC_DRAW)
        
        if sel_ind_vbo:
            GL.glBindBuffer(GL.GL_ELEMENT_ARRAY_BUFFER, self.sel_ind_vbo)
            GL.glBufferData(GL.GL_ELEMENT_ARRAY_BUFFER, self.indexes.nbytes, self.indexes, GL.GL_DYNAMIC_DRAW)
    
    def _load_color_vbo(self, colors = None):
        """ This function assigns the colors to
            be drawn by the function  draw_representation"""
        if colors is None:
            colors = self.colors
        GL.glBindBuffer(GL.GL_ARRAY_BUFFER, self.col_vbo)
        GL.glBufferData(GL.GL_ARRAY_BUFFER, colors.nbytes, colors, GL.GL_STATIC_DRAW)

            
    def _enable_anti_alias_to_lines(self):
        """ Function doc """
        GL.glEnable(GL.GL_DEPTH_TEST)
        GL.glEnable(GL.GL_BLEND)
        GL.glBlendFunc(GL.GL_SRC_ALPHA, GL.GL_ONE_MINUS_SRC_ALPHA)
        GL.glEnable(GL.GL_LINE_SMOOTH)
        GL.glHint(GL.GL_LINE_SMOOTH_HINT, GL.GL_NICEST)
    
    def _disable_anti_alias_to_lines(self):
        """ Function doc """
        GL.glDisable(GL.GL_LINE_SMOOTH)
        GL.glDisable(GL.GL_BLEND)
        GL.glDisable(GL.GL_DEPTH_TEST)
    
    def define_new_indexes_to_vbo(self, input_indexes):
        """ Function doc """
        self.indexes = np.array(input_indexes, dtype=np.uint32)
        self.elements = np.uint32(self.indexes.shape[0])


class PickingDotsRepresentation(Representation):
    """ Class doc """
    
    #def __init__(self, vismol_object, vismol_glcore, indexes=None, active=True, colors = None):
    def __init__(self, vismol_object, vismol_glcore, indexes=None, active=True, colors = None):
        """ Class initialiser """
        super(PickingDotsRepresentation, self).__init__(vismol_object, vismol_glcore, "picking_dots", active, indexes, colors)
        
        if colors:
            self.colors = colors
        else:
            self.colors = self.vm_session.vm_config.gl_parameters["picking_dots_color"]
        #print(self.colors)
        
    def _check_vao_and_vbos(self):
        self.shader_program = self.vm_glcore.core_shader_programs[self.name]
        if self.vao is None:
            self._make_gl_representation_vao_and_vbos()
    
    def _make_gl_representation_vao_and_vbos(self):
        """ Function doc """
        logger.debug("building '{}' representation VAO and VBOs".format(self.name))
        self.vao = self._make_gl_vao()
        self.ind_vbo = self._make_gl_index_buffer(self.indexes)
        self.coord_vbo = self._make_gl_coord_buffer(self.vm_object.frames[0], self.shader_program)
        #colors = self.vm_session.vm_config.gl_parameters["picking_dots_color"] * self.vm_object.frames.shape[1]
        colors = self.colors * self.vm_object.frames.shape[1]
        
        colors = np.array(colors, dtype=np.float32).reshape([self.vm_object.frames.shape[1], 3])
        self.col_vbo = self._make_gl_color_buffer(colors, self.shader_program)
    
    def draw_representation(self):
        """ Function doc """
        self._check_vao_and_vbos()
        _size = self.vm_session.vm_config.gl_parameters["dot_sel_size"]
        _size = _size * self.vm_glcore.height / (abs(self.vm_glcore.dist_cam_zrp)) / 2 
        

        #How to pass data to shader on the fly. It's not the most efficient mode, but it's an okay solution in this case.
        '''
        Using the "GL.glPointSize(_size)" function did not give 
        consistent results between different versions of openGL, 
        or linux distros. The solution, for now, was to pass the 
        pixel size value directly to the shader code.
        '''
        GL.glUseProgram(self.shader_program)
        custom_int_location = self.vm_glcore._get_uniform_location(self.shader_program, "_size")
        if custom_int_location != -1:
            GL.glUniform1i(custom_int_location, int(_size))  # Set the integer value to  _size
        else:
            #print("Uniform '_size' not found in shader program.")
            #print(int(_size))
            #GL.glPointSize(int(_size+256))
            GL.glPointSize(int(_size))
        
        GL.glEnable(GL.GL_VERTEX_PROGRAM_POINT_SIZE)
        self.vm_glcore.load_matrices(self.shader_program, self.vm_object.model_mat)
        GL.glBindVertexArray(self.vao)
        
        if self.was_rep_coord_modified:
            self._load_coord_vbo(coord_vbo=True)
            self.was_rep_coord_modified = False
        if self.was_rep_ind_modified:
            self._load_ind_vbo(ind_vbo=True)
            self.was_rep_ind_modified = False
        
        GL.glDrawElements(GL.GL_POINTS, self.elements, GL.GL_UNSIGNED_INT, None)
        
        GL.glBindVertexArray(0)
        GL.glDisable(GL.GL_VERTEX_PROGRAM_POINT_SIZE)
        #GL.glPointSize(1)
        GL.glUseProgram(0)
        GL.glDisable(GL.GL_DEPTH_TEST)
    
    def draw_background_sel_representation(self):
        pass


class OneColorDotsRepresentation(Representation):
    """ Class doc """
    
    def __init__(self, vismol_object, vismol_glcore, indexes, active=True, rgb=None):
        """ Class initialiser """
        super(OneColorDotsRepresentation, self).__init__(
            vismol_object,
            vismol_glcore,
            #"dots",
            'posdot_type',
            active,
            indexes
        )

        self.rgb = rgb

        # Se RGB foi fornecido, sobrescreve as cores e marca a
        # representacao como "cor fixa" (ver uses_uniform_color na
        # classe base) -- assim codigo generico de recoloracao por
        # atomo (ex: set_color_by_index) sabe que deve pular esta
        # representacao ao invés de sobrescrever sua cor.
        if rgb is not None:
            self._set_uniform_color(rgb)
            self.uses_uniform_color = True

    def _set_uniform_color(self, rgb):
        """
        rgb deve ser tupla/lista (R,G,B)
        Aceita valores em [0,255] ou [0,1]

        IMPORTANTE: o desenho e feito por indices (glDrawElements) sobre o
        array de coordenadas COMPLETO do objeto (self.vm_object.frames[0]),
        exatamente como coord_vbo/ind_vbo sao montados na classe base. Logo
        o buffer de cor precisa ter uma entrada POR ATOMO DO OBJETO (mesmo
        tamanho de vm_object.frames[0]/vm_object.colors), e nao apenas
        len(self.indexes) entradas -- caso contrario self.indexes (que
        contem os IDs reais dos atomos, normalmente >> numero de restricoes)
        acessa fora dos limites do buffer, gerando cores inconsistentes/
        indefinidas dependendo do driver/GPU.
        """
        rgb = np.array(rgb, dtype=np.float32)

        # Normaliza se necessário
        if np.any(rgb > 1.0):
            rgb /= 255.0

        n_atoms = self.vm_object.frames[0].shape[0]

        # Repete RGB para TODOS os atomos do objeto (nao so os N pontos
        # restringidos), para casar com o dimensionamento de coord_vbo.
        self.colors = np.tile(rgb, (n_atoms, 1)).astype(np.float32)

        self.was_col_modified = True

    def _make_gl_representation_vao_and_vbos(self):
        """ Sobrescreve a versao da classe base para nao usar
        vm_object.colors (cores por-atomo/por-elemento) quando uma cor
        uniforme foi definida. Sem este override, a primeira criacao do
        VAO usaria as cores normais dos atomos ate was_col_modified
        disparar uma recarga -- ou seja, o primeiro frame desenhado
        ficaria com a cor errada. """
        logger.debug("building '{}' representation VAO and VBOs".format(self.name))
        self.vao = self._make_gl_vao()
        self.ind_vbo = self._make_gl_index_buffer(self.indexes)
        self.coord_vbo = self._make_gl_coord_buffer(self.vm_object.frames[0], self.shader_program)

        if self.rgb is None:
            # Nenhuma cor uniforme foi passada: cai para o comportamento
            # padrao (cor por atomo), igual a classe base.
            self.col_vbo = self._make_gl_color_buffer(self.vm_object.colors, self.shader_program)
        else:
            # self.colors ja foi montado em _set_uniform_color com o
            # tamanho correto (n_atoms), entao so precisa subir pra GPU.
            self.col_vbo = self._make_gl_color_buffer(self.colors, self.shader_program)
    
    def draw_representation(self):
        """ Function doc """
        self._check_vao_and_vbos()
        self._enable_anti_alias_to_lines()
        _size = self.vm_glcore.vm_config.gl_parameters["dots_size"]
        _height = self.vm_glcore.height
        GL.glDisable(GL.GL_DEPTH_TEST)
        GL.glEnable(GL.GL_VERTEX_PROGRAM_POINT_SIZE)
        GL.glUseProgram(self.shader_program)
        point_size = _size * _height / (abs(self.vm_glcore.dist_cam_zrp)) / 2
        GL.glPointSize(point_size)
        # Fix: com GL_VERTEX_PROGRAM_POINT_SIZE habilitado, o glPointSize()
        # acima e ignorado pela especificacao OpenGL - o tamanho real do
        # ponto so vem do gl_PointSize escrito no vertex shader (ver
        # shaders/dots.py). Sem enviar esse uniform, o tamanho fica
        # indefinido pelo driver: em GPUs Intel/Mesa os dots ficavam
        # invisiveis (tamanho efetivamente zero), sem gerar nenhum erro.
        point_size_loc = self.vm_glcore._get_uniform_location(
            self.shader_program, "point_size"
        )
        GL.glUniform1f(point_size_loc, point_size)

        self.vm_glcore.load_matrices(
            self.shader_program,
            self.vm_object.model_mat
        )
        self.vm_glcore.load_fog(self.shader_program)
        self.vm_glcore.load_lights(self.shader_program)
        GL.glBindVertexArray(self.vao)

        if self.was_rep_coord_modified:
            self._load_coord_vbo(coord_vbo=True)
            self.was_rep_coord_modified = False

        if self.was_rep_ind_modified:
            self._load_ind_vbo(ind_vbo=True)
            self.was_rep_ind_modified = False

        if self.was_col_modified:
            self._load_color_vbo(None)
            self.was_col_modified = False

        GL.glDrawElements(
            GL.GL_POINTS,
            self.elements,
            GL.GL_UNSIGNED_INT,
            None
        )

        GL.glBindVertexArray(0)
        self._disable_anti_alias_to_lines()
        GL.glDisable(GL.GL_VERTEX_PROGRAM_POINT_SIZE)
        GL.glPointSize(1)
        GL.glUseProgram(0)
        
        
    def draw_background_sel_representation (self):
        """ Function doc """
        pass
    
    
    
class DotsRepresentation(Representation):
    """ Class doc """
    
    def __init__ (self, vismol_object, vismol_glcore, indexes, active=True):
        """ Class initialiser """
        super(DotsRepresentation, self).__init__(vismol_object, vismol_glcore, "dots", active, indexes)
    
    def draw_representation(self):
        """ Function doc """
        self._check_vao_and_vbos()
        self._enable_anti_alias_to_lines()
        _size = self.vm_glcore.vm_config.gl_parameters["dots_size"]
        _height = self.vm_glcore.height
        GL.glEnable(GL.GL_VERTEX_PROGRAM_POINT_SIZE)
        GL.glUseProgram(self.shader_program)
        GL.glPointSize(_size * _height / (abs(self.vm_glcore.dist_cam_zrp)) / 2)
        self.vm_glcore.load_matrices(self.shader_program, self.vm_object.model_mat)
        self.vm_glcore.load_fog(self.shader_program)
        self.vm_glcore.load_lights(self.shader_program)
        GL.glBindVertexArray(self.vao)
        
        if self.was_rep_coord_modified:
            self._load_coord_vbo(coord_vbo=True)
            self.was_rep_coord_modified = False
        if self.was_rep_ind_modified:
            self._load_ind_vbo(ind_vbo=True)
            self.was_rep_ind_modified = False
        if self.was_col_modified:
            self._load_color_vbo(None)
            self.was_col_modified = False
        
        GL.glDrawElements(GL.GL_POINTS, self.elements, GL.GL_UNSIGNED_INT, None)
        
        GL.glBindVertexArray(0)
        self._disable_anti_alias_to_lines()
        GL.glDisable(GL.GL_VERTEX_PROGRAM_POINT_SIZE)
        GL.glPointSize(1)
        GL.glUseProgram(0)
    
    def draw_background_sel_representation(self):
        """ Function doc """
        self._check_vao_and_vbos()
        _size = self.vm_glcore.vm_config.gl_parameters["dots_size"]
        _height = self.vm_glcore.height
        GL.glEnable(GL.GL_VERTEX_PROGRAM_POINT_SIZE)
        GL.glUseProgram(self.sel_shader_program)
        GL.glPointSize(_size * _height / (abs(self.vm_glcore.dist_cam_zrp)) / 2)
        self.vm_glcore.load_matrices(self.sel_shader_program, self.vm_object.model_mat)
        GL.glBindVertexArray(self.sel_vao)
        
        if self.was_sel_coord_modified:
            self._load_coord_vbo(sel_coord_vbo=True)
            self.was_sel_coord_modified = False
        if self.was_sel_ind_modified:
            self._load_ind_vbo(sel_ind_vbo=True)
            self.was_sel_ind_modified = False
        
        GL.glDrawElements(GL.GL_POINTS, self.elements, GL.GL_UNSIGNED_INT, None)
        
        GL.glBindVertexArray(0)
        GL.glDisable(GL.GL_DEPTH_TEST)
        GL.glPointSize(1)
        GL.glUseProgram(0)


LINES_MULTI_RING_DOUBLE = 0   # u_mode of the "lines_multi" geometry shader
LINES_MULTI_DOUBLE      = 1
LINES_MULTI_TRIPLE      = 2


def multiple_bond_references(pairs, orders, all_pairs, is_hydrogen=None, max_ring=8):
    """ Drawing mode and reference atom of each bond, for the double/triple
        bonds of "lines" and "sticks" (2026-10-14).

        pairs       (n, 2) bonds drawn by the representation
        orders      (n,)   bond orders (1, 2, 3; others are clamped)
        all_pairs   (m, 2) every bond of the object (neighbours, rings)
        is_hydrogen optional callable(atom index) -> bool

        Returns (mode, ref), int64 arrays of length n. mode: -1 single,
        LINES_MULTI_RING_DOUBLE, LINES_MULTI_DOUBLE or LINES_MULTI_TRIPLE.
        ref: for a double bond in a ring of up to max_ring atoms, the
        neighbour of i on the smallest ring (the inner line/cylinder goes on
        its side); for another double bond, a substituent of i or j, heavy
        atoms first (it fixes the plane of the two lines); -1 otherwise
        (triple bonds, no substituent: offset on the screen). """
    pairs = np.asarray(pairs, dtype=np.int64).reshape(-1, 2)
    orders = np.clip(np.asarray(orders, dtype=np.int64).ravel(), 1, 3)
    mode = np.full(pairs.shape[0], -1, dtype=np.int64)
    ref = np.full(pairs.shape[0], -1, dtype=np.int64)
    mode[orders >= 3] = LINES_MULTI_TRIPLE
    doubles = np.flatnonzero(orders == 2)
    if doubles.size == 0:
        return mode, ref
    # . neighbour lists (CSR) of the whole object
    allp = np.asarray(all_pairs, dtype=np.int64).reshape(-1, 2)
    n_atoms = int(max(allp.max(initial=-1), pairs.max(initial=-1))) + 1
    a = np.concatenate([allp[:, 0], allp[:, 1]])
    b = np.concatenate([allp[:, 1], allp[:, 0]])
    srt = np.argsort(a, kind="stable")
    a, b = a[srt], b[srt]
    starts = np.searchsorted(a, np.arange(n_atoms + 1))
    def neighbours(i):
        return b[starts[i]:starts[i + 1]]

    def ring_neighbour(i, j):
        """ First atom after i on the shortest path i -> j that avoids the bond i-j. """
        first = {}
        frontier = []
        for k in neighbours(i):
            k = int(k)
            if k == j or k in first:
                continue
            first[k] = k
            frontier.append(k)
        candidates = []
        for _ in range(max_ring - 2):
            nxt = []
            for x in frontier:
                for y in neighbours(x):
                    y = int(y)
                    if y == j:
                        candidates.append(first[x])
                    elif y != i and y not in first:
                        first[y] = first[x]
                        nxt.append(y)
            if candidates:
                # . several rings of the same size (fused systems): prefer a heavy atom
                if is_hydrogen is not None:
                    heavy = [k for k in candidates if not is_hydrogen(k)]
                    if heavy:
                        return heavy[0]
                return candidates[0]
            frontier = nxt
        return None

    for n in doubles:
        i, j = int(pairs[n, 0]), int(pairs[n, 1])
        k = ring_neighbour(i, j)
        if k is not None:
            mode[n], ref[n] = LINES_MULTI_RING_DOUBLE, k
            continue
        # . substituent plane (only the component of ref perpendicular to the
        #   bond is used, so a substituent of j gives the same plane)
        mode[n] = LINES_MULTI_DOUBLE
        subs = [int(x) for x in neighbours(i) if int(x) != j] + [int(x) for x in neighbours(j) if int(x) != i]
        if subs and is_hydrogen is not None:
            subs.sort(key=lambda s: bool(is_hydrogen(s)))
        if subs:
            ref[n] = subs[0]
    return mode, ref


def multiple_bond_split(pairs, orders, all_pairs, is_hydrogen=None, max_ring=8):
    """ Splits the bonds of a "lines" representation into single bonds and
        double/triple bonds (see multiple_bond_references).

        Returns (singles, groups): singles is the flat uint32 array of the
        single bonds (GL_LINES); groups maps LINES_MULTI_* to a flat uint32
        array of 4 indexes per bond (ref, i, j, j) for GL_LINES_ADJACENCY
        (ref == j when there is no reference atom). """
    pairs = np.asarray(pairs, dtype=np.int64).reshape(-1, 2)
    mode, ref = multiple_bond_references(pairs, orders, all_pairs, is_hydrogen, max_ring)
    singles = np.ascontiguousarray(pairs[mode < 0].ravel(), dtype=np.uint32)
    groups = {}
    for m in (LINES_MULTI_RING_DOUBLE, LINES_MULTI_DOUBLE, LINES_MULTI_TRIPLE):
        sel = np.flatnonzero(mode == m)
        if sel.size:
            r = np.where(ref[sel] >= 0, ref[sel], pairs[sel, 1])
            groups[m] = np.ascontiguousarray(np.stack([r, pairs[sel, 0], pairs[sel, 1], pairs[sel, 1]], axis=1).ravel(),
                                             dtype=np.uint32)
    return singles, groups


STICKS_TRIPLE_TRIANGLE = 3     # w of u_bond_side_tbo: triple bond as a triangle fixed to the molecule


def multiple_bond_sides(pairs, mode, ref, coords):
    """ Per-bond float32 (n, 4) array for the sticks geometry shader
        (u_bond_side_tbo), xyz in object space:
          double bonds: unit vector perpendicular to the bond in the plane
              bond/ref, w = mode (LINES_MULTI_RING_DOUBLE / LINES_MULTI_DOUBLE);
          triple bonds: a unit vector perpendicular to the bond (cross with the
              Cartesian axis least aligned with it), w = STICKS_TRIPLE_TRIANGLE:
              the three cylinders sit on a triangle that turns with the
              molecule, like the plane of the double bonds;
          anything else (single, no reference atom, reference collinear with
              the bond): w = 2, offset on the screen. """
    pairs = np.asarray(pairs, dtype=np.int64).reshape(-1, 2)
    out = np.zeros((pairs.shape[0], 4), dtype=np.float32)
    out[:, 3] = LINES_MULTI_TRIPLE
    xyz = np.asarray(coords, dtype=np.float64)
    tri = np.flatnonzero(mode == LINES_MULTI_TRIPLE)
    if tri.size:
        axis = xyz[pairs[tri, 1]] - xyz[pairs[tri, 0]]
        axis /= np.maximum(np.linalg.norm(axis, axis=1), 1e-9)[:, None]
        helper = np.eye(3)[np.argmin(np.abs(axis), axis=1)]
        perp = np.cross(axis, helper)
        norm = np.linalg.norm(perp, axis=1)
        ok = norm > 1e-6
        out[tri[ok], :3] = perp[ok] / norm[ok][:, None]
        out[tri[ok], 3] = STICKS_TRIPLE_TRIANGLE
    sel = np.flatnonzero(((mode == LINES_MULTI_RING_DOUBLE) | (mode == LINES_MULTI_DOUBLE)) & (ref >= 0))
    if sel.size == 0:
        return out
    p0 = xyz[pairs[sel, 0]]
    axis = xyz[pairs[sel, 1]] - p0
    axis /= np.maximum(np.linalg.norm(axis, axis=1), 1e-9)[:, None]
    r = xyz[ref[sel]] - p0
    perp = r - np.sum(r * axis, axis=1)[:, None] * axis
    norm = np.linalg.norm(perp, axis=1)
    ok = norm > 1e-3
    out[sel[ok], :3] = perp[ok] / norm[ok][:, None]
    out[sel[ok], 3] = mode[sel[ok]]
    return out


def _hydrogen_test(vm_object):
    """ callable(atom index) -> True for hydrogens of vm_object (reference-atom choice). """
    atoms = getattr(vm_object, "atoms", None) or {}
    def is_hydrogen(i):
        atom = atoms.get(i)
        return atom is not None and str(getattr(atom, "symbol", "")).upper() == "H"
    return is_hydrogen


class LinesRepresentation(Representation):
    """ Class doc """

    def __init__(self, vismol_object, vismol_glcore, indexes, active=True):
        """ Class initialiser """
        super(LinesRepresentation, self).__init__(vismol_object, vismol_glcore, "lines", active, indexes)

        #self.elements = np.uint32(self.indexes.shape[0])
        self.size_vbo  =  None
        self.ratio_vbo =  None
        #print(self.vm_object.cov_radii_array)
        # . double/triple bonds (2026-10-14): ind_vbo holds only the single
        #   bonds; the multiple ones go to multi_vao (program "lines_multi",
        #   GL_LINES_ADJACENCY). The selection pass still draws self.indexes.
        self.multi_vao = None
        self.multi_ind_vbo = None
        self._multi_key = None
        self._n_single_elems = int(self.elements)
        self._multi_ranges = []     # [(mode, byte offset, n elements)]

    def _make_gl_representation_vao_and_vbos(self):
        """ Base VAO plus the VAO of the multiple-bond pass, which shares the
            coordinate and color buffers. """
        super(LinesRepresentation, self)._make_gl_representation_vao_and_vbos()
        program = self.vm_glcore.shader_programs.get("lines_multi")
        if program is None:
            return
        self.multi_vao = self._make_gl_vao()
        self.multi_ind_vbo = self._make_gl_index_buffer(np.zeros(4, dtype=np.uint32))
        for vbo, name in ((self.coord_vbo, "vert_coord"), (self.col_vbo, "vert_color")):
            att = GL.glGetAttribLocation(program, name)
            if att == -1:
                continue
            GL.glBindBuffer(GL.GL_ARRAY_BUFFER, vbo)
            GL.glEnableVertexAttribArray(att)
            GL.glVertexAttribPointer(att, 3, GL.GL_FLOAT, GL.GL_FALSE, 12, ctypes.c_void_p(0))
        GL.glBindVertexArray(0)
        self._multi_key = None

    def _bond_orders(self):
        """ Order of each pair of self.indexes (1 when multiple_bonds is off). """
        n_bonds = int(self.indexes.shape[0] // 2)
        orders = np.ones(n_bonds, dtype=np.int64)
        try:
            if not self.vm_session.vm_config.gl_parameters.get("multiple_bonds", True):
                return orders
        except Exception:
            pass
        get_bond = getattr(self.vm_object, "get_bond", None)
        if get_bond is None or n_bonds == 0:
            return orders
        ib = self.indexes
        for k in range(n_bonds):
            bond = get_bond(int(ib[2 * k]), int(ib[2 * k + 1]))
            if bond is not None:
                orders[k] = int(bond.bond_order)
        return orders

    def _update_multiple_bonds(self):
        """ Re-splits single/multiple bonds when the drawn pairs, the bonds of
            the object or the multiple_bonds flag change (a Builder edit of a
            bond order recreates the representation). Call with self.vao bound. """
        try:
            mb_on = bool(self.vm_session.vm_config.gl_parameters.get("multiple_bonds", True))
        except Exception:
            mb_on = True
        key = (id(self.indexes), int(self.indexes.shape[0]), id(getattr(self.vm_object, "bonds", None)),
               id(getattr(self.vm_object, "bond_order_list", None)), mb_on, self.multi_vao is not None)
        if key == self._multi_key:
            return
        self._multi_key = key
        if self.multi_vao is None or self.indexes.shape[0] == 0:
            self._n_single_elems = int(self.elements)
            self._multi_ranges = []
            return
        orders = self._bond_orders()
        if not (orders >= 2).any():
            singles, groups = self.indexes, {}
        else:
            all_pairs = getattr(self.vm_object, "index_bonds", None)
            if all_pairs is None or len(all_pairs) == 0:
                all_pairs = self.indexes
            singles, groups = multiple_bond_split(self.indexes, orders, all_pairs, _hydrogen_test(self.vm_object))
        singles = np.ascontiguousarray(singles, dtype=np.uint32)
        # . self.vao is bound: this element buffer is the one the first pass draws
        GL.glBindBuffer(GL.GL_ELEMENT_ARRAY_BUFFER, self.ind_vbo)
        GL.glBufferData(GL.GL_ELEMENT_ARRAY_BUFFER, singles.nbytes, singles, GL.GL_DYNAMIC_DRAW)
        self._n_single_elems = int(singles.shape[0])
        self._multi_ranges = []
        if groups:
            flat = np.ascontiguousarray(np.concatenate([groups[m] for m in sorted(groups)]), dtype=np.uint32)
            offset = 0
            for m in sorted(groups):
                n = int(groups[m].shape[0])
                self._multi_ranges.append((m, offset * 4, n))
                offset += n
            GL.glBindVertexArray(self.multi_vao)
            GL.glBindBuffer(GL.GL_ELEMENT_ARRAY_BUFFER, self.multi_ind_vbo)
            GL.glBufferData(GL.GL_ELEMENT_ARRAY_BUFFER, flat.nbytes, flat, GL.GL_DYNAMIC_DRAW)
            GL.glBindVertexArray(self.vao)

    def _draw_multiple_bonds(self, line_width):
        """ Second pass: double/triple bonds (shaders/lines.py, lines_multi). """
        program = self.vm_glcore.shader_programs.get("lines_multi")
        if program is None or self.multi_vao is None or not self._multi_ranges:
            return
        GL.glUseProgram(program)
        GL.glLineWidth(line_width)
        self.vm_glcore.load_matrices(program, self.vm_object.model_mat)
        self.vm_glcore.load_fog(program)
        u_mode = self.vm_glcore._get_uniform_location(program, "u_mode")
        params = self.vm_session.vm_config.gl_parameters
        GL.glUniform1f(self.vm_glcore._get_uniform_location(program, "u_separation"),
                       float(params.get("lines_bond_separation", 0.20)) * float(params.get("lines_bond_spacing", 1.0)))
        GL.glUniform1f(self.vm_glcore._get_uniform_location(program, "u_shrink"), 0.15)
        GL.glBindVertexArray(self.multi_vao)
        for mode, offset, n in self._multi_ranges:
            GL.glUniform1i(u_mode, mode)
            GL.glDrawElements(GL.GL_LINES_ADJACENCY, n, GL.GL_UNSIGNED_INT, ctypes.c_void_p(offset))
        GL.glBindVertexArray(0)

    def draw_representation(self):
        """ Function doc """
        self._check_vao_and_vbos()
        self._enable_anti_alias_to_lines()
        GL.glUseProgram(self.shader_program)
        line_width = self.vm_session.vm_config.gl_parameters["line_width"]
        line_width = (line_width*200/abs(self.vm_glcore.dist_cam_zrp)/2)**0.5
        GL.glLineWidth(line_width)
        self.vm_glcore.load_matrices(self.shader_program, self.vm_object.model_mat)
        self.vm_glcore.load_fog(self.shader_program)
        GL.glBindVertexArray(self.vao)
        #print(self.vm_object.cov_radii_array)
        
        #'''simples bonds  or multiple bonds'''
        #if self.ratio_vbo == None:
        #    try:
        #        self.vm_object._get_covalent_radii()
        #        self.size_vbo, self.ratio_vbo = self._make_gl_impostor_buffer(self.vm_object.cov_radii_array, self.shader_program)
        #    except:
        #        pass
        #        #print('Failed: self.vm_object._get_covalent_radii()')
        #else:
        #    pass
        #    #self.size_vbo, self.ratio_vbo = self._make_gl_impostor_buffer(self.vm_object.cov_radii_array, self.shader_program)



        if self.was_rep_coord_modified:
            self._load_coord_vbo(coord_vbo=True)
            self.was_rep_coord_modified = False
            
            #
        if self.was_rep_ind_modified:
            self._load_ind_vbo(ind_vbo=True)
            self.was_rep_ind_modified = False
            self._multi_key = None      # ind_vbo got every pair again
        if self.was_col_modified:
            self._load_color_vbo(None)
            self.was_col_modified = False

        self._update_multiple_bonds()
        GL.glDrawElements(GL.GL_LINES, self._n_single_elems, GL.GL_UNSIGNED_INT, None)

        GL.glBindVertexArray(0)
        self._draw_multiple_bonds(line_width)
        self._disable_anti_alias_to_lines()
        GL.glLineWidth(1)
        GL.glUseProgram(0)

    def draw_background_sel_representation(self, line_width_factor=5):
        """ Function doc """
        self._check_vao_and_vbos()
        self._disable_anti_alias_to_lines()
        line_width = self.vm_session.vm_config.gl_parameters["line_width_selection"]
        GL.glUseProgram(self.sel_shader_program)
        GL.glLineWidth(line_width) # line_width_factor -> turn the lines thicker
        GL.glEnable(GL.GL_DEPTH_TEST)
        self.vm_glcore.load_matrices(self.sel_shader_program, self.vm_object.model_mat)
        GL.glBindVertexArray(self.sel_vao)
        
        if self.was_sel_coord_modified:
            self._load_coord_vbo(sel_coord_vbo=True)
            self.was_sel_coord_modified = False
        if self.was_sel_ind_modified:
            self._load_ind_vbo(sel_ind_vbo=True)
            self.was_sel_ind_modified = False
        
        GL.glDrawElements(GL.GL_LINES, self.elements, GL.GL_UNSIGNED_INT, None)
        
        GL.glBindVertexArray(0)
        GL.glDisable(GL.GL_DEPTH_TEST)
        GL.glLineWidth(1)
        GL.glUseProgram(0)


class NonBondedRepresentation(Representation):
    """ Class doc """
    
    def __init__ (self, vismol_object, vismol_glcore, indexes, active=True):
        """ Class initialiser """
        super(NonBondedRepresentation, self).__init__(vismol_object, vismol_glcore, "nonbonded", active, indexes)
    
    def draw_representation(self):
        """ Function doc """
        self._check_vao_and_vbos()
        self._enable_anti_alias_to_lines()
        line_width = self.vm_session.vm_config.gl_parameters["line_width"]
        GL.glUseProgram(self.shader_program)
        GL.glLineWidth(line_width*20/abs(self.vm_glcore.dist_cam_zrp))
        self.vm_glcore.load_matrices(self.shader_program, self.vm_object.model_mat)
        self.vm_glcore.load_fog(self.shader_program)
        GL.glBindVertexArray(self.vao)
        
        if self.was_rep_coord_modified:
            self._load_coord_vbo(coord_vbo=True)
            self.was_rep_coord_modified = False
        if self.was_rep_ind_modified:
            self._load_ind_vbo(ind_vbo=True)
            self.was_rep_ind_modified = False
        if self.was_col_modified:
            self._load_color_vbo(None)
            self.was_col_modified = False
        
        GL.glDrawElements(GL.GL_POINTS, self.elements, GL.GL_UNSIGNED_INT, None)
        
        GL.glBindVertexArray(0)
        self._disable_anti_alias_to_lines()
        GL.glLineWidth(1)
        GL.glUseProgram(0)
    
    def draw_background_sel_representation(self, line_width_factor=5):
        """ Function doc """
        self._check_vao_and_vbos()
        self._disable_anti_alias_to_lines()
        GL.glUseProgram(self.sel_shader_program)
        GL.glLineWidth(line_width_factor)
        GL.glEnable(GL.GL_DEPTH_TEST)
        self.vm_glcore.load_matrices(self.sel_shader_program, self.vm_object.model_mat)
        GL.glBindVertexArray(self.sel_vao)
        
        if self.was_sel_coord_modified:
            self._load_coord_vbo(sel_coord_vbo=True)
            self.was_sel_coord_modified = False
        if self.was_sel_ind_modified:
            self._load_ind_vbo(sel_ind_vbo=True)
            self.was_sel_ind_modified = False
        
        GL.glDrawElements(GL.GL_POINTS, self.elements, GL.GL_UNSIGNED_INT, None)
        
        GL.glBindVertexArray(0)
        GL.glDisable(GL.GL_DEPTH_TEST)
        GL.glLineWidth(1)
        GL.glUseProgram(0)


class SticksRepresentation(Representation):
    """ Class doc """
    
    def __init__(self, vismol_object, vismol_glcore, indexes, active=True, is_dynamic = False, name = "sticks", uniform_color = None):
        """ Class initialiser """
        super(SticksRepresentation, self).__init__(vismol_object, vismol_glcore, name, active, indexes, is_dynamic)
        # [EN] 2026-10-15: "sticks" follows gl_parameters['sticks_radius'] at
        # every draw (_sync_radius), so Preferences changes apply on the fly;
        # set_radius() fixes a radius for this representation instead.
        self._radius_fixed = False
        if  name == "sticks":
            self.radius = self.vm_session.vm_config.gl_parameters["sticks_radius"]
            self.spheres = None
        else:
            self.radius = self.vm_session.vm_config.gl_parameters["ribbon_width"]
            #print(set(indexes))
            self.spheres = SpheresRepresentation(vismol_object, vismol_glcore,
                                            active=True, indexes=list(vismol_object.atoms.keys()))

        # [NOVO] Cor unica opcional para as ligacoes (ex.: dynamic bonds em
        # branco). None => comportamento original (cor por atomo). Ver
        # set_uniform_color / _load_color_vbo / _make_gl_representation_vao_and_vbos.
        self.uniform_color = None
        if uniform_color is not None:
            self.set_uniform_color(uniform_color)

    def set_uniform_color(self, rgb):
        """Define uma cor solida para TODAS as ligacoes desta representacao.

        rgb: (R,G,B) ou (R,G,B,A); aceita [0,1] ou [0,255]; alfa e' ignorado.
        Marca uses_uniform_color=True (para recoloracao por atomo pular esta
        rep) e was_col_modified=True (para o buffer subir a GPU no proximo draw).
        """
        rgb = np.array(rgb, dtype=np.float32).ravel()[:3]
        if np.any(rgb > 1.0):
            rgb = rgb / 255.0
        self.uniform_color = rgb
        self.uses_uniform_color = True
        self.was_col_modified = True

    def clear_uniform_color(self):
        """Volta ao comportamento padrao (cor por atomo)."""
        self.uniform_color = None
        self.uses_uniform_color = False
        self.was_col_modified = True

    def _uniform_color_array(self):
        """Array de cor com UMA entrada por atomo do objeto (mesmo racional de
        OneColorDotsRepresentation._set_uniform_color: o desenho e' por indices
        sobre o array de coordenadas completo)."""
        n_atoms = self.vm_object.frames[0].shape[0]
        return np.tile(self.uniform_color, (n_atoms, 1)).astype(np.float32)
            

    def set_radius (self, radius):
        """ Fixed radius for this representation (no longer follows
            gl_parameters['sticks_radius']). """
        self.radius = radius
        self._radius_fixed = True

    def _multiple_bond_geometry(self):
        """ (radius double, radius triple, spacing double, spacing triple) in
            Angstrom from gl_parameters (defaults: the look before these
            parameters existed, 0.4 x and 1.0 x the stick radius). """
        if self.name != "sticks":
            return (0.0, 0.0, 0.0, 0.0)
        gp = self.vm_session.vm_config.gl_parameters
        r = float(self.radius)
        out = []
        for key, default in (("sticks_radius_double", 0.4 * r), ("sticks_radius_triple", 0.4 * r),
                             ("sticks_spacing_double", r), ("sticks_spacing_triple", r)):
            try:
                value = float(gp.get(key, default))
            except (TypeError, ValueError):
                value = default
            out.append(value if value > 0.0 else default)
        return out

    def _sync_radius(self):
        if self.name == "sticks" and not self._radius_fixed:
            try:
                self.radius = float(self.vm_session.vm_config.gl_parameters["sticks_radius"])
            except (KeyError, TypeError, ValueError):
                pass

    def _get_bond_order_per_bond(self):
        """ Retorna a ordem de ligacao PAREADA 1:1 com self.vm_object.index_bonds
            -- ou seja, o k-esimo elemento deste array e a ordem do k-esimo
            PAR (dois elementos) de index_bonds. Esta e a fonte de verdade
            usada pela Buffer Texture (u_bond_order_tbo) lida no geometry
            shader via gl_PrimitiveIDIn (ver shaders/sticks.py).

            [ATUALIZACAO] Substitui _compute_bond_order_per_vertex, que
            expandia isso para um array POR ATOMO tomando a MAIOR ordem entre
            as ligacoes de cada atomo -- correto apenas quando um atomo nao
            participa de ligacoes com ordens diferentes entre si. Agora nao
            ha essa perda de informacao: a ordem certa vai para a ligacao
            certa, nao para o atomo.
        """
        # [NOVO] Dynamic Bonds (regiao QC, is_dynamic=True): a conectividade
        # muda A CADA FRAME (self.vm_object.dynamic_bonds[f]), entao a ordem
        # tambem precisa ser recalculada por frame -- o dict estatico
        # self.vm_object.bonds (usado no branch abaixo) reflete so' a
        # conectividade INICIAL da molecula inteira e nao tem esses pares.
        # get_dynamic_bond_order_for_frame ja' cacheia por frame (so'
        # recalcula quando muda), entao chamar isso a cada draw e' barato.
        if self.is_dynamic:
            frame, f = self.vm_glcore._safe_frame_coords(self.vm_object)
            get_dyn = getattr(self.vm_object, "get_dynamic_bond_order_for_frame", None)
            dyn_pairs = getattr(self.vm_object, "dynamic_bonds", None)
            if get_dyn is not None and dyn_pairs is not None and 0 <= f < len(dyn_pairs):
                n_bonds_dyn = int(np.asarray(dyn_pairs[f]).ravel().shape[0] // 2)
                try:
                    mb_on = self.vm_session.vm_config.gl_parameters.get("multiple_bonds", True)
                except Exception:
                    mb_on = True
                if not mb_on:
                    return np.ones(max(n_bonds_dyn, 1), dtype=np.uint32)
                orders_dyn = get_dyn(f)
                if orders_dyn.shape[0] == n_bonds_dyn and n_bonds_dyn > 0:
                    return orders_dyn.astype(np.uint32)
                return np.ones(max(n_bonds_dyn, 1), dtype=np.uint32)
            # Sem dynamic_bonds populado ainda para este frame: cai no
            # comportamento antigo abaixo (estatico) como fallback seguro.
        # [BUG FIX / ROBUSTEZ] Antes, esta funcao pareava bond_order_list
        # com self.vm_object.index_bonds POR POSICAO. Isso quebra se
        # self.indexes (a copia usada de fato no GL_ELEMENT_ARRAY_BUFFER,
        # ver Representation.__init__ / define_new_indexes_to_vbo) ficar
        # fora de sincronia com self.vm_object.index_bonds -- por exemplo,
        # quando os bonds sao recalculados (nova percepcao de ordem) DEPOIS
        # que esta representacao ja foi criada, e ninguem chama
        # define_new_indexes_to_vbo(...) + was_rep_ind_modified=True para
        # atualizar o buffer da GPU. Nesse caso, self.indexes (o que a GPU
        # de fato desenha) e self.vm_object.index_bonds (usado aqui antes)
        # podem ter tamanhos/ordens diferentes, e o pareamento posicional
        # aplica a ordem ERRADA a ligacao ERRADA (o sintoma visto: um trecho
        # inteiro da cadeia parecendo "ligacao tripla", quando so' uma
        # ligacao no meio deveria ser).
        #
        # Agora buscamos a ordem PELO PAR DE ATOMOS (i,j), direto no dict
        # self.vm_object.bonds (fonte de verdade, chave = par normalizado),
        # usando os MESMOS pares que estao em self.indexes -- exatamente o
        # array que a GPU desenha. Isso e robusto a qualquer descompasso de
        # tamanho/ordem entre self.indexes e vm_object.index_bonds; o pior
        # caso e um par nao encontrado no dict, que cai em ordem 1 (simples).
        ib = np.asarray(self.indexes).ravel()
        n_bonds = int(ib.shape[0] // 2)
        if n_bonds == 0:
            return np.ones(1, dtype=np.uint32)
        # Flag global: se multiple_bonds estiver desligada, todas as ligacoes
        # sao desenhadas como simples (array fica todo 1). Default True se a
        # chave nao existir no config.
        try:
            if not self.vm_session.vm_config.gl_parameters.get("multiple_bonds", True):
                return np.ones(n_bonds, dtype=np.uint32)
        except Exception as _e:
            dprint("[multiple_bonds DEBUG] error reading flag:", _e)
        get_bond = getattr(self.vm_object, "get_bond", None)
        orders = np.ones(n_bonds, dtype=np.uint32)
        if get_bond is None:
            return orders
        n_not_found = 0
        for k in range(n_bonds):
            i = int(ib[2 * k]); j = int(ib[2 * k + 1])
            bond = get_bond(i, j)
            if bond is not None:
                orders[k] = int(bond.bond_order)
            else:
                n_not_found += 1
        # [DEBUG TEMPORARIO] print incondicional (nao usa dprint -- dprint so'
        # mostra algo com EASYHYBRID_DEBUG=1) para diagnosticar se os dados
        # chegam corretos ate aqui (lado Python) ou se o problema e' no
        # shader/GPU. Remova depois de confirmar.
        '''
        print("[DIAG _get_bond_order_per_bond] rep=%r n_bonds=%d nao_encontrados=%d orders=%s"
              % (getattr(self, "name", "?"), n_bonds, n_not_found, orders.tolist()))
        '''
        
        return orders

    def _make_gl_representation_vao_and_vbos(self):
        """ Same as base, plus a per-BOND (per-primitive) bond-order Buffer
            Texture feeding the geometry shader uniform 'u_bond_order_tbo'
            (used to draw double/triple bonds). Lida via gl_PrimitiveIDIn,
            entao nao depende do esquema de atributo por-vertice/por-atomo
            (ver _get_bond_order_per_bond acima e sticks.py). """
        super(SticksRepresentation, self)._make_gl_representation_vao_and_vbos()
        # [NOVO] Se uma cor unica foi definida, substitui o col_vbo (que a base
        # acabou de montar com as cores por atomo) por um buffer preenchido com
        # a cor solida. Reaproveita todo o resto (coord/ind/bond-order).
        if self.uniform_color is not None:
            self.col_vbo = self._make_gl_color_buffer(self._uniform_color_array(), self.shader_program)
        orders = self._get_bond_order_per_bond()
        self.bond_order_tbo_buf, self.bond_order_tex = self._make_gl_bond_order_texture_buffer(orders)
        self._bond_order_tbo_size = orders.shape[0]

    def _refresh_bond_order_tbo(self):
        """ Reenvia bond_order_list para a Buffer Texture. Chamado sempre que
            os indices (index_bonds) sao recarregados (ex.: Dynamic Bonds
            trocando de frame, ou uma nova selecao de atomos), ja que o
            numero/ordem das ligacoes pode ter mudado. Recria a textura do
            zero se o tamanho mudou (uma Buffer Texture nao pode so' 'crescer'
            sem realocar). """
        orders = self._get_bond_order_per_bond()
        self._last_bond_orders = orders
        if getattr(self, "bond_order_tbo_buf", None) is None:
            self.bond_order_tbo_buf, self.bond_order_tex = self._make_gl_bond_order_texture_buffer(orders)
            self._bond_order_tbo_size = orders.shape[0]
            return
        if orders.shape[0] != getattr(self, "_bond_order_tbo_size", -1):
            # Tamanho mudou: recria buffer+textura (mais simples e seguro do
            # que tentar redimensionar uma Buffer Texture existente).
            GL.glDeleteTextures([self.bond_order_tex])
            GL.glDeleteBuffers(1, [self.bond_order_tbo_buf])
            self.bond_order_tbo_buf, self.bond_order_tex = self._make_gl_bond_order_texture_buffer(orders)
            self._bond_order_tbo_size = orders.shape[0]
        else:
            self._update_gl_bond_order_texture_buffer(self.bond_order_tbo_buf, orders)

    def _refresh_bond_side_tbo(self, orders, coords_changed):
        """ [EN] 2026-10-14: Buffer Texture u_bond_side_tbo (RGBA32F, one texel
            per bond, see multiple_bond_sides and shaders/sticks.py): plane of
            each double bond, so that they look like the "lines" ones. The
            reference atoms (rings, substituents) are searched again only when
            the drawn pairs or their orders change; the vectors are recomputed
            when the frame or the coordinates change. Returns True when the
            texture holds valid data for this draw. """
        n_bonds = int(self.indexes.shape[0] // 2)
        if n_bonds == 0 or orders.shape[0] != n_bonds or not (orders >= 2).any():
            return False
        ref_key = (id(self.indexes), n_bonds, orders.tobytes())
        if ref_key != getattr(self, "_side_ref_key", None):
            if self.is_dynamic:
                all_pairs = self.indexes
            else:
                all_pairs = getattr(self.vm_object, "index_bonds", None)
                if all_pairs is None or len(all_pairs) == 0:
                    all_pairs = self.indexes
            self._side_mode, self._side_ref = multiple_bond_references(self.indexes, orders, all_pairs,
                                                                       _hydrogen_test(self.vm_object))
            self._side_ref_key = ref_key
            self._side_coord_key = None
        frame, f = self.vm_glcore._safe_frame_coords(self.vm_object)
        coord_key = (f, id(frame))
        if coords_changed or coord_key != getattr(self, "_side_coord_key", None):
            sides = np.ascontiguousarray(multiple_bond_sides(self.indexes, self._side_mode, self._side_ref, frame))
            if getattr(self, "bond_side_tex", None) is None:
                self.bond_side_tbo_buf = GL.glGenBuffers(1)
                self.bond_side_tex = GL.glGenTextures(1)
                GL.glBindTexture(GL.GL_TEXTURE_BUFFER, self.bond_side_tex)
                GL.glBindBuffer(GL.GL_TEXTURE_BUFFER, self.bond_side_tbo_buf)
                GL.glTexBuffer(GL.GL_TEXTURE_BUFFER, GL.GL_RGBA32F, self.bond_side_tbo_buf)
                GL.glBindTexture(GL.GL_TEXTURE_BUFFER, 0)
            GL.glBindBuffer(GL.GL_TEXTURE_BUFFER, self.bond_side_tbo_buf)
            GL.glBufferData(GL.GL_TEXTURE_BUFFER, sides.nbytes, sides, GL.GL_DYNAMIC_DRAW)
            GL.glBindBuffer(GL.GL_TEXTURE_BUFFER, 0)
            self._side_coord_key = coord_key
        return True

    # [ATUALIZACAO] O override de _make_gl_sel_representation_vao_and_vbos foi
    # removido: o picking/selecao nao precisa mais de nenhum buffer de ordem
    # de ligacao. A ordem agora e lida no geometry shader via
    # u_bond_order_tbo (Buffer Texture) + gl_PrimitiveIDIn; o programa de
    # selecao (sticks_sel) usa a variante sel_geometry_shader_sticks (mesmo
    # geometry shader compilado com PICKING_PASS), que nao declara essa
    # textura e fixa ordem 1 (ligacao simples) -> raio cheio, sem offset.
    # (Antes o sel usava o MESMO shader sem vincular a textura, supondo que
    # texelFetch retornaria 0 -- falso no core profile: um sampler nao
    # vinculado le a textura da unidade 0; no macOS isso gerava o aviso de
    # textura de inteiro sem sinal no picking.) Isso e exatamente o comportamento
    # ja documentado em draw_background_sel_representation ("Independe de
    # bond_order"), entao a classe base (sem bond-order) ja e suficiente.

    def _load_camera_pos(self, program):
        xyz_coords = self.vm_glcore.glcamera.get_modelview_position(self.vm_object.model_mat)
        u_campos = self.vm_glcore._get_uniform_location(program, "u_campos")
        GL.glUniform3fv(u_campos, 1, xyz_coords)
    
    def draw_representation(self):
        """ Function doc """
        self._check_vao_and_vbos ()
        self._sync_radius()
        self._enable_anti_alias_to_lines()
        GL.glEnable(GL.GL_CULL_FACE)
        GL.glCullFace(GL.GL_BACK)
        GL.glUseProgram(self.shader_program)
        self.vm_glcore.load_matrices(self.shader_program, self.vm_object.model_mat)
        self.vm_glcore.load_fog(self.shader_program)
        self.vm_glcore.load_lights(self.shader_program)
        self._load_camera_pos(self.shader_program)
        GL.glBindVertexArray(self.vao)
        
        #radius = self.vm_session.vm_config.gl_parameters["sticks_radius"]
        custom_int_location = self.vm_glcore._get_uniform_location(self.shader_program, "vert_rad")
        GL.glUniform1f(custom_int_location, self.radius)  # Set the integer value to  _size
        
        
        coords_changed = self.was_rep_coord_modified
        if self.was_rep_coord_modified:
            self._load_coord_vbo(coord_vbo=True)
            self.was_rep_coord_modified = False
        if self.was_rep_ind_modified:
            self._load_ind_vbo(ind_vbo=True)
            self.was_rep_ind_modified = False
        if self.was_col_modified:
            # [NOVO] com cor unica ativa, recarrega com a cor solida
            if self.uniform_color is not None:
                self._load_color_vbo(self._uniform_color_array())
            else:
                self._load_color_vbo(None)
            self.was_col_modified = False
        
        # Reenvia bond_order_list para a Buffer Texture sempre que desenha.
        # Cobre tanto o caso estatico (ordem mudou por edicao manual/Builder,
        # sem necessariamente marcar was_rep_ind_modified) quanto o dinamico
        # (Dynamic Bonds recalcula index_bonds/bond_order_list a cada frame
        # de trajetoria). Custo e uma unica glBufferData pequena (uma
        # ligacao = 4 bytes); _refresh_bond_order_tbo so recria buffer/textura
        # do zero se o NUMERO de ligacoes mudou, senao so' atualiza o conteudo.
        self._refresh_bond_order_tbo()
        has_side = self._refresh_bond_side_tbo(self._last_bond_orders, coords_changed)
        u_has_side_loc = self.vm_glcore._get_uniform_location(self.shader_program, "u_has_side")
        u_side_tbo_loc = self.vm_glcore._get_uniform_location(self.shader_program, "u_bond_side_tbo")
        has_side = has_side and u_has_side_loc != -1 and u_side_tbo_loc != -1
        if u_has_side_loc != -1:
            GL.glUniform1i(u_has_side_loc, 1 if has_side else 0)
        if has_side:
            GL.glActiveTexture(GL.GL_TEXTURE2)
            GL.glBindTexture(GL.GL_TEXTURE_BUFFER, self.bond_side_tex)
            GL.glUniform1i(u_side_tbo_loc, 2)
            GL.glActiveTexture(GL.GL_TEXTURE0)
        # [EN] 2026-10-15: radius and spacing of double/triple bonds (Angstrom,
        # Preferences > Sticks; shared by sticks and dynamic bonds, read at
        # every draw so changes apply on the fly). Other users of this class
        # (ribbons) send 0 = the shader's fallback.
        for uniform, value in zip(("u_rad_double", "u_rad_triple", "u_sep_double", "u_sep_triple"),
                                  self._multiple_bond_geometry()):
            loc = self.vm_glcore._get_uniform_location(self.shader_program, uniform)
            if loc != -1:
                GL.glUniform1f(loc, value)

        # Ligacoes multiplas: 3 passadas. A passada 0 desenha o cilindro
        # central/primeiro de TODAS as ligacoes; as passadas 1 e 2 desenham os
        # cilindros laterais apenas onde a ordem (>=2, >=3) exige -- o geometry
        # shader descarta as que nao se aplicam. Cada passada emite no maximo
        # um cilindro por ligacao, entao max_vertices=40 e o programa sempre
        # linka (nao estoura GL_MAX_GEOMETRY_TOTAL_OUTPUT_COMPONENTS).
        u_pass_loc = self.vm_glcore._get_uniform_location(self.shader_program, "u_pass")
        u_sep_loc = self.vm_glcore._get_uniform_location(self.shader_program, "u_separation")
        if u_sep_loc != -1:
            # Separacao entre cilindros de uma ligacao multipla. Reduzida para
            # acompanhar os tubos mais finos (eff_rad*0.6 no shader). Aumente
            # para afastar, diminua para aproximar.
            # [EN] fallback only: the spacings of double/triple bonds come
            # from u_sep_double / u_sep_triple (_multiple_bond_geometry).
            GL.glUniform1f(u_sep_loc, float(self.radius))
        
        # Vincula a Buffer Texture de ordem de ligacao numa unidade de textura
        # dedicada (GL_TEXTURE1, para nao colidir com a unidade 0 usada pelas
        # fontes/labels em outro lugar do codigo) e aponta o sampler do
        # geometry shader (u_bond_order_tbo) pra essa unidade.
        u_bond_tbo_loc = self.vm_glcore._get_uniform_location(self.shader_program, "u_bond_order_tbo")
        # [DIAG TEMPORARIO]
        '''
        if not getattr(self, "_diag_uniform_printed", False):
            print("[DIAG uniforms] rep=%r u_pass_loc=%d u_sep_loc=%d u_bond_tbo_loc=%d bond_order_tex=%r"
                  % (getattr(self, "name", "?"), u_pass_loc, u_sep_loc, u_bond_tbo_loc,
                     getattr(self, "bond_order_tex", None)))
            self._diag_uniform_printed = True
        
        '''
        if u_bond_tbo_loc != -1 and getattr(self, "bond_order_tex", None) is not None:
            GL.glActiveTexture(GL.GL_TEXTURE1)
            GL.glBindTexture(GL.GL_TEXTURE_BUFFER, self.bond_order_tex)
            GL.glUniform1i(u_bond_tbo_loc, 1)
        
        # [BUG FIX] Antes usava len(self.vm_object.index_bonds), que e o
        # array a nivel de OBJETO (pode ter sido recalculado -- ex.: nova
        # percepcao de ordem de ligacao -- com um NUMERO DIFERENTE de
        # ligacoes). O que realmente esta no GL_ELEMENT_ARRAY_BUFFER (e,
        # portanto, o que gl_PrimitiveIDIn de fato percorre) e self.indexes
        # (uma COPIA feita em __init__ / define_new_indexes_to_vbo). Se
        # os dois tamanhos divergirem, pedir mais elementos do que existem
        # no buffer faz o driver ler dados obsoletos/fora dos limites --
        # exatamente o tipo de desalinhamento que faz uma ligacao simples
        # aparecer com o offset/ordem de outra. Ver tambem
        # _get_bond_order_per_bond, que agora usa a MESMA fonte (self.indexes)
        # para montar a Buffer Texture, garantindo que os dois fiquem sempre
        # alinhados um com o outro.
        if self.is_dynamic:
            n_elem = int(len(self.indexes) + 4)
        else:
            n_elem = int(len(self.indexes))
        
        for _pass in range(3):
            if u_pass_loc != -1:
                GL.glUniform1i(u_pass_loc, _pass)
            GL.glDrawElements(GL.GL_LINES, n_elem, GL.GL_UNSIGNED_INT, None)
            if u_pass_loc == -1:
                break  # shader sem suporte a passadas: desenha so uma vez
        
        # Desvincula a Buffer Texture e volta pra unidade 0 (default usado em
        # outros pontos do codigo, ex.: fontes/labels) para nao deixar estado
        # de textura "vazando" para o proximo draw call de outra representacao.
        if u_bond_tbo_loc != -1 and getattr(self, "bond_order_tex", None) is not None:
            GL.glActiveTexture(GL.GL_TEXTURE1)
            GL.glBindTexture(GL.GL_TEXTURE_BUFFER, 0)
            GL.glActiveTexture(GL.GL_TEXTURE0)
        if has_side:
            GL.glActiveTexture(GL.GL_TEXTURE2)
            GL.glBindTexture(GL.GL_TEXTURE_BUFFER, 0)
            GL.glActiveTexture(GL.GL_TEXTURE0)
        
        GL.glBindVertexArray(0)
        self._disable_anti_alias_to_lines()
        GL.glDisable(GL.GL_CULL_FACE)
        GL.glUseProgram(0)
        GL.glLineWidth(1)
        #print(self.spheres)
        #if self.spheres:
        #    self.spheres.draw_representation()
        #    print(self.spheres.sphere_indexes)
            
    def draw_background_sel_representation(self):
        """ Function doc """
        self._check_vao_and_vbos()
        self._disable_anti_alias_to_lines()
        GL.glUseProgram(self.sel_shader_program)
        GL.glEnable(GL.GL_DEPTH_TEST)
        self.vm_glcore.load_matrices(self.sel_shader_program, self.vm_object.model_mat)
        GL.glBindVertexArray(self.sel_vao)
        
        self._sync_radius()
        custom_int_location = self.vm_glcore._get_uniform_location(self.sel_shader_program, "vert_rad")
        GL.glUniform1f(custom_int_location, self.radius)
        # Picking: um cilindro por ligacao (passada 0). Independe de bond_order,
        # entao o sel_vao nao precisa do VBO de ordem.
        u_pass_loc = self.vm_glcore._get_uniform_location(self.sel_shader_program, "u_pass")
        if u_pass_loc != -1:
            GL.glUniform1i(u_pass_loc, 0)
        u_sep_loc = self.vm_glcore._get_uniform_location(self.sel_shader_program, "u_separation")
        if u_sep_loc != -1:
            GL.glUniform1f(u_sep_loc, 0.0)
        
        
        if self.was_sel_coord_modified:
            self._load_coord_vbo(sel_coord_vbo=True)
            self.was_sel_coord_modified = False
        if self.was_sel_ind_modified:
            self._load_ind_vbo(sel_ind_vbo=True)
            self.was_sel_ind_modified = False
        
        GL.glDrawElements(GL.GL_LINES, self.elements, GL.GL_UNSIGNED_INT, None)
        
        GL.glBindVertexArray(0)
        GL.glDisable(GL.GL_DEPTH_TEST)
        GL.glLineWidth(1)
        GL.glUseProgram(0)


class SpheresRepresentation(Representation):
    """ Class doc """
    #            self, vismol_object, vismol_glcore, indexes, active=True, is_dynamic = False, name = "sticks"
    def __init__(self, vismol_object, vismol_glcore, indexes, active=True, vdw = False, mode = 0):
        """ Class initialiser """
        self.mode = mode
        if self.mode == 0 or self.mode == 2 or self.mode == 3:
            super(SpheresRepresentation, self).__init__(vismol_object, vismol_glcore, "spheres", active, indexes)
        elif self.mode ==1:
            super(SpheresRepresentation, self).__init__(vismol_object, vismol_glcore, "picking_spheres", active, indexes)
        
        
        else:
            pass
        import vismol.utils.sphere_data as sphd
        
        if self.mode == 2:
            #for ribbon spheres
            self.level = 1
            self.scale = 1
            self.rad   = 1.98
        
        if self.mode == 3:
            #for ribbon spheres
            self.level = 1
            self.scale = 1
            self.rad   = 0.16
        
        
        else:
            self.level = self.vm_session.vm_config.gl_parameters["sphere_quality"]
            self.scale = self.vm_session.vm_config.gl_parameters["sphere_scale"]
        
        if vdw:
            self.sphere_vertices = sphd.sphere_vertices[self.level] 
        else:
            self.sphere_vertices = sphd.sphere_vertices[self.level]* self.scale
        self.sphere_indexes = sphd.sphere_triangles[self.level]
        self.instances_elemns = self.sphere_indexes.shape[0]
    
    def _make_gl_representation_vao_and_vbos(self):
        """ Function doc """
        logger.debug("building '{}' representation VAO and VBOs".format(self.name))
        self.vao = self._make_gl_vao()
        self.ind_vbo = self._make_gl_index_buffer(self.sphere_indexes)
        self.coord_vbo = self._make_gl_coord_buffer(self.sphere_vertices, self.shader_program)
        self.col_vbo = self._make_gl_color_buffer(np.zeros(3, dtype=np.float32), self.shader_program, instances=True)
        self.rad_vbo = self._make_gl_radius_buffer(np.zeros(1, dtype=np.float32), self.shader_program, instances=True)
        self.insta_vbo = self._make_gl_instance_buffer(np.zeros(3, dtype=np.float32), self.shader_program)
    
    def _make_gl_sel_representation_vao_and_vbos(self):
        """ Function doc """
        logger.debug("building '{}' background selection VAO and VBOs".format(self.name))
        self.sel_vao = self._make_gl_vao()
        self.sel_ind_vbo = self._make_gl_index_buffer(self.sphere_indexes)
        self.sel_coord_vbo = self._make_gl_coord_buffer(self.sphere_vertices, self.sel_shader_program)
        self.sel_col_vbo = self._make_gl_color_buffer(np.zeros(3, dtype=np.float32), self.sel_shader_program, instances=True)
        self.sel_rad_vbo = self._make_gl_radius_buffer(np.zeros(1, dtype=np.float32), self.shader_program, instances=True)
        self.sel_insta_vbo = self._make_gl_instance_buffer(np.zeros(3, dtype=np.float32), self.shader_program)
    
    def _coords_colors_rads(self):
        """ Builds the per-instance coords/colors/rads arrays for the
            non-selection (normal view) spheres VBOs.

            coords change every frame during animation/trajectory playback
            or interactive geometry edits (dihedral rotation, dragging,
            etc). colors and rads only change when the selection (which
            atoms are shown) or the color scheme changes - they do NOT
            depend on the current frame. Splitting them avoids rebuilding
            colors/rads (and the underlying Python-level loop) on every
            single animation frame.
        """
        frame, f = self.vm_glcore._safe_frame_coords(self.vm_object)
        coords = np.ascontiguousarray(frame[self.indexes], dtype=np.float32)
        return coords

    def _colors_rads(self):
        """ Builds the per-instance colors/rads arrays. Vectorized via
            numpy fancy indexing into the object-level color array instead
            of a per-atom Python loop. ball_rad has no object-level array
            (it can be overridden per-atom, e.g. for picking highlights),
            so it still needs a comprehension, but this only runs when the
            atom selection/coloring changes - not on every frame.
        """
        if self.mode == 2 or self.mode == 3:
            colors = self.vm_object.colors[self.indexes]
            rads = np.full(len(self.indexes), self.rad, dtype=np.float32)
        else:
            colors = self.vm_object.colors[self.indexes]
            rads = np.fromiter(
                (self.vm_object.atoms[i].ball_rad for i in self.indexes),
                dtype=np.float32, count=len(self.indexes))
        return colors, rads

    def _sel_coords_colors_rads(self):
        """ Same split as _coords_colors_rads/_colors_rads but for the
            background-selection (picking) VBOs, which use color_indexes
            instead of colors.
        """
        frame, f = self.vm_glcore._safe_frame_coords(self.vm_object)
        coords = np.ascontiguousarray(frame[self.indexes], dtype=np.float32)
        return coords

    def _sel_colors_rads(self):
        colors = self.vm_object.color_indexes[self.indexes]
        rads = np.fromiter(
            (self.vm_object.atoms[i].ball_rad for i in self.indexes),
            dtype=np.float32, count=len(self.indexes))
        return colors, rads
    
    def draw_representation(self):
        """ Function doc """
        self._check_vao_and_vbos()
        GL.glEnable(GL.GL_DEPTH_TEST)
        GL.glEnable(GL.GL_CULL_FACE)
        GL.glCullFace(GL.GL_BACK)
        GL.glBlendFunc(GL.GL_SRC_ALPHA, GL.GL_ONE_MINUS_SRC_ALPHA)
        GL.glEnable(GL.GL_BLEND)
        GL.glUseProgram(self.shader_program)
        self.vm_glcore.load_matrices(self.shader_program, self.vm_object.model_mat)
        self.vm_glcore.load_lights(self.shader_program)
        self.vm_glcore.load_fog(self.shader_program)
        GL.glBindVertexArray(self.vao)
        
        # [EN] 2026-10-15: the joint spheres of the sticks (mode 3) follow
        # gl_parameters['sticks_radius'] (Preferences, on the fly).
        if self.mode == 3:
            try:
                r = float(self.vm_session.vm_config.gl_parameters.get("sticks_radius", self.rad))
            except (TypeError, ValueError):
                r = self.rad
            if r != self.rad:
                self.rad = r
                self.was_col_modified = True
        
        # colors/rads only depend on which atoms are shown and their color
        # scheme, not on the current frame, so they're only rebuilt when
        # the index/selection or coloring actually changed.
        if self.was_rep_ind_modified or self.was_col_modified:
            colors, rads = self._colors_rads()
            GL.glBindBuffer(GL.GL_ARRAY_BUFFER, self.col_vbo)
            GL.glBufferData(GL.GL_ARRAY_BUFFER, colors.nbytes, colors, GL.GL_STATIC_DRAW)
            GL.glBindBuffer(GL.GL_ARRAY_BUFFER, self.rad_vbo)
            GL.glBufferData(GL.GL_ARRAY_BUFFER, rads.nbytes, rads, GL.GL_STATIC_DRAW)
            self.elements = np.uint32(len(self.indexes))
            self.was_col_modified = False
        
        if self.was_rep_coord_modified or self.was_rep_ind_modified:
            coords = self._coords_colors_rads()
            GL.glBindBuffer(GL.GL_ARRAY_BUFFER, self.insta_vbo)
            GL.glBufferData(GL.GL_ARRAY_BUFFER, coords.nbytes, coords, GL.GL_DYNAMIC_DRAW)
            self.elements = np.uint32(coords.shape[0])
            self.was_rep_coord_modified = False
            self.was_rep_ind_modified = False
        
        GL.glDrawElementsInstanced(GL.GL_TRIANGLES, self.instances_elemns, GL.GL_UNSIGNED_INT, None, self.elements)
        
        GL.glBindVertexArray(0)
        GL.glUseProgram(0)
        GL.glDisable(GL.GL_CULL_FACE)
        GL.glDisable(GL.GL_DEPTH_TEST)
        #print('aqui', self.sphere_indexes)

    def draw_background_sel_representation(self):
        """ Function doc """
        
        if self.mode == 2 or self.mode == 3:
            #Selection is not necessary
            return None
        
        self._check_vao_and_vbos()
        GL.glEnable(GL.GL_DEPTH_TEST)
        GL.glUseProgram(self.sel_shader_program)
        self.vm_glcore.load_matrices(self.sel_shader_program, self.vm_object.model_mat)
        GL.glBindVertexArray(self.sel_vao)
        
        if self.was_sel_ind_modified:
            colors, rads = self._sel_colors_rads()
            GL.glBindBuffer(GL.GL_ARRAY_BUFFER, self.sel_col_vbo)
            GL.glBufferData(GL.GL_ARRAY_BUFFER, colors.nbytes, colors, GL.GL_STATIC_DRAW)
            GL.glBindBuffer(GL.GL_ARRAY_BUFFER, self.sel_rad_vbo)
            GL.glBufferData(GL.GL_ARRAY_BUFFER, rads.nbytes, rads, GL.GL_STATIC_DRAW)
            self.elements = np.uint32(len(self.indexes))
        
        if self.was_sel_coord_modified or self.was_sel_ind_modified:
            coords = self._sel_coords_colors_rads()
            GL.glBindBuffer(GL.GL_ARRAY_BUFFER, self.sel_insta_vbo)
            GL.glBufferData(GL.GL_ARRAY_BUFFER, coords.nbytes, coords, GL.GL_DYNAMIC_DRAW)
            self.elements = np.uint32(coords.shape[0])
            self.was_sel_coord_modified = False
            self.was_sel_ind_modified = False
        
        GL.glDrawElementsInstanced(GL.GL_TRIANGLES, self.instances_elemns, GL.GL_UNSIGNED_INT, None, self.elements)
        
        GL.glBindVertexArray(0)
        GL.glDisable(GL.GL_DEPTH_TEST)
        GL.glUseProgram(0)


class DashedLinesRepresentation(Representation):
    """ Class doc """
    
    def __init__(self, vismol_object, vismol_glcore, indexes, active=True, depth_test = False):
        """ Class initialiser """
        super(DashedLinesRepresentation, self).__init__(vismol_object, vismol_glcore, "dash", active, indexes)
        self.depth_test = depth_test
        self.color2 = [1.0 ,1.0, 0.0]
    
    def draw_representation(self):
        """ Function doc """
        self._check_vao_and_vbos()
        self._enable_anti_alias_to_lines()
        GL.glUseProgram(self.shader_program)
        line_width = self.vm_session.vm_config.gl_parameters["line_width"]
        line_width = (line_width*200/abs(self.vm_glcore.dist_cam_zrp)/2)**0.5
        GL.glLineWidth(line_width)
        self.vm_glcore.load_matrices(self.shader_program, self.vm_object.model_mat)
        self.vm_glcore.load_fog(self.shader_program)
        GL.glBindVertexArray(self.vao)
        
        # How to pass data to shader on the fly. It's not the most 
        # efficient mode, but it's an okay solution in this case.
        color = self.vm_glcore._get_uniform_location(self.shader_program, "uniform_color")
        #GL.glUniformMatrix4fv(proj, 1, GL.GL_FALSE, self.glcamera.projection_matrix)
        color2 = np.array(self.color2, dtype=np.float32)
        GL.glUniform3fv(color, 1, color2)
        
        '''
        #How to pass data to shader on the fly. It's not the most efficient wood, but it's an okay solution in this case.
        test = GL.glGetUniformLocation(self.shader_program, "test_int")
        test_int = 5
        GL.glUniform1i(test, test_int)
        '''
        
        if self.depth_test:
            pass
        else:
            GL.glDisable(GL.GL_DEPTH_TEST)
        
        
        
        if self.was_rep_coord_modified:
            self._load_coord_vbo(coord_vbo=True)
            self.was_rep_coord_modified = False
        if self.was_rep_ind_modified:
            self._load_ind_vbo(ind_vbo=True)
            self.was_rep_ind_modified = False
        if self.was_col_modified:
            self._load_color_vbo(None)
            self.was_col_modified = False
        
        GL.glDrawElements(GL.GL_LINES, self.elements, GL.GL_UNSIGNED_INT, None)
        
        GL.glBindVertexArray(0)
        self._disable_anti_alias_to_lines()
        GL.glLineWidth(1)
        GL.glUseProgram(0)
        
    def draw_background_sel_representation(self):
        """ Function doc """
        self._check_vao_and_vbos()
        self._disable_anti_alias_to_lines()
        line_width = self.vm_session.vm_config.gl_parameters["line_width_selection"]
        GL.glUseProgram(self.sel_shader_program)
        GL.glLineWidth(line_width)
        GL.glEnable(GL.GL_DEPTH_TEST)
        self.vm_glcore.load_matrices(self.sel_shader_program, self.vm_object.model_mat)
        GL.glBindVertexArray(self.sel_vao)
        
        if self.was_sel_coord_modified:
            self._load_coord_vbo(sel_coord_vbo=True)
            self.was_sel_coord_modified = False
        if self.was_sel_ind_modified:
            self._load_ind_vbo(sel_ind_vbo=True)
            self.was_sel_ind_modified = False
        
        GL.glDrawElements(GL.GL_LINES, self.elements, GL.GL_UNSIGNED_INT, None)
        
        GL.glBindVertexArray(0)
        GL.glDisable(GL.GL_DEPTH_TEST)
        GL.glLineWidth(1)
        GL.glUseProgram(0)


class ImpostorRepresentation(Representation):
    """ Class doc """
    
    def __init__ (self, vismol_object, vismol_glcore, indexes, active=True):
        """ Class initialiser """
        super(ImpostorRepresentation, self).__init__(vismol_object, vismol_glcore, "impostor", active, indexes)
    
    def _make_gl_representation_vao_and_vbos(self):
        """ Function doc """
        logger.debug("building '{}' representation VAO and VBOs".format(self.name))
        self.vao = self._make_gl_vao()
        self.ind_vbo = self._make_gl_index_buffer(self.indexes)
        self.coord_vbo = self._make_gl_coord_buffer(self.vm_object.frames[0], self.shader_program)
        self.col_vbo = self._make_gl_color_buffer(self.vm_object.colors, self.shader_program)
        self.size_vbo, self.ratio_vbo = self._make_gl_impostor_buffer(self.vm_object.cov_radii_array, self.shader_program)
    
    def _make_gl_sel_representation_vao_and_vbos(self):
        """ Function doc """
        logger.debug("building '{}' background selection VAO and VBOs".format(self.name))
        self.sel_vao = self._make_gl_vao()
        self.sel_ind_vbo = self._make_gl_index_buffer(self.indexes)
        self.sel_coord_vbo = self._make_gl_coord_buffer(self.vm_object.frames[0], self.sel_shader_program)
        self.sel_col_vbo = self._make_gl_color_buffer(self.vm_object.color_indexes, self.sel_shader_program)
        self.sel_size_vbo, self.sel_ratio_vbo = self._make_gl_impostor_buffer(self.vm_object.cov_radii_array, self.shader_program)
    
    # def _modified_window_size(self):
    #     ratio = self.vm_glcore.width / self.vm_glcore.height
    #     if self.ratio != ratio:
    #         return True
    #     return False
    
    def _load_camera_pos(self, program):
        xyz_coords = self.vm_glcore.glcamera.get_modelview_position(self.vm_object.model_mat)
        u_campos = self.vm_glcore._get_uniform_location(program, "u_campos")
        GL.glUniform3fv(u_campos, 1, xyz_coords)
    
    def draw_representation(self):
        """ Function doc """
        self._check_vao_and_vbos()
        # if self._modified_window_size():
        #     self._load_impostor_ratio_vbo(coord_vbo=True)
        self._enable_anti_alias_to_lines()
        GL.glUseProgram(self.shader_program)
        GL.glEnable(GL.GL_VERTEX_PROGRAM_POINT_SIZE)
        self.vm_glcore.load_matrices(self.shader_program, self.vm_object.model_mat)
        self.vm_glcore.load_fog(self.shader_program)
        self.vm_glcore.load_lights(self.shader_program)
        self._load_camera_pos(self.shader_program)
        GL.glBindVertexArray(self.vao)
        # pm = self.vm_glcore.glcamera.projection_matrix
        # print(pm[0,0] * pm[1,1], "fov")
        
        if self.was_rep_coord_modified:
            self._load_coord_vbo(coord_vbo=True)
            self.was_rep_coord_modified = False
        if self.was_rep_ind_modified:
            self._load_ind_vbo(ind_vbo=True)
            self.was_rep_ind_modified = False
        if self.was_col_modified:
            self._load_color_vbo(None)
            self.was_col_modified = False
        
        GL.glDrawElements(GL.GL_POINTS, self.elements, GL.GL_UNSIGNED_INT, None)
        
        GL.glBindVertexArray(0)
        self._disable_anti_alias_to_lines()
        GL.glPointSize(1)
        GL.glUseProgram(0)
    
    def draw_background_sel_representation(self):
        """ Function doc """
        self._check_vao_and_vbos()
        # if self._modified_window_size():
        #     self._load_impostor_ratio_vbo(sel_coord_vbo=True)
        GL.glUseProgram(self.sel_shader_program)
        GL.glEnable(GL.GL_VERTEX_PROGRAM_POINT_SIZE)
        self.vm_glcore.load_matrices(self.sel_shader_program, self.vm_object.model_mat)
        # self.vm_glcore.load_lights(self.shader_program)
        self._load_camera_pos(self.sel_shader_program)
        GL.glBindVertexArray(self.sel_vao)
        
        if self.was_sel_coord_modified:
            self._load_coord_vbo(sel_coord_vbo=True)
            self.was_sel_coord_modified = False
        if self.was_sel_ind_modified:
            self._load_ind_vbo(sel_ind_vbo=True)
            self.was_sel_ind_modified = False
        
        GL.glDrawElements(GL.GL_POINTS, self.elements, GL.GL_UNSIGNED_INT, None)
        
        GL.glBindVertexArray(0)
        GL.glDisable(GL.GL_DEPTH_TEST)
        GL.glPointSize(1)
        GL.glUseProgram(0)


class CellLineRepresentation:
    """ Class doc 
    
    This is a modified version of the Representations class. 
    The differences are in the methods:

    # _make_gl_representation_vao_and_vbos
    # _load_coord_vbo
    # _load_color_vbo

    In this case, there is a modification in the code so that 
    the coordinates, indices and colors transmitted to the GPU 
    refer to the cell parameters.

    In this class there is no selection presentation.
    """
    
    def __init__ (self, vismol_object, vismol_glcore, name, active, indexes, is_dynamic = False):
        self.vm_object = vismol_object
        self.vm_session = vismol_object.vm_session
        self.vm_glcore = vismol_glcore
        self.name = name
        self.active = active
        self.indexes = np.array(indexes, dtype=np.uint32)
        self.elements = np.uint32(self.indexes.shape[0])
        
        self.is_dynamic = is_dynamic
        
        self.was_rep_modified = False
        self.was_sel_modified = False
        self.was_col_modified = False
        self.was_rep_coord_modified = False
        self.was_sel_coord_modified = False
        self.was_rep_ind_modified = False
        self.was_sel_ind_modified = False
        self.was_rep_col_modified = False
        # representation
        self.vao = None
        self.ind_vbo = None
        self.coord_vbo = None
        self.col_vbo = None
        self.size_vbo = None
        # selection
        #self.sel_vao = None
        #self.sel_ind_vbo = None
        #self.sel_coord_vbo = None
        #self.sel_col_vbo = None
        #self.sel_size_vbo = None
        # shaders
        self.shader_program = None
        self.sel_shader_program = None
    
    def _check_vao_and_vbos(self):
        #print(self.name)
        self.shader_program = self.vm_glcore.shader_programs[self.name]
        #self.sel_shader_program = self.vm_glcore.shader_programs[self.name + "_sel"]
        if self.vao is None:
            self._make_gl_representation_vao_and_vbos()
        #if self.sel_vao is None:
        #    self._make_gl_sel_representation_vao_and_vbos()
    
    def _make_gl_representation_vao_and_vbos(self):
        """ Function doc """
        logger.debug("building '{}' representation VAO and VBOs".format(self.name))
        self.vao = self._make_gl_vao()
        self.ind_vbo   = self._make_gl_index_buffer(self.indexes)
        self.coord_vbo = self._make_gl_coord_buffer(self.vm_object.cell_coordinates[0], self.shader_program)
        self.col_vbo   = self._make_gl_color_buffer(self.vm_object.cell_colors, self.shader_program)
    
    #def _make_gl_sel_representation_vao_and_vbos(self):
    #    """ Function doc """
    #    logger.debug("building '{}' background selection VAO and VBOs".format(self.name))
    #    self.sel_vao = self._make_gl_vao()
    #    self.sel_ind_vbo = self._make_gl_index_buffer(self.indexes)
    #    self.sel_coord_vbo = self._make_gl_coord_buffer(self.vm_object.cell_coordinates[0], self.sel_shader_program)
    #    self.sel_col_vbo = self._make_gl_color_buffer(self.vm_object.color_indexes, self.sel_shader_program)
    
    def _make_gl_vao(self):
        """ Function doc """
        vao = GL.glGenVertexArrays(1)
        GL.glBindVertexArray(vao)
        return vao
    
    def _make_gl_index_buffer(self, indexes):
        """ Function doc """
        ind_vbo = GL.glGenBuffers(1)
        GL.glBindBuffer(GL.GL_ELEMENT_ARRAY_BUFFER, ind_vbo)
        GL.glBufferData(GL.GL_ELEMENT_ARRAY_BUFFER, indexes.nbytes, indexes, GL.GL_DYNAMIC_DRAW)
        return ind_vbo
    
    def _make_gl_coord_buffer(self, coords, program):
        """ Function doc """
        coord_vbo = GL.glGenBuffers(1)
        GL.glBindBuffer(GL.GL_ARRAY_BUFFER, coord_vbo)
        GL.glBufferData(GL.GL_ARRAY_BUFFER, coords.nbytes, coords, GL.GL_STATIC_DRAW)
        att_position = GL.glGetAttribLocation(program, "vert_coord")
        GL.glEnableVertexAttribArray(att_position)
        GL.glVertexAttribPointer(att_position, 3, GL.GL_FLOAT, GL.GL_FALSE, 3*coords.itemsize, ctypes.c_void_p(0))
        return coord_vbo
    
    def _make_gl_color_buffer(self, colors, program, instances=False):
        """ Function doc """
        col_vbo = GL.glGenBuffers(1)
        GL.glBindBuffer(GL.GL_ARRAY_BUFFER, col_vbo)
        GL.glBufferData(GL.GL_ARRAY_BUFFER, colors.nbytes, colors, GL.GL_STATIC_DRAW)
        att_colors = GL.glGetAttribLocation(program, "vert_color")
        GL.glEnableVertexAttribArray(att_colors)
        GL.glVertexAttribPointer(att_colors, 3, GL.GL_FLOAT, GL.GL_FALSE, 3*colors.itemsize, ctypes.c_void_p(0))
        if instances:
            GL.glVertexAttribDivisor(att_colors, 1)
        return col_vbo
    
    #def _make_gl_radius_buffer(self, radii, program, instances=False):
    #    """ Function doc """
    #    rad_vbo = GL.glGenBuffers(1)
    #    GL.glBindBuffer(GL.GL_ARRAY_BUFFER, rad_vbo)
    #    GL.glBufferData(GL.GL_ARRAY_BUFFER, radii.nbytes, radii, GL.GL_STATIC_DRAW)
    #    att_rads = GL.glGetAttribLocation(program, "vert_radius")
    #    GL.glEnableVertexAttribArray(att_rads)
    #    GL.glVertexAttribPointer(att_rads, 1, GL.GL_FLOAT, GL.GL_FALSE, radii.itemsize, ctypes.c_void_p(0))
    #    if instances:
    #        GL.glVertexAttribDivisor(att_rads, 1)
    #    return rad_vbo
    
    #def _make_gl_instance_buffer(self, instances, program):
    #    """ Function doc """
    #    insta_vbo = GL.glGenBuffers(1)
    #    GL.glBindBuffer(GL.GL_ARRAY_BUFFER, insta_vbo)
    #    GL.glBufferData(GL.GL_ARRAY_BUFFER, instances.nbytes, instances, GL.GL_STATIC_DRAW)
    #    gl_insta = GL.glGetAttribLocation(program, "vert_instance")
    #    GL.glEnableVertexAttribArray(gl_insta)
    #    GL.glVertexAttribPointer(gl_insta, 3, GL.GL_FLOAT, GL.GL_FALSE, 0, ctypes.c_void_p(0))
    #    GL.glVertexAttribDivisor(gl_insta, 1)
    #    return insta_vbo
    
    #def _make_gl_impostor_buffer(self, impostors_radii, program):
    #    """ Function doc """
    #    size_vbo = GL.glGenBuffers(1)
    #    GL.glBindBuffer(GL.GL_ARRAY_BUFFER, size_vbo)
    #    GL.glBufferData(GL.GL_ARRAY_BUFFER, impostors_radii.nbytes, impostors_radii, GL.GL_STATIC_DRAW)
    #    att_size = GL.glGetAttribLocation(program, "vert_dot_size")
    #    GL.glEnableVertexAttribArray(att_size)
    #    GL.glVertexAttribPointer(att_size, 1, GL.GL_FLOAT, GL.GL_FALSE, impostors_radii.itemsize, ctypes.c_void_p(0))
    #    
    #    self.ratio = self.vm_glcore.width / self.vm_glcore.height
    #    ratio_vbo = 1
    #    # ratio = np.repeat(self.ratio, impostors_radii.shape[0])
    #    # ratio_vbo = GL.glGenBuffers(1)
    #    # GL.glBindBuffer(GL.GL_ARRAY_BUFFER, ratio_vbo)
    #    # GL.glBufferData(GL.GL_ARRAY_BUFFER, ratio.nbytes, ratio, GL.GL_STATIC_DRAW)
    #    # att_ratio = GL.glGetAttribLocation(program, "hw_ratio")
    #    # GL.glEnableVertexAttribArray(att_ratio)
    #    # GL.glVertexAttribPointer(att_ratio, 1, GL.GL_FLOAT, GL.GL_FALSE, ratio.itemsize, ctypes.c_void_p(0))
    #    return size_vbo, ratio_vbo
    
    def _load_coord_vbo(self, coord_vbo=False, sel_coord_vbo=False):
        """ This function assigns the coordinates to 
        be drawn by the function  draw_representation"""
        _, f = self.vm_glcore._safe_frame_coords(self.vm_object)
        #frame = self.vm_object.cell_coordinates[0]

        # . "f" is clamped against the atom trajectory length
        #   (vm_object.frames), but the cell box may have fewer stored
        #   frames than that -- e.g. a single static cell shown for a
        #   multi-frame trajectory, or a cell recorded for only some frames
        #   of an NPT run. Clamp separately here so we never index past
        #   what cell_coordinates actually has (falls back to the last
        #   known cell frame instead of crashing).
        number_of_cell_frames = self.vm_object.cell_coordinates.shape[0]
        if f >= number_of_cell_frames:
            f = number_of_cell_frames - 1

        frame = self.vm_object.cell_coordinates[f]
        if coord_vbo:
            GL.glBindBuffer(GL.GL_ARRAY_BUFFER, self.coord_vbo)
            GL.glBufferData(GL.GL_ARRAY_BUFFER, frame.nbytes, frame, GL.GL_STATIC_DRAW)
        
        #if sel_coord_vbo:
        #    GL.glBindBuffer(GL.GL_ARRAY_BUFFER, self.sel_coord_vbo)
        #    GL.glBufferData(GL.GL_ARRAY_BUFFER, frame.nbytes, frame, GL.GL_STATIC_DRAW)
    
    def _load_ind_vbo(self, ind_vbo=False, sel_ind_vbo=False):
        """ Function doc """
        #if self.is_dynamic:
        #    frame, f = self.vm_glcore._safe_frame_coords(self.vm_object)
        #    self.define_new_indexes_to_vbo(input_indexes = self.vm_object.dynamic_bonds[f])
            
        if ind_vbo:
            GL.glBindBuffer(GL.GL_ELEMENT_ARRAY_BUFFER, self.ind_vbo)
            GL.glBufferData(GL.GL_ELEMENT_ARRAY_BUFFER, self.indexes.nbytes, self.indexes, GL.GL_DYNAMIC_DRAW)
        
        #if sel_ind_vbo:
        #    GL.glBindBuffer(GL.GL_ELEMENT_ARRAY_BUFFER, self.sel_ind_vbo)
        #    GL.glBufferData(GL.GL_ELEMENT_ARRAY_BUFFER, self.indexes.nbytes, self.indexes, GL.GL_DYNAMIC_DRAW)
    
    def _load_color_vbo(self, colors):
        """ This function assigns the colors to
            be drawn by the function  draw_representation"""
        GL.glBindBuffer(GL.GL_ARRAY_BUFFER, self.col_vbo)
        GL.glBufferData(GL.GL_ARRAY_BUFFER, self.vm_object.cell_colors.nbytes, self.vm_object.cell_colors, GL.GL_STATIC_DRAW)
    
    def _enable_anti_alias_to_lines(self):
        """ Function doc """
        GL.glEnable(GL.GL_DEPTH_TEST)
        GL.glEnable(GL.GL_BLEND)
        GL.glBlendFunc(GL.GL_SRC_ALPHA, GL.GL_ONE_MINUS_SRC_ALPHA)
        GL.glEnable(GL.GL_LINE_SMOOTH)
        GL.glHint(GL.GL_LINE_SMOOTH_HINT, GL.GL_NICEST)
    
    def _disable_anti_alias_to_lines(self):
        """ Function doc """
        GL.glDisable(GL.GL_LINE_SMOOTH)
        GL.glDisable(GL.GL_BLEND)
        GL.glDisable(GL.GL_DEPTH_TEST)
    
    def define_new_indexes_to_vbo(self, input_indexes):
        """ Function doc """
        self.indexes = np.array(input_indexes, dtype=np.uint32)
        self.elements = np.uint32(self.indexes.shape[0])


    def draw_representation(self):
        """ Function doc """
        self._check_vao_and_vbos()
        self._enable_anti_alias_to_lines()
        GL.glUseProgram(self.shader_program)
        line_width = self.vm_session.vm_config.gl_parameters["line_width"]
        line_width = (line_width*200/abs(self.vm_glcore.dist_cam_zrp)/2)**0.5
        GL.glLineWidth(line_width)
        self.vm_glcore.load_matrices(self.shader_program, self.vm_object.model_mat)
        self.vm_glcore.load_fog(self.shader_program)
        GL.glBindVertexArray(self.vao)
        
        if self.was_rep_coord_modified:
            self._load_coord_vbo(coord_vbo=True)
            self.was_rep_coord_modified = False
        if self.was_rep_ind_modified:
            self._load_ind_vbo(ind_vbo=True)
            self.was_rep_ind_modified = False
        #self.was_col_modified = True
        if self.was_col_modified:
            self._load_color_vbo(None)
            self.was_col_modified = False
        
        GL.glDrawElements(GL.GL_LINES, self.elements, GL.GL_UNSIGNED_INT, None)
        
        GL.glBindVertexArray(0)
        self._disable_anti_alias_to_lines()
        GL.glLineWidth(1)
        GL.glUseProgram(0)

    def draw_background_sel_representation(self, line_width_factor=5):
        """ Function doc """
        pass


class LabelRepresentation:
    
    def __init__ (self, vismol_object, vismol_glcore, indexes, labels, color = [1, 1, 1, 1]):
        self.vm_object = vismol_object
        self.vm_session = vismol_object.vm_session
        self.vm_glcore = vismol_glcore
        # Font family/size for THIS representation (atom index, MM charge,
        # residue name/index, chain, ...) are customizable independently
        # from the picking/distance labels, via the Preferences window
        # ("Atom Labels (glArea)"), and persisted in gl_parameters under
        # 'atom_label_font_file'/'atom_label_font_size'.
        gp = self.vm_session.vm_config.gl_parameters
        _font_file = gp.get("atom_label_font_file", DEFAULT_FONT_FILE)
        _font_size = gp.get("atom_label_font_size", DEFAULT_FONT_SIZE)
        self.vm_font = VismolFont(font_file=resolve_font_path(_font_file),
                                   char_width=_font_size, char_height=_font_size,
                                   color=color)
        # See VismolFont.zoom_sensitivity: 0.0 (default) keeps atom
        # labels a constant size on screen regardless of zoom; 1.0
        # restores the old behavior of shrinking/growing with camera
        # distance like the rest of the 3D scene. Shared by ALL glArea
        # labels via the single "Labels scale with zoom" checkbox on the
        # Viewer > General tab in Preferences (picking and distance
        # labels read the same 'labels_zoom_sensitivity' key -- see
        # VismolGLCore.__init__).
        self.vm_font.zoom_sensitivity = gp.get("labels_zoom_sensitivity", 1.0)
        self.indexes = indexes
        self.labels = labels
        self.active = True
        self.was_rep_ind_modified = False
        self.was_sel_ind_modified = False    
        self.is_dynamic = False
    def define_new_indexes_to_vbo(self, indexes):
        """ Function doc """
        self.indexes = indexes
    
    def draw_representation(self):
        """ Function doc """
        if self.vm_glcore.dragging:
            return False
        
        entries = []
        for index in self.indexes:
            atom = self.vm_object.atoms[index]
            # Falls back to the atom name if no explicit label_text was
            # set (e.g. via the "Show atom index/charge/residue.../" menu
            # in eSession.py), and to an empty string as a last resort so
            # this never raises on a None text.
            text = atom.label_text or atom.name or ''
            if not text:
                continue
            frame = self.vm_glcore._get_vismol_object_frame(atom.vm_object)
            point = atom.coords(frame)
            # x_shift = -len(text)/2.0 centers the label on the atom
            # instead of anchoring its first character there. Baked
            # per-entry (not into a shared string_shift uniform)
            # because this single draw call mixes many atoms with
            # DIFFERENT label lengths -- see _draw_text_labels().
            entries.append((text, point, -len(text) / 2.0))
        
        self.vm_glcore._draw_text_labels(self.vm_font, entries)


    def draw_background_sel_representation(self, line_width_factor=5):
        """ Function doc """
        pass


class SurfaceRepresentation(Representation):
    """ Class doc """
    

    def __init__(self, vismol_object, vismol_glcore, name, active, indexes, is_dynamic = False, iso_color = [0,1,1], surface_name = None):
    #def __init__(self, vismol_object = None, vismol_glcore = None, indexes = None, active=True, vdw = False, mode = 0):
        super(SurfaceRepresentation, self).__init__(vismol_object, vismol_glcore, name, active , indexes, is_dynamic)


        self.active    = True
        self.vm_glcore = vismol_glcore
        self.name      =  name   
        self.vismol_object = vismol_object
        self.surf_name     = surface_name
        # "surface" (malha preenchida, GL_FILL) ou "lines" (wireframe, GL_LINE).
        # Ambos usam o MESMO shader de triangulos/VAO -- so muda como o
        # rasterizador desenha os triangulos (glPolygonMode), entao nao ha
        # necessidade de manter um shader/VAO separado so para o wireframe.
        self.render_mode = "surface"
        self.alpha = 1.0   # 1.0 = opaco, 0.0 = totalmente transparente
        self.smooth_shading = False   # False = flat (normal por face), True = smooth (normal por vertice)

    def set_render_mode(self, mode):
        """ mode: 'surface' (preenchido) ou 'lines' (wireframe) """
        if mode not in ("surface", "lines"):
            raise ValueError("render_mode must be 'surface' or 'lines', got: {}".format(mode))
        self.render_mode = mode

    def set_alpha(self, alpha):
        """ alpha: 0.0 (totalmente transparente) a 1.0 (opaco) """
        self.alpha = max ( 0.0, min ( 1.0, float ( alpha ) ) )

    def set_shading_mode(self, mode):
        """ mode: 'flat' (normal constante por face -- mostra as facetas
        do marching cubes) ou 'smooth' (normal interpolada por vertice,
        media das faces adjacentes -- ver compute_smooth_normals em
        surface_analysis_window.py). """
        if mode not in ("flat", "smooth"):
            raise ValueError("shading mode must be 'flat' or 'smooth', got: {}".format(mode))
        self.smooth_shading = ( mode == "smooth" )

    def _make_gl_representation_vao_and_vbos(self, debug = False):
        """ Function doc """
        if debug:
            logger.debug("building '{}' representation VAO and VBOs".format(self.name))
        self.vao       = self._make_gl_vao()
        self.ind_vbo   = self._make_gl_index_buffer(self.indexes)
        self.coord_vbo = self._make_gl_coord_buffer(self.vertices, self.shader_program)
        # instances=True estava ligado aqui, mas draw_representation() usa
        # glDrawElements normal (nao instanciado) -- com o divisor de
        # atributo ligado, a GPU sempre lia colors[0] pra malha inteira,
        # nao a cor de cada vertice. Passava despercebido porque orbital/
        # potential/density usam uma cor uniforme (todo vertice igual);
        # quebrou visivelmente com o MEP, que tem cor diferente por vertice.
        self.col_vbo   = self._make_gl_color_buffer(self.colors, self.shader_program)
        self.norm_vbo  = self._make_gl_normal_buffer(self.normals, self.shader_program)

    def _make_gl_normal_buffer(self, normals, program):
        """ Buffer de normais por vertice, pro atributo 'vert_normal' que o
        vertex_shader_surface ja declarava mas nunca tinha dado real
        associado (normal ficava sempre (0,0,0) por padrao do OpenGL). """
        norm_vbo = GL.glGenBuffers(1)
        GL.glBindBuffer(GL.GL_ARRAY_BUFFER, norm_vbo)
        GL.glBufferData(GL.GL_ARRAY_BUFFER, normals.nbytes, normals, GL.GL_STATIC_DRAW)
        att_normal = GL.glGetAttribLocation(program, "vert_normal")
        GL.glEnableVertexAttribArray(att_normal)
        GL.glVertexAttribPointer(att_normal, 3, GL.GL_FLOAT, GL.GL_FALSE, 3*normals.itemsize, ctypes.c_void_p(0))
        return norm_vbo

    def _load_normal_vbo(self, normals):
        """ Reenvia as normais a cada frame, igual _load_color_vbo faz
        com as cores (a malha pode mudar de frame pra frame numa
        trajetoria). """
        GL.glBindBuffer(GL.GL_ARRAY_BUFFER, self.norm_vbo)
        GL.glBufferData(GL.GL_ARRAY_BUFFER, normals.nbytes, normals, GL.GL_STATIC_DRAW)
    
    def _load_coord_vbo(self, coord_vbo=False, vertices = None):
        """ This function assigns the coordinates to 
        be drawn by the function  draw_representation"""
        frame, f = self.vm_glcore._safe_frame_coords(self.vm_object)
        if coord_vbo:
            GL.glBindBuffer(GL.GL_ARRAY_BUFFER, self.coord_vbo)
            GL.glBufferData(GL.GL_ARRAY_BUFFER, vertices.nbytes, vertices, GL.GL_STATIC_DRAW)

    def _load_ind_vbo(self, ind_vbo=False, indexes = None ):
        """ Function doc """
 
        if ind_vbo:
            GL.glBindBuffer(GL.GL_ELEMENT_ARRAY_BUFFER, self.ind_vbo)
            GL.glBufferData(GL.GL_ELEMENT_ARRAY_BUFFER, indexes.nbytes, indexes, GL.GL_DYNAMIC_DRAW)
 
    def _load_color_vbo(self, colors):
        """ This function assigns the colors to
            be drawn by the function  draw_representation"""
        GL.glBindBuffer(GL.GL_ARRAY_BUFFER, self.col_vbo)
        GL.glBufferData(GL.GL_ARRAY_BUFFER, colors.nbytes, colors, GL.GL_STATIC_DRAW)

    def draw_background_sel_representation(self):
        pass
    
    def draw_representation(self):
        """ Function doc """
        frame_coords, frame = self.vm_glcore._safe_frame_coords(self.vismol_object)
        self.vertices, self.colors, self.indexes, self.normals = self.vismol_object.surface_trajectory[frame][self.surf_name]
        
        self._check_vao_and_vbos()
        GL.glEnable(GL.GL_DEPTH_TEST)
        #GL.glEnable(GL.GL_CULL_FACE)
        #GL.glCullFace(GL.GL_BACK)
        GL.glBlendFunc(GL.GL_SRC_ALPHA, GL.GL_ONE_MINUS_SRC_ALPHA)
        GL.glEnable(GL.GL_BLEND)
        GL.glUseProgram(self.shader_program)
        alpha_loc = GL.glGetUniformLocation(self.shader_program, "surf_alpha")
        GL.glUniform1f(alpha_loc, self.alpha)
        smooth_loc = GL.glGetUniformLocation(self.shader_program, "smooth_shading")
        GL.glUniform1i(smooth_loc, 1 if self.smooth_shading else 0)
        if self.alpha < 1.0:
            # com transparencia, desliga escrita no depth buffer (mas mantem
            # o teste de depth ligado) -- evita que a propria superficie se
            # auto-oclua de forma esquisita (frente vs fundo da mesma malha),
            # sem afetar como ela interage com objetos opacos na cena.
            GL.glDepthMask(GL.GL_FALSE)
        self.vm_glcore.load_matrices(self.shader_program, self.vm_object.model_mat)
        self.vm_glcore.load_lights(self.shader_program)
        self.vm_glcore.load_fog(self.shader_program)
        GL.glBindVertexArray(self.vao)
        
        #self.vismol_object
        

        self._load_coord_vbo( coord_vbo=True, vertices = self.vertices)
        self._load_ind_vbo  ( ind_vbo = True, indexes = self.indexes)
        self._load_color_vbo( self.colors)
        self._load_normal_vbo( self.normals)
        
        #if self.was_rep_coord_modified or self.was_rep_ind_modified:
        #    coords, colors, rads = self._coords_colors_rads()
        #    GL.glBindBuffer(GL.GL_ARRAY_BUFFER, self.insta_vbo)
        #    GL.glBufferData(GL.GL_ARRAY_BUFFER, coords.nbytes, coords, GL.GL_STATIC_DRAW)
        #    GL.glBindBuffer(GL.GL_ARRAY_BUFFER, self.col_vbo)
        #    GL.glBufferData(GL.GL_ARRAY_BUFFER, colors.nbytes, colors, GL.GL_STATIC_DRAW)
        #    GL.glBindBuffer(GL.GL_ARRAY_BUFFER, self.rad_vbo)
        #    GL.glBufferData(GL.GL_ARRAY_BUFFER, rads.nbytes, rads, GL.GL_STATIC_DRAW)
        #    self.elements = np.uint32(coords.shape[0])
        #    self.was_rep_coord_modified = False
        #    self.was_rep_ind_modified = False
        
        #GL.glDrawElementsInstanced(GL.GL_TRIANGLES, self.instances_elemns, GL.GL_UNSIGNED_INT, None, self.elements)
        if self.render_mode == "lines":
            GL.glPolygonMode(GL.GL_FRONT_AND_BACK, GL.GL_LINE)
        else:
            GL.glPolygonMode(GL.GL_FRONT_AND_BACK, GL.GL_FILL)
        GL.glDrawElements(GL.GL_TRIANGLES, int(len(self.indexes)), GL.GL_UNSIGNED_INT, None)
        # sempre restaura GL_FILL: glPolygonMode e estado global do OpenGL,
        # se nao resetar aqui, a proxima representacao desenhada (sticks,
        # cartoon, etc) tambem sairia em wireframe sem querer.
        GL.glPolygonMode(GL.GL_FRONT_AND_BACK, GL.GL_FILL)
        if self.alpha < 1.0:
            GL.glDepthMask(GL.GL_TRUE)   # restaura -- estado global do OpenGL
        GL.glBindVertexArray(0)
        GL.glUseProgram(0)
        GL.glDisable(GL.GL_CULL_FACE)
        GL.glDisable(GL.GL_DEPTH_TEST)

    def _check_vao_and_vbos(self):
        #print(self.name)
        self.shader_program = self.vm_glcore.shader_programs[self.name]
        #self.sel_shader_program = self.vm_glcore.shader_programs[self.name + "_sel"]
        if self.vao is None:
            self._make_gl_representation_vao_and_vbos()
        #if self.sel_vao is None:
        #    self._make_gl_sel_representation_vao_and_vbos()

def _cartoon_settings_from_gl_parameters ( gl_parameters ):
    """ [EN] Builds a ribbon_geometry-compatible settings dict (see
    DEFAULT_CARTOON_SETTINGS in vismol.utils.ribbon_geometry) from this
    session's own gl_parameters -- the "Preferences > Cartoon" tab (user
    request: adjustable ribbon thickness + a "variable criterion" for
    helix/strand assignment; see setup_interface.py/.glade's
    'cartoon_*' widgets and __apply_cartoon_parameters()) writes its
    values there under the 'cartoon_*' keys, read here with .get(...,
    default) against ribbon_geometry's OWN defaults -- so a
    .config.json that predates this feature (missing every 'cartoon_*'
    key) falls back to exactly the same look the app already shipped
    with, no migration needed. The 2 '*_tolerance' keys are stored as a
    PERCENTAGE in gl_parameters (100 = published P-SEA values, matching
    how the Preferences entry field reads "100%"), converted to the
    0-1-scale-factor form ribbon_geometry's classify_secondary_structure()
    expects. 'cartoon_coil_radius' is a single UI value applied to BOTH
    of ribbon_geometry's coil_half_width/coil_half_height (coil is
    rendered round, so only one thickness concept makes sense there). """
    from vismol.utils.ribbon_geometry import DEFAULT_CARTOON_SETTINGS
    d = DEFAULT_CARTOON_SETTINGS
    coil_radius = gl_parameters.get ( 'cartoon_coil_radius', d['coil_half_width'] )
    return {
        'helix_half_width':      gl_parameters.get ( 'cartoon_helix_width',   d['helix_half_width'] ),
        'helix_half_height':     gl_parameters.get ( 'cartoon_helix_height',  d['helix_half_height'] ),
        'strand_half_width':     gl_parameters.get ( 'cartoon_strand_width',  d['strand_half_width'] ),
        'strand_half_height':    gl_parameters.get ( 'cartoon_strand_height', d['strand_half_height'] ),
        'arrow_half_width':      gl_parameters.get ( 'cartoon_arrow_width',   d['arrow_half_width'] ),
        'arrow_tip_half_width':  d['arrow_tip_half_width'],
        'coil_half_width':       coil_radius,
        'coil_half_height':      coil_radius,
        'color_helix':  tuple ( gl_parameters.get ( 'cartoon_color_helix',  d['color_helix'] ) ),
        'color_strand': tuple ( gl_parameters.get ( 'cartoon_color_strand', d['color_strand'] ) ),
        'color_coil':   tuple ( gl_parameters.get ( 'cartoon_color_coil',   d['color_coil'] ) ),
        'helix_tolerance_scale':  gl_parameters.get ( 'cartoon_helix_tolerance',  d['helix_tolerance_scale']  * 100.0 ) / 100.0,
        'strand_tolerance_scale': gl_parameters.get ( 'cartoon_strand_tolerance', d['strand_tolerance_scale'] * 100.0 ) / 100.0,
        # [EN] 2026-09-26 user request ("um parametro minimo de residuos
        # para ser considerado uma fita e uma helice... faz sentido um
        # unico residuo?") -- replaced the single shared 'min_run' with
        # per-SS-type minimums; see ribbon_backbone.py's own
        # _smooth_short_runs() docstring for why the defaults differ.
        'min_run_helix':  int ( gl_parameters.get ( 'cartoon_min_run_helix',  d['min_run_helix']  ) ),
        'min_run_strand': int ( gl_parameters.get ( 'cartoon_min_run_strand', d['min_run_strand'] ) ),
        # [EN] 2026-09-26 aesthetic fix -- see ribbon_geometry.py's own
        # DEFAULT_CARTOON_SETTINGS comment for the full "shark fin" story.
        'arrow_taper_residues': int ( gl_parameters.get ( 'cartoon_arrow_taper_residues', d['arrow_taper_residues'] ) ),
        # [EN] 2026-09-26 aesthetic fix ("sinuosa" strand centerline vs a
        # PyMOL-style straighter "flecha") -- see ribbon_geometry.py's own
        # _smooth_strand_guide_points() docstring.
        'strand_smoothing_iterations': int ( gl_parameters.get ( 'cartoon_strand_smoothing', d['strand_smoothing_iterations'] ) ),
    }


class CartoonRepresentation(Representation):
    """ [EN] Independent Cartoon (ribbon) representation -- user's own
    explicit request: a from-scratch implementation, not influenced by
    the one this replaces, JMOL-like aesthetic (flat ribbons for helix/
    strand with an arrowhead on strands, a thin tube through coil/loop).

    ALL of the actual spline/secondary-structure/mesh-building work lives
    in vismol.utils.ribbon_backbone (chain-aware backbone extraction + a
    fresh P-SEA-style Ca-only H/S/C classifier) and vismol.utils.
    ribbon_geometry (Catmull-Rom spline, per-residue orientation frames
    with a flip-consistency pass, per-SS cross-section extrusion into a
    triangle mesh) -- this class owns ONLY the GL side (VAO/VBO upload +
    draw), same "fully generated geometry, nothing to do with vm_object's
    own per-atom frames/colors arrays" situation SurfaceRepresentation
    already has (its own _make_gl_representation_vao_and_vbos() override,
    just above in this same file, is the pattern this mirrors).

    Deliberately does NOT call Representation.__init__() (its positional
    `indexes` argument assumes atom/bond-derived indexes, which doesn't
    apply to a wholly generated mesh) -- sets the same attribute set by
    hand instead, matching this file's own established convention for
    "fully custom geometry" representations.

    Live trajectory/frame rebuild (user's own later follow-up request --
    this used to be a documented v1 limitation, mesh built ONCE at
    construction time): draw_representation() below now rebuilds the
    whole mesh from scratch whenever was_rep_coord_modified is set (the
    same dirty flag every frame-advance already sets on every active
    representation, see vm_glcore.render()'s own "updated_coords"
    block) -- see rebuild() and draw_representation()'s own comments for
    why a full rebuild, not a cheap coordinate-only update, is what a new
    frame actually needs here. """

    def __init__ ( self, vismol_object, vismol_glcore, name = 'cartoon', active = True, show_all = True ):
        self.name       = name
        self.active     = active
        self.vm_object  = vismol_object
        self.vm_session = vismol_object.vm_session
        self.vm_glcore  = vismol_glcore

        self.uses_uniform_color = True
        self.is_dynamic = False

        self.was_rep_modified       = False
        self.was_sel_modified       = False
        self.was_col_modified       = False
        self.was_rep_coord_modified = False
        self.was_sel_coord_modified = False
        self.was_rep_ind_modified   = False
        self.was_sel_ind_modified   = False
        self.was_rep_col_modified   = False
        self._last_uploaded_sel_ids = None

        # representation
        self.vao            = None
        self.ind_vbo        = None
        self.coord_vbo      = None
        self.norm_vbo       = None
        self.col_vbo        = None
        self.size_vbo       = None
           

        # bgrd selection   
        self.sel_vao        = None
        self.sel_ind_vbo    = None
        self.sel_coord_vbo  = None
        self.sel_col_vbo    = None
        self.sel_size_vbo   = None


        #     S H A D E R S
        self.shader_program     = None
        self.sel_shader_program = None

        # [EN] 2026-09-26 user request ("alterar a transparencia da
        # representacao de cartoon") -- same role/name/range as
        # SurfaceRepresentation's own self.alpha/set_alpha() just above in
        # this file; see draw_representation()'s own blend/depth-mask
        # handling below for how this is actually applied. Read from
        # gl_parameters at construction time (like every other cartoon_*
        # style setting) so a NEW cartoon rep created after the user has
        # already set a transparency preference starts with it applied,
        # not fully opaque until the next Preferences Apply.
        _gp_alpha = getattr ( vismol_object.vm_session.vm_config, 'gl_parameters', {} )
        self.alpha = 1.0 - float ( _gp_alpha.get ( 'cartoon_transparency', 0 ) ) / 100.0

        # [EN] User request (round 8 -- "Preferences: cartoon dinamico ou
        # nao, por padrao nao"): per-chain H/S/C classification cache, see
        # rebuild()'s own docstring below for the full story. None here
        # means "nothing cached yet" -- the very first build (right below)
        # always classifies fresh, same as before this cache existed.
        self._cached_ss_codes = None

        # [EN] 2026-09-26 user request ("tornar a representacao clicavel...
        # plot o c-alpha trace como referencia") -- separate picking-only
        # geometry (a thick Ca-Ca line trace, coloured by each residue's
        # own real picking-ID colour), built alongside the visible mesh
        # below and again on every rebuild(); see
        # _build_pick_trace_geometry()'s own docstring for the full story.
        self.pick_coords = None
        self.pick_colors = None

        # [EN] 2026-09-28 user request ("o usuario pode escolher o padrao
        # de cores para o cartoon (cada representacao pode ter um padrao
        # diferente)... colorir tipo rainbow"): the colour scheme belongs
        # to THIS representation instance, not to gl_parameters -- see
        # set_color_mode() and ribbon_geometry.residue_colors_for_records()
        # for the available modes. The SS colours themselves ("ss" mode)
        # still come from Preferences > Cartoon, shared by everyone.
        self.color_mode    = 'ss'
        self.uniform_color = ( 0.85, 0.85, 0.85 )

        # [EN] 2026-09-28 user request ("representar apenas parte da
        # estrutura na forma de cartoon, uma vez selecionada e atribuida
        # pelo usuario"): which residues get drawn is the per-atom
        # atom.cartoon flag on each residue's Ca (same per-atom-flag
        # convention as atom.sticks/atom.lines, driven by vm_session.
        # show_or_hide(rep_type='cartoon', ...)). show_all=True (the
        # treeview "Representation > Cartoon" toggle) means "the whole
        # object" -- every atom gets flagged here; show_or_hide() creates
        # the representation with show_all=False instead, AFTER flagging
        # only the selected atoms itself.
        if show_all:
            self.set_all_residues_shown ( True, rebuild = False )

        self._build_mesh ( cached_ss_codes = None )

    def _residue_is_shown ( self, ca_atom ):
        return bool ( getattr ( ca_atom, 'cartoon', False ) )

    def has_shown_residues ( self ):
        """ True if at least one Ca atom of this object is flagged for
        the cartoon (i.e. there is something to draw). """
        for atom in self.vm_object.atoms.values ( ):
            if atom.name == 'CA' and getattr ( atom, 'cartoon', False ):
                return True
        return False

    def set_all_residues_shown ( self, show = True, rebuild = True ):
        """ Flags (or clears) every atom of the object -- the "whole
        object" cartoon used by the treeview's Representation menu. """
        for atom in self.vm_object.atoms.values ( ):
            atom.cartoon = bool ( show )
        if rebuild:
            self.rebuild ( )

    def set_color_mode ( self, color_mode, uniform_color = None ):
        """ Per-representation colour scheme: 'ss' (secondary structure,
        Preferences colours), 'rainbow' (N->C blue->red per chain),
        'chain', 'atom' (each residue's Ca atom colour) or 'uniform'
        (`uniform_color`, RGB 0-1). Rebuilds the mesh -- the colour is
        baked into the vertex buffer, like everything else here -- but
        reuses the cached SS classification (colour never changes it). """
        from vismol.utils.ribbon_geometry import CARTOON_COLOR_MODES
        if color_mode not in CARTOON_COLOR_MODES:
            raise ValueError ( "unknown cartoon color mode: {}".format ( color_mode ) )
        self.color_mode = color_mode
        if uniform_color is not None:
            self.uniform_color = tuple ( float ( c ) for c in list ( uniform_color )[:3] )
        self.rebuild ( )

    def _build_mesh ( self, cached_ss_codes = None ):
        """ Shared by __init__ and rebuild(): runs build_ribbon_mesh_
        for_object() with the current gl_parameters settings + this
        representation's own colour scheme and shown-residue filter, and
        stores the resulting buffers (plus the pick trace). """
        from vismol.utils.ribbon_geometry import build_ribbon_mesh_for_object
        gl_parameters = getattr ( self.vm_session.vm_config, 'gl_parameters', {} )
        settings = _cartoon_settings_from_gl_parameters ( gl_parameters )
        # [EN] 2026-09-26: user-tunable via Preferences > Cartoon ("Ribbon
        # Detail"); default 5 UNCHANGED from before, since a higher value
        # visibly smooths out the ribbon (confirmed live, spline_detail=12
        # on a real protein) but roughly triples full-protein mesh-build
        # cost -- exactly the cost round 9 just finished cutting down for
        # large structures (6VXX-scale). Left low by default so that fix
        # isn't silently undone; a small/medium protein can afford to turn
        # this up for a noticeably smoother look.
        spline_detail = int ( gl_parameters.get ( 'cartoon_spline_detail', 5 ) )
        coords, normals, colors, indexes, ss_codes_per_chain = build_ribbon_mesh_for_object (
            self.vm_object, spline_detail = spline_detail, settings = settings,
            cached_ss_codes = cached_ss_codes, return_ss_codes = True,
            color_mode = self.color_mode, uniform_color = self.uniform_color,
            ca_filter = self._residue_is_shown )
        self._cached_ss_codes = ss_codes_per_chain

        self.coords  = np.ascontiguousarray ( coords,  dtype = np.float32 )
        self.colors  = np.ascontiguousarray ( colors,  dtype = np.float32 )
        self.normals = np.ascontiguousarray ( normals, dtype = np.float32 )
        self.indexes = np.ascontiguousarray ( indexes, dtype = np.uint32 )
        self.elements = np.uint32 ( self.indexes.shape[0] )

        # [EN] 2026-09-26: real atom positions may have moved (trajectory/
        # MD frame change) -- the pick trace must track the SAME current
        # frame the visible mesh above was just rebuilt from.
        self.pick_coords, self.pick_colors = self._build_pick_trace_geometry ( )

    def _build_pick_trace_geometry ( self ):
        """ [EN] 2026-09-26 user request: "vamos tornar a represetacao
        clicavel para a ferramenta de selecao / centralizacao, minha
        sugestao e que no background a gente plot o c-alpha trace como
        referencia para encontrar os residuos." Exactly this: the ribbon
        mesh itself has no 1:1 vertex-to-atom correspondence (it's a
        spline-extruded surface, nowhere near the same vertex count/order
        as vm_object's own atoms), so it can't get real per-vertex picking
        colours the way Lines/Sticks/Spheres do -- but the backbone Ca
        trace, ONE point per residue, genuinely does. Builds a plain list
        of Ca-Ca line SEGMENTS (GL_LINES topology, no index buffer needed
        -- see draw_background_sel_representation()'s glDrawArrays call),
        one segment per consecutive residue pair within each contiguous
        chain record (extract_backbone_chains() already refuses to bridge
        a chain change or a real sequence gap -- reused as-is here, so the
        pick trace never draws a segment where the ribbon itself has none
        either). BOTH endpoints of a given segment get the SAME colour
        (that segment's own first residue's real atom.color_id, the exact
        per-atom picking-ID encoding _pick() already knows how to decode
        via vm_session.atom_dic_id -- see Atom._generate_atom_unique_
        color_id()) -- deliberately, so there is no colour INTERPOLATION
        across a segment's length to corrupt the ID (a blended colour
        would decode to a bogus, nonexistent ID). Uses the RAW (un-
        smoothed) Ca positions, not round 11's smoothed spline centerline
        -- this is a real atom's real position, not a rendering guide. """
        from vismol.utils.ribbon_backbone import extract_backbone_chains
        records = extract_backbone_chains ( self.vm_object )
        coords = [ ]
        colors = [ ]
        for record in records:
            ca = record["ca"]
            ca_atoms = record["ca_atoms"]
            for i in range ( len ( ca ) - 1 ):
                # [EN] 2026-09-28 partial cartoon: only where the ribbon
                # itself is drawn (both ends of the segment shown).
                if not ( self._residue_is_shown ( ca_atoms[i] ) and self._residue_is_shown ( ca_atoms[i + 1] ) ):
                    continue
                color = ca_atoms[i].color_id
                coords.append ( ca[i] );     colors.append ( color )
                coords.append ( ca[i + 1] ); colors.append ( color )
        if not coords:
            return ( np.zeros ( ( 0, 3 ), dtype = np.float32 ),
                      np.zeros ( ( 0, 3 ), dtype = np.float32 ) )
        return ( np.ascontiguousarray ( coords, dtype = np.float32 ),
                  np.ascontiguousarray ( colors, dtype = np.float32 ) )

    def rebuild ( self, force_reclassify = False ):
        """ [EN] Re-runs the mesh build (current gl_parameters settings,
        current frame -- build_ribbon_mesh_for_object()'s own frame=None
        default resolves to "whichever frame is active for THIS object
        right now", so this naturally picks up trajectory playback too)
        and marks the GL buffers stale so the next draw goes through the
        base class's _check_vao_and_vbos() (triggers on `self.vao is
        None`) and rebuilds them from the fresh geometry. Two callers:
        setup_interface.py's __apply_cartoon_parameters()/reset (with
        force_reclassify=True -- a Preferences change, including to the
        classification tolerances themselves, must always be reflected
        immediately, never silently skipped by the cache below); and
        draw_representation() below (force_reclassify=False, the
        default), whenever was_rep_coord_modified says the session's
        current frame just changed (trajectory navigation/MD playback).
        Old VBO/VAO ids are simply replaced (not explicitly glDelete'd)
        -- the same lightweight "just overwrite the attribute" convention
        this codebase already uses for other forced-rebuild cases (e.g.
        picking_dots reset to None in setup_interface.py's
        color_set_viewing_selections()).

        [EN] User request (round 8 -- "cartoon dinamico ou nao, por
        padrao nao, ou seja, o que e helice ou fita nao muda com a troca
        de frames"): classify_secondary_structure() turned out to be
        ~45% of a rebuild's total cost once build_ribbon_mesh() itself
        got vectorized (round 7's own profiling) -- re-running it on
        EVERY trajectory frame is both wasted work (a single, rigid
        protein's SS assignment realistically doesn't change frame to
        frame in the vast majority of real use) and a source of visible
        flicker if it's borderline/noisy across frames. Controlled by
        gl_parameters['cartoon_dynamic_secondary_structure'] (Preferences
        > Cartoon, default OFF): OFF (the default) reuses this
        representation's own cached per-chain classification instead of
        recomputing it (build_ribbon_mesh_for_object()'s own
        `cached_ss_codes` parameter, which additionally guards against a
        length mismatch -- e.g. a genuine topology change mid-trajectory
        -- by falling back to a fresh classification for that one chain
        only); ON always reclassifies, matching round 6/7's original
        behaviour. Either way the cache is refreshed with whatever was
        actually used this call, so toggling the setting later picks up
        cleanly from wherever playback currently is. """
        gl_parameters = getattr ( self.vm_session.vm_config, 'gl_parameters', {} )
        dynamic_ss = bool ( gl_parameters.get ( 'cartoon_dynamic_secondary_structure', False ) )
        cached = None if ( force_reclassify or dynamic_ss ) else self._cached_ss_codes
        self._build_mesh ( cached_ss_codes = cached )
        self.vao     = None
        self.sel_vao = None

    def set_alpha ( self, alpha ):
        """ alpha: 0.0 (totalmente transparente) a 1.0 (opaco) -- same
        contract as SurfaceRepresentation.set_alpha(). """
        self.alpha = max ( 0.0, min ( 1.0, float ( alpha ) ) )

    def _make_gl_representation_vao_and_vbos ( self, debug = False ):
        """ [EN] Overrides the base Representation version (which reads
        self.vm_object.frames[0]/self.vm_object.colors -- the whole
        object's plain per-atom buffers) because Cartoon's geometry, like
        Surface's, is entirely generated (a spline-extruded ribbon mesh,
        built in __init__ above by vismol.utils.ribbon_geometry) and has
        nothing to do with vm_object's own per-atom arrays -- neither the
        vertex COUNT nor the ORDER match. Mirrors SurfaceRepresentation's
        own override of this same method, just above in this file. """
        if debug:
            logger.debug ( "building '{}' representation VAO and VBOs".format ( self.name ) )
        self.vao       = self._make_gl_vao ( )
        self.ind_vbo   = self._make_gl_index_buffer ( self.indexes )
        self.coord_vbo = self._make_gl_coord_buffer ( self.coords, self.shader_program )
        self.col_vbo   = self._make_gl_color_buffer ( self.colors, self.shader_program )
        self.norm_vbo  = self._make_gl_normal_buffer ( self.normals, self.shader_program )

    def _make_gl_normal_buffer ( self, normals, program ):
        """ Per-vertex normal buffer for the shader's own "vert_norm"
        attribute (shaders/cartoon.py) -- the base Representation class
        has no normal-buffer helper of its own (most representation
        types don't need per-vertex normals), so this mirrors the base
        class's own _make_gl_coord_buffer()/_make_gl_color_buffer()
        pattern for the one extra attribute Cartoon's shader needs. """
        norm_vbo = GL.glGenBuffers ( 1 )
        GL.glBindBuffer ( GL.GL_ARRAY_BUFFER, norm_vbo )
        GL.glBufferData ( GL.GL_ARRAY_BUFFER, normals.nbytes, normals, GL.GL_STATIC_DRAW )
        att_normal = GL.glGetAttribLocation ( program, "vert_norm" )
        GL.glEnableVertexAttribArray ( att_normal )
        GL.glVertexAttribPointer ( att_normal, 3, GL.GL_FLOAT, GL.GL_FALSE, 3 * normals.itemsize, ctypes.c_void_p ( 0 ) )
        return norm_vbo

    def _make_gl_sel_representation_vao_and_vbos ( self, debug = False ):
        """ [EN] 2026-09-26: background/picking VAO for the Ca-trace built
        by _build_pick_trace_geometry() -- NOT the visible ribbon mesh
        (self.coords/self.indexes), which has no 1:1 vertex-to-atom
        correspondence to give real per-vertex pick colours to. No index
        buffer: draw_background_sel_representation() below uses
        glDrawArrays(GL_LINES, ...) directly over pick_coords, which is
        already laid out as consecutive (start, end) pairs. """
        if debug:
            logger.debug ( "building '{}' background selection VAO and VBOs".format ( self.name ) )
        self.sel_vao       = self._make_gl_vao ( )
        self.sel_coord_vbo = self._make_gl_coord_buffer ( self.pick_coords, self.sel_shader_program )
        self.sel_col_vbo   = self._make_gl_color_buffer ( self.pick_colors, self.sel_shader_program )

    def draw_representation ( self ):
        # [EN] User request: "extender a representacao para todos os
        # frames de uma trajetoria" -- was previously a documented v1
        # limitation (mesh built ONCE, at construction). vm_glcore.render()
        # already marks was_rep_coord_modified=True on every ACTIVE
        # representation of every object whenever the session's current
        # frame advances (see the "updated_coords" block near the top of
        # render() in vismol_glcore.py -- the same flag every atom-based
        # representation already consumes to re-upload just its coord
        # VBO). Cartoon can't do the cheap "just update positions" thing
        # those do, since a new frame's geometry can change EVERYTHING a
        # rebuild computes from scratch -- the spline shape, each
        # residue's side/orientation vector, even its P-SEA secondary-
        # structure classification (a real MD trajectory can genuinely
        # unfold/refold) -- so this reuses the exact same rebuild() the
        # Preferences "Apply" flow already added: a full re-run of
        # build_ribbon_mesh_for_object(), reusing its own frame=None
        # default, which already resolves to "whichever frame is
        # currently active for THIS object" (Atom.coords()'s own
        # per-object frame-pin/clamping logic -- unchanged, nothing
        # cartoon-specific needed here). rebuild() sets self.vao = None,
        # so the _check_vao_and_vbos() call right below picks it straight
        # back up and re-uploads fresh GL buffers, same lazy-rebuild path
        # every other representation already uses.
        if self.was_rep_coord_modified:
            self.rebuild ( )
            self.was_rep_coord_modified = False
            self.was_sel_coord_modified = False

        self._check_vao_and_vbos ( )

        GL.glEnable ( GL.GL_DEPTH_TEST )
        # [EN] A ribbon is a THIN, open surface (helix/strand cross-
        # sections are only ~0.5 A tall) -- both sides are visible from
        # many camera angles, so face culling stays off; the shader's own
        # gl_FrontFacing-based normal flip (shaders/cartoon.py) is what
        # keeps lighting correct on both sides of the same surface.
        GL.glDisable ( GL.GL_CULL_FACE )

        GL.glUseProgram ( self.shader_program )
        self.vm_glcore.load_matrices ( self.shader_program, self.vm_object.model_mat )
        self.vm_glcore.load_lights   ( self.shader_program )
        self.vm_glcore.load_fog      ( self.shader_program )

        # [EN] 2026-09-26 user request ("alterar a transparencia da
        # representacao de cartoon") -- same glBlendFunc/glDepthMask
        # dance as SurfaceRepresentation.draw_representation() above:
        # depth TEST stays on (so an opaque object in front still hides
        # this ribbon), but depth WRITES are switched off while
        # translucent so the ribbon's own far/near walls don't self-
        # occlude in a way that looks wrong through the transparency.
        alpha_loc = GL.glGetUniformLocation ( self.shader_program, "cartoon_alpha" )
        GL.glUniform1f ( alpha_loc, self.alpha )
        if self.alpha < 1.0:
            GL.glBlendFunc ( GL.GL_SRC_ALPHA, GL.GL_ONE_MINUS_SRC_ALPHA )
            GL.glEnable ( GL.GL_BLEND )
            GL.glDepthMask ( GL.GL_FALSE )

        GL.glBindVertexArray ( self.vao )
        GL.glDrawElements ( GL.GL_TRIANGLES, int ( len ( self.indexes ) ), GL.GL_UNSIGNED_INT, None )

        if self.alpha < 1.0:
            GL.glDepthMask ( GL.GL_TRUE )
            GL.glDisable ( GL.GL_BLEND )
        GL.glDisable ( GL.GL_DEPTH_TEST )

    def draw_background_sel_representation ( self ):
        """ [EN] 2026-09-26 user request ("tornar a representacao clicavel
        para a ferramenta de selecao/centralizacao") -- draws the Ca-trace
        built by _build_pick_trace_geometry() into the shared picking
        framebuffer (see VismolGLCore._pick(), which calls this on every
        active representation of every active object, then reads back ONE
        pixel and decodes it via vm_session.atom_dic_id -- fully generic,
        NOT modified for this: any representation whose sel-pass draws a
        real atom's colour.color_id at a real 3D position already resolves
        correctly through that existing code, no changes needed there).
        Line-based, not point-based, so a click landing ANYWHERE along the
        ribbon's own length (not just exactly on a Ca) still resolves to
        the nearest residue -- reuses the same generously-wide
        'line_width_selection' every other representation's own sel-pass
        already uses for exactly this "make the click target easier to
        hit" reason (see e.g. LinesRepresentation.draw_background_sel_
        representation() above). """
        if self.pick_coords is None or len ( self.pick_coords ) == 0:
            return
        # [EN] 2026-09-28 BUG FIX: the picking pass can run BEFORE this
        # representation's first draw_representation() (e.g. the Builder's
        # hover pick, which render() runs ahead of the normal draw loop),
        # or right after a rebuild() reset sel_vao -- sel_shader_program/
        # sel_vao were still None there (ctypes "wrong type" crash in
        # glUseProgram). Same lazy check every other representation's own
        # sel-pass already starts with.
        self._check_vao_and_vbos ( )
        GL.glEnable ( GL.GL_DEPTH_TEST )
        GL.glUseProgram ( self.sel_shader_program )
        line_width = self.vm_session.vm_config.gl_parameters.get ( "line_width_selection", 10 )
        GL.glLineWidth ( line_width )
        self.vm_glcore.load_matrices ( self.sel_shader_program, self.vm_object.model_mat )
        GL.glBindVertexArray ( self.sel_vao )
        GL.glDrawArrays ( GL.GL_LINES, 0, int ( len ( self.pick_coords ) ) )
        GL.glBindVertexArray ( 0 )
        GL.glDisable ( GL.GL_DEPTH_TEST )
        GL.glLineWidth ( 1 )
        GL.glUseProgram ( 0 )







'''
class DynamicBonds(Representation):
    """ Class doc """
    
    def __init__ (self, vismol_object, vismol_glcore, name="sticks", indexes=None, active=True):
        """ Class initialiser """
        super(DynamicBonds, self).__init__(vismol_object, vismol_glcore, name, active, indexes)
    
    def draw_representation(self):
        """ Function doc """
        self._check_vao_and_vbos ()
        self._enable_anti_alias_to_lines()
        GL.glUseProgram(self.shader_program)
        GL.glLineWidth(40/abs(self.vm_glcore.dist_cam_zrp))
        self.vm_glcore.load_matrices(self.shader_program, self.vm_object.model_mat)
        self.vm_glcore.load_fog(self.shader_program)
        self.vm_glcore.load_lights(self.shader_program)
        GL.glBindVertexArray(self.vao)
        
        if self.vm_glcore.modified_view:
            pass
        else:
            # This function checks if the number of the called frame will not exceed 
            # the limit of frames that each object has. Allowing two objects with 
            # different trajectory sizes to be manipulated at the same time within the 
            # glArea
            frame = self.vm_glcore.frame
            if frame < len(self.vm_object.dynamic_bons):
                self.define_new_indexes_to_vbo(self.vm_object.dynamic_bons[frame])
                self._set_coordinates_to_buffer(coord_vbo=True, sel_coord_vbo=False)
                GL.glDrawElements(GL.GL_LINES, int(len(self.vm_object.dynamic_bons[frame])*2), GL.GL_UNSIGNED_INT, None)
            else:
                self.define_new_indexes_to_vbo(self.vm_object.dynamic_bons[-1])
                self._set_coordinates_to_buffer(coord_vbo=True, sel_coord_vbo=False)
                GL.glDrawElements(GL.GL_LINES, int(len(self.vm_object.dynamic_bons[-1])*2), GL.GL_UNSIGNED_INT, None)
        
        GL.glBindVertexArray(0)
        self._disable_anti_alias_to_lines()
        GL.glLineWidth(1)
        GL.glUseProgram(0)
    
    def draw_background_sel_representation(self):
        """ Function doc """
        self._check_vao_and_vbos()
        self._disable_anti_alias_to_lines()
        GL.glUseProgram(self.sel_shader_program)
        GL.glLineWidth(20)
        GL.glEnable(GL.GL_DEPTH_TEST)
        self.vm_glcore.load_matrices(self.sel_shader_program, self.vm_object.model_mat)
        GL.glBindVertexArray(self.sel_vao)
        
        if self.vm_glcore.modified_view:
            pass
        else:
            # This function checks if the number of the called frame will not exceed 
            # the limit of frames that each object has. Allowing two objects with 
            # different trajectory sizes to be manipulated at the same time within the 
            # glArea
            frame = self.vm_glcore.frame
            if frame < len(self.vm_object.dynamic_bons):
                self.define_new_indexes_to_vbo(self.vm_object.dynamic_bons[frame])
                self._set_coordinates_to_buffer(coord_vbo=True, sel_coord_vbo=False)
                GL.glDrawElements(GL.GL_LINES, int(len(self.vm_object.dynamic_bons[frame])*2), GL.GL_UNSIGNED_INT, None)
            else:
                self.define_new_indexes_to_vbo(self.vm_object.dynamic_bons[-1])
                self._set_coordinates_to_buffer(coord_vbo=True, sel_coord_vbo=False)
                GL.glDrawElements(GL.GL_LINES, int(len(self.vm_object.dynamic_bons[-1])*2), GL.GL_UNSIGNED_INT, None)
        
        GL.glBindVertexArray(0)
        GL.glDisable(GL.GL_DEPTH_TEST)
        GL.glLineWidth(1)
        GL.glUseProgram(0)



class RibbonsRepresentation(Representation):
    """ Class doc """
    
    def __init__(self, vismol_object, vismol_glcore, name="ribbon", indexes=None, active=True):
        """ Class initialiser """
        super(RibbonsRepresentation, self).__init__(vismol_object, vismol_glcore, name, active, indexes)
        
        if self.vm_object.c_alpha_bonds == []:
            self.vm_object.get_backbone_indexes()
        
        indexes = []
        for bond in self.vm_object.c_alpha_bonds:
            indexes.append(bond.atom_index_i)
            indexes.append(bond.atom_index_j)
        
        if indexes == []:
            self.activate = False
        else:
            self.indexes = np.array(indexes, dtype=np.uint32)
    
    def draw_representation(self):
        """ Function doc """
        self._check_vao_and_vbos ()
        self._enable_anti_alias_to_lines()
        GL.glUseProgram(self.shader_program)
        ribbon_width = self.vm_session.vm_config.gl_parameters["ribbon_width"]
        ribbon_width = (ribbon_width*20)/abs(self.vm_glcore.dist_cam_zrp)/2
        GL.glLineWidth(ribbon_width)
        self.vm_glcore.load_matrices(self.shader_program, self.vm_object.model_mat)
        self.vm_glcore.load_fog(self.shader_program)
        GL.glBindVertexArray(self.vao)
        
        if self.vm_glcore.modified_view:
            pass
        else:
            # This function checks if the number of the called frame will not exceed 
            # the limit of frames that each object has. Allowing two objects with 
            # different trajectory sizes to be manipulated at the same time within the 
            # glArea
            self._set_coordinates_to_buffer(coord_vbo=True, sel_coord_vbo=False)
            GL.glDrawElements(GL.GL_LINES, int(len(self.vm_object.index_bonds)*2), GL.GL_UNSIGNED_INT, None)
        
        GL.glBindVertexArray(0)
        self._disable_anti_alias_to_lines()
        GL.glLineWidth(1)
        GL.glUseProgram(0)
    
    def draw_background_sel_representation(self):
        """ Function doc """
        self._check_vao_and_vbos()
        self._disable_anti_alias_to_lines()
        line_width = self.vm_session.vm_config.gl_parameters["line_width_selection"] 
        GL.glUseProgram(self.sel_shader_program)
        GL.glLineWidth(line_width)
        GL.glEnable(GL.GL_DEPTH_TEST)
        self.vm_glcore.load_matrices(self.sel_shader_program, self.vm_object.model_mat)
        GL.glBindVertexArray(self.sel_vao)
        
        if self.vm_glcore.modified_view:
            pass
        else:
            # This function checks if the number of the called frame will not exceed 
            # the limit of frames that each object has. Allowing two objects with 
            # different trajectory sizes to be manipulated at the same time within the 
            # glArea
            self._set_coordinates_to_buffer(coord_vbo=False, sel_coord_vbo=True)
            GL.glDrawElements(GL.GL_LINES, int(len(self.vm_object.index_bonds)*2), GL.GL_UNSIGNED_INT, None)
        
        GL.glBindVertexArray(0)
        GL.glDisable(GL.GL_DEPTH_TEST)
        GL.glLineWidth(1)
        GL.glUseProgram(0)



class ImpostorRepresentation(Representation):
    """ Class doc """
    
    def __init__ (self, vismol_object, vismol_glcore, name = "impostor", indexes=None, active=True, scale=1.0):
        """ Class initialiser """
        super(ImpostorRepresentation, self).__init__(vismol_object, vismol_glcore, name, active, indexes)
        self.scale = scale
    
    def draw_representation(self):
        """ Function doc """
        self._check_vao_and_vbos()
        self._enable_anti_alias_to_lines()
        GL.glUseProgram(self.shader_program)
        height = self.vm_glcore.height
        dist_cam_zrp = self.vm_glcore.dist_cam_zrp
        xyz_coords = self.vm_glcore.glcamera.get_modelview_position(self.vm_object.model_mat)
        u_campos = GL.glGetUniformLocation(self.shader_program, "u_campos")
        GL.glUniform3fv(u_campos, 1, xyz_coords)
        self.vm_glcore.load_lights(self.shader_program)
        self.vm_glcore.load_matrices(self.shader_program, self.vm_object.model_mat)
        self.vm_glcore.load_fog(self.shader_program)
        GL.glBindVertexArray(self.vao)

        if self.vm_glcore.modified_view:
            pass
        else:
            # This function checks if the number of the called frame will not exceed 
            # the limit of frames that each object has. Allowing two objects with 
            # different trajectory sizes to be manipulated at the same time within the 
            # glArea
            self._set_coordinates_to_buffer(coord_vbo=True, sel_coord_vbo=False)
            GL.glDrawElements(GL.GL_POINTS, len(self.vm_object.atoms), GL.GL_UNSIGNED_INT, None)
        
        GL.glBindVertexArray(0)
        self._disable_anti_alias_to_lines()
        GL.glPointSize(1)
        GL.glUseProgram(0)
    
    def draw_background_sel_representation(self):
        """ Function doc """
        self._check_vao_and_vbos()
        self._disable_anti_alias_to_lines()
        GL.glUseProgram(self.sel_shader_program)
        self.vm_glcore.load_matrices(self.sel_shader_program, self.vm_object.model_mat)
        GL.glEnable(GL.GL_DEPTH_TEST)
        GL.glBindVertexArray(self.sel_vao)
        
        if self.vm_glcore.modified_view:
            pass
        else:
            # This function checks if the number of the called frame will not exceed 
            # the limit of frames that each object has. Allowing two objects with 
            # different trajectory sizes to be manipulated at the same time within the 
            # glArea
            self._set_coordinates_to_buffer(coord_vbo=False, sel_coord_vbo=True)
            GL.glDrawElements(GL.GL_POINTS, len(self.vm_object.atoms), GL.GL_UNSIGNED_INT, None)
        
        GL.glBindVertexArray(0)
        GL.glDisable(GL.GL_DEPTH_TEST)
        GL.glPointSize(1)
        GL.glUseProgram(0)
'''


'''
class SurfaceRepresentation(Representation):
    """ Class doc """
    
    def __init__ (self, name = "surface", active = True, rep_type = "mol", vismol_object = None, vm_glcore = None, indexes = []):
        """ Class initialiser """
        self.name               = name
        self.active             = active
        self.type               = rep_type

        self.vm_object             = vismol_object
        self.vm_glcore             = vm_glcore
        
        

        
        
        # representation 	
        self.vao            = None
        self.ind_vbo        = None
        self.coord_vbo      = None
        self.norm_vbo       = None
        self.col_vbo        = None
        self.size_vbo       = None
           

        # bgrd selection   
        self.sel_vao        = None
        self.sel_ind_vbo    = None
        self.sel_coord_vbo  = None
        self.sel_col_vbo    = None
        self.sel_size_vbo   = None


        #     S H A D E R S
        self.shader_program     = None
        self.sel_shader_program = None
        self.read_surface_data()
    
    
    ##### sub 2 vev3 vectors
    def sub_vec3(self, a, b):
        c = [ a[0] - b[0],
              a[1] - b[1],
              a[2] - b[2] ]

        return c

    ## add 2 vectors and take the avg
    ## if a vector is still 0 we just take b
    def avg_add_vec3(self, a, b):
        if a[0] == 0.0 and a[1] == 0.0 and a[2] == 0.0 :
            return b

        c = [ (a[0] + b[0]) * 0.5 ,
              (a[1] + b[1]) * 0.5 ,
              (a[2] + b[2]) * 0.5 ]

        return c    

    ## make the cross product of 2 vectors
    def cross_vec3(self, a, b):
        c = [a[1]*b[2] - a[2]*b[1],
             a[2]*b[0] - a[0]*b[2],
             a[0]*b[1] - a[1]*b[0]]

        return c
    #############################################
        
    
    
    def read_surface_data(self):
        """ Function doc """
        #from random import random 
        #
        #[verts, tris, verts_gpu, tris_gpu] = edtsurf.calc_surface("/home/fernando/programs/EasyHybrid3/Coords/pdbs/1bx4_H.pdb")
        #self.coords2  = verts_gpu
        #self.indexes2 = tris_gpu
        #self.colors2  = []
        #
        #
        #size = len( self.coords2 )
        #for i in range(size):
        #    self.colors2.append(float(i/size) + random())
        
        rawdata = open("../EasyHybrid3/Coords/pdbs/1bx4.ply", "r")
        lines  = rawdata.readlines()
        
        self.coords2 = []
        self.colors2 = []
        self.normals2 = []
        self.indexes2 = []
        avg_normals_indexes = []
        
        
        for line in lines:
            line2 = line.split()
            
            if len(line2) == 6:
                #print (line2)
                self.coords2.append(float(line2[0]))
                self.coords2.append(float(line2[1]))
                self.coords2.append(float(line2[2]))
                                                  
                self.colors2.append(float(line2[3])/255)
                self.colors2.append(float(line2[4])/255)
                self.colors2.append(float(line2[5])/255)
                
                self.normals2.append(float(line2[0]))
                self.normals2.append(float(line2[1]))
                self.normals2.append(float(line2[2]))                
                avg_normals_indexes.append( ( 0.0 , 0.0 , 0.0 ) )  ### NEW !!! 

            if len(line2) == 7:
                
                self.indexes2.append(int(line2[1]))
                self.indexes2.append(int(line2[2]))
                self.indexes2.append(int(line2[3]))
                
        
        ## calculate normals and interpolate them (thanks a lot Kai)
        for i in range( 0 , len(self.indexes2) , 3 ):

            index_1 = self.indexes2[i] * 3;
            index_2 = self.indexes2[i+1] * 3;
            index_3 = self.indexes2[i+2] * 3;
            vertex_1 = ( self.coords2[index_1] , self.coords2[index_1+1] , self.coords2[index_1+2] )
            vertex_2 = ( self.coords2[index_2] , self.coords2[index_2+1] , self.coords2[index_2+2] )
            vertex_3 = ( self.coords2[index_3] , self.coords2[index_3+1] , self.coords2[index_3+2] )

            vec_p0_p1 = self.sub_vec3( vertex_2 , vertex_1 )
            vec_p0_p2 = self.sub_vec3( vertex_3 , vertex_1 )
            norm_vec  = self.cross_vec3( vec_p0_p1, vec_p0_p2 )

            vert_index_1 = self.indexes2[i] ;
            vert_index_2 = self.indexes2[i+1] ;
            vert_index_3 = self.indexes2[i+2] ;
            
            avg_normals_indexes[vert_index_1] = self.avg_add_vec3( avg_normals_indexes[vert_index_1] , norm_vec )
            avg_normals_indexes[vert_index_2] = self.avg_add_vec3( avg_normals_indexes[vert_index_2] , norm_vec )
            avg_normals_indexes[vert_index_3] = self.avg_add_vec3( avg_normals_indexes[vert_index_3] , norm_vec )


        ## set all new interpolated normals   
        for i in range( 0 , len(self.indexes2) , 1 ):
            index_1 = self.indexes2[i] * 3;

            self.normals2[index_1]   = avg_normals_indexes[self.indexes2[i]][0]
            self.normals2[index_1+1] = avg_normals_indexes[self.indexes2[i]][1]
            self.normals2[index_1+2] = avg_normals_indexes[self.indexes2[i]][2]





               
                

    def _make_gl_vao_and_vbos (self, indexes = None):
        """ Function doc """
        #if indexes is not None:
        #    pass
        #else:
        
        #dot_qtty  = int(len(self.vm_object.frames[0])/3)
        #indexes = []
        #for i in range(dot_qtty):
        #    indexes.append(i)
        

        self.shader_program     = self.vm_glcore.shader_programs[self.name]
        self.sel_shader_program = self.vm_glcore.shader_programs[self.name+"_sel"]
        
        #indexes = np.array(self.vm_object.index_bonds, dtype=np.uint32)
        #indexes = np.array(self.vm_object.idex, dtype=np.uint32)

        coords  = np.array(self.coords2, dtype=np.float32)
        colors  = np.array(self.colors2, dtype=np.float32)
        normals = np.array(self.normals2, dtype=np.float32)
        #indexes = range(0, len(self.coords2))     
        #indexes = np.array(indexes, dtype=np.uint32)
        indexes = np.array(self.indexes2, dtype=np.uint32)
        #print (indexes)


        self._make_gl_representation_vao_and_vbos (indexes    = indexes,
                                                   coords     = coords ,
                                                   colors     = colors ,
                                                   dot_sizes  = None   ,
                                                   normals    = normals
                                                   )
        
        #self.centr_vbo = GL.glGenBuffers(1)
        #GL.glBindBuffer(GL.GL_ARRAY_BUFFER, coords)
        #GL.glBufferData(GL.GL_ARRAY_BUFFER, coords.itemsize*len(coords), coords, GL.GL_STATIC_DRAW)
        #gl_center = GL.glGetAttribLocation(self.shader_program , "vert_centr")
        #GL.glEnableVertexAttribArray(gl_center)
        #GL.glVertexAttribPointer(gl_center, 3, GL.GL_FLOAT, GL.GL_FALSE, 3*coords.itemsize, ctypes.c_void_p(0))
        
        
        
        colors_idx = self.vm_object.color_indexes
        self._make_gl_sel_representation_vao_and_vbos (indexes    = indexes    ,
                                                       coords     = coords     ,
                                                       colors     = colors_idx ,
                                                       dot_sizes  = None       ,
                                                       )

    def draw_representation (self):
        """ Function doc """
        self._check_vao_and_vbos ()
        #self._enable_anti_alias_to_lines()
        
        
        
        
        GL.glEnable(GL.GL_DEPTH_TEST)
        GL.glEnable(GL.GL_CULL_FACE)
        GL.glCullFace(GL.GL_BACK)
        view = self.vm_glcore.glcamera.view_matrix
        
        GL.glUseProgram(self.shader_program )
        
        #print (self.vm_object.model_mat,view)
        
        m_normal = np.array(np.matrix(np.dot(view, self.vm_object.model_mat)).I.T)
        
        self.vm_glcore.load_matrices(self.shader_program , self.vm_object.model_mat)
        self.vm_glcore.load_lights  (self.shader_program )
        self.vm_glcore.load_fog     (self.shader_program )
        GL.glBindVertexArray(self.vao)
        
        
        
        
        
        
        
        
        
        """
        #print ("DotsRepresentation")
        height = self.vm_glcore.height
        
        GL.glUseProgram(self.shader_program)
        #1*self.height dot_size
        #GL.glLineWidth(40/abs(self.vm_glcore.dist_cam_zrp))
        GL.glPointSize(0.1*height/abs(self.vm_glcore.dist_cam_zrp)) # dot size not included yet
        self.vm_glcore.load_matrices(self.shader_program, self.vm_object.model_mat)
        self.vm_glcore.load_fog(self.shader_program)
        GL.glBindVertexArray(self.vao)
        """
        if self.vm_glcore.modified_view:
            pass
        
        else:
            """
            This function checks if the number of the called frame will not exceed 
            the limit of frames that each object has. Allowing two objects with 
            different trajectory sizes to be manipulated at the same time within the 
            glArea"""
            # self._set_coordinates_to_buffer(coord_vbo = True, sel_coord_vbo = False)
            #GL.glDrawElements(GL.GL_POINTS, int(len(self.indexes2)), GL.GL_UNSIGNED_INT, None)
            #GL.glDrawElements(GL.GL_LINE_LOOP, int(len(self.coords2)), GL.GL_UNSIGNED_INT, None)
            #GL.glDrawElements(GL.GL_LINE_STRIP, int(len(self.indexes2)), GL.GL_UNSIGNED_INT, None)
            GL.glDrawElements(GL.GL_TRIANGLES, int(len(self.indexes2)), GL.GL_UNSIGNED_INT, None)
        
        #GL.glBindVertexArray(0)
        #GL.glLineWidth(1)
        #GL.glUseProgram(0)
        #GL.glDisable(GL.GL_LINE_SMOOTH)
        #GL.glDisable(GL.GL_BLEND)
        GL.glDisable(GL.GL_DEPTH_TEST)
        
            
    def draw_background_sel_representation  (self):
        """ Function doc """
        pass


class WiresRepresentation(Representation):
    """ Class doc """
    
    def __init__ (self, name = "wires", active = True, rep_type = "mol", vismol_object = None, vm_glcore = None, indexes = []):
        """ Class initialiser """
        self.name               = name
        self.active             = active
        self.type               = rep_type
        self.vm_object             = vismol_object
        self.vm_glcore             = vm_glcore
        
        # representation    
        self.vao            = None
        self.ind_vbo        = None
        self.coord_vbo      = None
        self.col_vbo        = None
        self.size_vbo       = None
        
        # bgrd selection   
        self.sel_vao        = None
        self.sel_ind_vbo    = None
        self.sel_coord_vbo  = None
        self.sel_col_vbo    = None
        self.sel_size_vbo   = None

        #     S H A D E R S
        self.shader_program     = None
        self.sel_shader_program = None
        self.read_surface_data()
    
    def read_surface_data(self):
        """ Function doc """
        rawdata = open("../EasyHybrid3/Coords/pdbs/1bx4.ply", "r")
        lines  = rawdata.readlines()
        
        self.coords2 = []
        self.colors2 = []
        self.indexes2 = []
        
        for line in lines:
            line2 = line.split()
            if len(line2) == 6:
                self.coords2.append(float(line2[0]))
                self.coords2.append(float(line2[1]))
                self.coords2.append(float(line2[2]))
                self.colors2.append(float(line2[3])/255)
                self.colors2.append(float(line2[4])/255)
                self.colors2.append(float(line2[5])/255)
            if len(line2) == 7:
                self.indexes2.append(int(line2[1]))
                self.indexes2.append(int(line2[2]))
                self.indexes2.append(int(line2[3]))

    def _make_gl_vao_and_vbos (self, indexes = None):
        """ Function doc """
        self.shader_program     = self.vm_glcore.shader_programs[self.name]
        self.sel_shader_program = self.vm_glcore.shader_programs[self.name+"_sel"]
        coords  = np.array(self.coords2, dtype=np.float32)
        colors  = np.zeros(len(self.colors2))
        indexes = np.array(self.indexes2, dtype=np.uint32)
        self._make_gl_representation_vao_and_vbos (indexes    = indexes,
                                                   coords     = coords ,
                                                   colors     = colors ,
                                                   dot_sizes  = None   ,
                                                   )
        colors_idx = self.vm_object.color_indexes
        self._make_gl_sel_representation_vao_and_vbos (indexes    = indexes    ,
                                                       coords     = coords     ,
                                                       colors     = colors_idx ,
                                                       dot_sizes  = None       ,
                                                       )

    def draw_representation (self):
        """ Function doc """
        self._check_vao_and_vbos ()
        pass
        #GL.glEnable(GL.GL_DEPTH_TEST)
        #GL.glEnable(GL.GL_CULL_FACE)
        #GL.glCullFace(GL.GL_BACK)
        #
        ##LineWidth = (80/abs(self.vm_glcore.dist_cam_zrp)/2)**0.5  #40/abs(self.vm_glcore.dist_cam_zrp)
        ##GL.glLineWidth(2)
        #
        #
        #view = self.vm_glcore.glcamera.view_matrix
        #GL.glUseProgram(self.shader_program )
        #m_normal = np.array(np.matrix(np.dot(view, self.vm_object.model_mat)).I.T)
        #self.vm_glcore.load_matrices(self.shader_program , self.vm_object.model_mat)
        ##self.vm_glcore.load_lights  (self.shader_program )
        #self.vm_glcore.load_fog     (self.shader_program )
        #GL.glBindVertexArray(self.vao)
        #if self.vm_glcore.modified_view:
        #    pass
        #
        #else:
        #    """
        #    This function checks if the number of the called frame will not exceed 
        #    the limit of frames that each object has. Allowing two objects with 
        #    different trajectory sizes to be manipulated at the same time within the 
        #    glArea"""
        #    # self._set_coordinates_to_buffer(coord_vbo = True, sel_coord_vbo = False)
        #    GL.glDrawElements(GL.GL_TRIANGLES, int(len(self.indexes2)), GL.GL_UNSIGNED_INT, None)
        #GL.glDisable(GL.GL_DEPTH_TEST)
        
    def draw_background_sel_representation  (self):
        """ Function doc """
        pass


class LabelRepresentation:
    """ Class doc """
    
    def __init__ (self, name = "labels", active = True, rep_type = "mol", vismol_object = None, vm_glcore = None, indexes = []):
        """ Class initialiser """
        self.vm_object = vismol_object
        self.name   = name
        self.active = True
        self.vm_glcore = vm_glcore
        
        self.chars     = 0 
        #self._check_vao_and_vbos()
        
    def _check_vao_and_vbos (self, indexes = None):
        """ Function doc """
        if self.vm_object.vm_font.vao is None:
            self.vm_object.vm_font.make_freetype_font()
            self.vm_object.vm_font.make_freetype_texture(self.vm_glcore.freetype_program)
        
        if self.chars == 0:
            dprint("self._build_buffer()")
            self._build_buffer()
        GL.glBindBuffer(GL.GL_ARRAY_BUFFER, self.vm_object.vm_font.vbos[0])
        GL.glBufferData(GL.GL_ARRAY_BUFFER, self.xyz_pos.itemsize*len(self.xyz_pos), self.xyz_pos, GL.GL_DYNAMIC_DRAW)
        GL.glBindBuffer(GL.GL_ARRAY_BUFFER, self.vm_object.vm_font.vbos[1])
        GL.glBufferData(GL.GL_ARRAY_BUFFER, self.uv_coords.itemsize*len(self.uv_coords), self.uv_coords, GL.GL_DYNAMIC_DRAW)
        GL.glBindBuffer(GL.GL_ARRAY_BUFFER, 0)


    def _build_buffer (self, indexes = None):
        self.chars     = 0
        self.xyz_pos   = []
        self.uv_coords = []
        for atom in self.vm_object.atoms:
            
            texto = atom.name
            point = np.array(atom.coords (self.vm_glcore.frame),np.float32)
            point = np.array((point[0],point[1],point[2],1),np.float32)
            point = np.dot(point, self.vm_object.model_mat)

            GL.glBindTexture(GL.GL_TEXTURE_2D, self.vm_object.vm_font.texture_id)
            for i,c in enumerate(texto):
                self.chars += 1
                c_id = ord(c)
                x = c_id%16
                y = c_id//16-2
                self.xyz_pos.append(point[0]+i*self.vm_object.vm_font.char_width)
                self.xyz_pos.append(point[1])
                self.xyz_pos.append(point[2])

                self.uv_coords.append(x*self.vm_object.vm_font.text_u)
                self.uv_coords.append(y*self.vm_object.vm_font.text_v)
                self.uv_coords.append((x+1)*self.vm_object.vm_font.text_u)
                self.uv_coords.append((y+1)*self.vm_object.vm_font.text_v)
            #print(texto)
        #print("xyz_pos  ",len(self.xyz_pos))
        #print("uv_coords",len(self.uv_coords))
        #print("atoms    ",len(self.vm_object.atoms))
        #print("chars    ",self.chars)
        
        self.xyz_pos   = np.array(self.xyz_pos  , np.float32)
        self.uv_coords = np.array(self.uv_coords, np.float32)
        


    
    
    def draw_representation (self):
        """ Function doc """
        self._check_vao_and_vbos()
        
        GL.glDisable(GL.GL_DEPTH_TEST)
        GL.glEnable(GL.GL_BLEND)
        GL.glBlendFunc(GL.GL_SRC_ALPHA, GL.GL_ONE_MINUS_SRC_ALPHA)
        GL.glUseProgram(self.vm_glcore.freetype_program)
        
        self.vm_object.vm_font.load_matrices(self.vm_glcore.freetype_program, self.vm_glcore.glcamera.view_matrix, self.vm_glcore.glcamera.projection_matrix)
        self.vm_object.vm_font.load_font_params(self.vm_glcore.freetype_program)
        
        GL.glBindVertexArray(self.vm_object.vm_font.vao)
        GL.glDrawArrays(GL.GL_POINTS, 0, self.chars)
        GL.glDisable(GL.GL_BLEND)
        GL.glBindVertexArray(0)
        GL.glUseProgram(0)

        

    def draw_background_sel_representation  (self):
        """ Function doc """
        pass

'''
