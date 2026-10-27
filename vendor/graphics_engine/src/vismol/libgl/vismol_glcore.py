#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
#  vismol_glcore.py
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
import time
import numpy as np
from OpenGL import GL
from logging import getLogger
from vismol.libgl.glaxis import GLAxis
from vismol.libgl.glcamera import GLCamera
from vismol.libgl.vismol_font import VismolFont
from vismol.libgl.vismol_font import resolve_font_path, list_available_fonts
from vismol.libgl.vismol_font import DEFAULT_FONT_FILE, DEFAULT_FONT_SIZE
from vismol.libgl.selection_box import SelectionBox
import vismol.libgl.shapes as shapes
import vismol.libgl.shaders.pick as shaders_pick
import vismol.libgl.shaders.lighting as shaders_lighting
import vismol.libgl.shaders.dots as shaders_dots
import vismol.libgl.shaders.lines as shaders_lines
import vismol.libgl.shaders.wires as shaders_wires
import vismol.libgl.shaders.sticks as shaders_sticks
import vismol.libgl.shaders.cartoon as shaders_cartoon
import vismol.libgl.shaders.surface as shaders_surface
import vismol.libgl.shaders.spheres as shaders_spheres
import vismol.libgl.shaders.impostor as shaders_impostor
import vismol.libgl.shaders.nonbonded as shaders_nonbonded
import vismol.libgl.shaders.vm_freetype as shaders_vm_freetype
import vismol.libgl.shaders.dashed_lines as shaders_dashed_lines
import vismol.utils.matrix_operations as mop

from vismol.core.vismol_object import VismolObject
from vismol.model.atom import Atom
import math

# [EN] BUG FIXED (reported by the user: real crash on app startup --
# "ImportError: cannot import name 'VismolGLCore' from
# 'vismol.libgl.vismol_glcore'"). This USED to be a top-level
# `from gui.windows.builder import click_mode/atom_ops/empty_object`
# right here -- reasoned (at length, in this project's own history) to
# be safe based on reading every file in the dependency chain and
# finding no import that led back to this module. That static reading
# turned out to be incomplete somehow (the exact live import order
# vismol_gtkwidget.py -> vismol_glcore.py -> gui.windows.builder.* takes
# during actual app startup isn't fully visible from source alone, and
# apparently loops back to THIS module before the VismolGLCore class
# below is defined, in a way the static check missed). Reverted to
# local/deferred imports (see mouse_pressed()/mouse_motion()/
# mouse_released()/render(), the only 4 methods that actually use
# click_mode/atom_ops/empty_object) -- imported lazily, inside each of
# those methods, the FIRST time any of them actually runs, which is
# necessarily well after the whole application (including
# vismol_gtkwidget.py) has already finished loading, so this can no
# longer participate in any load-time cycle, regardless of the exact
# mechanism of the one that broke here. Costs nothing extra after the
# first call either (Python caches already-imported modules in
# sys.modules; a repeated `import` is just a fast dict lookup, not a
# re-execution) -- same reasoning that justified the top-level import
# in the first place, just without the load-time risk.

logger = getLogger(__name__)


class _OptionalBuilderFallback:
    """No-op fallback for EasyHybrid-only builder hooks in standalone VisMol."""

    def _read_depth_and_atom_at_pixel(self, *_args, **_kwargs):
        return None, None

    def __getattr__(self, _name):
        def _noop(*_args, **_kwargs):
            return False
        return _noop


_BUILDER_FALLBACK = _OptionalBuilderFallback()


def _builder_hooks():
    try:
        from gui.windows.builder import click_mode, atom_ops, empty_object
        return click_mode, atom_ops, empty_object
    except ImportError:
        return _BUILDER_FALLBACK, _BUILDER_FALLBACK, _BUILDER_FALLBACK


class VismolGLCore:
    
    def __init__(self, widget, vismol_session=None, width=640.0, height=420.0):
        """ Constructor of the class.
            
            Keyword arguments:
            vismol_session - 
        """
        self.parent_widget = widget
        self.vm_session = vismol_session
        self.vm_config = self.vm_session.vm_config
        self.width = np.float32(width)
        self.height = np.float32(height)
        self.shader_programs = {}
        self.core_shader_programs = {}
        self.representations_available = self.vm_config.representations_available
        self.sphere_selection = None
        # Cache for glGetUniformLocation results, keyed by (program, name).
        # Uniform locations are fixed once a program is linked, so looking
        # them up by string every frame (the previous behaviour) is wasted
        # driver/CPU work. This cache is invalidated whenever shaders are
        # (re)compiled - see create_gl_programs().
        self._uniform_loc_cache = {}
        # Cache do ULTIMO valor enviado para cada (program, uniform_name).
        # fog/light/antialias sao IDENTICOS para todos os objetos do frame e
        # raramente mudam entre frames, mas antes eram re-enviados via
        # glUniform* dentro de CADA draw_representation (dezenas de chamadas
        # GL por objeto, todo frame). Aqui guardamos o valor ja residente no
        # programa e so re-enviamos no diff. Invalidado em create_gl_programs
        # (relink zera os uniforms no driver). Ver _uniform_changed().
        self._uniform_value_cache = {}
        # --- Medidor de frame-time (Gargalo: instrumentacao) ---------------
        # Leve e opcional. Quando self.show_fps == True, render() acumula o
        # tempo de cada frame e imprime FPS medio + ms/frame a cada
        # self._fps_report_every frames. Default desligado: nenhum custo alem
        # de um 'if' por frame quando False. Para medir, basta setar
        # glcore.show_fps = True em runtime (ou aqui). Totalmente reversivel.
        self.show_fps = False
        self._fps_report_every = 60      # reporta a cada N frames
        self._fps_frame_count = 0        # frames desde o ultimo report
        self._fps_accum_time = 0.0       # tempo acumulado (s) desde o report
        self._fps_last_t = None          # timestamp do frame anterior
        # [EN] 2026-10-13 multisample antialiasing: Gtk.GLArea gives us a
        # single-sample framebuffer (GL_SAMPLES == 0) and has no API to ask
        # for more, so the scene is drawn into our own multisampled
        # framebuffer and resolved into the widget's one at the end of
        # render() -- see _msaa_begin()/_msaa_end().
        self._msaa_fbo = None
        self._msaa_color_rbo = None
        self._msaa_depth_rbo = None
        self._msaa_state = None          # (width, height, samples) of the storage
        self._msaa_target = None         # (framebuffer, width, height) being resolved into
        self._msaa_failed = False        # a driver refused it: do not retry every frame
        # -------------------------------------------------------------------
        
        
    def initialize(self):
        """ Enables the buffers and other charasteristics of the OpenGL context.
            sets the initial projection, view and model matrices
            
            self.flag -- Needed to only create one OpenGL program, otherwise a bunch of
                         programs will be created and use system resources. If the OpenGL
                         program will be changed change this value to True
        """
        self.model_mat = np.identity(4, dtype=np.float32) # Not sure if this is used :S
        # --- UBO de matrizes de camera (Gargalo 2) -------------------------
        # view_mat e proj_mat sao IGUAIS para todos os objetos do frame, mas
        # antes eram re-enviados por objeto via glUniformMatrix4fv (2 uploads
        # x N objetos por frame). Agora vivem num Uniform Buffer Object unico,
        # atualizado 1x por frame e compartilhado por todos os shaders via o
        # binding point CAMERA_UBO_BINDING. model_mat segue como uniform
        # por-objeto (muda por objeto, nao entra aqui).
        # Layout std140: dois mat4 contiguos = 2 * 64 = 128 bytes.
        #   offset   0: mat4 view_mat
        #   offset  64: mat4 proj_mat
        self.CAMERA_UBO_BINDING = 0          # binding point fixo
        self._camera_ubo = None              # id do buffer GL (lazy init)
        self._camera_ubo_size = 128          # 2 mat4 float32
        # -------------------------------------------------------------------
        self.zero_reference_point = np.zeros(3, dtype=np.float32)
        self.glcamera = GLCamera(self.vm_config.gl_parameters["field_of_view"],
                                 self.width / self.height,
                                 np.array([0,0,10], dtype=np.float32),
                                 self.zero_reference_point)
        # [EN] 2026-09-26 user request ("um parametro para a extensao do
        # fog") -- see GLCamera.fog_range_scale's own docstring.
        # 2026-09-28: plus the Preferences > General > Fog settings
        # (enable/intensity/custom start-end), see GLCamera.apply_fog_parameters().
        self.glcamera.apply_fog_parameters(self.vm_config.gl_parameters)

        # Font family/size for the labels drawn in the glArea can be
        # customized by the user in the Preferences window ("Labels Font
        # (glArea)"). 'label_font_file'/'label_font_size' drive the
        # PICKING labels (#1 #2 #3 #4, self.vm_font below) -- this is the
        # original preference key, kept as-is for backward compatibility.
        # Distance labels (self.vm_font_dist) now have their OWN,
        # independent size via 'pk_dist_label_font_size' (falling back to
        # the picking size when unset), but share the SAME font family as
        # picking labels (one combo box in Preferences).
        #
        # ALL glArea labels (picking, distance, atom labels) share a
        # single "zoom_sensitivity" (0.0..1.0) set on the Viewer >
        # General tab in Preferences ('labels_zoom_sensitivity'): 0.0
        # keeps every label a constant size on screen no matter how far
        # the camera dollies in/out (the default, see
        # VismolFont.zoom_sensitivity); 1.0 restores the pre-billboard-
        # refactor behavior where a label's on-screen size shrinks as the
        # camera moves away and grows as it gets closer, exactly like the
        # rest of the 3D scene. Values in between blend the two.
        _label_font_file = self.vm_config.gl_parameters.get("label_font_file", DEFAULT_FONT_FILE)
        _label_font_size = self.vm_config.gl_parameters.get("label_font_size", DEFAULT_FONT_SIZE)
        _label_font_path = resolve_font_path(_label_font_file)
        _labels_zoom_sensitivity = self.vm_config.gl_parameters.get("labels_zoom_sensitivity", 1.0)
        
        _dist_font_size = self.vm_config.gl_parameters.get("pk_dist_label_font_size", _label_font_size)
        
        self.vm_font        = VismolFont(font_file=_label_font_path, char_width=_label_font_size,
                                          char_height=_label_font_size, color=[1, 1, 1, 1])
        self.vm_font.zoom_sensitivity = _labels_zoom_sensitivity
        self.vm_font_static = VismolFont(font_file=_label_font_path, char_width=_label_font_size,
                                          char_height=_label_font_size, color=[1, 1, 1, 1])
        self.vm_font_dist   = VismolFont(char_res=264, font_file=_label_font_path,
                                          char_width=_dist_font_size, char_height=_dist_font_size,
                                          color = [1, 1, 1, 1])
        self.vm_font_dist.zoom_sensitivity = _labels_zoom_sensitivity
        # [EN] 2026-10-04: atom index labels of the Builder's DYFF Parameters
        # window (builder_atom_types_labeled_atoms) -- same family, size and
        # zoom sensitivity as the picking labels, its own (white) colour
        self.vm_font_index  = VismolFont(font_file=_label_font_path, char_width=_label_font_size,
                                          char_height=_label_font_size, color=[1, 1, 1, 1])
        self.vm_font_index.zoom_sensitivity = _labels_zoom_sensitivity
        
        self.axis = GLAxis(vm_glcore = self)
        self.selection_box = SelectionBox()
        self.parent_widget.set_has_depth_buffer(True)
        self.parent_widget.set_has_alpha(True)
        self.scroll = self.vm_config.gl_parameters["scroll_step"]
        # per-tick exponential rate applied to the current slab WIDTH
        # (z_far - z_near) -- see mouse_scroll(); kept separate from
        # self.scroll, which Ctrl+scroll still uses unchanged for its own
        # fixed-step model-matrix translation.
        self.scroll_step_fraction = self.vm_config.gl_parameters.get("scroll_step_fraction", 0.138)
        self.bckgrnd_color = self.vm_config.gl_parameters["background_color"]
        #                       Light Parameters                                
        self.apply_light_parameters()
        #                              Variables                                
        self.right = self.width / self.height
        self.left = -self.right
        self.top = np.float32(1.0)
        self.bottom = np.float32(-1.0)
        self.button = None
        self.dist_cam_zrp = np.linalg.norm(self.glcamera.get_position())
        self.shader_flag = True
        self.modified_data = False # Used somewhere?
        self.updated_coords = False
        self.dragging = False
        self.editing_mols = False
        self.show_axis = True
        self.ctrl = False
        self.shift = False
        # [EN] macOS trackpad fix: set True only by the Meta_L/R and
        # Super_L/R key handlers in vismol_gtkwidget.py, which only exist
        # on macOS -- stays permanently False on Linux/Windows, so the
        # "and not self.cmd"/"or (... self.cmd ...)" gating below is a
        # no-op there and mouse_pressed's behavior is unchanged.
        self.cmd = False
        self.atom_picked = None
        self.selection_box_picking = False
        self.picking = False
        self.picking_x, self.picking_y = None, None
        self.show_selection_box = False
        self.show_selection_box_x, self.show_selection_box_y = None, None
        self.mouse_x, self.mouse_y = np.float32(0.0), np.float32(0.0)
        self.mouse_rotate, self.mouse_zoom, self.mouse_pan = False, False, False
        self._macos_buttons_down = 0
        self.drag_pos_x, self.drag_pos_y, self.drag_pos_z = None, None, None
    
    def resize_window(self, width, height):
        """ Resizing function, takes the widht and height of the widget
            and modifies the view in the camera acording to the new values
        
            Keyword arguments:
            width -- Actual width of the window
            height -- Actual height of the window
        """
        self.width = np.float32(width)
        self.height = np.float32(height)
        self.right = self.width / self.height
        self.left = -self.right
        self.center_x = self.width / 2.0
        self.center_y = self.height / 2.0
        if not hasattr(self, "glcamera"):
            return
        self.glcamera.viewport_aspect_ratio = self.width / self.height
        # [EN] 2026-09-28 BUG FIX (macOS report, not macOS-specific): the
        # slab/clipping zoom deliberately lets the STORED z_near go below
        # zero (it keeps tracking the slab), and every other place that
        # builds the projection clamps it with max(z_near, min_znear) --
        # this one didn't, so a resize after widening the clipping built
        # an invalid perspective matrix (near <= 0).
        _proj_mat = mop.my_glPerspectivef(self.glcamera.field_of_view,
                                          self.glcamera.viewport_aspect_ratio,
                                          max(self.glcamera.z_near, self.glcamera.min_znear),
                                          self.glcamera.z_far)
        self.glcamera.set_projection_matrix(_proj_mat)

    def _macos_modes_from_buttons(self, left, middle, right):
        """ Derives rotate/zoom/pan from a set of *currently held* button
            booleans. Used for both the macOS press and motion handlers so
            the same live bitmask is always the single source of truth --
            see the note in mouse_motion() for why that matters.
        """
        if left and right and not middle:
            # Both trackpad sides held at once: treat as zoom, anchored
            # to the viewport center, regardless of self.cmd.
            return False, True, False
        rotate = left   and not (middle or right)
        zoom   = right  and not (middle or left) and not self.cmd
        pan    = (middle and not (right or left)) \
                 or (right and self.cmd and not (middle or left))
        return rotate, zoom, pan

    def mouse_pressed(self, button_number, mouse_x, mouse_y, button_state=None):
        """ Function doc
        """
        left   = np.int32(button_number) == 1
        middle = np.int32(button_number) == 2
        right  = np.int32(button_number) == 3
        # [EN] macOS trackpad fix: MacBook trackpads have no middle
        # button and no default gesture emulating one, so plain
        # middle-drag pan was unreachable there. Cmd+right-drag (Cmd +
        # two-finger-drag on a trackpad) now also pans, matching PyMOL's
        # own Cmd+Right convention on macOS -- self.cmd is only ever set
        # True by macOS-only key handlers (see vismol_gtkwidget.py), so
        # this is a no-op everywhere else: plain right-drag still zooms,
        # plain middle-drag still pans, unchanged.
        if button_state is not None:
            # [EN] macOS reports the very start of a two-finger
            # secondary-click gesture ambiguously -- it can deliver this
            # press as a plain left button (1) for the first event or
            # two, before "upgrading" it to the right button (3) once it
            # finishes disambiguating one- vs two-finger contact. Deciding
            # the gesture mode ONCE here, from just this press's button
            # number, used to latch onto whatever button macOS reported
            # first (often rotate) and never re-evaluate -- so a zoom
            # gesture that started ambiguous stayed stuck rotating for
            # its whole drag, which is what read as "zoom drifting
            # sideways". mouse_motion() now recomputes the mode from the
            # live button bitmask on every single motion event instead,
            # so this initial guess self-corrects within one frame as
            # soon as macOS reports the real button. Only used to seed
            # the mode for the (rare) case a press is immediately
            # followed by a release with no motion in between.
            button_mask = {1: 1, 2: 2, 3: 4}.get(int(button_number), 0)
            self._macos_buttons_down |= button_mask
            self.mouse_rotate, self.mouse_zoom, self.mouse_pan = \
                self._macos_modes_from_buttons(
                    bool(self._macos_buttons_down & 1),
                    bool(self._macos_buttons_down & 2),
                    bool(self._macos_buttons_down & 4))
        else:
            self.mouse_rotate = left   and not (middle or right)
            self.mouse_zoom   = right  and not (middle or left) and not self.cmd
            self.mouse_pan    = (middle and not (right or left)) \
                                 or (right and self.cmd and not (middle or left))
        self.mouse_x = np.float32(mouse_x)
        self.mouse_y = np.float32(mouse_y)
        self.drag_pos_x, self.drag_pos_y, self.drag_pos_z = self._mouse_pos(self.mouse_x, self.mouse_y)
        self.dragging = False
        if left:
            
            if  self.ctrl:
                pass
                #print(self.shift, self.ctrl)
            
            if self.shift and not self.ctrl:
                #print('picking_selection_mode:',self.vm_session.picking_selection_mode)
                #
                if self.vm_session.picking_selection_mode: # bachega 06 / 18 /2025
                    #print(mouse_x, mouse_y, button_number)
                    pass 
                else: # bachega 06 / 18 /2025
                    self.show_selection_box = True
                    self.selection_box.start = self.get_viewport_pos(mouse_x, mouse_y)
                    self.selection_box.end = self.get_viewport_pos(mouse_x, mouse_y)
                    self.selection_box.update_points()
                    self.selection_box_x = mouse_x
                    self.selection_box_y = self.height - mouse_y

            # [EN] Builder "click-and-drag to create a bonded atom" -- see
            # gui/windows/builder/click_mode.py (click_mode.start_bond_drag() /
            # click_mode.update_bond_drag() / click_mode.finish_bond_drag()). For the "add"
            # tool, plain press (not shift, same reasoning as the
            # click-to-place/click-to-replace hook in mouse_released():
            # shift always falls through to normal multi-select so 'b'
            # add-bond-between-two-selected-atoms keeps working
            # unchanged; not ctrl either -- Ctrl+click is reserved for
            # setting a BOND's order (to the sidebar's current
            # selection), see the mouse_released() hook below and
            # click_mode.apply_selected_bond_order()).
            #
            # ALSO covers the "move" tool -- a plain drag there means
            # "reposition this existing atom" instead of "create a new
            # bonded one" (see click_mode.start_atom_drag(), and the
            # tool branch added to the candidate-resolution block in
            # mouse_motion() below that decides which of the two this
            # turns into). Sharing this SAME press-recording/candidate
            # mechanism between both tools is safe because they're
            # mutually exclusive -- only one tool is ever active at a
            # time -- so there's no ambiguity about which action a given
            # drag should resolve to.
            #
            # Deliberately does NOT touch mouse_rotate here (was an
            # earlier version's bug, fixed after the user reported the
            # camera not rotating AT ALL anymore while the "add" tool was
            # selected, even over empty space): we can't know yet whether
            # this press landed on an atom or on empty space (the actual
            # pick has to happen inside render() -- see the GL-context
            # constraint documented next to builder_placing_atom in
            # render()), so suppressing rotation unconditionally here
            # would kill normal camera-rotate-by-drag for every plain
            # empty-space drag too, not just the bond-drag case. Instead:
            # just remember the press position; mouse_motion() only
            # actually starts the bond-drag (and only THEN turns rotation
            # off) once real mouse movement confirms the press resolved
            # to an atom AND the user is genuinely dragging, not just
            # clicking -- see mouse_motion() and the render() hook below.
            if ( getattr ( self.vm_session, "builder_atom_mode", False )
                 and getattr ( self.vm_session, "builder_tool", "add" ) in ( "add", "move" )
                 and not self.shift and not self.ctrl ):
                self.builder_press_x = np.float32 ( mouse_x )
                self.builder_press_y = np.float32 ( mouse_y )
                self.builder_checking_press = True

            # [EN] Builder "Ctrl+drag to reposition an existing atom" --
            # see gui/windows/builder/click_mode.py (click_mode.start_atom_drag() /
            # click_mode.update_atom_drag() / click_mode.finish_atom_drag()). Same reasoning as
            # the plain-drag press recording just above (deferred pick,
            # no mouse_rotate touched here) -- just gated on self.ctrl
            # instead of its absence. A plain Ctrl+click (no real drag)
            # on a BOND still means "cycle its order" (see the
            # mouse_released() hook below) -- these two never conflict,
            # since one only ever resolves to an ATOM and the other only
            # ever resolves to a BOND's on-screen line.
            if ( getattr ( self.vm_session, "builder_atom_mode", False )
                 and getattr ( self.vm_session, "builder_tool", "add" ) == "add"
                 and self.ctrl and not self.shift ):
                self.builder_ctrl_press_x = np.float32 ( mouse_x )
                self.builder_ctrl_press_y = np.float32 ( mouse_y )
                self.builder_ctrl_checking_press = True

        if middle:
            self.picking_x = np.float32(mouse_x)
            self.picking_y = np.float32(mouse_y)
            self.picking = True
            self.parent_widget.queue_draw()
        if right:
            self.picking_x = np.float32(mouse_x)
            self.picking_y = np.float32(mouse_y)
            self.picking = True
            self.parent_widget.queue_draw()
    
    def mouse_released(self, button_number, mouse_x, mouse_y, button_state=None):
        """ Function doc
        int(event.button)
        
        info      = menu header info
        menu_type = "pick_menu" / "bg_menu" / "sele_menu" / "ob_menu"
        
        """
        # [EN] Local/deferred import -- see the note at the top of this
        # file (where this used to be a module-level import) for why:
        # importing here, the first time this method actually runs
        # (necessarily well after the whole app has finished loading),
        # avoids the load-time circular import that broke when this was
        # a top-level import instead. Cheap on every subsequent call --
        # Python caches already-imported modules in sys.modules.
        click_mode, atom_ops, empty_object = _builder_hooks()
        left   = np.int32(button_number) == 1
        middle = np.int32(button_number) == 2
        right  = np.int32(button_number) == 3
        if button_state is not None:
            button_mask = {1: 1, 2: 2, 3: 4}.get(int(button_number), 0)
            self._macos_buttons_down &= ~button_mask
            if self._macos_buttons_down == 0:
                self.mouse_rotate = False
                self.mouse_zoom = False
                self.mouse_pan = False
            else:
                # One button released but another is still held (e.g. a
                # both-sides zoom where only one side lifted) -- fall
                # back to whatever mode the still-held button(s) imply,
                # same live-bitmask source of truth as mouse_motion().
                self.mouse_rotate, self.mouse_zoom, self.mouse_pan = \
                    self._macos_modes_from_buttons(
                        bool(self._macos_buttons_down & 1),
                        bool(self._macos_buttons_down & 2),
                        bool(self._macos_buttons_down & 4))
        else:
            self.mouse_rotate = False
            self.mouse_zoom = False
            self.mouse_pan = False
        # [EN] Builder "click-and-drag to create a bonded atom" -- checked
        # FIRST and returns immediately, same "early-exit, don't touch
        # anything below" structure as the builder_placing_atom hook
        # further down. Checking builder_bond_drag_active explicitly
        # here, rather than relying on self.dragging's value, avoids ANY
        # ambiguity with the plain-click builder_placing_atom path below
        # (which would otherwise ALSO fire on this same release and
        # place/replace a second atom on top of the one just finished
        # being dragged).
        if getattr ( self.vm_session, "builder_bond_drag_active", False ):
            click_mode.finish_bond_drag ( self )
            self.dragging = False
            self.parent_widget.queue_draw ( )
            return
        # [EN] Builder "Ctrl+drag to reposition an existing atom" --
        # checked next, same early-exit reasoning as builder_bond_drag_
        # active just above. MUST come before the Ctrl+click-a-bond
        # handling right below: without this, releasing after a genuine
        # Ctrl+drag would fall through into click_mode.find_bond_at_pixel() at the
        # RELEASE position and could cycle an unrelated bond's order as
        # an unintended side effect of what was actually just a
        # reposition gesture.
        if getattr ( self.vm_session, "builder_ctrl_drag_active", False ):
            click_mode.finish_atom_drag ( self )
            self.dragging = False
            self.parent_widget.queue_draw ( )
            return
        # [EN] Builder "Ctrl+click on a bond -> set its order to the
        # sidebar's current selection" -- see click_mode.
        # apply_selected_bond_order()/click_mode.find_bond_at_pixel(). Pure Python/
        # numpy (a 2D-projection distance check, not a GPU pick -- see
        # click_mode.find_bond_at_pixel()'s own docstring for why), so unlike
        # builder_placing_atom/builder_checking_press below, this needs
        # NO deferral to render() for a GL context -- runs directly here.
        # Scoped to ONLY vm_session.builder_target_object, matching the
        # "only one object is editable at a time" design. Checked before
        # self.dragging/builder_placing_atom below and returns immediately
        # so Ctrl+click never ALSO tries to place/replace an atom on the
        # same click (mouse_pressed() already keeps Ctrl+click out of the
        # builder_checking_press flow that would otherwise trigger that).
        if ( left and self.ctrl and not self.shift
             and getattr ( self.vm_session, "builder_atom_mode", False )
             and getattr ( self.vm_session, "builder_tool", "add" ) == "add" ):
            target_object = getattr ( self.vm_session, "builder_target_object", None )
            if target_object is not None:
                bond = click_mode.find_bond_at_pixel ( self, target_object, mouse_x, mouse_y )
                if bond is not None:
                    click_mode.apply_selected_bond_order ( self, target_object, bond )
            # limpa qualquer candidato de atom-drag que tenha ficado
            # pendente (um Ctrl+clique sem arrasto real nunca chega a
            # consumir esse candidato em mouse_motion())
            self.vm_session.builder_ctrl_press_candidate_atom  = None
            self.vm_session.builder_ctrl_press_candidate_depth = None
            # [EN] BUG FIX: same root cause as builder_checking_press just
            # below in the plain-click path (see that fix's comment for the
            # full mechanism) -- a fast Ctrl+click (press+release with no
            # render() in between) left builder_ctrl_checking_press True,
            # so the NEXT render() call would re-pick whatever's at this
            # click's position and hand it to mouse_motion() as a phantom
            # Ctrl-drag candidate, even though this click already fully
            # resolved as a bond-order-cycle with no drag involved.
            self.builder_ctrl_checking_press = False
            self.dragging = False
            self.parent_widget.queue_draw ( )
            return
        # A plain (non-Ctrl) press landed on an atom but never turned into a real drag (a
        # plain click) -- the candidate mouse_motion() would have
        # consumed on the first real movement never got claimed, so
        # clear it here rather than letting it leak into the NEXT
        # press/drag. This is exactly the "click on an atom -> replace"
        # case, handled below via the normal builder_placing_atom path.
        self.vm_session.builder_press_candidate_atom  = None
        self.vm_session.builder_press_candidate_depth = None
        # [EN] BUG FIX (reported by the user: placing a new atom left a
        # "residual" bond-drag, as if a second click had immediately tried
        # to connect two atoms). Root cause: mouse_pressed() sets
        # self.builder_checking_press = True on EVERY plain press while the
        # "add"/"move" tool is active, BEFORE it's known whether the press
        # will resolve to a drag or to a plain click -- that flag is only
        # ever CONSUMED (turned into builder_press_candidate_atom) inside
        # render(), never cleared here. For a fast click (press+release
        # with no render() in between), it was still True when this method
        # returned early via builder_placing_atom below -- so on the VERY
        # NEXT render() call (the one builder_placing_atom itself needs to
        # actually place the atom), builder_checking_press ALSO fired,
        # re-picking whatever is now at the click's pixel position -- which
        # is the atom that render() had JUST placed one line earlier (see
        # the order of those two checks in render()). That atom became a
        # dangling "candidate", and the very next mouse_motion() (no button
        # even held anymore) grabbed it and called start_bond_drag() on it,
        # exactly matching the reported symptom. Clearing the flag here,
        # alongside the candidate it would have produced, closes the race:
        # by the time a genuine drag is actually in progress, render() has
        # already had a chance to run at least once and already reset this
        # to False itself, so this is a no-op for that path (see
        # builder_bond_drag_active/builder_ctrl_drag_active's early returns
        # above, both already resolved before ever reaching this line).
        self.builder_checking_press = False
        if self.dragging:
            if left:
                if self.shift:
                    self.selection_box_picking = True
                    self.show_selection_box = False
                    self.selection_box.start = None
                    self.selection_box.end = None
                    self.parent_widget.queue_draw()
        else:
            if left:
                # [EN] Builder "click to place atom" mode -- see
                # gui/windows/builder/click_mode.py. Checked FIRST, before
                # any of the normal picking/selection state below, and
                # returns immediately if handled -- deliberately kept as
                # a single, minimal-footprint early-exit so normal
                # click-to-select behaviour (the rest of this method) is
                # completely unaffected when this mode is off (the
                # default / existing behaviour for every non-Builder use
                # of the app).
                #
                # IMPORTANT (bug fixed after a live GLError report):
                # mouse_released() is a plain GTK event handler, NOT the
                # GLArea "render" callback -- the OpenGL context is only
                # GUARANTEED current inside render() (that's exactly why
                # self.picking below is just a flag set here and actually
                # acted on, via self._pick(), from inside render() --
                # never called directly from this handler). Calling
                # click_mode.handle_click_to_place_atom() (which does real GL calls:
                # a whole extra draw pass to read the depth buffer)
                # directly, right here, raised
                # "GLError: invalid operation" on glBindVertexArray in
                # the user's real environment -- exactly the class of
                # error you get issuing GL calls with no current context.
                # Fixed the same way self.picking already works: just
                # record the click position and a flag here; the actual
                # click_mode.handle_click_to_place_atom() call now happens inside
                # render(), right next to "if self.picking: self._pick()".
                # [EN] Only intercepts for the "add" tool, and only for a
                # PLAIN click (not self.shift) -- holding shift always
                # falls through to normal click-to-select below,
                # regardless of Builder mode, so the user can shift-click
                # two atoms (the app's existing, already-accumulating
                # multi-select -- see VismolViewingSelection.
                # selecting_by_atom() in vismol_selections.py) to pick a
                # pair for the 'b' add-bond shortcut without leaving
                # Builder mode. The "delete" tool ALSO falls through here
                # on purpose: it needs _pick() to actually run (to know
                # WHICH atom was clicked) before builder_tool=='delete'
                # is acted on -- see the new check added next to
                # "if self.picking: self._pick()" in render().
                # [EN] "Add Structure" tool (vm_session.builder_tool ==
                # "add_structure") shares this EXACT same deferred-pick
                # mechanism -- placing a whole structure also needs a real
                # depth-buffer read, same GL-context constraint as placing
                # a single atom (see the comment block above). Which of the
                # two actually runs is decided in render(), by re-reading
                # builder_tool at that point -- see the hook next to
                # "if self.picking: self._pick()".
                if ( getattr ( self.vm_session, "builder_atom_mode", False )
                     and getattr ( self.vm_session, "builder_tool", "add" ) in ( "add", "add_structure" )
                     and not self.shift ):
                    self.builder_click_x = np.float32(mouse_x)
                    self.builder_click_y = np.float32(mouse_y)
                    self.builder_placing_atom = True
                    self.dragging = False
                    self.parent_widget.queue_draw ( )
                    return
                self.picking_x = np.float32(mouse_x)
                self.picking_y = np.float32(mouse_y)
                self.picking = True
                self.button = 1
                #dragging is set to false here
                self.dragging = False
                self.parent_widget.queue_draw()
            if middle:
                if self.atom_picked is not None:
                    self.dragging = True
                    self.button = 2
                    self.center_on_atom(self.atom_picked)
                    self.atom_picked = None
            if right:
                # [EN] Builder "right-click to delete" (Avogadro-style) --
                # checked FIRST, before the normal context-menu logic
                # below. Only active while builder_atom_mode + tool==
                # "add", and ONLY for vm_session.builder_target_object
                # (same "one editable object at a time" scoping as
                # Ctrl+click/bond-order-cycling). Right-click an ATOM
                # (self.atom_picked, from the normal _pick() pass that
                # already ran for this release -- see mouse_pressed(),
                # which already sets self.picking=True for the right
                # button, same as it always did) deletes that atom;
                # right-click a BOND (found via the same projection-
                # distance check Ctrl+click uses, see
                # click_mode.find_bond_at_pixel()) deletes just that
                # bond. If NEITHER is under the cursor, falls through to
                # the normal context menu below UNCHANGED -- this never
                # suppresses the menu on a background/other-object click,
                # only when something was actually deleted.
                if ( getattr ( self.vm_session, "builder_atom_mode", False )
                     and getattr ( self.vm_session, "builder_tool", "add" ) == "add" ):
                    target_object = getattr ( self.vm_session, "builder_target_object", None )
                    if target_object is not None:
                        if self.atom_picked is not None and self.atom_picked.vm_object is target_object:
                            dprint ( "DEBUG vismol_glcore: builder right-click delete -- atom #{} ('{}')".format (
                                    self.atom_picked.atom_id, self.atom_picked.symbol ) )
                            atom_ops.push_undo_snapshot ( target_object )

                            deleted_symbol = self.atom_picked.symbol
                            deleted_id     = self.atom_picked.atom_id

                            # [EN] Capture the (still-live) Atom OBJECT
                            # references of every former neighbour BEFORE
                            # removing -- atom_ops.remove_atom() renumbers every
                            # atom_id above the removed one, but mutates
                            # each SURVIVING atom's .atom_id IN PLACE, so
                            # reading it fresh off the object afterwards
                            # (below) is always correct regardless. Skips
                            # neighbours that are themselves hydrogens --
                            # nothing to "adjust" on an H (its own target
                            # valence is always satisfied by its one bond).
                            neighbor_objs = [ ]
                            for bond in target_object.bonds.values ( ):
                                if bond.atom_index_i != deleted_id and bond.atom_index_j != deleted_id:
                                    continue
                                other_id = bond.atom_index_j if bond.atom_index_i == deleted_id else bond.atom_index_i
                                neighbor = target_object.atoms[other_id]
                                if neighbor.symbol != 'H':
                                    neighbor_objs.append ( neighbor )

                            atom_ops.remove_atom ( target_object, deleted_id )

                            # [EN] Deliberately SKIPPED when the deleted
                            # atom was itself a hydrogen: adjusting its
                            # (former) parent right after would just
                            # immediately re-add a replacement H, making
                            # "delete this H" a complete no-op from the
                            # user's point of view -- clearly not what a
                            # deliberate right-click delete on an H means.
                            if deleted_symbol != 'H':
                                for neighbor in neighbor_objs:
                                    atom_ops.adjust_hydrogens ( target_object, neighbor.atom_id )

                            empty_object.sync_pdynamo_system ( target_object )

                            self.atom_picked = None
                            self.dragging = False
                            self.parent_widget.queue_draw ( )
                            return
                        elif self.atom_picked is None:
                            # so tenta achar um BOND se nenhum atomo foi
                            # pego pelo picking normal (um atomo do objeto-
                            # alvo sempre tem prioridade sobre uma ligacao)
                            bond = click_mode.find_bond_at_pixel ( self, target_object, mouse_x, mouse_y )
                            if bond is not None:
                                dprint ( "DEBUG vismol_glcore: builder right-click delete -- bond #{} <-> #{}".format (
                                        bond.atom_index_i, bond.atom_index_j ) )
                                atom_ops.push_undo_snapshot ( target_object )
                                atom_a_obj = target_object.atoms[bond.atom_index_i]
                                atom_b_obj = target_object.atoms[bond.atom_index_j]
                                atom_ops.remove_bond ( target_object, bond.atom_index_i, bond.atom_index_j )
                                # [EN] Removing a bond FREES UP valence on
                                # both sides -- unlike the atom-deletion
                                # case above, there's no "undo-my-own-
                                # click" ambiguity here (deleting a BOND,
                                # not an H atom directly), so both atoms
                                # always get adjusted, hydrogens included.
                                atom_ops.adjust_hydrogens ( target_object, atom_a_obj.atom_id )
                                atom_ops.adjust_hydrogens ( target_object, atom_b_obj.atom_id )

                                empty_object.sync_pdynamo_system ( target_object )

                                self.dragging = False
                                self.parent_widget.queue_draw ( )
                                return
                # The right button (button = 3) always opens one of the available menus.
                self.button = 3
                menu_type = None
                # Check if there is anything in the selection list
                # If {} means that there are no selection points on the screen
                # Checks if vismol_session.current_selection has any selection.
                # Also needs to check whether "picking" mode is enabled.
                if not bool(self.vm_session.selections[self.vm_session.current_selection].selected_objects) \
                   or self.vm_session.picking_selection_mode:
                    # Checks if the list of atoms selected by the picking function has any elements. 
                    # If the list is empty, the pick menu is not shown.
                    if self.vm_session.picking_selection_mode \
                       and self.vm_session.picking_selections.picking_selections_list != [None,None,None,None]:
                        info = None
                        menu_type = "pick_menu"
                    else:
                        # Here the obj_menu is activated based on the atom that
                        # was identified by the picking function.
                        # The picking function detects the selected pixel and
                        # associates it with the respective object.
                        # There is no selection (blue dots) but an atom was
                        # identified in the click with the right button
                        if self.atom_picked is not None:
                            # Getting the info about the atom that was identified in the click
                            label = "{} / {} / {}({}) / {}({} / {})".format(self.atom_picked.vm_object.name,
                                                                        self.atom_picked.chain.name,
                                                                        self.atom_picked.residue.name,
                                                                        self.atom_picked.residue.index,
                                                                        self.atom_picked.name,
                                                                        self.atom_picked.index,
                                                                        self.atom_picked.symbol)
                            self.info_atom = self.atom_picked
                            self.atom_picked = None
                            menu_type = "obj_menu"
                            info = label
                        else:
                            # When no atom is identified in the click (user clicked
                            # on a point in the background)
                            menu_type = "bg_menu"
                            info = None
                else:
                    # When a selection (viewing selection) is active, the
                    # selection menu is passed as an option.
                    info = self.vm_session.selections[self.vm_session.current_selection].get_selection_info()
                    menu_type = "sele_menu"
                # The right button (button = 3) always opens one of the available menus.
                self.parent_widget.show_gl_menu(menu_type=menu_type, info=info)
        
        
        
        self.dragging = False
        self.parent_widget.queue_draw()
    
    def mouse_motion(self, mouse_x, mouse_y, button_state=None):
        """ Function doc
        """
        # [EN] Local/deferred import -- see the note near mouse_released()
        # (and at the top of this file) for why.
        click_mode, _, _ = _builder_hooks()
        x = np.float32(mouse_x)
        y = np.float32(mouse_y)
        dx = x - self.mouse_x
        dy = y - self.mouse_y
        if (dx == 0) and (dy == 0):
            return
        self.mouse_x, self.mouse_y = x, y
        if button_state is not None:
            # [EN] macOS trackpad fix: re-derive rotate/zoom/pan from the
            # ACTUAL buttons GTK reports as held on THIS motion event,
            # every single event, rather than trusting whatever mode was
            # decided once back in mouse_pressed(). macOS can report the
            # start of a two-finger secondary-click ambiguously (see the
            # long comment in mouse_pressed()), so a mode picked once at
            # press time could get permanently stuck wrong (typically
            # rotating) for an entire gesture the user intended as a
            # zoom -- which is what read as "zoom drifting sideways".
            # Recomputing on every motion event means a wrong initial
            # guess self-corrects within a single frame, as soon as
            # macOS reports the real button combination.
            from gi.repository import Gdk
            self._macos_buttons_down = (
                (1 if button_state & Gdk.ModifierType.BUTTON1_MASK else 0) |
                (2 if button_state & Gdk.ModifierType.BUTTON2_MASK else 0) |
                (4 if button_state & Gdk.ModifierType.BUTTON3_MASK else 0))
            self.mouse_rotate, self.mouse_zoom, self.mouse_pan = \
                self._macos_modes_from_buttons(
                    bool(self._macos_buttons_down & 1),
                    bool(self._macos_buttons_down & 2),
                    bool(self._macos_buttons_down & 4))
        # [EN] Builder "click-and-drag to create a bonded atom" -- checked
        # FIRST, before the normal rotate/pan/zoom branches below.
        #
        # Two stages:
        #  1) A drag is ALREADY active (builder_bond_drag_active) --
        #     just reposition the dragged atom every motion event.
        #     click_mode.update_bond_drag() does pure math (no GL calls:
        #     world_pos_from_mouse() only touches the GPU when depth is
        #     None, and here it's always the FIXED depth captured once
        #     in click_mode.start_bond_drag()), so unlike click_mode.handle_click_to_place_atom()
        #     this is safe to call directly from a plain GTK handler,
        #     every single motion event, without deferring to render().
        #  2) No drag active yet, but render() left a CANDIDATE atom
        #     waiting (a press that landed on an atom, but we don't yet
        #     know if the user is dragging or just about to release a
        #     plain click) -- THIS motion event is the proof that real
        #     dragging is happening, so the bond-drag actually starts
        #     now (lazily). Only from this point on is mouse_rotate
        #     turned off -- a press that never turns into a real drag
        #     (a plain click) never touches mouse_rotate at all, so
        #     normal camera-rotate-by-drag keeps working exactly as
        #     before for every press that ISN'T on an atom (bug fixed
        #     after the user reported the camera not rotating AT ALL
        #     anymore while the "add" tool was active, even over empty
        #     space -- see the note in mouse_pressed()).
        if getattr ( self.vm_session, "builder_bond_drag_active", False ):
            click_mode.update_bond_drag ( self, mouse_x, mouse_y )
            self.dragging = True
            return

        candidate_atom = getattr ( self.vm_session, "builder_press_candidate_atom", None )
        if candidate_atom is not None:
            candidate_depth = getattr ( self.vm_session, "builder_press_candidate_depth", None )
            self.vm_session.builder_press_candidate_atom  = None
            self.vm_session.builder_press_candidate_depth = None

            # [EN] "move" tool -> plain drag repositions the existing
            # atom (click_mode.start_atom_drag(), the SAME functions
            # Ctrl+drag uses in the "add" tool -- see that feature's own
            # comments); any other tool ("add") -> plain drag creates a
            # new bonded atom, as before.
            if getattr ( self.vm_session, "builder_tool", "add" ) == "move":
                click_mode.start_atom_drag ( self, candidate_atom, candidate_depth )
                click_mode.update_atom_drag ( self, mouse_x, mouse_y )
            else:
                click_mode.start_bond_drag ( self, candidate_atom, candidate_depth )
                click_mode.update_bond_drag ( self, mouse_x, mouse_y )   # move it to THIS event's position right away, not just the press position
            self.mouse_rotate = False
            self.mouse_pan    = False
            self.mouse_zoom   = False
            self.dragging = True
            return

        # [EN] Builder "Ctrl+drag to reposition an existing atom" -- same
        # "already active" / "lazy start from a candidate" two-stage
        # structure as the plain bond-drag handling just above (see that
        # block's own comments for the full reasoning; identical here).
        if getattr ( self.vm_session, "builder_ctrl_drag_active", False ):
            click_mode.update_atom_drag ( self, mouse_x, mouse_y )
            self.dragging = True
            return

        ctrl_candidate_atom = getattr ( self.vm_session, "builder_ctrl_press_candidate_atom", None )
        if ctrl_candidate_atom is not None:
            ctrl_candidate_depth = getattr ( self.vm_session, "builder_ctrl_press_candidate_depth", None )
            self.vm_session.builder_ctrl_press_candidate_atom  = None
            self.vm_session.builder_ctrl_press_candidate_depth = None

            click_mode.start_atom_drag ( self, ctrl_candidate_atom, ctrl_candidate_depth )
            click_mode.update_atom_drag ( self, mouse_x, mouse_y )
            self.mouse_rotate = False
            self.mouse_pan    = False
            self.mouse_zoom   = False
            self.dragging = True
            return

        changed = False
        if self.mouse_rotate:
            changed = self._rotate_view(dx, dy, x, y)
        elif self.mouse_pan:
            changed = self._pan_view(x, y)
        elif self.mouse_zoom:
            changed = self._zoom_view(dy)
        if changed:
            self.dragging = True
            self.parent_widget.queue_draw()

        # [EN] "Hover -> print which atom + draw a highlight ring" -- used
        # to be a pure-CPU 2D-projection distance check (see click_mode.
        # find_atom_at_pixel_2d_any_object(), still defined there but no
        # longer called from here). BUG FIXED (reported by the user after
        # actually testing this): that approach ignored OCCLUSION -- it
        # only compared screen-space DISTANCE to each atom's projected
        # centre, with no idea whether some OTHER atom/geometry was
        # actually in front, blocking the view -- so it disagreed with
        # real click-picking (which reads the real depth/colour buffer,
        # respecting whatever's actually visible at that pixel) any time
        # something occluded the "closest in 2D" atom. Fixed by using the
        # EXACT SAME GPU colour-ID pick real clicks already use
        # (click_mode._read_depth_and_atom_at_pixel()) -- guaranteed
        # consistent BY CONSTRUCTION, since it's literally the same code
        # path, not a separate approximation that could disagree again in
        # some other way later.
        #
        # Real GPU picks need a render pass + a synchronous glReadPixels
        # (see that function's own docstring for why that's a genuine,
        # well-documented stall, not just "one more call") -- too
        # expensive to do on every single motion event, so THROTTLED by
        # time here (builder_hover_last_check_time) to at most ~12
        # checks/second, however fast the mouse is actually moving.
        # Recording the position and setting builder_hover_checking=True
        # is all that happens HERE, though -- the actual GL calls are
        # deferred to render() (see the hook added there), same GL-
        # context-only-guaranteed-inside-render() constraint documented
        # at length next to builder_placing_atom/builder_checking_press.
        elif not self.mouse_rotate and not self.mouse_pan and not self.mouse_zoom:
            now = time.time ( )
            last_check = getattr ( self, "builder_hover_last_check_time", 0.0 )
            if now - last_check >= 0.08:
                self.builder_hover_last_check_time = now
                self.builder_hover_check_x = np.float32 ( mouse_x )
                self.builder_hover_check_y = np.float32 ( mouse_y )
                self.builder_hover_checking = True
                self.parent_widget.queue_draw ( )
        
        #else:
        #    self.dragging = False
        #    self.parent_widget.queue_draw()
            
    
    def mouse_scroll(self, direction):
        """ Function doc
        """
        up = np.int32(direction) == 1
        down = np.int32(direction) == -1
        if self.ctrl:
            if self.editing_mols:
                for index, vm_object in self.vm_session.vm_objects_dic.items():
                    if vm_object.editing:
                        if up:
                            vm_object.model_mat = mop.my_glTranslatef(vm_object.model_mat, np.array([0.0, 0.0, -self.scroll]))
                        if down:
                            vm_object.model_mat = mop.my_glTranslatef(vm_object.model_mat, np.array([0.0, 0.0, self.scroll]))
                
                for key in self.vm_session.vm_geometric_object_dic.keys():
                    vm_object = self.vm_session.vm_geometric_object_dic[key]
                    if vm_object:
                        if vm_object.editing:
                            if up:
                                vm_object.model_mat = mop.my_glTranslatef(vm_object.model_mat, np.array([0.0, 0.0, -self.scroll]))
                            if down:
                                vm_object.model_mat = mop.my_glTranslatef(vm_object.model_mat, np.array([0.0, 0.0, self.scroll]))
            else:
                if up:
                    self.model_mat = mop.my_glTranslatef(self.model_mat, np.array([0.0, 0.0, -self.scroll]))
                if down:
                    self.model_mat = mop.my_glTranslatef(self.model_mat, np.array([0.0, 0.0, self.scroll]))
                
                for index, vm_object in self.vm_session.vm_objects_dic.items():
                    if up:
                        vm_object.model_mat = mop.my_glTranslatef(vm_object.model_mat, np.array([0.0, 0.0, -self.scroll]))
                    if down:
                        vm_object.model_mat = mop.my_glTranslatef(vm_object.model_mat, np.array([0.0, 0.0, self.scroll]))
                
                for key in self.vm_session.vm_geometric_object_dic.keys():
                    vm_object = self.vm_session.vm_geometric_object_dic[key]
                    if vm_object:
                        if up:
                            vm_object.model_mat = mop.my_glTranslatef(vm_object.model_mat, np.array([0.0, 0.0, -self.scroll]))
                        if down:
                            vm_object.model_mat = mop.my_glTranslatef(vm_object.model_mat, np.array([0.0, 0.0, self.scroll]))
        
        elif self.shift:
            #print(direction)
            #theta = direction/10
            
            # Editing Dihedral Rotamer / moving group os atoms
            if self.vm_session.picking_selection_mode:
                theta = direction/10
                atom1 = self.vm_session.picking_selections.picking_selections_list[0]
                atom2 = self.vm_session.picking_selections.picking_selections_list[1]
                atom3 = self.vm_session.picking_selections.picking_selections_list[2]
                atom4 = self.vm_session.picking_selections.picking_selections_list[3]
                
                if atom1 and atom2 and atom3 and atom4:
                                       
                    subgroup = self.vm_session.picking_selections.subgroup
                    index1   = atom2.index-1
                    index2   = atom3.index-1
                    vobject  = atom2.vm_object
                    #print(subgroup)
                    #self.vm_session.get_dihedral (theta)
                    
                    self.vm_session.rotate_dihedral ( 
                                       vobject   = vobject, 
                                       index1    = index1, 
                                       index2    = index2, 
                                       subgroup  = subgroup,
                                       theta     = theta)
                
                elif atom1 and atom2 and not atom3 and not atom4:
                    #print (atom1 ,atom2)
                    index1   = atom1.index-1
                    index2   = atom2.index-1
                    subgroup = self.vm_session.picking_selections.subgroup
                    vobject  = atom2.vm_object
                    #print(direction)
                    self.vm_session.move_subgroup  (
                                       vobject   = vobject,
                                       index1    = index1, 
                                       index2    = index2, 
                                       subgroup  = subgroup,
                                       direction = direction)
                                       
                    
                    
                    
        else:
            # [EN] 2026-09-25 slab-zoom fix (user's own report: the plain
            # scroll wheel widened/narrowed z_near/z_far by a FIXED constant
            # (self.scroll) per tick, so a step near the min_znear floor felt
            # exactly as coarse as one taken from a wide-open slab; user
            # wants finer increments specifically as z_near/z_far get close
            # together, and suggested a quadratic or exponential model).
            # Went through 2 iterations before this one, both self-caught by
            # a live repeated-tick test before shipping:
            #   1. Step proportional to the slab-MIDPOINT distance (the
            #      formula _zoom_view() already uses for its own right-drag
            #      dolly) -- wrong reference: that distance stays constant
            #      here, since every plain tick narrows/widens symmetrically
            #      about a fixed centre, so the step never actually shrank.
            #   2. Step = slab WIDTH (z_far - z_near) * a fixed fraction,
            #      i.e. w -> w*(1 +/- 2*fraction) per tick -- shrinks/grows
            #      as wanted, but is only an approximately-exponential
            #      progression: since (1+2f) and (1-2f) are not exact
            #      multiplicative inverses of each other, N ticks up
            #      followed by N ticks down does NOT return to the original
            #      width (a real, measured net drift over repeated ticks).
            # Fixed by using a genuinely exponential (multiplicative) model:
            # width *= exp(+RATE) on the way out, width *= exp(-RATE) on the
            # way in. exp(RATE) and exp(-RATE) are EXACT reciprocals of each
            # other at every scale, so this is perfectly reversible (equal
            # up/down ticks always return to the exact starting width, not
            # just approximately) while still giving finer absolute
            # increments the narrower the slab already is and coarser ones
            # the wider it is -- exactly the model the user asked for.
            # `max(..., floor)` only keeps the width used for this
            # calculation from being exactly zero/negative; the existing
            # crossing-guard below (`z_near + step + 0.005`) is what
            # actually stops further narrowing once the slab gets very thin.
            # scroll_step_fraction (a per-tick exponential RATE, not a plain
            # fraction) defaults to 0.138 -- half the per-tick step of the
            # original 0.26 default (user's own follow-up: "podemos deixar o
            # passo mais fino, talvez metade do atual" -- 0.138 is the exact
            # rate solving exp(rate)-1 = (exp(0.26)-1)/2, not a rough guess).
            slab_width = max ( self.glcamera.z_far - self.glcamera.z_near, 0.05 )
            rate = self.scroll_step_fraction
            if up:
                new_width = slab_width * np.exp ( rate )
                step = ( new_width - slab_width ) / 2.0
                self.glcamera.z_near -= step
                self.glcamera.z_far += step
            if down:
                new_width = slab_width * np.exp ( -rate )
                step = ( slab_width - new_width ) / 2.0
                if (self.glcamera.z_far-step) >= (self.glcamera.min_zfar):
                    if (self.glcamera.z_far-step) > (self.glcamera.z_near + step + 0.005):
                        self.glcamera.z_near += step
                        self.glcamera.z_far -= step

            if (self.glcamera.z_near >= self.glcamera.min_znear):
                self.glcamera.set_projection_matrix(mop.my_glPerspectivef(self.glcamera.field_of_view,
                        self.glcamera.viewport_aspect_ratio, self.glcamera.z_near, self.glcamera.z_far))
            else:
                if self.glcamera.z_far < (self.glcamera.min_zfar + self.glcamera.min_znear):
                    self.glcamera.z_near -= step
                    # [EN] 2026-09-25 bug fix: was `self.glcamera.min_clip`,
                    # an attribute that does not exist anywhere on GLCamera
                    # (only min_znear/min_zfar do) -- would raise
                    # AttributeError the first time this recovery branch was
                    # ever reached. The floor being tested two lines above
                    # (`min_zfar + min_znear`) is clearly the intended
                    # recovery value for z_far here, matching the same
                    # "if below floor: snap to floor" idiom used elsewhere
                    # in this function (e.g. the `down` branch above).
                    self.glcamera.z_far = self.glcamera.min_zfar + self.glcamera.min_znear
                self.glcamera.set_projection_matrix(mop.my_glPerspectivef(self.glcamera.field_of_view,
                                                    self.glcamera.viewport_aspect_ratio,
                                                    self.glcamera.min_znear, self.glcamera.z_far))
            self.glcamera.update_fog()
        self.parent_widget.queue_draw()

    def _rotate_view(self, dx, dy, x, y):
        """
        Rotate the scene or selected objects based on mouse movement.
        """

        sens = self.vm_config.gl_parameters['mouse_rotation_sensibility']
        dx *= sens
        dy *= sens

        angle = np.hypot(dx, dy) / max(self.width, self.height) * 180.0

        # SHIFT → selection box update
        if self.shift:
            if self.vm_session.picking_selection_mode:
                return False

            if self.selection_box.start is None:
                self.selection_box.start = self.get_viewport_pos(self.mouse_x, self.mouse_y)
            self.selection_box.end = self.get_viewport_pos(self.mouse_x, self.mouse_y)
            self.selection_box.update_points()
            return True

        # --- Compute rotation axis ---
        if self.ctrl:
            if abs(dx) >= abs(dy):
                sign = -1 if (y - self.height / 2.0) >= 0 else 1
                axis = np.array([0.0, 0.0, sign * dx], dtype=np.float32)
            else:
                sign = -1 if (x - self.width / 2.0) >= 0 else 1
                axis = np.array([0.0, 0.0, sign * dy], dtype=np.float32)
        else:
            axis = np.array([-dy, -dx, 0.0], dtype=np.float32)

        # Normalize axis
        norm = np.linalg.norm(axis)
        if norm == 0:
            return False
        axis /= norm

        rotation_matrix = mop.my_glRotatef(np.identity(4), angle, axis)

        # --- Apply rotation ---
        def iter_all_objects():
            yield from self.vm_session.vm_objects_dic.values()
            for obj in self.vm_session.vm_geometric_object_dic.values():
                if obj:
                    yield obj

        if self.editing_mols:
            for obj in iter_all_objects():
                if obj.editing:
                    obj.model_mat = mop.my_glMultiplyMatricesf(obj.model_mat, rotation_matrix)
        else:
            self.model_mat = mop.my_glMultiplyMatricesf(self.model_mat, rotation_matrix)
            for obj in iter_all_objects():
                obj.model_mat = mop.my_glMultiplyMatricesf(obj.model_mat, rotation_matrix)

        # --- Gizmo axis ---
        if not self.editing_mols:
            self.axis.model_mat = mop.my_glTranslatef(self.axis.model_mat, -self.axis.zrp)
            self.axis.model_mat = mop.my_glRotatef(self.axis.model_mat, angle, axis)
            self.axis.model_mat = mop.my_glTranslatef(self.axis.model_mat, self.axis.zrp)

        return True                 

    def _rotate_view_old(self, dx, dy, x, y):
        """ Function doc """
        dx =  dx*self.vm_config.gl_parameters['mouse_rotation_sensibility']
        dy =  dy*self.vm_config.gl_parameters['mouse_rotation_sensibility']
        angle = np.sqrt(dx**2 + dy**2) / (self.width + 1) * 180.0
        if self.shift:
            if self.vm_session.picking_selection_mode:
                return False
            else:
                self.selection_box.end = self.get_viewport_pos(self.mouse_x, self.mouse_y)
                self.selection_box.update_points()
        
        else:
            if self.ctrl:
                if abs(dx) >= abs(dy):
                    if (y - self.height / 2.0) < 0:
                        rot_mat = mop.my_glRotatef(np.identity(4), angle, np.array([0.0, 0.0, dx]))
                    else:
                        rot_mat = mop.my_glRotatef(np.identity(4), angle, np.array([0.0, 0.0, -dx]))
                else:
                    if (x - self.width / 2.0) < 0:
                        rot_mat = mop.my_glRotatef(np.identity(4), angle, np.array([0.0, 0.0, -dy]))
                    else:
                        rot_mat = mop.my_glRotatef(np.identity(4), angle, np.array([0.0, 0.0, dy]))
            else:
                rot_mat = mop.my_glRotatef(np.identity(4), angle, np.array([-dy, -dx, 0.0]))
            
            if self.editing_mols:
                for index, vm_object in self.vm_session.vm_objects_dic.items():
                    if vm_object.editing:
                        vm_object.model_mat = mop.my_glMultiplyMatricesf(vm_object.model_mat, rot_mat)
                
                for key in self.vm_session.vm_geometric_object_dic.keys():
                    vm_object = self.vm_session.vm_geometric_object_dic[key]
                    if vm_object:
                        if vm_object.editing:
                            vm_object.model_mat = mop.my_glMultiplyMatricesf(vm_object.model_mat, rot_mat)
            else:
                self.model_mat = mop.my_glMultiplyMatricesf(self.model_mat, rot_mat)
                for index, vm_object in self.vm_session.vm_objects_dic.items():
                    vm_object.model_mat = mop.my_glMultiplyMatricesf(vm_object.model_mat, rot_mat)
                
                for key in self.vm_session.vm_geometric_object_dic.keys():
                    vm_object = self.vm_session.vm_geometric_object_dic[key]
                    if vm_object:
                        vm_object.model_mat = mop.my_glMultiplyMatricesf(vm_object.model_mat, rot_mat)
            
            # Axis operations, this code only affects the gizmo axis
            if not self.editing_mols:
                self.axis.model_mat = mop.my_glTranslatef(self.axis.model_mat, -self.axis.zrp)
                if self.ctrl:
                    if abs(dx) >= abs(dy):
                        if (y - self.height / 2.0) < 0:
                            self.axis.model_mat = mop.my_glRotatef(self.axis.model_mat, angle, np.array([0.0, 0.0, dx]))
                        else:
                            self.axis.model_mat = mop.my_glRotatef(self.axis.model_mat, angle, np.array([0.0, 0.0, -dx]))
                    else:
                        if (x - self.width / 2.0) < 0:
                            self.axis.model_mat = mop.my_glRotatef(self.axis.model_mat, angle, np.array([0.0, 0.0, -dy]))
                        else:
                            self.axis.model_mat = mop.my_glRotatef(self.axis.model_mat, angle, np.array([0.0, 0.0, dy]))
                else:
                    self.axis.model_mat = mop.my_glRotatef(self.axis.model_mat, angle, np.array([dy, dx, 0.0]))
                self.axis.model_mat = mop.my_glTranslatef(self.axis.model_mat, self.axis.zrp)
                #self.axis.model_mat = mop.my_glScalef(self.axis.model_mat, self.axis.zrp)
                
            # Axis operations, this code only affects the gizmo axis
        return True
    
    def _pan_view(self, x, y):
        """ Function doc """
        px, py, pz = self._mouse_pos(x, y)
        pan_mat = mop.my_glTranslatef(np.identity(4, dtype=np.float32), np.array(
                                    [(px - self.drag_pos_x) * self.glcamera.z_far / 10.0,
                                     (py - self.drag_pos_y) * self.glcamera.z_far / 10.0,
                                     (pz - self.drag_pos_z) * self.glcamera.z_far / 10.0]))
        if self.editing_mols:
            for index, vm_object in self.vm_session.vm_objects_dic.items():
                if vm_object.editing:
                    vm_object.model_mat = mop.my_glMultiplyMatricesf(vm_object.model_mat, pan_mat)
            
            for key in self.vm_session.vm_geometric_object_dic.keys():
                vm_object = self.vm_session.vm_geometric_object_dic[key]
                if vm_object:
                    if vm_object.editing:
                        vm_object.model_mat = mop.my_glMultiplyMatricesf(vm_object.model_mat, pan_mat)
        
        else:
            self.model_mat = mop.my_glMultiplyMatricesf(self.model_mat, pan_mat)
            for index, vm_object in self.vm_session.vm_objects_dic.items():
                vm_object.model_mat = mop.my_glMultiplyMatricesf(vm_object.model_mat, pan_mat)
            
            for key in self.vm_session.vm_geometric_object_dic.keys():
                vm_object = self.vm_session.vm_geometric_object_dic[key]
                if vm_object:
                    vm_object.model_mat = mop.my_glMultiplyMatricesf(vm_object.model_mat, pan_mat)
            self.zero_reference_point = mop.get_xyz_coords(self.model_mat)
        
        self.drag_pos_x = px
        self.drag_pos_y = py
        self.drag_pos_z = pz
        return True

    def _zoom_view(self, dy):
        """ Function doc """
        delta = (((self.glcamera.z_far - self.glcamera.z_near) / 2.0) + self.glcamera.z_near) / 200.0
        move_z = dy * delta
        # Dolly the camera along its own Z axis by writing straight into
        # the view matrix's translation row -- mathematically identical to
        # the previous my_glTranslatef(view_matrix, [0,0,move_z]) call
        # (verified: both produce byte-identical matrices), just cheaper
        # than a full 4x4 multiply for a single-component update. This is
        # NOT a fix for sideways zoom drift -- that bug is not in this
        # matrix update; see the macOS gesture-detection code in
        # mouse_pressed()/mouse_motion() instead.
        moved_mat = np.array(self.glcamera.view_matrix, dtype=np.float32, copy=True)
        moved_mat[3, 2] += move_z
        moved_pos = mop.get_xyz_coords(moved_mat)
        if moved_pos[2] > 0.101:
            new_z_near = self.glcamera.z_near - move_z
            new_z_far = self.glcamera.z_far - move_z

            if new_z_far < (self.glcamera.min_zfar + self.glcamera.min_znear):
                return False
            proj_z_near = max(new_z_near, self.glcamera.min_znear)
            if new_z_far <= (proj_z_near + 0.005):
                return False

            self.glcamera.set_view_matrix(moved_mat)
            self.glcamera.z_near = new_z_near
            self.glcamera.z_far = new_z_far
            self.glcamera.set_projection_matrix(mop.my_glPerspectivef(self.glcamera.field_of_view,
                                                self.glcamera.viewport_aspect_ratio,
                                                proj_z_near, self.glcamera.z_far))
            self.glcamera.update_fog()
            self.dist_cam_zrp += -move_z
            return True
        return False
    
    def _print_hover_atom_info_idle(self, atom):
        """ One-shot idle callback: prints the hovered atom's details, if
        it is still the hovered atom (skips stale ones when the cursor kept
        moving). """
        if atom is getattr(self, "builder_hover_atom", None):
            click_mode, _, _ = _builder_hooks()
            click_mode.print_hover_atom_info(self, atom)
        return False

    def _hover_atom_callback_idle(self, atom):
        """ One-shot idle callback: hands the hovered atom (or None) to
        vm_session.hover_atom_callback, skipping stale ones. """
        callback = getattr(self.vm_session, "hover_atom_callback", None)
        if callback is not None and atom is getattr(self, "builder_hover_atom", None):
            try:
                callback(atom)
            except Exception as error:
                dprint("hover_atom_callback failed: {}".format(error))
        return False

    def render(self, rotating = False):
        """ This is the function that will be called everytime the window
            needs to be re-drawed.
        """
        # [EN] Local/deferred import -- see the note near mouse_released()
        # (and at the top of this file) for why.
        click_mode, _, _ = _builder_hooks()
        # Medidor opcional: marca o inicio do corpo do render. Mede o custo
        # de CPU+submissao GL deste frame (nao inclui o tempo ocioso entre
        # frames), que e exatamente o numero relevante para avaliar se o
        # refactor de UBO (Gargalo 2) compensa.
        if self.show_fps:
            _t_start = time.perf_counter()
        if self.shader_flag:
            self.create_gl_programs()
            self.selection_box.initialize_gl()
            self.axis.initialize_gl()
            self.shader_flag = False
        # [EN] BUG FIXED (reported by the user after testing frame
        # navigation): this coordinate-dirty-flagging block used to sit
        # much further down in render(), AFTER self._selection_box_pick(),
        # self._pick(), AND every Builder deferred picking check
        # (builder_checking_press/builder_ctrl_checking_press/
        # builder_hover_checking) -- all of which call
        # draw_background_sel_representation() themselves, which only
        # refreshes a representation's sel_coord_vbo when THIS block has
        # already set was_sel_coord_modified=True for it. On the exact
        # render() call where a trajectory frame just advanced
        # (vm_session.forward_frame()/reverse_frame() setting
        # self.updated_coords=True), every one of those picking passes
        # would run BEFORE this block ever got a chance to mark anything
        # dirty -- reading whatever STALE positions were left over from
        # the PREVIOUS frame instead. Moved to the very top of render(),
        # before any picking pass whatsoever, so a frame change is always
        # fully reflected before anything tries to pick against it, no
        # matter which of those picking paths happens to also be pending
        # on that same render() call. Same reasoning _pick() already
        # applies to itself for the camera UBO specifically (see its own
        # comment: "O picking roda ANTES do update_camera_ubo() do render()
        # principal... logo atualizamos agora") -- this fixes the
        # equivalent problem for coordinates.
        if self.updated_coords:
            for vm_object in self.vm_session.vm_objects_dic.values():
                if vm_object.core_representations["picking_dots"] is None:
                    vm_object.build_core_representations()
                vm_object.core_representations["picking_dots"].was_rep_coord_modified = True
                for rep in vm_object.representations.values():
                    if rep is not None:
                        # Only flag coordinates as dirty for representations
                        # that are actually drawn. Inactive representations
                        # (e.g. a hidden "spheres" rep) would otherwise have
                        # their full coordinate VBO re-uploaded to the GPU
                        # every single frame of trajectory/MD playback for
                        # nothing, since draw_representation() never runs
                        # for them while inactive. was_sel_coord_modified is
                        # safe to skip here too: draw_background_sel_representation()
                        # (which consumes that flag) is itself only called
                        # for representations where rep.active is True - see
                        # _selection_box_pick(). Re-activating a representation
                        # already re-marks its indexes as modified elsewhere
                        # (vismol_session.py), and we mark coords there too
                        # below so a reactivated rep never draws stale coords.
                        if rep.active:
                            rep.was_rep_coord_modified = True
                            rep.was_sel_coord_modified = True
                        if rep.is_dynamic:
                            rep.was_rep_ind_modified = True
                            rep.was_sel_ind_modified = True
            # [EN] 2026-09-29 PERFORMANCE BUG FIX (user report: "a
            # representacao de cartoon esta reduzindo muito a performance"):
            # this flag was never reset here -- EasyHybridSession.set_frame()
            # (called on load, by trajectory windows, ...) sets it and it
            # stayed True FOREVER, so every later redraw (each mouse event
            # while rotating) re-marked every active representation dirty.
            # Atom-based reps "only" re-uploaded their coordinate VBOs every
            # frame; Cartoon rebuilds its whole mesh on that flag (measured
            # on 6VXX: 21 full rebuilds in 21 redraws, 254 ms/frame). It is
            # a one-shot "the frame changed" signal: consumed right here.
            # (The Builder's live atom drags mark their own reps directly --
            # see atom_ops.move_atom() -- they don't rely on it staying set.)
            self.updated_coords = False
        if self.selection_box_picking:
            self._selection_box_pick()
        if self.picking:
            self._pick()
        # [EN] Builder "click to place atom" mode -- see the note in
        # mouse_released() for why this can't be called from there
        # directly (no guaranteed-current GL context outside render()).
        # Mirrors the self.picking / self._pick() pattern immediately
        # above: mouse_released() only sets a flag + the click
        # coordinates; the actual GL work (click_mode.handle_click_to_place_atom(),
        # which reads the depth buffer via a real draw pass) only runs
        # here, where GTK guarantees the context is current.
        if getattr ( self, "builder_placing_atom", False ):
            # [EN] "Add Structure" tool re-uses this exact deferred flag
            # (see the mouse_released() hook that sets it) -- branch on
            # builder_tool here, the one place both paths actually run the
            # real GL depth-buffer read.
            if getattr ( self.vm_session, "builder_tool", "add" ) == "add_structure":
                msg = click_mode.handle_click_to_add_structure ( self, self.builder_click_x, self.builder_click_y )
                if msg is not None:
                    dprint ( "DEBUG click_mode: {}".format ( msg ) )
            else:
                click_mode.handle_click_to_place_atom ( self, self.builder_click_x, self.builder_click_y )
            self.builder_placing_atom = False
        # [EN] Builder "click-and-drag to create a bonded atom" -- see
        # the note added to mouse_pressed() for why this has to be
        # resolved here (GL context only guaranteed current inside
        # render(), same reasoning as builder_placing_atom right above).
        # Only records a CANDIDATE atom (on vm_session, consumed by
        # mouse_motion()) if the press actually landed on an atom
        # belonging to the CURRENT Builder target object -- pressing on
        # empty space, or on an atom from a different object, leaves the
        # candidate at None, and the interaction falls back to whatever
        # mouse_released() ends up doing with that press/release pair.
        # Deliberately does NOT start the drag itself here (that used to
        # be a bug: it made EVERY press-on-an-atom immediately create a
        # new bonded atom, even a plain click with no real drag,
        # breaking the "click on an atom -> replace its element"
        # interaction, which depends on NO drag having happened) -- see
        # mouse_motion() for where the candidate actually turns into a
        # live drag, only once genuine mouse movement confirms the user
        # is dragging, not just clicking.
        if getattr ( self, "builder_checking_press", False ):
            self.builder_checking_press = False
            if ( getattr ( self.vm_session, "builder_atom_mode", False )
                 and getattr ( self.vm_session, "builder_tool", "add" ) in ( "add", "move" ) ):
                pressed_atom, depth = click_mode._read_depth_and_atom_at_pixel ( self, self.builder_press_x, self.builder_press_y )
                target_object = getattr ( self.vm_session, "builder_target_object", None )
                if pressed_atom is not None and depth is not None and pressed_atom.vm_object is target_object:
                    self.vm_session.builder_press_candidate_atom  = pressed_atom
                    self.vm_session.builder_press_candidate_depth = depth
        # [EN] Builder "Ctrl+drag to reposition an existing atom" -- same
        # deferred-pick reasoning as builder_checking_press just above,
        # just for the Ctrl-modified press instead of the plain one.
        if getattr ( self, "builder_ctrl_checking_press", False ):
            self.builder_ctrl_checking_press = False
            if ( getattr ( self.vm_session, "builder_atom_mode", False )
                 and getattr ( self.vm_session, "builder_tool", "add" ) == "add" ):
                pressed_atom, depth = click_mode._read_depth_and_atom_at_pixel ( self, self.builder_ctrl_press_x, self.builder_ctrl_press_y )
                target_object = getattr ( self.vm_session, "builder_target_object", None )
                if pressed_atom is not None and depth is not None and pressed_atom.vm_object is target_object:
                    self.vm_session.builder_ctrl_press_candidate_atom  = pressed_atom
                    self.vm_session.builder_ctrl_press_candidate_depth = depth
        # [EN] "Hover -> print which atom + draw a highlight ring" --
        # resolves the throttled request from mouse_motion() (see that
        # hook's own comment for the full reasoning: a real GPU pick,
        # deferred here for the GL-context reason, throttled by time
        # there so this doesn't run on every motion event). Works ANY
        # TIME (not gated on builder_atom_mode), over EVERY active object
        # in the session -- _read_depth_and_atom_at_pixel() already
        # iterates vm_session.vm_objects_dic.values() itself, so no
        # per-object loop is needed here.
        if getattr ( self, "builder_hover_checking", False ):
            self.builder_hover_checking = False
            hovered_atom, _hover_depth = click_mode._read_depth_and_atom_at_pixel (
                    self, self.builder_hover_check_x, self.builder_hover_check_y )
            if hovered_atom is not getattr ( self, "builder_hover_atom", None ):
                self.builder_hover_atom = hovered_atom
                if hovered_atom is not None:
                    dprint ( "DEBUG click_mode: hovering atom #{} ('{}') of object '{}'".format (
                            hovered_atom.atom_id, hovered_atom.symbol, hovered_atom.vm_object.name ) )
                    # [EN] 2026-10-03: atom details -> terminal (system, chain,
                    # residue, element, name, DYFF type, FF type). Deferred to
                    # idle so a (cached, but sometimes recomputed) DYFF typing
                    # never stalls this frame.
                    # 2026-10-04: disabled -- superseded by the main window's
                    # hover statusbar (hover_atom_callback below).
                    #from gi.repository import GLib
                    #GLib.idle_add ( self._print_hover_atom_info_idle, hovered_atom )
                # [EN] 2026-10-04: GUI hook (main window hover statusbar) --
                # called on EVERY change, including None (cursor moved onto
                # the background), so the GUI can clear its readout.
                if getattr ( self.vm_session, "hover_atom_callback", None ) is not None:
                    from gi.repository import GLib
                    GLib.idle_add ( self._hover_atom_callback_idle, hovered_atom )
        # [EN] Builder "delete atom" tool -- unlike "add" (handled
        # entirely above via builder_placing_atom, no atom identification
        # needed), "delete" needs to know WHICH atom was clicked, so it
        # deliberately let the click fall through to normal self.picking
        # / self._pick() above (see the comment in mouse_released()) --
        # self.atom_picked is already set correctly by the time we get
        # here. Only acts when the tool is actually 'delete', so a normal
        # (non-Builder) click-to-select still works exactly as before.
        if ( getattr ( self.vm_session, "builder_atom_mode", False )
             and getattr ( self.vm_session, "builder_tool", "add" ) == "delete"
             and self.atom_picked is not None ):
            # [EN] 2026-09-28: an exception here used to abort render()
            # right after the picking pass, so the screen showed the
            # picking-ID colours instead of the scene. Report it and keep
            # drawing the normal frame.
            try:
                click_mode.handle_click_to_delete_atom ( self )
            except Exception:
                import traceback
                traceback.print_exc ( )
                self.atom_picked = None

        # [EN] Builder "Bond Order" tool -- same "let the click fall
        # through to normal picking first" reasoning as "delete" just
        # above (click_mode.handle_click_to_set_bond_order() needs to
        # know WHICH atom was clicked). See that function's own docstring
        # for the full two-click "arm first atom, click a second to
        # apply" interaction -- added after the user found both Ctrl+
        # click-on-the-bond-line and the pk1/pk2 button too fiddly for
        # changing an EXISTING bond's order.
        if ( getattr ( self.vm_session, "builder_atom_mode", False )
             and getattr ( self.vm_session, "builder_tool", "add" ) == "bond_order"
             and self.atom_picked is not None ):
            msg = click_mode.handle_click_to_set_bond_order ( self )
            if msg is not None:
                dprint ( "DEBUG click_mode: {}".format ( msg ) )

        # [EN] Builder "Attach Fragment" tool -- same "let the click fall
        # through to normal picking first" wiring as "delete"/"Bond
        # Order" above. See click_mode.handle_click_to_attach_fragment()
        # for the full behaviour (click a hydrogen atom to replace it
        # with whatever fragment is currently selected in the Fragment
        # Library window).
        if ( getattr ( self.vm_session, "builder_atom_mode", False )
             and getattr ( self.vm_session, "builder_tool", "add" ) == "attach_fragment"
             and self.atom_picked is not None ):
            msg = click_mode.handle_click_to_attach_fragment ( self )
            if msg is not None:
                dprint ( "DEBUG click_mode: {}".format ( msg ) )

        # [EN] Builder "Rotate Dihedral" tool -- same "let the click fall
        # through to normal picking first" wiring as "delete"/"Bond
        # Order"/"Attach Fragment" above. See click_mode.handle_click_to_
        # pick_dihedral_atom() for the full behaviour (click 4 atoms in
        # sequence, pk1->pk4, to arm the sidebar's rotation slider).
        if ( getattr ( self.vm_session, "builder_atom_mode", False )
             and getattr ( self.vm_session, "builder_tool", "add" ) == "rotate_dihedral"
             and self.atom_picked is not None ):
            msg = click_mode.handle_click_to_pick_dihedral_atom ( self )
            if msg is not None:
                dprint ( "DEBUG click_mode: {}".format ( msg ) )

        #print('self.dragging', self.dragging)
        
        # from here on only the visible scene is drawn (every picking pass
        # above read its ID colours back from the single-sample framebuffer)
        self._msaa_begin()
        GL.glClearColor(self.bckgrnd_color[0], self.bckgrnd_color[1],
                        self.bckgrnd_color[2], self.bckgrnd_color[3])
        GL.glClear(GL.GL_COLOR_BUFFER_BIT | GL.GL_DEPTH_BUFFER_BIT)
        # Atualiza view/proj no UBO uma unica vez por frame (Gargalo 2).
        # Todos os shaders convertidos leem deste buffer compartilhado.
        self.update_camera_ubo()
        
        '''
                           - - -  R E P R E S E N T A T I O N S - - -  
        '''
        for vm_object in self.vm_session.vm_objects_dic.values():
            if vm_object.active:
                if vm_object.frames.shape[0] > 0:
                    #print(vm_object.representations)
                    for representation in vm_object.representations.values():
                        if representation is not None:
                            # Only shows the representation if
                            # representations[rep_name].active = True
                            if representation.active:
                                representation.draw_representation()
                                
        
        if self.vm_session.picking_selection_mode:
            '''
                               - - -  P I C K I N G   S E L E C T I O N S - - -  
                                Here is where we will draw the dashed lines
                 Drawing labels on the fly is not very efficient; however, since there are only a few labels 
                 to be generated, the impact on performance is minimal.
            '''
            #'''
            '''Two loops are necessary so that the labels are always drawn above the dash lines.'''
            for geo_object in ["pk1pk2", "pk2pk3", "pk3pk4", "pk_angle", "pk_dihedral"]:
                if self.vm_session.vm_geometric_object_dic.get(geo_object):
                    if self.vm_session.vm_geometric_object_dic[geo_object].representations["dash"].active:
                        self.vm_session.vm_geometric_object_dic[geo_object].representations["dash"].draw_representation()
            
            '''                              drawing distance labels                     '''
            for geo_object in ["pk1pk2", "pk2pk3", "pk3pk4"]:
                if self.vm_session.vm_geometric_object_dic[geo_object]:
                    if self.vm_session.vm_geometric_object_dic[geo_object].representations["dash"].active:
                        self._draw_distance_labels(self.vm_session.vm_geometric_object_dic[geo_object])
            for geo_object in ["pk_angle", "pk_dihedral"]:
                arc = self.vm_session.vm_geometric_object_dic.get(geo_object)
                if arc and arc.representations["dash"].active:
                    self._draw_arc_labels(arc)
            
            '''                               drawing picking spheres                     '''
            for geo_object in ["pk1", "pk2", "pk3", "pk4"]:
                if self.vm_session.vm_geometric_object_dic[geo_object]:
                    if self.vm_session.vm_geometric_object_dic[geo_object].representations["picking_spheres"].active:
                        self.vm_session.vm_geometric_object_dic[geo_object].representations["picking_spheres"].draw_representation()
                        #print(geo_object,self.vm_session.vm_geometric_object_dic[geo_object].representations["picking_spheres"].active)
            '''                               drawing picking labels                     '''
            self._draw_picking_label()
            #------------------------------------------------------------------
            #'''
        
        
        elif getattr(self.vm_session, "builder_atom_mode", False):
            # [EN] 2026-09-25, user's own explicit request: "enquanto o
            # modo 'editing' estiver 'on', nao mostrar nenhuma selecao"
            # -- suppresses the viewing-selection dots (and the drag-box
            # preview below) while the Builder's own "Editing: ON/OFF"
            # toggle is on, EVEN IF a selection was already made earlier
            # while paused (Editing off) -- this is a pure RENDERING
            # suppression; the underlying selected_atoms set is left
            # untouched, so it reappears exactly as it was the moment
            # Editing is switched back off (see vismol_session.py's own
            # _selection_function_set(), which separately blocks NEW
            # selection changes from being made at all while this is on).
            pass
        else:
            '''                    - - -  V I E W I N G   S E L E C T I O N S - - -                  '''
            '''                               drawing blues dot selections                           '''
            for vm_object in self.vm_session.selections[self.vm_session.current_selection].selected_objects:
                # Here are represented the blue dots referring to the atom's selections
                if vm_object.core_representations["picking_dots"] is None:
                    vm_object.build_core_representations()
                pdots = vm_object.core_representations["picking_dots"]
                # Gargalo 1: a reconstrucao da lista de indices + re-upload do
                # VBO so precisa acontecer quando o conjunto de atomos
                # selecionados muda. Antes isso era feito INCONDICIONALMENTE
                # a cada frame (rebuild de lista Python, np.array novo e
                # glBufferData), desperdicando trabalho durante rotacao/zoom
                # com a selecao parada. Comparamos com um snapshot do ultimo
                # conjunto subido; so re-subimos no diff. O draw continua
                # todo frame, pois a cena e redesenhada normalmente.
                sel_ids = vm_object.selected_atom_ids
                if pdots._last_uploaded_sel_ids != sel_ids:
                    pdots.was_rep_ind_modified = True
                    pdots.define_new_indexes_to_vbo(list(sel_ids))
                    # snapshot por copia: selected_atom_ids e mutado in-place
                    # (add/discard/clear), entao guardar a referencia nao
                    # detectaria mudancas. set(...) congela o estado atual.
                    pdots._last_uploaded_sel_ids = set(sel_ids)
                pdots.draw_representation()
        
        if ( not self.vm_session.picking_selection_mode and self.show_selection_box and self.shift
             and not getattr(self.vm_session, "builder_atom_mode", False) ):
            if self.selection_box.vao is None:
                self.selection_box._make_gl_selection_box()
            else:
                self.selection_box._draw()
        
        if self.show_axis:
            self.axis._draw(True)
            self.axis._draw(False)

        # [EN] Builder "hover -> highlight ring" -- see click_mode.
        # draw_hover_highlight()'s own docstring. Drawn last (after every
        # other representation, axis, selection box) so it isn't
        # accidentally covered by anything -- it still respects normal
        # depth testing against the actual scene geometry (GL_DEPTH_TEST
        # stays enabled inside draw_hover_highlight() itself), it's just
        # not competing with any OTHER overlay for draw order here.
        # [EN] 2026-10-17: atoms deleted/renumbered since they were stored
        # (Delete by molecule, Undo...) must not reach the highlights below.
        click_mode.drop_stale_builder_atoms(self)
        hovered_atom = getattr(self, "builder_hover_atom", None)
        if hovered_atom is not None:
            world_center, hover_up, hover_radius = click_mode.draw_hover_highlight(self, hovered_atom)
            click_mode.draw_hover_info_text(self, hovered_atom, world_center, hover_up, hover_radius)

        # [EN] Builder "drag-to-bond -- live preview" (roadmap section 4,
        # plus the "light up which EXISTING atom this will bond to" the
        # user asked for right after seeing section 4's ring in action):
        # reuses the SAME hover-ring drawing code above (click_mode.
        # draw_hover_highlight(), already proven live/working here),
        # colour-coded by the bond's CURRENT effective order/aromaticity
        # (click_mode.update_bond_drag() keeps builder_bond_drag_preview_
        # order/_preview_aromatic in sync every motion event -- ALWAYS,
        # whether that order comes from the live drag distance or from an
        # explicit sidebar pick, see that function's own comment for why
        # both cases need a colour here, not just the distance-decided
        # one).
        #
        # WHICH atom gets the ring: builder_bond_drag_snap_target (also
        # kept in sync by update_bond_drag(), via the SAME _current_snap_
        # target() finish_bond_drag() itself uses, so this never lights
        # up an atom that release would then not actually snap onto) --
        # the EXISTING atom the drag is currently close enough to connect
        # to, if any; otherwise the ring stays on the dragged (temporary)
        # atom itself, exactly as before this atom-highlighting addition.
        #
        # Deliberately does NOT reuse draw_hover_info_text() (that one is
        # dead code -- returns False immediately, see its own docstring)
        # for the NUMERIC LABEL drawn on the ring below -- draw_ring_label()
        # is a new, small, CORRECT replacement (see its own docstring for
        # exactly why draw_hover_info_text()'s old approach can't just be
        # revived as-is).
        drag_preview_order = getattr(self.vm_session, "builder_bond_drag_preview_order", None)
        if drag_preview_order is not None:
            snap_target  = getattr(self.vm_session, "builder_bond_drag_snap_target", None)
            dragged_atom = getattr(self.vm_session, "builder_bond_drag_new_atom", None)
            highlight_atom = snap_target if snap_target is not None else dragged_atom
            if highlight_atom is not None:
                preview_color = click_mode.bond_order_preview_color(
                        drag_preview_order, getattr(self.vm_session, "builder_bond_drag_preview_aromatic", False))
                world_center, _hl_up, _hl_radius = click_mode.draw_hover_highlight(self, highlight_atom, color=preview_color, alpha=0.35)
                # [EN] Numeric label (user's own request, after seeing the
                # coloured ring: "quero colocar label nos circulos... '1',
                # '2', '3'") -- drawn dead-centred on the SAME ring, in
                # plain white so it stays readable against any of the 3
                # order colours (and the distinct aromatic purple) alike.
                click_mode.draw_ring_label(self, world_center, str(drag_preview_order))

        # [EN] Builder "Bond Order" tool -- persistent highlight on the
        # FIRST atom armed by click_mode.handle_click_to_set_bond_order()
        # (see its own docstring), waiting for a second click to complete
        # the pair. Same ring+label mechanism as the drag preview just
        # above, coloured/labelled by the sidebar's CURRENT selection
        # (there's no distance involved here to override it -- this tool
        # always uses exactly what the sidebar says, unlike drag-to-bond).
        armed_atom = getattr(self.vm_session, "builder_bond_pick_first_atom", None)
        if armed_atom is not None:
            armed_order    = getattr(self.vm_session, "builder_bond_order", 1)
            armed_aromatic = getattr(self.vm_session, "builder_bond_aromatic", False)
            armed_color = click_mode.bond_order_preview_color(armed_order, armed_aromatic)
            world_center, _hl_up, _hl_radius = click_mode.draw_hover_highlight(self, armed_atom, color=armed_color, alpha=0.35)
            click_mode.draw_ring_label(self, world_center, str(armed_order))

        # [EN] Builder "Rotate Dihedral" tool -- highlights each of the up
        # to 4 atoms picked so far (click_mode.handle_click_to_pick_
        # dihedral_atom()), numbered "1".."4" in click order so it's
        # visually obvious which atom is which once the dihedral is armed
        # (atoms #2/#3 are the rotation axis, #1/#4 are only there to
        # define/measure the angle). Same ring+label mechanism as the
        # other tools above; a single neutral color (not order-coded --
        # there's no bond order concept here) keeps it visually distinct
        # from the Bond Order tool's own rings.
        dihedral_pick_atoms = getattr(self.vm_session, "builder_dihedral_pick_atoms", None)
        if dihedral_pick_atoms:
            for pick_index, picked_atom in enumerate(dihedral_pick_atoms):
                world_center, _hl_up, _hl_radius = click_mode.draw_hover_highlight(
                        self, picked_atom, color=(0.2, 0.6, 1.0), alpha=0.35)
                click_mode.draw_ring_label(self, world_center, str(pick_index + 1))

        # [EN] Builder "Atom Types" window -- user's own follow-up
        # request ("colocar o label dos indexes no vobject, [para]
        # serem facilmente identificados"): while that window is open,
        # numbers each atom it's currently listing with its own atom_id
        # (matching column 0 of that window's own treeview exactly), so
        # a row there and the atom it refers to on screen are trivially
        # easy to match up. Text-ONLY (no highlight ring, unlike the
        # tools above -- user's own follow-up request, see the comment
        # right at the draw call below for why);
        # vm_session.builder_atom_types_labeled_atoms is set/cleared by
        # atom_types_window.py itself (populated on open/refresh, wiped
        # on close -- also wiped by builder_sidebar.py's own
        # close_window() when the WHOLE Builder session ends).
        #
        # Guards against a STALE reference (an atom deleted by some
        # OTHER Builder action since this window's own last refresh, not
        # yet caught up) the same way the Bond Order tool's own armed-
        # atom guard already does -- drawing a highlight for an atom_id
        # that no longer exists in its object would crash the same way
        # a stale representation's own indexes already did before that
        # was fixed (project_edit_builder_stale_representation_crash).
        #
        # BUG FIX (caught by this session's own test, not a live report):
        # checking ONLY "does this atom_id still exist" is NOT enough --
        # remove_atom() RENUMBERS every surviving atom above the removed
        # one DOWN by one, mutating .atom_id IN PLACE on the SAME Atom
        # objects (see remove_atom()'s own docstring); the just-deleted
        # atom's OWN .atom_id is never touched, so it can end up
        # numerically EQUAL to whatever atom_id a surviving atom just got
        # renumbered into. Comparing by IDENTITY (is this atom_id slot
        # still pointing back at THIS SAME object) is what actually
        # detects staleness correctly.
        #
        # [EN] BUG FIX/DESIGN CHANGE (user's own explicit follow-up --
        # "nao precisamos das esferas amarelas de realce ... vamos
        # deixar apenas os labels"): no longer calls draw_hover_
        # highlight() here at all -- just the number, positioned via the
        # new, highlight-free click_mode.atom_world_position() instead.
        #
        # [EN] 2026-10-04 (user's request: the index labels were too large
        # while editing parameters): drawn like the picking labels pk1-pk4
        # (_draw_text_labels, atom.coords(frame)) with self.vm_font_index,
        # so they follow the Preferences label font/size and ZOOM
        # SENSITIVITY instead of draw_ring_label()'s constant on-screen size
        # (kept for the bond-order ring, where a constant size is wanted).
        labeled_atoms = getattr(self.vm_session, "builder_atom_types_labeled_atoms", None)
        if labeled_atoms:
            entries = []
            for labeled_atom in labeled_atoms:
                if labeled_atom.vm_object.atoms.get(labeled_atom.atom_id) is not labeled_atom:
                    continue
                text = str(labeled_atom.atom_id)
                frame = self._get_vismol_object_frame(labeled_atom.vm_object)
                # centred on the atom: glyph centres at i * spacing, i = 0..n-1
                entries.append((text, labeled_atom.coords(frame), -(len(text) - 1) * 0.65 / 2.0))
            if entries:
                self._draw_text_labels(self.vm_font_index, entries, char_spacing=0.65)

        self._msaa_end()
        # Fecha o cronometro do frame e reporta periodicamente. Medimos so o
        # tempo de CPU/submissao (sem glFinish, que serializaria a GPU e
        # distorceria o numero). 'ms/frame' aqui = custo da thread Python por
        # render; e o indicador para decidir sobre o UBO. O 'FPS' derivado e
        # o teto teorico se o render fosse o unico limite.
        if self.show_fps:
            dt = time.perf_counter() - _t_start
            self._fps_accum_time += dt
            self._fps_frame_count += 1
            if self._fps_frame_count >= self._fps_report_every:
                avg_s = self._fps_accum_time / self._fps_frame_count
                avg_ms = avg_s * 1000.0
                fps = (1.0 / avg_s) if avg_s > 0 else float("inf")
                n_obj = len(self.vm_session.vm_objects_dic)
                dprint("[FPS] {:.1f} fps | {:.3f} ms/render | {} objects | "
                      "average of {} frames".format(fps, avg_ms, n_obj,
                                                   self._fps_frame_count))
                self._fps_frame_count = 0
                self._fps_accum_time = 0.0
        return True
    
    def antialias_samples(self):
        """ Number of samples per pixel asked for by the preferences
            (gl_parameters "antialias" on/off and "antialias_samples").
            Qt's QOpenGLWidget owns its framebuffer, so do not redirect it
            through the GTK-oriented multisample framebuffer path. """
        if getattr(self.vm_session, "toolkit", None) == "Qt5":
            return 0
        gl_parameters = self.vm_config.gl_parameters
        if not gl_parameters.get("antialias", True) or self._msaa_failed:
            return 0
        try:
            samples = int(gl_parameters.get("antialias_samples", 4))
        except (TypeError, ValueError):
            samples = 4
        return samples if samples > 1 else 0
    
    def _msaa_begin(self):
        """ Redirects the drawing that follows to a multisampled
            framebuffer of the size of the current viewport; _msaa_end()
            resolves it into the framebuffer that was bound here (the
            Gtk.GLArea's, the macOS offscreen one, or render_to_image()'s).
            Does nothing when antialiasing is off or not available. """
        self._msaa_target = None
        samples = self.antialias_samples()
        if samples == 0:
            return
        try:
            target = int(GL.glGetIntegerv(GL.GL_FRAMEBUFFER_BINDING))
            if self._msaa_fbo is not None and target == int(self._msaa_fbo):
                # a previous render() raised before _msaa_end()
                return
            viewport = GL.glGetIntegerv(GL.GL_VIEWPORT)
            width, height = int(viewport[2]), int(viewport[3])
            if width <= 0 or height <= 0:
                return
            samples = min(samples, int(GL.glGetIntegerv(GL.GL_MAX_SAMPLES)))
            if samples < 2:
                return
            if self._msaa_fbo is None:
                self._msaa_fbo = GL.glGenFramebuffers(1)
                self._msaa_color_rbo = GL.glGenRenderbuffers(1)
                self._msaa_depth_rbo = GL.glGenRenderbuffers(1)
            GL.glBindFramebuffer(GL.GL_FRAMEBUFFER, self._msaa_fbo)
            if self._msaa_state != (width, height, samples):
                self._msaa_state = None
                GL.glBindRenderbuffer(GL.GL_RENDERBUFFER, self._msaa_color_rbo)
                GL.glRenderbufferStorageMultisample(GL.GL_RENDERBUFFER, samples, GL.GL_RGBA8,
                                                    width, height)
                GL.glBindRenderbuffer(GL.GL_RENDERBUFFER, self._msaa_depth_rbo)
                GL.glRenderbufferStorageMultisample(GL.GL_RENDERBUFFER, samples,
                                                    GL.GL_DEPTH_COMPONENT24, width, height)
                GL.glBindRenderbuffer(GL.GL_RENDERBUFFER, 0)
                GL.glFramebufferRenderbuffer(GL.GL_FRAMEBUFFER, GL.GL_COLOR_ATTACHMENT0,
                                             GL.GL_RENDERBUFFER, self._msaa_color_rbo)
                GL.glFramebufferRenderbuffer(GL.GL_FRAMEBUFFER, GL.GL_DEPTH_ATTACHMENT,
                                             GL.GL_RENDERBUFFER, self._msaa_depth_rbo)
                status = GL.glCheckFramebufferStatus(GL.GL_FRAMEBUFFER)
                if status != GL.GL_FRAMEBUFFER_COMPLETE:
                    raise RuntimeError("multisampled framebuffer incomplete (status 0x{:x})".format(int(status)))
                self._msaa_state = (width, height, samples)
            GL.glEnable(GL.GL_MULTISAMPLE)
            self._msaa_target = (target, width, height)
        except Exception as error:
            # never lose the frame because of the antialiasing: fall back
            # to drawing straight into the widget's framebuffer, for good
            logger.warning("antialiasing disabled, multisampled framebuffer not available: {}".format(error))
            self._msaa_failed = True
            self._msaa_target = None
            try:
                GL.glBindFramebuffer(GL.GL_FRAMEBUFFER, target)
            except Exception:
                pass
    
    def _msaa_end(self):
        """ Resolves the multisampled framebuffer into the one render()
            was called with, and leaves that one bound. """
        if self._msaa_target is None:
            return
        target, width, height = self._msaa_target
        self._msaa_target = None
        GL.glBindFramebuffer(GL.GL_READ_FRAMEBUFFER, self._msaa_fbo)
        GL.glBindFramebuffer(GL.GL_DRAW_FRAMEBUFFER, target)
        GL.glBlitFramebuffer(0, 0, width, height, 0, 0, width, height,
                             GL.GL_COLOR_BUFFER_BIT, GL.GL_NEAREST)
        GL.glBindFramebuffer(GL.GL_FRAMEBUFFER, target)

    
    # Representations left out of the depth map used for exported outlines
    # (see render_to_image(return_depth=True)): thin/flat overlays that
    # should not get a contour of their own.
    _NO_OUTLINE_REPRESENTATIONS = {"LinesRepresentation", "NonBondedRepresentation",
                                   "DotsRepresentation", "OneColorDotsRepresentation",
                                   "PickingDotsRepresentation", "DashedLinesRepresentation",
                                   "CellLineRepresentation", "LabelRepresentation",
                                   "WiresRepresentation"}

    def _read_solid_depth(self, width, height):
        """ Draws only the solid representations (spheres, sticks, cartoon,
            surfaces...) into the currently bound framebuffer and returns
            their depth as a (height, width) float32 array of distances from
            the camera in Angstrom, top row first, np.inf on the background.
            Also returns the size of one pixel, in Angstrom per unit of
            distance, and the surface normals as a (height, width, 3) uint8
            array (128 + 127 * n; the shaders' own smooth normals, so the
            facets of the stick/sphere meshes do not show up as creases).
        """
        # the shared lighting block writes normals instead of colours in
        # this mode; fog off so they are not blended with the background
        shading_mode, fog_enabled = self.shading_mode, self.glcamera.fog_enabled
        self.shading_mode = shaders_lighting.NORMAL_MAP_MODE
        self.glcamera.fog_enabled = False
        try:
            GL.glClearColor(0.5, 0.5, 0.5, 0.0)
            GL.glClear(GL.GL_COLOR_BUFFER_BIT | GL.GL_DEPTH_BUFFER_BIT)
            self.update_camera_ubo()
            for vm_object in self.vm_session.vm_objects_dic.values():
                if vm_object.active and vm_object.frames.shape[0] > 0:
                    for representation in vm_object.representations.values():
                        if (representation is not None and representation.active
                                and type(representation).__name__ not in self._NO_OUTLINE_REPRESENTATIONS):
                            representation.draw_representation()
            GL.glFinish()
        finally:
            self.shading_mode = shading_mode
            self.glcamera.fog_enabled = fog_enabled
        GL.glPixelStorei(GL.GL_PACK_ALIGNMENT, 1)
        data = GL.glReadPixels(0, 0, width, height, GL.GL_RGB, GL.GL_UNSIGNED_BYTE)
        normals = np.frombuffer(data, dtype=np.uint8).reshape((height, width, 3))
        normals = np.ascontiguousarray(np.flip(normals, axis=0))
        data = GL.glReadPixels(0, 0, width, height, GL.GL_DEPTH_COMPONENT, GL.GL_FLOAT)
        window_z = np.frombuffer(data, dtype=np.float32).reshape((height, width))
        window_z = np.flip(window_z, axis=0)
        # invert the projection (see mop.my_glPerspectivef): z_ndc = -A - B/z_eye;
        # abs() because only the distance matters here
        proj = self.glcamera.projection_matrix
        with np.errstate(divide="ignore"):
            distance = proj[3, 2] / (2.0 * window_z - 1.0 + proj[2, 2])
        distance = np.where(window_z >= 1.0, np.inf, np.abs(distance)).astype(np.float32)
        pixel_size = 2.0 / (float(proj[1, 1]) * height)
        return distance, pixel_size, normals

    def render_to_image(self, scale_factor=1, return_depth=False):
        """ Renders the current scene into an offscreen framebuffer at
            scale_factor times the widget's current resolution, then
            reads it back as a numpy RGBA array (height, width, 4) with
            origin at the top-left (already flipped from OpenGL's
            bottom-left convention, ready to hand to PIL/Image.fromarray).

            Rationale for *rendering* at a higher resolution instead of
            just upscaling the captured screenshot afterwards: naive
            image upscaling only blurs/interpolates existing pixels, it
            doesn't add detail. Rendering at scale_factor x resolution
            actually rasterizes thinner/more precise lines, smoother
            circles for spheres/dots, and crisper text, because the GPU
            recomputes every fragment at the higher pixel density.

            scale_factor must be uniform (same factor for width and
            height) so the aspect ratio - and therefore the projection
            matrix - stays exactly the same as the on-screen view; only
            the pixel density changes. This also means dot/point sizes
            (which several representations compute in screen pixels via
            vm_glcore.height - see representations.py) come out scaled
            consistently with line/stick thickness (which is computed in
            world units via the projection matrix), instead of looking
            disproportionately small at the higher resolution.

            Keyword arguments:
            scale_factor -- integer or float multiplier applied to both
                            width and height (e.g. 2 for 2x, 3 for 3x).

            return_depth -- also return the depth of the solid
                            representations, for outlines drawn from the
                            real geometry (see _read_solid_depth() and
                            vismol_gtkwidget.apply_contour_filter()).

            Returns:
            image -- np.uint8 array of shape (height*scale_factor,
                     width*scale_factor, 4), or None if the offscreen
                     framebuffer could not be created/completed.
            With return_depth: (image, depth) where depth is a dict with
            "distance", "pixel_size", "normals", "scale_factor" and
            "background" (its rgb colour); None on failure.
        """
        if scale_factor == 1 and not return_depth:
            # No offscreen framebuffer needed - just read whatever is
            # already in the on-screen framebuffer (caller is expected
            # to have rendered a fresh frame into it already).
            width = int(self.width)
            height = int(self.height)
            GL.glFinish()  # mesma garantia de sincronizacao do caminho off-screen abaixo
            data = GL.glReadPixels(0, 0, width, height, GL.GL_RGBA, GL.GL_UNSIGNED_BYTE)
            image = np.frombuffer(data, dtype=np.uint8).reshape((height, width, 4))
            return np.ascontiguousarray(np.flip(image, axis=0))

        # --- Save every bit of state render_to_image is about to touch ---
        orig_width  = self.width
        orig_height = self.height
        orig_right  = self.right
        orig_left   = self.left
        orig_proj_mat = np.copy(self.glcamera.projection_matrix)
        orig_aspect   = self.glcamera.viewport_aspect_ratio
        prev_fbo      = GL.glGetIntegerv(GL.GL_FRAMEBUFFER_BINDING)
        prev_viewport = GL.glGetIntegerv(GL.GL_VIEWPORT)

        new_width  = int(orig_width  * scale_factor)
        new_height = int(orig_height * scale_factor)

        fbo = color_tex = depth_rbo = None
        try:
            # --- Build the offscreen framebuffer at the higher resolution ---
            color_tex = GL.glGenTextures(1)
            GL.glBindTexture(GL.GL_TEXTURE_2D, color_tex)
            GL.glTexImage2D(GL.GL_TEXTURE_2D, 0, GL.GL_RGBA8, new_width, new_height,
                             0, GL.GL_RGBA, GL.GL_UNSIGNED_BYTE, None)
            GL.glTexParameteri(GL.GL_TEXTURE_2D, GL.GL_TEXTURE_MIN_FILTER, GL.GL_LINEAR)
            GL.glTexParameteri(GL.GL_TEXTURE_2D, GL.GL_TEXTURE_MAG_FILTER, GL.GL_LINEAR)

            depth_rbo = GL.glGenRenderbuffers(1)
            GL.glBindRenderbuffer(GL.GL_RENDERBUFFER, depth_rbo)
            GL.glRenderbufferStorage(GL.GL_RENDERBUFFER, GL.GL_DEPTH_COMPONENT24,
                                      new_width, new_height)

            fbo = GL.glGenFramebuffers(1)
            GL.glBindFramebuffer(GL.GL_FRAMEBUFFER, fbo)
            GL.glFramebufferTexture2D(GL.GL_FRAMEBUFFER, GL.GL_COLOR_ATTACHMENT0,
                                       GL.GL_TEXTURE_2D, color_tex, 0)
            GL.glFramebufferRenderbuffer(GL.GL_FRAMEBUFFER, GL.GL_DEPTH_ATTACHMENT,
                                          GL.GL_RENDERBUFFER, depth_rbo)

            status = GL.glCheckFramebufferStatus(GL.GL_FRAMEBUFFER)
            if status != GL.GL_FRAMEBUFFER_COMPLETE:
                logger.error("render_to_image: offscreen framebuffer incomplete (status=%s)", status)
                return None

            # --- Point width/height/projection at the higher resolution ---
            # Same aspect ratio (width and height scale by the same
            # factor), so the projection matrix comes out geometrically
            # identical to the on-screen one - this only changes pixel
            # density, not framing/zoom/distortion.
            self.width  = np.float32(new_width)
            self.height = np.float32(new_height)
            self.right  = self.width / self.height
            self.left   = -self.right
            self.glcamera.viewport_aspect_ratio = self.width / self.height
            # same effective-near clamp as resize_window() (the stored
            # z_near can be negative after widening the clipping slab)
            self.glcamera.set_projection_matrix(
                mop.my_glPerspectivef(self.glcamera.field_of_view,
                                      self.glcamera.viewport_aspect_ratio,
                                      max(self.glcamera.z_near, self.glcamera.min_znear),
                                      self.glcamera.z_far))

            GL.glViewport(0, 0, new_width, new_height)
            self.render()

            # [BUG FIX] Sem isso, glReadPixels podia rodar antes do comando
            # de desenho estar de fato completo na GPU -- normalmente
            # mascarado por uma superficie PREENCHIDA (rasterizacao mais
            # pesada, geralmente ja termina a tempo), mas exposto com
            # wireframe (poucos pixels efetivos, o comando de desenho volta
            # rapido demais do lado da CPU, driver/GPU-dependente). Reportado
            # como "PNG exportado sai completamente preto (so o fundo)" ao
            # exportar em 2x/3x/4x com uma superficie so' em modo wireframe
            # -- funcionava normalmente em 1x (sem framebuffer off-screen,
            # sem essa race). glFinish() bloqueia ate a GPU terminar TUDO
            # que foi submetido, garantindo que o proprio framebuffer
            # off-screen esta com o conteudo final antes de ler os pixels.
            GL.glFinish()

            data = GL.glReadPixels(0, 0, new_width, new_height, GL.GL_RGBA, GL.GL_UNSIGNED_BYTE)
            image = np.frombuffer(data, dtype=np.uint8).reshape((new_height, new_width, 4))
            image = np.ascontiguousarray(np.flip(image, axis=0))
            if return_depth:
                distance, pixel_size, normals = self._read_solid_depth(new_width, new_height)
                return image, {"distance": distance, "pixel_size": pixel_size,
                               "normals": normals, "scale_factor": scale_factor,
                               "background": [float(c) for c in self.bckgrnd_color[:3]]}
            return image

        finally:
            # --- Restore everything, even if something above raised ---
            self.width  = orig_width
            self.height = orig_height
            self.right  = orig_right
            self.left   = orig_left
            self.glcamera.viewport_aspect_ratio = orig_aspect
            self.glcamera.set_projection_matrix(orig_proj_mat)

            GL.glBindFramebuffer(GL.GL_FRAMEBUFFER, prev_fbo)
            GL.glViewport(int(prev_viewport[0]), int(prev_viewport[1]),
                          int(prev_viewport[2]), int(prev_viewport[3]))

            if fbo is not None:
                GL.glDeleteFramebuffers(1, [fbo])
            if color_tex is not None:
                GL.glDeleteTextures(1, [color_tex])
            if depth_rbo is not None:
                GL.glDeleteRenderbuffers(1, [depth_rbo])

            # [BUG FIX] Antes so' chamava self.queue_draw() aqui -- mas
            # isso e' ASSINCRONO (so' AGENDA um redraw pro GTK processar no
            # proximo ciclo do loop principal, nao garante que aconteca
            # antes do usuario ver/capturar o conteudo de novo). Como
            # render_to_image() e' chamado de fora do callback "render" da
            # propria GLArea (o clique de Export/Refresh acontece numa
            # janela SEPARADA, PreviewWindow), o binding de FBO que o GTK/
            # GDK considera "correto" pro proximo repaint natural pode nao
            # bater com prev_fbo (capturado via GL_FRAMEBUFFER_BINDING fora
            # do ciclo de render da propria GLArea -- nao confiavel em
            # todo backend/compositor). Resultado pratico: a GLArea principal
            # ficava em branco ate ALGO MAIS disparar um redraw de verdade
            # (ex.: mudar o zoom) -- exatamente o sintoma relatado.
            #
            # Agora, alem de agendar o proximo redraw natural do GTK
            # (queue_draw, mantido como rede de seguranca), forcamos um
            # RE-RENDER SINCRONO aqui mesmo -- usando a largura/altura/
            # projecao/FBO JA restaurados acima -- entao o conteudo on-screen
            # ja sai correto imediatamente, sem depender do agendamento
            # assincrono do GTK nem de nenhuma acao adicional do usuario.
            self.render()
            self.queue_draw()
    
    def _create_sphere_selection (self):
        """ Function doc """
        self.sphere_selection = VismolObject(self.vm_session, len(self.vm_session.vm_objects_dic), name='test')
        self.sphere_selection.set_model_matrix(self.model_mat)
        self.sphere_selection.create_representation(rep_type="picking_spheres")
        coords    = np.empty([1, 4, 3], dtype=np.float32)
        self.sphere_selection.frames = coords
        
        for index in range(0,4):
            atom = Atom(vismol_object = self.sphere_selection, 
                    name          ='Br'  , 
                    index         = index,
                    residue       = None , 
                    chain         = None, 
                    atom_id       = index,
                    occupancy     = 0, 
                    bfactor       = 0 ,
                    charge        = 0 )
            
            atom.vdw_rad  = 2.3 
            atom.cov_rad  = 2.3 
            atom.ball_rad = 2.3 
            color = [0.0,0.0,0.2]
            atom.color = np.array(color, dtype=np.float32)
            # unique_id and color_id (picking color) are None by default on
            # a freshly constructed Atom - they're normally set while
            # walking the parsed file in load_molecule(). This object is
            # built by hand, so set them explicitly; _generate_color_vectors
            # below needs a real atom.color_id to fill color_indexes.
            atom.unique_id = index
            atom._generate_atom_unique_color_id()
            self.sphere_selection.atoms[index] = atom
        
        # Build the object-level "colors"/"color_indexes" arrays from the
        # per-atom .color/.color_id values set above. Every normally loaded
        # VismolObject gets these via load_molecule() -> _generate_color_vectors();
        # this synthetic picking object skipped that step, so
        # SpheresRepresentation._colors_rads()/_sel_colors_rads() (which index
        # into vm_object.colors) would raise AttributeError as soon as
        # picking_spheres tried to draw.
        self.sphere_selection._generate_color_vectors(-1)
        
        self.sphere_selection.representations["picking_spheres"].define_new_indexes_to_vbo(range(0,4))
        self.sphere_selection.representations["picking_spheres"].active = True
        self.sphere_selection.representations["picking_spheres"].was_rep_ind_modified = True
        self.vm_session.vm_geometric_object_dic['picking_spheres'] = self.sphere_selection
        
        '''
        for index in range(1,5):
            atoms.append({"name"     : 'Au',
                          "resn"     : 'SEL',
                          "chain"    : 'A',
                          "resi"     : '1',
                          "occupancy": '1',
                          "bfactor"  : 0,
                          "charge"   : 0,
                          "index"    : index})
        
        for _atom in atoms:
            if _atom["chain"] not in vm_object.chains.keys():
                vm_object.chains[_atom["chain"]] = Chain(vm_object, name=_atom["chain"])
            _chain = vm_object.chains[_atom["chain"]]
            
            if _atom["resi"] not in _chain.residues.keys():
                _r = Residue(vm_object, name=_atom["resn"], index=_atom["resi"], chain=_chain)
                vm_object.residues[_atom["resi"]] = _r
                _chain.residues[_atom["resi"]] = _r
            _residue = _chain.residues[_atom["resi"]]
            
            atom = Atom(vm_object, name=_atom["name"], index=_atom["index"],
                        residue=_residue, chain=_chain, atom_id=atom_id,
                        occupancy=_atom["occupancy"], bfactor=_atom["bfactor"],
                        charge=_atom["charge"])
            atom.unique_id = unique_id
            atom._generate_atom_unique_color_id()
            _residue.atoms[atom_id] = atom
            vm_object.atoms[atom_id] = atom
            atom_id += 1
            unique_id += 1
        #logger.debug("Time used to build the tree: {:>8.5f} secs".format(time.time() - initial))
        vm_object.frames = PDBFiles.get_coords_from_raw_frames(rawframes, atom_id, vismol_session.vm_config.n_proc)
        '''

    def create_gl_programs(self):
        """ Function doc
        """
        logger.info("OpenGL version: {}".format(GL.glGetString(GL.GL_VERSION)))
        logger.info("OpenGL major version: {}".format(GL.glGetDoublev(GL.GL_MAJOR_VERSION)))
        logger.info("OpenGL minor version: {}".format(GL.glGetDoublev(GL.GL_MINOR_VERSION)))
        # Programs are about to be (re)linked, so any previously cached
        # uniform locations are no longer valid.
        self._uniform_loc_cache.clear()
        # Relink reinicia todos os uniforms do programa para o default no
        # driver, entao o cache de VALORES tambem precisa ser zerado: caso
        # contrario acreditariamos que um valor antigo ainda esta residente
        # e pulariamos o reenvio, deixando o shader com lixo.
        self._uniform_value_cache.clear()
        self._compile_shader_picking_dots()
        self._compile_shader_freetype()
        for rep in self.representations_available:
            func = getattr(self, "_compile_shader_" + rep)
            try:
                func()
            except AttributeError as ae:
                logger.error("Representation of type '{}' not implemented".format(rep))
                logger.error(ae)
    
    def load_shaders(self, vertex, fragment, geometry=None):
        """ Here the shaders are loaded and compiled to an OpenGL program. By default
            the constructor shaders will be used, if you want to change the shaders
            use this function. The flag is used to create only one OpenGL program.
            
            Keyword arguments:
            vertex -- The vertex shader to be used
            fragment -- The fragment shader to be used
        """
        my_vertex_shader = self.create_shader(vertex, GL.GL_VERTEX_SHADER)
        my_fragment_shader = self.create_shader(fragment, GL.GL_FRAGMENT_SHADER)
        if geometry is not None:
            my_geometry_shader = self.create_shader(geometry, GL.GL_GEOMETRY_SHADER)
        program = GL.glCreateProgram()
        GL.glAttachShader(program, my_vertex_shader)
        GL.glAttachShader(program, my_fragment_shader)
        if geometry is not None:
            GL.glAttachShader(program, my_geometry_shader)
        GL.glLinkProgram(program)
        # Checa o status do LINK. Sem isto, um programa que nao linka (ex.: GS
        # excedendo GL_MAX_GEOMETRY_TOTAL_OUTPUT_COMPONENTS) e devolvido como
        # objeto invalido e so estoura mais tarde, de forma confusa, num
        # glUseProgram com GL_INVALID_OPERATION. Aqui o erro aparece na hora,
        # com o info log do driver.
        if GL.glGetProgramiv(program, GL.GL_LINK_STATUS) != GL.GL_TRUE:
            info = GL.glGetProgramInfoLog(program)
            try:
                info = info.decode("utf-8", "replace")
            except (AttributeError, UnicodeDecodeError):
                pass
            logger.critical("Shader program link FAILED:\n{}".format(info))
            raise RuntimeError("Shader program link failed: {}".format(info))
        # Liga o bloco de camera (se o shader o declara) ao binding point
        # compartilhado. No-op para shaders ainda nao convertidos.
        self._bind_camera_ubo_to_program(program)
        return program
    
    def create_shader(self, shader_prog, shader_type):
        """ Creates, links to a source, compiles and returns a shader.
            
            Keyword arguments:
            shader -- The shader text to use
            shader_type -- The OpenGL enum type of shader, it can be:
                           GL.GL_VERTEX_SHADER, GL.GL_GEOMETRY_SHADER or GL.GL_FRAGMENT_SHADER
            
            Returns:
            A shader object identifier or pops out an error
        """
        shader = GL.glCreateShader(shader_type)
        # lit fragment shaders pull in the shared lighting block
        # (shaders/lighting.py); a no-op for every other shader
        shader_prog = shaders_lighting.inject(shader_prog)
        GL.glShaderSource(shader, shader_prog)
        GL.glCompileShader(shader)
        if GL.glGetShaderiv(shader, GL.GL_COMPILE_STATUS) != GL.GL_TRUE:
            logger.critical("Error compiling the shader: {}".format(shader_type))
            raise RuntimeError(logger.critical(GL.glGetShaderInfoLog(shader)))
        return shader

    def _selection_box_pick_new(self):
        """
        Select atoms using a rectangular selection box based on color picking.

        This method renders a special background where each atom is encoded as a unique color.
        Then, it reads pixels from the framebuffer inside the selection rectangle and decodes
        them back into atom IDs.
        """
        
        # Ensure mouse coordinates exist
        if not hasattr(self, "mouse_x") or not hasattr(self, "mouse_y"):
            return False
        
        # Original selection box corner (mouse press)
        x1, y1 = self.selection_box_x, self.selection_box_y
        # Current mouse position (convert Y from GTK to OpenGL coordinates)
        x2, y2 = self.mouse_x, self.height - self.mouse_y


        # Compute normalized rectangle (bottom-left corner + width/height)
        pos_x = int(max(0, min(x1, x2)))
        pos_y = int(max(0, min(y1, y2)))
        width = int(abs(x2 - x1))
        height = int(abs(y2 - y1))

        # Ignore empty selections
        if width == 0 or height == 0:
            return False

        # Clear buffers before rendering picking scene
        GL.glClear(GL.GL_COLOR_BUFFER_BIT | GL.GL_DEPTH_BUFFER_BIT)

        # Render all active objects using their selection representation
        for vm_object in self.vm_session.vm_objects_dic.values():
            if not vm_object.active:
                continue
            for rep in vm_object.representations.values():
                if rep and rep.active:
                    rep.draw_background_sel_representation()

        # Ensure tight packing of pixel data (no row alignment padding)
        GL.glPixelStorei(GL.GL_PACK_ALIGNMENT, 1)
        
        # Read pixels from framebuffer within the selection rectangle
        data = GL.glReadPixels(pos_x, pos_y, width, height,
                               GL.GL_RGBA, GL.GL_UNSIGNED_BYTE)

        # Decode RGB values into unique IDs
        # Each pixel encodes an ID as: R + G*256 + B*256^2
        picked_set = {
            data[i] + data[i+1]*256 + data[i+2]*256*256
            for i in range(0, len(data), 4)
        }

        atom_dic = self.vm_session.atom_dic_id
        
        # Convert IDs into atom objects, ignoring background (white = 16777215)
        selected = {
            atom_dic[pid] for pid in picked_set
            if pid != 16777215 and pid in atom_dic
        }

        # Apply selection (disable=False means additive selection)
        self.vm_session._selection_function_set(selected, disable=False)
        
        # Disable selection box picking mode
        self.selection_box_picking = False

    def _selection_box_pick(self):
        """ Selects a set of atoms from pixels obtained by the rectangle selection.  
            This function (method) is called in the render method, when the 
            "self.selection_box_picking" attribute is active. 
        
        """
        
        # glReadPixels and glReadnPixels return pixel data from the frame buffer, 
        # starting with the pixel whose lower left corner is at location (x, y), 
        # into client memory starting at location data.
        #
        # In GTK, x=0 and y=0 set to upper left corner (unlike openGL input data, 
        # the following lines do the coordinate conversion) 
        
        #try: # bachega 06 / 18 /2025
        #    selection_box_x2 = self.mouse_x
        #    selection_box_y2 = self.height - self.mouse_y
        #    selection_box_width  = selection_box_x2 - self.selection_box_x
        #    selection_box_height = selection_box_y2 - self.selection_box_y
        #except: # bachega 06 / 18 /2025
        #    return False
        
        
        # Bachega 03/22/2026 
        mouse_x = getattr(self, "mouse_x", None)
        mouse_y = getattr(self, "mouse_y", None)

        if mouse_x is None or mouse_y is None:
            return False
        x2 = mouse_x
        y2 = self.height - mouse_y  # convert to OpenGL coordinates
        
        try:
            width = x2 - self.selection_box_x
            height = y2 - self.selection_box_y
        except:
            return None

        '''
        #Simplificar cálculo do retângulo (grande melhoria)
        #Looking for the lower left corner of the checkbox
        if selection_box_width > 0 and selection_box_height > 0:
            pos_x = self.selection_box_x
            pos_y = self.selection_box_y
            width = selection_box_width
            height = selection_box_height
        
        elif selection_box_width < 0 and selection_box_height > 0:
            pos_x = selection_box_x2
            pos_y = self.selection_box_y
            width = -selection_box_width
            height = selection_box_height
        
        elif selection_box_width < 0 and selection_box_height < 0:
            pos_x = selection_box_x2
            pos_y = selection_box_y2
            width =  -selection_box_width
            height = -selection_box_height
        else:
            pos_x = self.selection_box_x
            pos_y = selection_box_y2
            width = selection_box_width
            height = -selection_box_height
        
        # taking the module from the width and height values 
        if pos_x < 0:
            pos_x = 0.0
        if pos_y < 0:
            pos_y = 0.0
        #'''
        
        #       cleaner and shoter bachega 03/22/2026
        # -------------------------------------------------
        x1 = self.selection_box_x
        y1 = self.selection_box_y
        x2 = self.mouse_x
        y2 = self.height - self.mouse_y

        pos_x = int(min(x1, x2))
        pos_y = int(min(y1, y2))
        width = int(abs(x2 - x1))
        height = int(abs(y2 - y1))
        # -------------------------------------------------

        # -------------------------------------------------
        GL.glClearColor(1, 1, 1, 1)
        GL.glClear(GL.GL_COLOR_BUFFER_BIT | GL.GL_DEPTH_BUFFER_BIT)
        # Mesma correcao do _pick: garante UBO de camera atualizado antes de
        # desenhar a cena de selecao (shaders sel_* convertidos leem daqui).
        self.update_camera_ubo()
        for index, vm_object in self.vm_session.vm_objects_dic.items():
            if vm_object.active:
                #vismol_object has few different types of representations
                '''
                for rep_name in vm_object.representations:
                    # checking all the representations in vismol_object.representations dictionary
                    if vm_object.representations[rep_name] is not None:
                        #  vismol_object.representations[rep_name] may be active or not  True/False
                        if vm_object.representations[rep_name].active:
                            vm_object.representations[rep_name].draw_background_sel_representation()
                '''
                # bachega 03/22/2026
                for rep in vm_object.representations.values():
                    if rep and rep.active:
                        rep.draw_background_sel_representation()
        # -------------------------------------------------

        '''
        # this is an old version
        GL.glPixelStorei(GL.GL_PACK_ALIGNMENT, 1)
        #data = GL.glReadPixels(float(pos_x), float(pos_y), width, height, GL.GL_RGBA, GL.GL_UNSIGNED_BYTE)
        # Bachega 03/22/2026
        data = GL.glReadPixels(int(pos_x), int(pos_y), int(width), int(height),
                       GL.GL_RGBA, GL.GL_UNSIGNED_BYTE)
        data = list(data)
        picked_set = set()
        for i in range(0, len(data), 4):
            #converting RGB values to atoms address (unique id)
            pickedID = data[i] + data[i+1] * 256 + data[i+2] * 256 * 256;
            picked_set.add(pickedID)
        _selected = set()
        for pickedID in picked_set:
            if pickedID == 16777215:
                pass
            else:
                self.atom_picked = self.vm_session.atom_dic_id[pickedID]
                _selected.add(self.vm_session.atom_dic_id[pickedID])
                # The disable variable does not allow, if the selected 
                # atom is already in the selected list, to be removed.
                # The disable variable is "False" for when we use 
                # selection by area (selection box)
                # self.vm_session._selection_function(selected=self.atom_picked, disable=False)
        self.vm_session._selection_function_set(_selected, disable=False)
        self.selection_box_picking = False
        '''


        GL.glPixelStorei(GL.GL_PACK_ALIGNMENT, 1)

        # Read pixels (RGBA, unsigned byte)
        data = GL.glReadPixels(
            int(pos_x), int(pos_y),
            int(width), int(height),
            GL.GL_RGBA, GL.GL_UNSIGNED_BYTE
        )

        BACKGROUND_ID = 16777215
        atom_dic = self.vm_session.atom_dic_id

        # Decode RGB → ID (ignore alpha)
        picked_ids = {
            data[i] + (data[i+1] << 8) + (data[i+2] << 16)
            for i in range(0, len(data), 4)
        }

        # Map IDs → atoms (filter background + missing IDs)
        selected = {
            atom_dic[pid]
            for pid in picked_ids
            if pid != BACKGROUND_ID and pid in atom_dic
        }

        # Apply selection
        self.vm_session._selection_function_set(selected, disable=False)

        # Disable selection mode
        self.selection_box_picking = False


        
    def _pick(self):
        """
        Perform single-pixel picking using color-encoded IDs.
        """

        BACKGROUND_ID = 16777215

        # Clear buffers before rendering picking scene
        GL.glClearColor(1, 1, 1, 1)
        GL.glClear(GL.GL_COLOR_BUFFER_BIT | GL.GL_DEPTH_BUFFER_BIT)
        # O picking roda ANTES do update_camera_ubo() do render() principal,
        # entao o UBO ainda teria as matrizes do frame anterior. Os shaders
        # sel_* ja convertidos ao UBO leem daqui, logo atualizamos agora para
        # garantir que a cena de picking use a camera ATUAL (corrige erro de
        # selecao logo apos rotacao/zoom).
        self.update_camera_ubo()

        # Render selection representations
        for vm_object in self.vm_session.vm_objects_dic.values():
            if not vm_object.active:
                continue

            for rep in vm_object.representations.values():
                if rep and rep.active:
                    rep.draw_background_sel_representation()

        # Ensure tight packing
        GL.glPixelStorei(GL.GL_PACK_ALIGNMENT, 1)

        # Convert coordinates (GTK → OpenGL)
        x = int(self.picking_x)
        y = int(self.height - self.picking_y)

        # Read single pixel
        data = GL.glReadPixels(x, y, 1, 1, GL.GL_RGBA, GL.GL_UNSIGNED_BYTE)

        # Decode RGB → ID
        pickedID = data[0] + data[1] * 256 + data[2] * 256 * 256

        atom_dic = self.vm_session.atom_dic_id

        if pickedID == BACKGROUND_ID:
            self.atom_picked = None

            if self.button == 1:
                self.vm_session._selection_function_set(None)
                self.button = None
        else:
            atom = atom_dic.get(pickedID)

            if atom is not None:
                self.atom_picked = atom

                if self.button == 1:
                    dprint(atom)
                    self.vm_session._selection_function_set({atom})
                    self.button = None
            else:
                logger.debug(f"pickedID {pickedID} not found")
                self.button = None

        self.picking = False
        return True

    def _pick_old(self):
        """ Function doc """
        GL.glClearColor(1, 1, 1, 1)
        GL.glClear(GL.GL_COLOR_BUFFER_BIT | GL.GL_DEPTH_BUFFER_BIT)
        for index, vm_object in self.vm_session.vm_objects_dic.items():
            if vm_object.active:
                #vismol_object has few different types of representations
                for rep_name in vm_object.representations:
                    # checking all the representations in vismol_object.representations dictionary
                    if vm_object.representations[rep_name] is None:
                        pass
                    else:
                        #  vismol_object.representations[rep_name] may be active or not  True/False
                        if vm_object.representations[rep_name].active:
                            vm_object.representations[rep_name].draw_background_sel_representation()
        
        GL.glPixelStorei(GL.GL_PACK_ALIGNMENT, 1)
        pos = [self.picking_x, self.height - self.picking_y]
        data = GL.glReadPixels(float(pos[0]), float(pos[1]), 1, 1, GL.GL_RGBA, GL.GL_UNSIGNED_BYTE)
        
        #converting RGB values to atoms address (unique id)
        pickedID = data[0] + data[1] * 256 + data[2] * 256 * 256;
        if pickedID == 16777215:
            self.atom_picked = None
            if self.button == 1:
                self.vm_session._selection_function_set(None)
                self.button = None
        else:
            try:
                # Using antialias, in some rare cases, the pick function is not 
                # identifying the right color of the selected atom. This event is 
                # rare, but can impair viewing if it is not properly ignored
                self.atom_picked = self.vm_session.atom_dic_id[pickedID]
                if self.button == 1:
                    self.vm_session._selection_function_set({self.atom_picked})
                    self.button = None
            except KeyError as ke:
                logger.debug("pickedID {} not found".format(pickedID))
                logger.debug(ke)
                self.button = None
        self.picking = False
        return True
    
    def _get_uniform_location(self, program, name):
        """ Cached wrapper around glGetUniformLocation.

            glGetUniformLocation does a string lookup against the linked
            program every time it is called. Locations don't change once a
            program is linked, so we look each one up once and reuse it for
            the lifetime of the program (cache is cleared in
            create_gl_programs when shaders are recompiled).
        """
        key = (program, name)
        loc = self._uniform_loc_cache.get(key)
        if loc is None:
            loc = GL.glGetUniformLocation(program, name)
            self._uniform_loc_cache[key] = loc
        return loc

    def _uniform_changed(self, program, name, value):
        """ Decide se o uniform 'name' precisa ser reenviado ao 'program'.

            Retorna (location, True) se o valor mudou (ou nunca foi enviado)
            desde a ultima vez, registrando o novo valor no cache. Retorna
            (location, False) se o valor identico ja esta residente no
            programa -- nesse caso o chamador deve PULAR o glUniform*.

            Motivacao (Gargalo: custo por-objeto): fog, luz e antialias sao
            os mesmos para todos os objetos do frame e mudam raramente. Antes,
            cada draw_representation reenviava todos eles incondicionalmente,
            gerando dezenas de chamadas GL por objeto * N objetos por frame.
            Com este diff, a partir do 2o objeto (e nos frames seguintes) o
            valor ja esta no cache e nada e reenviado.

            'value' precisa ser hashable de forma estavel. Escalares e tuplas
            servem direto; arrays numpy sao convertidos pelo chamador via
            _hashable() antes de chegar aqui.
        """
        loc = self._get_uniform_location(program, name)
        if loc == -1:
            # uniform inexistente/otimizado para fora: nada a enviar.
            return loc, False
        key = (program, name)
        if self._uniform_value_cache.get(key) == value:
            return loc, False
        self._uniform_value_cache[key] = value
        return loc, True

    @staticmethod
    def _hashable(arr):
        """ Converte um valor de uniform (escalar np, lista, ndarray) numa
            chave hashable e comparavel por igualdade para o cache. Mantem a
            ordem dos componentes. """
        try:
            return tuple(np.asarray(arr, dtype=np.float32).ravel().tolist())
        except (TypeError, ValueError):
            return float(arr)

    def load_fog(self, program):
        """ Load the fog parameters in the specified program
            
            fog_start -- The coordinates where the fog will begin (always
                         positive)
            fog_end -- The coordinates where the fog will begin (always positive
                       and greater than fog_start)
            fog_color -- The color for the fog (same as background)

            NOTA (cache de uniforms): so reenvia cada parametro quando ele
            muda em relacao ao ultimo valor residente no programa. Em regime
            permanente (camera/fundo parados) isto vira tres lookups de dict
            e zero chamadas GL.
        """
        val = self._hashable(self.glcamera.fog_start)
        loc, changed = self._uniform_changed(program, "fog_start", val)
        if changed:
            GL.glUniform1fv(loc, 1, self.glcamera.fog_start)
        val = self._hashable(self.glcamera.fog_end)
        loc, changed = self._uniform_changed(program, "fog_end", val)
        if changed:
            GL.glUniform1fv(loc, 1, self.glcamera.fog_end)
        val = self._hashable(self.bckgrnd_color)
        loc, changed = self._uniform_changed(program, "fog_color", val)
        if changed:
            GL.glUniform4fv(loc, 1, self.bckgrnd_color)
        # 2026-09-28 (Preferences > General > Fog): 0 = no fog, 1 = the
        # original full-strength fog. Every shader that has fog declares
        # it; programs that don't simply get loc == -1 and skip it.
        intensity = self.glcamera.get_effective_fog_intensity()
        loc, changed = self._uniform_changed(program, "fog_intensity", intensity)
        if changed:
            GL.glUniform1f(loc, intensity)

    def load_fog_legacy(self, program):
        """ Load the fog parameters in the specified program
            
            fog_start -- The coordinates where the fog will begin (always
                         positive)
            fog_end -- The coordinates where the fog will begin (always positive
                       and greater than fog_start)
            fog_color -- The color for the fog (same as background)
        """
        fog_s = self._get_uniform_location(program, "fog_start")
        GL.glUniform1fv(fog_s, 1, self.glcamera.fog_start)
        fog_e = self._get_uniform_location(program, "fog_end")
        GL.glUniform1fv(fog_e, 1, self.glcamera.fog_end)
        fog_c = self._get_uniform_location(program, "fog_color")
        GL.glUniform4fv(fog_c, 1, self.bckgrnd_color)
    
    def _ensure_camera_ubo(self):
        """ Cria o UBO de camera uma unica vez (lazy). Idempotente: chamadas
            seguintes apos a criacao nao fazem nada. Seguro chamar dentro do
            contexto GL (ex.: inicio de render). """
        if self._camera_ubo is not None:
            return
        self._camera_ubo = GL.glGenBuffers(1)
        GL.glBindBuffer(GL.GL_UNIFORM_BUFFER, self._camera_ubo)
        # Aloca o buffer vazio (DYNAMIC: atualizado todo frame)
        GL.glBufferData(GL.GL_UNIFORM_BUFFER, self._camera_ubo_size, None,
                        GL.GL_DYNAMIC_DRAW)
        # Liga o buffer inteiro ao binding point compartilhado
        GL.glBindBufferBase(GL.GL_UNIFORM_BUFFER, self.CAMERA_UBO_BINDING,
                            self._camera_ubo)
        GL.glBindBuffer(GL.GL_UNIFORM_BUFFER, 0)

    def update_camera_ubo(self):
        """ Sobe view_mat e proj_mat para o UBO. Chamar UMA vez por frame,
            antes de desenhar as representacoes. Substitui os dois
            glUniformMatrix4fv por-objeto que existiam em load_matrices. """
        self._ensure_camera_ubo()
        view = np.ascontiguousarray(self.glcamera.view_matrix, dtype=np.float32)
        proj = np.ascontiguousarray(self.glcamera.projection_matrix,
                                    dtype=np.float32)
        GL.glBindBuffer(GL.GL_UNIFORM_BUFFER, self._camera_ubo)
        # view em offset 0, proj em offset 64 (std140: mat4 ocupa 64 bytes)
        GL.glBufferSubData(GL.GL_UNIFORM_BUFFER, 0,  64, view)
        GL.glBufferSubData(GL.GL_UNIFORM_BUFFER, 64, 64, proj)
        GL.glBindBuffer(GL.GL_UNIFORM_BUFFER, 0)

    def _bind_camera_ubo_to_program(self, program):
        """ Liga o bloco 'CameraMatrices' de um programa ao binding point
            compartilhado. Chamar UMA vez por programa, logo apos o link.
            Se o programa nao declara o bloco (ex.: shader ainda nao
            convertido, ou que nao usa camera), glGetUniformBlockIndex
            devolve GL_INVALID_INDEX e simplesmente ignoramos -- isso e o
            que permite a convivencia com shaders nao-convertidos. """
        block_index = GL.glGetUniformBlockIndex(program, "CameraMatrices")
        if block_index != GL.GL_INVALID_INDEX:
            GL.glUniformBlockBinding(program, block_index,
                                     self.CAMERA_UBO_BINDING)

    def load_matrices(self, program=None, model_mat=None):
        """ Load the matrices to OpenGL.
            
            model_mat -- transformation matrix for the objects rendered
            view_mat -- transformation matrix for the camera used
            proj_mat -- matrix for the space to be visualized in the scene

            NOTA (Gargalo 2): para shaders JA CONVERTIDOS ao UBO de camera,
            view_mat e proj_mat vem do UBO (atualizado 1x/frame por
            update_camera_ubo) e NAO sao enviados aqui. Para shaders LEGADOS
            (ainda com 'uniform mat4 view_mat/proj_mat'), mantemos o envio
            individual como antes -- e isso que permite converter um shader
            de cada vez sem quebrar os demais.

            A deteccao usa o cache de uniform location: se 'view_mat' existe
            como uniform individual no programa (loc != -1), e legado.
        """
        model = self._get_uniform_location(program, "model_mat")
        GL.glUniformMatrix4fv(model, 1, GL.GL_FALSE, model_mat)
        # Fallback para shaders nao convertidos: se ainda tem view_mat como
        # uniform individual, envia view+proj como antes.
        view = self._get_uniform_location(program, "view_mat")
        if view != -1:
            GL.glUniformMatrix4fv(view, 1, GL.GL_FALSE,
                                  self.glcamera.view_matrix)
            proj = self._get_uniform_location(program, "proj_mat")
            if proj != -1:
                GL.glUniformMatrix4fv(proj, 1, GL.GL_FALSE,
                                      self.glcamera.projection_matrix)

    def load_dot_params(self, program):
        """ Function doc
        """
        # Extern line
        linewidth = np.float32(80 / abs(self.dist_cam_zrp))
        if linewidth > 3.73:
            linewidth = 3.73
        # Intern line
        antialias = np.float32(80 / abs(self.dist_cam_zrp))
        if antialias > 3.73:
            antialias = 3.73
        # Dot size factor
        dot_factor = np.float32(500 / abs(self.dist_cam_zrp))
        if dot_factor > 150.0:
            dot_factor = 150.0
        uni_vext_linewidth = self._get_uniform_location(program, "vert_ext_linewidth")
        GL.glUniform1fv(uni_vext_linewidth, 1, linewidth)
        uni_vint_antialias = self._get_uniform_location(program, "vert_int_antialias")
        GL.glUniform1fv(uni_vint_antialias, 1, antialias)
        uni_dot_size = self._get_uniform_location(program, "vert_dot_factor")
        GL.glUniform1fv(uni_dot_size, 1, dot_factor)
    
    def apply_light_parameters(self):
        """ (Re)reads every light/shading setting from gl_parameters. Call
            it again after changing them (Preferences) and redraw: the
            uniform cache in load_lights() picks up whatever changed.

            The key light is gl_parameters "light_position"/"light_intensity";
            "light_preset" may move it and adds the fill/rim lights (see
            shaders/lighting.py), and "extra_lights" appends user-defined
            ones ({"direction": [x, y, z], "intensity": [r, g, b]}), up to
            MAX_LIGHTS in total. Directions are fixed relative to the camera.
        """
        params = self.vm_config.gl_parameters
        self.light_position = params["light_position"]
        self.light_ambient_coef = params["light_ambient_coef"]
        self.light_shininess = params["light_shininess"]
        self.light_intensity = params["light_intensity"]
        preset = params.get("light_preset", "headlight")
        if preset not in shaders_lighting.LIGHT_PRESETS:
            preset = "headlight"
        key_dir = shaders_lighting.PRESET_KEY_DIRECTION.get(preset, self.light_position)
        lights = [{"direction": key_dir, "intensity": self.light_intensity}]
        lights += shaders_lighting.LIGHT_PRESETS[preset]
        lights += list(params.get("extra_lights", []))
        lights = lights[:shaders_lighting.MAX_LIGHTS]
        dirs = np.zeros((shaders_lighting.MAX_LIGHTS, 3), dtype=np.float32)
        ints = np.zeros((shaders_lighting.MAX_LIGHTS, 3), dtype=np.float32)
        for i, light in enumerate(lights):
            direction = np.array(light["direction"], dtype=np.float32)
            norm = np.linalg.norm(direction)
            dirs[i] = direction / norm if norm > 0 else (0.0, 0.0, 1.0)
            ints[i] = light["intensity"]
        self.n_lights = len(lights)
        self.light_directions = dirs
        self.light_intensities = ints
        self.light_ambient = (np.float32(self.light_ambient_coef)
                              * np.array(self.light_intensity, dtype=np.float32))
        mode = params.get("shading_mode", "standard")
        self.shading_mode = shaders_lighting.SHADING_MODE_IDS.get(mode, 0)
        self.flat_softness = float(params.get("flat_softness", 0.25))
        self.toon_bands = max(1.0, float(params.get("toon_bands", 3)))

    def load_lights(self, program):
        """ Uploads the shared lighting uniforms (shaders/lighting.py).

            Cache de uniforms: a posicao/intensidade da luz e constante entre
            objetos e quase sempre entre frames, entao so reenviamos no diff.
        """
        loc, changed = self._uniform_changed(program, "vm_n_lights", self.n_lights)
        if changed:
            GL.glUniform1i(loc, self.n_lights)
        val = self._hashable(self.light_directions)
        loc, changed = self._uniform_changed(program, "vm_light_dir", val)
        if changed:
            GL.glUniform3fv(loc, shaders_lighting.MAX_LIGHTS, self.light_directions)
        val = self._hashable(self.light_intensities)
        loc, changed = self._uniform_changed(program, "vm_light_int", val)
        if changed:
            GL.glUniform3fv(loc, shaders_lighting.MAX_LIGHTS, self.light_intensities)
        val = self._hashable(self.light_ambient)
        loc, changed = self._uniform_changed(program, "vm_ambient", val)
        if changed:
            GL.glUniform3fv(loc, 1, self.light_ambient)
        val = self._hashable(self.light_shininess)
        loc, changed = self._uniform_changed(program, "vm_shininess", val)
        if changed:
            GL.glUniform1f(loc, self.light_shininess)
        loc, changed = self._uniform_changed(program, "vm_shading_mode", self.shading_mode)
        if changed:
            GL.glUniform1i(loc, self.shading_mode)
        loc, changed = self._uniform_changed(program, "vm_flat_softness", self.flat_softness)
        if changed:
            GL.glUniform1f(loc, self.flat_softness)
        loc, changed = self._uniform_changed(program, "vm_toon_bands", self.toon_bands)
        if changed:
            GL.glUniform1f(loc, self.toon_bands)
    
    def load_antialias_params(self, program):
        """ Function doc

            Cache de uniforms: antialias_length e fixo e alias_color so muda
            com o fundo; reenvio apenas no diff. """
        val = self._hashable(0.05)
        loc, changed = self._uniform_changed(program, "antialias_length", val)
        if changed:
            GL.glUniform1fv(loc, 1, 0.05)
        val = self._hashable(self.bckgrnd_color[:3])
        loc, changed = self._uniform_changed(program, "alias_color", val)
        if changed:
            GL.glUniform3fv(loc, 1, self.bckgrnd_color[:3])
    
    def _draw_text_labels(self, vm_font, entries, string_shift=(0.0, 0.0), char_spacing=1.0):
        """ [EN] Unified label/text drawing routine, shared by
            _draw_labels(), _draw_distance_labels(), _draw_picking_label()
            and representations.LabelRepresentation.draw_representation()
            (the actual, active "Atom Labels" feature) -- previously each
            of those duplicated almost the exact same "build
            xyz_pos/uv_coords -> upload VBOs -> set GL state -> draw"
            sequence.

            More importantly, this is also where the billboard fix
            lives: earlier, each caller baked the per-character
            horizontal advance directly into the WORLD-space X
            coordinate on the CPU (`point[0] + i * char_width`), which
            only looked correct when the camera happened to be looking
            straight down -Z -- any other camera orientation sheared
            the text, because the advance direction was the world X
            axis, not the camera's screen-right direction. Here we
            instead upload the SAME anchor point (in world/model space)
            for every character of a string, plus that character's slot
            index (0, 1, 2, ...) as a separate per-vertex attribute.
            The geometry shader (shaders/vm_freetype.py) transforms the
            anchor to view space ONCE and then computes the advance
            along the view-space X/Y axes -- which are always the
            camera's screen-right/screen-up, regardless of camera
            rotation -- and scales both the advance and the glyph size
            by the point's distance from the camera (see
            VismolFont.zoom_sensitivity for how much/little that
            scaling actually happens).

            Input parameters:
                vm_font -- the VismolFont instance to draw with (already
                           carrying color/size/zoom_sensitivity/etc).
                entries -- iterable of (text, (x, y, z)) or
                           (text, (x, y, z), x_shift) tuples, where
                           (x, y, z) is the *world/model-space* anchor
                           point for that string (e.g. an atom's
                           coordinates, or a distance label's midpoint),
                           and the optional x_shift is a PER-STRING
                           horizontal nudge in character-size units
                           (e.g. -len(text)/2.0 to center a label of
                           that length on its anchor point). This is
                           baked directly into that string's char_idx
                           values, so -- unlike string_shift below --
                           different entries in the SAME call can each
                           have a different x_shift (needed for Atom
                           Labels, where many differently-sized strings
                           are drawn together in one glDrawArrays call).
                string_shift -- small constant (x, y) nudge applied to
                                EVERY string in this call alike, in
                                character-size units (e.g. to offset a
                                picking label like "#1" away from the
                                atom it names). Use this instead of a
                                per-entry x_shift when the nudge doesn't
                                depend on the text itself.
                char_spacing -- glyph advance as a fraction of the font's
                                char width (1.0 = the default spacing;
                                < 1 packs the glyphs, e.g. multi-digit
                                atom indices).

            Returns the number of glyphs (GL_POINTS) drawn.
        """
        if vm_font.vao is None:
            vm_font.make_freetype_font()
            vm_font.make_freetype_texture(self.core_shader_programs["freetype"])
        
        xyz_pos = []
        uv_coords = []
        char_idx = []
        
        GL.glBindTexture(GL.GL_TEXTURE_2D, vm_font.texture_id)
        for entry in entries:
            if len(entry) == 3:
                text, (x, y, z), x_shift = entry
            else:
                text, (x, y, z) = entry
                x_shift = 0.0
            point = np.array([x, y, z, 1], dtype=np.float32)
            point = np.dot(point, self.model_mat)
            for i, c in enumerate(text):
                c_id = ord(c)
                cx = c_id % 16
                cy = c_id // 16 - 2
                # Every glyph of this string shares the SAME anchor --
                # the advance (i + x_shift) is applied later, on the
                # GPU, in screen-aligned space (see the geometry
                # shader).
                xyz_pos.append(point[0])
                xyz_pos.append(point[1])
                xyz_pos.append(point[2])
                uv_coords.append(cx * vm_font.text_u)
                uv_coords.append(cy * vm_font.text_v)
                uv_coords.append((cx + 1) * vm_font.text_u)
                uv_coords.append((cy + 1) * vm_font.text_v)
                char_idx.append(float(i) * char_spacing + x_shift)
        
        chars = len(char_idx)
        if chars == 0:
            return 0
        
        xyz_pos = np.array(xyz_pos, dtype=np.float32)
        uv_coords = np.array(uv_coords, dtype=np.float32)
        char_idx = np.array(char_idx, dtype=np.float32)
        
        GL.glBindBuffer(GL.GL_ARRAY_BUFFER, vm_font.coord_vbo)
        GL.glBufferData(GL.GL_ARRAY_BUFFER, xyz_pos.itemsize * len(xyz_pos),
                        xyz_pos, GL.GL_DYNAMIC_DRAW)
        GL.glBindBuffer(GL.GL_ARRAY_BUFFER, vm_font.text_vbo)
        GL.glBufferData(GL.GL_ARRAY_BUFFER, uv_coords.itemsize * len(uv_coords),
                        uv_coords, GL.GL_DYNAMIC_DRAW)
        GL.glBindBuffer(GL.GL_ARRAY_BUFFER, vm_font.char_idx_vbo)
        GL.glBufferData(GL.GL_ARRAY_BUFFER, char_idx.itemsize * len(char_idx),
                        char_idx, GL.GL_DYNAMIC_DRAW)
        GL.glBindBuffer(GL.GL_ARRAY_BUFFER, 0)
        GL.glDisable(GL.GL_DEPTH_TEST)
        GL.glEnable(GL.GL_BLEND)
        GL.glBlendFunc(GL.GL_SRC_ALPHA, GL.GL_ONE_MINUS_SRC_ALPHA)
        GL.glUseProgram(self.core_shader_programs["freetype"])
        
        vm_font.string_shift = np.array(string_shift, dtype=np.float32)
        vm_font.load_matrices(self.core_shader_programs["freetype"],
                              self.glcamera.view_matrix,
                              self.glcamera.projection_matrix)
        vm_font.load_font_params(self.core_shader_programs["freetype"])
        
        GL.glBindVertexArray(vm_font.vao)
        GL.glDrawArrays(GL.GL_POINTS, 0, chars)
        GL.glDisable(GL.GL_BLEND)
        GL.glBindVertexArray(0)
        GL.glUseProgram(0)
        return chars

    def _draw_labels(self):
        """ Draws one "residue/atom_name/index" label per atom of every
            vm_object in the session. (Not currently wired into
            render() -- kept working/consistent with the other two
            label drawers in case it gets turned back on, but nothing
            in this refactor changes whether it's called.)
        """
        entries = []
        for vm_object in self.vm_session.vm_objects_dic.values():
            for index, atom in vm_object.atoms.items():
                text = atom.residue.name + '/' + atom.name + '/' + str(atom.index)
                frame = self._get_vismol_object_frame(atom.vm_object)
                entries.append((text, atom.coords(frame)))
        self._draw_text_labels(self.vm_font, entries)

    def _draw_distance_labels(self, vm_object):
        """ Draws the numeric distance label ("1.23") at the midpoint of
            a picking distance/dash line.
        """
        if self.vm_font_dist.vao is None:
            self.vm_font_dist.color = np.array(self.vm_session.vm_config.gl_parameters["pk_dist_label_color"], dtype=np.float32)
        self.vm_font_dist.char_res = 15
        text = '{:.2f}'.format(vm_object.dist)
        midpoint = (vm_object.midpoint[0], vm_object.midpoint[1], vm_object.midpoint[2])
        # x_shift = -len(text)/2.0 centers the string on the midpoint
        # (half the string's width to the left) instead of anchoring
        # its first character there -- replaces the old "-0.2"
        # world-space fudge factor, now applied in screen-aligned space
        # by the shader so it no longer skews with camera rotation.
        self._draw_text_labels(self.vm_font_dist, [(text, midpoint, -len(text) / 2.0)])
    
    def _draw_arc_labels(self, vm_object):
        """ The value of an angle or dihedral arc of the picked atoms (degrees). """
        text = getattr(vm_object, "label_text", None)
        pos  = getattr(vm_object, "label_pos", None)
        if not text or pos is None:
            return
        if self.vm_font_dist.vao is None:
            self.vm_font_dist.color = np.array(self.vm_session.vm_config.gl_parameters["pk_dist_label_color"], dtype=np.float32)
        self.vm_font_dist.char_res = 15
        self._draw_text_labels(self.vm_font_dist, [(text, (pos[0], pos[1], pos[2]), -len(text) / 2.0)])

    def _draw_picking_label(self):
        """ This function draws the labels of the atoms selected by the
            function picking #1 #2 #3 #4
        """
        if self.vm_font.vao is None:
            self.vm_font.color = np.array(self.vm_session.vm_config.gl_parameters["pk_label_color"], dtype=np.float32)
        
        entries = []
        for number, atom in enumerate(self.vm_session.picking_selections.picking_selections_list, start=1):
            if atom:
                text = "#" + str(number)
                frame = self._get_vismol_object_frame(atom.vm_object)
                entries.append((text, atom.coords(frame)))
        # Small constant nudge (in char-size units) so the "#N" label
        # doesn't sit exactly on top of the picking sphere/dot -- same
        # role as the old "-0.075"/"-0.05" world-space offsets, now
        # applied in screen-aligned space.
        self._draw_text_labels(self.vm_font, entries, string_shift=(-0.2, -0.15))
        
        '''
        atomlist = self.vm_session.picking_selections.picking_selections_list
        
        text_d1 = ""
        text_d2 = ""
        text_d3 = ""
        angle_1 = ""
        # this should be a new function later
        if atomlist[0] and atomlist[1]:
            
            crd1 = atomlist[0].coords(frame)
            crd2 = atomlist[1].coords(frame)
            d1 = ((crd2[0]-crd1[0])**2+ 
                 (crd2[1]-crd1[1])**2+ 
                 (crd2[2]-crd1[2])**2)**0.5
            
            text_d1 = "dist 1-2: {:7.5f} ".format(d1)
            
            
            if atomlist[2]:
                #print('distance #1 - #2: ', d)
                crd3 = atomlist[2].coords(frame)
                
                d2 = ((crd2[0]-crd3[0])**2+ 
                     ( crd2[1]-crd3[1])**2+ 
                     ( crd2[2]-crd3[2])**2)**0.5
                
                text_d2 = "dist 2-3: {:7.5f} ".format(d2)
                
                
                v1 = [crd2[0]-crd1[0], 
                      crd2[1]-crd1[1], 
                      crd2[2]-crd1[2]]
                v2 = [crd2[0]-crd3[0], 
                      crd2[1]-crd3[1], 
                      crd2[2]-crd3[2]]
                dot_product = sum(v1[i] * v2[i] for i in range(len(v1)))  
                
                magnitude_a = math.sqrt(sum(x**2 for x in v1))
                magnitude_b = math.sqrt(sum(x**2 for x in v2))
                cosine_theta = dot_product / (magnitude_a * magnitude_b)
                angle_rad = math.acos(cosine_theta)
                angle_deg = math.degrees(angle_rad)
                
                angle_1 = "angle 1-2-3: {:7.5f} ".format(angle_deg)
                
        dprint(text_d1, text_d2, text_d3, angle_1)
                #print(angle_deg)
        '''

        '''
        self._draw_distance_labels( distance = d)
        '''

    def _compile_shader_picking_dots(self):
        """ Function doc """
        
        safe =  self.vm_config.gl_parameters["picking_dots_safe"]
        if safe:
            self.core_shader_programs["picking_dots"] = self.load_shaders(shaders_pick.vertex_shader_picking_dots_safe,
                                                                          shaders_pick.fragment_shader_picking_dots_safe)
        else:
            self.core_shader_programs["picking_dots"] = self.load_shaders(shaders_pick.vertex_shader_picking_dots,
                                                                      shaders_pick.fragment_shader_picking_dots)
        
    def _compile_shader_dots(self):
        """ Function doc """
        dot_type = self.vm_config.gl_parameters["dot_type"]
        self.shader_programs["dots"] = self.load_shaders(shaders_dots.shader_type[dot_type]["vertex_shader"],
                                                shaders_dots.shader_type[dot_type]["fragment_shader"])
        self.shader_programs["dots_sel"] = self.load_shaders(shaders_dots.shader_type[dot_type]["sel_vertex_shader"],
                                                    shaders_dots.shader_type[dot_type]["sel_fragment_shader"])
    
    def _compile_shader_posdot_type(self):
        """ Function doc """
        dot_type = 1 #self.vm_config.gl_parameters["posdot_type"]
        self.shader_programs["posdot_type"] = self.load_shaders(shaders_dots.shader_type[dot_type]["vertex_shader"],
                                                shaders_dots.shader_type[dot_type]["fragment_shader"])
        self.shader_programs["posdot_type_sel"] = self.load_shaders(shaders_dots.shader_type[dot_type]["sel_vertex_shader"],
                                                    shaders_dots.shader_type[dot_type]["sel_fragment_shader"])
    
    def _compile_shader_lines(self):
        """ Function doc """
        line_type = self.vm_config.gl_parameters["line_type"]
        self.shader_programs["lines"] = self.load_shaders(shaders_lines.shader_type[line_type]["vertex_shader"],
                                                  shaders_lines.shader_type[line_type]["fragment_shader"],
                                                  shaders_lines.shader_type[line_type]["geometry_shader"])
        self.shader_programs["lines_sel"] = self.load_shaders(shaders_lines.shader_type[line_type]["sel_vertex_shader"],
                                                      shaders_lines.shader_type[line_type]["sel_fragment_shader"],
                                                      shaders_lines.shader_type[line_type]["sel_geometry_shader"])
        # second pass of LinesRepresentation: double/triple bonds (GL_LINES_ADJACENCY)
        self.shader_programs["lines_multi"] = self.load_shaders(shaders_lines.vertex_shader_lines_multi,
                                                                shaders_lines.fragment_shader_lines,
                                                                shaders_lines.geometry_shader_lines_multi)
        '''
        #sticks_type = self.vm_config.gl_parameters["sticks_type"]
        #self.shader_programs["lines"] = self.load_shaders(shaders_sticks.shader_type[sticks_type]["vertex_shader"],
        #                                           shaders_sticks.shader_type[sticks_type]["fragment_shader"],
        #                                           shaders_sticks.shader_type[sticks_type]["geometry_shader"])
        #self.shader_programs["lines_sel"] = self.load_shaders(shaders_sticks.shader_type[sticks_type]["sel_vertex_shader"],
        #                                               shaders_sticks.shader_type[sticks_type]["sel_fragment_shader"],
        #                                               shaders_sticks.shader_type[sticks_type]["sel_geometry_shader"])
        '''
        
    def _compile_shader_nonbonded(self):
        """ Function doc """
        self.shader_programs["nonbonded"] = self.load_shaders(shaders_nonbonded.vertex_shader_non_bonded,
                                                      shaders_nonbonded.fragment_shader_non_bonded,
                                                      shaders_nonbonded.geometry_shader_non_bonded)
        self.shader_programs["nonbonded_sel"] = self.load_shaders(shaders_nonbonded.sel_vertex_shader_non_bonded,
                                                          shaders_nonbonded.sel_fragment_shader_non_bonded,
                                                          shaders_nonbonded.sel_geometry_shader_non_bonded)
    
    def _compile_shader_sticks(self):
        """ Function doc """
        sticks_type = self.vm_config.gl_parameters["sticks_type"]
        self.shader_programs["sticks"] = self.load_shaders(shaders_sticks.shader_type[sticks_type]["vertex_shader"],
                                                   shaders_sticks.shader_type[sticks_type]["fragment_shader"],
                                                   shaders_sticks.shader_type[sticks_type]["geometry_shader"])
        self.shader_programs["sticks_sel"] = self.load_shaders(shaders_sticks.shader_type[sticks_type]["sel_vertex_shader"],
                                                       shaders_sticks.shader_type[sticks_type]["sel_fragment_shader"],
                                                       shaders_sticks.shader_type[sticks_type]["sel_geometry_shader"])
    
    def _compile_shader_ribbons(self):
        """ Function doc """
        #sticks_type = self.vm_config.gl_parameters["ribbon_type"]
        sticks_type = self.vm_config.gl_parameters["sticks_type"]
        sticks_type = 3
        self.shader_programs["ribbons"] = self.load_shaders(shaders_sticks.shader_type[sticks_type]["vertex_shader"],
                                                   shaders_sticks.shader_type[sticks_type]["fragment_shader"],
                                                   shaders_sticks.shader_type[sticks_type]["geometry_shader"])
        
        self.shader_programs["ribbons_sel"] = self.load_shaders(shaders_sticks.shader_type[sticks_type]["sel_vertex_shader"],
                                                       shaders_sticks.shader_type[sticks_type]["sel_fragment_shader"],
                                                       shaders_sticks.shader_type[sticks_type]["sel_geometry_shader"])
    
    def _compile_shader_dynamic(self):
        """ Function doc """
        sticks_type = self.vm_config.gl_parameters["sticks_type"]
        self.shader_programs["dynamic"] = self.load_shaders(shaders_sticks.shader_type[sticks_type]["vertex_shader"],
                                                   shaders_sticks.shader_type[sticks_type]["fragment_shader"],
                                                   shaders_sticks.shader_type[sticks_type]["geometry_shader"])
        self.shader_programs["dynamic_sel"] = self.load_shaders(shaders_sticks.shader_type[sticks_type]["sel_vertex_shader"],
                                                       shaders_sticks.shader_type[sticks_type]["sel_fragment_shader"],
                                                       shaders_sticks.shader_type[sticks_type]["sel_geometry_shader"])
    
    def _compile_shader_spheres(self):
        """ Function doc """
        self.shader_programs["spheres"] = self.load_shaders(shaders_spheres.vertex_shader_spheres,
                                                    shaders_spheres.fragment_shader_spheres)
        self.shader_programs["spheres_sel"] = self.load_shaders(shaders_spheres.sel_vertex_shader_spheres,
                                                        shaders_spheres.sel_fragment_shader_spheres)
    
    def _compile_shader_picking_spheres(self):
        """ Function doc """
        #print('\npicking_spheres'*10)
        #self.shader_programs["picking_spheres"] = self.load_shaders(shaders_spheres.vertex_shader_spheres,
        #                                                            shaders_spheres.fragment_shader_spheres)
        #self.shader_programs["picking_spheres_sel"] = self.load_shaders(shaders_spheres.vertex_shader_picking_spheres,
        #                                                                shaders_spheres.fragment_shader_picking_spheres)
        #
        self.shader_programs["picking_spheres"] = self.load_shaders(shaders_spheres.vertex_shader_picking_spheres,
                                                                    shaders_spheres.fragment_shader_picking_spheres)
        self.shader_programs["picking_spheres_sel"] = self.load_shaders(shaders_spheres.vertex_shader_picking_spheres,
                                                                        shaders_spheres.fragment_shader_picking_spheres)
    
    def _compile_shader_vdw_spheres(self):
        """ Function doc """
        self.shader_programs["vdw_spheres"] = self.load_shaders(shaders_spheres.sel_vertex_shader_spheres,
                                                    shaders_spheres.sel_fragment_shader_spheres)
        self.shader_programs["vdw_spheres_sel"] = self.load_shaders(shaders_spheres.sel_vertex_shader_spheres,
                                                        shaders_spheres.sel_fragment_shader_spheres)
    
    def _compile_shader_dash(self):
        """ Function doc """
        self.shader_programs["dash"] = self.load_shaders(shaders_dashed_lines.vertex_shader_dashed_lines,
                                                         shaders_dashed_lines.fragment_shader_dashed_lines,
                                                         shaders_dashed_lines.geometry_shader_dashed_lines)
        self.shader_programs["dash_sel"] = self.load_shaders( shaders_dashed_lines.sel_vertex_shader_dashed_lines,
                                                              shaders_dashed_lines.sel_fragment_shader_dashed_lines,
                                                              shaders_dashed_lines.sel_geometry_shader_dashed_lines)
    
    def _compile_shader_freetype(self):
        """ Function doc """
        self.core_shader_programs["freetype"] = self.load_shaders(shaders_vm_freetype.vertex_shader_freetype,
                                                         shaders_vm_freetype.fragment_shader_freetype,
                                                         shaders_vm_freetype.geometry_shader_freetype)
    
    def _compile_shader_static_freetype(self):
        """ Function doc """
        self.core_shader_programs["static_freetype"] = self.load_shaders(shaders_vm_freetype.static_vertex_shader_freetype,
                                                                         shaders_vm_freetype.static_fragment_shader_freetype,
                                                                         shaders_vm_freetype.static_geometry_shader_freetype)

    def _compile_shader_impostor(self):
        """ Function doc """
        im_type = self.vm_config.gl_parameters["impostor_type"]
        self.shader_programs["impostor"] = self.load_shaders(shaders_impostor.shader_type[im_type]["vertex_shader"],
                                                     shaders_impostor.shader_type[im_type]["fragment_shader"],
                                                     shaders_impostor.shader_type[im_type]["geometry_shader"])
        self.shader_programs["impostor_sel"] = self.load_shaders(shaders_impostor.shader_type[im_type]["sel_vertex_shader"],
                                                        shaders_impostor.shader_type[im_type]["sel_fragment_shader"],
                                                        shaders_impostor.shader_type[im_type]["sel_geometry_shader"])
    
    def _compile_shader_surface(self):
        """ Function doc """
        # Switched from the "lines" (wireframe) shader pair to the real
        # triangle + Phong "surface" shader pair. See geometry_shader_surface
        # in shaders/surface.py: it now derives a flat per-face normal from
        # the triangle's own edges instead of relying on a "vert_normal"
        # vertex attribute the VAO never populated.
        self.shader_programs["surface"] = self.load_shaders(shaders_surface.vertex_shader_surface,
                                                    shaders_surface.fragment_shader_surface,
                                                    shaders_surface.geometry_shader_surface)
        self.shader_programs["surface_sel"] = self.load_shaders(shaders_spheres.vertex_shader_spheres,
                                                        shaders_spheres.fragment_shader_spheres)
    
    # [EN] BUG FIX: this method used to sit inside a '''...''' dead-code
    # block (immediately followed by another one literally commented
    # "NOT IMPLEMENTED YET") -- so it was never actually a real method
    # on this class at all, despite looking like normal, live code on a
    # quick read. initialize()'s shader-compile loop
    # (for rep in self.representations_available: getattr(self,
    # "_compile_shader_" + rep)()) would hit AttributeError for
    # "cartoon" every single time, which its own except AttributeError
    # handler swallows silently (just an error-level log line) -- so
    # "cartoon" never actually landed in self.shader_programs, and the
    # very first attempt to draw a Cartoon representation crashed with
    # KeyError('cartoon') in representations.py's _check_vao_and_vbos().
    # Restored to live code once the underlying secondary-structure-
    # classification bug (originally in the now-replaced cartoon_BCK.py --
    # see the independent CartoonRepresentation/ribbon_backbone.py that
    # replaced it) was fixed -- that was almost certainly the original
    # reason this got commented out and marked "not implemented yet" to
    # begin with. Untouched by that replacement: still just loads
    # shaders/cartoon.py's own v_shader_triangles/f_shader_triangles
    # (now an independently-written shader pair, same export names).
    def _compile_shader_cartoon(self):
        """ Function doc """
        self.shader_programs["cartoon"] = self.load_shaders(shaders_cartoon.v_shader_triangles,
                                                    shaders_cartoon.f_shader_triangles)
        # [EN] 2026-09-26 fix: this used to compile "cartoon_sel" from the
        # SAME lit shader pair as the visible mesh above -- wrong for
        # picking (see shaders/cartoon.py's own sel_vertex_shader_cartoon/
        # sel_fragment_shader_cartoon docstring for why). Now a real flat/
        # unlit shader, matching every other representation's own "_sel"
        # convention (e.g. _compile_shader_lines() just below).
        self.shader_programs["cartoon_sel"] = self.load_shaders(shaders_cartoon.sel_vertex_shader_cartoon,
                                                    shaders_cartoon.sel_fragment_shader_cartoon)
    
    '''
    #----------------------------NOT IMPLEMENTED YET---------------------------#
    # def _dynamic_bonds_shaders(self):
    #     """ Function doc """
    #     self.shader_programs["dynamic"] = self.load_shaders(shaders_sticks.vertex_shader_sticks,
    #                                                 shaders_sticks.fragment_shader_sticks,
    #                                                 shaders_sticks.geometry_shader_sticks)
    #     self.shader_programs["dynamic_sel"] = self.load_shaders(shaders_sticks.sel_vertex_shader_sticks,
    #                                                    shaders_sticks.sel_fragment_shader_sticks,
    #                                                    shaders_sticks.sel_geometry_shader_sticks)
    
    # def _wires_dot_shaders(self):
    #     """ Function doc """
    #     self.shader_programs["wires"] = self.load_shaders(shaders_wires.vertex_shader_wires,
    #                                               shaders_wires.fragment_shader_wires,
    #                                               shaders_wires.geometry_shader_wires)
    #     self.shader_programs["wires_sel"] = self.load_shaders(shaders_spheres.vertex_shader_spheres,
    #                                                     shaders_spheres.fragment_shader_spheres)
    #----------------------------NOT IMPLEMENTED YET---------------------------#
    #'''
    
    def _safe_frame_coords(self, vismol_object):
        """ Function doc
        This function checks if the number of the called frame will not exceed
        the limit of frames that each object has. Allowing two objects with
        different trajectory sizes to be manipulated at the same time within the
        glArea
        """
        # [EN] BUG FIX (reported by the user: visualizing a Normal Modes
        # object crashed the render loop with "IndexError: index -1 is
        # out of bounds for axis 0 with size 0"). A vismol object can
        # legitimately have zero frames for a while -- e.g. Normal Modes
        # import creates the object with vismol_object.frames =
        # np.empty([0, natoms, 3]) BEFORE looping over the selected mode
        # folders to append real frames -- and any render/pick pass that
        # touches it in that window used to hit frames[-1] on an empty
        # array. Returning the (already empty) frames array itself keeps
        # the return type/dtype the caller expects (frame_coords.nbytes
        # == 0, so the resulting glBufferData call is a harmless no-op)
        # instead of raising.
        if vismol_object.frames.shape[0] == 0:
            return vismol_object.frames, 0
        # [EN] Per-object frame pin (VismolObject.pin_frame()): when
        # set, this object's own frame overrides the session-global
        # scrubber, letting e.g. a docking receptor stay on frame 17
        # while a ligand object elsewhere in the scene is on frame 4.
        # An object that never calls pin_frame() falls through to the
        # original global-frame behavior unchanged.
        if vismol_object._frame_override is not None:
            base_frame = vismol_object._frame_override
            if base_frame < 0:
                base_frame = 0
        else:
            if self.vm_session.frame < 0:
                self.vm_session.frame = 0
            base_frame = self.vm_session.frame
        if base_frame >= vismol_object.frames.shape[0] - 1:
            frame_coords = vismol_object.frames[vismol_object.frames.shape[0] - 1]
            frame = vismol_object.frames.shape[0] - 1
        else:
            frame_coords = vismol_object.frames[base_frame]
            frame = base_frame
        return frame_coords, frame

    def _get_vismol_object_frame(self, vismol_object):
        """ Function doc """
        # [EN] Same zero-frames guard as _safe_frame_coords() above.
        if vismol_object.frames.shape[0] == 0:
            return 0
        # [EN] Same per-object frame-pin override as _safe_frame_coords().
        if vismol_object._frame_override is not None:
            base_frame = vismol_object._frame_override
            if base_frame < 0:
                base_frame = 0
        else:
            if self.vm_session.frame < 0:
                self.vm_session.frame = 0
            base_frame = self.vm_session.frame
        if base_frame >= vismol_object.frames.shape[0] - 1:
            frame = vismol_object.frames.shape[0] - 1
        else:
            frame = base_frame
        return frame
    
    def get_viewport_pos(self, x, y):
        """ Function doc """
        px = (2.0 * x - self.width) / self.width
        py = (2.0 * y - self.height) / self.height
        return np.array([px, -py], dtype=np.float32)
    
    def _mouse_pos_old(self, x, y):
        """
        Use the ortho projection and viewport information
        to map from mouse co-ordinates back into world
        co-ordinates
        """
        px = x / self.width
        py = y / self.height
        px = self.left + px * (self.right - self.left)
        py = self.top + py * (self.bottom - self.top)
        pz = self.glcamera.z_near
        return px, py, pz

    def _mouse_pos(self, x, y):
        """
        Convert screen (mouse) coordinates to world coordinates
        using orthographic projection.

        Assumes:
        - Screen origin is top-left (e.g., GTK)
        - OpenGL origin is bottom-left
        """

        # Avoid division by zero
        if self.width == 0 or self.height == 0:
            return None

        # Normalize to [0, 1]
        nx = float(x) / float(self.width)
        ny = 1.0 - (float(y) / float(self.height))  # invert Y axis

        # Map to world coordinates
        world_x = self.left + nx * (self.right - self.left)
        world_y = self.bottom + ny * (self.top - self.bottom)

        # Use near plane as default Z
        world_z = self.glcamera.z_near

        return world_x, world_y, world_z

    def center_on_atom(self, atom):
        """ Function doc
        """
        frame_index = self._get_vismol_object_frame(atom.vm_object)
        self.center_on_coordinates(atom.vm_object, atom.coords(frame_index))

    def fit_to_object(self, vismol_object, margin=1.15, atom_indices=None):
        """Adjust the camera distance and clipping planes to show an object."""
        frames = np.asarray(vismol_object.frames)
        if frames.ndim != 3 or frames.shape[0] == 0 or frames.shape[1] == 0:
            raise ValueError("Cannot fit the camera to an object without coordinates")

        frame_index = int(getattr(self.vm_session, "frame", 0))
        frame_index = min(max(frame_index, 0), frames.shape[0] - 1)
        coordinates = np.asarray(frames[frame_index], dtype=np.float32)
        if atom_indices is not None:
            indices = np.asarray(atom_indices, dtype=np.intp)
            if indices.size == 0:
                raise ValueError("Cannot fit the camera to an empty atom selection")
            coordinates = coordinates[indices]
        if not np.isfinite(coordinates).all():
            raise ValueError("Cannot fit the camera to non-finite coordinates")

        coordinates_h = np.concatenate(
            (coordinates, np.ones((len(coordinates), 1), dtype=np.float32)),
            axis=1,
        )
        world_coordinates = coordinates_h @ vismol_object.model_mat
        radius = float(np.max(np.linalg.norm(world_coordinates[:, :3], axis=1)))
        if not math.isfinite(radius):
            raise ValueError("Cannot fit the camera to non-finite coordinates")
        if radius == 0.0:
            radius = 1.0

        camera = self.glcamera
        aspect = float(camera.viewport_aspect_ratio)
        fov = math.radians(float(camera.field_of_view))
        half_angle = math.atan(math.tan(fov / 2.0) * min(aspect, 1.0))
        fit_radius = radius * margin
        distance = fit_radius / math.sin(half_angle)

        direction = np.asarray(camera.get_position(), dtype=np.float32)
        direction_length = float(np.linalg.norm(direction))
        if direction_length == 0.0:
            direction = np.array([0.0, 0.0, 1.0], dtype=np.float32)
        else:
            direction /= direction_length

        camera.set_view_matrix(camera._get_view_matrix(direction * distance))
        camera.z_near = max(float(camera.min_znear), distance - fit_radius)
        camera.z_far = distance + fit_radius
        camera.update_projection()
        camera.update_fog()
        self.dist_cam_zrp = distance
        self.parent_widget.queue_draw()

    def center_on_coordinates_new(self, vismol_object, target):
        """
        Smoothly translate all objects so that the target coordinate becomes centered.
        
        Takes the coordinates of an atom in absolute coordinates and first
        transforms them in 4D world coordinates, then takes the unit vector
        of that atom position to generate the loop animation. To generate
        the animation, first obtains the distance from the zero reference
        point (always 0,0,0) to the atom, then divides this distance in a
        defined number of cycles, this result will be the step for
        translation. For the translation, the world will move a number of
        steps defined, and every new point will be finded by multiplying the
        unit vector by the step. As a final step, to avoid biases, the world
        will be translated to the atom position in world coordinates.
        The effects will be applied on the model matrices of every VisMol
        object and the model matrix of the window.
        
        """

        if not np.allclose(self.zero_reference_point, target):
            self.zero_reference_point[:] = target

            target_pos_4d = np.array([*target, 1.0], dtype=np.float32)

            # Transform to world space
            target_in_world = vismol_object.model_mat.T.dot(target_pos_4d)[:3]

            norm = np.linalg.norm(target_in_world)
            if norm == 0:
                return

            unit_vec = target_in_world / norm

            steps = self.vm_config.gl_parameters.get("center_steps", 15)
            step_size = norm / steps

            def iter_all_objects():
                yield from self.vm_session.vm_objects_dic.values()
                for obj in self.vm_session.vm_geometric_object_dic.values():
                    if obj:
                        yield obj

            # Animate movement
            for _ in range(steps):
                delta = unit_vec * step_size

                for obj in iter_all_objects():
                    obj.model_mat = mop.my_glTranslatef(obj.model_mat, -delta)

                # Trigger redraw (GTK)
                if self.vm_session.toolkit == "Gtk_3.0":
                    win = self.parent_widget.get_window()
                    win.invalidate_rect(None, False)
                    win.process_updates(False)
                elif self.vm_session.toolkit == "Qt5":
                    self.parent_widget.process_updates()
                else:
                    raise RuntimeError(
                        "Centering is not implemented for toolkit {!r}".format(
                            self.vm_session.toolkit
                        )
                    )

                time.sleep(self.vm_config.gl_parameters["center_on_coord_sleep_time"])

            # Final correction to eliminate accumulated error
            for obj in iter_all_objects():
                final_pos = obj.model_mat.T.dot(target_pos_4d)[:3]
                obj.model_mat = mop.my_glTranslatef(obj.model_mat, -final_pos)

            self.parent_widget.queue_draw()

    def center_on_coordinates(self, vismol_object, target):
        """ Takes the coordinates of an atom in absolute coordinates and first
            transforms them in 4D world coordinates, then takes the unit vector
            of that atom position to generate the loop animation. To generate
            the animation, first obtains the distance from the zero reference
            point (always 0,0,0) to the atom, then divides this distance in a
            defined number of cycles, this result will be the step for
            translation. For the translation, the world will move a number of
            steps defined, and every new point will be finded by multiplying the
            unit vector by the step. As a final step, to avoid biases, the world
            will be translated to the atom position in world coordinates.
            The effects will be applied on the model matrices of every VisMol
            object and the model matrix of the window.
        """
        if (self.zero_reference_point[0] != target[0]) or \
           (self.zero_reference_point[1] != target[1]) or \
           (self.zero_reference_point[2] != target[2]):
            self.zero_reference_point[:] = target
            pos = np.array([target[0],target[1],target[2],1], dtype=np.float32)
            model_pos = vismol_object.model_mat.T.dot(pos)[:3]
            self.model_mat = mop.my_glTranslatef(self.model_mat, -model_pos)
            unit_vec = model_pos / np.linalg.norm(model_pos)
            step = np.linalg.norm(model_pos)/15.0
            for i in range(15):
                to_move = unit_vec * step
                
                for index, vm_object in self.vm_session.vm_objects_dic.items():
                    vm_object.model_mat = mop.my_glTranslatef(vm_object.model_mat, -to_move)
                
                for key in self.vm_session.vm_geometric_object_dic.keys():
                    vm_object = self.vm_session.vm_geometric_object_dic[key]
                    if vm_object:
                        vm_object.model_mat = mop.my_glTranslatef(vm_object.model_mat, -to_move)
                
                # WARNING: Method only works with GTK!!!
                if self.vm_session.toolkit == "Gtk_3.0":
                    self.parent_widget.get_window().invalidate_rect(None, False)
                    self.parent_widget.get_window().process_updates(False)
                elif self.vm_session.toolkit == "Qt5":
                    self.parent_widget.process_updates()
                else:
                    raise RuntimeError(
                        "Centering is not implemented for toolkit {!r}".format(
                            self.vm_session.toolkit
                        )
                    )
                time.sleep(self.vm_config.gl_parameters["center_on_coord_sleep_time"])
            
            for index, vm_object in self.vm_session.vm_objects_dic.items():
                model_pos = vm_object.model_mat.T.dot(pos)[:3]
                vm_object.model_mat = mop.my_glTranslatef(vm_object.model_mat, -model_pos)
            
            for key in self.vm_session.vm_geometric_object_dic.keys():
                vm_object = self.vm_session.vm_geometric_object_dic[key]
                if vm_object:
                    model_pos = vm_object.model_mat.T.dot(pos)[:3]
                    vm_object.model_mat = mop.my_glTranslatef(vm_object.model_mat, -model_pos)
            #self.dragging = True
            self.parent_widget.queue_draw()
    
    def queue_draw(self):
        """ Function doc """
        self.parent_widget.queue_draw()
