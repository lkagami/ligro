#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
#  vismol_config.py
#  
#  Copyright 2022 Fernando <fernando@winter>
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

import os
import json

class VismolConfig:
    """ Class doc """
    
    def __init__ (self, vismol_session):
        """ Class initialiser """
        self.vismol_session = vismol_session
        self.gl_parameters = {"background_color": [0.0, 0.0, 0.0, 1.0],#[1.0, 1.0, 1.0, 1.0],#"background_color": [0.0, 0.0, 0.0, 1.0],
                              
                              "color_type": 0,
                              "dot_size": 2,
                              "dots_size": 2,
                              "dot_type": 1,
                              "dot_sel_size": 1.5,
                              "dashed_dist_lines_color" : [0.4, 0.4, 0.4, 1.0],
                              "line_width": 3,
                              "line_width_selection": 80,
                              "line_type": 0,
                              "line_color": 0,
                              "ribbon_width": 1000,
                              "ribbon_width_selection": 100,
                              "ribbon_type": 2,
                              "ribbon_color": 0,
                              "sphere_type": 0,
                              "sphere_scale": 0.20,
                              #"sphere_scale": 0.25,
                              "sphere_quality": 2,
                              "impostor_type": 0,
                              #"sticks_radius": 2.5,
                              "sticks_radius": 0.010,
                              "sticks_color": 0,
                              "sticks_type": 0,
                              # Liga/desliga o desenho de ligacoes duplas/triplas
                              # nos sticks. True = percebe e desenha as ordens;
                              # False = tudo como ligacao simples (cilindro unico).
                              # Pode ser alterado em runtime (ver SticksRepresentation).
                              "multiple_bonds": True,
                              "lines_bond_spacing": 1.0,
                              # Liga/desliga a representacao pontilhada das ligacoes
                              # de coordenacao metalica (metal-ligante). True =
                              # ligacoes com metal viram linha pontilhada e saem
                              # dos sticks; False = tratadas como ligacao normal.
                              "metal_dashed_bonds": True,
                              "antialias": True,
                              "antialias_samples": 4,
                              "mouse_rotation_sensibility" : 1.5,

                              "scroll_step": 0.9,
                              # per-tick exponential rate applied to the
                              # current slab WIDTH (z_far - z_near) -- a
                              # genuinely exponential (multiplicatively
                              # reversible) step, not a fixed one -- see
                              # VismolGLCore.mouse_scroll(); 0.138 gives half
                              # the per-tick step of an initial 0.26 default
                              # (user's own follow-up request for a finer
                              # scroll).
                              "scroll_step_fraction": 0.138,
                              # Fog settings -- see GLCamera.apply_fog_parameters().
                              "fog_enabled": True,
                              "fog_intensity": 100.0,
                              "fog_range_scale": 1.0,
                              "fog_custom_range": False,
                              "fog_start_percent": 50.0,
                              "fog_end_percent": 100.0,
                              "fog_start_reference": "far_plane",
                              "fog_center_offset": 5.0,
                              "field_of_view": 10,
                              "light_position": [0, 0, 10.0],
                              #"light_position": [-2.5, 2.5, 3.0],
                              "light_color": [ 1.0, 1.0, 1.0, 1.0],
                              "light_ambient_coef": 0.4,
                              "light_shininess": 5.5,
                              "light_intensity": [0.6, 0.6, 0.6],
                              "light_specular_color": [1.0, 1.0, 1.0],
                              # Lighting / shading -- see libgl/shaders/lighting.py
                              # and VismolGLCore.apply_light_parameters().
                              #   light_preset -- "headlight" (the single light
                              #       above) or "three_point" (key + fill + rim)
                              #   extra_lights -- more lights, each
                              #       {"direction": [x,y,z], "intensity": [r,g,b]},
                              #       fixed relative to the camera (max 4 in total)
                              #   shading_mode -- "standard", "flat", "toon", "matte",
                              #       "glossy" or "rim"
                              "light_preset": "headlight",
                              "extra_lights": [],
                              "shading_mode": "standard",
                              "flat_softness": 0.25,
                              "toon_bands": 3,
                              "center_on_coord_sleep_time": 0.01,
                              "gridsize": 0.8,
                              "maxbond": 2.4,
                              "bond_tolerance": 1.4,
                              #"dynamic_bond_tolerance": 2.0,
                              "picking_dots_color": [0.0, 1.0, 1.0],
                              "picking_dots_safe"          : True,
                              "pk_label_color"             : [1.0, 1.0, 1.0, 1.0],
                              "pk_dist_label_color"        : [1.0, 1.0, 0.0, 1.0],
                              # Font family (bundled .ttf filename, see
                              # vismol/libgl/fonts/) and size used to draw
                              # every text label in the glArea: atom
                              # labels, picking labels (#1, #2, ...) and
                              # distance labels. Customizable via the
                              # Preferences window ("Labels Font (glArea)").
                              "label_font_file"            : "Amiko-SemiBold.ttf",
                              "label_font_size"            : 0.35,
                              # Separate font family/size for the "Labels"
                              # representation (atom index, MM charge,
                              # residue name/index, chain -- see
                              # 'label_content' below), customizable
                              # independently from the picking/distance
                              # labels above via the Preferences window
                              # ("Atom Labels (glArea)").
                              "atom_label_font_file"       : "Amiko-SemiBold.ttf",
                              "atom_label_font_size"       : 0.35,
                              # What text each atom's persistent "Labels"
                              # representation shows: one of 'name',
                              # 'symbol', 'index', 'mm_charge',
                              # 'residue_name', 'residue_index', 'chain'.
                              # Set/applied from the Preferences window
                              # ("Atom Labels (glArea)"). The same options
                              # are also available per-selection via the
                              # glArea right-click menu (Show > labels).
                              "label_content"              : "name",
                              "label_show_all"             : False,

                              # [EN] Cartoon (ribbon) representation --
                              # "Preferences > Cartoon" tab. Every value
                              # here matches vismol.utils.ribbon_geometry.
                              # DEFAULT_CARTOON_SETTINGS exactly (the
                              # single source of truth for the actual
                              # numeric defaults) so a config missing
                              # these keys entirely (any .config.json
                              # saved before this feature existed) and a
                              # brand-new default config behave
                              # identically -- see representations.py's
                              # _cartoon_settings_from_gl_parameters().
                              # The 2 "*_tolerance" values are a
                              # PERCENTAGE (100 = published P-SEA
                              # thresholds as-is), not the 0-1 scale
                              # factor ribbon_geometry itself uses
                              # internally -- conversion happens in
                              # _cartoon_settings_from_gl_parameters().
                              "cartoon_helix_width"        : 1.60,
                              "cartoon_helix_height"       : 0.18,
                              "cartoon_strand_width"       : 1.40,
                              "cartoon_strand_height"      : 0.18,
                              "cartoon_arrow_width"        : 2.20,
                              "cartoon_coil_radius"        : 0.32,
                              "cartoon_helix_tolerance"    : 100.0,
                              "cartoon_strand_tolerance"   : 150.0,
                              "cartoon_min_run_helix"      : 3,
                              "cartoon_min_run_strand"     : 2,
                              "cartoon_transparency"       : 0,
                              "cartoon_arrow_taper_residues" : 2,
                              "cartoon_strand_smoothing"      : 4,
                              "cartoon_spline_detail"        : 5,
                              "cartoon_dynamic_secondary_structure" : False,
                              "surface_type"                : "ses",
                              "surface_probe_radius"        : 1.4,
                              "surface_grid_spacing"        : 0.8,
                              "surface_taubin_enabled"      : False,
                              "surface_taubin_iterations"   : 15,
                              "surface_color_mode"          : "atom",
                              "cartoon_color_helix"        : [0.90, 0.20, 0.55],
                              "cartoon_color_strand"       : [0.95, 0.85, 0.15],
                              "cartoon_color_coil"         : [0.85, 0.85, 0.85],
                              }
        self.n_proc = 2
        # self.representations_available = {"dots", "lines", "nonbonded", "dotted_lines",
        #                                   "ribbon", "sticks", "spheres", "impostor",
        #                                   "surface", "cartoon", "freetype",
        #                                   "picking_dots"}
        # [EN] "cartoon" re-added (was only in the commented-out version
        # above, not the active set below -- meaning self.representations
        # (built from this set in VismolObject.__init__) never had a
        # "cartoon" key pre-declared at all). Without it,
        # vm_session.show_or_hide() -- used by the terminal's "show"/
        # "hide" commands -- would raise KeyError('cartoon') on
        # `vm_object.representations[rep_type]`, since that line indexes
        # the dict directly rather than using .get(). Safe to re-enable
        # now that calculate_secondary_structure()'s vector-math bug
        # (see cartoon_BCK.py) is fixed -- see the changelog there.
        self.representations_available = {"dots", "lines", "nonbonded", "impostor",'dash', "posdot_type",
                                          "sticks", "spheres", 'ribbons', 'cartoon', #'ribbon_sphere', 
                                          'dynamic','vdw_spheres', 'picking_spheres', 'static_freetype'}
    
    
    def save_easyhybrid_config(self):
        """ Function doc """
        config_path = os.path.join(os.environ["HOME"], ".VisMol", "VismolConfig.json")
        with open(config_path, "w") as config_file:
            json.dump(self.gl_parameters, config_file, indent=2)
    
    def load_easyhybrid_config(self, config_path):
        """ Function doc """
        if not os.path.isfile(config_path):
          config_path = os.path.join(os.environ["HOME"], ".VisMol", "VismolConfig.json")
        # Keep a copy of the built-in defaults (set in __init__) so that
        # keys added in newer versions (e.g. 'label_font_file',
        # 'label_font_size') are still present even when loading an older
        # config file saved to disk before those keys existed. Without
        # this, code that indexes gl_parameters[key] directly (instead of
        # .get(key, default)) could raise KeyError after an update.
        defaults = self.gl_parameters
        with open(config_path, "r") as config_file:
            loaded = json.load(config_file)
        merged = dict(defaults)
        merged.update(loaded)
        self.gl_parameters = merged
    
