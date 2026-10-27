#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
#  vismol_object.py
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
from logging import getLogger
from vismol.model.atom import Atom
from vismol.model.bond import Bond
from vismol.model.chain import Chain
from vismol.model.residue import Residue
from vismol.model.molecule import Molecule
#from vismol.model.molecular_properties import COLOR_PALETTE
from vismol.libgl.vismol_font import VismolFont
from vismol.libgl.vismol_font import resolve_font_path, DEFAULT_FONT_FILE, DEFAULT_FONT_SIZE
from vismol.libgl.representations import DotsRepresentation
from vismol.libgl.representations import LinesRepresentation
from vismol.libgl.representations import NonBondedRepresentation
from vismol.libgl.representations import PickingDotsRepresentation
from vismol.libgl.representations import ImpostorRepresentation
from vismol.libgl.representations import SticksRepresentation
from vismol.libgl.representations import SpheresRepresentation
from vismol.libgl.representations import DashedLinesRepresentation
from vismol.libgl.representations import LabelRepresentation
from vismol.libgl.representations import SurfaceRepresentation
from vismol.libgl.representations import CartoonRepresentation
# from vismol.libgl.representations import WiresRepresentation
# from vismol.libgl.representations import RibbonsRepresentation
import vismol.utils.c_distances as cdist

logger = getLogger(__name__)

# [EN] BUG FIX: as chaves deste dicionario eram TODAS MAIUSCULAS ('CL', 'NA',
# 'FE', ...), mas self.atoms[i].symbol usa a notacao quimica padrao ('Cl',
# 'Na', 'Fe', ...) -- ver Atom._get_symbol() / periodic_table.elements_by_symbol.
# Resultado: para TODO elemento de duas letras (Cl, Na, Fe, Ca, Mg, Br, Zn,
# Mn, Ni, Cu, Co, Si, Al, Ar, Ne, Be, Li, He) o .get(symbol, 4) nunca batia
# e sempre caia no default 4 -- ou seja, a valencia maxima desses elementos
# era ignorada silenciosamente. Corrigido usando as chaves na mesma notacao
# de self.atoms[i].symbol (assim nao precisa de .upper()/.lower() espalhado
# pelo codigo que consome este dicionario).
GABEDIT_MAX_VALENCE = {
    'H' : 1, 'He': 0,
    'Li': 1, 'Be': 2, 'B' : 3, 'C' : 4, 'N' : 3, 'O' : 2, 'F' : 1, 'Ne': 0,
    'Na': 1, 'Mg': 2, 'Al': 3, 'Si': 4, 'P' : 3, 'S' : 2, 'Cl': 1, 'Ar': 0,
    'K' : 1, 'Ca': 2, 'Br': 1, 'I' : 1,
    'Fe': 2, 'Zn': 2, 'Cu': 2, 'Mn': 2, 'Ni': 2, 'Co': 2,
}

class VismolObject:
    """ Visual Object contains the information necessary for openGL to draw 
        a model on the screen. Everything that is represented in the graphical 
        form is stored in the form of a VismolObject.
        
        
        Arguments:
        - vismol_session: Vismol Session - Necessary to build the "atomtree_structure".
                          The vismol_session contains the atom_id_counter (self.vm_session.atom_id_counter).
        - index: Unique index for the VismolObject to find it in self.vm_session.vismol_objects_dic.
        - name: Label that describes the object (default: "UNK").
        - active: Boolean flag to enable/disable the object (default: False).
        - trajectory: A list of coordinates representing the trajectory frames.
                      Each frame is represented as a list of [x, y, z] coordinates (optional).
        - color_palette: Integer number to access the color pick for carbon atoms (0 = green, 1 = purple, ...) (optional).
        - bonds_pair_of_indexes: Pair of atoms used to define bonds, like [1, 3, 1, 17, 3, 4, 4, 20] (optional).
   
        
        Arguments
        
        name       = string  - Label that describes the object  
        atoms      = list of atoms  - [index, at_name, cov_rad,  at_pos, at_res_i, at_res_n, at_ch]
        vismol_session  = Vismol Session - Necessary to build the "atomtree_structure"
                     vismol_session contains the atom_id_counter (self.vm_session.atom_id_counter)
        trajectory = A list of coordinates - eg [ [x1,y1,z1, x2,y2,z2...], [x1,y1,z1, x2,y2,z2...]...]
                     One frame is is required at last.
    """
    def __init__(self, vismol_session, index = 0, name="UNK", active=False, trajectory=None,
                 color_palette=None, bonds_pair_of_indexes=None):
        """ Class initialiser """
        # References to vismol_session and its configuration
        self.vm_session = vismol_session
        self.vm_config = vismol_session.vm_config
                                   
                                   # Unique index for the VismolObject to find it in self.vm_session.vismol_objects_dic
        self.index = index         # import to find vboject in self.vm_session.vismol_objects_dic
        
        self.name = name
        self.active = active       # Boolean flag to enable/disable the object (show and hide) # for "show and hide"   enable/disable
        
        self.frames = trajectory  # A list of coordinates representing the trajectory frames
                                  # Each frame is represented as a list of [x, y, z] coordinates
                                  # (set as None if not provided)
        
        if color_palette is None:
            self.color_palette = self.vm_session.periodic_table.get_color_palette() #None #COLOR_PALETTE[0] # Default color palette for carbon atoms (0 = green)
        else:
            self.color_palette = color_palette #this is an integer num to access the color pick for carbon atoms (0 = green, 1 = purple, ...) 
                                               # Integer number to access the specified color palette
        
                                   
        self.editing = False       # for translate and rotate  xyz coords
                                   # Boolean flag for translate and rotate XYZ coordinates 
        
        self.mass_center = np.zeros(3, dtype=np.float32)
        # NOTE: kept as a plain default-font VismolFont. It is NOT the
        # font actually used to draw the "labels" representation (atom
        # index/charge/residue name/etc.) -- that one lives on
        # LabelRepresentation.vm_font (see vismol/libgl/representations.py),
        # created with the family/size configured in gl_parameters
        # ('atom_label_font_file'/'atom_label_font_size') and adjustable
        # live from the Preferences window ("Atom Labels (glArea)").
        self.vm_font = VismolFont() # Font object used for visualization
        self.coords = None          # Set as None for now (not defined in the provided code)
        
        
        # Dictionaries to store atoms, residues, chains, and molecules with their respective indexes as keys
        self.atoms = {}
        #self.residues = {}
        self.chains = {}
        self.molecules = {}

        # [EN] Free-form dict for callers to attach structured info to
        # this object (e.g. docking-ensemble frame->receptor_id maps,
        # per-pose score/rank -- see src/util/vismol_pdbqt_builder.py).
        # Empty by default; nothing in the render/pick path reads it.
        self.metadata = {}

        # [EN] Per-object frame pin: when set (an int), every frame
        # lookup for THIS object (see VismolGLCore._safe_frame_coords/
        # _get_vismol_object_frame and Atom.coords()) uses this value
        # instead of the session-global vm_session.frame, letting one
        # object be shown at a fixed frame while every other object
        # keeps following the global trajectory scrubber. None (the
        # default) means "follow the global frame", so an object that
        # never calls pin_frame() behaves exactly as before this
        # attribute existed. Set/cleared via pin_frame()/unpin_frame().
        self._frame_override = None

        self.atom_unique_id_dic = {}   # Dictionary to store atom unique identifiers (not yet defined in the code)
        self.selected_atom_ids = set() # Set to store the IDs of selected atoms (not yet defined in the code)
        
        
        self.bonds = None       # dict {(i,j): Bond} -- ligacoes indexadas pelo par
                                # de indices de atomos NORMALIZADO (menor, maior).
                                # Itere com self.bonds.values(); busque uma ligacao
                                # com self.get_bond(i, j).
        self.index_bonds = None # Pair of atoms, something like: [1, 3, 1, 17, 3, 4, 4, 20]
                                # Pair of atoms used to define bonds (set as None if not provided)
        
        self.bond_order_list = None  # A bond order list [1,1,1,2,1,1, and so on...]        
        self.bond_order_per_atom = None # Ordem de ligacao por atomo (VBO vert_bond_order)
        '''
        self.index_bonds is a flattened list of atom pairs: [i0, j0, i1, j1, i2, j2, ...] — two elements per bond (atom i, atom j).
        self.bond_order_list contains one bond order value per bond: [order0, order1, order2, ...].
        '''
        
        
        
        # Ligacoes de coordenacao metalica (>=1 metal) e covalentes (sem metal).
        # Arrays achatados [i0,j0, i1,j1, ...]. metal_bonds e desenhado como
        # linha pontilhada; covalent_bonds alimenta sticks/lines. As metalicas
        # CONTINUAM em index_bonds/bonds (o grafo nao muda) -- isto e so um
        # roteamento de VISUALIZACAO.
        self.metal_bonds = None
        self.covalent_bonds = None
        self.non_bonded_atoms = None # Array of indexes
                                     # Array of indexes for non-bonded atoms (not yet defined in the code)
                                     
        self.cov_radii_array = None  # a list of covalent radius values for all  --> will be used to calculate de bonds
                                     # List of covalent radius values for all atoms (not yet defined in the code)
        
        self.electronegativity_array = None # a list Electronegativity from UFF force field
                                            # see ../utils/elements.py , colunm 8 - En_UFF
          
        
        self.topology = {} # {92: [93, 99], 93: [92, 94, 96, 100], 99: [92],...}
                           # important to define molecules
        
        self.rings = [] # [[1,2,3,4], [5,5,6,7], ...] each element is list of indexes and a ring
                           
        self.representations = {}# Dictionary to store different visualization representations for the object
        # (initialized with None values for each representation type)
        
        for rep_type in self.vm_config.representations_available:
            self.representations[rep_type] = None
        self.representations['ribbon_sphere'] = None
        self.representations['stick_spheres'] = None
        self.representations['metal_dash'] = None  # ligacoes de coordenacao metalica
        
        # Additional transformation matrices for the object's visual representation
        self.model_mat = np.identity(4, dtype=np.float32)
        self.trans_mat = np.identity(4, dtype=np.float32)
        
        self.core_representations = {"picking_dots":None, "picking_text":None} # Core representation objects
        self.selection_dots_vao = None
        self.selection_dot_buffers = None
        self.picking_dots_vao = None
        self.picking_dot_buffers = None
        
        self.dynamic_bonds  = [] # Pair of atoms, something like: [[0,1,1,2,3,4] , [0,1,1,2], ...]
        # Cache de ordem de ligacao POR FRAME, paralelo a self.dynamic_bonds
        # (self.dynamic_bond_orders[f] e' um array numpy com uma ordem por
        # par de self.dynamic_bonds[f], na mesma ordem). None = ainda nao
        # calculado para aquele frame. Ver get_dynamic_bond_order_for_frame
        # e perceive_bond_order_for_pairs.
        self.dynamic_bond_orders = []
        # [NOVO] Overrides manuais de ordem de ligacao para Dynamic Bonds,
        # POR FRAME -- dict {frame_idx: {(i,j): ordem, ...}, ...}. Usado
        # pelo comando de terminal 'bond ... frame=...' (ver atom_ops.
        # set_dynamic_bond_order()) para forcar uma ordem especifica num
        # frame especifico, sem que a percepcao automatica (rodada a cada
        # frame, ver get_dynamic_bond_order_for_frame abaixo) a sobrescreva.
        # Esparso de proposito: so' tem entrada para os pares/frames que o
        # usuario realmente editou -- todo o resto continua vindo da
        # percepcao automatica normalmente.
        self.dynamic_manual_bond_orders = {}
                                # Like self.index_bonds but for each frame
        
        self.c_alpha_bonds = [] # List of pair of atoms defining dynamic bonds for each frame
        self.c_alpha_atoms = [] # List of pair of atoms defining C-alpha bonds
    
        ''' Cell and Symmetry'''
        # Cell and Symmetry attributes (not yet defined in the code)
        self.cell_parameters  = None
        self.cell_coordinates = None
        self.cell_indexes     = None
        self.cell_colors      = None
        self.cell_bonds       = None
        self.representations['labels'] = None

        # ------------------------------------------------------------------
        # [EN] AUDIT (requested explicitly): every attribute found anywhere
        # in the codebase being SET on a VismolObject instance (self.X in
        # this file's own methods, or vismol_object.X/vm_object.X elsewhere)
        # but that was never declared here in __init__ -- meaning a FRESH
        # object (one that hasn't yet gone through whichever code path first
        # assigns each of these) would raise AttributeError on read before
        # that. Declared here as None (exactly as requested -- not a
        # type-appropriate default like False/{}/[] for the ones that would
        # normally suggest one, e.g. is_builder_only), so every VismolObject
        # has a well-defined value for all of these from the moment it's
        # constructed, regardless of which code path (if any) touches it
        # afterwards. Grouped by where each one actually gets assigned, for
        # anyone maintaining this list later.

        # -- Builder / undo (gui/windows/builder/{atom_ops,empty_object}.py) --
        self.manual_bonds       = None  # atom_ops.py: set of (atom_id_a, atom_id_b) explicit bonds
        self.manual_bond_orders = None  # atom_ops.py: {(atom_id_a, atom_id_b): order}
        self.undo_stack         = None  # atom_ops.py: list of snapshots (push_undo_snapshot/undo)
        self.e_id               = None  # empty_object.py (also pDynamo2EasyHybrid/session.py for normal objects): pDynamo psystem index, None for builder-only objects
        self.is_builder_only    = None  # empty_object.py: True for objects created by the Builder (no pDynamo system backing them)
        self.key6               = None  # empty_object.py: short random tag used as part of this object's unique display name

        # -- eSession / treeview / GUI bookkeeping (gui/main/*.py, pdynamo/pDynamo2EasyHybrid/session.py) --
        self.e_treeview_iter            = None  # main_treeview.py: Gtk.TreeIter for this object's row in the main treeview
        self.e_treeview_iter_parent_key = None  # main_treeview.py: key of the parent row, for nested/grouped objects
        self.liststore_iter             = None  # main_window.py: Gtk.TreeIter for this object's row in a liststore-based widget
        self.is_surface                 = None  # main_treeview.py: True for surface/grid objects (as opposed to atomistic ones)
        self.e_sequence                 = None  # util/sequence_plot.py: cached sequence-plot data for this object

        # -- Surface/MEP analysis (gui/windows/analysis/surface_analysis_window.py) --
        self.mep_cmap_name     = None  # colormap name used for the last Molecular Electrostatic Potential rendering
        self.mep_vmin          = None  # MEP colour-scale lower bound
        self.mep_vmax          = None  # MEP colour-scale upper bound
        self.surface_trajectory = None # surface/grid frames, analogous to `frames` for atomistic trajectories

        # -- pDynamo integration (pdynamo/pDynamo2EasyHybrid/*.py) --
        self.normal_modes_dict = None  # import_trajectory.py: parsed normal-mode data, keyed by mode index
        self.results           = None  # simulations_mixin.py: last simulation's result payload

        # -- Rendering: per-atom colour arrays (this file's own
        #    _generate_color_vectors(), called by e.g. atom_ops.add_atom()) --
        self.colors        = None  # np.float32[n_atoms, 3] -- current display colour per atom
        self.color_indexes = None  # np.float32[n_atoms, 3] -- picking-ID colour per atom (see Atom._generate_atom_unique_color_id())
        self.color_rainbow = None  # np.float32[n_atoms, 3] -- rainbow/spectrum colour-by-index variant
        self.cov_dot_sizes = None  # np.float32[n_atoms]    -- covalent-radius-derived point size, per atom
        self.vdw_dot_sizes = None  # np.float32[n_atoms]    -- van-der-Waals-radius-derived point size, per atom

        # -- Rendering: ribbons/dot-surface GPU buffers (libgl/shapes.py) --
        self.ribbons_vao             = None
        self.ribbons_buffers         = None
        self.dots_surface_vao        = None
        self.dots_surface_buffers    = None
        self.sel_dots_surface_vao     = None
        self.sel_dots_surface_buffers = None
        # ------------------------------------------------------------------

    
    
    def build_core_representations(self):
        """
        Function doc: Builds core representations for the VismolObject.

        Note: The function creates core representation objects for the VismolObject and
        updates the 'core_representations' dictionary with the created representations.
        """
        
        # Create a PickingDotsRepresentation object and add it to the 'core_representations' dictionary
        self.core_representations["picking_dots"] = PickingDotsRepresentation(self,
                                                    self.vm_session.vm_glcore, active=True,
                                                    indexes=list(self.atoms.keys()))
        
        # Create a DashedLinesRepresentation object and add it to the 'core_representations' dictionary
        self.core_representations["dash"] = DashedLinesRepresentation(self, self.vm_session.vm_glcore,
                                                active=True, indexes=self.index_bonds)
    
    def _ensure_metal_dash(self):
        """ Cria e ativa a representacao pontilhada das ligacoes metalicas se
            houver ligacoes metalicas e ela ainda nao existir. Idempotente:
            chamar varias vezes nao recria. Respeita a flag de config
            'metal_dashed_bonds' (default True). """
        return False
        
        try:
            if not self.vm_config.gl_parameters.get("metal_dashed_bonds", True):
                return
        except Exception:
            pass
        met = getattr(self, "metal_bonds", None)
        if met is None or len(met) == 0:
            return
        if self.representations.get("metal_dash") is not None:
            self.representations["metal_dash"].active = True
            return
        self.create_representation(rep_type="metal_dash")

    def _covalent_indexes(self):
        """ Indices de ligacao para sticks/lines: covalentes apenas (as
            metalicas vao para a representacao pontilhada). Cai em index_bonds
            se a deteccao nao tiver rodado OU se a flag metal_dashed_bonds
            estiver desligada (nesse caso metais voltam a ser sticks normais).
        """
        try:
            if not self.vm_config.gl_parameters.get("metal_dashed_bonds", True):
                return self.index_bonds
        except Exception:
            pass
        cov = getattr(self, "covalent_bonds", None)
        if cov is not None:
            return cov
        return self.index_bonds

    def create_representation(self, rep_type="lines", indexes=None, show_all=True):
        """   
        Function doc: Creates a visualization representation for the VismolObject.

        Parameters:
        - rep_type: The type of representation to create.
        - indexes: List of atom indexes or bond indexes based on the representation type (optional).
    
        Note: The function creates a new visualization representation object for the VismolObject
        and updates its internal 'representations' dictionary. It initializes the representation
        according to the provided 'rep_type' and other relevant parameters.
        """
        # Create a DotsRepresentation object and add it to the 'representations' dictionary
        if rep_type == "dots":
            dot_indexes = list(self.atoms.keys()) if indexes is None else indexes
            self.representations["dots"] = DotsRepresentation(self, self.vm_session.vm_glcore,
                                                active=True, indexes=dot_indexes)
        elif rep_type == "lines":
            # Create a LinesRepresentation object and add it to the 'representations' dictionary
            # Usa apenas ligacoes covalentes; as metalicas viram pontilhado.
            self.representations["lines"] = LinesRepresentation(self, self.vm_session.vm_glcore,
                                                active=True, indexes=self._covalent_indexes())
            self._ensure_metal_dash()
        elif rep_type == "nonbonded":
            # Create a NonBondedRepresentation object and add it to the 'representations' dictionary
            self.representations["nonbonded"] = NonBondedRepresentation(self, self.vm_session.vm_glcore,
                                                    active=True, indexes=self.non_bonded_atoms)
        elif rep_type == "impostor":
            # Create an ImpostorRepresentation object and add it to the 'representations' dictionary
            self.representations["impostor"] = ImpostorRepresentation(self, self.vm_session.vm_glcore,
                                                active=True, indexes=list(self.atoms.keys()))
        elif rep_type == "sticks":
            # Create a SticksRepresentation object and add it to the 'representations' dictionary
            # Usa apenas ligacoes covalentes; as metalicas viram pontilhado.
            self.representations["sticks"] = SticksRepresentation(self, self.vm_session.vm_glcore,
                                                                  active=True, indexes=self._covalent_indexes())
            # Deteccao automatica (3a): se ha ligacoes metalicas, cria/ativa a
            # representacao pontilhada dedicada junto com os sticks.
            self._ensure_metal_dash()
        elif rep_type == 'stick_spheres':
            # Create a SpheresRepresentation object and add it to the 'representations' dictionary
            self.representations['stick_spheres'] = SpheresRepresentation(self, self.vm_session.vm_glcore,
                                                    active=True, indexes=list(self.atoms.keys()), mode = 3)
                           
                                                                  
        elif rep_type == "spheres":
            # Create a SpheresRepresentation object and add it to the 'representations' dictionary
            self.representations["spheres"] = SpheresRepresentation(self, self.vm_session.vm_glcore,
                                                    active=True, indexes=list(self.atoms.keys()) )
        
        elif rep_type == "picking_spheres":
            # Create a SpheresRepresentation object with picking mode and add it to the 'representations' dictionary
            self.representations["picking_spheres"] = SpheresRepresentation(self, self.vm_session.vm_glcore,
                                                    active=True, indexes=list(self.atoms.keys()), mode =1 )
        
        elif rep_type == "vdw_spheres":
            # Create a SpheresRepresentation object with Van der Waals radius and add it to the 'representations' dictionary
            self.representations["vdw_spheres"] = SpheresRepresentation(self, self.vm_session.vm_glcore,
                                                    active=True, indexes=list(self.atoms.keys()), vdw =True)

        elif rep_type == "dash":
            # Create a DashedLinesRepresentation object and add it to the 'representations' dictionary
            self.representations["dash"] = DashedLinesRepresentation(self, self.vm_session.vm_glcore,
                                                active=True, indexes=self.index_bonds)
        elif rep_type == "metal_dash":
            # Representacao pontilhada DEDICADA as ligacoes de coordenacao
            # metalica. Usa o mesmo shader 'dash', mas so com os indices das
            # ligacoes que envolvem metal (self.metal_bonds). Cor distinta para
            # diferenciar das medicoes/dash genericas.
            met = getattr(self, "metal_bonds", None)
            if met is None:
                met = np.array([], dtype=np.uint32)
            rep = DashedLinesRepresentation(self, self.vm_session.vm_glcore,
                                            active=True, indexes=met)
            # cor das ligacoes metalicas (cinza-claro); ajuste a gosto.
            rep.color2 = [0.6, 0.6, 0.65]
            self.representations["metal_dash"] = rep
        elif rep_type == "ribbons":
            # Create a SticksRepresentation object with 'ribbons' name and add it to the 'representations' dictionary
            self.representations["ribbons"] = SticksRepresentation(self, self.vm_session.vm_glcore,
                                                                   active=True, indexes=self.index_bonds, name  = 'ribbons')
        elif rep_type == 'ribbon_sphere':
            # Create a SpheresRepresentation object and add it to the 'representations' dictionary
            self.representations['ribbon_sphere'] = SpheresRepresentation(self, self.vm_session.vm_glcore,
                                                    active=True, indexes=list(self.atoms.keys()), mode = 2  )
                                                                    
        elif rep_type == "dynamic":
            # Create a SticksRepresentation object with 'dynamic' flag and add it to the 'representations' dictionary
            #print(self.dynamic_bonds)
            # [NOVO] Cor unica opcional para as dynamic bonds, lida das
            # preferencias (EasyHybrid Preferences). Se 'dynamic_bonds_single_color'
            # estiver ligado, usa 'dynamic_bonds_color' (RGB/RGBA, branco por
            # padrao); senao mantem a cor por atomo (uniform_color = None).
            gp = self.vm_session.vm_config.gl_parameters
            dyn_uniform_color = None
            if gp.get("dynamic_bonds_single_color", False):
                dyn_uniform_color = gp.get("dynamic_bonds_color", [1.0, 1.0, 1.0, 1.0])
            self.representations["dynamic"] = SticksRepresentation(self, self.vm_session.vm_glcore,
                                                                  active=True, indexes=self.index_bonds, is_dynamic = True,
                                                                  uniform_color = dyn_uniform_color)
        
        elif rep_type == "labels":
            # Create a LabelRepresentation object and add it to the 'representations' dictionary
            self.representations["labels"] = LabelRepresentation(vismol_object  = self  ,  
                                                                  vismol_glcore = self.vm_session.vm_glcore , 
                                                                  indexes       = [0,1,2] , 
                                                                  labels        = None     , 
                                                                  color         = [1, 1, 0, 1])
        elif rep_type == "surface":
            # [EN] User request: a general (VDW/SAS/SES) whole-molecule
            # surface, see vismol.utils... wait, util/molecular_surface.py
            # (util/molecular_surface.build_surface_mesh()) -- this branch
            # already existed (unused by the QC surface tool, which always
            # uses its own "surface1"/"surface2" keys instead) but never
            # passed `surface_name`, leaving self.surf_name at its own
            # None default; explicit now so draw_representation()'s own
            # `surface_trajectory[frame][self.surf_name]` lookup has a
            # readable key ('surface') instead of an implicit None one.
            self.representations["surface"] = SurfaceRepresentation(vismol_object = self                      ,
                                                                    vismol_glcore = self.vm_session.vm_glcore ,
                                                                    name          = 'surface'                 ,
                                                                    active        = True                      ,
                                                                    indexes       = []                        ,
                                                                    is_dynamic    = False                     ,
                                                                    surface_name  = 'surface'                 )


            #self.representations["surface"] = SurfaceRepresentation(vismol_object = self  ,
            #                                                        vismol_glcore = self.vm_session.vm_glcore)
                                                                    #rep_type      = "mol", 
                                                                    #vismol_object = self, 
                                                                    #vm_glcore     = self.vm_session.vm_glcore, 
                                                                    #indexes       = [])
            
            #self.representations["surface"] = SurfaceRepresentation(vismol_object  = self  ,  
            #                                                      vismol_glcore = self.vm_session.vm_glcore , 
            #                                                      indexes       = [0,1,2] , 
            #                                                      labels        = None     , 
            #                                                      color         = [1, 1, 0, 1])
        
        
        elif rep_type == 'cartoon':
            # [EN] Independent Cartoon representation -- see
            # CartoonRepresentation's own docstring (representations.py)
            # and vismol.utils.ribbon_backbone/ribbon_geometry for the
            # actual spline/secondary-structure/mesh work. Constructor no
            # longer takes rep_type='mol'/indexes=[] -- neither was ever
            # used by the mesh-building path, only by the old atom-frames-
            # based Representation.__init__() this class still doesn't call.
            # [EN] 2026-09-28: show_all=False is used by vm_session.
            # show_or_hide() (partial cartoon of a selection), which sets
            # the per-atom atom.cartoon flags itself; the default (True)
            # keeps the old "whole object" behaviour for every other caller.
            self.representations["cartoon"] = CartoonRepresentation (
                    vismol_object = self, vismol_glcore = self.vm_session.vm_glcore,
                    show_all = show_all )
        
        # elif rep_type == "dotted_lines":
        #     self.representations["dotted_lines"] = LinesRepresentation(self, self.vm_session.vm_glcore,
        #                                                                active=True, indexes=indexes)
        else:
            # Handle error when an unsupported representation type is provided
            logger.error("Representation {} not implemented".format(rep_type))
            raise NotImplementedError("Representation {} not implemented".format(rep_type))
    
    def _generate_color_vectors(self, colors_id_start, do_colors_raindow=True,
                                do_vdw_dot_sizes=True, do_cov_dot_sizes=True):
        """ (1) This method assigns to each atom of the system a 
            unique identifier based on the RGB color standard. 
            This identifier will be used in the selection function. 
            There are no two atoms with the same color ID in  
            
            (2) This method builds the "colors" np array that will 
            be sent to the GPU and which contains the RGB values 
            for each atom of the system.


            Parameters:
            - colors_id_start: The starting value for color IDs, which determines the color scheme for atoms.
            
            - do_colors_rainbow: If True, generates a rainbow color scheme based on the atom's position in the atom list.
                                 If False, the atom colors will be based on the atom's individual color attribute.
            
            - do_vdw_dot_sizes: If True, generates Van der Waals dot sizes based on the atom's Van der Waals radius.
            
            - do_cov_dot_sizes: If True, generates covalent dot sizes based on the atom's covalent radius.

            Note: The function updates several attributes of the VismolObject, including 'colors', 'color_indexes',
                  'color_rainbow', 'vdw_dot_sizes', and 'cov_dot_sizes', based on the provided parameters.
        """
        
        atom_qtty = len(self.atoms)
        
        if atom_qtty ==0:
            return False
        
        half = int(atom_qtty/2)
        quarter = int(atom_qtty/4)
        color_step = 1.0/(atom_qtty/4.0)
        red = 0.0
        green = 0.0
        blue = 1.0
        
        self.colors = np.empty([len(self.atoms), 3], dtype=np.float32)
        self.color_indexes = np.empty([len(self.atoms), 3], dtype=np.float32)
        if do_colors_raindow:
            self.color_rainbow = np.empty([len(self.atoms), 3], dtype=np.float32)
        if do_vdw_dot_sizes:
            self.vdw_dot_sizes = np.empty(len(self.atoms), dtype=np.float32)
        if do_cov_dot_sizes:
            self.cov_dot_sizes = np.empty(len(self.atoms), dtype=np.float32)
        
        for i, atom in self.atoms.items():
            self.colors[i] = atom.color
            self.color_indexes[i] = atom.color_id
            if do_vdw_dot_sizes: 
                self.vdw_dot_sizes[i] = atom.vdw_rad * 3
            if do_cov_dot_sizes: 
                self.cov_dot_sizes[i] = atom.cov_rad
            if do_colors_raindow:
                if i <= 1*quarter:
                    self.color_rainbow[i,:] = red, green, blue
                    green += color_step
                
                if (i >= 1*quarter) and (i <= 2*quarter):
                    self.color_rainbow[i,:] = red, green, blue
                    blue -= color_step
                
                if (i >= 2*quarter) and (i <= 3*quarter):
                    self.color_rainbow[i,:] = red, green, blue
                    red += color_step
                
                if (i >= 3*quarter) and (i <= 4*quarter):
                    self.color_rainbow[i,:] = red, green, blue
                    green -= color_step
    
    def define_bonds_from_external (self, index_bonds = [], bond_orders=None, internal = True):
        """ Function doc """
        if internal is True:
            self.index_bonds = index_bonds

            self._bonds_from_pair_of_indexes_list(external_orders=bond_orders)
            self._get_non_bonded_from_bonded_list()
            
            self._generate_topology_from_index_bonds()
            #self.find_rings(self.topology)
            self.define_molecules()
            self.define_Calpha_backbone()
        else:
            return index_bonds
    
    def find_bonded_and_nonbonded_atoms(self, selection=None, frame=0, gridsize=1.2,
                                         maxbond=2.4, tolerance=1.2, internal = True, debug = False):
        """
        Function doc: Determines bonded and nonbonded atoms based on selection in a VismolObject.
        
        Parameters:
        - selection: A dictionary containing the selected atoms (optional).
        - frame: Frame number for the atom coordinates (default is 0).
        - gridsize: Grid size for bond calculation (usually should not be changed, default is 0.8).
        - maxbond: Size of the maximum bond to be monitored (bonds greater than maxbond may be disregarded, default is 2.4 Å).
        - tolerance: Safety factor that multiplies the term (ra_cov + rb_cov)**2 (default is 1.4).
        - internal: If True, calculates the bindings for the object itself (used in the object's genesis).
                    If False, returns a list of atom_ids representing the bonds between atoms.

        Returns:
        If internal is True, the function updates the VismolObject's internal attributes:
        - index_bonds: List of pairs of atom indexes representing the bonded atoms.
        - bonds: List of Bond objects representing the bonds between atoms.
        - topology: Dictionary representing the atom topology and connectivity.
        - Non-bonded interactions and other attributes are also updated.

        If internal is False, the function returns a list of atom_ids representing the bonds between atoms.

        Note: This function calculates bonds between atoms based on their positions and covalent radii.
        
            Receives a dictionary as a selection:

            selection = {'atom_id': atom_object, ...} ()
            
            frame - (frame number)
            
            grid_size - (usually should not be changed. Default is: 0.8)
            
            maxbond - (Size of the maximum bond to be monitored. Bonds greater than "maxbond" 
            may be disregarded. Default is: 2.4 A)

            tolerance (Safety factor that multiplies the term (ra_cov + rb_cov)**2. Default is: 1.4)

            internal (If True, calculates the bindings for the object itself. Used in the object's 
            
            genesis. If False, returns a list of atom_ids like  [0,1, # bond  between 0 and 1 
                                                                 1,2, # bond  between 1 and 2
                                                                 0,3] # bond  between 0 and 3 )
        

        """
        # Check if internal is True and there is already information about contacts.
        if internal:
            if self.index_bonds is not None:
                logger.critical("It seems that there is already information about "\
                    "the contacts in this VismolObject, trying to override the data "\
                    "can produce serious problems :(")
        
        
        # Initialize variables and data structures
        initial = time.time()
        atoms_frame_mask = np.zeros(len(self.atoms), bool)
        if self.cov_radii_array is None:
            
            self.cov_radii_array = np.empty(len(self.atoms), dtype=np.float32)
            self.electronegativity_array = np.empty(len(self.atoms), dtype=np.float32)
            
            for i, atom in self.atoms.items():
                self.cov_radii_array[i] = atom.cov_rad
                self.electronegativity_array[i] = atom.electronegativity
                
        # Create a mask to identify atoms in the frame (all atoms if selection is None)
        if selection is None:
            selection = self.atoms
            atoms_frame_mask[:] = True
        else:
            atoms_frame_mask[:] = False
            for atom in selection.values():
                atoms_frame_mask[atom.atom_id] = True
        
        
        # Extract relevant data for bond calculation
        #cov_rads = self.cov_radii_array[atoms_frame_mask]
        #coords = self.frames[frame][atoms_frame_mask]
        cov_rads = self.cov_radii_array
        coords = self.frames[frame]#[atoms_frame_mask]
        
        
        # Generate indexes and grid positions for each atom in the selection
        indexes = []
        gridpos_list = []
        for atom in selection.values():
            indexes.append(atom.atom_id)
            gridpos_list.append(atom.get_grid_position(gridsize=gridsize, frame=frame))
        #print(  'grid elementes:',gridpos_list)
        if debug:
            logger.debug("Time used for preparing the atom mask, covalent radii list "\
                         "and grid positions: {}".format(time.time() - initial))
        
        
        # Calculate bonds based on grid positions and covalent radii
        if internal is True:
            initial = time.time()
            self.index_bonds = cdist.get_atomic_bonds_from_grid(indexes, coords,
                                            cov_rads, gridpos_list, gridsize, maxbond, tolerance)
            msg = """Building grid elements  :
        Total number of Atoms   : {}
        Gridsize                : {}
        Bonds                   : {}
        Bonds calcultation time : {} seconds""".format(len(selection), gridsize,
                                    len(self.index_bonds), time.time() - initial)
            logger.info(msg)
            
            # Create Bond objects and update atom bond lists
            self._bonds_from_pair_of_indexes_list()
            
            # Generate non-bonded interactions from bonded atoms
            self._get_non_bonded_from_bonded_list()
            
            # Generate atom topology from index_bonds
            #bachega's code
            initial = time.time()
            self._generate_topology_from_index_bonds()
            
            # Define molecules based on the atom topology
            self.define_molecules()
            
            
            final = time.time()
            dprint('        Defining molecule indexes: ', final - initial)
            
            # Define Calpha backbone atoms
            self.define_Calpha_backbone()

        else:            
            # Calculate bonds based on grid positions and covalent radii and return the results
            index_bonds = cdist.get_atomic_bonds_from_grid(indexes, coords,
                                            cov_rads, gridpos_list, gridsize, maxbond, tolerance)
            
#            msg = """Building grid elements  :
#        Total number of Atoms   : {}
#        Gridsize                : {}
#        Bonds                   : {}
#        Bonds calcultation time : {} seconds""".format(len(selection), gridsize,
#                                    len(index_bonds), time.time() - initial)
#            logger.info(msg)
            
                       
            # Dynamic bonds is not working properly for pure QC system
            #'''
            bonds = index_bonds
            new_bonds = []
            for index in range (0,len(index_bonds), 2):
                if [bonds[index], bonds[index+1]] in new_bonds or [bonds[index+1], bonds[index]] in new_bonds:
                    pass
                else:
                    new_bonds.append([bonds[index], bonds[index+1]])
            #print(new_bonds)
            index_bonds = [item for sublista in new_bonds for item in sublista]
            #'''
            return index_bonds
    
    def _get_covalent_radii (self):
        """ Function doc """
        if self.cov_radii_array is None:
            self.cov_radii_array = np.empty(len(self.atoms), dtype=np.float32)
            for i, atom in self.atoms.items():
                self.cov_radii_array[i] = atom.cov_rad

    def _bonds_from_pair_of_indexes_list_old(self, exclude_list = [['H','H']], external_orders=None ):
        """ 
        Creates Bond objects based on pairs of indexes in self.index_bonds list.
        The bonds list is populated with the created Bond objects, and each
        atom involved in a bond is updated with the respective Bond object.
        
        self.index_bonds = [0,1  , 0,4  ,  1,3  , ...]
        self.bonds = {(i,j): bond, ...}  # dict, chave = par de indices 
        """
        assert self.bonds is None # Ensure the bonds list is not already initialized
        self.bonds = {} # dict {(i,j): Bond}, chave = par de indices normalizado
        self.bond_order_list =[] # Initialize an empty list to store the Bond orders
        new_index_bonds = []
        # Loop through the self.index_bonds list in pairs
        for i in range(0, len(self.index_bonds)-1, 2):
            
            index_i = self.index_bonds[i]    # Get the first atom's index of the bond
            index_j = self.index_bonds[i+1]  # Get the second atom's index of the bond
            dprint(index_i, index_j)
            is_excluded = False
            
            for excluded_bond in exclude_list: 
                if self.atoms[index_i].symbol in excluded_bond and self.atoms[index_j].symbol in excluded_bond:
                    is_excluded = True
            
            
            if is_excluded:
                pass
            else:
                new_index_bonds.append(index_i)
                new_index_bonds.append(index_j)
            
                #index_i = self.index_bonds[i]    # Get the first atom's index of the bond
                #index_j = self.index_bonds[i+1]  # Get the second atom's index of the bond
                
                
                # Create a Bond object with the atoms and their indexes
                bond = Bond(atom_i=self.atoms[index_i], atom_index_i=index_i,
                            atom_j=self.atoms[index_j], atom_index_j=index_j)
                

                # Define a ordem de ligação:
                #  - se veio ordem externa, usa a posição bond_pair_idx
                #    (Convenção B: lista já filtrada, então o k-ésimo
                #     sobrevivente casa com external_orders[k]);
                #  - senão, cai no palpite geométrico (UFF).
                if external_orders is not None:
                    #bond.bond_order = int(external_orders[bond_pair_idx])
                    dprint('here')
                    bond.bond_order = int(external_orders[i])+1
                else:
                    bond.get_bond_order()
                
                
                #bond.get_bond_order()       # fallback UFF
                
                
                #bond.get_bond_order()
                self.bond_order_list.append(bond.bond_order)
                
                
                # Add the created Bond object to the bonds dict, keyed by the
                # normalized (smaller, larger) atom-index pair so that a lookup
                # works regardless of the order the two atoms are given.
                bkey = (index_i, index_j) if index_i <= index_j else (index_j, index_i)
                self.bonds[bkey] = bond
                
                # Update the atoms with the created Bond object, indicating their bond connections
                self.atoms[index_i].bonds.append(bond)
                self.atoms[index_j].bonds.append(bond)
        
        # Convert the index_bonds list to a numpy array of unsigned 32-bit integers
        self.index_bonds = new_index_bonds
        #print(self.index_bonds)
        self.index_bonds     = np.array(self.index_bonds, dtype=np.uint32)
        self.bond_order_list = np.array(self.bond_order_list, dtype=np.uint32)

        # Percepcao robusta de ordem de ligacao (estilo Antechamber/Wang&Case).
        # So roda quando NAO ha ordens externas confiaveis (ex.: XYZ puro): nesse
        # caso substitui o palpite geometrico (UFF) por uma atribuicao baseada
        # em valencia + modelos de oxianion. Se o arquivo trouxe ordens
        # (MOL2/pDynamo via external_orders), respeita-as.
        if external_orders is None:
            self._perceive_bond_orders()

        self._build_bond_order_per_atom() # bachega 2026/Jun/24
        self._detect_metal_bonds()


    def _detect_metal_bonds(self):
        """ Classifica as ligacoes em metalicas (>=1 metal) e covalentes,
            preenchendo self.metal_bonds e self.covalent_bonds (arrays achatados,
            mesmo formato de index_bonds). NAO altera index_bonds nem bonds: o
            grafo permanece intacto; isto e so roteamento de visualizacao.

            Em caso de falha (modulo ausente etc.), trata tudo como covalente
            para nunca quebrar o carregamento.
        """
        
        try:
            try:
                from vismol.core.metal_bonds import split_metal_bonds
            except Exception:
                try:
                    from core.metal_bonds import split_metal_bonds
                except Exception:
                    from metal_bonds import split_metal_bonds
            symbols = [self.atoms[i].symbol for i in range(len(self.atoms))]
            cov, met = split_metal_bonds(symbols, self.index_bonds)
            self.covalent_bonds = cov
            self.metal_bonds = met
            if len(met):
                logger.info("Metal bonds detected: %d", len(met) // 2)
        except Exception as e:
            logger.warning("Metal bond detection failed (%s); "
                           "treating all bonds as covalent." % e)
            self.covalent_bonds = self.index_bonds
            self.metal_bonds = np.array([], dtype=np.uint32)

    def _perceive_bond_orders(self):
        """ Preenche self.bond_order_list (alinhado a self.index_bonds) usando
            percepcao robusta estilo Antechamber (modulo bond_order_perception).

            Usado quando NAO ha ordens externas confiaveis (ex.: XYZ puro). Em
            caso de qualquer falha (modulo ausente, elemento desconhecido, etc.)
            mantem o bond_order_list atual e apenas avisa, para nunca quebrar o
            carregamento da estrutura.
        """
        # Flag global: se multiple_bonds estiver desligada, nao percebe ordens
        # (deixa tudo simples). Evita ate o custo de CPU do algoritmo.
        try:
            if not self.vm_config.gl_parameters.get("multiple_bonds", True):
                return
        except Exception:
            pass
        try:
            from vismol.core.bond_order_perception import perceive_for_vismol
        except Exception:
            try:
                from core.bond_order_perception import perceive_for_vismol
            except Exception:
                try:
                    from bond_order_perception import perceive_for_vismol
                except Exception as e:
                    logger.warning("bond_order_perception not available (%s); "
                                   "keeping the current orders." % e)
                    return
        try:
            symbols = [self.atoms[i].symbol for i in range(len(self.atoms))]
            orders = perceive_for_vismol(symbols, self.index_bonds)
            self.bond_order_list = np.array(orders, dtype=np.uint32)
        except Exception as e:
            logger.warning("Bond order perception failed (%s); "
                           "keeping the current orders." % e)

    def _build_bond_order_per_atom(self):
        """
        Constrói self.bond_order_per_atom: um inteiro por átomo, alinhado com
        self.frames[x] (mesma ordem das coordenadas). É ESSE array que vira VBO
        e alimenta `vert_bond_order` nos shaders de sticks/linhas.

        Regra: cada átomo recebe a maior ordem dentre os bonds em que aparece.
        Átomos sem bond ficam com 1.
        """
        n_atoms = len(self.atoms)
        order_per_atom = np.ones(n_atoms, dtype=np.int32)  # default = 1
        if self.bond_order_list is None or self.index_bonds is None:
            self.bond_order_per_atom = order_per_atom
            return
        for k, bond_order in enumerate(self.bond_order_list):
            idx_i = int(self.index_bonds[2 * k])
            idx_j = int(self.index_bonds[2 * k + 1])
            o = int(bond_order)
            if idx_i < n_atoms and o > order_per_atom[idx_i]:
                order_per_atom[idx_i] = o
            if idx_j < n_atoms and o > order_per_atom[idx_j]:
                order_per_atom[idx_j] = o
        self.bond_order_per_atom = order_per_atom

    def get_bond(self, index_i, index_j):
        """ Retorna o objeto Bond entre os atomos index_i e index_j em O(1),
            ou None se nao houver ligacao. A ordem dos argumentos nao importa
            (a chave e normalizada para (menor, maior)). """
        if self.bonds is None:
            return None
        key = (index_i, index_j) if index_i <= index_j else (index_j, index_i)
        return self.bonds.get(key)

    def perceive_bond_order_for_pairs(self, flat_pairs):
        """
        Heuristica de valencia (GABEDIT_MAX_VALENCE), igual a usada em
        _bonds_from_pair_of_indexes_list, mas como funcao PURA: dado um
        array achatado de pares de atomos [i0,j0,i1,j1,...], devolve um
        array numpy (uint32) com uma ordem de ligacao por par, na MESMA
        ordem dos pares recebidos.

        Diferenca-chave em relacao ao codigo original: o grau de cada atomo
        usado no teste de valencia e' calculado LOCALMENTE a partir apenas
        dos pares recebidos (dict 'degree' interno) -- NAO mexe em
        self.atoms[i].nbonds nem em nenhum outro estado persistente do
        objeto. Isso permite chamar esta funcao repetidamente com listas de
        pares DIFERENTES (ex.: a conectividade de cada frame das Dynamic
        Bonds, que pode mudar a cada passo da trajetoria) sem que uma
        chamada "contamine" a proxima com graus inflados de uma chamada
        anterior.

        [EN] HISTORICO / BUG FIX (aneis conjugados): a atribuicao de duplas
        era feita com um GULOSO de uma passada so' (primeira ligacao
        "candidata" -- ambos atomos com folga -- na ordem em que aparece no
        array e' promovida). Isso e', na pratica, uma selecao gulosa de
        arestas para um problema de CASAMENTO MAXIMO -- e selecao gulosa de
        casamento e' conhecida por depender da ordem de visita. Na pratica:
        o MESMO 1,3-butadieno (CH2=CH-CH=CH2) dava a estrutura certa (duplas
        nas pontas) se a ligacao central aparecesse por ultimo no arquivo,
        mas dava uma estrutura QUIMICAMENTE INVALIDA (dupla so' no meio,
        atomo de ponta com valencia faltando) se a ligacao central aparecesse
        primeiro -- mesma molecula, mesma topologia, resposta diferente e as
        vezes errada, so' por causa da ordem de escrita do arquivo/parser.
        Aneis fundidos/aromaticos com heteroatomo (ex. imidazol em
        histidina/purinas) tinham o mesmo problema.

        Corrigido delegando a decisao de quais ligacoes promover a dupla
        para um CASAMENTO MAXIMO EXATO por componente conexo (ver
        bond_order_perception.perceive_bond_order_for_pairs_pure) -- da'
        sempre a mesma resposta independente da ordem de entrada, e resolve
        aneis conjugados corretamente. Continua rapido porque o problema so'
        e' resolvido dentro dos poucos atomos realmente conjugados (o resto
        da molecula, ex. cadeia principal sp3 de uma proteina, nem entra no
        grafo de candidatas).

        [EN] HISTORICO / BUG FIX (fronteira QC/MM): quando flat_pairs e' um
        SUBCONJUNTO das ligacoes reais de um atomo -- caso das Dynamic Bonds,
        onde find_bonded_and_nonbonded_atoms monta o grid so' com os atomos
        da selecao/regiao QC -- um atomo de fronteira (ligado tambem a um
        atomo da regiao MM) tem essa ligacao para a MM "invisivel" aqui: ela
        nunca aparece em flat_pairs. Sem correcao, o grau local desse atomo
        fica subestimado em 1 (ou mais), e o algoritmo acha que ele ainda
        tem folga de valencia que na verdade ja foi consumida pela ligacao
        para a regiao MM -- promovendo indevidamente uma ligacao vizinha
        (dentro da regiao QC) a dupla perto da fronteira.

        Corrigido comparando, para cada atomo presente em flat_pairs, o grau
        local (dentro deste subconjunto) com self.atoms[i].nbonds -- a
        contagem REAL de ligacoes desse atomo, calculada uma unica vez sobre
        a estrutura estatica completa (QC+MM juntos) no carregamento do
        sistema, antes de qualquer selecao Dynamic Bonds existir. A diferenca
        (nbonds real - grau local aqui) e' passada para
        perceive_bond_order_for_pairs_pure via extra_degree, que "pre-
        consome" essa valencia antes da primeira passada.

        Mesma regra de duas passadas de sempre (dupla primeiro, tripla so'
        promove quem ja e' dupla).
        """
        ib = np.asarray(flat_pairs).ravel()
        n_bonds = int(ib.shape[0] // 2)
        if n_bonds == 0:
            return np.ones(0, dtype=np.uint32)

        symbols = [self.atoms[i].symbol for i in range(len(self.atoms))]

        # Grau local dentro APENAS deste subconjunto de pares (pode ser
        # menor que o grau real do atomo, se flat_pairs for um recorte --
        # ex.: Dynamic Bonds restritas a' regiao QC).
        local_degree = {}
        for k in range(n_bonds):
            i = int(ib[2 * k]); j = int(ib[2 * k + 1])
            local_degree[i] = local_degree.get(i, 0) + 1
            local_degree[j] = local_degree.get(j, 0) + 1

        # Diferenca entre o grau REAL (estrutura estatica completa,
        # self.atoms[i].nbonds) e o grau visivel neste subconjunto -- so'
        # atomos de fronteira (com ligacoes fora do subconjunto, tipicamente
        # para a regiao MM) tem diferenca > 0 aqui.
        extra_degree = {}
        for atom_idx in local_degree:
            true_n = getattr(self.atoms[atom_idx], "nbonds", None)
            if true_n is not None:
                missing = int(true_n) - local_degree[atom_idx]
                if missing > 0:
                    extra_degree[atom_idx] = missing

        try:
            from vismol.core.bond_order_perception import perceive_bond_order_for_pairs_pure
        except Exception:
            try:
                from core.bond_order_perception import perceive_bond_order_for_pairs_pure
            except Exception:
                try:
                    from bond_order_perception import perceive_bond_order_for_pairs_pure
                except Exception as e:
                    logger.warning("bond_order_perception (maximum matching) "
                                   "not available (%s); using the local greedy "
                                   "fallback." % e)
                    return self._perceive_bond_order_for_pairs_greedy_fallback(
                        ib, symbols, extra_degree)

        order = perceive_bond_order_for_pairs_pure(symbols, ib, GABEDIT_MAX_VALENCE,
                                                     extra_degree=extra_degree)
        return np.asarray(order, dtype=np.uint32)

    def _perceive_bond_order_for_pairs_greedy_fallback(self, ib, symbols, extra_degree=None):
        """ Fallback local -- o algoritmo guloso de uma passada so' (com o
            bug de case do GABEDIT_MAX_VALENCE ja corrigido e a correcao de
            fronteira QC/MM via extra_degree), usado APENAS se o modulo
            bond_order_perception nao puder ser importado por algum motivo.
            Mantido para nunca quebrar o carregamento da estrutura; sabidamente
            sujeito ao problema de ordem descrito em
            perceive_bond_order_for_pairs (nao usar como caminho principal).
        """
        n_bonds = int(ib.shape[0] // 2)
        degree = {}
        for k in range(n_bonds):
            i = int(ib[2 * k]); j = int(ib[2 * k + 1])
            degree[i] = degree.get(i, 0) + 1
            degree[j] = degree.get(j, 0) + 1

        if extra_degree:
            for atom, extra in extra_degree.items():
                if atom in degree and extra > 0:
                    degree[atom] += extra

        order = np.ones(n_bonds, dtype=np.uint32)

        for k in range(n_bonds):
            i = int(ib[2 * k]); j = int(ib[2 * k + 1])
            max_i = GABEDIT_MAX_VALENCE.get(symbols[i], 4)
            max_j = GABEDIT_MAX_VALENCE.get(symbols[j], 4)
            if degree[i] < max_i and degree[j] < max_j:
                order[k] = 2
                degree[i] += 1
                degree[j] += 1

        for k in range(n_bonds):
            if order[k] != 2:
                continue
            i = int(ib[2 * k]); j = int(ib[2 * k + 1])
            max_i = GABEDIT_MAX_VALENCE.get(symbols[i], 4)
            max_j = GABEDIT_MAX_VALENCE.get(symbols[j], 4)
            if degree[i] < max_i and degree[j] < max_j:
                order[k] = 3
                degree[i] += 1
                degree[j] += 1

        return order

    def get_dynamic_bond_order_for_frame(self, f):
        """
        Retorna a ordem de ligacao para o frame 'f' das Dynamic Bonds,
        pareada 1:1 com self.dynamic_bonds[f] (mesma ordem de pares).
        Calcula sob demanda e cacheia em self.dynamic_bond_orders[f]; so'
        recalcula se o conteudo dos pares mudou (nao so' o tamanho -- em
        uma regiao QC com ligacoes quebrando/formando, o NUMERO de ligacoes
        pode ficar igual enquanto os PARES mudam, ex.: A-B quebra e C-D
        forma ao mesmo tempo; comparar so' o tamanho deixaria passar esse
        caso e reusaria a ordem errada).

        [NOVO] Depois da percepcao automatica, aplica
        self.dynamic_manual_bond_orders[f] (se houver) por cima -- os
        pares que o usuario forcou explicitamente (via 'bond ...
        frame=...', ver atom_ops.set_dynamic_bond_order()) sempre vencem
        a percepcao automatica para aquele par NAQUELE frame especifico.
        Os setters (set_dynamic_bond_order/unset_dynamic_bond) sao
        responsaveis por invalidar self.dynamic_bond_orders[f] (setar
        para None) sempre que mexerem no override ou nos pares daquele
        frame, para que este metodo recalcule (aplicando o override de
        novo) na proxima chamada -- ver essas funcoes em atom_ops.py.
        """
        if self.dynamic_bonds is None or f < 0 or f >= len(self.dynamic_bonds):
            return np.ones(0, dtype=np.uint32)
        while len(self.dynamic_bond_orders) <= f:
            self.dynamic_bond_orders.append(None)
        pairs = np.asarray(self.dynamic_bonds[f]).ravel()
        cached = self.dynamic_bond_orders[f]  # None ou (pairs_snapshot, order)
        if cached is not None:
            cached_pairs, cached_order = cached
            if cached_pairs.shape == pairs.shape and np.array_equal(cached_pairs, pairs):
                return cached_order
        order = self.perceive_bond_order_for_pairs(pairs)

        dmb = getattr(self, "dynamic_manual_bond_orders", None)
        overrides = dmb.get(f) if dmb else None
        if overrides:
            n_bonds = int(pairs.shape[0] // 2)
            for k in range(n_bonds):
                key = (int(pairs[2 * k]), int(pairs[2 * k + 1]))
                key = (min(key), max(key))
                if key in overrides:
                    order[k] = overrides[key]

        self.dynamic_bond_orders[f] = (pairs.copy(), order)
        return order

    def _bonds_from_pair_of_indexes_list(self, exclude_list=[['H', 'H']],
                                         external_orders=None):
        """
        Creates Bond objects based on pairs of indexes in self.index_bonds.

        self.index_bonds = [0,1 , 0,4 , 1,3 , ...]   (achatado, pares)
        self.bonds       = {(i,j): bond, ...}  # dict, chave = par de indices

        external_orders (Convenção B):
            Lista de ordens de ligação alinhada APENAS aos bonds que
            SOBREVIVEM ao filtro de exclusão (exclude_list), na mesma ordem
            em que eles aparecem. Ou seja, external_orders[k] é a ordem do
            k-ésimo bond NÃO-excluído. Quem monta essa lista (ex.: o parser
            de MOL2/SDF ou o wrapper do pDynamo) deve aplicar a MESMA regra
            de exclusão para manter o alinhamento.

            Se None, a ordem é estimada por distância (UFF) via
            bond.get_bond_order(), preservando o comportamento antigo.
        """
        assert self.bonds is None  # garante que a lista ainda não foi inicializada
        self.bonds = {}            # dict {(i,j): Bond}, chave = par normalizado
        self.bond_order_list = []  # ordem de cada bond (paralela a self.bonds)
        new_index_bonds = []

        # [EN] BUG FIX (found via the Builder tool, but pre-existing in
        # this method regardless -- normal single-shot file loading never
        # triggered it, since bonds were only ever computed once per
        # object's lifetime; the Builder's add_atom()/remove_atom()/
        # add_bond() call this repeatedly as atoms are added/removed one
        # at a time, which exposed it): the loop below does
        # self.atoms[index_i].bonds.append(bond) for every bond, but
        # nothing was ever resetting each ATOM's OWN .bonds list first --
        # only the OBJECT-level self.bonds dict got reset just above.
        # Confirmed live: after 3 add_atom() calls plus one add_bond()
        # call (4 total recomputations) on the same small object, atom 0
        # had accumulated 6 stale/duplicate Bond objects in its .bonds
        # list instead of the correct 2. Fixed by clearing every atom's
        # .bonds list here, once, before the loop that repopulates it.
        for atom in self.atoms.values ( ):
            atom.bonds = []

        # Contador de bonds NÃO-excluídos. Só avança quando um bond é
        # efetivamente criado, mantendo o alinhamento com external_orders
        # (Convenção B). Diferente de i//2, que contaria também os excluídos.
        bond_pair_idx = 0

        # Percorre self.index_bonds de 2 em 2 (cada par = um bond)
        n = 0
        
        
        #''' este loop é para montar o self.bonds dict com a ligaçãoes'''
        for i in range(0, len(self.index_bonds) - 1, 2):
            index_i = self.index_bonds[i]      # índice do primeiro átomo
            index_j = self.index_bonds[i + 1]  # índice do segundo átomo

            
            # Verifica se o par está na lista de exclusão (ex.: H–H)
            is_excluded = False
            for excluded_bond in exclude_list:
                if (self.atoms[index_i].symbol in excluded_bond and
                        self.atoms[index_j].symbol in excluded_bond):
                    is_excluded = True

            if is_excluded:
                # Bond excluído: não cria objeto, não consome posição em
                # external_orders e NÃO incrementa bond_pair_idx.
                pass
            
            else:
                new_index_bonds.append(index_i)
                new_index_bonds.append(index_j)

                # Cria o objeto Bond com os átomos e seus índices
                bond = Bond(atom_i=self.atoms[index_i], atom_index_i=index_i,
                            atom_j=self.atoms[index_j], atom_index_j=index_j)

                # [EN] BUG FIX: external_orders was accepted as a parameter
                # and documented above ("Convenção B" -- external_orders[k]
                # is the order of the k-th NON-EXCLUDED bond) but was NEVER
                # ACTUALLY CONSULTED anywhere in this method: every bond
                # silently got bond.bond_order = 1 here, then got
                # UNCONDITIONALLY overwritten again a few lines below by
                # perceive_bond_order_for_pairs() -- regardless of what the
                # caller passed in external_orders. This silently discarded
                # every manual bond-order override (Builder's terminal
                # 'bond order=N' / Ctrl+click cycle_bond_order()): the
                # command ran with no error, manual_bond_orders was
                # correctly updated, external_orders was correctly built by
                # _reapply_manual_bonds() in atom_ops.py -- but by the time
                # a Bond object actually existed, its order had already been
                # replaced by whatever the automatic valence heuristic
                # decided. Fixed by applying external_orders[bond_pair_idx]
                # HERE, 1:1 with the k-th non-excluded bond exactly as the
                # docstring already specified, and by only running the
                # automatic perception below when external_orders is None
                # (see that block's own updated comment). The caller is
                # responsible for MERGING automatic perception with any
                # manual overrides before calling this method with
                # external_orders -- see _reapply_manual_bonds()'s own
                # updated comment for how it now does that.
                if external_orders is not None and bond_pair_idx < len(external_orders):
                    bond.bond_order = int(external_orders[bond_pair_idx])
                else:
                    bond.bond_order = 1  # provisorio -- sobrescrito abaixo se external_orders for None

                # Avança o contador só agora, após criar um bond válido.
                bond_pair_idx += 1

                # Registra o bond no dict, chave = par de indices normalizado
                # (menor, maior), para busca independente da ordem dos atomos.
                bkey = (index_i, index_j) if index_i <= index_j else (index_j, index_i)
                
                if bkey in self.bonds.keys():
                    pass
                else:
                    self.bonds[bkey] = bond
                    self.atoms[index_i].bonds.append(bond)
                    self.atoms[index_j].bonds.append(bond)
        
        
        #atribuindo aos átomos o numero de parceiros
        for atom in self.atoms.values ( ):
            atom.nbonds = len(atom.bonds)
            #print(atom.nbonds)
        
    
        #print(new_index_bonds)
        self._index_bonds_from_bonds_dict()
        
        #print(self.index_bonds)
        self.index_bonds = np.array(self.index_bonds, dtype=np.uint32)
       
        #print(self.bond_order_list)
        #print(self.bonds.keys())
        
        
        
        
        
        #- - - - - - - - - percepcao de ordem (dupla/tripla) - - - - - - - -
        # Reusa a mesma heuristica de valencia usada pelas Dynamic Bonds
        # (ver perceive_bond_order_for_pairs), agora como funcao pura: nao
        # muta mais self.atoms[i].nbonds como efeito colateral (antes os
        # dois loops separados incrementavam isso durante o calculo, o que
        # inflava .nbonds alem do numero real de vizinhos do atomo).
        #
        # [EN] BUG FIX: this used to run UNCONDITIONALLY, overwriting
        # every bond.bond_order set above from external_orders (see that
        # loop's own updated comment). Now only runs when external_orders
        # is None -- i.e. the normal file-loading path (no externally
        # known/trusted orders), exactly what the docstring already said
        # was supposed to happen ("Se None, a ordem é estimada..."). When
        # external_orders IS given, bond.bond_order was already set
        # correctly per-bond in the loop above -- just read it back into
        # bond_order_list here for anything that consults that list
        # separately.
        if external_orders is None:
            computed_orders = self.perceive_bond_order_for_pairs(self.index_bonds)
            # [EN] 2026-10-15: residues WITHOUT hydrogens (raw PDB files) --
            # the valence heuristic above assumes every H is present and puts
            # spurious double bonds there. Standard residues get their
            # template, other residues EasyHybrid's native perception (from
            # the heavy-atom geometry); see residue_bond_orders.py.
            try:
                from vismol.core.residue_bond_orders import residue_aware_orders
                computed_orders = residue_aware_orders(self.atoms, self.index_bonds, computed_orders,
                                                       self.frames[0])
            except Exception as e:
                logger.warning("residue-aware bond orders failed (%s); valence orders kept" % e)
            self.bond_order_list = computed_orders.tolist()
            #print(len(self.bond_order_list))
            #print(self.bond_order_list)


            # Sincroniza bond_order_list de volta para os objetos Bond em
            # self.bonds (representations.py le' bond.bond_order direto do
            # objeto via self.get_bond(i,j), nao posicionalmente da lista).
            # Mesma ordem de iteracao usada em _index_bonds_from_bonds_dict
            # (dict nao foi mutado desde entao, entao a correspondencia por
            # posicao continua valida).
            for k, bond in enumerate(self.bonds.values()):
                bond.bond_order = int(self.bond_order_list[k])
        else:
            self.bond_order_list = [int(bond.bond_order) for bond in self.bonds.values()]
        
        
        # Constroi o array de ordem-de-ligacao por atomo (alinhado com as
        # coordenadas) usado pelo VBO 'vert_bond_order' nos sticks/lines.
        
        #self._build_bond_order_per_atom()
        
        
        
        '''----------------------------------------------------------'''
        #self._detect_metal_bonds()
        '''----------------------------------------------------------'''
        
        
        #print(self.index_bonds)
        #print(external_orders)
        
        # Propaga a ordem (por aresta) para um array por átomo, alinhado
        # com as coordenadas — é esse array que alimenta o VBO do shader.
        #self._build_bond_order_per_atom()

    def _index_bonds_from_bonds_dict (self):
        """ 
        Esta função  gera  as listas  self.index_bonds e a self.bonds a 
        partir do dicionário self.bonds. Ela garante que não há ligações 
        duplicadas no self.index_bonds
        
        
        self.index_bonds = [i0, j0, i1, j1, ...] — two elements per bond (atom i, atom j).
        self.bond_order_list = [order0, order1, ...]. - all int s 
        
        self.bonds = {(i0,j0):Bond, (i1,j1):Bond,, (i2,j2):Bond, } / Bond is a bond object
        
        """
        self.index_bonds = [x for par in self.bonds.keys() for x in par]
        self.bond_order_list = []
        for bond in self.bonds.values():
            self.bond_order_list.append(bond.bond_order)
        #self.bond_order_list(list(self.bonds.values()))
    
    
    def _get_non_bonded_from_bonded_list(self):
        """ Function doc """
        assert self.non_bonded_atoms is None
        bonded_set = set(self.index_bonds)
        self.non_bonded_atoms = []
        for i, atom in self.atoms.items():
            if i in bonded_set:
                atom.nonbonded = False
            else:
                atom.nonbonded = True
                self.non_bonded_atoms.append(i)
        self.non_bonded_atoms.sort()
        self.non_bonded_atoms = np.array(self.non_bonded_atoms, dtype=np.int32)
    
    def _get_center_of_mass(self, frame=0):
        """ Function doc """
        if frame >= self.frames.shape[0]:
            logger.info("Frame {} is out of range for trajectory of size {}. \
                Using the last frame.".format(frame, self.frames.shape[0] - 1))
            frame = self.frames.shape[0] - 1
        return np.mean(self.frames[frame], axis=0)
    
    def set_model_matrix(self, mat):
        """ Function doc """
        self.model_mat = np.copy(mat)

    def pin_frame(self, frame_index):
        """ Pins this object's rendering/coordinate lookups to a fixed
            frame, overriding the session-global vm_session.frame for
            THIS object only (see self._frame_override, set in
            __init__, and its use in VismolGLCore._safe_frame_coords/
            _get_vismol_object_frame and Atom.coords()). Every other
            object in the scene keeps following the global trajectory
            scrubber normally. Call unpin_frame() to release it.
        """
        self._frame_override = int(frame_index)

    def unpin_frame(self):
        """ Releases a pin set by pin_frame(), returning this object to
            following the session-global vm_session.frame again.
        """
        self._frame_override = None
    
    def _generate_topology_from_index_bonds(self, bonds = None):
        """ 
        bonds = [92,93  ,  92,99  ,  ...]
        
        Returns a graph in dictionary form (this may in turn be 
        needed to determine which objects are molecules).
        
        defines: self.topology = {92: [93, 99], 93: [92, 94, 96, 100], 99: [92],...}
        
        """
        if bonds is None:
            bonds =  self.index_bonds
        else:
            pass

        
        bonds_pairs = []
        topology    = {}
        
        for i  in range(0, len(bonds),2):
            bonds_pairs.append([bonds[i], bonds[i+1]])
        
        
        for bond in bonds_pairs:
            if bond[0] in topology.keys():
                topology[bond[0]].append(bond[1])
            else:
                topology[bond[0]] = []
                topology[bond[0]].append(bond[1])
            
            
            if bond[1] in topology.keys():
                topology[bond[1]].append(bond[0])
            else:
                topology[bond[1]] = []
                topology[bond[1]].append(bond[0])
        self.topology = topology
        #self.find_rings(topology)
        
        
    def find_rings(self, graph):
        '''
        The algorithm traverses all paths using DFS (depth-first search).
        A ring is identified when we return to the starting node after passing through at least two other nodes (len(path) > 2).
        After the search, it removes duplicates based on the set of nodes in each ring.
        Expected output:

        csharp
        Copy
        Edit
        Cycles found:
        [92, 93, 99]
        [93, 96, 100]
        '''
        rings = []

        def dfs(current, start, path, visited):
            path.append(current)
            visited.add(current)

            for neighbor in graph[current]:
                if neighbor == start and len(path) > 2:
                    ring = sorted(path)
                    if ring not in rings:
                        rings.append(ring[:])
                elif neighbor not in visited:
                    dfs(neighbor, start, path, visited)

            path.pop()
            visited.remove(current)

        for node in graph:
            dfs(node, node, [], set())

        # Eliminate duplicates (same cycle in different orders)
        unique_rings = []
        for ring in rings:
            if not any(set(ring) == set(c) for c in unique_rings):
                unique_rings.append(ring)
        self.rings = unique_rings
        return unique_rings

    def define_molecules (self):
        """ Function doc 
        self.topology  = It is a graph written in the form of a dictionary:
                        {
                         index1 : [index2 , index3, ...], 
                         index2 : [index1 , index4, ...]
                         ...}
        
        groups = [{0, 1, 2}, {3, 4, 5, 6}, {8, 9, 7}]
        
        self.atoms =  {
                       0: <pdynamo.pDynamo2EasyHybrid.Atom object at 0x7f3323cbcdf0>, 
                       1: <pdynamo.pDynamo2EasyHybrid.Atom object at 0x7f3323cbcfd0>, 
                       2: <pdynamo.pDynamo2EasyHybrid.Atom object at 0x7f3323cbcf40>, 
                       3: <pdynamo.pDynamo2EasyHybrid.Atom object at 0x7f3323e38250>, 
                       4: <pdynamo.pDynamo2EasyHybrid.Atom object at 0x7f3323e38040>, 
                       5: ...
                       }
        
        This function populates the molecule dictionary (self.molecules), each 
        molecule object is created very similarly to the residue object
        
        """
        #try:
        #groups = find_groups(self.topology)
        groups = find_connected_components(self.topology)
        mol_index = 0
        for mol_index, group in enumerate( groups ):
            atoms = {}
            molecule = Molecule(self, name="UNK", index = mol_index)
            
            for atom_index in list(group):
                molecule.atoms[atom_index] = self.atoms[atom_index]
                self.atoms[atom_index].molecule = molecule
                
            self.molecules[mol_index] = molecule
        
        #-------------------------------------------------------------
        # non_bonded_atoms
        # Should be here, ohterwise we will have selection problems
        #-------------------------------------------------------------

        for atom_index in self.non_bonded_atoms:
            #try:
            mol_index += 1
            molecule = Molecule(self, name="UNK", index = mol_index)
            molecule.atoms[atom_index] = self.atoms[atom_index]
            self.atoms[atom_index].molecule = molecule
            self.molecules[mol_index] = molecule

            #print(self.molecules)
        #except:
        #    print('Failure to determine the list of molecules!')
            
    def define_Calpha_backbone (self):
        """ Function doc 
        Verifica quais conexões entre c_alphas são válidas.
        """
        
        
        self.c_alpha_bonds = []
        self.c_alpha_atoms = []
        
        #
        # Building the self.c_alpha_atoms dict
        # {atom_id : atom_object, ...}
        #
        for c_index, chain in self.chains.items():
            for r_index, residue in chain.residues.items():
                residue._is_protein()
                if residue.is_protein:
                    for a_index, atom in residue.atoms.items():
                        if atom.name == "CA":
                            self.c_alpha_atoms.append(atom)

        
        for i in range(1, len(self.c_alpha_atoms)):

            atom_before  = self.c_alpha_atoms[i-1]
            resi_before  = atom_before.residue.index
            index_before = atom_before.atom_id
            

            atom  = self.c_alpha_atoms[i]
            resi  = atom.residue.index
            index = atom.atom_id
            
            '''Checks whether the two residues are in sequence 
            (otherwise there is a break in the backbone structure)'''
            if resi == resi_before + 1:
                
                bond = Bond(atom_i=atom_before, atom_index_i=index_before,
                            atom_j=atom, atom_index_j=index)
                
                distance = bond.distance()
                if distance < 4.0:
                    self.c_alpha_bonds.append(bond)
        
        #print(self.c_alpha_bonds)
        #print(self.c_alpha_atoms)
        

    def _calculate_unit_cell_vertices(self, a, b, c, alpha, beta, gamma):
        '''
        Returns the list of vertices positions of 
        a box with parameters a, b, c, alpha, beta, gamma.
        '''
        
        # Convert angles to radians
        alpha_rad = np.radians(alpha)
        beta_rad  = np.radians(beta)
        gamma_rad = np.radians(gamma)

        # Calculate unit cell vectors
        v1 = np.array([a, 0, 0])
        v2 = np.array([b * np.cos(gamma_rad), b * np.sin(gamma_rad), 0])
        v3_x = c * np.cos(beta_rad)
        v3_y = c * (np.cos(alpha_rad) - np.cos(beta_rad) * np.cos(gamma_rad)) / np.sin(gamma_rad)
        v3_z = np.sqrt(c**2 - v3_x**2 - v3_y**2)
        v3 = np.array([v3_x, v3_y, v3_z])

        # Define the eight vertices of the unit cell
        vertices = [
            np.array([0, 0, 0]),
            v1,
            v2,
            v2 + v1,
            v3,
            v3 + v1,
            v3 + v2,
            v3 + v2 + v1
        ]

        return vertices
    
    def set_cell (self, a, b, c, alpha, beta, gamma, color = [0.5, 0.5, 0.5]):
        """
        Assign the cell parameters to the 
        "vismol_object" object.
        
        Arguments a, b, c, alpha, beta, gamma 
        must be obtained externally.
        
        """
        
        self.cell_parameters = {'a'     : a, 
                                'b'     : b, 
                                'c'     : c, 
                                'alpha' : alpha, 
                                'beta'  : beta, 
                                'gamma' : gamma
                               }
         
        vertices = self._calculate_unit_cell_vertices(a, b, c, alpha, beta, gamma)


        '''
        The coordinates follow the same structure as the atoms, 
        they are organized in a trajectory structure with only one 
        frame.
        '''
        self.cell_coordinates = np.empty([1, 8, 3], dtype=np.float32)
        self.cell_indexes     = []
       
        
        '''
        The connections between the vertices of the boxes 
        are pre-established.
        '''
        self.cell_bonds       = [0,1, 0,2, 2,3, 0,4 , 1,3 , 1,5 , 2,3 , 2,6 , 3,7 , 4,6 , 4,5, 5,7, 6,7]       
        self.cell_bonds       = np.array(self.cell_bonds, dtype=np.uint32)
        
        
        '''
        The colors of the vertices follow the same 
        structure used in atoms
        '''
        self.cell_colors      = np.empty([8, 3], dtype=np.float32)

        for i, vertex in enumerate(vertices, 0):
        
            dprint("Vertex {}: {:7.3f} {:7.3f} {:7.3f}".format(i, vertex[0] ,vertex[1] ,vertex[2]))
            self.cell_indexes.append(i)
            
            self.cell_colors[i] = np.array(color, dtype=np.float32)
            self.cell_coordinates[0,i,:] = vertex[0] ,vertex[1] ,vertex[2]  
            

        dprint('\n\n')
        #print (self.cell_colors)









def find_connected_components(graph):
    """
    
    graph = topology (self.topology)
    
    eg: 

        self.topology = {0: [1, 2, 1], 1: [0, 0], 2: [0], 3: [6, 4, 5], 6: [3], 4: [3], 5: [3], 7: [8, 9], 8: [7], 9: [7]}
    
    
    
    This version of the function also takes a graph represented as a dictionary and 
    returns a list of connected components, where each connected component is a 
    list of nodes.

    The function starts by initializing an empty set called visited to keep track 
    of the nodes that have been visited and an empty list called components to 
    store the connected components.

    The function then iterates over each node in the graph and checks if it has 
    been visited yet. If the node has not been visited, the function starts a DFS 
    from that node.

    The DFS is implemented using a stack-based approach. The function initializes 
    a stack with the starting node and an empty list called component to store the 
    nodes in the connected component.

    The function then enters a loop that continues as long as the stack is not 
    empty. In each iteration of the loop, the function pops a node from the stack 
    and checks if it has been visited yet. If the node has not been visited, the 
    function adds it to the visited set, appends it to the component list, and 
    adds its unvisited neighbors to the stack.

    After the DFS has completed, the function appends the component list to the 
    components list.

    Finally, the function returns the components list, which contains the connected 
    components of the graph. Each connected component is represented as a list of 
    nodes.
    
    eg:
        components = [[0, 1, 2], [3, 5, 4, 6], [7, 9, 8]]
    
    """
    
    
    visited = set()
    components = []

    for start_node in graph:
        if start_node not in visited:
            stack = [start_node]
            component = []
            while stack:
                node = stack.pop()
                if node not in visited:
                    visited.add(node)
                    component.append(node)
                    stack.extend([neighbor for neighbor in graph[node] if neighbor not in visited])
            components.append(component)
    return components

def DFS(graph, node, visited):
    ''' 
        The DFS function takes a graph, a node, and a set of 
        visited nodes as inputs, and performs a depth-first 
        search starting from the node. 
        
        is not being used

    '''
    visited.add(node)
    for neighbor in graph[node]:
        if neighbor not in visited:
            DFS(graph, neighbor, visited)

def find_groups(graph):
    '''
    
    is not being used due: 
    
    Python: maximum recursion depth exceeded while calling a Python object

    for more info access: 
    https://stackoverflow.com/questions/6809402/python-maximum-recursion-depth-exceeded-while-calling-a-python-object
        
    '''
    visited = set()
    groups = []
    for node in graph:
        if node not in visited:
            group = set()
            DFS(graph, node, group)
            groups.append(group)
            visited |= group
    return groups
