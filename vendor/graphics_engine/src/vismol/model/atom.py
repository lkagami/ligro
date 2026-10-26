#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
#  atom.py
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

import numpy as np
#from vismol.model.molecular_properties import ATOM_TYPES


class Atom:
    """ Class doc """
    
    def __init__(self, vismol_object, name="Xx", index=None, residue=None,
                 chain=None, pos=None, symbol=None, atom_id=None, color=None,
                 vdw_rad=None, cov_rad=None, ball_rad=None, electronegativity=None,
                 occupancy=0.0, bfactor=0.0, charge=0.0, formal_charge=0, bonds_indexes=None):
        """ Class initializer """
        self.vm_object   = vismol_object
        self.vm_session  = vismol_object.vm_session
        
        
        self.name     = name
        self.index    = index   # - Remember that the "index" attribute refers to the numbering of atoms (it is not a zero base, it starts at 1 for the first atom)
        self.residue  = residue
        self.chain    = chain
        self.molecule = None

        self.pos = pos     # - coordinates of the first frame
        self.unique_id = None
        
        # self.symbol = self._get_symbol(self.name) if (symbol is None) else symbol
        if symbol is None:
            self.symbol = self._get_symbol()
        else:
            self.symbol = symbol
        self.atom_id = atom_id
        
        if color is None:
            self.color = self._init_color()
        else:
            self.color = color
        
        if vdw_rad is None:
            self.vdw_rad = self._init_vdw_rad()
        else:
            self.vdw_rad = vdw_rad
        
        if cov_rad is None:
            self.cov_rad = self._init_cov_rad()
        else:
            self.cov_rad = cov_rad
        
        if ball_rad is None:
            self.ball_rad = self._init_ball_rad()
        else:
            self.ball_rad = ball_rad
        
        #Bachega April/05/2026
        if electronegativity is None:
            self.electronegativity = self._init_electronegativity()
        else:
            self.electronegativity = electronegativity           
        
        self.color_id = None
        self.occupancy = occupancy
        self.bfactor = bfactor
        self.charge = charge
        # [EN] 2026-09-25: a SEPARATE field from `charge` above -- `charge`
        # is the MM/QC PARTIAL charge (a float, populated from `system.
        # mmState.charges` once a force field is assigned; see
        # pDynamo2EasyHybrid/session.py's define_a_new_QCModel()).
        # `formal_charge` is the INTEGER formal/protonation charge (0 by
        # default) the Builder's own atom_ops.protonate_atom()/
        # deprotonate_atom() read and write -- kept distinct so assigning
        # an MM force field never silently overwrites a user-set
        # protonation state, and vice versa.
        self.formal_charge = formal_charge
        if bonds_indexes is None:
            self.bonds_indexes = []
        else:
            self.bonds_indexes = bonds_indexes
        
        self.selected       = False
        self.lines          = True
        self.dots           = False
        self.nonbonded      = False
        self.impostor       = False
        self.ribbons        = False
        self.ribbon_sphere  = False
        self.dynamic        = False
        self.ball_and_stick = False
        self.sticks         = False
        self.stick_spheres  = False
        self.spheres        = False
        self.vdw_spheres    = False
        self.dash           = False
        self.surface        = False
        # [EN] 2026-09-28 partial cartoon: residue is drawn by the Cartoon
        # representation when its CA atom has this set -- see
        # CartoonRepresentation._residue_is_shown() (representations.py).
        self.cartoon        = False
        self.bonds          = []
        self.isfree         = True
        self.labels         = False
        self.label_text     = None
        
    def _get_symbol(self):
        """ Function doc """
        ATOM_TYPES = self.vm_session.periodic_table.elements_by_symbol
        
        name = self.name.strip()
        if name == "":
            return ""
        
        _n = name
        for char in name:
            if char.isnumeric():
                _n = _n.replace(char, "")
        name = _n
        
        if len(name) >= 3:
            name = name[:2]
        
        # This can't happen since before all numbers have been converted to strings
        # if len(name) == 2:
        #     if name[1].isnumeric():
        #         symbol = name[0]
        
        # The capitalization of hte name can solve all the next ifs
        # name = name.lower().capitalize()
        if name in ATOM_TYPES.keys():
            return name
        else:
            if name[0] == "H":
                if name[1] == "g":
                    symbol =  "Hg"
                elif name[1] =="e":
                    symbol = "He"
                else:
                    symbol =  "H"
            
            elif name[0] == "C":
                if name[1] == "a":
                    symbol = "Ca"
                elif name[1] =="l":
                    symbol = "Cl"
                elif name[1] =="L":
                    symbol = "Cl"
                elif name[1] =="d":
                    symbol = "Cd"
                elif name[1] =="u":
                    symbol = "Cu"
                elif name[1] =="U":
                    symbol = "Cu"
                else:
                    symbol = "C"
            
            elif name[0] == "N":
                if name[1] == "i" or name[1] == "I":
                    symbol = "Ni"
                elif name[1] == "a":
                    symbol = "Na"
                elif name[1] == "e":
                    symbol = "Ne"
                elif name[1] == "b":
                    symbol = "Nb"
                else:
                    symbol = "N"
            
            elif name[0] == "O":
                if name[1] == "s":
                    symbol = "Os"
                else:
                    symbol = "O"
            
            elif name[0] == "S":
                if name[1] == "I":
                    symbol = "Si"
                elif name[1] == "e":
                    symbol = "Se"
                else:
                    symbol = "S"
            
            elif name[0] == "P":
                if name[1] == "d":
                    symbol = "Pd"
                elif  name[1] == "b":
                    symbol = "Pb"
                elif  name[1] == "o":
                    symbol = "Po"
                else:
                    symbol = "P" 
            
            elif name[0] == "Z":
                if name[1] == "r":
                    symbol = "Zr"
                elif  name[1] == "N":
                    symbol = "Zn"
                else:
                    symbol = "Zn" 
            
            elif name[0] == "F":
                if name[1] == "E":
                    symbol = "Fe"
                elif  name[1] == "e":
                    symbol = "Fe"
                else:
                    symbol = "F" 
            
            elif name[0] == "M":
                if name[1] == "n":
                    symbol = "Mn"
                elif name[1] == "N":
                    symbol = "Mn"
                elif name[1] == "o":
                    symbol = "Mo"
                elif name[1] == "G":
                    symbol = "Mg"
                else:
                    symbol = "X"
            else:
                symbol = "X"
        return symbol
    
    def _init_color(self):
        """ Return the color of an atom in RGB. Note that the returned
            value is in scale of 0 to 1, but you can change this in the
            index. If the atomname does not match any of the names
            given, it returns the default dummy value of atom X.
        """
        try:
            color = self.vm_object.color_palette[self.name]
        except KeyError:
            # #print(self.symbol, "Atom")
            color = self.vm_object.color_palette[self.symbol]
        return np.array(color, dtype=np.float32)
    
    def _generate_atom_unique_color_id(self):
        ATOM_TYPES = self.vm_session.periodic_table.elements_by_symbol
        """ Function doc """
        r = (self.unique_id & 0x000000FF) >> 0
        g = (self.unique_id & 0x0000FF00) >> 8
        b = (self.unique_id & 0x00FF0000) >> 16
        # pickedID = r + g * 256 + b * 256*256
        self.color_id = np.array([r/255.0, g/255.0, b/255.0], dtype=np.float32)
    
    def _init_vdw_rad(self):
        """ Function doc """
        ATOM_TYPES = self.vm_session.periodic_table.elements_by_symbol
        try:
            vdw = ATOM_TYPES[self.name][6]
        except KeyError:
            vdw = ATOM_TYPES[self.symbol][6]
        return vdw
    
    def _init_electronegativity(self):
        """ Function doc """
        ATOM_TYPES = self.vm_session.periodic_table.elements_by_symbol
        try:
            en_UFF = ATOM_TYPES[self.name][8]
        except KeyError:
            en_UFF = ATOM_TYPES[self.symbol][8]
        return en_UFF
    
    def _init_cov_rad(self):
        """ Function doc """
        ATOM_TYPES = self.vm_session.periodic_table.elements_by_symbol
        try:
            #cov = ATOM_TYPES[self.name][5]
            cov = ATOM_TYPES[self.name][7]
        except KeyError:
            #cov = ATOM_TYPES[self.symbol][5]
            cov = ATOM_TYPES[self.symbol][7]
        return cov
    
    def _init_ball_rad(self):
        """ Function doc """
        ATOM_TYPES = self.vm_session.periodic_table.elements_by_symbol
        try:
            ball = ATOM_TYPES[self.name][6]
        except KeyError:
            ball = ATOM_TYPES[self.symbol][6]
        return ball
    
    def coords(self, frame=None):
        """ 
        frame = int
        
        Returns the coordinates of an atom according to the specified frame. 
        If no frame is specified, the frame set by easyhybrid (probably by the 
        scale bar of the trajectory manipulation window) is used. If the object 
        (vobject) has a smaller number of frames than the one set by the 
        interface, the last frame of the object is used.
        
        return  xyz 
        """
        if frame is None:
            # [EN] Per-object frame pin (VismolObject.pin_frame()) --
            # same override VismolGLCore._safe_frame_coords()/
            # _get_vismol_object_frame() apply, kept in sync here for
            # non-render callers (distance/dihedral calc, hover labels)
            # that read coordinates directly through this method.
            if self.vm_object._frame_override is not None:
                frame = self.vm_object._frame_override
            else:
                frame = self.vm_object.vm_session.frame
            #print (frame, len(self.vm_object.frames))
            if len(self.vm_object.frames)-1 <= frame:
                frame = len(self.vm_object.frames)-1
            else:
                pass
        return self.vm_object.frames[frame, self.atom_id]
    
    def get_grid_position(self, gridsize=3, frame=None):
        """ Function doc """
        coords = self.coords(frame)
        gridpos = (int(coords[0]/gridsize), int(coords[1]/gridsize), int(coords[2]/gridsize))
        return gridpos
    
    def get_cov_rad(self):
        """ Function doc """
        return self.cov_rad
    
    def init_color_rgb(self, name):
        """ Return the color of an atom in RGB. Note that the returned
            value is in scale of 0 to 1, but you can change this in the
            index. If the atomname does not match any of the names
            given, it returns the default dummy value of atom X.
        """
        
        try:
            color = color =self.vm_object.color_palette[name]
            #color = ATOM_TYPES[name][1]
        except KeyError:
            if name[0] == "H":# or name in self.hydrogen:
                #color = ATOM_TYPES["H"][1]
                color = self.vm_object.color_palette["H"]
            
            elif name[0] == "C":
                #color = ATOM_TYPES["C"][1]
                color = self.vm_object.color_palette["C"]
            
            elif name[0] == "O":
                #color = ATOM_TYPES["O"][1]
                color = self.vm_object.color_palette["O"]
            
            elif name[0] == "N":
                #color = ATOM_TYPES["N"][1]
                color = self.vm_object.color_palette["N"]
                
            elif name[0] == "S":
                #color = ATOM_TYPES["S"][1]
                color = self.vm_object.color_palette["S"]
            else:
                #color = ATOM_TYPES["X"][1]
                color = self.vm_object.color_palette["X"]
                
        color = [int(color[0]*250), int(color[1]*250), int(color[2]*250)]
        return color
    
    def init_radius(self, name):
        """
        """
        ATOM_TYPES = self.vm_session.periodic_table.elements_by_symbol
        try:
            rad = ATOM_TYPES[name][6]/5.0
        except KeyError:
            if name[0] == "H" or name in self.hydrogen:
                rad = ATOM_TYPES["H"][6]/5.0
            elif name[0] == "C":
                rad = ATOM_TYPES["C"][6]/5.0
            elif name[0] == "O":
                rad = ATOM_TYPES["O"][6]/5.0
            elif name[0] == "N":
                rad = ATOM_TYPES["N"][6]/5.0
            elif name[0] == "S":
                rad = ATOM_TYPES["S"][6]/5.0
            else:
                rad = 0.30
        return rad







#class Atom:
#    """ Class doc """
#    
#    def __init__(self, vismol_object=None, name="Xx", index=None, residue=None,
#                 chain=None, pos=None, symbol=None, atom_id=None, color=None,
#                 vdw_rad=None, cov_rad=None, ball_rad=None,
#                 occupancy=0.0, bfactor=0.0, charge=0.0, bonds_indexes=None):
#        
#        """ Class initializer """
#        self.vm_object = vismol_object
#        self.name = name
#        self.index = index   # - Remember that the "index" attribute refers to the numbering of atoms (it is not a zero base, it starts at 1 for the first atom)
#        self.residue = residue
#        self.chain = chain
#        self.molecule = None
#        
#        self.pos = pos     # - coordinates of the first frame
#        self.unique_id = None
#        
#        # self.symbol = self._get_symbol(self.name) if (symbol is None) else symbol
#        if symbol is None:
#            self.symbol = self._get_symbol()
#        else:
#            self.symbol = symbol
#        self.atom_id = atom_id
#        
#        if color is None:
#            self.color = self._init_color()
#        else:
#            self.color = color
#        
#        if vdw_rad is None:
#            self.vdw_rad = self._init_vdw_rad()
#        else:
#            self.vdw_rad = vdw_rad
#        
#        if cov_rad is None:
#            self.cov_rad = self._init_cov_rad()
#        else:
#            self.cov_rad = cov_rad
#        
#        if ball_rad is None:
#            self.ball_rad = self._init_ball_rad()
#        else:
#            self.ball_rad = ball_rad
#        
#        self.color_id = None
#        self.occupancy = occupancy
#        self.bfactor = bfactor
#        self.charge = charge
#        if bonds_indexes is None:
#            self.bonds_indexes = []
#        else:
#            self.bonds_indexes = bonds_indexes
#        
#        self.selected       = False
#        self.lines          = True
#        self.dots           = False
#        self.nonbonded      = False
#        self.impostor       = False
#        self.ribbons        = False
#        self.ribbon_sphere  = False
#        self.ball_and_stick = False
#        self.sticks         = False
#        self.stick_spheres  = False
#        self.spheres        = False
#        self.vdw_spheres    = False
#        self.dash           = False
#        self.surface        = False
#        self.bonds          = []
#        self.isfree         = True
#        self.labels         = False
#        self.label_text     = ''
#                              #{
#                              #'index'     : False,
#                              #'name'      : False,
#                              #'symbol'    : False,
#                              #'MM_type'   : False,
#                              #'MM_charge' : False,
#                              #'res_type'  : False,
#                              #'res_number': False,
#                              #'chain'     : False,
#                              #}
#        
#    def _get_symbol(self):
#        """ Function doc """
#        name = self.name.strip()
#        if name == "":
#            return ""
#        
#        _n = name
#        for char in name:
#            if char.isnumeric():
#                _n = _n.replace(char, "")
#        name = _n
#        
#        if len(name) >= 3:
#            name = name[:2]
#        
#        # This can't happen since before all numbers have been converted to strings
#        # if len(name) == 2:
#        #     if name[1].isnumeric():
#        #         symbol = name[0]
#        
#        # The capitalization of hte name can solve all the next ifs
#        # name = name.lower().capitalize()
#        if name in ATOM_TYPES.keys():
#            return name
#        else:
#            if name[0] == "H":
#                if name[1] == "g":
#                    symbol =  "Hg"
#                elif name[1] =="e":
#                    symbol = "He"
#                else:
#                    symbol =  "H"
#            
#            elif name[0] == "C":
#                if name[1] == "a":
#                    symbol = "Ca"
#                elif name[1] =="l":
#                    symbol = "Cl"
#                elif name[1] =="L":
#                    symbol = "Cl"
#                elif name[1] =="d":
#                    symbol = "Cd"
#                elif name[1] =="u":
#                    symbol = "Cu"
#                elif name[1] =="U":
#                    symbol = "Cu"
#                else:
#                    symbol = "C"
#            
#            elif name[0] == "N":
#                if name[1] == "i" or name[1] == "I":
#                    symbol = "Ni"
#                elif name[1] == "a":
#                    symbol = "Na"
#                elif name[1] == "e":
#                    symbol = "Ne"
#                elif name[1] == "b":
#                    symbol = "Nb"
#                else:
#                    symbol = "N"
#            
#            elif name[0] == "O":
#                if name[1] == "s":
#                    symbol = "Os"
#                else:
#                    symbol = "O"
#            
#            elif name[0] == "S":
#                if name[1] == "I":
#                    symbol = "Si"
#                elif name[1] == "e":
#                    symbol = "Se"
#                else:
#                    symbol = "S"
#            
#            elif name[0] == "P":
#                if name[1] == "d":
#                    symbol = "Pd"
#                elif  name[1] == "b":
#                    symbol = "Pb"
#                elif  name[1] == "o":
#                    symbol = "Po"
#                else:
#                    symbol = "P" 
#            
#            elif name[0] == "Z":
#                if name[1] == "r":
#                    symbol = "Zr"
#                elif  name[1] == "N":
#                    symbol = "Zn"
#                else:
#                    symbol = "Zn" 
#            
#            elif name[0] == "F":
#                if name[1] == "E":
#                    symbol = "Fe"
#                elif  name[1] == "e":
#                    symbol = "Fe"
#                else:
#                    symbol = "F" 
#            
#            elif name[0] == "M":
#                if name[1] == "n":
#                    symbol = "Mn"
#                elif name[1] == "N":
#                    symbol = "Mn"
#                elif name[1] == "o":
#                    symbol = "Mo"
#                elif name[1] == "G":
#                    symbol = "Mg"
#                else:
#                    symbol = "X"
#            else:
#                symbol = "X"
#        return symbol
#    
#    def _init_color(self):
#        """ Return the color of an atom in RGB. Note that the returned
#            value is in scale of 0 to 1, but you can change this in the
#            index. If the atomname does not match any of the names
#            given, it returns the default dummy value of atom X.
#        """
#        try:
#            color = self.vm_object.color_palette[self.name]
#        except KeyError:
#            # print(self.symbol, "Atom")
#            color = self.vm_object.color_palette[self.symbol]
#        return np.array(color, dtype=np.float32)
#    
#    def _generate_atom_unique_color_id(self):
#        """ Function doc """
#        r = (self.unique_id & 0x000000FF) >> 0
#        g = (self.unique_id & 0x0000FF00) >> 8
#        b = (self.unique_id & 0x00FF0000) >> 16
#        # pickedID = r + g * 256 + b * 256*256
#        self.color_id = np.array([r/255.0, g/255.0, b/255.0], dtype=np.float32)
#    
#    def _init_vdw_rad(self):
#        """ Function doc """
#        try:
#            vdw = ATOM_TYPES[self.name][6]
#        except KeyError:
#            vdw = ATOM_TYPES[self.symbol][6]
#        return vdw
#    
#    def _init_cov_rad(self):
#        """ Function doc """
#        try:
#            cov = ATOM_TYPES[self.name][5]
#        except KeyError:
#            cov = ATOM_TYPES[self.symbol][5]
#        return cov
#    
#    def _init_ball_rad(self):
#        """ Function doc """
#        try:
#            ball = ATOM_TYPES[self.name][6]
#        except KeyError:
#            ball = ATOM_TYPES[self.symbol][6]
#        return ball
#    
#    
#    def coords(self, frame=None):
#        """ 
#        frame = int
#        
#        Returns the coordinates of an atom according to the specified frame. 
#        If no frame is specified, the frame set by easyhybrid (probably by the 
#        scale bar of the trajectory manipulation window) is used. If the object 
#        (vobject) has a smaller number of frames than the one set by the 
#        interface, the last frame of the object is used.
#        
#        return  xyz 
#        """
#        if frame is None:
#            frame  = self.vm_object.vm_session.frame
#            print (frame, len(self.vm_object.frames))
#            
#            if len(self.vm_object.frames)-1 <= frame:
#                frame = len(self.vm_object.frames)-1
#            else:
#                pass
#        
#        return self.vm_object.frames[frame, self.atom_id]
#    
#    
#    #def coords(self, frame=None):
#    #    """ Function doc """
#    #    if frame is None:
#    #        frame  = self.vm_object.vm_session.frame
#    #    # x = self.vm_object.frames[frame][(self.index-1)*3  ]
#    #    # y = self.vm_object.frames[frame][(self.index-1)*3+1]
#    #    # z = self.vm_object.frames[frame][(self.index-1)*3+2]
#    #    # return np.array([x, y, z])
#    #    return self.vm_object.frames[frame, self.atom_id]
#    
#    def get_grid_position(self, gridsize=3, frame=None):
#        """ Function doc """
#        coords = self.coords(frame)
#        gridpos = (int(coords[0]/gridsize), int(coords[1]/gridsize), int(coords[2]/gridsize))
#        return gridpos
#    
#    def get_cov_rad(self):
#        """ Function doc """
#        return self.cov_rad
#    
#    def init_color_rgb(self, name):
#        """ Return the color of an atom in RGB. Note that the returned
#            value is in scale of 0 to 1, but you can change this in the
#            index. If the atomname does not match any of the names
#            given, it returns the default dummy value of atom X.
#        """
#        try:
#            color = color =self.vm_object.color_palette[name]
#            #color = ATOM_TYPES[name][1]
#        except KeyError:
#            if name[0] == "H":# or name in self.hydrogen:
#                #color = ATOM_TYPES["H"][1]
#                color = self.vm_object.color_palette["H"]
#            
#            elif name[0] == "C":
#                #color = ATOM_TYPES["C"][1]
#                color = self.vm_object.color_palette["C"]
#            
#            elif name[0] == "O":
#                #color = ATOM_TYPES["O"][1]
#                color = self.vm_object.color_palette["O"]
#            
#            elif name[0] == "N":
#                #color = ATOM_TYPES["N"][1]
#                color = self.vm_object.color_palette["N"]
#                
#            elif name[0] == "S":
#                #color = ATOM_TYPES["S"][1]
#                color = self.vm_object.color_palette["S"]
#            else:
#                #color = ATOM_TYPES["X"][1]
#                color = self.vm_object.color_palette["X"]
#                
#        color = [int(color[0]*250), int(color[1]*250), int(color[2]*250)]
#        return color
#    
#    def init_radius(self, name):
#        """
#        """
#        try:
#            rad = ATOM_TYPES[name][6]/5.0
#        except KeyError:
#            if name[0] == "H" or name in self.hydrogen:
#                rad = ATOM_TYPES["H"][6]/5.0
#            elif name[0] == "C":
#                rad = ATOM_TYPES["C"][6]/5.0
#            elif name[0] == "O":
#                rad = ATOM_TYPES["O"][6]/5.0
#            elif name[0] == "N":
#                rad = ATOM_TYPES["N"][6]/5.0
#            elif name[0] == "S":
#                rad = ATOM_TYPES["S"][6]/5.0
#            else:
#                rad = 0.30
#        return rad
#
