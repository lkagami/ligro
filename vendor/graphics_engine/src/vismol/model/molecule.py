#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
#  residue.py
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
from vismol.model.molecular_properties import residues_dictionary, solvent_dictionary



class Molecule:
    """ Class doc """
    
    def __init__(self, vismol_object, name="UNK", index=None):
        """ Class initialiser """
        self.vm_object = vismol_object
        self.atoms = {}
        #self.name = name
        #self.index = index
        #if atoms is None:
        #    self.atoms = {}
        #else:
        #    self.vm_object.atoms = atoms
        self.topology = {}
    
    def _is_protein(self):
        """ Function doc """
        # is it a protein residue?
        if self.name in residues_dictionary.keys():
            self.is_protein = True
        # is it a salvent molecule?
        if self.name in solvent_dictionary.keys():
            self.is_solvent = True
    
    def geometry_center(self, frame=0):
        """ Function doc """
        if frame > len(self.vm_object.frames)-1:
            frame = len(self.vm_object.frames)-1
        gc = np.zeros(3, dtype=np.float32)
        for atom in self.atoms.values():
            gc += atom.coords(frame)
        gc /= len(self.atoms.values())
        return gc
    
    def get_center_of_mass(self, mass=False, frame=0):
        """ Function doc """
        frame_size = len(self.vm_object.frames)-1
        
        if frame <= frame_size:
            pass
        else:
            frame = frame_size
        
        total = len(self.atoms)
        sum_x = 0.0
        sum_y = 0.0
        sum_z = 0.0
        
        for atom in self.atoms.values():
            coord = atom.coords (frame)
            sum_x += coord[0]
            sum_y += coord[1]
            sum_z += coord[2]
        
        self.mass_center = np.array([sum_x / total,
                                     sum_y / total, 
                                     sum_z / total])
    





