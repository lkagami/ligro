#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
#  chain.py
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


class Chain:
    """ Class doc """
    
    def __init__ (self, vismol_object, name=None, residues=None, label=None):
        """ Class initialiser """
        self.vm_object = vismol_object
        if residues is None:
            self.residues = {}
        else:
            self.residues = residues
        self.residues_by_index = {}
        self.backbone = []
        self.name = name
        self.backbone_pair_indexes_full = []
        self.backbone_pair_indexes_show = []
    
    def get_CA_list(self):
        """ Function doc """
        self.backbone = []
        
        for residue in self.residues:
            if residue.isProtein:
                for atom in residue.atoms:
                    if atom.name == "CA":
                        self.backbone.append(atom)
                    else:
                        pass
        
        return self.backbone
    
    def get_secundary_structure(_):
        """ Function doc """
        return None
    
    def return_name (self):
        """ Function doc """
        return self.name
