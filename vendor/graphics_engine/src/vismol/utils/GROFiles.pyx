#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#

import os
import time
import multiprocessing
import numpy as np
cimport numpy as np
from vismol.core.vismol_object import VismolObject
#from vismol.model.molecular_properties import AtomTypes

from pprint import pprint

cpdef load_gro_file(infile=None, gridsize=3, vismol_session=None):
    """ Function doc 

    gridsize =

    The longest covalent bond I can find is the bismuth-iodine single bond.
    The order of bond lengths is single > double > triple.
    The largest atoms should form the longest covalent bonds. 
    So we look at atoms in the lower right corner of the Periodic Table.
    The most likely candidates are Pb, Bi, and I.
    The experimental bond lengths are:
    Bi-I = 281 pm; Pb-I = 279 pm; I-I = 266.5 pm.
    So the polar covalent Bi-I bond is the longest covalent measured so far.

    """
    #-------------------------------------------------------------------------------------------
    #                                P D B     P A R S E R 
    #-------------------------------------------------------------------------------------------
    # [EN] BUG FIX: 'at = AtomTypes()' -- AtomTypes nao existe mais em
    # vismol.model.molecular_properties (o import ja' estava comentado,
    # entao esta linha quebrava com NameError sempre que load_gro_file()
    # era chamada -- carregar um .gro estava completamente quebrado).
    # PDBFiles.pyx (mesmo parser, ja' modernizado) resolve isso de outro
    # jeito: nao calcula mais o simbolo do atomo aqui -- Atom._get_symbol()
    # ja' deriva o simbolo automaticamente a partir do nome do atomo (ver
    # vismol/model/atom.py), entao nem precisa ser passado no dict abaixo.
    with open(infile, 'r') as gro_file:
        
        grotext = gro_file.readlines()
        
        atoms, frame = get_atom_info_from_raw_line(grotext, gridsize = 3)
        frame = np.array(frame, dtype=np.float32)
        frames = [frame]
        #frames = []
    #    n    = 0 
    #    atoms = []
    #    while atoms == []:
    #        atoms     = get_list_of_atoms_from_rawframe(rawframes[n], gridsize, at =  at)
    #        n += 1
    #    
    #    
    ##-------------------------------------------------------------------------------------------
    #
    name = os.path.basename(infile)
    vismol_object  = VismolObject(name        = name, 
                                  atoms       = atoms, 
                                  vismol_session   = vismol_session, 
                                  trajectory  = frames,
                                  auto_find_bonded_and_nonbonded = True)
    '''
    #-------------------------------------------------------------------------------------------
    #                                Bonded and NB lists 
    #-------------------------------------------------------------------------------------------
    
    #atoms, bonds_full_indices, bonds_pair_of_indices, NB_indices_list = cdist.generete_full_NB_and_Bonded_lists(atoms)
    
    #-------------------------------------------------------------------------------------------
    #print (bonds_pair_of_indices, NB_indices_list )
    #for atom in atoms:
    #    pprint (atom[8])
    #-------------------------------------------------------------------------------------------
    #                         Building   V I S M O L    O B J
    #-------------------------------------------------------------------------------------------

    
    #vismol_object.non_bonded_atoms  = NB_indices_list
    #vismol_object._generate_atomtree_structure()
    #vismol_object._generate_atom_unique_color_id()
    #vismol_object.index_bonds       = bonds_full_indices
    #vismol_object.import_bonds(bonds_pair_of_indices)
	'''

    #-------------------------------------------------------------------------------------------
    return vismol_object
    
    
    
cpdef get_atom_info_from_raw_line(lines, gridsize = 3, at =  None):
    #try:
    atoms           = []
    index           = 0
    size            = int(lines[1])
    frame           = []
    for line in lines[2:size+2]:
        at_resi    = int(line[0:5])
        
        at_resn    = line[5:10].strip()

        at_name    = line[10:15].strip()
        
        #index      = int(line[15:20])
        
        x          =float(line[20:28])*10
        y          =float(line[28:36])*10
        z          =float(line[36:44])*10
        frame.append(x)
        frame.append(y)
        frame.append(z)
        at_pos     = np.array([x,y,z])
        
        at_ch      = 'X'          

        at_occup   = 0.0   #occupancy
        at_bfactor = 0.0
        at_charge  = 0.0

        gridpos  = [int(at_pos[0]/gridsize), int(at_pos[1]/gridsize), int(at_pos[2]/gridsize)]
        #ocupan   = float(line[54:60])
        #bfactor  = float(line[60:66])

                        #0      1        2        3       4        5        6       9       10          11        12      
        #atoms.append([index, at_name, at_pos, at_resi, at_resn, at_ch, [], gridpos, at_occup, at_bfactor, at_charge ])
        atoms.append({
                      'index'      : index      , 
                      'name'       : at_name    , 
                      'resi'       : at_resi    , 
                      'resn'       : at_resn    , 
                      'chain'      : at_ch      , 
                      'occupancy'  : at_occup   , 
                      'bfactor'    : at_bfactor , 
                      'charge'     : at_charge   
                      })
        
        
        #print (index, at_name, cov_rad,  at_pos, at_resi, at_resn, at_ch, at_symbol, [], gridpos, at_occup, at_bfactor, at_charge )
        index += 1

    return atoms, frame
