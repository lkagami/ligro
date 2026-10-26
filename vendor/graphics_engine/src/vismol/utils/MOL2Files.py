#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#

from vismol.utils.debug import dprint
import os
import time
import multiprocessing
import numpy as np
import vismol.utils.c_distances as cdist
from vismol.core.vismol_object import VismolObject


def load_mol2_files(infile=None, vismol_session=None, gridsize=3):
    """ Function doc """
    dprint ('\nstarting: parse_mol2')

    #at  = MolecularProperties.AtomTypes()
    with open(infile, 'r') as mol2_file:
        filetext = mol2_file.read()

        molecules     =  filetext.split('@<TRIPOS>MOLECULE')
        firstmolecule =  molecules[1].split('@<TRIPOS>ATOM')
        header        =  firstmolecule[0]
        firstmolecule =  firstmolecule[1].split('@<TRIPOS>BOND')
        raw_atoms     =  firstmolecule[0]
        bonds         =  firstmolecule[1]

    header    = header.split('\n')
    raw_atoms = raw_atoms.split('\n')
    bonds     = bonds.split('\n')

    dprint (raw_atoms)
    atoms, frames = get_atom_list_from_mol2_frame(raw_atoms = raw_atoms, frame = True)#,  gridsize = gridsize,  at = at)

    #-------------------------------------------------------------------------------------------
    #                         Building   V I S M O L    O B J
    #-------------------------------------------------------------------------------------------
    name = os.path.basename(infile)
    vismol_object  = VismolObject.VismolObject(name        = name, 
                                               atoms       = atoms, 
                                               vismol_session   = vismol_session, 
                                               trajectory  = frames)
    #-------------------------------------------------------------------------------------------
    return vismol_object

def get_atom_list_from_mol2_frame(raw_atoms, frame=True):#, gridsize = 3, at =  None):
    """ Function doc """
    #nCPUs =  multiprocessing.cpu_count()
    #pool  = multiprocessing.Pool(nCPUs)
    #pdb_file_lines  = frame.split('\n')   
    #atoms = (pool.map(parse_pdb_line, pdb_file_lines))
    
    atoms  = []
    frames = []
    frame_coordinates = []
    #print (raw_atoms)
    for line in raw_atoms:
        line = line.split()
        if len(line) > 1:
            #print (line) 
            index    = int(line[0])-1
            
            at_name  = line[1]
            
            at_pos   = np.array([float(line[2]), float(line[3]), float(line[4])])
            
            at_resi = int(line[6])
            
            at_resn = line[7]


            at_ch   = 'X'          

            at_occup   = 0.0     #occupancy
            at_bfactor = 0.0
            at_charge  = float(line[8])


            at_symbol  = None# at.get_symbol(at_name)
            
            #at_symbol= line[5].split('.')
            #at_symbol= at_symbol[0]
            #cov_rad  = at.get_cov_rad (at_symbol)



            #gridpos  = [int(at_pos[0]/gridsize), int(at_pos[1]/gridsize), int(at_pos[2]/gridsize)]
            #atoms.append([index, at_name, cov_rad,  at_pos, at_resi, at_resn, at_ch, at_symbol, [], gridpos, at_occup, at_bfactor, at_charge ])
            #atoms.append([index, at_name, cov_rad,  at_pos, at_resi, at_resn, at_ch, at_symbol, [], gridpos, at_occup, at_bfactor, at_charge ])
            
            atoms.append({
                          'index'      : index      , 
                          'name'       : at_name    , 
                          'resi'       : at_resi    , 
                          'resn'       : at_resn    , 
                          'chain'      : at_ch      , 
                          'symbol'     : at_symbol  , 
                          'occupancy'  : at_occup   , 
                          'bfactor'    : at_bfactor , 
                          'charge'     : at_charge   
                          })


            #atoms.append([index, at_name, cov_rad,  at_pos, at_res_i, at_res_n, at_ch])
            
            #atom     = Atom(name      =  at_name, 
            #                index     =  index, 
            #                pos       =  at_pos, 
            #                resi      =  at_res_i, 
            #                resn      =  at_res_n, 
            #                chain     =  at_ch, 
            #                )
            #atoms.append(atom)
            frame_coordinates.append(float(line[2]))
            frame_coordinates.append(float(line[3]))
            frame_coordinates.append(float(line[4]))
    frame_coordinates = np.array(frame_coordinates, dtype=np.float32)
    frames.append(frame_coordinates)
    #print (frames)
    #print (atoms)
    return atoms, frames#, coords

def get_bonds (raw_bonds):
    """ Function doc """
    index_bonds              = []
    index_bonds_pairs        = []
    index_bonds_pairs_orders = []
    
    #print (raw_bonds)
    #print ('Obtain bonds from original MOL2 file')
    for line in raw_atoms:
        line = line.split()
        if len(line) == 4:
            index    = int(line[0])            
            atom1    = int(line[1]-1)
            atom2    = int(line[2]-1)
            order    = line[3]

            index_bonds      .append(atom1)
            index_bonds      .append(atom2)
            index_bonds_pairs.append([atom1,atom2])
            
            index_bonds_pairs_orders.append(order)

    return [index_bonds, index_bonds_pairs]
