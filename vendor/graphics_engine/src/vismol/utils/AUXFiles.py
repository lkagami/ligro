#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#

from vismol.utils.debug import dprint
import os
import time
import multiprocessing
import numpy as np
from vismol.core.vismol_object import VismolObject


def load_aux_file(infile=None, vismol_session=None, gridsize=3):
    """ Function doc """
    dprint ("\nstarting: parse_aux")
    #at  =  MolecularProperties.AtomTypes()
    #initial = time.time()

    with open(infile, "r") as aux_file:
        

        frames = []
        
        aux_lines        = aux_file.read()
        text1 = aux_lines.split("HEAT_OF_FORM_UPDATED")
        
        
        # atoms names
        h1 = text1[0].split("ATOM_EL")
        h2 = h1[-1].split("ATOM_CORE")
        
        atoms_temp = []
        atoms      = []
        atoms_raw  = h2[0].split()

 
 
 
        for atom in atoms_raw:
            if "\n" in atom or "[" in atom:
                pass
            else:
                atoms_temp.append(atom)

        #first frame
        coord = text1[0].split("ATOM_X:ANGSTROMS")
        frame_coordinates = coord[-1].split("AO_ATOMINDEX")
        frame_coordinates = frame_coordinates[0].replace("\n","")
        frame_coordinates =frame_coordinates.split()
        frame  = []
        
        for coord in frame_coordinates[1:]:
        #frame_coordinates.pop[0]
            frame.append(float(coord))
        

        atom_coord = []
        
        i = 0
        index = 0
        for atom in atoms_temp:
            
            at_name = atom
            
            at_pos  = np.array([frame[i], frame[i+1], frame[i+2]])
            
            at_resi = 1
            
            at_resn = "UNK"
            
            at_ch   = "X"          
            
            at_occup   = 0.0     #occupancy
            at_bfactor = 0.0
            at_charge  = 0.0
            

            #at_symbol = line[5].split(".")
            at_symbol = at_name
            #cov_rad   = at.get_cov_rad (at_name)
            #gridpos  = [int(at_pos[0]/gridsize), int(at_pos[1]/gridsize), int(at_pos[2]/gridsize)]     
            #atoms.append([index, at_name, cov_rad,  at_pos, at_resi, at_resn, at_ch, at_symbol, [], gridpos, at_occup, at_bfactor, at_charge ])
            
            #print ([index, at_name, cov_rad,  at_pos, at_resi, at_resn, at_ch, at_symbol, [], gridpos, at_occup, at_bfactor, at_charge ])
            atoms.append({
                          "index"      : index      , 
                          "name"       : at_name    , 
                          "resi"       : at_resi    , 
                          "resn"       : at_resn    , 
                          "chain"      : at_ch      , 
                          "symbol"     : at_symbol  , 
                          "occupancy"  : at_occup   , 
                          "bfactor"    : at_bfactor , 
                          "charge"     : at_charge   
                          })
            
            
            
            index += 1
            #print (atom, frame[i], frame[i+1], frame[i+2])
            i+=3 
        
        
        
        
        frame = np.array(frame, dtype=np.float32)
        frames.append(frame)
        
        for frame_raw in text1[1:-2]:
            frame_raw = frame_raw.split("\n")
            frame = []
            
            for data in frame_raw: 
                data2 = data.split()
                if len(data2) == 3:
                    frame.append(data2[0])
                    frame.append(data2[1])
                    frame.append(data2[2])
            frame = np.array(frame, dtype=np.float32)
            frames.append(frame)
    
    name = os.path.basename(infile)
    #print (name, atoms, frames)     

    vismol_object  = VismolObject.VismolObject(name        = name, 
                                               atoms       = atoms, 
                                               vismol_session   = vismol_session, 
                                               trajectory  = frames)
    
    return vismol_object











def get_atom_list_from_xyz_frame (raw_atoms, frame = True, gridsize = 3, at = None):
    """ Function doc """
    #nCPUs =  multiprocessing.cpu_count()
    #pool  = multiprocessing.Pool(nCPUs)
    #pdb_file_lines  = frame.split("\n")   
    #atoms = (pool.map(parse_pdb_line, pdb_file_lines))
    
    atoms  = []
    frames = []
    frame_coordinates = []
    #print (raw_atoms)
    index  = 1
    for line in raw_atoms:
        line = line.split()
        if len(line) > 1:
            #print (line) 
            #index    = int(line[0])-1
            
            at_name = line[0]
            
            at_pos  = np.array([float(line[1]), float(line[2]), float(line[3])])
            
            at_resi = 1
            
            at_resn = "UNK"
            
            at_ch   = "X"          
            
            at_occup   = 0.0     #occupancy
            at_bfactor = 0.0
            at_charge  = 0.0
            

            #at_symbol = line[5].split(".")
            at_symbol = at_name
            cov_rad   = at.get_cov_rad (at_name)
            gridpos  = [int(at_pos[0]/gridsize), int(at_pos[1]/gridsize), int(at_pos[2]/gridsize)]
            
            #atoms.append([index, at_name, cov_rad,  at_pos, at_resi, at_resn, at_ch, at_symbol, [], gridpos, at_occup, at_bfactor, at_charge ])

            atoms.append([index, at_name, cov_rad,  at_pos, at_resi, at_resn, at_ch, at_symbol, [], gridpos, at_occup, at_bfactor, at_charge ])
            
            index += 1

            frame_coordinates.append(float(line[1]))
            frame_coordinates.append(float(line[2]))
            frame_coordinates.append(float(line[3]))
    
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
    #print ("Obtain bonds from original MOL2 file")
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
