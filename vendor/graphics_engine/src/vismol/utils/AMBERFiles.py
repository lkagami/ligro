#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#

from vismol.utils.debug import dprint
import os
import time
import numpy as np
import vismol.utils.c_distances as cdist
from vismol.core.vismol_object import VismolObject


def load_netcdf4_file (filein = None, visObj = None):
    """ Function doc """
    
    #from netCDF4 import Dataset
    #f = Dataset(filein)
    
    #try:
    from netCDF4 import MFDataset
    from netCDF4 import Dataset
    f = Dataset(filein)
    #except:
    #    import netCDF4 as nc
    #    f= nc.Dataset(filein)
    size = len(visObj.atoms)*3
    
    frames = np.array(f.variables["coordinates"][:], dtype=np.float32)
    frames = frames.reshape((frames.shape[0], frames.shape[1]*frames.shape[2]))


    #frames = []
    #
    ##for netcdf_frame in f.variables["coordinates"]:
    ##    frames.append(netcdf_frame)
    #
    #for netcdf_frame in f.variables["coordinates"][:]:
    #    try:
    #        frame = []
    #        for coords in netcdf_frame:
    #            frame.append(coords[0])
    #            frame.append(coords[1])
    #            frame.append(coords[2])
    #            #print (coords)
    #        frame  = np.array(frame, dtype=np.float32)
    #        frames.append(frame)
    #        frame  = []
    #    except:
    #        print('netcdf error')
    #print ('frames netcdf:',frames[0])
    return frames

    
    
def load_amber_crd_file (filein = None, visObj = None):
    """ Function doc """
    
    size = len(visObj.atoms)*3
    #size = 158*3
    filein =  open(filein, 'r')
    
    data   = filein.readlines()
    line1  = data[1].split()
    dprint(line1)
    #print(data)
    if int(float(line1[0])) == size/3:
        start = 2
    else:
        start = 1


    frames = []
    frame  = []
    frame_counter = 0
    counter       = 0
    for line in data[start:]:
        
        line2 = line.split()
        
        for coord in line2:
            #print (coord)
            frame.append(float(coord))
            counter +=1
            
            if counter == size:
                frame = np.array(frame, dtype=np.float32)
                frames.append(frame)
                
                #print(len(frames), len(frame))
                frame   = []
                counter = 0
    #for i in range(0 , len(frames[0]), 6):
    #    print (frames[0][i:i+6])
    return frames


def load_amber_topology_file (infile = None, vismolSession =  None, gridsize = 3):
    """ Function doc """
    #at  =  MolecularProperties.AtomTypes()
    filename = infile
    infile = open(infile, 'r')
    text = infile.read()
    #print(text)
    text2 = text.split("%FLAG")
    
    atom_names      = None
    atomic_numbers  = None
    bonds_H         = None
    bonds_noH       = None
    residues_names  = None  
    pointers        = None
    
    at_names    = []
    at_numbers  = []
    bonds       = []
    bonds2      = []
    res_names   = []
    res_indexes = []
    
    total_bonds = []
    for block in text2:
        
        block2 = block.split('\n')
        string = ''
        for block3  in block2[2:]:
            string += block3
        
        #print(string)
        if 'ATOM_NAME' in block:
            atom_names = string
            #print(atom_names)
            for i in range(0, len(atom_names),4):
                at_name = atom_names[i:i+4]
                at_names.append(at_name.strip())
        
        if 'ATOMIC_NUMBER' in block:
            atomic_numbers = string
            #atomic_numbers = atomic_numbers.split()
            
            for i in range(0, len(atomic_numbers),8):
                atomic_number = atomic_numbers[i:i+8]
                at_numbers.append(int(atomic_number.strip()))
        
        if 'BONDS_INC_HYDROGEN' in block:
            bonds_H = string
            
            for i in range(0, len(bonds_H),8):
                bond = bonds_H[i:i+8]
                #bond = int((int(bond.strip())/3) + 1) 
                bond = int((int(bond.strip())/3)) 
                bonds.append(bond)
            for i in range(0, len(bonds),3):
                #print(bonds[i:i+2])
                total_bonds.append(bonds[i:i+2])
                
        if 'BONDS_WITHOUT_HYDROGEN' in block:
            bonds_noH = string
            for i in range(0, len(bonds_noH),8):
                bond = bonds_noH[i:i+8]
                #bond = int((int(bond.strip())/3) + 1) 
                bond = int((int(bond.strip())/3)) 
                bonds2.append(bond)
            
            for i in range(0, len(bonds2),3):
                #print(bonds2[i:i+2])
                total_bonds.append(bonds2[i:i+2])




        if 'RESIDUE_LABEL' in block:
            labels = string
            res_names_short = []
            #print (labels)
            for i in range(0, len(labels),4):
                label = labels[i:i+4]
                res_names_short.append(label.strip())
                #bonds2.append(bond)
            
            
        if 'RESIDUE_POINTER' in block:
            
            pointers = string.split()
            #print (pointers)
            pointers2 = []
            
            for i in pointers:
                pointers2.append(int(i))
            
            pointers2.append(len(at_names)+1)
            
            #print (pointers2)
            
            
            
            idx = 0
            residue_counter = 1
            
            for atom in at_names:
                if idx < pointers2[residue_counter]-1:
                    #print (idx, atom, residue_counter, pointers2, res_names_short[residue_counter-1])
                    
                    res_indexes.append(residue_counter)
                    res_names.append(res_names_short[residue_counter-1])
                else:
                    residue_counter += 1 
                    #print (idx, atom, residue_counter, pointers2, res_names_short[residue_counter-1])
                    res_indexes.append(residue_counter)
                    res_names.append(res_names_short[residue_counter-1])
                idx+=1
            
            
            
            '''
            residue_counter = 1
            for i in range(0, len(at_names)):
                
                if i <= pointers2[residue_counter]-1:
                    dprint (pointers2[residue_counter]-1, residue_counter, pointers2, at_names[i],i )
                    pass

                    #print (i+1, pointers, residue_counter, at_names[i], res_names_short[residue_counter-1])
                    res_indexes.append(residue_counter)
                    res_names.append(res_names_short[residue_counter-1])
                else:
                    residue_counter += 1 
                #print (len(res_indexes), len(res_names))
            '''
            #for i in range(0, len(pointers),8):
            #    pointers = pointers[i:i+8]
            #    print (pointers ) 
            #    res_indexes.append(pointers)
            #    #bonds2.append(bond)
 



    #print (len(at_names), len(at_numbers), len(res_names),len(res_indexes))
    #print (at_numbers)
    #print (at_names, at_numbers)

    #print (bonds)
    #print (bonds2)
    #print (res_names_short)
    

    #print (res_names)
    #print (res_indexes)
    
    resi_counter  = 1 
    atoms = []
    
    for index in range(len(at_names)):
        
        at_name  = at_names[index]
        at_pos   = np.array([0.0,0.0,0.0])
        
        
        
        at_resn  = res_names[index]
        at_resi  = res_indexes[index]
        at_ch    = 'X'
        
        at_occup   = 0.0     #occupancy
        at_bfactor = 0.0
        at_charge  = 0.0
                
        #print (index, at_name, at_resn, at_resi,  at_chain, at_pos, at_numbers[index] )
        #print(at.ATOM_TYPES_BY_ATOMICNUMBER, at_numbers[index])
        #if at_numbers != []:
        #    at_symbol = at.ATOM_TYPES_BY_ATOMICNUMBER[at_numbers[index]]
        #else:
        at_symbol = None #at.get_symbol(at_name)
        #cov_rad   = at.get_cov_rad (at_symbol)
        #gridpos   = [0,0,0]
        

        #atoms.append([index, at_name, cov_rad,  at_pos, at_resi, at_resn, at_chain, at_symbol, [], gridpos, at_occup, at_bfactor, at_charge ])
        #print([index, at_name, cov_rad,  at_pos, at_resi, at_resn, at_chain, at_symbol, [], gridpos, at_occup, at_bfactor, at_charge ])
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
            
    #print (total_bonds)

    
    bonds_pair_of_indexes = []
    for bond in total_bonds:
        bonds_pair_of_indexes.append(bond[0])
        bonds_pair_of_indexes.append(bond[1])
        #print(bond, atoms[bond[0]], atoms[bond[1]])


    name = os.path.basename(filename)
    vismol_object  = VismolObject.VismolObject(name                           = name       , 
                                               atoms                          = atoms      , 
                                               vismolSession                  = vismolSession  , 
                                               bonds_pair_of_indexes          = bonds_pair_of_indexes,
                                               trajectory                     = []         ,
                                               auto_find_bonded_and_nonbonded = False      )
        
    return   vismol_object  
    
        

        
        
