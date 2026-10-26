#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#

from vismol.utils.debug import dprint
import os
import numpy as np
from vismol.core.vismol_object import VismolObject


def load_PSF_topology_file(infile=None, vismol_session=None, gridsize=3):
    """ Function doc """
    filename = infile
    infile = open(infile, "r")
    text = infile.read()
    
    text2 = text.split("!")
    
    atoms = []
    bonds = []
    
    for block in text2: 
        text3 = block.split("\n")
        if "NATOM" in  text3[0]:
         
            for line in text3:
                line2 = line.split()
                if len(line2)> 6:
                    #print (line2)
        
        
        
                    index      = int(line2[0]) -1
                    at_resi    = int(line2[2])
                    at_resn    = line2[3]
                    at_name    = line2[4]
                    #index      = int(line[15:20])
                    

                    at_pos     = np.array([0.0,0.0,0.0])
                    try:
                        at_ch      = line2[1][0]
                    except:
                        at_ch      = "X"          
                    
                    #at_symbol  = "H"
                    
                    at_symbol  = None#at.get_symbol(at_name)

                    at_occup   = 0.0   #occupancy
                    at_bfactor = 0.0
                    at_charge  = float(line2[6])
                    gridpos    = [0,0,0]
                    #cov_rad  = 0.00 #at.get_cov_rad (at_symbol)
                    #gridpos  = 0.00 #[int(at_pos[0]/gridsize), int(at_pos[1]/gridsize), int(at_pos[2]/gridsize)]
                    
                    #cov_rad  = at.get_cov_rad (at_symbol)
                    #gridpos  = [int(at_pos[0]/gridsize), int(at_pos[1]/gridsize), int(at_pos[2]/gridsize)]
                    
                    #ocupan   = float(line[54:60])
                    #bfactor  = float(line[60:66])
                                    #0      1        2        3       4        5        6       7       8       9       10          11        12      
                    #atoms.append([index, at_name, cov_rad,  at_pos, at_resi, at_resn, at_ch, at_symbol, [], gridpos, at_occup, at_bfactor, at_charge ])
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
                    
                    #print([index, at_name, cov_rad,  at_pos, at_resi, at_resn, at_ch, at_symbol, [], gridpos, at_occup, at_bfactor, at_charge ])
                    
        
        
        if "NBOND" in  text3[0]:
            for line in text3:
                line2 = line.split()
                if len(line2)> 6:
                    
                    for i in range(0, len(line2),2):
                        dprint(int(line2[i]), int(line2[i+1]))
                        #print(line2[i:i+2])
                        bonds.append(int(line2[i])-1)
                        bonds.append(int(line2[i+1])-1)
        
        else:
            pass
        
    #print (bonds)
        
        
        
    name = os.path.basename(filename)
    vismol_object  = VismolObject.VismolObject(name                           = name       ,    
                                               atoms                          = atoms      ,    
                                               vismol_session                      = vismol_session  ,    
                                               bonds_pair_of_indexes          = bonds      ,    
                                               trajectory                     = []         ,    
                                               auto_find_bonded_and_nonbonded = False      )    
            
    return   vismol_object
