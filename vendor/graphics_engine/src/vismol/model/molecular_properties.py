#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
import yaml
import os


'''
CARBON_COLOR_PICK = {
                     0 : [0.      , 1.0     , 0.      ],
                     1 : [0.3     , 0.0     , 0.46    ],
                     2 : [0.8     , 0.5     , 0.1     ],
                     3 : [0.0     , 0.5     , 0.5     ],
                     4 : [0.6     , 0.4     , 0.4     ]
                     }
'''

solvent_dictionary  = {
                       'WAT': 'O', 
                       'SOL': 'O', 
                       'HOH': 'O', 
                       'H2O': 'O', 
                       }


residues_dictionary = {'CYS': 'C', 
                       'ASP': 'D', 
                       'SER': 'S', 
                       'GLN': 'Q', 
                       'LYS': 'K',
                       'ILE': 'I', 
                       'PRO': 'P', 
                       'THR': 'T', 
                       'PHE': 'F', 
                       'ASN': 'N', 
                       'GLY': 'G', 
                       'HIS': 'H', 
                       
                       # amber
                       "HID": "H",
                       "HIE": "H",
                       "HIP": "H",
                       "ASH": "D",
                       "GLH": "E",
                       "CYX": "C",
                       
                       # charmm
                       "HSD": "H", 
                       "HSE": "H", 
                       "HSP": "H", 
                       
                       
                       'LEU': 'L', 
                       'ARG': 'R', 
                       'TRP': 'W', 
                       'ALA': 'A', 
                       'VAL': 'V', 
                       'GLU': 'E', 
                       'TYR': 'Y', 
                       'MET': 'M'}

