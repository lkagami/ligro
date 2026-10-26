#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#

import multiprocessing
import numpy as np
cimport numpy as np


'''
cpdef tuple selection_spherical_expansion (set selected_indexes, set selectable_indexes,  coordinates, float radius ):
        cdef float radius_sqr
        cdef float dx
        cdef float dy
        cdef float dz
        
        radius_sqr = radius**2
        new_selected_indexes   = set()
        
        for i in selected_indexes:
            #print ('len(selectable_indexes)', len(selectable_indexes))
            for j in selectable_indexes:
                
                i_xyz  = coordinates[i]
                j_xyz  = coordinates[j]
                
                dx = (i_xyz[0] - j_xyz[0])**2
                dy = (i_xyz[1] - j_xyz[1])**2
                dz = (i_xyz[2] - j_xyz[2])**2
                
                if dx + dy + dz <= radius_sqr:
                    new_selected_indexes.add(j)
                else:
                    pass
            
            selectable_indexes = selectable_indexes - new_selected_indexes
        
        for index in selected_indexes:
            new_selected_indexes.add(index)

        return new_selected_indexes,  selectable_indexes
'''

# cython: boundscheck=False
# cython: wraparound=False
# cython: cdivision=True

import numpy as np
cimport numpy as cnp
cimport cython

ctypedef cnp.float64_t DTYPE_t


@cython.boundscheck(False)
@cython.wraparound(False)
cpdef tuple selection_spherical_expansion(
        set selected_indexes,
        set selectable_indexes,
        coordinates,
        float radius):

    cdef:
        double radius_sqr
        double dx, dy, dz
        double dist2
        double ix, iy, iz
        Py_ssize_t i, j
        set new_selected_indexes

        cnp.ndarray[DTYPE_t, ndim=2] coords

    # Aceita qualquer ndarray compatível sem mudar API
    coords = np.asarray(coordinates, dtype=np.float64)

    radius_sqr = radius * radius
    new_selected_indexes = set()

    for i in selected_indexes:

        # lê coordenadas uma vez por átomo selecionado
        ix = coords[i, 0]
        iy = coords[i, 1]
        iz = coords[i, 2]

        # iterar sobre snapshot evita problemas ao alterar o set depois
        for j in selectable_indexes:

            dx = ix - coords[j, 0]
            dy = iy - coords[j, 1]
            dz = iz - coords[j, 2]

            dist2 = dx * dx + dy * dy + dz * dz

            if dist2 <= radius_sqr:
                new_selected_indexes.add(j)

        # evita criar um novo set
        selectable_indexes.difference_update(new_selected_indexes)

        # early exit: nada mais para selecionar
        if not selectable_indexes:
            break

    # preserva seleção original
    new_selected_indexes.update(selected_indexes)

    return new_selected_indexes, selectable_indexes









'''


cpdef tuple selection_spherical_expansion(set selected_indexes, set selectable_indexes, coordinates, float radius):
    # Expands a selection of points in 3D space by adding all candidate points
    # that fall within a spherical radius from the already selected points.
    # Typical use: computational chemistry / molecular visualization
    # (e.g., "select all atoms within a 5 Å radius around this residue").
    #
    # Parameters:
    #   selected_indexes   - set of already selected indices (expansion seeds)
    #   selectable_indexes - set of candidate indices that may be included
    #   coordinates        - xyz coordinates for each index
    #   radius             - selection sphere radius
    #
    # Returns:
    #   tuple (new_expanded_selection, remaining_candidates)

    # Cython typing (cdef float) -> compiled code for better performance
    # in this double loop, which can be computationally expensive.
    cdef float radius_sqr
    cdef float dx
    cdef float dy
    cdef float dz

    # Work with the squared radius to avoid the square root in distance
    # calculations, since sqrt is computationally expensive.
    radius_sqr = radius**2
    new_selected_indexes = set()

    # For each already selected point...
    for i in selected_indexes:
        #print ('len(selectable_indexes)', len(selectable_indexes))
        i_xyz = coordinates[i]
        # ...iterate through all candidate points.
        for j in selectable_indexes:
            j_xyz = coordinates[j]

            # Squared distance between the two points (dx^2 + dy^2 + dz^2).
            #dx = (i_xyz[0] - j_xyz[0])**2
            #dy = (i_xyz[1] - j_xyz[1])**2
            #dz = (i_xyz[2] - j_xyz[2])**2
            
            dx = i_xyz[0] - j_xyz[0]
            dy = i_xyz[1] - j_xyz[1]
            dz = i_xyz[2] - j_xyz[2]

            dist2 = dx*dx + dy*dy + dz*dz
            
            
            # If the candidate lies inside the sphere, add it to the new selection.
            if dist2 <= radius_sqr:
            #if dx + dy + dz <= radius_sqr:
                new_selected_indexes.add(j)
            else:
                pass

        # Remove newly selected points from the candidate set so they are not
        # processed again. Optimization: starting from the second selected point,
        # the candidate set is already reduced. Since a 'set' has no guaranteed
        # iteration order, small differences may occur depending on iteration order.
        #selectable_indexes = selectable_indexes - new_selected_indexes
        selectable_indexes.difference_update(new_selected_indexes)
    
    # Ensure the original selection is preserved in the final result.
    for index in selected_indexes:
        new_selected_indexes.add(index)

    # Return the expanded selection and the remaining candidates.
    return new_selected_indexes, selectable_indexes
'''
