#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
cimport numpy as np


#cython: boundscheck=False
#cython: wraparound=False
#cython: initializedcheck=False
#cython: cdivision=True

from collections import deque
cimport cython

cpdef list find_subgroup(int atom1, int atom2, dict top):
    """
    Iterative, efficient, and cycle-safe version.
    Returns all atoms accessible from atom2, excluding atom1 and avoiding cyclic loops.

    Parameters:

        atom1 : int
            Initial “parent” atom — used to prevent traversal back.

        atom2 : int
            Starting atom for the search (the traversal origin).

        top : dict[int, list[int]]
            System connectivity (graph). Example:
            top[i] = [list of atoms connected to atom i]

    Return:
        list[int] — list of atoms belonging to the subgroup.
    """

    cdef list subgroup = []
    cdef object visited = set([atom1])        # Python set
    cdef object queue = deque([atom2])        # Python deque

    cdef int current, neighbor
    cdef list neighbors

    while queue:
        current = queue.pop()  # ou popleft() se quiser BFS

        if current not in top:
            continue  # evita quebrar e mantém consistência

        neighbors = top[current]

        for neighbor in neighbors:
            if neighbor in visited:
                continue

            visited.add(neighbor)
            subgroup.append(neighbor)

            if neighbor in top and len(top[neighbor]) > 1:
                queue.append(neighbor)

    return subgroup





# =====================================================================
# GRID OFFSET
# =====================================================================

cpdef list calculate_grid_offset(double gridsize, double maxbond=2.6):
    """
    Compute grid offsets within maxbond distance.

 
    This function calculates the offset for each grid element in the atomic grid. 
    The grid elements define regions in 3D space within a certain distance from 
    each atom. The gridsize parameter controls the size of each grid element, 
    and the maxbond parameter determines the maximum distance within which bonds 
    are considered. The function returns a list of grid offsets representing 
    grid elements.
 
    
    """
    cdef int border = int(maxbond / gridsize)
    cdef list offsets = []
    cdef int i, j, k

    #      -------- First floor (k=0) --------
    #   |-------|-------|-------|-------|-------| 
    #   |\\\\\\\|\\\\\\\|\\\\\\\|\\\\\\\|\\\\\\\| 
    #   |\\\\\\\|\\\\\\\|\\\\\\\|\\\\\\\|\\\\\\\| 
    #   |-2,2,0 |-1,2,0 | 0,2,0 | 1,2,0 | 2,2,0 | 
    #   |-------|-------|-------|-------|-------| 
    #   |\\\\\\\|\\\\\\\|\\\\\\\|\\\\\\\|\\\\\\\| 
    #   |\\\\\\\|\\\\\\\|\\\\\\\|\\\\\\\|\\\\\\\| 
    #   |-1,1,0 |-1,1,0 | 0,1,0 | 1,1,0 | 2,1,0 | 
    #   |-------|-------|-------|-------|-------| 
    #   |       |       |XXXXXXX|\\\\\\\|\\\\\\\|
    #   |       |       |XXXXXXX|\\\\\\\|\\\\\\\|
    #   |-1,0,0 |-1,0,0 | 0,0,0 | 1,0,0 | 2,0,0 |
    #   |-------|-------|-------|-------|-------|
    #   |       |       |       |       |       |
    #   |       |       |       |       |       |
    #   |-1,-1,0|-1,-1,0| 0,-1,0| 1,-1,0| 2,-1,0|
    #   |-------|-------|-------|-------|-------|
    
    for i in range(-border, border + 1):
        for j in range(0, border + 1):

            # Skip unused regions
            if i < -1 and j == 0:
                continue

            # Skip origin
            if i == 0 and j == 0:
                continue

            offsets.append([i, j, 0])

    #                    -------- Floors k > 0 --------
    #--------------------- floors above---------------------------------
    #                                     |-------|-------|-------|-------| 
    #                                     |\\\\\\\|\\\\\\\|\\\\\\\|\\\\\\\| 
    #                                     |\\\\\\\|\\\\\\\|\\\\\\\|\\\\\\\| 
    #                                     |-1,2,2 | 0,2,2 | 1,2,2 | 2,2,2 | 
    #                                     |-------|-------|-------|-------| 
    # |-------|-------|-------|-------|   |\\\\\\\|\\\\\\\|\\\\\\\|\\\\\\\| 
    # |\\\\\\\|\\\\\\\|\\\\\\\|\\\\\\\|   |\\\\\\\|\\\\\\\|\\\\\\\|\\\\\\\| 
    # |\\\\\\\|\\\\\\\|\\\\\\\|\\\\\\\|   |-1,1,2 | 0,1,2 | 1,1,2 | 2,1,2 | 
    # |-1,2,1 | 0,2,1 | 1,2,1 | 2,2,1 |   |-------|-------|-------|-------| 
    # |-------|-------|-------|-------|   |\\\\\\\|\\\\\\\|\\\\\\\|\\\\\\\|
    # |\\\\\\\|\\\\\\\|\\\\\\\|\\\\\\\|   |\\\\\\\|\\\\\\\|\\\\\\\|\\\\\\\|
    # |\\\\\\\|\\\\\\\|\\\\\\\|\\\\\\\|   |-1,0,2 | 0,0,2 | 1,0,2 | 2,0,2 |
    # |-1,1,1 | 0,1,1 | 1,1,1 | 2,1,1 |   |-------|-------|-------|-------|
    # |-------|-------|-------|-------|   |\\\\\\\|\\\\\\\|\\\\\\\|\\\\\\\|
    # |\\\\\\\|\\\\\\\|\\\\\\\|\\\\\\\|   |\\\\\\\|\\\\\\\|\\\\\\\|\\\\\\\|
    # |\\\\\\\|\\\\\\\|\\\\\\\|\\\\\\\|   |-1,-1,2| 0,-1,2| 1,-1,2| 2,-1,2|
    # |-1,0,1 | 0,0,1 | 1,0,1 | 2,0,1 |   |-------|-------|-------|-------|  
    # |-------|-------|-------|-------|
    # |\\\\\\\|\\\\\\\|\\\\\\\|\\\\\\\|
    # |\\\\\\\|\\\\\\\|\\\\\\\|\\\\\\\|
    # |-1,-1,1| 0,-1,1| 1,-1,1| 2,-1,1|
    # |-------|-------|-------|-------|

    for i in range(-border, border + 1):
        for j in range(-border, border + 1):
            for k in range(1, border + 1):
                offsets.append([i, j, k])

    return offsets


# =====================================================================
# DISTANCE
# =====================================================================

cpdef double calculate_sqrt_distance(int i, int j, coords):
    """
    Compute squared distance.
    
    This Cython function calculates the squared distance between two atoms 
    with indices i and j using their coordinates coords. It is used to 
    calculate the distance between atoms and determine if a bond exists 
    between them. 
     
    """
    cdef double dx = coords[i][0] - coords[j][0]
    cdef double dy = coords[i][1] - coords[j][1]
    cdef double dz = coords[i][2] - coords[j][2]
    return dx*dx + dy*dy + dz*dz


# =====================================================================
# CONNECTIONS WITHIN ONE GRID ELEMENT
# =====================================================================

cpdef list get_connections_within_grid_element(
        list atoms,
        coords,
        cov_rad,
        double tolerance,
        double gridsize):

    """
    This function calculates the bonds between atoms within a single 
    element of the atomic grid. It takes a list of atoms, their coordinates, 
    covalent radii (cov_rad), a tolerance factor, and the grid size. 
    It returns a list of pairs of atom indices representing the bonds 
    within the grid element.
    
        Calculate the distances and bonds 
        between atoms within a single element 
        of the atomic grid
        
                  |-------|-------|-------|
                  |       |       |       |
                  |       |       |       |
                  |       |       |       |
                  |-------|-atoms-|-------|
                  |       |       |       |
                  |       | i<->j |       |
                  |       |       |       |
                  |-------|-------|-------|
                  |       |       |       |
                  |       |       |       |
                  |       |       |       |
                  |-------|-------|-------|
    
    
    
    atoms = [[index, at_name, cov_rad,  at_pos, at_res_i, at_res_n, at_ch], ...]
            each elemte is a list contain required data.
    
    
    bonds_pair_of_indexes [[a,b],[b,c], ...] where a and b are indices. 
    returns a list of pair of indices "bonds_pair_of_indexes"
    
    """

    cdef list bonds = []
    cdef int idx_i, idx_j, i, j
    cdef double r2, cutoff

    for i in range(len(atoms) - 1):
        idx_i = atoms[i]
        for j in range(i + 1, len(atoms)):
            idx_j = atoms[j]

            r2 = calculate_sqrt_distance(idx_i, idx_j, coords)
            cutoff = (cov_rad[idx_i] + cov_rad[idx_j]) ** 2 * tolerance

            if r2 <= cutoff:
                bonds.append(idx_i)
                bonds.append(idx_j)

    return bonds


# =====================================================================
# CONNECTIONS BETWEEN GRID ELEMENTS
# =====================================================================

cpdef list get_connections_between_grid_elements(
        list grid1,
        list grid2,
        coords,
        cov_rad,
        double tolerance,
        double gridsize):

    '''
    This function calculates the bonds between atoms in two different grid
    elements. It takes two atomic grids, their coordinates, covalent radii, 
    a tolerance factor, and the grid size. It returns a list of pairs of atom
    indices representing the bonds between the two grid elements.
    '''

    cdef list bonds = []
    cdef int i, j
    cdef double r2, cutoff

    if grid1 is grid2:
        return bonds

    for i in grid1:
        for j in grid2:
            if i == j:
                continue

            r2 = calculate_sqrt_distance(i, j, coords)
            cutoff = (cov_rad[i] + cov_rad[j]) ** 2 * tolerance

            if r2 <= cutoff:
                bonds.append(i)
                bonds.append(j)

    return bonds


# =====================================================================
# BUILD ATOMIC GRID
# =====================================================================

cpdef dict build_the_atomic_grid(list indexes, list positions):
    """
    Build dict: grid_pos → list of atom indices
 
    This function builds the atomic grid by grouping atoms based on their 
    grid positions. It takes a list of atom indices and their corresponding 
    grid positions. It returns a dictionary where each key represents a grid 
    element, and the value is a list of atom indices within that grid element.
     
    """
    cdef dict grid = {}
    cdef int a

    for a, pos in enumerate(positions):
        if pos not in grid:
            grid[pos] = []
        grid[pos].append(indexes[a])

    return grid


# =====================================================================
# MAIN BOND FINDER
# =====================================================================

cpdef list get_atomic_bonds_from_grid(
        list indexes,
        coords,
        cov_rad,
        list gridpos_list,
        double gridsize,
        double maxbond,
        double tolerance):
    '''
    function is responsible for calculating the bonds between atoms within a given grid. It takes the following parameters:

    indexes: A list of atom indices.
    coords:  A list of atom coordinates.
    cov_rad: A list of covalent radii for each atom.
    
    gridpos_list: A list of grid positions that represent the atomic grid.
    gridsize:     A double representing the grid size.
    maxbond:      A double representing the maximum bond length.
    
    
    The function first builds an atomic grid using the build_the_atomic_grid 
    function, which organizes atoms into grid elements based on their positions. 
    It then calculates the grid offset using the calculate_grid_offset function. 
    The grid offset is used to define neighboring grid elements around each 
    element in the atomic grid.

    The function iterates over each grid element and calculates the bonds 
    between atoms within that element. If there is more than one atom in 
    the element, it calls the get_connections_within_grid_element function 
    to calculate the intra-element bonds. The result is added to the 
    bonds_pair_of_indexes list.

    Next, the function iterates over neighboring grid elements based on the 
    grid offset and calculates the bonds between atoms in the current element 
    and those in the neighboring element. The result is again added to the 
    bonds_pair_of_indexes list.

    Finally, the function returns the bonds_pair_of_indexes list containing 
    the pairs of atom indices representing the bonds within the atomic grid.
    '''
    cdef dict atomic_grid = build_the_atomic_grid(indexes, gridpos_list)
    cdef list offsets = calculate_grid_offset(gridsize, maxbond)
    cdef list bonds = []
    cdef tuple element, nei
    cdef list atoms1, atoms2

    for element in atomic_grid:

        atoms1 = atomic_grid[element]

        # Intra–grid bonds
        if len(atoms1) > 1:
            bonds += get_connections_within_grid_element(
                atoms1, coords, cov_rad, tolerance, gridsize
            )

        # Inter–grid bonds
        for off in offsets:
            nei = (element[0] + off[0], element[1] + off[1], element[2] + off[2])

            if nei in atomic_grid:
                atoms2 = atomic_grid[nei]
                bonds += get_connections_between_grid_elements(
                    atoms1, atoms2, coords, cov_rad, tolerance, gridsize
                )

    return bonds
