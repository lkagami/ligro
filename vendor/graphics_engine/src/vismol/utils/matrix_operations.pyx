#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
#  matrix_operations.pyx
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
cimport numpy as np

cpdef my_glSideVectorAbs(in_matrix):
    """ Uses a numpy array as input to obtain the side vector untransformed,
        i.e. when you apply a modification to a vector using a matrix,
        its coordinates are changed, so this function returns the
        "original" vector or the vector in tha ABSOLUTE coordinates
        system. The entry numpy array has shape (4, 4), size of 16 and ndim=2.
        
        Keyword arguments:
        in_matrix -- a numpy array that can be the PROJECTION, VIEW or 
                     MODEL matrix.
        
        Returns:
            side -- Numpy array of 3 elements corresponding to the Side vector.
    """
    inv_mat = np.matrix(in_matrix, dtype=np.float32).I
    side = np.array([[1],[0],[0],[1]],dtype=np.float32)
    side = np.array(inv_mat*side).T
    return side[0,:3]

cpdef my_glUpVectorAbs(in_matrix):
    """ Uses a 4x4 matrix as input to obtain the up vector untransformed,
        i.e. when you apply a modification to a vector using a matrix,
        its coordinates are changed, so this function returns the
        "original" vector or the vector in tha ABSOLUTE coordinates
        system.
        
        Keyword arguments:
        in_matrix -- a 4x4 matrix that can be the PROJECTION, VIEW or 
                     MODEL matrix.
        
        Returns:
            Numpy array of 3 elements corresponding to the Up vector.
    """
    inv_mat = np.matrix(in_matrix, dtype=np.float32).I
    side = np.array([[0],[1],[0],[1]],dtype=np.float32)
    side = np.array(inv_mat*side).T
    return side[0,:3]

cpdef my_glForwardVectorAbs(in_matrix):
    """ Uses a 4x4 matrix as input to obtain the forward vector untransformed,
        i.e. when you apply a modification to a vector using a matrix,
        its coordinates are changed, so this function returns the
        "original" vector or the vector in tha ABSOLUTE coordinates
        system.
        
        Keyword arguments:
        in_matrix -- a 4x4 matrix that can be the PROJECTION, VIEW or 
                     MODEL matrix.
        
        Returns:
            Numpy array of 3 elements corresponding to the Forward vector.
    """
    inv_mat = np.matrix(in_matrix, dtype=np.float32).I
    side = np.array([[0],[0],[-1],[1]],dtype=np.float32)
    side = np.array(inv_mat*side).T
    return side[0,:3]

cpdef my_glTranslatef(np.ndarray in_matrix, np.ndarray position):
    """ Creates a translation matrix using an identity matrix and the
        position's coordinates at the last row, but only the xyz components.
        As initial result, you get the following matrix:
        
                        |1 0 0 0|       |1 0 0 0|
                        |0 1 0 0|  ==>  |0 1 0 0|
                        |0 0 1 0|  ==>  |0 0 1 0|
                        |0 0 0 1|       |x y z 1|
        
        then, this translation matrix is multiplied to the input matrix to
        obtain the final result, a translated matrix to position.
        
        Keyword arguments:
            in_matrix -- a 4x4 matrix that you want to apply a translation.
            position -- a 3 elements numpy ndarray with the coordinates of
                        your desired translation.
        
        Returns:
            A multiplied 4x4 matrix of [in_matrix] x [trans_matrix].
    """
    cdef np.ndarray[np.float32_t, ndim=2] trans_matrix = np.identity(4, dtype=np.float32)
    trans_matrix[3,:3] = position
    return my_glMultiplyMatricesf(in_matrix, trans_matrix)

cpdef my_glMultiplyMatricesf(np.ndarray mat1, np.ndarray mat2):
    """ Multiplication of matrices in the order [mat1] x [mat2].
        
        Keyword arguments:
            mat1 -- a 4x4 matrix.
            mat2 -- a 4x4 matrix.
        
        Returns:
            result -- a multiplied 4x4 matrix of [mat1] x [mat2].
    """
    cdef np.ndarray[np.float32_t, ndim=2] result = np.zeros((4,4), dtype=np.float32)
    result[0,0] = mat1[0,0]*mat2[0,0]+mat1[0,1]*mat2[1,0]+mat1[0,2]*mat2[2,0]+mat1[0,3]*mat2[3,0]
    result[1,0] = mat1[1,0]*mat2[0,0]+mat1[1,1]*mat2[1,0]+mat1[1,2]*mat2[2,0]+mat1[1,3]*mat2[3,0]
    result[2,0] = mat1[2,0]*mat2[0,0]+mat1[2,1]*mat2[1,0]+mat1[2,2]*mat2[2,0]+mat1[2,3]*mat2[3,0]
    result[3,0] = mat1[3,0]*mat2[0,0]+mat1[3,1]*mat2[1,0]+mat1[3,2]*mat2[2,0]+mat1[3,3]*mat2[3,0]
    
    result[0,1] = mat1[0,0]*mat2[0,1]+mat1[0,1]*mat2[1,1]+mat1[0,2]*mat2[2,1]+mat1[0,3]*mat2[3,1]
    result[1,1] = mat1[1,0]*mat2[0,1]+mat1[1,1]*mat2[1,1]+mat1[1,2]*mat2[2,1]+mat1[1,3]*mat2[3,1]
    result[2,1] = mat1[2,0]*mat2[0,1]+mat1[2,1]*mat2[1,1]+mat1[2,2]*mat2[2,1]+mat1[2,3]*mat2[3,1]
    result[3,1] = mat1[3,0]*mat2[0,1]+mat1[3,1]*mat2[1,1]+mat1[3,2]*mat2[2,1]+mat1[3,3]*mat2[3,1]
    
    result[0,2] = mat1[0,0]*mat2[0,2]+mat1[0,1]*mat2[1,2]+mat1[0,2]*mat2[2,2]+mat1[0,3]*mat2[3,2]
    result[1,2] = mat1[1,0]*mat2[0,2]+mat1[1,1]*mat2[1,2]+mat1[1,2]*mat2[2,2]+mat1[1,3]*mat2[3,2]
    result[2,2] = mat1[2,0]*mat2[0,2]+mat1[2,1]*mat2[1,2]+mat1[2,2]*mat2[2,2]+mat1[2,3]*mat2[3,2]
    result[3,2] = mat1[3,0]*mat2[0,2]+mat1[3,1]*mat2[1,2]+mat1[3,2]*mat2[2,2]+mat1[3,3]*mat2[3,2]
    
    result[0,3] = mat1[0,0]*mat2[0,3]+mat1[0,1]*mat2[1,3]+mat1[0,2]*mat2[2,3]+mat1[0,3]*mat2[3,3]
    result[1,3] = mat1[1,0]*mat2[0,3]+mat1[1,1]*mat2[1,3]+mat1[1,2]*mat2[2,3]+mat1[1,3]*mat2[3,3]
    result[2,3] = mat1[2,0]*mat2[0,3]+mat1[2,1]*mat2[1,3]+mat1[2,2]*mat2[2,3]+mat1[2,3]*mat2[3,3]
    result[3,3] = mat1[3,0]*mat2[0,3]+mat1[3,1]*mat2[1,3]+mat1[3,2]*mat2[2,3]+mat1[3,3]*mat2[3,3]
    return result

cpdef my_glRotatef(np.ndarray in_matrix, np.float32_t angle, np.ndarray dir_vec):
    """ Produces a rotation matrix of "angle" degrees around the vector 
        "dir_vec", then multiply the input matrix "in_matrix" with this 
        rotation matrix "rot_matrix" in the order [in_matrix] x [trans_matrix].
        The rotation matrix takes the form of:
        
        | x*x*(1-c)+c   x*y*(1-c)-z*s   x*z*(1-c)+y*s   0|
        |y*x*(1-c)+z*s   y*y*(1-c)+c    y*z*(1-c)-x*s   0|
        |x*z*(1-c)-y*s  y*z*(1-c)+x*s    z*z*(1-c)+c    0|
        |      0              0               0         1|
        
        Keyword arguments:
            in_matrix -- a 4x4 matrix that you want to apply a rotation.
            angle -- the angle of rotation in radians.
            dir_vec -- a 3 elements vector with the xyz components to use as
                       axis of rotation.
        
        Returns:
            A multiplied 4x4 matrix of [in_matrix] x [rot_matrix].
    """
    cdef np.ndarray[np.float32_t, ndim=1] vector = np.array(dir_vec, dtype=np.float32)
    assert(np.linalg.norm(vector)>0.0)
    #assert(angle>=0.0)
    angle = angle*np.pi/180.0
    cdef np.float32_t x, y, z
    x, y, z = vector/np.linalg.norm(vector)
    cdef np.float32_t c = np.cos(angle)
    cdef np.float32_t s = np.sin(angle)
    cdef np.ndarray[np.float32_t, ndim=2] rot_matrix = np.identity(4, dtype=np.float32)
    rot_matrix[0,0] = x*x*(1-c)+c
    rot_matrix[1,0] = y*x*(1-c)+z*s
    rot_matrix[2,0] = x*z*(1-c)-y*s
    rot_matrix[0,1] = x*y*(1-c)-z*s
    rot_matrix[1,1] = y*y*(1-c)+c
    rot_matrix[2,1] = y*z*(1-c)+x*s
    rot_matrix[0,2] = x*z*(1-c)+y*s
    rot_matrix[1,2] = y*z*(1-c)-x*s
    rot_matrix[2,2] = z*z*(1-c)+c
    return my_glMultiplyMatricesf(in_matrix, rot_matrix)

cpdef my_glPerspectivef(np.float32_t fovy, np.float32_t aspect, np.float32_t z_near, np.float32_t z_far):
    """ Creates a perspective matrix with "fovy" as field of view, "aspect" as
        viewport aspect ratio, "z_near" as the near clipping plane and "z_far"
        as far clipping plane. The perpective matrix is constructed in the form:
        
        |f/aspect  0                 0                   0|
        |   0      f                 0                   0|
        |   0      0    (z_near+z_far)/(z_near-z_far)   -1|
        |   0      0   2.0*z_near*z_far/(z_near-z_far)   0|
        
        With f = cotangent(fovy/2)
        
    """
    assert(aspect>0.0)
    assert(z_far>z_near)
    cdef np.float32_t f = np.float32(1.0/(np.tan(fovy*np.pi/180.0)))
    cdef np.ndarray[np.float32_t, ndim=2] pers_matrix = np.zeros((4,4), dtype=np.float32)
    pers_matrix[0,0] = f/aspect
    pers_matrix[1,1] = f
    pers_matrix[2,2] = (z_near+z_far)/(z_near-z_far)
    pers_matrix[3,2] = 2.0*z_near*z_far/(z_near-z_far)
    pers_matrix[2,3] = -1.0
    return pers_matrix

cpdef my_glFrustumf(np.float32_t left, np.float32_t rigth, np.float32_t bottom, np.float32_t top, np.float32_t near, np.float32_t far):
    """ Creates a frustrum to use as perspective matrix. Uses the "left" and 
        "right" as leaft and rigth vertical clipping planes of the frustrum.
        The "top" and "bottom" as the top and bottom horizontal clipping planes.
        The "near" and "far" especify the near and far depth clipping planes.
        The resulting matrix is constructed in the form:
        
        |  2.0*near/(rigth-left)                0                         0               0|
        |            0                2.0*near/(top-bottom)               0               0|
        |(rigth+left)/(rigth-left)  (top+bottom)/(top-bottom)   (far+near)/(near-far)    -1|
        |            0                          0              -2.0*near*far/(far-near)   0|
        
    """
    cdef np.ndarray[np.float32_t, ndim=2] frust = np.zeros((4,4), dtype=np.float32)
    frust[0,0] = 2.0*near/(rigth-left)
    frust[1,1] = 2.0*near/(top-bottom)
    frust[2,0] = (rigth+left)/(rigth-left)
    frust[2,1] = (top+bottom)/(top-bottom)
    frust[2,2] = -(far+near)/(far-near)
    frust[2,3] = -1.0
    frust[3,2] = -2.0*near*far/(far-near)
    return frust

cpdef my_glOrthof(np.float32_t left, np.float32_t rigth, np.float32_t bottom, np.float32_t top, np.float32_t near, np.float32_t far):
    """ 
    """
    cdef np.ndarray[np.float32_t, ndim=2] ortho = np.zeros((4,4), dtype=np.float32)
    ortho[0,0] = 2.0/(rigth-left)
    ortho[1,1] = 2.0/(top-bottom)
    ortho[2,2] = -2.0/(far-near)
    ortho[3,0] = (rigth+left)/(left-rigth)
    ortho[3,1] = (top+bottom)/(bottom-top)
    ortho[3,2] = (far+near)/(near-far)
    ortho[3,3] = 1.0
    return ortho

cpdef get_xyz_coords(np.ndarray xyz_mat):
    """ Returns the x, y, z position contained in the xyz_mat matrix. The
        input matrix needs to be a 4x4 matrix.
    """
    assert(xyz_mat.ndim==2)
    assert(xyz_mat.size==16)
    cdef np.ndarray[np.float32_t, ndim=2] rot_mat = xyz_mat[:3,:3]
    cdef np.ndarray[np.float32_t, ndim=1] pos = -xyz_mat[3,:3]
    cdef np.ndarray[np.float32_t, ndim=1] position = pos.dot(rot_mat)
    return position

cpdef get_inverse_matrix(mat):
    """ Function doc
    """
    mat_o = np.matrix(np.copy(mat))
    assert(mat.shape == (4,4))
    return np.array(mat_o.I)

cpdef angle(v0, v1):
    """ Return angle [0..pi] between two vectors."""
    try:
        v0 = np.array(v0, dtype=np.float32)
        v1 = np.array(v1, dtype=np.float32)
        v0 /= np.linalg.norm(v0)
        v1 /= np.linalg.norm(v1)
        return np.acos(np.dot(v0, v1))
    except:
        return False

cpdef dihedral(p0, p1, p2, p3):
    """ Return angle [0..2*pi] formed by vertices p0-p1-p2-p3."""
    p0 = np.array(p0, dtype=np.float32)
    p1 = np.array(p1, dtype=np.float32)
    p2 = np.array(p2, dtype=np.float32)
    p3 = np.array(p3, dtype=np.float32)
    v01 = p0 - p1
    v32 = p3 - p2
    v12 = p1 - p2
    v0 = np.cross(v12, v01)
    v3 = np.cross(v12, v32)
    # The cross product vectors are both normal to the axis
    # vector v12, so the angle between them is the dihedral
    # angle that we are looking for.  However, since "angle"
    # only returns values between 0 and pi, we need to make
    # sure we get the right sign relative to the rotation axis
    a = angle(v0, v3)
    if np.dot(np.cross(v0, v3), v12) > 0:
        a = -a
    return a

