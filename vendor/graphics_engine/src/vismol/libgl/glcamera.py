#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
#  glcamera.py
#  
#  Copyright 2016 Carlos Eduardo Sequeiros Borja <casebor@gmail.com>
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

import math
import numpy as np
import vismol.utils.matrix_operations as mop
from logging import getLogger

logger = getLogger(__name__)


class GLCamera():
    """ The GLCamera object creates a "camera" to be used in OpenGL.
        It automatically creates a viewing and projection matrices with the
        values defined in the constructor.
    """
    
    def __init__ (self, fov=30.0, var=(4.0 / 3.0), pos=np.array([0,0,10], dtype=np.float32),
                  zrp=np.array([0,0,0], dtype=np.float32)):
        """ Depending of the distance from the camera position to a defined
            reference point, it creates different clipping planes (this function
            could be improved, but will work for now).
            
            Input parameters:
                fov -- Specifies the field of view angle, in degrees, in the
                       Y direction.
                var -- Specifies the aspect ratio that determines the field of
                       view in the x direction. The aspect ratio is the ratio
                       of x (width) to y (height).
                pos -- The position in world coordinates of the camera.
                zrp -- The zero reference point defined for rotation functions.
            
            Automatically generates:
                z_near -- Specifies the distance from the viewer to the near
                          clipping plane (always positive).
                z_far -- Specifies the distance from the viewer to the far
                         clipping plane (always positive).
                fog_start -- Specifies the beggining of the fog effect (always
                             positive and lower than z_far).
                fog_end -- Specifies the end of the fog effect (always positive
                           and preferable equal to z_far).
                view_matrix -- The viewing matrix used in the render by shaders.
                projection_matrix -- The projection matrix used in the render
                                     by shaders.
        """
        self.field_of_view = np.float32(fov)
        self.viewport_aspect_ratio = np.float32(var)
        self.zero_reference_point = np.array(zrp, dtype=np.float32)
        self.max_vertical_angle = 85.0  # must be less than 90 to avoid gimbal lock
        self.horizontal_angle = 0.0
        self.vertical_angle = 0.0
        self.min_znear = 0.1
        self.min_zfar = 9.0
        # [EN] 2026-09-26 user request ("um parametro para a extensao do
        # fog") -- multiplies min_zfar ONLY for the fog-range computation
        # below (see update_fog()), decoupled from min_zfar's OTHER role
        # as the hard zoom-out floor elsewhere (mouse_scroll()'s slab
        # guard, _zoom_view()'s own clamp) -- changing min_zfar itself, as
        # this class's own (now outdated) update_fog() docstring used to
        # suggest, would have side-affected zoom behaviour too. 1.0 =
        # original behaviour (fog span exactly min_zfar); >1 = fog starts
        # closer to the camera and reaches further into the scene (more
        # pronounced depth cue); <1 = fog stays confined closer to z_far.
        self.fog_range_scale = 1.0
        # [EN] 2026-09-28 user request (Preferences > General > Fog: "mexer
        # na intensidade, onde comeca e onde termina"). Defaults reproduce
        # the previous look exactly (fog on, full intensity, automatic
        # range) -- see update_fog() and VismolGLCore.load_fog().
        #   fog_enabled       -- False = no fog at all (intensity sent as 0)
        #   fog_intensity     -- 0..1, how much fog is reached at fog_end
        #                        (1.0 = fully fog-coloured, the old behaviour)
        #   fog_custom_range  -- False = automatic range (old formula);
        #                        True = fog_start/end_percent below
        #   fog_start_percent / fog_end_percent -- 0..100 % of the VISIBLE
        #                        depth (0 = near clipping plane, 100 = far
        #                        clipping plane), so the fog follows zoom/
        #                        slab changes instead of fixed Angstroms.
        self.fog_enabled = True
        self.fog_intensity = 1.0
        self.fog_custom_range = False
        self.fog_start_percent = 50.0
        self.fog_end_percent = 100.0
        # [EN] 2026-09-29 (user request: "seria interessante o fog comecar
        # antes do ponto onde a camera esta centrada"): reference for the
        # AUTOMATIC range (custom range unaffected) --
        #   "far_plane" (default, the original behaviour): starts
        #       min_zfar * fog_range_scale in front of the far plane;
        #   "center": starts fog_center_offset A in front of the rotation
        #       center and ends at the far plane.
        self.fog_start_reference = "far_plane"
        self.fog_center_offset = 5.0
        dist = np.linalg.norm(pos - zrp)
        if dist <= 10.0:
            self.z_near = dist - 3.0
            self.z_far = dist + 3.0
        elif dist <= 20.0:
            self.z_near = dist - 6.0
            self.z_far = dist + 6.0
        elif dist <= 40.0:
            self.z_near = dist - 9.0
            self.z_far = dist + 9.0
        elif dist <= 80.0:
            self.z_near = dist - 12.0
            self.z_far = dist + 12.0
        else:
            self.z_near = dist - 15.0
            self.z_far = dist + 15.0
        
        self.fog_end = self.z_far
        self.fog_start = self.fog_end - self.min_zfar
        self.view_matrix = self._get_view_matrix(pos)
        self.projection_matrix = self._get_projection_matrix()
    
    def _get_view_matrix(self, position):
        """ Creates the view matrix, i.e. the matrix with the position and
            orientation of the camera used in OpenGL.
            First creates a translation matrix at the position defined by the
            "position" argument, then creates a orientation matrix initially
            with the vertical angle, then the horizontal angle is applied to
            this matrix. Finally the multiplication of translation x orientation
            is returned as the view matrix.
            
            Input parameters:
                position -- a numpy array of size 3 and ndim of 1 containing
                            the camera's position in XYZ coordinates.
            
            Returns:
                view -- a numpy array of size 16 and ndim of 2 containing the
                        information for the camera's position and orientation.
        """
        trans = mop.my_glTranslatef(np.identity(4,dtype=np.float32), -position)
        orient = mop.my_glRotatef(np.identity(4,dtype=np.float32), self.vertical_angle, np.array([1,0,0]))
        orient = mop.my_glRotatef(orient, self.horizontal_angle, np.array([0,1,0]))
        view = mop.my_glMultiplyMatricesf(trans, orient)
        return view
    
    def _get_projection_matrix(self):
        """ Creates the projection matrix using the parameters supplied in the
            constructor of the class. This function creates a perspective with
            a defined field of view, viewport aspect ratio, near and far depth
            clipping planes.
            
            Returns:
                persp -- a numpy array of size 16 and ndim of 2, containing the
                         information for the camera's field and depth view.
        """
        assert(self.field_of_view>0.0 and self.field_of_view<180.0)
        assert(self.z_near>0.0)
        assert(self.z_near<self.z_far)
        assert(self.viewport_aspect_ratio>0.0)
        persp = mop.my_glPerspectivef(self.field_of_view,self.viewport_aspect_ratio,self.z_near,self.z_far)
        return persp
    
    def get_position(self):
        """ Returns the x, y, z position of the camera in 
            absolute coordinates.
        """
        return mop.get_xyz_coords(self.view_matrix)
    
    def get_modelview_position(self, model_matrix):
        #mais atual        
        modelview = mop.my_glMultiplyMatricesf(model_matrix, self.view_matrix)
        #crd_xyz = -1 * np.mat(modelview[:3,:3]) * np.mat(modelview[3,:3]).T
        #return crd_xyz.A1
        A = modelview[:3,:3]
        T = modelview[3,:3]
        crd_xyz = A @ T
        #crd_xyz = modelview[:3,:3] @ modelview[3,:3]
        return crd_xyz 
        
    def _normalize_angles(self):
        """ DEPRECATED FUNCTION??? SEEMS TO NOT BE USED ANYWHERE
        """
        self.horizontal_angle = self.horizontal_angle % 360.0
        if self.horizontal_angle<0:
            self.horizontal_angle += 360.0
        if self.vertical_angle>self.max_vertical_angle:
            self.vertical_angle = self.max_vertical_angle
        elif self.vertical_angle<-self.max_vertical_angle:
            self.vertical_angle = -self.max_vertical_angle
    
    def add_orientation_angles(self, h_angle, v_angle):
        """ DEPRECATED FUNCTION??? SEEMS TO NOT BE USED ANYWHERE
        """
        self.horizontal_angle += h_angle
        self.vertical_angle += v_angle
        self._normalize_angles()
        return True
    
    def look_at(self, target):
        """ DEPRECATED FUNCTION??? SEEMS TO NOT BE USED ANYWHERE
        """
        position = self.get_position()
        assert(position[0]!=target[0] and
               position[1]!=target[1] and
               position[2]!=target[2])
        direction = target - position
        direction /= np.linalg.norm(direction)
        self.vertical_angle = -math.asin(direction[1])*180/math.pi
        self.horizontal_angle = -(math.atan2(-direction[0], -direction[2])*180/math.pi)
        self._normalize_angles()
        return True
    
    def get_proj_view_matrix(self):
        """ Returns:
                proview -- a numpy array of size 16 and ndim of 2, containing
                           the information for the camera's projection-view
                           matrix.
        """
        proview = mop.my_glMultiplyMatricesf(self.get_projection_matrix(), self.get_view_matrix())
        return proview
    
    def set_view_matrix(self, new_view_matrix):
        """ Sets a new matrix as view matrix for the camera.
            
            Input parameters:
                new_view_matrix -- a numpy array of size 16 and ndim of 2
                                   containing the new view matrix.
            
            Returns:
                True
        """
        self.view_matrix = new_view_matrix
        return True
    
    def set_projection_matrix(self, new_proj_matrix):
        """ Sets a new matrix as projection matrix for the camera.
            
            Input parameters:
                new_proj_matrix -- a numpy array of size 16 and ndim of 2
                                   containing the new projection matrix.
            
            Returns:
                True
        """
        self.projection_matrix = new_proj_matrix
        return True
    
    def update_projection(self):
        """ Updates the projection matrix. If you change any parameter of the
            camera, you should use this function right after so the changes can
            be applied.
            
            Returns:
                True
        """
        self.projection_matrix = self._get_projection_matrix()
        return True
    
    def update_fog(self):
        """ Updates automatically the fog. This function was created to avoid
            errors when the fog values were changed manually, this way, the fog
            start and end will always have constant values. To expand or
            decrease the fog's own extent without also affecting min_zfar's
            other role as the zoom-out floor, set self.fog_range_scale
            instead (see its own docstring in __init__).

            Returns:
                True
        """
        if getattr(self, "fog_custom_range", False):
            # effective near plane: the stored z_near may be negative after
            # widening the clipping slab; the projection clamps it to
            # min_znear, so the visible depth starts there too.
            near = max(float(self.z_near), float(self.min_znear))
            depth = max(float(self.z_far) - near, 1e-3)
            start = min(max(float(self.fog_start_percent), 0.0), 100.0)
            end = min(max(float(self.fog_end_percent), 0.0), 100.0)
            # every shader divides by (fog_end - fog_start): keep a
            # strictly positive span even if start >= end was typed in.
            if end <= start:
                end = min(start + 1.0, 100.0)
                start = end - 1.0
            self.fog_start = near + depth * start / 100.0
            self.fog_end = near + depth * end / 100.0
        elif getattr(self, "fog_start_reference", "far_plane") == "center":
            # The visible slab is kept centered on the rotation center by
            # every zoom/slab operation (verified: its midpoint always
            # equals VismolGLCore.dist_cam_zrp), so the center's distance
            # from the camera is simply the slab midpoint.
            center = 0.5 * (float(self.z_near) + float(self.z_far))
            self.fog_end = self.z_far
            self.fog_start = min(center - float(self.fog_center_offset), self.fog_end - 0.01)
        else:
            self.fog_end = self.z_far
            self.fog_start = self.fog_end - self.min_zfar * self.fog_range_scale
        return True

    def get_effective_fog_intensity(self):
        """ Value uploaded as the shaders' fog_intensity uniform: 0 when
            the fog is disabled, otherwise fog_intensity clamped to 0..1. """
        if not getattr(self, "fog_enabled", True):
            return 0.0
        return float(min(max(getattr(self, "fog_intensity", 1.0), 0.0), 1.0))

    def apply_fog_parameters(self, gl_parameters):
        """ Reads the fog_* keys of gl_parameters (Preferences > General
            > Fog) into this camera and recomputes the fog range. Missing
            keys keep the defaults, i.e. the original look. """
        self.fog_range_scale = float(gl_parameters.get("fog_range_scale", 1.0))
        self.fog_enabled = bool(gl_parameters.get("fog_enabled", True))
        self.fog_intensity = float(gl_parameters.get("fog_intensity", 100.0)) / 100.0
        self.fog_custom_range = bool(gl_parameters.get("fog_custom_range", False))
        self.fog_start_percent = float(gl_parameters.get("fog_start_percent", 50.0))
        self.fog_end_percent = float(gl_parameters.get("fog_end_percent", 100.0))
        self.fog_start_reference = gl_parameters.get("fog_start_reference", "far_plane")
        self.fog_center_offset = float(gl_parameters.get("fog_center_offset", 5.0))
        self.update_fog()
    
    def print_parms(self):
        """ Prints camera parameters in the terminal. Method created only for
            debugging purposes. It will come out in the final distribution?
        """
        logger.debug("######## GLCAMERA PARAMETERS ########")
        logger.debug("<= z_near    => {}".format(self.z_near))
        logger.debug("<= z_far     => {}".format(self.z_far))
        logger.debug("<= fog_start => {}".format(self.fog_start))
        logger.debug("<= fog_end   => {}".format(self.fog_end))
        logger.debug("<= position  => {}".format(self.get_position()))
        logger.debug("######## GLCAMERA PARAMETERS ########")
        return True
    
    def print_matrices(self):
        """ Prints camera matrices in the terminal. Method created only for
            debugging purposes. It will come out in the final distribution?
        """
        logger.debug("######## GLCAMERA MATRICES ########")
        logger.debug("<= view_matrix => {}".format(self.view_matrix))
        logger.debug("<= projection_matrix => {}".format(self.projection_matrix))
        logger.debug("######## GLCAMERA MATRICES ########")
        return True
    
