#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
#  vismol_session.py
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

from vismol.utils.debug import dprint
import os
import numpy as np
from vismol.utils import parser
from vismol.utils.elements import PeriodicTable
from logging import getLogger
from vismol.core.vismol_config import VismolConfig
from vismol.core.vismol_selections import VismolPickingSelection as VMPick
from vismol.core.vismol_selections import VismolViewingSelection as VMSele

logger = getLogger(__name__)

class Config:
    """ Class doc """
    
    def __init__ (self):
        """ Class initialiser """
        pass
    def set_pk_sphr_selection_color (self, color, pk = 1):
        """ Function doc """
        
        if pk == 1:
            self.picking_selections.pk_scolor['pk1'] = color
        elif pk == 2:
            self.picking_selections.pk_scolor['pk2'] = color
        elif pk == 3:
            self.picking_selections.pk_scolor['pk3'] = color
        elif pk == 4:
            self.picking_selections.pk_scolor['pk4'] = color
        else:
            return False
        
        
    
    
class VismolSession(Config):
    """ Class doc """
    
    def __init__(self, toolkit, widget=None, main_session=None, vm_config = None):
        """ Class initialiser """
        self.main_session = None
        self.toolkit = toolkit
        self.frame = 0
        
        if vm_config:
            self.vm_config = vm_config
        else:
            self.vm_config = VismolConfig(self)
        
        self.vm_objects_dic = {}
        self.atom_dic_id = {}
        self.atom_id_counter = np.uint32(0) # This variable could be replaced by len(self.atom_dic_id)?
        self.picking_selection_mode = False # True/False  - interchange between viewing  and picking mode
        self.selections = {"sel_00": VMSele(self)}
        self.current_selection = "sel_00"
        self.picking_selections = VMPick(self)
        
        self.periodic_table = PeriodicTable()
        #print(self.periodic_table.elements_by_symbol)
        
        '''
        vm_geometric_object_dic contains items with specific names that 
        are used in selection representations
        
        see vismol_glcore.py  ->  render function
        '''
        self.vm_geometric_object_dic = {"pk1pk2":None, "pk2pk3":None, "pk3pk4":None,
                                        "pk1": None  , "pk2": None  , "pk3": None  , "pk4": None,
                                        "pk_angle": None, "pk_dihedral": None}
        
        if toolkit == "Gtk_3.0":
            from vismol.gui.vismol_gtkwidget import VismolGTKWidget
            self.selection_box_frame = None
            if widget is None:
                self.vm_widget = VismolGTKWidget(self)
            else:
                self.vm_widget = widget
            self.vm_glcore = self.vm_widget.vm_glcore
            self.vm_glcore.queue_draw()
        
        elif toolkit == "Qt5":
            from vismol.gui.vismol_qtwidget import VismolQtWidget
            self.selection_box_frame = None
            if widget is None:
                self.vm_widget = VismolQtWidget(self)
            else:
                self.vm_widget = widget
                self.vm_widget.vm_session = self
                # The widget owns the OpenGL core, created with this session.
                if getattr(self.vm_widget, "vm_glcore", None) is None:
                    from vismol.libgl.vismol_glcore import VismolGLCore
                    self.vm_widget.vm_glcore = VismolGLCore(self.vm_widget, self)
            self.vm_glcore = self.vm_widget.vm_glcore
        
        else:
            self.vm_widget = None
            logger.error("Toolkit not defined or syntax error, try 'Gtk_3.0'. Quitting.")
            raise RuntimeError("Toolkit not defined or syntax error, try 'Gtk_3.0'. Quitting.")
            quit()
        
        
        
    def _add_vismol_object(self, vismol_object, show_molecule=True, autocenter=True):
        """ Function doc """

        # Check if the vismol_object's index already exists in the dictionary
        if vismol_object.index in self.vm_objects_dic.keys():
            logger.warning("The VismolObject with id {} already exists. \
                The data will be overwritten.".format(vismol_object.index))

        # Add the vismol_object to the dictionary with a new index
        self.vm_objects_dic[len(self.vm_objects_dic)] = vismol_object

        # Update the atom_id_counter by adding the number of atoms in the new vismol_object
        self.atom_id_counter += len(vismol_object.atoms)

        # Iterate through the atoms in the vismol_object and add them to the atom_dic_id dictionary
        # using their unique_id as the key and the atom object as the value.
        for atom in vismol_object.atoms.values():
            self.atom_dic_id[atom.unique_id] = atom

        # If show_molecule is True, create representations for "lines" and "nonbonded" to display the object.
        if show_molecule:
            vismol_object.create_representation(rep_type="lines")
            vismol_object.create_representation(rep_type="nonbonded")
            #vismol_object.create_representation(rep_type="dash")
            #vismol_object.create_representation(rep_type="surface")

            # If autocenter is True, center the view on the mass center of the vismol_object.
            if autocenter:
                self.vm_glcore.center_on_coordinates(vismol_object, vismol_object.mass_center)
            else:
                self.vm_glcore.queue_draw()
    
    def load_molecule(self, infile):
        """ Probably would be better to join this with _add_vismol_object
        """
        vismol_object, show_molecule = parser.parse_file(self, infile)
        vismol_object._generate_color_vectors(self.atom_id_counter)
        vismol_object.active = True
        self._add_vismol_object(vismol_object, show_molecule=show_molecule)
    
    def _change_attributes_for_atoms(self, atoms, rep_type, show):
        """ Function doc """
        for atom in atoms:
            try:
                if show:
                    setattr(atom, rep_type, True)
                else:
                    setattr(atom, rep_type, False)
            except AttributeError as ae:
                logger.error("Representation of type {} not implemented".format(rep_type))
                logger.error(ae)
    
    def show_or_hide(self, rep_type="lines", selection=None, show=True):
        """ Function doc """
        if selection is None:
            selection = self.selections[self.current_selection]
            #print('\n\n\selection.get_selection_info()')
            #print(selection.get_selection_info())
            #print('\n\n')
        self._change_attributes_for_atoms(selection.selected_atoms, rep_type, show)
        for vm_object in selection.selected_objects:
            if rep_type == "cartoon":
                # [EN] 2026-09-28 user request: partial cartoon, only for
                # the residues of the selection -- handled apart from the
                # index-list reps below (the ribbon is a generated mesh,
                # see _show_or_hide_cartoon()).
                self._show_or_hide_cartoon(vm_object)
                continue
            if vm_object.representations[rep_type] is None:
                vm_object.create_representation(rep_type = rep_type)
            show_hide_indexes = []
            if rep_type == "lines":
                for bond in vm_object.bonds.values():
                    if bond.atom_i.lines and bond.atom_j.lines:
                        show_hide_indexes.append(bond.atom_index_i)
                        show_hide_indexes.append(bond.atom_index_j)
            
            elif rep_type == "sticks":
                for bond in vm_object.bonds.values():
                    if bond.atom_i.sticks and bond.atom_j.sticks:
                        show_hide_indexes.append(bond.atom_index_i)
                        show_hide_indexes.append(bond.atom_index_j)
            
            elif rep_type == "dash":
                for bond in vm_object.bonds.values():
                    if bond.atom_i.dash and bond.atom_j.dash:
                        show_hide_indexes.append(bond.atom_index_i)
                        show_hide_indexes.append(bond.atom_index_j)
            
            elif rep_type == "dynamic":
                self.define_dynamic_bonds(selection = selection)
                for bond in vm_object.bonds.values():
                    if bond.atom_i.dynamic and bond.atom_j.dynamic:
                        show_hide_indexes.append(bond.atom_index_i)
                        show_hide_indexes.append(bond.atom_index_j)
            
            elif rep_type == "dots":
                for atom in vm_object.atoms.values():
                    if atom.dots:
                        show_hide_indexes.append(atom.atom_id)
            
            elif rep_type == "nonbonded":
                for atom in vm_object.atoms.values():
                    if atom.nonbonded:
                        show_hide_indexes.append(atom.atom_id)
            
            elif rep_type == "impostor":
                for atom in vm_object.atoms.values():
                    if atom.impostor:
                        show_hide_indexes.append(atom.atom_id)
            
            elif rep_type == "spheres":
                for atom in vm_object.atoms.values():
                    if atom.spheres:
                        show_hide_indexes.append(atom.atom_id)
            
            elif rep_type == "stick_spheres":
                for bond in vm_object.bonds.values():
                    if bond.atom_i.stick_spheres and bond.atom_j.stick_spheres:
                        show_hide_indexes.append(bond.atom_index_i)
                        show_hide_indexes.append(bond.atom_index_j)
                
                
                #for atom in vm_object.atoms.values():
                #    if atom.stick_spheres:
                #        show_hide_indexes.append(atom.atom_id)
            
            elif rep_type == "vdw_spheres":
                for atom in vm_object.atoms.values():
                    if atom.vdw_spheres:
                        show_hide_indexes.append(atom.atom_id)
            
            elif rep_type == "ribbons": # add by bachega at 02/02/2023
                for bond in vm_object.c_alpha_bonds:
                    if bond.atom_i.ribbons and bond.atom_j.ribbons:
                        show_hide_indexes.append(bond.atom_index_i)
                        show_hide_indexes.append(bond.atom_index_j)
            
            elif rep_type == "ribbon_sphere":
                for atom in vm_object.atoms.values():
                    if atom.ribbon_sphere and atom.name=='CA' and atom.residue.is_protein:
                        show_hide_indexes.append(atom.atom_id)
                #print('show_hide_indexes',show_hide_indexes)
                #logger.error("Not implementer for 'ribbon' yet.")
                #raise NotImplementedError("Not implementer for 'ribbon' yet.")
            
            
            elif rep_type == "surface":
                logger.error("Not implementer for 'surface' yet.")
                raise NotImplementedError("Not implementer for 'surface' yet.")
            
            elif rep_type == "labels":
                for atom in vm_object.atoms.values():
                    if atom.labels:
                        show_hide_indexes.append(atom.atom_id)
            
            if len(show_hide_indexes) > 0:
                #print('vm_object.representations[rep_type].active = True')
                vm_object.representations[rep_type].define_new_indexes_to_vbo(show_hide_indexes)
                vm_object.representations[rep_type].active = True
                vm_object.representations[rep_type].was_rep_ind_modified = True
                vm_object.representations[rep_type].was_sel_ind_modified = True
                # Coordinates may be stale if this representation was
                # inactive while the trajectory frame advanced (render()
                # only refreshes coords for representations that are
                # active - see vismol_glcore.py), so force a refresh now
                # that it's visible again.
                vm_object.representations[rep_type].was_rep_coord_modified = True
                vm_object.representations[rep_type].was_sel_coord_modified = True
            else:
                vm_object.representations[rep_type].active = False
                #print('vm_object.representations[rep_type].active = False')

        self.vm_widget.queue_draw()
        return selection
        
    def _show_or_hide_cartoon(self, vm_object):
        """ Cartoon counterpart of show_or_hide()'s per-object block. The
            atom.cartoon flags were already set by
            _change_attributes_for_atoms(); a residue is drawn when its CA
            is flagged (see CartoonRepresentation._residue_is_shown()).
            Since the ribbon is a generated mesh (no index list to swap),
            showing/hiding part of it means rebuilding it -- reusing the
            cached secondary structure, which is always computed over the
            FULL chain, so hiding a piece never changes how the rest is
            classified. Creates the representation on first use with
            show_all=False, so only the selection appears. """
        rep = vm_object.representations.get("cartoon")
        if rep is None:
            if not any(getattr(atom, "cartoon", False) for atom in vm_object.atoms.values()):
                return
            vm_object.create_representation(rep_type="cartoon", show_all=False)
            rep = vm_object.representations["cartoon"]
        else:
            rep.rebuild()
        rep.active = rep.has_shown_residues()

    def forward_frame(self):
        """ Function doc """
        #return 0
        frame = self.frame + 1
        for i, vm_object in enumerate(self.vm_objects_dic.values()):
            if frame < vm_object.frames.shape[0]:
                self.frame += 1
                self.vm_glcore.updated_coords = True
                break
            else:
                pass
        else:
            self.vm_glcore.updated_coords = False
        
        #if self.picking_selection_mode:
        #    self.picking_selections.print_pk_distances()
            
    def reverse_frame(self):
        """ Function doc """
        if self.frame - 1 >= 0:
            self.frame -= 1
            self.vm_glcore.updated_coords = True
        else:
            self.vm_glcore.updated_coords = False
        
        #if self.picking_selection_mode:
        #    self.picking_selections.print_pk_distances()
            
    def set_frame(self, frame=0):
        """ Function doc """
        assert frame >= 0
        self.frame = np.uint32(frame)
        
        #self.picking_selections
        #for 
        self.picking_selections.update_pki_pkj_rep_coordinates()
        self.vm_widget.queue_draw()
        #print('\n\n\nuhuuu')
        #if self.picking_selection_mode:
        #    self.picking_selections.print_pk_distances()
        
    def get_frame(self):
        """ Function doc """
        return self.frame
    
    def _selection_function_set(self, selected, _type=None, disable=True):
        """ Function doc """
        #print('selected', selected)
        if self.picking_selection_mode: # True for picking mode
            if selected:
                assert len(selected) == 1
                selected = list(selected)[0]
                self.picking_selections.selection_function_picking(selected)
            else:
                self.picking_selections.selection_function_picking(None)
        else: # False for viewing mode
            # [EN] 2026-09-25, user's own explicit request: "enquanto o
            # modo 'editing' estiver 'on', nao mostrar nenhuma selecao.
            # Quando o builder estiver aberto e o modo editing estiver
            # desligado, selecoes sao permitidas." -- builder_atom_mode
            # is only ever True while the Builder sidebar is open AND
            # its own "Editing: ON/OFF" toggle is on (see builder_
            # sidebar.py's open_window()/close_window()/on_placemode_
            # toggled()) -- False both while the Builder is closed and
            # while it's open but paused, so a single check here already
            # covers exactly the 2 cases the request describes. No-ops
            # instead of calling selection_function_viewing_set() at all,
            # so shift+click/shift+drag genuinely can't CHANGE the
            # viewing selection while actively editing -- not just a
            # rendering suppression (see vm_glcore.py's render(), which
            # ALSO hides whatever selection already existed from before
            # Editing was turned back on).
            if getattr(self, "builder_atom_mode", False):
                return
            self.selections[self.current_selection].selection_function_viewing_set(selected, _type, disable)
    
    def viewing_selection_mode(self, sel_type="atom"):
        """ Function doc """
        if self.selection_box_frame:
            self.selection_box_frame.change_sel_type_in_combobox(sel_type)
        self.selections[self.current_selection].selection_mode = sel_type
    
    def define_dynamic_bonds (self, selection = None):
        """ Function doc """
        if selection:
            pass
        else:
            selection = self.selections[self.current_selection]
        
        selection_dict = {}
        vobject = None
        for atom in selection.selected_atoms:
            selection_dict[atom.atom_id] = atom
            vobject = atom.vm_object
        
        #tolerance = self.vm_config.gl_parameters["dynamic_bond_tolerance"]
        tolerance = self.vm_config.gl_parameters['bond_tolerance']
        
        #tolerance = 1.05
        #tolerance = 1.2
        
        vobject.dynamic_bonds = []
        for frame in range(len(vobject.frames)):            
            bonds = vobject.find_bonded_and_nonbonded_atoms(selection=selection_dict, 
                                                            frame=frame, internal = False, 
                                                            tolerance = tolerance)
            
            '''       
            new_bonds = []

            for index in range (0,len(bonds), 2):
                dprint(bonds[index], bonds[index+1])
                
                if [bonds[index], bonds[index+1]] in new_bonds or [bonds[index+1], bonds[index]] in new_bonds:
                    pass
                else:
                    new_bonds.append([bonds[index], bonds[index+1]])
            dprint(new_bonds)
            #bonds = []
            bonds = [item for sublista in new_bonds for item in sublista]
            #for b in new_bonds:
            #    bonds.append(b[0], b[1])
            '''
            #bonds.append(0)
            #bonds.append(0)
            #bonds.append(0)
            #bonds.append(0)
            #
            #print( bonds)

            vobject.dynamic_bonds.append(bonds)
            #print(len(bonds), bonds)
        #print(vobject.dynamic_bonds)

    def hide_axes (self):
        """ Function doc """
        self.vm_glcore.show_axis = False
    def show_axes (self):
        """ Function doc """
        self.vm_glcore.show_axis = True

        
    
    # def delete_vismol_object_by_index(self, index):
    #     """ Function doc
    #     """
    #     try:
    #         vm_object = self.vm_objects_dic.pop(index)
    #     except KeyError:
    #         logger.warning("VismolObject with index {} not found".format(index))
    #         return None
    #     return vm_object
    
    # def select(self, vismol_object=None, indexes=None, sele=None):
    #     """ Function doc """
    #     if vismol_object is None:
    #         vismol_object = self.vm_objects[-1]
        
    #     if sele is None:
    #         sele = self.current_selection
        
    #     if indexes == "all":
    #         self.selections[sele].selecting_by_indexes(vismol_object=vismol_object,
    #                                                    indexes=range(0, int(len(vismol_object.atoms)/2)))
    #     self.vm_widget.queue_draw()
    
    # def disable_by_index(self, index):
    #     """ When the variable "dictionary" is active, the function accesses 
    #         a vismol object through the dictionary "self.vm_objects_dic". 
    #         Each vismol object has a unique access key (int), which, in 
    #         easyhybrid, is generated in the method: add_vismol_object.
            
    #         In the vismol interface the enable_by_index/disable_by_index methods
    #         access the vismol objects by their position in the "self.vm_objects" 
    #         list (this is because when an object is deleted in the vismol 
    #         interface, the treeview"s liststore is rewritten)
    #     """
    #     try:
    #         self.vm_objects_dic[index].active = False
    #         self.vm_glcore.queue_draw()
    #     except KeyError:
    #         logger.error("VismolObject with index {} not found".format(index))
    #         return False
    #     return True
    
    # def enable_by_index(self, index):
    #     """ When the variable "dictionary" is active, the function accesses 
    #         a vismol object through the dictionary "self.vm_objects_dic". 
    #         Each vismol object has a unique access key (int), which, in 
    #         easyhybrid, is generated in the method: add_vismol_object.
            
    #         In the vismol interface the enable_by_index/disable_by_index methods
    #         access the vismol objects by their position in the "self.vm_objects" 
    #         list (this is because when an object is deleted in the vismol 
    #         interface, the treeview"s liststore is rewritten)
    #     """
    #     try:
    #         self.vm_objects_dic[index].active = True
    #         self.vm_glcore.queue_draw()
    #     except KeyError:
    #         logger.error("VismolObject with index {} not found".format(index))
    #         return False
    #     return True
    
    # def edit_by_index(self, index):
    #     """ Function doc
    #     """
    #     try:
    #         self.vm_objects_dic[index].editing = not self.vm_objects_dic[index].editing
    #         self.vm_glcore.queue_draw()
    #     except KeyError:
    #         logger.error("VismolObject with index {} not found".format(index))
    #         return False
    #     return True
    
    # def set_color_by_index(self, vismol_object, indexes=None, color=None):
    #     """ NOT SURE WHAT THIS FUNCTION DOES
    #     """
    #     if indexes is None:
    #         indexes = []
    #     if color is None:
    #         color = np.array([0.9, 0.9, 0.9], dtype=np.float32)
        
    #     for atom_index in indexes:
    #         vismol_object.atoms[atom_index].color = color
    #     vismol_object._generate_color_vectors(do_colors=True, do_colors_idx=False,
    #                                           do_colors_raindow=False, do_vdw_dot_sizes=False,
    #                                           do_cov_dot_sizes=False)
    #     self.vm_widget.queue_draw()
    #     for rep  in vismol_object.representations.keys():
    #         if vismol_object.representations[rep]:
    #             vismol_object.representations[rep]._set_colors_to_buffer()
    #             # try:
    #             #     vismol_object.representations[rep]._set_colors_to_buffer()
    #             # except:
    #             #     print('"VisMol/vModel/Representations.py, line 123, in _set_colors_to_buffer GL.glBindBuffer(GL.GL_ARRAY_BUFFER, ctypes.ArgumentError: argument 2: <class "TypeError">: wrong type"')
    #     return True
    
