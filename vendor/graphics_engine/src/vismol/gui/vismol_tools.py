#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#

from vismol.utils.debug import dprint
import gi
gi.require_version("Gtk", "3.0")
from gi.repository import Gtk, Gdk

class VismolStatusBar(Gtk.Statusbar):
    """ Class doc """
    
    def __init__ (self, vismol_session = None):
        """ Class initialiser """
        
        pass
        self.vm_session = vismol_session
        self.statusbar = Gtk.Statusbar()

class VismolGoToAtomWindow2(Gtk.Window):
    def OpenWindow (self):
        """ Function doc """
        if self.Visible  ==  False:
            
            self.vm_session.vm_objects_list_store
            
            #------------------------------------------------------------------#
            #                  VISOBJ combobox and Label
            #------------------------------------------------------------------#
            self.box_vertical    = Gtk.Box(orientation = Gtk.Orientation.VERTICAL,   spacing = 6)
            self.box_horizontal1 = Gtk.Box(orientation = Gtk.Orientation.HORIZONTAL, spacing = 6)
            
            
            self.label1  = Gtk.Label()
            self.label1.set_text('Object:')
            self.box_horizontal1.pack_start(self.label1, False, False, 0)

            combobox_vobjects = Gtk.ComboBox.new_with_model(self.vm_session.vm_objects_list_store)
            combobox_vobjects.connect("changed", self.on_combobox_vobjects_changed)
            renderer_text = Gtk.CellRendererText()
            combobox_vobjects.pack_start(renderer_text, True)
            combobox_vobjects.add_attribute(renderer_text, "text", 2)
            #
            #vbox.pack_start(combobox_vobjects, False, False, True)

            self.box_horizontal1.pack_start(combobox_vobjects, False, False, 0)
            #------------------------------------------------------------------#
            
            
            
            #------------------------------------------------------------------#
            #                  CHAIN combobox and Label
            #------------------------------------------------------------------#
            self.box_horizontal2 = Gtk.Box(orientation = Gtk.Orientation.HORIZONTAL, spacing = 6)
            
            
            self.label2  = Gtk.Label()
            self.label2.set_text('Chain:')
            self.box_horizontal2.pack_start(self.label2, False, False, 0)

            self.liststore_chains = Gtk.ListStore(str)
            
            self.combobox_chains = Gtk.ComboBox.new_with_model(self.liststore_chains)
            self.combobox_chains.connect("changed", self.on_combobox_chains_changed)
            renderer_text = Gtk.CellRendererText()
            self.combobox_chains.pack_start(renderer_text, True)
            self.combobox_chains.add_attribute(renderer_text, "text", 0)
            #vbox.pack_start(self.combobox_chains, False, False, True)
            self.box_horizontal2.pack_start(self.combobox_chains, False, False, 0)
            #------------------------------------------------------------------#
            #                  RESIDUES combobox and Label
            #------------------------------------------------------------------#
            
            self.label3  = Gtk.Label()
            self.label3.set_text('Residue type:')
            self.box_horizontal2.pack_start(self.label3, False, False, 0)

            self.liststore_residues = Gtk.ListStore(str)
            
            self.combobox_residues = Gtk.ComboBox.new_with_model(self.liststore_residues)
            self.combobox_residues.connect("changed", self.on_combobox_residues_changed)
            renderer_text = Gtk.CellRendererText()
            self.combobox_residues.pack_start(renderer_text, True)
            self.combobox_residues.add_attribute(renderer_text, "text", 0)
            #vbox.pack_start(self.combobox_chains, False, False, True)
            self.box_horizontal2.pack_start(self.combobox_residues, False, False, 0)
            
            
            #------------------------------------------------------------------#
            
            
            
            
            
            self.treeviewbox_horizontal = Gtk.Box(orientation = Gtk.Orientation.HORIZONTAL, spacing = 6)
            
            #------------------------------------------------------------------------------------------
            #self.treeview = Gtk.TreeView(model =self.residue_liststore)
            
            
            #-----------------------------------------------------------------------------------------
            #                                      Chain filter
            #-----------------------------------------------------------------------------------------
            self.current_filter_chain = None
            # Creating the filter, feeding it with the liststore model
            self.chain_filter = self.residue_liststore.filter_new()
            # setting the filter function, note that we're not using the
            self.chain_filter.set_visible_func(self.chain_filter_func)
            #-----------------------------------------------------------------------------------------
            
            
            self.treeview = Gtk.TreeView(model = self.chain_filter)
            self.treeview.get_selection().set_mode(Gtk.SelectionMode.MULTIPLE)
            
            self.treeview.connect("button-release-event", self.on_treeview_Objects_button_release_event)
            self.treeview.connect("row-activated", self.on_treeview_row_activated_event)
            
            for i, column_title in enumerate(
                ['', "index", "Residue",  "Chain", 'size']
            ):
                if i == 0:
                    cell = Gtk.CellRendererToggle()
                    cell.set_property('activatable', True)
                    cell.connect('toggled', self.on_chk_renderer_toggled, self.residue_liststore)
                    column = Gtk.TreeViewColumn(column_title, cell )
                    column.add_attribute(cell, 'active', 0)
                    self.treeview.append_column(column)
                    #print ('aqui')
                else:
                    renderer = Gtk.CellRendererText()
                    #renderer.connect('toggled', self.on_chk_renderer_toggled, self.residue_liststore)
                    column = Gtk.TreeViewColumn(column_title, renderer, text=i)
                    self.treeview.append_column(column)

            
            self.current_filter_chain = None
            

            self.scrollable_treelist = Gtk.ScrolledWindow()
            self.scrollable_treelist.set_vexpand(True)
            self.scrollable_treelist.add(self.treeview)
            #------------------------------------------------------------------------------------------
            


            
            
            #------------------------------------------------------------------------------------------
            self.treeview_atom = Gtk.TreeView(model =self.atom_liststore)
            self.treeview_atom.connect("button-release-event", self.on_treeview_atom_button_release_event)
            self.treeview_atom.connect("row-activated", self.on_treeview_atom_row_activated_event)

            for i, column_title in enumerate(
                ['', "index", "name",  "type", 'mass']
            ):
                if i == 0:
                    cell = Gtk.CellRendererToggle()
                    cell.set_property('activatable', True)
                    cell.connect('toggled', self.on_chk_renderer_toggled, self.atom_liststore)
                    column = Gtk.TreeViewColumn(column_title, cell )
                    column.add_attribute(cell, 'active', 0)
                    self.treeview_atom.append_column(column)
                    #print ('aqui')
                else:
                    renderer = Gtk.CellRendererText()
                    column = Gtk.TreeViewColumn(column_title, renderer, text=i)
                    self.treeview_atom.append_column(column)

            self.scrollable_treelist2 = Gtk.ScrolledWindow()
            self.scrollable_treelist2.set_vexpand(True)
            self.scrollable_treelist2.add(self.treeview_atom)
            #------------------------------------------------------------------------------------------
            
            
            
            
            
            
            
            
            
            
            
            self.box_vertical.pack_start(self.box_horizontal1, False, True, 0)
            self.box_vertical.pack_start(self.box_horizontal2, False, True, 0)
            self.treeviewbox_horizontal.pack_start(self.scrollable_treelist, True, True, 0)
            self.treeviewbox_horizontal.pack_start(self.scrollable_treelist2, True, True, 0)
            
            self.box_vertical.pack_start(self.treeviewbox_horizontal, False, True, 0)


            combobox_vobjects.set_active(0)
            self.window =  Gtk.Window()
            self.window.set_border_width(10)
            self.window.set_default_size(600, 600)  
            self.window.add(self.box_vertical)
            self.window.connect("destroy", self.CloseWindow)
            self.window.show_all()                                               
            #                                                                
            #self.builder.connect_signals(self)                                   
            
            self.Visible  =  True
            #self.PutBackUpWindowData()
            #gtk.main()
            #----------------------------------------------------------------

    def CloseWindow (self, button):
        """ Function doc """
        #self.BackUpWindowData()
        self.window.destroy()
        self.Visible    =  False
        dprint('self.Visible',self.Visible)
    
    def __init__(self, vismol_session = None):
        """ Class initialiser """
        self.vm_session = vismol_session
        #if EasyHybridSession != None:
        #    self.project          = EasyHybridSession.project
        #    self.main_builder     = EasyHybridSession.builder
        #    self.EasyHybridSession = EasyHybridSession        
        #    self.window_control   = EasyHybridSession.window_control
        #
        #self.atom1_index = ''
        #self.name1       = ''
        #self.atom2_index = ''
        #self.name2       = ''
        #self.atom3_index = ''
        #self.name3       = ''
        #
        #
        #self.distanceType  = 0
        #self.DMINIMUM      = ''
        #self.minitype      = 0
        #
        #self.DINCREMENT    = '0.1'
        #self.NWINDOWS      = '10'
        #self.FORCECONSTANT = '4000'
        #self.max_int       = '500'
        #self.rms_grad      = '0.1'
        #
        #self._mass_weight_check = False
        #
        #
        #self.project   =  project
        self.Visible    =  False
        
        self.residue_liststore = Gtk.ListStore(bool, int, str, str, int)
        self.atom_liststore    = Gtk.ListStore(bool, int, str, str, int, int, )
        self.residue_filter    = False


    def on_combobox_residues_changed (self, widget):
        """ Function doc """
        tree_iter = widget.get_active_iter()
        if tree_iter is not None:
            model = widget.get_model()
            residue = model[tree_iter][0]
            #print("Selected: country=%s" % country)
        
            self.current_filter_residue = residue
            dprint("%s Chain selected!" % self.current_filter_residue)
            # we update the filter, which updates in turn the view
            if self.residue_filter:
                self.residue_filter.refilter()
        
        
    def on_combobox_chains_changed (self, widget):
        """ Function doc """
        ##---------------------------------------------------------------
        #self.current_filter_chain = None
        ## Creating the filter, feeding it with the liststore model
        #self.chain_filter = self.residue_liststore.filter_new()
        ## setting the filter function, note that we're not using the
        #self.chain_filter.set_visible_func(self.chain_filter_func)
        ##---------------------------------------------------------------
        
        tree_iter = self.combobox_chains.get_active_iter()
        if tree_iter is not None:
            model = self.combobox_chains.get_model()
            chain = model[tree_iter][0]
            #print("Selected: country=%s" % country)
        
        self.current_filter_chain = chain
        dprint("%s Chain selected!" % self.current_filter_chain)
        # we update the filter, which updates in turn the view
        self.chain_filter.refilter()
    
    
    def on_combobox_vobjects_changed (self, widget):
        """ Function doc """
        #print(widget)
        #print(widget.get_active())
        
        #print(widget.get_active_id())
        #print(widget.get_active_iter())
        
        #self.vismol_session.vm_objects_dic.items()
        
        #self.vismol_session.vm_objects_dic.items()
        self.VObj = self.vm_session.vm_objects_dic[widget.get_active()]
        #self.VObj = self.vismol_session.vismol_objects[widget.get_active()]
        
        
        self.liststore_chains = Gtk.ListStore(str)
        self.liststore_chains.append(['all'])
        chains = self.VObj.chains.keys()
        
        #self.chain_liststore = Gtk.ListStore(str)
        
        for chain in chains:
            self.liststore_chains.append([chain])
        self.combobox_chains.set_model(self.liststore_chains)
        self.combobox_chains.set_active(0)
        
        
        self.residue_liststore = Gtk.ListStore(bool, int, str, str, int)
        for chain in self.VObj.chains:
            for res in self.VObj.chains[chain].residues:
                #print(res.resi, res.resn, chain,  len(res.atoms) ) 
                
                self.residue_liststore.append(list([True, res.resi, res.resn, chain,  len(res.atoms)]))
        #-----------------------------------------------------------------------------------------
        #                                      Chain filter
        #-----------------------------------------------------------------------------------------
        self.current_filter_chain = None
        # Creating the filter, feeding it with the liststore model
        self.chain_filter = self.residue_liststore.filter_new()
        # setting the filter function, note that we're not using the
        self.chain_filter.set_visible_func(self.chain_filter_func)
        #-----------------------------------------------------------------------------------------
        
        
        
        
        
        
        #-----------------------------------------------------------------------------------------
        #                                      Residue combobox
        #-----------------------------------------------------------------------------------------
        self.liststore_residues = Gtk.ListStore(str)
        self.liststore_residues.append(['all'])
        
        resn_labels = {}
        
        for residue in self.VObj.residues:
            resn_labels[residue.resn] = True
        
        for resn in resn_labels.keys():
            #print (resn)
            self.liststore_residues.append([resn])
        
        self.combobox_residues.set_model(self.liststore_residues)
        self.combobox_residues.set_active(0)
        
        #-----------------------------------------------------------------------------------------
        #                                      Residue filter
        #-----------------------------------------------------------------------------------------
        self.current_filter_residue = None
        # Creating the filter, feeding it with the liststore model
        self.residue_filter = self.chain_filter.filter_new()
        # setting the filter function, note that we're not using the
        self.residue_filter.set_visible_func(self.residue_filter_func)
        #-----------------------------------------------------------------------------------------        
        
        
        #self.treeview.set_model(self.residue_liststore)
        #self.treeview.set_model(self.chain_filter)
        self.treeview.set_model(self.residue_filter)
        
        
    def on_treeview_atom_button_release_event(self, tree, event):
        if event.button == 2:
            selection     = tree.get_selection()
            model         = tree.get_model()
            (model, iter) = selection.get_selected()
            
            
            if iter != None:
                self.selectedID  = int(model.get_value(iter, 1))-1  # @+
                atom = self.VObj.atoms[self.selectedID]
                self.vm_session.glwidget.vm_widget.center_on_atom(atom)
       
    def on_treeview_atom_row_activated_event (self, tree, rowline , column):
        """ Function doc """
        #print (A,B,C)
        selection     = tree.get_selection()
        model         = tree.get_model()
        
        #print(model)
        #print(rowline, list(model[rowline]))
        
        data  = list(model[rowline])
        #print(data)
        #self.selectedID  = int(data[1])  # @+
        #self.selectedObj = str(data[2])
        pickedID = data[-1]
        
        atom_picked = self.vm_session.atom_dic_id[pickedID]
        
        
        #res = self.VObj.chains[self.selectedChn].residues_by_index[self.selectedID]
        #
        ##print('Selecting by doble click')
        ##self.vismol_session.selections[self.vismol_session.current_selection].selecting_by_residue (res.atoms[0])
        ##self.vismol_session.selections[self.vismol_session.current_selection].selection_function_viewing (res.atoms[0])
       
        self.vm_session._selection_function (atom_picked, _type = 'atom')
        self.vm_session.glwidget.queue_draw()
        #
        #self.atom_liststore.clear()
        #for atom in res.atoms:
        #    self.atom_liststore.append(list([True, int(atom.index), atom.name, atom.symbol, atom.charge]))
    
    
    def on_treeview_row_activated_event(self, tree, rowline , column ):
        #print (A,B,C)
        selection     = tree.get_selection()
        model         = tree.get_model()
        
        #print(model)
        #print(rowline, list(model[rowline]))
        
        data  = list(model[rowline])
        self.selectedID  = int(data[1])  # @+
        self.selectedObj = str(data[2])
        self.selectedChn = str(data[3])
        res = self.VObj.chains[self.selectedChn].residues_by_index[self.selectedID]
        
        
        '''centering and selecting'''
        frame = self.vm_session.get_frame ()
        res.get_center_of_mass(frame = frame)
        self.vm_session.glwidget.vm_widget.center_on_coordinates(res.Vobject, res.mass_center)
        
        self.vm_session._selection_function (res.atoms[0], _type = 'residue')
        self.vm_session.glwidget.queue_draw()
        
        self.atom_liststore.clear()
        for atom in res.atoms:
            self.atom_liststore.append(list([True, int(atom.index), atom.name, atom.symbol, atom.charge, atom.atom_id ]))

        #self.treeview_atom.set_model(self.atom_liststore)
  
    
    def on_treeview_Objects_button_release_event(self, tree, event):
        #print ( tree, event)
        
        if event.button == 3:
            dprint (3)
            #selection     = tree.get_selection()
            #model         = tree.get_model()
            #(model, iter) = selection.get_selected()
            #
            #
            #
            #if iter != None:
            #    self.selectedID  = int(model.get_value(iter, 1))  # @+
            #    self.selectedObj = str(model.get_value(iter, 2))
            #    print(self.selectedID, self.selectedObj, self.VObj.residues[self.selectedID])
            #    #self.builder.get_object('TreeViewObjLabel').set_label('- ' +self.selectedObj+' -' )
            #
            #    #widget = self.builder.get_object('treeview_menu')
            #    #widget.popup(None, None, None, None, event.button, event.time)
            #    #print ('button == 3')


        if event.button == 2:
            #print ('button == 2')
            self.treeview.get_selection().set_mode(Gtk.SelectionMode.SINGLE)

            
            selection     = tree.get_selection()
            model         = tree.get_model()
            (model, iter) = selection.get_selected()

            if iter != None:
                self.selectedID  = int(model.get_value(iter, 1))  # @+
                self.selectedObj = str(model.get_value(iter, 2))
                self.selectedChn = str(model.get_value(iter, 3))
                res = self.VObj.chains[self.selectedChn].residues_by_index[self.selectedID]
                frame = self.vm_session.get_frame ()
                res.get_center_of_mass(frame = frame)
                
                self.vm_session.glwidget.vm_widget.center_on_coordinates(res.Vobject, res.mass_center)
        
                self.atom_liststore.clear()
                for atom in res.atoms:
                     self.atom_liststore.append(list([True, int(atom.index), atom.name, atom.symbol, atom.charge, atom.atom_id ]))
            
            
            
            self.treeview.get_selection().set_mode(Gtk.SelectionMode.MULTIPLE)

        
        if event.button == 1:
            self.treeview.get_selection().set_mode(Gtk.SelectionMode.SINGLE)

            selection     = tree.get_selection()
            model         = tree.get_model()
            (model, iter) = selection.get_selected()
            
            
            if iter != None:
                self.selectedID  = int(model.get_value(iter, 1))  # @+
                self.selectedObj = str(model.get_value(iter, 2))
                self.selectedChn = str(model.get_value(iter, 3))
                res = self.VObj.chains[self.selectedChn].residues_by_index[self.selectedID]
                
                
                self.atom_liststore.clear()
                #self.atom_liststore = Gtk.ListStore(bool, int, str, str, float)
                for atom in res.atoms:
                     self.atom_liststore.append(list([True, int(atom.index), atom.name, atom.symbol, atom.charge, int(atom.atom_id)]))
                
                #self.treeview_atom.set_model(self.atom_liststore)
            self.treeview.get_selection().set_mode(Gtk.SelectionMode.MULTIPLE)

        
      
    
    def on_chk_renderer_toggled(self, cell, path, model):
        dprint(model[path][0])

   
   
    def residue_filter_func(self, model, iter, data):
        """Tests if the language in the row is the one in the filter"""
        if (
            self.current_filter_residue is None
            or self.current_filter_residue == "all"
        ):
            return True
        else:
            return model[iter][2] == self.current_filter_residue
   
            
    def chain_filter_func(self, model, iter, data):
        """Tests if the language in the row is the one in the filter"""
        if (
            self.current_filter_chain is None
            or self.current_filter_chain == "all"
        ):
            return True
        else:
            return model[iter][3] == self.current_filter_chain

    #def on_selection_button_clicked(self, widget):
    #    """Called on any of the button clicks"""
    #    # we set the current language filter to the button's label
    #    self.current_filter_language = widget.get_label()
    #    print("%s language selected!" % self.current_filter_language)
    #    # we update the filter, which updates in turn the view
    #    self.language_filter.refilter()
    #
    #def on_button1_clicked(self, widget):
    #    print("Hello")
    #
    #def on_button2_clicked(self, widget):
    #    print("Goodbye")
    def update (self):
        """ Function doc """
        dprint('VismolGoToAtomWindow2 update')
        pass
        #self.combobox_vobjects.set_active(-1)

class VismolSelectionTypeBox(Gtk.Box):
    """ Class doc """
    
    def __init__ (self, vismol_session=None):
        """ Class initialiser """
        Gtk.Box.__init__(self)
        #self.set_orientation(Gtk.Orientation.VERTICAL)
        #self.set_orientation(Gtk.Orientation.HORIZONTAL)
        #self.set_spacing(5)
        self.box           = Gtk.Box(orientation = Gtk.Orientation.HORIZONTAL, spacing = 6)
        self.vm_session = vismol_session
        #combobox
        self.combobox_selection_type = Gtk.ComboBox.new_with_model(self.vm_session.vm_selection_modes_list_store)
        self.combobox_selection_type.set_model(self.vm_session.vm_selection_modes_list_store)
        
        self.renderer_text = Gtk.CellRendererText()
        self.combobox_selection_type.pack_start(self.renderer_text, True)
        self.combobox_selection_type.add_attribute(self.renderer_text, "text", 0)
        
        self.combobox_selection_type.connect('changed', self.on_combobox_selection_type)
        #
        
        #labels        
        self.label_selecting_by = Gtk.Label('selecting by: ')
        
        #toggle_button
        self.toggle_button_selecting_mode = Gtk.ToggleButton('Viewing')
        self.toggle_button_selecting_mode.connect('clicked', self.on_toggle_button_selecting_mode)
        

        # Packing 
        self.box.pack_start(self.toggle_button_selecting_mode, False, False, 0)
        self.box.pack_start(self.label_selecting_by          , False, False, 0)
        self.box.pack_start(self.combobox_selection_type     , False, False, 0)
        self.combobox_selection_type.set_active(1)

    def on_combobox_selection_type (self, combobox):
        """ Function doc """
        self.active = combobox.get_active()
        
        if self.active == 0:
            self.vm_session.viewing_selection_mode(sel_type = 'atom')
        if self.active == 1:
            self.vm_session.viewing_selection_mode(sel_type = 'residue')
        if self.active == 2:
            self.vm_session.viewing_selection_mode(sel_type = 'chain')
        
        if self.active == 3:
            self.vm_session.viewing_selection_mode(sel_type = 'molecule')
            
    def change_sel_type_in_combobox (self, sel_type):
        """ Function doc """
        
        if sel_type == 'atom':
            self.combobox_selection_type.set_active(0)
        if sel_type == 'residue':
            self.combobox_selection_type.set_active(1)
        if sel_type == 'chain':
            self.combobox_selection_type.set_active(2)
        if sel_type == 'molecule':
            self.combobox_selection_type.set_active(3)
    
    def on_toggle_button_selecting_mode (self, button):
        """ Function doc """
        if button.get_active():
            state = "on"
            self.vm_session.picking_selection_mode = True
            button.set_label('Picking')
            dprint(self.combobox_selection_type.get_active())
            self.vm_session._selection_function (None)
            self.vm_session.glwidget.vm_widget.queue_draw()
            
            self.combobox_selection_type.set_sensitive(False)
            self.label_selecting_by.set_sensitive(False)
            
            circle = Gdk.Cursor(Gdk.CursorType.CROSSHAIR)
            button.get_window().set_cursor(circle)
            
        else:
            state = "off"
            self.vm_session.picking_selection_mode = False
            button.set_label('Viewing')
            self.vm_session.glwidget.vm_widget.queue_draw()
            
            self.combobox_selection_type.set_sensitive(True)
            self.label_selecting_by.set_sensitive(True)
            
            
            circle = Gdk.Cursor(Gdk.CursorType.ARROW)
            button.get_window().set_cursor(circle)
            
    def change_toggle_button_selecting_mode_status (self, status = False):
        """ Function doc """
        self.toggle_button_selecting_mode.set_active(status)

    
    
    def update (self):
        """ Function doc """
        dprint('VismolSelectionTypeBox update')

class VismolTrajectoryFrame(Gtk.Frame):
    """ Class doc """
    
    def __init__ (self, vismol_session = None):
        """ Class initialiser """
        self.vm_session = vismol_session 
        
        self.frame      =Gtk.Frame()
        #self.frame.set_shadow_type(Gtk.SHADOW_IN)
        self.frame.set_border_width(4)
        
        self.box        = Gtk.Box(orientation = Gtk.Orientation.VERTICAL, spacing = 6) 
        
        self.box.set_margin_top    (3)
        self.box.set_margin_bottom (3)
        self.box.set_margin_left   (3)
        self.box.set_margin_right  (3)
        
        self.value      = 1
        self.scale      = Gtk.Scale()
        
        #self.adjustment = Gtk.Adjustment(self.value, 1, 1, 0, 1, 0)
        self.adjustment     = Gtk.Adjustment(value         = self.value,
                                             lower         = 0,
                                             upper         = 100,
                                             step_increment= 1,
                                             page_increment= 1,
                                             page_size     = 1)
        
        
        
        self.scale.set_adjustment ( self.adjustment)
        self.scale.set_digits(0)
        self.scale.connect("change_value", self.on_scaler_frame_change_change_value)
        self.box.pack_start(self.scale, True, True, 0)
        
        self.vbox =  Gtk.Box(orientation = Gtk.Orientation.HORIZONTAL, spacing = 6)
        
        self.functions = [self.reverse,None, None, self.forward ]
        c = 0
        for label in ['<<','#', '>','>>']:
            button = Gtk.Button(label)
            self.vbox.pack_start(button, True, True, 0)
            
            if self.functions[c]:
                button.connect("clicked", self.functions[c])
            c += 1 
            
        self.box.pack_start(self.vbox, True, True, 0)
        
        
        #----------------------------------------------------------------------------
        self.label2 =  Gtk.Label('Object:')
        self.combobox_vobjects = Gtk.ComboBox.new_with_model(self.vm_session.vm_objects_list_store)
        self.combobox_vobjects.connect("changed", self.on_combobox_vobjects_changed)
        self.renderer_text = Gtk.CellRendererText()
        self.combobox_vobjects.pack_start(self.renderer_text, True)
        self.combobox_vobjects.add_attribute(self.renderer_text, "text", 1)
        
        #self.box.pack_start(self.combobox, True, True, 0)
        #----------------------------------------------------------------------------
        
        #----------------------------------------------------------------------------
        self.vbox2 = Gtk.Box(orientation = Gtk.Orientation.HORIZONTAL, spacing = 6)
        self.label = Gtk.Label('FPS:')
        
        self.entry = Gtk.Entry()
        self.entry.set_text(str(25))
        
        self.fps_adjustment = Gtk.Adjustment(value          = 24 , 
                                             upper          = 100, 
                                             step_increment = 1  , 
                                             page_increment = 10 )
        #self.fps_adjustment = Gtk.Adjustment(value         = float, 
        #                                     lower         = float, 
        #                                     upper         = float, 
        #                                     step_increment= float, 
        #                                     page_increment= float, 
        #                                     page_size     = float)
        self.entry = Gtk.SpinButton()
        self.entry.set_adjustment ( self.fps_adjustment)
        
        self.vbox2.pack_start(self.label2, True, True, 0)
        self.vbox2.pack_start(self.combobox_vobjects, True, True, 0)
        self.vbox2.pack_start(self.label, True, True, 0)
        self.vbox2.pack_start(self.entry, True, True, 0)
        self.box.pack_start(self.vbox2, True, True, 0)
        #----------------------------------------------------------------------------


    def forward (self, button):
        """ Function doc """
        value =  int(self.scale.get_value())
        value = value+1
        self.scale.set_value(int(value))
        self.vm_session.set_frame(int(value))
        dprint(value)

    def reverse (self, button):
        """ Function doc """
        value = int(self.scale.get_value())
        
        if value == 0:
            pass
        else:
           value = value-1
        
        self.vm_session.set_frame(int(value))
        self.scale.set_value(value)
        dprint(value)
    
    def get_box (self):
        """ Function doc """
        #self.add(self.box)
        return self.box
        #return self.frame
        
    def on_combobox_vobjects_changed (self, widget):
        """ Function doc """
        dprint('\n\n',widget)
        dprint('\n\n',widget.get_active())
        
        cb_index = widget.get_active()
        if cb_index in self.vm_session.vm_objects_dic:
            self.VObj = self.vm_session.vm_objects_dic[widget.get_active()]
            #self.VObj = self.vismol_session.vismol_objects[widget.get_active()]
            number_of_frames = len(self.VObj.frames)
            self.scale.set_range(0, int(number_of_frames))
            self.scale.set_value(self.vm_session.get_frame())
        else:
            pass

    def on_scaler_frame_change_change_value (self, hscale, text= None,  data=None):
        """ Function doc """
        value = hscale.get_value()
        self.vm_session.set_frame(int(value)) 
        #print(value)

    def update (self):
        """ Function doc """
        dprint('VismolTrajectoryFrame update')
        #for index , visObj in self.vismol_session.vm_objects_dic.items():
        #last_obj = len(self.vismol_session.vismol_objects) -1
        last_obj = len(self.vm_session.vm_objects_dic.items()) -1
        self.combobox_vobjects.set_active(last_obj)


#VismolTrajectoryFrame()
#frame = VismolTrajectoryFrame()
#win = Gtk.Window()
#win.add(frame.get_box())
#win.show_all()        
#Gtk.main()



