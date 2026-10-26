#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
#  gtk3.py
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


#import VISMOL.glCore.shapes as shapes
#import VISMOL.glCore.glaxis as glaxis
#import VISMOL.glCore.glcamera as cam
#import VISMOL.glCore.operations as op
#import VISMOL.glCore.sphere_data as sph_d
#import VISMOL.glCore.vismol_shaders as vm_shader
#import VISMOL.glCore.matrix_operations as mop

from vismol.utils.debug import dprint
import gi, sys
gi.require_version('Gtk', '3.0')
from gi.repository import Gtk, Gdk
import time

class PlayerFrame:
    """ Class doc """
    
    def __init__ (self, vismol_session):
        """ Class initialiser """
        pass
        xml ='''
<?xml version="1.0" encoding="UTF-8"?>
<!-- Generated with glade 3.22.2 -->
<interface>
  <requires lib="gtk+" version="3.20"/>
  <object class="GtkWindow" id="main_window">
    <property name="can_focus">False</property>
    <child type="titlebar">
      <placeholder/>
    </child>
    <child>
      <object class="GtkBox" id="main_frame">
        <property name="visible">True</property>
        <property name="can_focus">False</property>
        <property name="orientation">vertical</property>
        <property name="spacing">7</property>
        <child>
          <object class="GtkScale" id="scaler_frame_change">
            <property name="visible">True</property>
            <property name="can_focus">True</property>
            <property name="round_digits">2</property>
          </object>
          <packing>
            <property name="expand">False</property>
            <property name="fill">True</property>
            <property name="position">0</property>
          </packing>
        </child>
        <child>
          <object class="GtkBox">
            <property name="visible">True</property>
            <property name="can_focus">False</property>
            <property name="spacing">7</property>
            <child>
              <placeholder/>
            </child>
            <child>
              <object class="GtkButton">
                <property name="visible">True</property>
                <property name="can_focus">True</property>
                <property name="receives_default">True</property>
                <child>
                  <object class="GtkImage">
                    <property name="visible">True</property>
                    <property name="can_focus">False</property>
                    <property name="pixbuf">../old/icon_player_backward.png</property>
                  </object>
                </child>
              </object>
              <packing>
                <property name="expand">False</property>
                <property name="fill">True</property>
                <property name="position">1</property>
              </packing>
            </child>
            <child>
              <object class="GtkButton" id="button_stop">
                <property name="visible">True</property>
                <property name="can_focus">True</property>
                <property name="receives_default">True</property>
                <child>
                  <object class="GtkImage">
                    <property name="visible">True</property>
                    <property name="can_focus">False</property>
                    <property name="pixbuf">../old/icon_player_stop.png</property>
                  </object>
                </child>
              </object>
              <packing>
                <property name="expand">False</property>
                <property name="fill">True</property>
                <property name="position">2</property>
              </packing>
            </child>
            <child>
              <object class="GtkButton" id="button_play">
                <property name="visible">True</property>
                <property name="can_focus">True</property>
                <property name="receives_default">True</property>
                <child>
                  <object class="GtkImage">
                    <property name="visible">True</property>
                    <property name="can_focus">False</property>
                    <property name="pixbuf">../old/icon_player_play.png</property>
                  </object>
                </child>
              </object>
              <packing>
                <property name="expand">False</property>
                <property name="fill">True</property>
                <property name="position">3</property>
              </packing>
            </child>
            <child>
              <object class="GtkButton">
                <property name="visible">True</property>
                <property name="can_focus">True</property>
                <property name="receives_default">True</property>
                <child>
                  <object class="GtkImage">
                    <property name="visible">True</property>
                    <property name="can_focus">False</property>
                    <property name="pixbuf">../old/icon_player_forward.png</property>
                  </object>
                </child>
              </object>
              <packing>
                <property name="expand">False</property>
                <property name="fill">True</property>
                <property name="position">4</property>
              </packing>
            </child>
            <child>
              <placeholder/>
            </child>
          </object>
          <packing>
            <property name="expand">False</property>
            <property name="fill">True</property>
            <property name="position">1</property>
          </packing>
        </child>
        <child>
          <placeholder/>
        </child>
        <child>
          <placeholder/>
        </child>
      </object>
    </child>
  </object>
</interface>


        '''
        
        self.vismol_session = vismol_session
        
        self.builder = Gtk.Builder()
        #self.builder.add_from_string(xml)
        self.builder.add_from_file('GTK3VisMol/GTKGUI/gtkWidgets/VisMolPlayer.glade')
        self.builder.connect_signals(self)
        self.value = 0
        self.main  = self.builder.get_object('main_window')
        self.main_frame = self.builder.get_object('main_frame')
        self.scale = self.builder.get_object('scaler_frame_change')
        self.button_play = self.builder.get_object('button_play')
        self.scale.set_digits(0)
        
        self.adjustment = Gtk.Adjustment(self.value, 0, 100, 0, 10, 0)
        
        self.scale.set_adjustment ( self.adjustment)
        
        #for 
        
        
        # combo box  - using the same Vismol_Objects_ListStore
        self.combo_box_objs = self.builder.get_object('combo_box_objects')
        self.combo_box_objs.set_model(self.vismol_session.Vismol_Objects_ListStore)
        renderer_text = Gtk.CellRendererText()
        self.combo_box_objs.pack_start(renderer_text, True)
        self.combo_box_objs.add_attribute(renderer_text, "text", 2)
        
        #print ('print', self.main)
    
    def show_player_main_window (self):
        """ Function doc """
        self.main.show()
        
    
        #Gtk.main()
    def on_combobox_change (self, combo):
        """ Function doc """
        #print (combobox)
        
        tree_iter = combo.get_active_iter()
        
        if tree_iter is not None:
            model = combo.get_model()
            row_id, name = model[tree_iter][:2]
            dprint(model[tree_iter][:])
            number_of_frames = model[tree_iter][4]
            self.scale.set_range(0, int(number_of_frames)-1)
            #print("Selected: ID=%d, name=%s" % (row_id, name))
        else:
            entry = combo.get_child()
            #print("Entered: %s" % entry.get_text())




    def forward (self, button):
        """ Function doc """
        value =  int(self.scale.get_value())
        value = value+1
        self.scale.set_value(int(value))
        self.vismol_session.set_frame(int(value))
        dprint(value)

    def reverse (self, button):
        """ Function doc """
        value = int(self.scale.get_value())
        
        if value == 0:
            pass
        else:
           value = value-1

        self.vismol_session.set_frame(int(value))
        self.scale.set_value(value)
        dprint(value)

    def on_button_print (self, button):
        """ Function doc """
        #self.adjustment = Gtk.Adjustment(self.value+1, 0, 90, 1, 10, 0)
        value = self.scale.get_value()
        
        while value != 100:
            value += 1
            self.vismol_session.set_frame(int(value))
            
            self.scale.set_value(int(value))
            time.sleep(0.1)
            
            dprint(value)
        #self.scale.set_range(10  , 90)
        #self.scale.set_digits(0)

    def on_scaler_frame_change_change_value (self, hscale, text= None,  data=None):
        """ Function doc """
        value = hscale.get_value()
        #cmd.frame( int (valor) )
        #BondTable = self.project.BondTable
        
        
        #MAX  = int(self.builder.get_object('trajectory_max_entrey').get_text())
        #MIN  = int(self.builder.get_object('trajectory_min_entrey').get_text())

        #scale = self.builder.get_object("trajectory_hscale")
        #scale.set_increments(1, 10)
        #scale.set_digits(0)
        self.vismol_session.set_frame(int(value)) 
        #print(value)
        
    def show_player (self, show = True):
        """ Function doc """
        self.main.show()
        Gtk.main()
        
        
#PlayerFrame = PlayerFrame()
