#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
#  simple.py
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

import logging
import gi, sys, os
gi.require_version("Gtk", "3.0")
from gi.repository import Gtk, Gdk

#               Installation is not necessary anymore.
#.This retrieves the absolute path of the script file that is currently being executed.
main_file = os.path.abspath(__file__)

#.This extracts the directory (folder) from the absolute path obtained in the previous step.
HOME = os.path.dirname(main_file)
HOME = os.path.split(HOME)
HOME = HOME[0]
#.Adding GRAPHIC ENGINE LIB
sys.path.append(os.path.join(HOME,"src/graphics_engine/src"))
sys.path.append(os.path.join(HOME,"src/"))



from vismol.core.vismol_session import VismolSession

logger = logging.getLogger(__name__)





def main():
    logging.basicConfig(format="%(asctime)s,%(msecs)d %(levelname)-8s [%(filename)s:%(lineno)d] %(message)s",
                        datefmt="%Y-%m-%d:%H:%M:%S", level=logging.DEBUG)
    # logging.basicConfig(level=logging.DEBUG)
    vm_session = VismolSession(toolkit="Gtk_3.0")
    vm_session.vm_widget.insert_glmenu()
    window = Gtk.Window(title="Vismol window")
    container = Gtk.Box(orientation=Gtk.Orientation.VERTICAL)
    container.pack_start(vm_session.vm_widget, True, True, 0)
    window.connect("key-press-event", vm_session.vm_widget.key_pressed)
    window.connect("key-release-event", vm_session.vm_widget.key_released)
    window.add(container)
    window.connect("delete-event", Gtk.main_quit)
    window.show_all()
    try:
        filein = sys.argv[-1]
        vm_session.load_molecule(filein)
    except:
        pass
    Gtk.main()
    return 0

if __name__ == "__main__":
    main()

