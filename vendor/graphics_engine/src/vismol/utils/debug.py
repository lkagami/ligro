#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
#  VisMol / graphics_engine - Debug/verbose-output toggle
#
#  Description:
#      Centralises the ad-hoc "print(...)" debug statements scattered
#      across this package behind a single ON/OFF switch, so that by
#      default it runs quietly, but the exact same messages can be brought
#      back at any time for troubleshooting.
#
#      Shares the same "EASYHYBRID_DEBUG" environment variable used by
#      EasyHybrid3's own util.debug module, so a single switch controls
#      both. Enable by setting it before starting the application, e.g.:
#
#          EASYHYBRID_DEBUG=1 python easyhybrid.py
#
#      or at runtime:
#
#          from vismol.utils.debug import set_debug
#          set_debug(True)
#
import os

_TRUE_VALUES = ( '1', 'true', 'yes', 'on' )

DEBUG = os.environ.get ( 'EASYHYBRID_DEBUG', '0' ).strip ( ).lower ( ) in _TRUE_VALUES


def set_debug ( value ):
    """ Turns debug output on/off at runtime. """
    global DEBUG
    DEBUG = bool ( value )


def is_debug ( ):
    """ Returns whether debug output is currently enabled. """
    return DEBUG


def dprint ( *args, **kwargs ):
    """
    Drop-in replacement for the builtin print(), silenced unless debug
    output is enabled (see set_debug() / the EASYHYBRID_DEBUG environment
    variable).
    """
    if DEBUG:
        print ( *args, **kwargs )
