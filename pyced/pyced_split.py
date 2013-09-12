'''
Python module 'pyced', for processing events for the CED.
Version 1.0
Created on Sep 04, 2013
Authors: Frank Gaede (DESY), Bjoern Klaas (Uni Goettingen)

Uses ctypes to call the C Event Display (CED) and the 
PyLCIO bindings to display various collections.

Make sure PyLCIO and PyROOT are in your Python path, eg:

    export PYTHONPATH=${LCIO}/src/python:${ROOTSYS}/lib

and that libCED and libMarlinUtil are in your dynamica linker path,
    export LD_LIBRARY_PATH=$ILCSOFT/CED/HEAD/lib:$ILCSOFT/MarlinUtil/HEAD/lib:$LD_LIBRARY_PATH

and that glced is running.

Note: CED (incl. glced) can be downloaded from:
    svn co https://svnsrv.desy.de/public/marlinreco/CED/trunk

'''

# --- Python dependencies ---
import os, sys, select, platform
from math import sin, cos, fabs, sqrt, atan2
from ctypes import *
import time

# --- LCIO dependencies ---
from pyLCIO.io.LcioReader import LcioReader
from pyLCIO import UTIL, EVENT #, IMPL

# --- Local/Custom dependencies ---
from printPDG import printPdg

####################################################################

# Definition of class PYCED
execfile( "pycedClass.py" )
# For loading libraries, wrappers
execfile( "pycedSupport.py" )
# Drawing functions for individual types (Hits, Tracks, MCParticles...)
execfile( "pycedDraw.py" )
# readEvent, drawEvent, drawEvents and idleLoop
execfile( "pycedEvent.py" )
# Functions for userfriendly interactive mode
execfile( "pycedInteractive.py" )

####################################################################

# Load config file
execfile( "pyced.cfg.py" )

# Build parameter maps
buildDm()
buildTm()
buildNm()

####################################################################

if __name__ == "__main__":

    if len( sys.argv ) < 2:
        print 'Usage:\n\tpython %s <fileName>' % ( os.path.split( sys.argv[0] )[1] )
        sys.exit( 0 )
        
#   Read the file name from the command line input
    g.fileName = sys.argv[1]


    init()
    drawEvents()

####################################################################
