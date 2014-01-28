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

from ROOT import TVector3
####################################################################
####################################################################

# Definition of class PYCED
class PYCED :
    """Main class for pyced holding constants and color codes.

    The variables set here are to be left unchanged. For modifications turn to
    the separate 'pyced.cfg.py'. The class is dynamically expanded with
    attributes for customizing the CED.
    """

    def __init__(self):
        """Sets constants and color codes. Leave unchanged."""
#       some enum like constants 
#       --- for drawHelix:
        self.Default=0
        self.AtIP=1
        self.AtFirstHit=2
        self.AtLastHit=3
        self.AtCalo=4
#       --- for color schemes
        self.nCol       = 20 
        self.nScheme    = 10 
        self.Red         =   0 
        self.Orange      =   1 
        self.Plum        =   2 
        self.Violet      =   3 
        self.Blue        =   4 
        self.LightBlue   =   5 
        self.Aquamarine  =   6 
        self.Green       =   7 
        self.Olive       =   8 
        self.Yellow      =   9 
        self.Dark        =   10 
        self.Light       =   11 
        self.Classic     =   12 
#       List containing colorcodes
        self.cols = [  

            0x510000,0x660202,0x720202,0x840202,0x960303,0xa80303,0xbc0303,
            0xce0404,0xe00404,0xef0404,0xff0000,0xfc0505,0xf91111,0xf92222,
            0xfc4949,0xf97777,0xf98e8e,0xf9aeae,0xf7b2b2,0xf7cdcd,
            
            0x512401,0x662d03,0x7a3602,0x934204,0xa54a04,0xb75103,0xc65803,
            0xd35e04,0xe56604,0xff6e00,0xf4710c,0xf4791a,0xf9842a,0xf98f3e,
            0xf99d57,0xf9a768,0xf9b47f,0xf9bf93,0xf9c9a4,0xf9d2b3,
            
            0x48014c,0x610266,0x700277,0x93039b,0xb103ba,0xc904d3,0xda04e5,
            0xe604f2,0xfa00ff,0xf00ffc,0xec1df7,0xef2ff9,0xeb42f4,0xec58f4,
            0xed6bf4,0xf486f9,0xf59af9,0xf8b5fc,0xf4c3f7,0xf8d6f9,
            
            0x2d0251,0x3d026d,0x4a0284,0x550399,0x5f03aa,0x6903bc,0x7102cc,
            0x8004e5,0x9800ff,0x8e0ef7,0x9922f9,0xa134f9,0xa845f9,0xb057f9,
            0xbc70f9,0xbf77f9,0xc98ef9,0xd3a4f9,0xddbbf9,0xecd6ff,
            
            0x00004f,0x020268,0x03037f,0x030399,0x0303b2,0x0404cc,0x0404e0,
            0x0404ef,0x0000ff,0x0c0cf9,0x1b1bf9,0x2a2af9,0x3939f9,0x4d4df9,
            0x6363f9,0x7272f9,0x8484f9,0x9898f9,0xb3b3f9,0xcacaf9,
            
            0x01434c,0x025c68,0x027382,0x04899b,0x0397aa,0x03abc1,0x04b9d1,
            0x04cbe5,0x00d8ff,0x0cdef9,0x1bddf7,0x2ae1f9,0x3ddff4,0x59e4f7,
            0x6be9f9,0x7cebf9,0x94f0fc,0xa3edf7,0xb6f2f9,0xc7f4f9,
            
            0x014443,0x01605f,0x027f7d,0x039996,0x03b2af,0x01c6c3,0x04ddda,
            0x04efeb,0x00ffff,0x09f9f5,0x0ef9f5,0x20f9f6,0x32fcf9,0x3ef9f6,
            0x52f9f7,0x6bf9f7,0x7ff9f7,0x95f9f8,0xb1f9f8,0xcaf9f9,
            
            0x016001,0x027002,0x027f02,0x029102,0x05aa05,0x05bf05,0x06d306,
            0x04e504,0x00ff00,0x09f909,0x18f918,0x2cf92c,0x43f943,0x52f952,
            0x63f963,0x77f977,0x8bf98b,0x9ff99f,0xb3f9b3,0xcaf9ca,
            
            0x344701,0x4b6603,0x608202,0x739b04,0x83b203,0x96cc04,0xa7e204,
            0xb1ef04,0xb6ff00,0xbaf713,0xbff725,0xc5f73b,0xcbf751,0xd3f968,
            0xd7f97a,0xd8f48b,0xe2f9a2,0xe1f4ad,0xe7f7bb,0xe9f4c8,
            
            0x565501,0x727002,0x898702,0xa5a303,0xb7b403,0xd1cd04,0xe2df04,
            0xefeb04,0xffff00,0xf9f509,0xf9f618,0xf9f62a,0xf7f43b,0xf9f64a,
            0xf9f759,0xf9f76b,0xf9f77c,0xf9f88b,0xf9f89f,0xfcfbbd ]
        
#   cs == ColorScheme
    def  setColorScheme( self, cs ):
        """Fills 'colors' with the color codes from the selected scheme"""

        self.colors = [] 

#       If a specific color was chosen load the corresponding block of colors
        if( cs >= self.Red and cs <= self.Yellow ):
            for i in xrange( self.nCol ):
                self.colors.append(  self.cols[ (cs*self.nCol) + i ]  ) 

#       If scheme dark or light are chosen load one color from each colorblock
        elif( cs == self.Dark):
            for i in xrange( self.nCol ):
                for j in xrange( self.nScheme ):
                    self.colors.append(  self.cols[ ( j*self.nCol) + i ]  ) 

        elif( cs == self.Light):
            for i in xrange( self.nCol-1, 0 , -1 ):
                for j in xrange( self.nScheme ):
                    self.colors.append(  self.cols[ ( j*self.nCol) + i ]  ) 

#       If no scheme was chosen load all colors
        else:
            self.colors = self.cols
        self.lenColors = len(self.colors)
    
    ### END OF setColorSchme #######################################

#--- dymmy function for userDraw method (overwrite in pyced_user.py)
    def userDraw(self):
        nothing=42

### END OF CLASS 'PYCED' ###########################################

# Create a global instance of the class. Will be extended to 
# contain all settings.
g = PYCED()

####################################################################
####################################################################

# For small or supporting tasks

def fixedColor( colorCode ) :
    """Sets a fixed color using a callable"""

    def getFixedColor( hit ) :
        return colorCode
    return getFixedColor

def SimTrackerHitColor(hit) :
    """Sets color of a SimTrackerHit depending on hitID"""

    if hit.getMCParticle() :
        return g.colors[ hit.getMCParticle().id() % g.lenColors ]
    else :
        return g.colors[ 0 ]

def SimCalorimeterHitColor(hit) :
    """Sets color of a SimCalorimeterHit depending on hitID"""

    if hit.getNMCContributions() != 0 and hit.getParticleCont(0) :
        return g.colors[ hit.getParticleCont(0).id() % g.lenColors ]
    else :
        return g.colors[ 0 ]

### END OF ColorSet Functions ######################################

def clearLayerDescription() :
    """Clears the dict holding the layer descriptions"""

    g.layerDescriptions = {}


def applyLayerDescription() :
    """Sends the layer descriptions to the CED"""

    for layer in g.layerDescriptions :
        describe_layer( g.layerDescriptions[layer] , layer)


def addLayerDescription( description, layer ) :
    """Adds a layer description to the dict"""

    if g.layerDescriptions.get(layer) == None :
        g.layerDescriptions[layer] = description
    else :
        g.layerDescriptions[layer] = ( g.layerDescriptions[layer] + ", " +
                                       description )

### END OF LayerDescription Functions ##############################

def getHitsFromTracks(track) :
    """Collects hits from subtracks of tracks"""

    hits = []
    if not len( track.getTracks() ) :

#       extend hits with all elements from track
        for hit in track.getTrackerHits() :
            hits.append( hit )
        for subtrack in track.getTracks() :
            for hit in subtrack.getTrackerHits() :
                hits.append(hit)

    if not len(hits) : 
        hits = track.getTrackerHits() 
    return hits

### END OF getHitsFromTracks #######################################

def addDrawingParameters( col ) :
    """Attaches the drawing parameters from the maps to the collections"""

    options = {}
    options.update( g.dm )
    if col.type in g.tm :
        options.update( g.tm[col.type] )
    if col.name in g.nm :
        options.update( g.nm[col.name] )

    for option in options :
        setattr( col, option, options[ option ] )

### END OF addDrawingParameters ####################################

def init():
    """Loads the required libraries and connects to glced."""

    g.ced = load_library( "libCED" ) 
    g.marlinutil = load_library( "libMarlinUtil" )
    g.reader = LcioReader( g.fileName )

    g.ced.ced_register_elements()
    g.ced.ced_client_init( g.host , g.port ) 
    g.setColorScheme( g.colorScheme ) 

    g.objects = dict() 

### END OF init ####################################################

def eCut( cut, type = "greater" ) :
    """Performs an energy cut.

    Arguments:
    cut specifies the energy threshold,
    type specifies if energies "greater", "smaller", or between cut
        and type should be displayed.
    """

    def make_ecut( obj ) :
        try :
            if type == "greater" :
                return obj.getEnergy() > cut
            elif type == "smaller" :
                return obj.getEnergy() < cut
            else :
                return obj.getEnergy() > cut and obj.getEnergy() < type
        except AttributeError :
#           True or False, draw if no getEnergy or not?
            return True
    return make_ecut

def ptCut( cut, type = "greater" ) :
    """Performs a p_t cut.

    Arguments:
    cut specifies the p_t threshold,
    type specifies if p_t "greater", "smaller", or between cut
        and type should be displayed.
    """

    def make_ptcut( obj ) :
#       TODO: Is better with px, py or getMom = obj.getMomentum?
        try :
            pt = sqrt( obj.getMomentum()[0]*obj.getMomentum()[0] + 
                       obj.getMomentum()[1]*obj.getMomentum()[1]
                       )
        except AttributeError :
            return True

        if type == "greater" :
            return pt > cut
        elif type == "smaller" :
            return pt < cut
        else :
            return pt > cut and pt < type
    return make_ptcut

def cosThetaCut( cut, type = "greater" ) :
    """Performs a cos Theta cut.

    Arguments:
    cut specifies the cos Theta threshold,
    type specifies if cos Theta "greater", "smaller", or between cut
        and type should be displayed.
    """

    def make_costhetacut( obj ) :
#       TODO: Test and optimize
        try : 
            p = sqrt( obj.getMomentum()[0]*obj.getMomentum()[0] + 
                      obj.getMomentum()[1]*obj.getMomentum()[1] +
                      obj.getMomentum()[2]*obj.getMomentum()[2]
                      )
            cosTheta = obj.getMomentum()[2] / p
        except AttributeError :
            return True

        if type == "greater" :
            return cosTheta > cut
        elif type == "smaller" :
            return cosTheta < cut
        else :
            return cosTheta > cut and cosTheta < type
    return make_costhetacut

### END OF cuts ####################################################

def load_library( fname ):
    """Loads the specified C library"""

    if (sys.platform == 'darwin') :
        fname = fname + '.dylib'
    else:
        fname = fname + '.so'

# if not( os.path.isfile(fname) ) :
#     print "ERROR: cannot load CED library: ",  fname 
#     print "       please set LD_LIBRARY_PATH properly ! " 
#     exit( 1 )

    return CDLL( fname )

### END OF load_library ############################################

# --- Wrappers --------------------------
def pced_hit_ID( x, y, z, type, layer, size, color, lcioID):
    """Wraps the CED function ced_hit_ID using ctypes"""

    g.ced.ced_hit_ID( c_float(x), c_float(y), c_float(z), type, layer, size, color, lcioID )

def pced_line_ID( x, y, z, px, py, pz, layer, size, color, mcpid ) :
    """Wraps the CED function ced_line_ID using ctypes"""

    g.ced.ced_line_ID( c_float(x), c_float(y), c_float(z), 
        c_float(px), c_float(py), c_float(pz), layer, size, color, mcpid )

def drawHelix( bField, charge, xs, ys, zs, px, py, pz, marker, layer , size ,
        color, helix_minr, _helix_max_r, _helix_max_z, trkid ):
    """Wraps the MarlinUtil function draw_helix using ctypes"""

    ml = marker | ( layer << 8 )
    g.marlinutil.draw_helix( c_float(bField), c_float(charge), c_float(xs),
            c_float(ys), c_float(zs), c_float(px), c_float(py), c_float(pz),
            ml, size, color, c_float(helix_minr),
            c_float(_helix_max_r), c_float(_helix_max_z), trkid )

def describe_layer( description , layerNumber) :
    """Wraps the CED function ced_describe_layer using ctypes"""

    g.ced.ced_describe_layer( c_char_p( description) , layerNumber )

### END OF Wrappers ################################################

####################################################################
####################################################################

# Drawing functions for individual types (Hits, Tracks, MCParticles...)

def drawHits( col ):
    """Draws collections containing only hits.

    Argument:
    Collection of type '(Sim)Tracker-/CalorimeterHit' and similar.
    """
#   Setting this increases speed (avoid dots in Python).
    cut = col.cut
    marker = col.marker
    layer = col.layer
    size = col.size
    color = col.color
    objects = g.objects
    
    for hit in col:
        if cut( hit ) :
#           Long notation for increased speed (Root vectors are
#           very expensive).
            x = hit.getPosition()[0]
            y = hit.getPosition()[1]
            z = hit.getPosition()[2]

            pced_hit_ID( x, y, z, marker, layer, 
                         size, color(hit), hit.id() )

#           Stores individual number for picking
            objects[ hit.id() ] = hit

### END OF drawHits ################################################


def drawTracks( col ):
    """Draws collections of type Track.

    Argument:
    Collection of type 'Track'.
    """

#   For increased speed (avoid dots in Python).
    cut = col.cut
    colors = g.colors
    lenColors = g.lenColors
    marker = col.marker
    layer = col.layer
    size = col.size
    hitsize = col.hitsize
    drawHelixForTracks = g.drawHelixForTracks
    bField = g.BField
    helix_max_r = g.helix_max_r
    helix_max_z = g.helix_max_z

    for i, track in enumerate(col) :
        if cut( track ) :
#           collect hits from all track segments
            hits = getHitsFromTracks(track)
            color = colors[ i % lenColors ]

            for hit in hits :
#               Long notation for increased speed (Root 
#               vectors are very expensive).
                x = hit.getPosition()[0]
                y = hit.getPosition()[1]
                z = hit.getPosition()[2]
                pced_hit_ID( x, y, z, marker, layer, hitsize,
                             color, track.id() )

#           Stores individual number for picking
            g.objects[ track.id() ] = track

            ts = 0
#           Case 0 is the default value which takes the whatever first track state
#           (e.g. InteractionPoint for ILD, or the simple track for test beam data)
            if drawHelixForTracks == g.Default :
                if track.getTrackStates().size() > 0 :
                    ts = track.getTrackStates().at(0)

#           The rest states are pre-defined
            elif drawHelixForTracks == g.AtIP :
                ts = track.getTrackState( EVENT.TrackState.AtIP )
            elif drawHelixForTracks == g.AtFirstHit :
                ts = track.getTrackState( EVENT.TrackState.AtFirstHit )
            elif drawHelixForTracks == g.AtLastHit :
                if track.getTracks().empty() :
                    ts = track.getTrackState( EVENT.TrackState.AtLastHit )
                else :
                    ts = track.getTracks().back().getTrackState( EVENT.TrackState.AtLastHit )
            elif drawHelixForTracks == g.AtCalo :
                if track.getTracks().empty() :
                    ts = track.getTrackState( EVENT.TrackState.AtCalorimeter ) 
                else :
                    ts = track.getTracks().back().getTrackState( EVENT.TrackState.AtCalorimeter )
                pced_hit_ID( ts.getReferencePoint()[0],
                             ts.getReferencePoint()[1],
                             ts.getReferencePoint()[2],
                             1 , layer , col.hitsize*10 , color, track.id() 
                             )

            if ts != 0 :
#               bField = Global::GEAR->getBField().at(  gear::Vector3D(0,0,0)  ).z() ; 

                if bField != 0.0 and fabs( ts.getOmega() ) > 0.00001 :
                    pt = bField * 3e-4 / fabs( ts.getOmega() )
                else :
                    pt = 1.e10

                if ts.getOmega() > 0 : charge = 1.
                else : charge = -1.

                px = pt * cos( ts.getPhi() )
                py = pt * sin( ts.getPhi() )
                pz = pt * ts.getTanLambda()

#               start point for drawing ( PCA to reference point )
                xs = ts.getReferencePoint()[0] -  ts.getD0() * sin( ts.getPhi() )
                ys = ts.getReferencePoint()[1] +  ts.getD0() * cos( ts.getPhi() )
                zs = ts.getReferencePoint()[2] +  ts.getZ0()

                if drawHelixForTracks >= 0 and pt > 0.01 :
#                   TODO: Remove hard-coded coloring?
                    drawHelix( bField, charge, xs, ys, zs, px, py, pz, 
                               marker, layer, size, 0xdddddd, 0.0,
                               helix_max_r, helix_max_z, track.id() 
                               )

### END OF drawTracks ##############################################

def drawMCParticles( col ) :
    """Draws collections of type MCParticle.

    Argument:
    Collection of type 'MCParticle'.
    """

#   For increased speed (avoid dots in Python).
    cut = col.cut
    marker = col.marker
    layer = col.layer
    size = col.size
    bField = g.BField
    helix_max_r = g.helix_max_r
    helix_max_z = g.helix_max_z
    usingParticleGun = g.usingParticleGun

    for mcp in col :
        if cut( mcp ) :
            charge = mcp.getCharge()

#           stable particles only   
            if mcp.getGeneratorStatus() != 1 and usingParticleGun == False : 
                continue

#           Long notation for increased speed (Root vectors are
#           very expensive).
            px = mcp.getMomentum()[0]
            py = mcp.getMomentum()[1]
            pz = mcp.getMomentum()[2]
            x = mcp.getVertex()[0]
            y = mcp.getVertex()[1]
            z = mcp.getVertex()[2]

#           Stores individual number for picking
            g.objects[ mcp.id() ] = mcp

#           Charged particles
            if fabs( charge ) > 0.0001 :
#                bField = Global::GEAR->getBField().at(  gear::Vector3D(0,0,0)  ).z() ; 

#               TODO: Remove hard-coded coloring?
                drawHelix( bField , charge, x, y, z, px, py, pz, 
                           marker, layer, size, 0x7af774, 0.0,
                           helix_max_r, helix_max_z, mcp.id() 
                           )

#           neutral particles
            else :
                if fabs( mcp.getPDG() ) == 22 :         # photon
                    color = 0xf9f920
                    r_max = g.ecalR
                    z_max = g.ecalZ
                elif (fabs( mcp.getPDG() ) == 12 or 
                      fabs( mcp.getPDG() ) == 14 or 
                      fabs( mcp.getPDG() ) == 16
                      ) :    # neutrino
                    color = 0xdddddd
                    r_max = g.hcalR * 2.
                    z_max = g.hcalZ * 2.
                else :                                  # neutral hadron
                    color = 0xb900de
                    r_max = g.hcalR
                    z_max = g.hcalZ

                pt = sqrt( px*px + py*py )
                p  = sqrt( pt*pt + pz*pz )

#               Hit in barrel or endcap
                if fabs( pt/pz ) > r_max/z_max :
                    length = r_max * p / pt
                else :
                    length = fabs( z_max * p / pz )

#               TODO: Remove hard-coded coloring?
                pced_line_ID( x, y, z, 
                              length*px/p, length*py/p, length*pz/p,
                              layer, size, color, mcp.id() )

### END OF drawMCParticles #########################################

def drawClusters( col ) :
    """Draws collections of type Cluster.

    Argument:
    Collection of type 'Cluster'.
    """

#   For increased speed (avoid dots in Python).
    cut = col.cut
    colors = g.colors
    lenColors = g.lenColors
    marker = col.marker
    layer = col.layer
    size = col.size
    hitsize = col.hitsize

#   TODO: Achieve using min() max()
    emin, emax = 1.e99, 0.
    for clu in col :
        e = clu.getEnergy()
        if   e > emax : emax = e
        elif e < emin : emin = e
      
    for i, clu in enumerate(col) :
        if cut( clu ) :
            hits = clu.getCalorimeterHits()

            color =  colors[ i % lenColors ]
            for hit in hits :
#               Long notation for increased speed (Root vectors
#               are very expensive).
                x = hit.getPosition()[0]
                y = hit.getPosition()[1]
                z = hit.getPosition()[2]
                pced_hit_ID( x, y, z, marker, layer, hitsize, color, clu.id() )
              
#           Long notation for increased speed (Root vectors
#           are very expensive).
            x = clu.getPosition()[0]
            y = clu.getPosition()[1] 
            z = clu.getPosition()[2]
            pced_hit_ID( x, y, z, marker, layer, hitsize*3, color, clu.id() ) ;
            
#           Implementation of functions from the LCVector3D (identical 
#           to Hep3Vector) class

#           Contructor and length
            d = [1., 1., 1.]
            length = 100 + 500 * ( clu.getEnergy() - emin ) / ( emax - emin )

#           setMag(length)
            mag = sqrt(d[0]*d[0] + d[1]*d[1] + d[2]*d[2])
            d = [d[0] * length/mag, d[1] * length/mag, d[2] * length/mag]

#           setTheta(clu.getITheta())
            mag = sqrt(d[0]*d[0] + d[1]*d[1] + d[2]*d[2])
            if d[0] == 0.0 and d[1] == 0.0 :
                phi = 0.0 
            else : 
                phi = atan2(d[1], d[0])
            d = [mag*sin(clu.getITheta())*cos(phi), 
                 mag*sin(clu.getITheta())*sin(phi),
                 mag*cos(clu.getITheta())
                 ]

#           setPhi(clu.getIPhi())
            xy = sqrt(d[0]*d[0] + d[1]*d[1])
            d = [xy*cos(clu.getIPhi()), xy*sin(clu.getIPhi()), d[2]]

#           Constructor
            dp = [x+d[0], y+d[1], z+d[2]]
            dm = [x-d[0], y-d[1], z-d[2]]

            pced_line_ID( dp[0], dp[1], dp[2],  
                          dm[0], dm[1], dm[2],
                          layer, 1, color, clu.id() );	 

### END OF drawClusters ############################################

def drawReconstructedParticle(col) :
    """Draws collections of type ReconstructedParticle.

    Argument:
    Collection of type 'ReconstructedParticle'.
    """

#   For increased speed (avoid dots in Python).
    cut = col.cut
    colors = g.colors
    lenColors = g.lenColors
    marker = col.marker
    layer = col.layer
    size = col.size
    hitsize = col.hitsize
    clustersize = col.clustersize
    bField = g.BField
    helix_max_r = g.helix_max_r
    helix_max_z = g.helix_max_z

    for i, particle in enumerate(col) :
        if cut( particle ) :
            color = colors[ i % lenColors ]

            tracks = particle.getTracks()
            clusters = particle.getClusters()

#           Long notation for increased speed
            px  = particle.getMomentum()[0]
            py  = particle.getMomentum()[1]
            pz  = particle.getMomentum()[2]
            ene = particle.getEnergy()
            type = particle.getType()

#           Stores individual number for picking
            g.objects[ particle.id() ] = particle


#           Is the if necessary? Re-add if there are problems
#            if len(clusters) :
            for cluster in clusters :
                hits = cluster.getCalorimeterHits()
                for hit in hits :
#                   Long notation for increased speed
                    x = hit.getPosition()[0]
                    y = hit.getPosition()[1]
                    z = hit.getPosition()[2]
                    pced_hit_ID(x, y, z, marker, layer ,
                                clustersize, color, particle.id())

#           Is the if necessary? Re-add if there are problems
#            if len(tracks) :
            for track in tracks :
#               collect hits from all track segments
                hits = getHitsFromTracks(track)

                if len(hits) :
                    for hit in hits :
#                       Long notation for increased speed
                        x = hit.getPosition()[0]
                        y = hit.getPosition()[1]
                        z = hit.getPosition()[2]
                        pced_hit_ID(x, y, z, marker, layer, 
                                    hitsize, color, particle.id())
#               no hits
                else :
#                   Duplicate of line 308, any reason not to delete it?
#                    px = particle.getMomentum()[0]
#                    py = particle.getMomentum()[1]
#                    pz = particle.getMomentum()[2]

#                   start point
                    refx = 0.0
                    refy = 0.0
                    refz = 0.0

#                   helix
                    charge = particle.getCharge()
#                    bField = Global::GEAR->getBField().at(  gear::Vector3D(0,0,0)  ).z() 
                    drawHelix(bField, charge, refx, refy, refz, px, py, pz, 
                              marker, layer, size, color, 0.0, 
                              helix_max_r, helix_max_z, part.id() )

### END OF drawReconstructedParticle ###############################

####################################################################
####################################################################

# readEvent, drawEvent, drawEvents and idleLoop

def drawEvent(evt):
    """Loops over all collections within an event and initiates drawing.

    Argument:
    Event, loaded from lcio file

    Loops over all collections in the event, attaches the 
    drawing parameters, tests if drawing the collection is 
    desired, invokes drawing, and applies the description
    to the layer.
    """

    totTime = time.time()

    clearLayerDescription()
    names = evt.getCollectionNames()

    for name in names :
#        colTime = time.time()

        col = evt.getCollection( name )
        col.name = name
        col.type = col.getTypeName()
        addDrawingParameters( col )

#       Draw only collections enabled for drawing
        if col.draw == True :
            try :
                col.callDraw(col)
                addLayerDescription( name, col.layer )
#                print ("Drawing " + col.type + " " + name + 
#                       " on layer " + str(col.layer))
#           If the attribute collDraw was not set properly in the Maps
            except AttributeError :
                print ( col.type + " " + col.name + 
                        " has no drawFunction associated.")
#        else :
#            print col.type + " " + col.name + " is disabled, skipping."

#        print "\tExecution time " + "%.9f" % (time.time()-colTime)
            
    g.userDraw()

    applyLayerDescription()
    print "Execution time drawEvent(): " + "%.9f" % (time.time()-totTime)

### END OF drawEvent ###############################################

def readEvent(evt = -1, run = -1):
    """Reads the event from the lcio file"""

    if run != -1 and evt != -1 :
        g.event = g.reader.reader.readEvent( run, evt )
    else:
        g.event = g.reader.next()
    return g.event

### END OF readEvent ###############################################

def drawEvents(prev = False):
    """Initiates reading and drawing of the events.
    
    Argument:
    Draws previous event if True.
    """

    if not prev:
        event = readEvent() 
    else:
        event = readEvent(g.event.getEventNumber()-1,
                          g.event.getRunNumber()
                          )
#    event = readEvent() 

    while( event != 0 ) :
        g.ced.ced_new_event()
        g.marlinutil.drawDetectorFromGearFile( g.gearFile )
        drawEvent( event )
        g.ced.ced_send_event()

        cmd = idleLoop() 

        if( cmd == '' ):
            print " reading next event ..." 
        if( cmd == 'q' ):
            break ;

        g.objects.clear()
        event = readEvent() 


### END OF drawEvents ##############################################
       
# TODO: Could be prettier, hopefully (raw_input?)
def idleLoop():
    """Watches user input, next, quit or picking"""

    text = "event " + str(g.event.getEventNumber()) + " " +  str(g.event.getRunNumber())  
    text = text + " - hit enter to draw next - [q] to quit - double click on objects for picking" 

#  input loop (keyboard or mouse)
    print text
    while True:
        i, o, e = select.select( [sys.stdin], [], [], 0.5 )
        if(i):
            c = sys.stdin.readline().strip()
        else:
            c = -1
        if( c == 'q' or c == '' ):
            break ;
        if( c != -1 ):
            print text
        pid = g.ced.ced_selected_id_noblock();
        if pid != -1:
            print ' picked object with Id: ' , pid 
            if( pid in g.objects ):
                print UTIL.toString( g.objects[ pid ] ) 
    return c

### END OF idleloop ################################################

####################################################################
####################################################################

# Functions for userfriendly interactive mode

def redraw() :
    """Redraws the current Event."""

    g.ced.ced_new_event()
    drawEvent( g.event )
    g.ced.ced_send_event()

def next() :
    """Draws the next Event."""

    drawEvents()

# TODO: Does this work?
def prev() :
    """Draws the previous Event."""

    drawEvents(True)

def set( col, atr, val ) :
    """Sets the attribute of a collection to value.

    The collection can either be a name, a type, "all" or "default".
    """
    if col == "default" :
        g.dm[ atr ] = val
        print "Set g.dm[ " + atr + " ] to " + str(val)
    elif col == "all" :
        for collection in g.nm :
            g.nm[ collection ][ atr ] = val
            print "Set g.nm[ " + collection + " ][ " + atr + " ] to " + str(val)
    elif col in g.tm :
        g.tm[ col ][ atr ] = val
        print "Set g.tm[ " + col + " ][ " + atr + " ] to " + str(val)
    elif col in g.nm :
        g.nm[ col ][ atr ] = val
        print "Set g.nm[ " + col + " ][ " + atr + " ] to " + str(val)
    else :
        print "Entry not found, existing keys: "
        print sorted(g.nm) + sorted(g.tm)
        print "'all', 'default'"

def redrawWith( col, atr, val ) :
    """Redraws the Event with the modified settings."""

    set( col, atr, val )
    redraw()

def reset( *maps ) :
    """Resets one or multiple maps.

    Argument can be "all" or a combination of "dm", "tm", "nm".
    """
    if "all" in maps :
        buildDm()
        buildTm()
        buildNm()
    else :
        if "dm" in maps :
            buildDm()
        if "tm" in maps :
            buildTm()
        if "nm" in maps :
            buildNm()

def enable( col ) :
    """Includes the collection when drawing."""
    set( col, "draw", True )

def disable( col ) :
    """Excludes the collection when drawing. Improves performance."""
    set( col, "draw", False )

def picking() :
    """Re-enters the command loop for picking."""
    idleLoop()

####################################################################
####################################################################

# Load config file
execfile( "pyced.cfg.py" )

# Build parameter maps
buildDm()
buildTm()
buildNm()


execfile( "pyced_user.py")

####################################################################

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
