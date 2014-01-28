####################################################
#
# user defined (cut) functions for pyced
#
# 
###################################################
from math import *

## cut function for reconstruced particles
def rcislepton(rc):
    pdg = abs( rc.getType() )
    return ( pdg==15 or pdg==13 or pdg==11 )

## cut function for reconstruced particles
def mcislepton(mc):
    pdg = abs( mc.getPDG() )
    return ( pdg==15 or pdg==13 or pdg==11 )


## example cut for displaying tracks
def trackCut( trk ):

    bField = 3.5
    alpha =  2.99792458e-4 * bField

    pt = abs( alpha / trk.getOmega()  )

    tl = trk.getTanLambda()

    costh = tl / sqrt( 1 + tl * tl )

    tsfH = trk.getTrackState( EVENT.TrackState.AtFirstHit ) 

    zfH = tsfH.getReferencePoint()[2]

    deltaZ = zfH -  trk.getZ0()

#DEBUG    return ( pt < 1. and abs( costh )<0.05  and abs( deltaZ ) < ( pi * tl / trk.getOmega()  ))

    #---- remove "curler" segments ------
    return ( abs( deltaZ ) < ( pi * tl / trk.getOmega()  ) )


###################################################################
#
# apply cuts - per type or per collection name
##################################################################

## draw only reconstructed particles that have a lepton id
#set( "ReconstructedParticle", "cut" , rcislepton )

## draw only monte carlo  particles that have a lepton id
#set( "MCParticle", "cut" , mcislepton )

set( "Track", "cut" ,  trackCut )
