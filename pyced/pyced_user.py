####################################################
#
# user defined (cut) functions for pyced
#
# 
###################################################

## cut function for reconstruced particles
def rcislepton(rc):
    pdg = abs( rc.getType() )
    return ( pdg==15 or pdg==13 or pdg==11 )

## cut function for reconstruced particles
def mcislepton(mc):
    pdg = abs( mc.getPDG() )
    return ( pdg==15 or pdg==13 or pdg==11 )


## draw only reconstructed particles that have a lepton id
#set( "ReconstructedParticle", "cut" , rcislepton )

## draw only monte carlo  particles that have a lepton id

#set( "MCParticle", "cut" , mcislepton )

