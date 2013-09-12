# This is the config file for pyced, a Python class to control the CED.
# All global variables are set here, change them if you know what you are doing.

g.gearFile              = "gear_ILD_o1_v05.xml"
g.host                  = "localhost"
g.port                  = 7286 
g.colorScheme           = -1
g.BField                = 3.5


# --- draw a helix for Track objects: 
#     0 none, 1: atIP, 2: atFirstHit, 3: atLastHit, 4: atCalorimeter
g.drawHelixForTracks    = 0
g.helix_max_r           = 2000.
g.helix_max_z           = 2500.


# --- For drawing MCParticles ---
# (Global::GEAR->getEcalBarrelParameters().getExtent()[0] + 
#   Global::GEAR->getEcalEndcapParameters().getExtent()[1])/2.;
g.ecalR = 3000.
# std::abs(Global::GEAR->getEcalEndcapParameters().getExtent()[2] + 
#   Global::GEAR->getEcalEndcapParameters().getExtent()[3])/2.;
g.ecalZ = 3000.
# (Global::GEAR->getHcalBarrelParameters().getExtent()[0] + 
#   Global::GEAR->getHcalEndcapParameters().getExtent()[1])/2.;
g.hcalR = 3000.
# std::abs(Global::GEAR->getHcalEndcapParameters().getExtent()[2] + 
#   Global::GEAR->getHcalEndcapParameters().getExtent()[3])/2.;
g.hcalZ = 3000.

# If set to True the generator status is ignored for MCParticles to be drawn
g.usingParticleGun = False


# Set default values for drawing, dm: default map, tm: type map, nm: name map
def buildDm() :
    """Builds the DefaultMap for adding drawing attributes to the Collections.
    
    Required attributes:
    marker (int), layer (int), size (int), draw (bool), cut (callable).
    """
    g.dm = dict(
            marker  = 0,
            layer   = 1,
            size    = 1,
#           To filter out LCRelations etc
            draw    = False,
#           Used to apply specific cuts, implemented are eCut, ptCut,
#           cosThetaCut. Second argument can be "greater" (than cut), 
#           "smaller" (than cut) or a number (range: cut < energy < 
#           cutType), default is "greater".
#           Only apply to Collections with the desired attribute defined.
            cut     = eCut(0)
            )

def buildTm() :
    """Builds the TypeMap for adding drawing attributes to the Collections.

    The attributes of collections of types matching a dictionary entry will be
    updated with the values entered here. This extends and possibly overwrites 
    the values from the DefaultMap.
    
    Possible attributes:
    marker (int), layer (int), size (int), draw (bool), cut (callable).

    Required attributes:
    callDraw (drawFunction).

    Required for certain collection types:
    color (callable), hitsize (int), clustersize (int).
    """
    g.tm = dict(
            TrackerHit = dict(
                    size    = 5,
                    color   = fixedColor("0xee0044"),
                    draw    = False,
                    callDraw = drawHits
                    ),
            TrackerHitZCylinder = dict(
                    size    = 5,
                    color   = fixedColor("0xee0044"),
                    draw    = False,
                    callDraw = drawHits
                    ),
            SimTrackerHit = dict(
                    size    = 5,
                    color   = SimTrackerHitColor,
                    draw    = False,
                    callDraw = drawHits
                    ),
            CalorimeterHit = dict(
                    layer   = 2,
                    size    = 3,
                    color   = fixedColor("0xee0000"),
                    draw    = False,
                    callDraw = drawHits
                    ),
            SimCalorimeterHit = dict(
                    layer   = 2,
                    size    = 3,
                    color   = SimCalorimeterHitColor,
                    draw    = False,
                    callDraw = drawHits
                    ),
            Track = dict(
                    layer   = 3,
                    size    = 3,
                    hitsize = 5,
                    draw    = False,
                    callDraw = drawTracks
                    ),
            Cluster = dict(
                    layer   = 12,
                    hitsize = 3,
                    draw    = False,
                    callDraw = drawClusters
                    ),
            MCParticle = dict(
                    layer   = 0,
                    hitsize = 3,
                    draw    = True,
                    callDraw = drawMCParticles
                    ),
            ReconstructedParticle = dict(
                    size    = 3,
                    layer   = 11,
                    hitsize = 5,
                    clustersize = 5,
                    draw    = True,
                    callDraw = drawReconstructedParticle
                    )
            )




def buildNm() :
    """Builds the NameMap for adding drawing attributes to the Collections.

    The attributes of the collection with a name matching a dictionary 
    entry will be updated with the values entered here. This exetends and 
    possibly overwrites the values from the DefaultMap and the TypeMap.
    
    Possible attributes:
    marker (int), layer (int), size (int), draw (bool), cut (callable),
    callDraw (drawFunction), color (callable), hitsize (int), clustersize (int).

    Empty collections could be deleted but are kept for convenience.
    """
    g.nm = dict(
            VXDCollection               = dict(  ) , 
            SITCollection               = dict(  ) ,
            FTD_PIXELCollection         = dict(  ) ,
            FTD_STRIPCollection         = dict(  ) ,
            TPCCollection               = dict(  ) ,
            SETCollection               = dict(  ) ,
            ETDCollection               = dict(  ) ,

            VXDTrackerHits              = dict( marker = 1, size = 10, layer = 11 ) , 
            SITTrackerHits              = dict( marker = 1, size = 10, layer = 11 ) ,
            TPCTrackerHits              = dict( marker = 1, size =  5, layer = 11 ) ,
            FTDTrackerHits              = dict( marker = 1, size = 10, layer = 11 ) ,

            HcalEndCapsCollection       = dict(  ) , 
            LHcalCollection             = dict(  ) ,
            LumiCalCollection           = dict(  ) ,
            MuonBarrelCollection        = dict(  ) ,
            MuonEndCapCollection        = dict(  ) ,
            EcalBarrelSiliconCollection = dict(  ) ,  
            EcalBarrelSiliconPreShower  = dict(  ) ,  
            EcalEndcapRingCollection    = dict(  ) ,  
            EcalEndcapRingPreShower     = dict(  ) ,  
            EcalEndcapSiliconCollection = dict(  ) ,  
            EcalEndcapSiliconPreShower  = dict(  ) ,  
            HcalBarrelRegCollection     = dict(  ) ,
            HcalEndCapRingsCollection   = dict(  ) ,

            HCALEndcap                  = dict( size =  5, layer = 12 ) ,     
            HCALOther                   = dict( size =  5, layer = 12 ) ,    
            MUON                        = dict( size =  2, layer = 12 ) ,    
            LHCAL                       = dict( size =  3, layer = 12 ) ,      
            LCAL                        = dict( size =  3, layer = 12 ) ,    
            BCAL                        = dict( size =  3, layer = 12 ) ,    
            ECALBarrel                  = dict( size =  2, layer = 12 ) ,    
            ECALEndcap                  = dict( size =  2, layer = 12 ) ,    
            ECALOther                   = dict( size =  2, layer = 12 ) ,    
            HCALBarrel                  = dict( size =  5, layer = 12 ) ,    

            TruthTracks                 = dict( layer =  3 ) ,    
            ForwardTracks               = dict( layer =  4 ) ,    
            SiTracks                    = dict( layer =  5 ) ,    
            ClupatraTracks              = dict( layer =  6 ) ,    
            ClupatraTrackSegments       = dict( size =  4, layer =  6 ) ,
            CluTrackSegments            = dict( size =  3, layer =  6 ) , 
            MarlinTrkTracks             = dict( layer =  7 ) ,    

            PandoraClusters             = dict( size =  3, layer =  8 ) ,
            SeedCluster                 = dict( size =  4, layer =  8 ) , 
            LeftoverClusters            = dict( size =  3, layer =  8 ) ,
            PandoraPFOs                 = dict( size =  3, layer =  9 ) ,
            BCALParticles               = dict( size =  3, layer = 19 ) ,    

            MCParticle                  = dict( size =  3, layer =  0 )    
    )
