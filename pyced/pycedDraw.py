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
