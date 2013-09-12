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
