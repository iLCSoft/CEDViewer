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
