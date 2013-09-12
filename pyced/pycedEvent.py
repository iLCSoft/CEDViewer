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

    event = readEvent() 

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
#       TODO: Does this work?
        if not prev:
            event = readEvent() 
        else:
            event = readEvent(g.event.getEventNumber()-1,
                              g.event.getRunNumber()
                              )

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
