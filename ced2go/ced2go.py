#!/usr/bin/python
import sys
import re
import subprocess
import os
import signal
import time
import getopt


#Help output
def help():
    print "Usage: " + str(sys.argv[0]) + " [OPTION]  <LCIO Input File>"
    print "Options: "
    print "     -v            Viewer1"
    print "     -v            Viewer2"
    print "     -d            Gearfile"
    print "     -h            Help"
    print "All options are optional"
    print "Viewers are e.g. CEDViewer, DSTViewer, GenericViewer"

#import shlex
#from Tkinter import *
#import threading

#The global variables, paths and commands
viewer   = "CEDViewer" #default viewer

#templateXMLFile = "/afs/desy.de/user/h/hhoelbe/workspace/ilcsoft/v01-07/Wrapper/template.xml"
templateXMLFile = "./template.xml"

if not os.path.isfile(templateXMLFile):
    print "Error: Template File \"" + templateXMLFile + "\" not found!"
    sys.exit()


#extractDetector="/afs/desy.de/user/h/hhoelbe/workspace/ilcsoft/v01-07/Wrapper/build/extractdetector"
extractDetector="./build/extractdetector"

if not os.path.isfile(extractDetector):
    print "Error: Helper Tool \"" + extractDetector + "\" not found!"
    sys.exit()

steeringFile="steering.xml"
CED= "glced"
Marlin="Marlin"

try:
    standartconfig=os.environ["STANDARDCONFIG"]
except: 
    print "Error: Please source the init ilcsoft script first"
    sys.exit(1)

gearfiles={'ILD_00':standartconfig + '/mc2008/gear_ILD_00.xml',
        'LDC01_06Sc':standartconfig+ '/mc2008/gear_LDC01_06Sc.xml',
        'LDCPrime_02Sc':standartconfig+'/mc2008/gear_LDCPrime_02Sc.xml',
        'LDC_GLD_01Sc':standartconfig+ '/mc2008/gear_LDC_GLD_01Sc.xml'
}

keys = {'$LCIOInputFiles$':'', 
        '$GearXMLFile$':'',
        '$Viewer$':''
       }

#Parse the command line arguments
try: 
    opts, args = getopt.getopt(sys.argv[1:], "hv:d:", "help") 
    #print 'Opts:',opts
    #print 'Extra parameters:',args
except getopt.GetoptError, err:
    print str(err)
    help()
    sys.exit(2)

detector = ""
for key in opts:
    #print key[0] + str(key[1])
    if key[0] == "-d":
        detector=key[1]
        if not os.path.isfile(detector):
            print "Error: Gearfile \"" + detector + "\" not found!"
            sys.exit()
    if key[0]  == "-v":
          keys['$Viewer$'] += "\n<processor name=\"My" + str(key[1]) + "\"/>"
    if key[0] in ("-h", "--help"):
        help()
        sys.exit(0)

if keys['$Viewer$'] == '': #default
    keys['$Viewer$'] = "<processor name=\"My" + viewer + "\"/>" 


#class Thread1(threading.Thread):
#    def run(self):
#        print "thread 1 startet and block"
#        while True:
#            pass

#Thread1().start()
#print "main sleep 5sec"
#time.sleep(5)
#sys.exit()

#check if lcio file is given
if not len(args) == 1:
    help()
    sys.exit()


#check if lcio file exists
if not os.path.isfile(args[0]):
    print "Error: LCIO File \"" + sys.argv[1] + "\" not found!"
    help()
    sys.exit()
    
keys['$LCIOInputFiles$']=args[0]



if detector != "": #user defined gearfile
     keys['$GearXMLFile$'] = detector
else:  #find out the detector name used in the lcio file
    detector= str(subprocess.Popen([extractDetector,keys['$LCIOInputFiles$']], stdout=subprocess.PIPE).communicate()[0])
    detector=detector.rstrip() #remove newline

    #print '"' + detector + '"'

    tmp = detector.split("\n")
    if len(tmp) != 1:
        print "Warning: LCIO File contains more than one detector model, try the first one"
        detector=tmp[0]
    try: 
        keys['$GearXMLFile$'] = gearfiles[detector] 
    except(KeyError):
        print "Error: Unknown detector in LCIO file"
        sys.exit(1)

#root = Tk()

#w=Label(root, text="hello world")
#w.pack()

#root.mainloop()


#replace the words from dict with the keywords in the template file
templateFH = open(templateXMLFile,'r')
xmlFH = open(steeringFile,'w')

for line in templateFH:
    for key in keys: 
        line= re.sub(re.escape(key),keys[key],line);
    xmlFH.write(line) 
xmlFH.close()
templateFH.close()

#start ced and marlin
CEDobj = subprocess.Popen([CED,'-geometry', '600x400'])
Marlinobj = subprocess.Popen([Marlin, steeringFile])

#terminate the both processes when one is closed
while True: 
    if Marlinobj.poll() != None:
        os.kill(CEDobj.pid,15)
        break

    if CEDobj.poll() != None: 
        os.kill(Marlinobj.pid,15)
        break

    time.sleep(0.5) #sleep for 0.5 seconds

print "Done"
