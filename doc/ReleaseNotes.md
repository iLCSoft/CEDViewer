# v01-20

* 2024-07-23 tmadlener ([PR#28](https://github.com/iLCSoft/CEDViewer/pull/28))
  - Add key4hep based CI workflows

* 2024-07-23 Leonhard Reichenbach ([PR#27](https://github.com/iLCSoft/CEDViewer/pull/27))
  - Added the processor parameter `DrawMCParticlesCreatedInSimulation` to enable drawing of the MCParticles that were created in the simulation.

# v01-19-01

* 2022-06-22 Thomas Madlener ([PR#23](https://github.com/iLCSoft/CEDViewer/pull/23))
  - Make sure that `ced2go` can be run without the `-n` argument even if `CED_HOST` is set. Fixes #21
  - Remove some compatibility code for python < 2.4

* 2022-03-11 Bohdan Dudar ([PR#20](https://github.com/iLCSoft/CEDViewer/pull/20))
  - Fix typo bug in the ced2go local host argument option

# v01-19

* 2021-11-01 Thomas Madlener ([PR#19](https://github.com/iLCSoft/CEDViewer/pull/19))
  - Fix coverity CI workflow

* 2021-11-01 Bohdan Dudar ([PR#17](https://github.com/iLCSoft/CEDViewer/pull/17))
  - Fix position of the tracks displated by the DSTViewer (fixes #16)
  - Fix CEDViewer `Weff-c++` warnings

* 2021-10-29 Thomas Madlener ([PR#18](https://github.com/iLCSoft/CEDViewer/pull/18))
  - Migrate CI to github actions

# v01-18

* 2021-06-30 scott snyder ([PR#15](https://github.com/iLCSoft/CEDViewer/pull/15))
  - Fix miscellaneous compilation warnings.

* 2021-06-30 scott snyder ([PR#14](https://github.com/iLCSoft/CEDViewer/pull/14))
  - Updated ced2go for python3 compatibility.

# v01-17-01

* 2020-04-13 Frank Gaede ([PR#13](https://github.com/iLCSoft/CEDViewer/pull/13))
  - make compatible w/ -std=c++17 
         - needed on macos w/ clang

# v01-17

* 2018-11-06 Marko Petric ([PR#9](https://github.com/iLCSoft/CEDViewer/pull/9))
  - DDCEDViewer
    - Add option to color reconstructed particles by energy
      - the color range is from lowest blue to highest red, with possibility of the user to change brightens and saturation of the color pallet.
      - the color range can be fixed from a minimal energy to a maximal for all events or auto adopted event by event to span from min to max energy of a given/current event
    - Add option for scaling line thickness of helices and marker size, this is needed for good image dumps, since the displayed sizes differ from the ones in the dumped image

# v01-16

* 2018-08-21 Frank Gaede ([PR#7](https://github.com/ilcsoft/CEDViewer/pull/7))
  - improved CEDViewer/DDCEDViewer (also used by ced2go)
        - print run and event number at end of event
        - print only shown collections 
        - added new collections ( mostly calo hits) from ILD mass production

* 2018-08-20 Akiya Miyamoto ([PR#6](https://github.com/ilcsoft/CEDViewer/pull/6))
  - Added a command line option, "-n 1".   If this option is specified, glced is not executed. This option is useful for a case to run glced at local client and run CEDViewer at remote host.

# v01-15

* 2017-09-14 Shaojun Lu ([PR#5](https://github.com/iLCSoft/CEDViewer/pull/5))
  -   Update the  ced2go template for ILD_l4_v02.
      - add the Ecal and Hcal collection names from ILD_l4_v02.
      - remove duplicated  LumiCalCollection.

# v01-14-01

* 2017-07-06 Frank Gaede ([PR#4](https://github.com/iLCSoft/CEDViewer/pull/4))
  - add "EcalXXXCollection" to ced2go  template for drawing SimCalorimeterHits
  - fix the drawing of non-prompt neutral MCParticles in CEDViewer

# v01-14

* 2017-06-20 Andre Sailer ([PR#3](https://github.com/iLCSoft/CEDViewer/pull/3))
  - Adapt to changes in namespaces and LCDD -->  Detector
  - Ignore warnings from external header files

# v01-13

# v01-12

F. Gaede
   - made compatible with c++11
   - removed -ansi -pedantic -Wno-long-long
   - fixed narrowing in initializer lists

# v01-11

  - ced2go
    - increased size of tracker hits drawn for tracks
    - added optional EventSelector to ced2go/ced2go-template-DD4.xml
     ( commented out )
  
 - DDCEDViewer.cc
   - added printout of shown collections and event number to DDCEDViewer
   - fixed duplicate drawing of collections in DDCEDViewer
   - fixed some log levels
   - introduced flag DrawEllipsoidForPFOClusters
     (default is Disabled)
   - draw calorimeter hits for PFO clusters ...


T.Quast, RWTH Aachen
   
  - ced2go
   - Allow for Marlin options '--', '-n' option removed due to default 'global.SkipNEvents=...' option
   - Extend ced2go's help message
  - DDCEDViewer
   - Jet collections to be defined in DrawInLayer, jets are drawn if the collection contains the substring 'Jet'
   - use logarithmic (ROOT) color map for reconstructed clusters' energy visualization
   - energy-weighted track lengths for jet entities
   - energy-weighted track lengths for jet entities
   - Draw energy clusters as ellipsoids. Energy scaling with colors, (energy weighted) spatial extension <--> size of according ellipsoid
   - Improve jet visualization for CEDViewer
   - removed GEAR dependencies, merge of processors, adjusted ced2go and add skipNEvents option, 
     implement jet visualization  (commit requires tip version of MarlinUtil



# v01-10
    T.Quast, RWTH Aachen

    - add DrawDetectorDD4hep
    - general detector (CLIC & ILD) drawing for the event display without GEAR
    - adapted ced2go to transparently use DD4hep if a compact file is given
      with -d instead of a Gear file    
    - add option to ced2go to draw tracking surfaces from DDRec (-s)

    => depends now on DD4hep/Root


# v01-09-01 patch release

   - M.Petric:
     - updated ced2go template files for new CLIC simulation model


# v01-09

       - adapted CED event display for CLIC-Sid detector ( M.Petric/FG)
          - added new processor DrawCLIC
       	     - overwrites the drawDetectorFromGEAR method in MarlinUtil
               ( eventually this code should go to MarlinUtil)
	     - uses different colors for detector (hardcoded)
          - added option DrawHelixForPFOs to CEDViewer
            -> draws the helix in the same colors as the PFO hits

          - added gear file and template to ced2go, use e.g.:
            ced2go -d ./ced2go/gear_CLIC_SiD.xml -t ./ced2go/ced2go-CLIC-template.xml h_nunu_rec_4019_20994.slcio


v01-08 - does not exist 

# v01-07-02
      - some minor improvements to pyced 
         - added userDraw() method for
           drawing additional objects
         - added example for track cut function

# v01-07-01
   -  added pyced python event display
      (see: ./pyced/README )
    
   - no change to the rest of the code

# v01-07

   - changed calling attributes of ced_hit_ID to newer version
     -> needed for CED v01-09 which has no function overloading 
        in library anymore (c-style)
 
  - changed the default conventions for the track state to be drawn. Now the following processor options are valid
    -1: do not draw helix
     0: default track state (the 0-th element in the track state vector)
      (the rest options as before)
     1: AtIP
     2: AtFirstHit
     3: AtLastHit
     4: AtCalorimeter
    If no track state were found, exit with an error

 - avoid wrong values for the helix parameters for the case of small curvature
   or zero magnetic field by assigning large value of Pt to the track

 - update for the multi-module support via generic GEAR interface


# v01-06-01

 - display all tracker hits from all segments
    for Tracks and ReconstructedParticles

  - use last track segment (if any) to get
    track states at last hit and at calorimeter

# v01-06

  - removed TrackerRawViewer 
     - lives now in MarlinTPC


# v01-05-02

 - fixed FTD hit collection names 
 - draw all tracks and clusters for recoparticles.


# v01-05-01

    - added drawing of the strip measurment for SIT TrackerHitPlane
    - fixed ced2go in case variable STANDARDCONFIG is not set
    - made compatible with clang++


# v01-05

    CEDViewer:
      added parameter WaitForKeyboard
    ced2go:
      modifed collections names to by the same as in bbudsc_stdreco.xml
      allow to specifiy steering file template on the command line with -t


# v01-04-01
     - fixed collection names for ced2go
       (EcalBarrel/EndcapSiliconCollection)


# v01-04

    - Added option to supply limits to helix drawing, other than the extent of the TPC
    - Added option for drawing MC truth info for particle gun events
    - switched to TPC drawing with surfaces and cuts

    - adopted steering example steering files:
       - renamed ChannelMapperProcessor to ChannelMappingProcessor
       - changed raw data file pathes on afs

    - changed parameter DrawHelixForTrack to int:
       draw a helix for Track objects:
       0 none, 1: atIP, 2: atFirstHit, 3: atLastHit, 4: atCalorimeter

    - changes in DSTViewer processor:
      - code cleanup
      - use streamlog instead of cout
      - added parameter DrawDetectorID 
      - this processor needs to be revisited ....

    - color code SimTrackerHits and SimCalorimeterHits by the MCParticle 
    - changed drawing of track helix color to light grey

    - added parameter MCParticleEnergyCut 
    - fixed drawing of neutrals 
         - neg. pz particles where drawn in wrong direction
         - start drawing particles at their vertex
    - added lightgrey color for neutrinos

    - added drawing of TrackerHitPlanae and TrackerHitZCylinder

    - added various color schemes for displaying particles, tracks and clusters
      -> steering parameter ColorScheme 
    - print evt/run numbers at the end of the event

    - cmake changes:
        - added cmake policy CMP0008
        - removed obsolete install of header files



# v01-03

   -  added new TrackerRawViewer

       - displays raw data from LCTPC large prototype

       - needs CEDViewer to be build with MarlinTPC and LCCD
       - see example/viewTrackerRawData.xml and example/reconstructLCTPCRawData.xml 
       

   - CEDViewer
      - added new steering parameters 
        "DrawHelixForTrack" and "DrawDetectorID"

      - fixed display of ReconstructParticles (did not respect
	the marker flags - patch provided by T. Tanabe)


# v01-02

 - cmake changes:
    - simplified CMakeLists.txt
    - improvements in dependency handling:
        - removed dependencies GEAR, LCIO and streamlog (now exported through Marlin)
        - removed dependencies CED and CLHEP (now exported thorough MarlinUtil)
        - removed dependency to GSL (not used anywhere in CEDViewer). Dependency was propagated 
            from a "bug" in MarlinUtil header files ( #include <gsl/...> statements )
        - exchanged CMakeModules dependency with new package ILCUTIL

    - removed CEDViewerConfig.cmake.in
    - removed BuildSetup.cmake
 - removed old makefile
 - removed old steering file


# v01-01
    new features: (H. Hoelbe, DESY)
    -  CEDViewer: 
        -  Adds now his layer description to CED
        -  Now also displays Reconstructed Particles. If there are no hits,
           than draw the helix. 
    -  DSTViewer:         
        -  Adds now his layer description to CED
        -  Cones have now a central line for picking    
    -  GenericViewer: 
        -  Adds now his layer description to CED

    [ -  ced2go:
        -  A script which simplify the use of Marlin with CED in some cases.
           ced2go needs the name of the LCIO file as an argument. Than it 
           look to this file find out which gearfile is needed, create an 
           steering file with this parameters. Than it starts CED and Marlin.  
    ]

    changes/bug fixes: (H. Hoelbe, DESY)
    -  CEDViewer:
        - Helix part of tracks was not pickable. fixed

# v01-00
     - new release of CEDViewer with 'picking' functionality 
        - implemented  in CEDViewer and GenericViewer
        - see ./doc/CEDPicking.pdf for details   


# v00-07-02

    - bug fix: incorrect library version numbers

# v00-07-01

    - replaced "const double Pi = acos(-1.)" with M_PI from math.h

# v00-07
    - implemented a DSTViewer (S. Darasz):
        draws ReconstructedParticles as helices and straight lines;
        Clusters as Cylinders (Ellipsoids) and different jet hypotheses as cones
    - made cmake 2.6 compliant
    - added 32 bit compatibility build option
    
# v00-06
    - made dependent on CED directly 
      -> need CED  >= v00-04-01  
      -> need MarlinUtil >= v00-10 

# v00-05
   - VertexViewer processor added   (A.Raspereza)
     allows detailed vertex views in r-phi plane
 
   - impoved cmake file (J.Engels)

   - minor fixes

# v00-04
 
 - removed HepPDT dependency
 

# v00-03

      - cmake is now default build tool:  (J.Engels)
        # edit BuildSetup.cmake as needed
        mkdir build ; cd build
        cmake -C ../BuildSetup.cmake ..
        make install
     -> creates plugin library $CEDViewer/lib/libCEDViewer.so


# v00-02
 
  - added cmake support (epxerimental)
  - renamed Physical_Geometrical_database to PGdb (see MarlinUtil)
