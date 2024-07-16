/***********************************************************************************************
Single processor for drawing detector geometry and LCIO data; removed dependence on GEAR.

author: Thorben Quast (CERN Summer Student 2015)
date: 10 August 2015
***********************************************************************************************/

#ifndef CEDViewer_h
#define CEDViewer_h 1

#include "marlin/Processor.h"
#include "lcio.h"

#include <DD4hep/Detector.h>

#include <string>
#include <vector>

using namespace lcio ;
using namespace marlin ;



/**Helper struct for drawing collections*/
struct DrawParameters{
    DrawParameters(const std::string& colName, int size, int marker, int layer ) :
    ColName( colName ),
    Size( size ),
    Marker( marker ),
    Layer( layer ) {
    }
    std::string ColName ;
    int Size ;
    int Marker ;
    int Layer ;
};


class DDCEDViewer : public Processor {

 public:
  virtual Processor*  newProcessor() { return new DDCEDViewer ; }

  DDCEDViewer() ;
  /** Called at the begin of the job before anything is read.
   * Use to initialize the processor, e.g. book histograms.
   */
  virtual void init() ;

  /** Called for every run.
   */
  virtual void processRunHeader( LCRunHeader* run ) ;

  /** Called for every event - the working horse.
   */
  virtual void processEvent( LCEvent * evt ) ;


  virtual void check( LCEvent * evt ) ;


  /** Called after data processing for clean up.
   */
  virtual void end() ;

  //helpers
  void printParticle(int id, LCEvent * evt);
  bool detailledDrawing(std::string detName);


  /** LCIO collection based viewer function for CED event dislplay (A. Zhelezov).
 *  Define the collection name, marker, size and optionally the layer number in the steering file.
 *  For 'DrawCollection' a default layer number in CED is chosen.
 *
 *  @param DrawCollection - collection to be displayed ( ColName, marker type[0-2] )
 *  @param DrawInLayer    - collection to be displayed ( ColName, marker type[0-2], size)
 *
 *  @author F.Gaede, DESY
 *  @version $Id: CEDViewer.h 4773 2015-02-27 09:01:23Z gaede $
 *
 * Modified original code by Thorben Quast (CERN Summer Student 2015) on 07 August 2015.
 * - removed GEAR dependence
 * - refactored drawing of straight lines up to calorimeters
 */
  void drawDD4LCIO(LCEvent * evt, dd4hep::Detector& theDetector);
  //convenience function to improve readibility of main code
  void drawCluster(dd4hep::Detector& theDetector, int& layer, unsigned& np, std::string colName, int& marker, LCCollection* col, int& size);
  void drawTrack(dd4hep::Detector& theDetector, int& layer, unsigned& np, std::string colName, int& marker, LCCollection* col, int& size);
  void drawMCParticle(dd4hep::Detector& theDetector, int& layer, unsigned& np, std::string colName, int& marker, LCCollection* col, int& size);
  void drawSIMTrackerHit(int& layer, unsigned& np, std::string colName, int& marker, LCCollection* col, std::vector<int>& _colors, int& size);
  void drawSIMCalorimeterHit(int& layer, unsigned& np, std::string colName, int& marker, LCCollection* col, std::vector<int>& _colors, int& size);
  void drawTrackerHit(int& layer, unsigned& np, std::string colName, int& marker, LCCollection* col, int& size);
  void drawCalorimeterHit(int& layer, unsigned& np, std::string colName, int& marker, LCCollection* col, int& size);
  void drawReconstructedParticle(dd4hep::Detector& theDetector, int& layer, unsigned& np, std::string colName, int& marker, LCCollection* col, int& size);
  void drawJets(dd4hep::Detector& theDetector, int layer, std::string colName, LCCollection* col);

 protected:

  constexpr static int ncol = 20 ;
  constexpr static int nscheme = 10 ;
  constexpr static int Red         =   0 ;
  constexpr static int Orange      =   1 ;
  constexpr static int Plum        =   2 ;
  constexpr static int Violet      =   3 ;
  constexpr static int Blue        =   4 ;
  constexpr static int LightBlue   =   5 ;
  constexpr static int Aquamarine  =   6 ;
  constexpr static int Green       =   7 ;
  constexpr static int Olive       =   8 ;
  constexpr static int Yellow      =   9 ;
  constexpr static int Dark        =   10 ;
  constexpr static int Light       =   11 ;
  constexpr static int Classic     =   12 ;

  StringVec _drawCollections{} ;
  StringVec _drawCollectionsLayer{} ;
  std::vector< DrawParameters > drawParameters{} ;

  //general options
  int _nRun = 0 ;
  int _nEvt = 0;
  //lcio options
  bool      _usingParticleGun            = false ;
  bool      _drawMCParticlesCreatedInSimulation = false;
  int       _drawHelixForTracks          = 0 ;
  int       _colorScheme                 = 10 ;
  bool      _colorEnergy                 = false ;
  double    _colorEnergyMin              = 0.0 ;
  double    _colorEnergyMax              = 35.0 ;
  double    _colorEnergySaturation       = 0.8 ;
  double    _colorEnergyValue            = 0.8 ;
  bool      _colorEnergyAuto             = false ;
  double    _scaleLineThickness          = 1 ;
  double    _scaleMarkerSize             = 1 ;
  float     _mcpECut                     = 0.001 ;
  float     _helix_max_r                 = 2000 ;
  float     _helix_max_z                 = 2500 ;
  bool      _useTrackerForLimitsOfHelix  = true ;
  int       _waitForKeyboard             = 1 ;
  int       _drawHelixForPFOs            = 0 ;
  int       _useColorForHelixTracks      = 0 ;
  int       _drawEllipsoidForPFOClusters = 0 ;
  IntVec    _colors{} ;

  //detector options
  bool _begin          = false;
  StringVec _detailled{};
  StringVec _jets{};
  bool _surfaces       = false;
} ;


/***********************************************************************************************
Helper or 'utility' functions & structures for DDCEDViewer processor.

author: Thorben Quast (CERN Summer Student 2015)
date: 12 August 2015
***********************************************************************************************/

#include "marlin/Processor.h"
#include <string>

//Includes for detector drawing
#include "DD4hep/Detector.h"
#include "DDRec/DetectorData.h"

//structure for calculating the track length of given particles
struct CalorimeterDrawParams {
    double r_inner, delta_r, z_0, delta_z;
};

/***lcio draw helpers***/
//get the outer extents of the tracker
double* getTrackerExtent(dd4hep::Detector& theDetector);

//get the outer extents of the yoke
double* getYokeExtent(dd4hep::Detector& theDetector);

//calculates and returns the relevant calorimeter parameters for track length calculations.
CalorimeterDrawParams getCalorimeterParameters(dd4hep::Detector& theDetector, std::string name, bool selfCall = false);

//It suffices to perform the calculations in the first quadrant due to the detector's symmetry.
//The signs of the tracks' directions are ultimately determined by the momenta.
double calculateTrackLength(std::string type, dd4hep::Detector& theDetector, double x, double y, double z, double px, double py, double pz);

int returnRGBClusterColor(float eneCluster, float cutoff_min, float cutoff_max, int color_steps, char scale, int colorMap);

#endif



