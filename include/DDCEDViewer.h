/***********************************************************************************************
Single processor for drawing detector geometry and LCIO data; removed dependence on GEAR.

author: Thorben Quast (CERN Summer Student 2015)
date: 10 August 2015
***********************************************************************************************/

#ifndef CEDViewer_h
#define CEDViewer_h 1

#include "marlin/Processor.h"
#include "lcio.h"
#include <string>
#include <vector>

using namespace lcio ;
using namespace marlin ;
namespace DD4hep{ 
  namespace Geometry{
    class LCDD ;
  }
}



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
  void drawDD4LCIO(LCEvent * evt, DD4hep::Geometry::LCDD& lcdd);
  //convenience function to improve readibility of main code
  void drawCluster(DD4hep::Geometry::LCDD& lcdd, int& layer, unsigned& np, std::string colName, int& marker, LCCollection* col, int& size);
  void drawTrack(DD4hep::Geometry::LCDD& lcdd, int& layer, unsigned& np, std::string colName, int& marker, LCCollection* col, int& size);
  void drawMCParticle(DD4hep::Geometry::LCDD& lcdd, int& layer, unsigned& np, std::string colName, int& marker, LCCollection* col, int& size);
  void drawSIMTrackerHit(int& layer, unsigned& np, std::string colName, int& marker, LCCollection* col, std::vector<int>& _colors, int& size);
  void drawSIMCalorimeterHit(int& layer, unsigned& np, std::string colName, int& marker, LCCollection* col, std::vector<int>& _colors, int& size);
  void drawTrackerHit(int& layer, unsigned& np, std::string colName, int& marker, LCCollection* col, int& size);
  void drawCalorimeterHit(int& layer, unsigned& np, std::string colName, int& marker, LCCollection* col, int& size);
  void drawReconstructedParticle(DD4hep::Geometry::LCDD& lcdd, int& layer, unsigned& np, std::string colName, int& marker, LCCollection* col, int& size);
  void drawJets(DD4hep::Geometry::LCDD& lcdd, int layer, std::string colName, LCCollection* col);
 
 protected:

  static const int ncol = 20 ;
  static const int nscheme = 10 ;
  static const int Red         =   0 ;
  static const int Orange      =   1 ;
  static const int Plum        =   2 ;
  static const int Violet      =   3 ;
  static const int Blue        =   4 ;
  static const int LightBlue   =   5 ;
  static const int Aquamarine  =   6 ;
  static const int Green       =   7 ;
  static const int Olive       =   8 ;
  static const int Yellow      =   9 ;
  
  static const int Dark        =   10 ;
  static const int Light       =   11 ;
  static const int Classic     =   12 ;
  
  StringVec _drawCollections ;
  StringVec _drawCollectionsLayer ;
  std::vector< DrawParameters > drawParameters ;

  //general options
  int _nRun ;
  int _nEvt ;
  //lcio options
  bool      _usingParticleGun ;
  int       _drawHelixForTracks ;
  int       _colorScheme ;
  float     _mcpECut ;
  float     _helix_max_r;
  float     _helix_max_z;
  bool      _useTrackerForLimitsOfHelix;
  int       _waitForKeyboard ;
  int       _drawHelixForPFOs;
  int       _useColorForHelixTracks ;
  IntVec    _colors ;

  //detector options
  bool _begin;
  StringVec _detailled;
  StringVec _jets;
  bool _surfaces;
} ;





/***********************************************************************************************
Helper or 'utility' functions & structures for DDCEDViewer processor.

author: Thorben Quast (CERN Summer Student 2015)
date: 12 August 2015
***********************************************************************************************/

#include "marlin/Processor.h"
#include <string>

//Includes for detector drawing
#include "DD4hep/LCDD.h"
#include "DDRec/DetectorData.h"
using namespace DD4hep::Geometry ;

//structure for calculating the track length of given particles
struct CalorimeterDrawParams {
    double r_inner, delta_r, z_0, delta_z;
};

/***lcio draw helpers***/
//get the outer extents of the tracker
double* getTrackerExtent(DD4hep::Geometry::LCDD& lcdd);

//get the outer extents of the yoke
double* getYokeExtent(DD4hep::Geometry::LCDD& lcdd);

//calculates and returns the relevant calorimeter parameters for track length calculations.
CalorimeterDrawParams getCalorimeterParameters(DD4hep::Geometry::LCDD& lcdd, std::string name, bool selfCall = false);

//It suffices to perform the calculations in the first quadrant due to the detector's symmetry.
//The signs of the tracks' directions are ultimately determined by the momenta.
double calculateTrackLength(std::string type, DD4hep::Geometry::LCDD& lcdd, double x, double y, double z, double px, double py, double pz);

int returnClusterColor(float eneCluster, float cutoff_min, float cutoff_max);


#endif



