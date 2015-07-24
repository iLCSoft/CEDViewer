#ifndef DrawDetectorDD4hep_h
#define DrawDetectorDD4hep_h 1

#include "marlin/Processor.h"
#include "lcio.h"
#include <string>


using namespace lcio ;
using namespace marlin ;


/**  Processor for marlin.
 *   Draw a CED detector from a DD4hep detector model
 * 
 * @author ...
 * @version $Id: $
 */

namespace DD4hep{ 
  namespace Geometry{
    class LCDD ;
  }
}

class DrawDetectorDD4hep : public Processor {
  
 public:
  
  virtual Processor*  newProcessor() { return new DrawDetectorDD4hep ; }
  
  
  DrawDetectorDD4hep() ;
  
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

  //remove 'static' to access new defined parameter
  bool detailledDrawing(std::string detName);
  void drawDD4hepDetector( DD4hep::Geometry::LCDD& lcdd) ;
  
  
 protected:
  bool _begin;
  int _nRun ;
  int _nEvt ;
  StringVec _detailled;

} ;


//Set of geometric parameters for initialization of a CEDGeoTube class object
struct CEDGeoTubeParams {
  double Rmax; double Rmin; double inner_symmetry; double outer_symmetry; double phi0; double delta_phi; double delta_z; double z0; 
  //boolean that decides if the GeoTube is drawn twice at two different zPositions
  bool isBarrel;
};
//Set of geometric parameters for initialization of a CEDGeoBox class object
struct CEDGeoBox {
  double  sizes[3] ;
  double  center[3] ;
  double rotate[3];
};
//Convenient summary of both parameter sets above as (tracker) layers may be drawn as one tube or as a sequence of staves (-->GeoBox)
struct LayerGeometry {
  CEDGeoTubeParams tube;
  std::vector<CEDGeoBox> staves;
};

#endif



