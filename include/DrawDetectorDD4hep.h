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

  static void drawDD4hepDetector( DD4hep::Geometry::LCDD& lcdd) ;
  
  
 protected:
  bool _begin;
  int _nRun ;
  int _nEvt ;
} ;

#endif



