#ifndef CEDViewer_h
#define CEDViewer_h 1

#include "marlin/Processor.h"

#include "lcio.h"
#include <string>


using namespace lcio ;
using namespace marlin ;

/** Cluster viewer based on CED by A. Zhelezov.
 *  @author F.Gaede, DESY
 *  @version $Id: CEDViewer.h,v 1.2 2005-08-04 13:45:21 gaede Exp $ 
 */

class CEDViewer : public Processor {
  
 public:
  
  virtual Processor*  newProcessor() { return new CEDViewer ; }
  
  
  CEDViewer() ;
  
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
  
  
 protected:


  /** Input collection name.
   */
  StringVec _drawCollections ;
  
  IntVec _colors ;

  int _nRun ;
  int _nEvt ;
} ;

#endif



