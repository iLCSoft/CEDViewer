#ifndef CEDViewer_h
#define CEDViewer_h 1

#include "marlin/Processor.h"

#include "lcio.h"
#include <string>


using namespace lcio ;
using namespace marlin ;

/** LCIO collection based viewer for CED event dislplay (A. Zhelezov).
 *  Define the collection name, marker, size and optionally the layer number in the steering file.
 *  For 'DrawCollection' a default layer number in CED is chosen.
 * 
 *  @param DrawCollection - collection to be displayed ( ColName, marker type[0-2] )
 *  @param DrawInLayer    - collection to be displayed ( ColName, marker type[0-2], size) 
 *
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
  void printParticle(int id, LCEvent * evt);

  
  
 protected:


  /** Input collection name.
   */
  StringVec _drawCollections ;
  StringVec _drawCollectionsLayer ;
  
  IntVec _colors ;

  int _nRun ;
  int _nEvt ;
} ;

#endif



