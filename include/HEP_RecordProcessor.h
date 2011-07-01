#ifndef HEP_RecordProcessor_h
#define HEP_RecordProcessor_h 1

#include "marlin/Processor.h"
#include "lcio.h"

using namespace lcio ;
using namespace marlin ;


namespace marlin {
  /**  Displays the MCParticles in CED. <br>
   *  Color scheme:<br>
   *  neutral hadrons:  pink 'clusters' + magenta 'tracks' <br>
   *  photons: yellow 'clusters' <br>
   *  charged tracks:  green <br>
   *  electrons:  magenta-pink  <br>        
   *  muons:  light green <br>
   * <p> Toggle display of neutrals with [1] and charged Tracks with [Shift-1].
   * <p> F.Gaede: changed to use MarlinCED
   * 
   * @author V.Morgunov, DESY
   * @version $Id$
   */
  class HEP_RecordProcessor : public Processor {

  public:
    virtual Processor*  newProcessor() { return new HEP_RecordProcessor ; }
    HEP_RecordProcessor() ;
    virtual void init() ; //the begin of the job
    virtual void processRunHeader( LCRunHeader* run ) ; //  every run.
    virtual void processEvent( LCEvent * evt ) ;   // very event
    virtual void check( LCEvent * evt ) ; 
    virtual void end() ; //  after data processing 
  protected:
    int _nRun ;
    int _nEvt ;
  } ;

}// namespace marlin
#endif
