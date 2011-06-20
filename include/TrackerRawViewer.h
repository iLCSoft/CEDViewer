#ifndef TrackerRawViewer_h
#define TrackerRawViewer_h 1

#include "marlin/Processor.h"
#include "lcio.h"
#include <string>

// MarlinTPC
#define USE_LCCD
#include "ADCChannelMapping.h"

//---- LCCD ----
#include "lccd/ConditionsMap.hh"
typedef lccd::ConditionsMap< tpcconddata::ADCChannelMapping::key_type , tpcconddata::ADCChannelMapping > ChannelMap ;

// -- Gear
#include "gear/GEAR.h"
#include "gear/TPCParameters.h"

using namespace lcio ;
using namespace marlin ;

class ChannelPosMap ;

/** Simple event display that visualizes TrackerRawData for the LCTPC large prototype with CED.
 *  
 *  <h4>Input </h4>
 *  A collection of TrackerRawData (and optionally a collection of TrackerHits).
 *
 *  <h4>Output</h4> 
 *  collection of Tracks 
 * 
 * @param CollectionName:            Name of the TrackerRawData collection (default: AltroRawData )
 * @param ChannelMappingCollection:  Name of the LCCD collection with channel mapping (default: ADCChannelMapping)
 * @param ChannelPositionTextFile:   Optionally use a text file for the hardware channel to position mapping - overwrites mapping from LCCD and GEAR 
 * @param DriftVelocity:             Drift velocity in micron/ns (default: 80.)
 * @param ColorScaleMaxADC:          ADC value used for the maximum of the color scale (default: 32)
 * @param ColorScheme:               Color scheme - 1: hot , 2 : cold    (default: 1)
 * @param HitCollectionName:         Name of the TrackerHit collection   (default: TPCTrackerHits)
 *   
 * 
 *  @author F. Gaede, DESY
 *  @version $Id:$
 */

class TrackerRawViewer : public Processor {
  
 public:
  
  virtual Processor*  newProcessor() { return new TrackerRawViewer ; }
  
  
  TrackerRawViewer() ;
  
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

  /** Computes the pad center for given hardware ID in padCenter - returns false if nothing found */
  bool getPadCenter( double *padCenter, int cellID0, int cellID1 ) ;


  /** Input collection name.
   */
  std::string _colName ;
  std::string _hitColName ;
  std::string _chMapCollection ;
  std::string _chPosTextFile ;

  const gear::TPCParameters* _tpcParams  ;
  ChannelPosMap* _posMap ;
  ChannelMap* _chMap ;

  float _driftVelocity ;
  int _colorScaleMaxADC ;
  int _colorScheme ;

  int _nRun ;
  int _nEvt ;
} ;

#endif



