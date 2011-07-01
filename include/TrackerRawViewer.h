#ifndef TrackerRawViewer_h
#define TrackerRawViewer_h 1

#include "marlin/Processor.h"
#include "lcio.h"
#include <string>
#include <vector>
#include <map>

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
 *  Draws all ADC channels with a color code on layer 1. Optionally, hits and tracks are displayed
 *  on layers 2 and 3. On Layer 4 a color coded map of the mean number of integrated ADC counts 
 *  per pad is drawn with a color code. 
 *  
 *  <h4>Input </h4>
 *  A collection of TrackerRawData (and optionally collections of TrackerHits and Tracks).
 *
 * 
 * @param CollectionName:            Name of the TrackerRawData collection (default: AltroRawData )
 * @param ChannelMappingCollection:  Name of the LCCD collection with channel mapping (default: ADCChannelMapping)
 * @param ChannelPositionTextFile:   Optionally use a text file for the hardware channel to position mapping - overwrites mapping from LCCD and GEAR 
 * @param DriftVelocity:             Drift velocity in micron/ns (default: 80.)
 * @param ColorScaleMaxADC:          ADC value used for the maximum of the color scales (default: 32)
 * @param WaitForKeyBoard:           true: display one event and wait for 'return' - false: continiously display events
 * @param ColorScheme:               Color scheme - 1: hot , 2 : cold    (default: 1)
 * @param HitCollectionName:         Name of the TrackerHit collection   (optional, default: TPCTrackerHits)
 * @param TrackCollectionName:       Name of the Track collection        (optional, default: TPCTracks)
 *   
 *   
 * 
 *  @author F. Gaede, DESY
 *  @version $Id$
 */

class TrackerRawViewer : public Processor {
  
  /** Helper struct for running mean of ADC counts on pads */
  struct PadMeanADC{
    
    PadMeanADC() :Total(0.),N(0){}
    
    float Total ;
    int N ;
    
    void operator+=(float adc){
      Total += adc ;
      //      ++N ; 
    }
  };
  
  typedef std::map< std::pair<int, int>,  PadMeanADC > ADCMap ; 
  
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

  /** Input collection name.
   */
  std::string _colName ;
  std::string _hitColName ;
  std::string _trkColName ;
  std::string _chMapCollection ;
  std::string _chPosTextFile ;

  const gear::TPCParameters* _tpcParams  ;
  ChannelPosMap* _posMap ;
  ChannelMap* _chMap ;

  float _driftVelocity ;
  int _colorScaleMaxADC ;
  int _colorScheme ;
  bool _waitForKeyboard ;

  std::vector<int> _colors ;

  int _nRun ;
  int _nEvt ;

  ADCMap _adcMap ;
} ;

#endif



