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
 *  @version $Id$ 
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
  
  /** Input collection name.
   */
  StringVec _drawCollections {};
  StringVec _drawCollectionsLayer {};

  bool      _usingParticleGun            = false ;
  int       _drawHelixForTracks          = 0 ;
  bool      _colorEnergy                 = false ;
  double    _colorEnergyMin              = 0.0 ;
  double    _colorEnergyMax              = 35.0 ;
  double    _colorEnergySaturation       = 0.8 ;
  double    _colorEnergyValue            = 0.8 ;
  bool      _colorEnergyAuto             = false ;
  double    _scaleLineThickness          = 1 ;
  double    _scaleMarkerSize             = 1 ;
  int       _drawDetectorID              = 0  ;
  int       _colorScheme                 = 10;
  float     _mcpECut                     = 0.001;

  float     _helix_max_r                 = 2000 ;
  float     _helix_max_z                 = 2500 ;
  bool      _useTPCForLimitsOfHelix      = false ;
  int       _waitForKeyboard             = 1 ;
  int       _drawHelixForPFOs            = 0 ;
  int       _useColorForHelixTracks      = 0;
  IntVec    _colors {};


  int _nRun = 0 ;
  int _nEvt = 0 ;
} ;

#endif



