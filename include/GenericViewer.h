#ifndef VIEWERAR_H
#define VIEWERAR_H 1

#include "marlin/Processor.h"
#include "EVENT/MCParticle.h"
#include "lcio.h"
#include <string>
#include <vector>

using namespace lcio ;
using namespace marlin ;


/** Viewer Processor <br>
 *  This processor displays collections of LCIO objects, <br>
 *  including SimCalorimeterHits, SimTrackerHits, <br>
 *  CalorimeterHits, TrackerHits, reconstructed Tracks, <br>
 *  reconstructed Clusters and reconstructed particle flow objects. <br>
 *  Display of a given collection is manipulated with numeric keys <br>
 *  associated with toggle layers. <br>
 *  Processor also draws wire-frame sketch of the detector. <br>
 *  A current event is displayed in the window panel activated by C Event <br>
 *  Display (CED) program glced. You have to launch this program before <br>
 *  starting Marlin. <br>
 *  The user should specify the name of collections to be displayed <br>
 *  via processor parameters described below. <br>
 *  @param SimCaloHitCollections : names of SimCalorimeterHit 
 *  collections to be displayed <br>
 *  @param LayerSimCaloHit : toggle layer for SimCalorimeterHits <br>
 *  @param SimTrackerHitCollections : SimTrackerHit collections
 *  to be displayed <br>
 *  @param LayerSimTrackerHit : toggle layer for SimTrackerHits <br>
 *  @param CaloHitCollections : names of CalorimeterHit 
 *  collections to be displayed <br>
 *  @param LayerCaloHit : toggle layer for CalorimeterHits <br>
 *  @param TrackerHitCollections : names of TrackerHit collections
 *  to be displayed <br>
 *  @param LayerTrackerHit : toggle layer for TrackerHits <br>
 *  @param TrueTrackCollection : True Monte Carlo track collection name  
 *  to be displayed <br>
 *  @param LayerTrueTracks : toggle layer for True Monte Carlo tracks <br>
 *  @param TrueClusterCollection : name of true MC cluster collection
 *  to be displayed <br>
 *  @param LayerTrueClusters : toggle layer for true MC clusters <br>
 *  @param TrackCollection : name of reconstructed track collection 
 *  collections to be displayed <br>
 *  @param LayerTracks : toggle layer for reconstructed tracks <br>
 *  @param ClusterCollection : name of reconstructed Cluster collection
 *  to be displayed <br>
 *  @param LayerClusters : toggle layer for reconstructed clusters <br>
 *  @param ParticleCollection : name of reconstructed Particle collection
 *  to be displayed <br>
 *  @param LayerReco : toggle layer for reconstructed Particles <br>
 *  <br>
 *  If a layer toggle is set to -1, the corresponding collection is not displayed <br>
 *  <br>
 *  @param DetectorModel : detector model <br>
 *   0 - LDC model (Mokka program, models D09, D10, D12) <br>
 *   1 - SiD model (SLIC program) <br>
 *   2 - GLD model (Jupiter program) <br>
 *  @param MagneticField : Magnetic field in Tesla <br>
 *  @param WaitForKeyboard : toggle to wait for keyboard input or not <br>
 *                           0 - the GenericViewer processor does not wait after 
 *                               the display of the current event <br>
 *                           1 - the GenericViewer processor waits for a keyboard 
 *                               (std) input after each event. This is the default.
 *  <br>
 *
 *  @authors A.Raspereza (DESY), O.Wendt (DESY)
 *  @version $Id$ 
 */
class GenericViewer : public Processor {
  
 public:
  
  virtual Processor*  newProcessor() { return new GenericViewer ; }
  
  
  GenericViewer() ;
  
  virtual void init() ;
  
  virtual void processRunHeader( LCRunHeader* run ) ;
  
  virtual void processEvent( LCEvent * evt ) ; 
  
  virtual void check( LCEvent * evt ) ; 
  
  virtual void end() ;
  
  
 protected:

  int _nRun ;
  int _nEvt ;
  
  std::vector<std::string> _caloHitCollections;
  std::vector<std::string> _simCaloHitCollections;
  std::vector<std::string> _trackerHitCollections;
  std::vector<std::string> _simTrackerHitCollections;
  std::string _trueClustersCollection;
  std::string _trueTracksCollection;
  std::string _clustersCollection;
  std::string _tracksCollection;
  std::string _particleCollection;

  int _layerCaloHit;
  int _layerSimCaloHit;
  int _layerTrackerHit;
  int _layerSimTrackerHit;
  int _layerTrueClusters;
  int _layerTrueTracks;
  int _layerClusters;
  int _layerTracks;
  int _layerMCP;
  int _layerBosons;
  int _layerReco;

  int _detModel;

  int _waitForKeyboard;

  std::map<MCParticle *, int > _mcpList;
  
  int returnColor(int counter);

  float _bField;


} ;

#endif



