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
 *                           0 - the DSTViewer processor does not wait after 
 *                               the display of the current event <br>
 *                           1 - the DSTViewer processor waits for a keyboard 
 *                               (std) input after each event. This is the default.
 *  <br>
 *
 *  @authors A.Raspereza (DESY), O.Wendt (DESY)
 *  @version $Id: DSTViewer.h,v 1.1 2008-09-09 14:38:53 darasz Exp $ 
 */
class DSTViewer : public Processor {
  
 public:
  
  virtual Processor*  newProcessor() { return new DSTViewer ; }
  
  DSTViewer() ;
  
  virtual void init() ;
  
  virtual void processRunHeader( LCRunHeader* run ) ;
  
  virtual void processEvent( LCEvent * evt ) ; 
  
  virtual void check( LCEvent * evt ) ; 
  
  virtual void end() ;
  
 protected:

  int _nRun ;
  int _nEvt ;
  // bool dispHelix;
  
//  std::vector<std::string> _caloHitCollections;
//  std::vector<std::string> _simCaloHitCollections;
//  std::vector<std::string> _trackerHitCollections;
//  std::vector<std::string> _simTrackerHitCollections;
//  std::string _trueClustersCollection;
//  std::string _trueTracksCollection;
//  std::string _clustersCollection;
//  std::string _tracksCollection;
  std::string _particleCollection;
  
  StringVec _jetCollections;

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
  
  int returnTrackColor(int type);
  
  /**
   * @author: S.Daraszewicz (UoE)
   * @date: 21.08.08
   * 
   * @param eneCluster : Energy of the cluster
   * @param cutoff_min : Min cut for the cluster energy
   * @param cutoff_max : Max cut for the cluster energy
   * 
   * Returns the colour for the visualisation of a cluster (in log scale). Colours run from red (0x660000) to blue (0x6600FF) 
   * for high and low energy particles respectively.
   */
  int returnClusterColor(float eneCluster, float cutoff_min, float cutoff_max);

  /**
   * @author: S.Daraszewicz (UoE)
   * @date: 22.08.08
   * 
   * @param eneCluster : Energy of the cluster
   * @param cutoff_min : Min cut for the cluster energy
   * @param cutoff_max : Max cut for the cluster energy
   * @param color_steps: Number of colours to be used in the scale
   * @param scale	   : Linear or Logarithmic colour scale
   * @param colorMap   : Color map type to be used
   * 
   * Returns the size for the visualisation of a cluster (in log scale) within user defined sizes.
   */
   
    void showLegendSpectrum(const unsigned int color_steps, char scale, int colorMap, float ene_min, float ene_max, unsigned int ticks);

  	int returnRGBClusterColor(float eneCluster, float cutoff_min, float cutoff_max, int color_steps, char scale, int colorMap);

  	int returnClusterSize(float eneCluster, float cutoff_min, float cutoff_max);
    
  	int returnTrackSize(float type);
  	
  	int returnJetLayer(std::string jetColName);
  	
  	int returnJetColor(std::string jetColName, int colNumber);
  	
  	float * returnConeColor(std::string jetColName);

  	float _bField;

} ;

#endif



