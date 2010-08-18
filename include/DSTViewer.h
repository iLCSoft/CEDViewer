#ifndef VIEWERAR_H
#define VIEWERAR_H 1

#include "marlin/Processor.h"
#include "EVENT/MCParticle.h"
#include "lcio.h"
#include <string>
#include <vector>

using namespace lcio ;
using namespace marlin ;


/** DSTViewer <br>
 *  This processor displays .... <br>
 *  @param ParticleCollection : name of reconstructed Particle collection to be displayed <br>
     .....
 *  <br>
 *  @authors Szymon Daraszewicz (DESY/UOE)
 *  @version $Id: DSTViewer.h,v 1.3 2008-09-15 10:13:48 darasz Exp $ 
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
  	
  	int returnIpLayer(std::string jetColName);
  	
  	int returnJetColor(std::string jetColName, int colNumber);
  	
  	int addAlphaChannelToColor(int color, int alphaChannel);
  	
  	float * returnConeColor(std::string jetColName);

  	float _bField;
    
    void writeLayerDescription(void);

} ;

#endif



