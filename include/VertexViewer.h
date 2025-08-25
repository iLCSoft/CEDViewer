#ifndef VERTEXVIEWER_H
#define VERTEXVIEWER_H 1

#include "marlin/Processor.h"
#include "EVENT/MCParticle.h"
#include "lcio.h"
#include <string>
#include <vector>

using namespace lcio ;
using namespace marlin ;


/** Vertex Viewer Processor <br>
 *  @author A.Raspereza, DESY
 *  @version $Id$ 
 */
class VertexViewer : public Processor {
  
 public:
  
  virtual Processor*  newProcessor() { return new VertexViewer ; }
  
  
  VertexViewer() ;
  
  virtual void init() ;
  
  virtual void processRunHeader( LCRunHeader* run ) ;
  
  virtual void processEvent( LCEvent * evt ) ; 
  
  virtual void check( LCEvent * evt ) ; 
  
  virtual void end() ;
  
  
 protected:

  int _nRun{-1} ;
  int _nEvt{0} ;
  
  std::vector<std::string> _trackerHitCollection{};
  std::string _trueTracksCollection{};
  std::string _tracksCollection{};
  std::string _trueTracksMCPCollection{};
  std::string _tracksMCPCollection{};
  std::vector<std::string> _simTrackerHitCollection{};
  
  int _layerTrackerHits{-1};
  int _layerTrueTracks{-1};
  int _layerTracks{-1};
  int _layerSimTrackerHits{-1};
  int returnColor(int counter);
  float _cutOnD0{1e+20f}, _cutOnZ0{1e+20f};
  float _bField{};

  // int _nTPCCut;
  int _detModel{0};

} ;

#endif



