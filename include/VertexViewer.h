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

  int _nRun ;
  int _nEvt ;
  
  std::vector<std::string> _trackerHitCollection;
  std::string _trueTracksCollection;
  std::string _tracksCollection;
  std::string _trueTracksMCPCollection;
  std::string _tracksMCPCollection;
  std::vector<std::string> _simTrackerHitCollection;
  
  int _layerTrackerHits;
  int _layerTrueTracks;
  int _layerTracks;
  int _layerSimTrackerHits;
  int returnColor(int counter);
  float _cutOnD0,_cutOnZ0;
  float _bField;

  int _nTPCCut;
  int _detModel;

} ;

#endif



