/* -*- Mode: C++; tab-width: 2; indent-tabs-mode: nil; c-basic-offset: 2 -*- */
#include "CEDViewer.h"
#include <iostream>

#include <EVENT/LCCollection.h>
#include <EVENT/Cluster.h>
#include <EVENT/Track.h>
#include <EVENT/MCParticle.h>
#include <EVENT/SimCalorimeterHit.h>
#include <EVENT/SimTrackerHit.h>

#include <EVENT/ReconstructedParticle.h>

#include <UTIL/LCTypedVector.h>
#include <UTIL/CellIDDecoder.h>
#include <UTIL/ILDConf.h>

#include "MarlinCED.h"

#include <math.h>
#include <cmath>
#include <cstdlib>
#include <algorithm>

#include <marlin/Global.h>
#include <gear/GEAR.h>
#include <gear/BField.h>
#include <gear/TPCParameters.h>
#include <gear/TPCModule.h>
#include <gear/CalorimeterParameters.h>

using namespace lcio ;
using namespace marlin ;

#define SIMTRACKERHIT_LAYER 1
#define TRACKERHIT_LAYER 2
#define TRACK_LAYER 3
#define TRACKHELIX_LAYER 3

#define SIMCALORIMETERHIT_LAYER 4
#define CALORIMETERHIT_LAYER 5
#define CLUSTER_LAYER 6

#define MCPARTICLE_LAYER 7
#define RECOPARTICLE_LAYER 8


CEDViewer aCEDViewer ;


/**Helper struct for drawing collections*/
struct DrawParameters{
    DrawParameters(const std::string& colName, int size, int marker, int layer ) :
    ColName( colName ),
    Size( size ),
    Marker( marker ),
    Layer( layer ) {
    }
    std::string ColName ;
    int Size ;
    int Marker ;
    int Layer ;
};



CEDViewer::CEDViewer() : Processor("CEDViewer") {
    
    // modify processor description
    _description = "CEDViewer: event display of LCIO objects  - based on CED by A.Zhelezov." ;
    
    StringVec cluExample ;
    cluExample.push_back("CalorimeterHits" ) ;
    cluExample.push_back("0" ) ;
    cluExample.push_back("3" ) ;
    
    registerOptionalParameter( "DrawCollection" ,
                              "collection to be displayed ( ColName, marker type[0-2], size)"  ,
                              _drawCollections ,
                              cluExample ,
                              cluExample.size() ) ;
    
    StringVec layerExample ;
    layerExample.push_back("TPCTracks" ) ;
    layerExample.push_back("0" ) ;
    layerExample.push_back("3" ) ;
    layerExample.push_back("4" ) ;
    
    registerOptionalParameter( "DrawInLayer" ,
                              "collection to be displayed ( ColName, marker type[0-2], size, layer)",
                              _drawCollectionsLayer ,
                              layerExample ,
                              layerExample.size() ) ;
    
    
    registerProcessorParameter( "DrawHelixForTrack" ,
                               "draw a helix for Track objects: -1: none, 0: default, 1: atIP, 2: atFirstHit, 3: atLastHit, 4: atCalorimeter",
                               _drawHelixForTracks ,
                               0  ) ;
    
    registerProcessorParameter( "DrawDetectorID" ,
                               "draw detector from GEAR file with given ID (see MarlinCED::newEvent() ) : 0 ILD, 1 CLIC, -1 none",
                               _drawDetectorID ,
                               0 ) ;
    
    registerProcessorParameter( "ColorScheme" ,
                               "color scheme to be used for drawing - see startup log MESSAGEs for options",
                               _colorScheme,
                               10 ) ;
    
    registerProcessorParameter( "MCParticleEnergyCut" ,
                               "minimum energy of MCParticles to be drawn",
                               _mcpECut ,
                               float(0.001) ) ;
    
    
    registerProcessorParameter( "UsingParticleGun" ,
                               "If set true generator status is ignored for MCParticles to be drawn",
                               _usingParticleGun ,
                               bool(false) ) ;
    
    registerProcessorParameter( "UseTPCForLimitsOfHelix" ,
                               "Use the gear parameters to define the max extent of drawing a helix",
                               _useTPCForLimitsOfHelix ,
                               bool(true) ) ;
    
    registerProcessorParameter( "HelixMaxR" ,
                               "Max R (mm) Extent for drawing Helix if UseTPCForLimitsOfHelix false",
                               _helix_max_r ,
                               float(2000.0) ) ;
    
    registerProcessorParameter( "HelixMaxZ" ,
                               "Max Z (mm) Extent for drawing Helix if UseTPCForLimitsOfHelix false",
                               _helix_max_z ,
                               float(2500.0) ) ;
    
    registerProcessorParameter("WaitForKeyboard",
                               "Wait for Keyboard before proceed",
                               _waitForKeyboard,
                               (int)1);
    
    registerProcessorParameter("DrawHelixForPFOs",
                               "draw a helix for PFO objects (usefull for Si tracker): 0: disabled, 1: enabled",
                               _drawHelixForPFOs,
                               (int)0);

    registerProcessorParameter("UseColorForHelixTracks",
                               "draw helices in the color of the track/PFO: 0: disabled (lightgrey), 1: enabled",
                               _useColorForHelixTracks,
                               (int)0);
 
}

void CEDViewer::init() {
    
    // usually a good idea to
    printParameters() ;
    MarlinCED::init(this) ;
    
    
    static unsigned bojum_colors[ncol*nscheme]={
        0x510000,0x660202,0x720202,0x840202,0x960303,0xa80303,0xbc0303,
        0xce0404,0xe00404,0xef0404,0xff0000,0xfc0505,0xf91111,0xf92222,
        0xfc4949,0xf97777,0xf98e8e,0xf9aeae,0xf7b2b2,0xf7cdcd,
        
        0x512401,0x662d03,0x7a3602,0x934204,0xa54a04,0xb75103,0xc65803,
        0xd35e04,0xe56604,0xff6e00,0xf4710c,0xf4791a,0xf9842a,0xf98f3e,
        0xf99d57,0xf9a768,0xf9b47f,0xf9bf93,0xf9c9a4,0xf9d2b3,
        
        0x48014c,0x610266,0x700277,0x93039b,0xb103ba,0xc904d3,0xda04e5,
        0xe604f2,0xfa00ff,0xf00ffc,0xec1df7,0xef2ff9,0xeb42f4,0xec58f4,
        0xed6bf4,0xf486f9,0xf59af9,0xf8b5fc,0xf4c3f7,0xf8d6f9,
        
        0x2d0251,0x3d026d,0x4a0284,0x550399,0x5f03aa,0x6903bc,0x7102cc,
        0x8004e5,0x9800ff,0x8e0ef7,0x9922f9,0xa134f9,0xa845f9,0xb057f9,
        0xbc70f9,0xbf77f9,0xc98ef9,0xd3a4f9,0xddbbf9,0xecd6ff,
        
        0x00004f,0x020268,0x03037f,0x030399,0x0303b2,0x0404cc,0x0404e0,
        0x0404ef,0x0000ff,0x0c0cf9,0x1b1bf9,0x2a2af9,0x3939f9,0x4d4df9,
        0x6363f9,0x7272f9,0x8484f9,0x9898f9,0xb3b3f9,0xcacaf9,
        
        0x01434c,0x025c68,0x027382,0x04899b,0x0397aa,0x03abc1,0x04b9d1,
        0x04cbe5,0x00d8ff,0x0cdef9,0x1bddf7,0x2ae1f9,0x3ddff4,0x59e4f7,
        0x6be9f9,0x7cebf9,0x94f0fc,0xa3edf7,0xb6f2f9,0xc7f4f9,
        
        0x014443,0x01605f,0x027f7d,0x039996,0x03b2af,0x01c6c3,0x04ddda,
        0x04efeb,0x00ffff,0x09f9f5,0x0ef9f5,0x20f9f6,0x32fcf9,0x3ef9f6,
        0x52f9f7,0x6bf9f7,0x7ff9f7,0x95f9f8,0xb1f9f8,0xcaf9f9,
        
        0x016001,0x027002,0x027f02,0x029102,0x05aa05,0x05bf05,0x06d306,
        0x04e504,0x00ff00,0x09f909,0x18f918,0x2cf92c,0x43f943,0x52f952,
        0x63f963,0x77f977,0x8bf98b,0x9ff99f,0xb3f9b3,0xcaf9ca,
        
        0x344701,0x4b6603,0x608202,0x739b04,0x83b203,0x96cc04,0xa7e204,
        0xb1ef04,0xb6ff00,0xbaf713,0xbff725,0xc5f73b,0xcbf751,0xd3f968,
        0xd7f97a,0xd8f48b,0xe2f9a2,0xe1f4ad,0xe7f7bb,0xe9f4c8,
        
        0x565501,0x727002,0x898702,0xa5a303,0xb7b403,0xd1cd04,0xe2df04,
        0xefeb04,0xffff00,0xf9f509,0xf9f618,0xf9f62a,0xf7f43b,0xf9f64a,
        0xf9f759,0xf9f76b,0xf9f77c,0xf9f88b,0xf9f89f,0xfcfbbd
    };
    
    
    
    switch (_colorScheme) {
            
        case Red         :
        case Orange      :
        case Plum        :
        case Violet      :
        case Blue        :
        case LightBlue   :
        case Aquamarine  :
        case Green       :
        case Olive       :
        case Yellow      :
            
            for(int i = 0 ; i<ncol ; ++i)
                _colors.push_back( bojum_colors[  ( _colorScheme * ncol  ) + i ] ) ;
            
            break;
            
        case Dark:
            
            for(int i = 0 ; i<ncol ; ++i )
                for(int j = 0 ; j<nscheme ; ++j )
                    _colors.push_back( bojum_colors[  ( j * ncol  ) + i ] ) ;
            break ;
            
        case Light:
            
            for(int i = ncol-1 ; i>=0 ; --i )
                for(int j = 0 ; j<nscheme ; ++j )
                    _colors.push_back( bojum_colors[  ( j * ncol  ) + i ] ) ;
            break ;
            
        case Classic:
        default:
            _colors.push_back( 0xff00ff );
            _colors.push_back( 0x00ff00 );
            _colors.push_back( 0xffff00 );
            _colors.push_back( 0x0000ff );
            _colors.push_back( 0xff00ff );
            _colors.push_back( 0x00ffff );
            _colors.push_back( 0xffffff );
            
            _colors.push_back( 0xff88ff );
            _colors.push_back( 0x88ff88 );
            _colors.push_back( 0xffff88 );
            _colors.push_back( 0x8888ff );
            _colors.push_back( 0xff88ff );
            _colors.push_back( 0x88ffff );
            _colors.push_back( 0xffffff );
            break ;
            
    }
    
    streamlog_out(MESSAGE) << " ---- selected color scheme:  "  << _colorScheme << "\n"
    << "    possible options: "  <<   "\n"
    << "      Red          "     << Red          <<   "\n"
    << "      Orange       "     << Orange       <<   "\n"
    << "      Plum         "     << Plum         <<   "\n"
    << "      Violet       "     << Violet       <<   "\n"
    << "      Blue         "     << Blue         <<   "\n"
    << "      LightBlue    "     << LightBlue    <<   "\n"
    << "      Aquamarine   "     << Aquamarine   <<   "\n"
    << "      Green        "     << Green        <<   "\n"
    << "      Olive        "     << Olive        <<   "\n"
    << "      Yellow       "     << Yellow       <<   "\n"
    << "      Dark         "     << Dark         <<   "\n"
    << "      Light        "     << Light        <<   "\n"
    << "      Classic      "     << Classic      <<   "\n"
    << " -------- " << std::endl ;
    
    
    _nRun = 0 ;
    _nEvt = 0 ;
    
}

void CEDViewer::processRunHeader( LCRunHeader* run) {
    _nRun++ ;
}


void CEDViewer::processEvent( LCEvent * evt ) {
    
    //-----------------------------------------------------------------------
    // Reset drawing buffer and START drawing collection

    MarlinCED::newEvent(this , _drawDetectorID ) ;
    
    CEDPickingHandler &pHandler=CEDPickingHandler::getInstance();
    
    pHandler.update(evt);
    
    //-----------------------------------------------------------------------
    
    if( _useTPCForLimitsOfHelix ){
        try{
            // Support for multi-module TPCs
            const gear::TPCParameters& gearTPC = Global::GEAR->getTPCParameters() ;
            
            _helix_max_r = gearTPC.getPlaneExtent()[1] + 300.0;
            _helix_max_z = gearTPC.getMaxDriftLength()+600.0;
        }
        catch(gear::UnknownParameterException& e) {}
    }
    
    //add DrawInLayer to
    std::vector< DrawParameters > drawParameters ;
    drawParameters.reserve( _drawCollections.size() + _drawCollectionsLayer.size()  ) ;
    
    if( parameterSet( "DrawCollection" ) ) {
        
        unsigned index = 0 ;
        while( index < _drawCollections.size() ){
            
            const std::string & colName = _drawCollections[ index++ ] ;
            int marker = std::atoi( _drawCollections[ index++ ].c_str() ) ;
            int size = std::atoi( _drawCollections[ index++ ].c_str() ) ;
            //std::cout << "#######################      SIZE: " << size << " #######################################" << std::endl;
            int layer = -1 ;
            
            drawParameters.push_back(DrawParameters( colName,size,marker,layer ) );
        }
    }
    if( parameterSet( "DrawInLayer" ) ) {
        
        unsigned index = 0 ;
        while( index < _drawCollectionsLayer.size() ){
            
            const std::string & colName = _drawCollectionsLayer[ index++ ] ;
            int marker = std::atoi( _drawCollectionsLayer[ index++ ].c_str() ) ;
            int size   = std::atoi( _drawCollectionsLayer[ index++ ].c_str() ) ;
            int layer  = std::atoi( _drawCollectionsLayer[ index++ ].c_str() ) ;
            
            drawParameters.push_back(DrawParameters( colName,size,marker,layer ) );
            //std::cout << "layer: " << layer << " description: " << colName << std::endl; //hauke
        }
    }
    
    unsigned nCols = drawParameters.size() ;
    for(unsigned np=0 ; np < nCols ; ++np){
        
        const std::string & colName = drawParameters[np].ColName ;
        int size =   drawParameters[np].Size ;
        int marker = drawParameters[np].Marker ;
        int layer =  drawParameters[np].Layer ;
        
        LCCollection* col = 0 ;
        try{
            
            col = evt->getCollection( colName ) ;
            
        }catch(DataNotAvailableException &e){
            
            streamlog_out( WARNING ) << " collection " << colName <<  " not found ! "   << std::endl ;
            continue ;
            
        }
        
        if( col->getTypeName() == LCIO::CLUSTER ){
            
            // find Emin and Emax of cluster collection for drawing
            float emin=1.e99, emax=0. ;
            for( int i=0 ; i< col->getNumberOfElements() ; i++ ){
                Cluster* clu = dynamic_cast<Cluster*>( col->getElementAt(i) ) ;
                float e = clu->getEnergy() ;
                if( e > emax) emax = e  ;
                else if( e < emin )  emin = e ;
            }
            streamlog_out( DEBUG )  << " nClu: " << col->getNumberOfElements()
            << ", Emin: " << emin << " GeV , Emax: " << emax << " GeV " << std::endl ;
            
            for( int i=0 ; i< col->getNumberOfElements() ; i++ ){
                
                Cluster* clu = dynamic_cast<Cluster*>( col->getElementAt(i) ) ;
                const CalorimeterHitVec& hits = clu->getCalorimeterHits() ;
                
                layer = ( layer > -1 ? layer : CLUSTER_LAYER ) ;
                drawParameters[np].Layer = layer ;
                
                
                //ced_describe_layer( colName.c_str() ,layer);
                
                
                int color =  _colors[ i % _colors.size() ] ;
                for( CalorimeterHitVec::const_iterator it = hits.begin();  it != hits.end() ; it++ ) {
                    
                    //hauke hoelbe: add id for picking!
                    ced_hit_ID( (*it)->getPosition()[0],
                               (*it)->getPosition()[1],
                               (*it)->getPosition()[2],
                               marker, layer, size , color, clu->id() ) ;
                    
                } // hits
                
                float x = clu->getPosition()[0] ;
                float y = clu->getPosition()[1] ;
                float z = clu->getPosition()[2] ;
                
                ced_hit_ID( x,y,z, marker, layer , size*3 , color, clu->id() ) ;
                
                LCVector3D v(x,y,z) ;
                
                LCVector3D d(1.,1.,1.) ;
                
                
                float length = 100 + 500 * ( clu->getEnergy() - emin )  / ( emax - emin )  ;
                
                // 	  const FloatVec& pv = clu->getShape() ;
                // 	  if( pv.size() > 7 ){
                // 	      length =  3 * clu->getShape()[3] ;   //FIXME : this is not generic !!!!
                // 	  } else {
                //       std::cout <<   "  ---- no cluster shape[3] :( " << std::endl ;
                //     }
                
                d.setMag( length  ) ;
                
                d.setTheta( clu->getITheta() ) ;
                d.setPhi( clu->getIPhi() ) ;
                
                LCVector3D dp( v + d ) , dm( v - d )   ;
                
                ced_line_ID( dp.x() , dp.y() , dp.z(),
                            dm.x() , dm.y() , dm.z(),
                            layer  , 1 , color, clu->id() );
                
                
            } // cluster
            
        } else if( col->getTypeName() == LCIO::TRACK ){
            
            for( int i=0 ; i< col->getNumberOfElements() ; i++ ){
                
                Track* trk = dynamic_cast<Track*>( col->getElementAt(i) ) ;
                
                // -- collect hits from all track segments
                TrackerHitVec tHV ;
                if( !trk->getTracks().empty() ){
                    
                    streamlog_out( DEBUG ) << " -- track has "<< trk->getTracks().size() << " subtracks - will use these for displaying hits "
                    << std::endl ;
                    
                    std::copy( trk->getTrackerHits().begin() , trk->getTrackerHits().end() , std::back_inserter(  tHV ) ) ;
                    
                    for( unsigned j=0 ,N = trk->getTracks().size() ; j<N ; ++j ){
                        
                        const Track* t = trk->getTracks()[j] ;
                        
                        streamlog_out( DEBUG ) << " -- track j= "<< j << " has " <<  t->getTrackerHits().size() << " hits  "
                        << std::endl ;
                        
                        std::copy( t->getTrackerHits().begin() , t->getTrackerHits().end() , std::back_inserter(  tHV ) ) ;
                    }
                    
                }
                const TrackerHitVec& hits = (  tHV.empty() ?  trk->getTrackerHits()  :  tHV  ) ;
                
                layer = ( layer > -1 ? layer : TRACK_LAYER ) ;
                drawParameters[np].Layer = layer ;
                
                int color = _colors[ i % _colors.size() ] ;

                MarlinCED::add_layer_description(colName, layer);
                
                for( TrackerHitVec::const_iterator it = hits.begin();  it != hits.end() ; it++ ) {
                    
                    //hauke hoelbe: add id for picking
                    ced_hit_ID( (*it)->getPosition()[0],
                               (*it)->getPosition()[1],
                               (*it)->getPosition()[2],
                               marker, layer , size , color , trk->id() ) ;
                    
                } // hits
                
                const TrackState* ts = 0 ;
                
                switch( _drawHelixForTracks ){
                        
                        // Case 0 is the default value which takes the whatever first track state
                        // (e.g. InteractionPoint for ILD, or the simple track for test beam data)
                    case 0:
                        if (trk->getTrackStates().size() > 0) {
                            ts = trk->getTrackStates().at(0) ;
                        }
                        break ;
                        
                        // The rest states are pre-defined
                    case 1:  ts = trk->getTrackState( TrackState::AtIP          ) ; break ;
                    case 2:  ts = trk->getTrackState( TrackState::AtFirstHit    ) ; break ;
                        
                    case 3: // AtLastHit
                        
                        ts = (  trk->getTracks().empty() ?  trk->getTrackState( TrackState::AtLastHit )
                              :       trk->getTracks().back()->getTrackState( TrackState::AtLastHit ) )   ;
                        break ;
                        
                    case 4:  // AtCalo
                        
                        ts = (  trk->getTracks().empty() ?  trk->getTrackState( TrackState::AtCalorimeter )
                              :       trk->getTracks().back()->getTrackState( TrackState::AtCalorimeter ) )   ;
                        
                        ced_hit_ID( ts->getReferencePoint()[0],
                                   ts->getReferencePoint()[1],
                                   ts->getReferencePoint()[2],
                                   1 , layer , size*10 , color , trk->id() ) ;
                        break ;
                }
                
                if( ts !=0 ){
                    
                    double bField = Global::GEAR->getBField().at(  gear::Vector3D(0,0,0)  ).z() ;
                    
                    double pt;
                    if (bField != 0.0 && std::abs(ts->getOmega()) > 0.00001 )
                        pt = bField * 3e-4 / std::abs( ts->getOmega() ) ;
                    else
                        pt = 1.e10;
                    double charge = ( ts->getOmega() > 0. ?  1. : -1. ) ;
                    
                    double px = pt * std::cos(  ts->getPhi() ) ;
                    double py = pt * std::sin(  ts->getPhi() ) ;
                    double pz = pt * ts->getTanLambda() ;
                    
                    
#define Correct_Track_Params 1
#ifdef Correct_Track_Params
                    // start point for drawing ( PCA to reference point )
                    
                    double xs = ts->getReferencePoint()[0] -  ts->getD0() * sin( ts->getPhi() ) ;
                    double ys = ts->getReferencePoint()[1] +  ts->getD0() * cos( ts->getPhi() ) ;
                    double zs = ts->getReferencePoint()[2] +  ts->getZ0() ;
                    
#else  // assume track params have reference point at origin ...
                    
                    double xs = 0 -  ts->getD0() * sin( ts->getPhi() ) ;
                    double ys = 0 +  ts->getD0() * cos( ts->getPhi() ) ;
                    double zs = 0 +  ts->getZ0() ;
#endif
                    
                    
                    layer = ( layer > -1 ? layer : TRACKHELIX_LAYER ) ;
                    drawParameters[np].Layer = layer ;
                    
                    if( _drawHelixForTracks >= 0 && pt > 0.01 ) {
                        
                        const int ml = marker | ( layer << CED_LAYER_SHIFT ) ;
                        
                        int helixColor = ( _useColorForHelixTracks ? color : 0xdddddd ) ;
                        
                        MarlinCED::drawHelix( bField , charge, xs, ys, zs ,
                                             px, py, pz, ml ,  2 , helixColor  ,
                                             0.0, _helix_max_r,
                                             _helix_max_z, trk->id() ) ;

                    }
                    
                }
                
            } // track
            
        } else if( col->getTypeName() == LCIO::MCPARTICLE ){
            
            streamlog_out( DEBUG ) << "  drawing MCParticle collection " << std::endl ;
            
            
            double ecalR =  ( Global::GEAR->getEcalBarrelParameters().getExtent()[0] +
                             Global::GEAR->getEcalEndcapParameters().getExtent()[1]  ) / 2. ;
            
            double ecalZ = std::abs ( Global::GEAR->getEcalEndcapParameters().getExtent()[2] +
                                     Global::GEAR->getEcalEndcapParameters().getExtent()[3]  ) / 2. ;
            
            double hcalR =  ( Global::GEAR->getHcalBarrelParameters().getExtent()[0] +
                             Global::GEAR->getHcalEndcapParameters().getExtent()[1]  ) / 2. ;
            
            double hcalZ = std::abs ( Global::GEAR->getHcalEndcapParameters().getExtent()[2] +
                                     Global::GEAR->getHcalEndcapParameters().getExtent()[3]  ) / 2. ;
            
            
            
            for(int i=0; i<col->getNumberOfElements() ; i++){
                
                MCParticle* mcp = dynamic_cast<MCParticle*> ( col->getElementAt( i ) ) ;
                
                float charge = mcp->getCharge ();
                
                if( mcp-> getGeneratorStatus() != 1 && _usingParticleGun == false ) continue ; // stable particles only
                //if( mcp-> getNumberOfDaughters() != 0 ) continue ; // stable particles only
                //if( charge == 0.0 ) continue ;
                // 	  if( mcp-> getSimulatorStatus() != 0 ) continue ; // stable particles only
                //if( mcp->getDaughters().size() > 0  ) continue ;    // stable particles only
                // FIXME: need definition of stable particles (partons, decays in flight,...)
                
                if ( mcp->getEnergy() < _mcpECut ) continue ;           // ECut ?
                
                streamlog_out( DEBUG ) << "  drawing MCParticle pdg "
                << mcp->getPDG()
                << " genstat: " << mcp->getGeneratorStatus()
                << std::endl ;
                
                
                
                double px = mcp->getMomentum()[0];
                double py = mcp->getMomentum()[1];
                double pz = mcp->getMomentum()[2];
                
                
                double x = mcp->getVertex()[0] ;
                double y = mcp->getVertex()[1] ;
                double z = mcp->getVertex()[2] ;
                
                
                layer = ( layer > -1 ? layer : MCPARTICLE_LAYER ) ;
                drawParameters[np].Layer = layer ;
                
                //ced_describe_layer( colName.c_str() ,layer);
                MarlinCED::add_layer_description(colName, layer);
                
                
                
                if( std::fabs( charge ) > 0.0001  ) {
                    
                    
                    double bField = Global::GEAR->getBField().at(  gear::Vector3D(0,0,0)  ).z() ;
                    
                    streamlog_out( DEBUG ) << "  drawing MCParticle helix for p_t "
                    << sqrt(px*px+py*py)
                    << std::endl ;
                    
                    const int ml = marker | ( layer << CED_LAYER_SHIFT ) ;
                    
                    MarlinCED::drawHelix( bField , charge, x, y, z,
                                         px, py, pz, ml , size , 0x7af774  ,
                                         0.0,  _helix_max_r ,
                                         _helix_max_z, mcp->id() ) ;
                    
                } else { // neutral
                    
                    int color  ;
                    //          double r_min = 300 ;
                    double z_max ;
                    double r_max ;
                    
                    
                    switch(  std::abs(mcp->getPDG() )  ){
                            
                        case 22:
                            
                            color = 0xf9f920;          // photon
                            r_max = ecalR ;
                            z_max = ecalZ ;
                            
                            break ;
                            
                        case 12:  case 14: case 16: // neutrino
                            
                            color =  0xdddddd  ;
                            r_max = hcalR * 2. ;
                            z_max = hcalZ * 2. ;
                            break ;
                            
                        default:
                            color = 0xb900de  ;        // neutral hadron
                            r_max = hcalR ;
                            z_max = hcalZ ;
                    }
                    
                    
                    
                    double pt = hypot(px,py);
                    double p = hypot(pt,pz);
                    
                    double length = ( std::abs( pt/pz) > r_max/z_max ) ?  // hit barrel or endcap ?
                    r_max * p / pt  :  std::abs( z_max * p / pz ) ;
                    
                    ced_line_ID( x , y , z ,
                                length*px/p ,  length*py/p ,  length*pz/p ,
                                layer  , size, color, mcp->id() );
                    
                }
            }
        } else if( col->getTypeName() == LCIO::SIMTRACKERHIT ){
            
            
            layer = ( layer > -1 ? layer : SIMTRACKERHIT_LAYER ) ;
            drawParameters[np].Layer = layer ;
            
            //ced_describe_layer( colName.c_str() ,layer);
            MarlinCED::add_layer_description(colName, layer);
            
            
            
            for( int i=0, n=col->getNumberOfElements(); i<n ; i++ ){
                
                SimTrackerHit* h = dynamic_cast<SimTrackerHit*>( col->getElementAt(i) ) ;
                
                // color code by MCParticle
                const int mci = ( h->getMCParticle() ? h->getMCParticle()->id() :  0 ) ;
                int color = _colors[  mci % _colors.size() ] ;
                
                int id =    h->id();
                ced_hit_ID( h->getPosition()[0],
                           h->getPosition()[1],
                           h->getPosition()[2],
                           marker,layer, size , color, id ) ;
                
            }
            
        } else if( col->getTypeName() == LCIO::SIMCALORIMETERHIT ){
            
            layer = ( layer > -1 ? layer : SIMCALORIMETERHIT_LAYER ) ;
            drawParameters[np].Layer = layer ;
            
            MarlinCED::add_layer_description(colName, layer);
            
            for( int i=0, n=col->getNumberOfElements(); i<n ; i++ ){
                
                SimCalorimeterHit* h = dynamic_cast<SimCalorimeterHit*>( col->getElementAt(i) ) ;
                
                // color code by MCParticle
                //        const int mci = (  h->getNMCContributions() !=0  ?  h->getParticleCont(0)->id()  :  0 ) ;
                const int mci = (  h->getNMCContributions() !=0  && h->getParticleCont(0) ?  h->getParticleCont(0)->id()  :  0 ) ;
                int color = _colors[  mci % _colors.size() ] ;
                
                int id =    h->id();
                ced_hit_ID( h->getPosition()[0],
                           h->getPosition()[1],
                           h->getPosition()[2],
                           marker,layer, size , color, id ) ;
                
            }
            
        } else if(  col->getTypeName() == LCIO::TRACKERHIT       ||
                  col->getTypeName() == LCIO::TRACKERHITPLANE  ||
                  col->getTypeName() == LCIO::TRACKERHITZCYLINDER   ){
            
            int color = 0xee0044 ;
            
            layer = ( layer > -1 ? layer : TRACKERHIT_LAYER ) ;
            drawParameters[np].Layer = layer ;
            
            //ced_describe_layer( colName.c_str() ,layer);
            MarlinCED::add_layer_description(colName, layer);
            
            
            
            if( col->getTypeName() != LCIO::TRACKERHITPLANE ){
                
                // draw a marker at hit position
                
                LCTypedVector<TrackerHit> v( col ) ;
                MarlinCED::drawObjectsWithPosition( v.begin(), v.end() , marker, size , color, layer) ;
                
            } else { // LCIO::TRACKERHITPLANE
                
                lcio::CellIDDecoder<TrackerHitPlane> dec( col ) ;
                
                LCTypedVector<TrackerHitPlane> hits( col ) ;
                
                for( unsigned i=0,N=hits.size() ; i<N ; ++i){
                    
                    TrackerHitPlane* h = hits[i] ;
                    
                    gear::Vector3D p( h->getPosition()[0] ,  h->getPosition()[1] ,  h->getPosition()[2] ) ;
                    
                        
                    ced_hit_ID( p[0], p[1], p[2], marker, layer , size , color, h->id() );
                        
                    
                    // draw an additional line for strip hits 
                    if(  BitSet32( h->getType() )[ ILDTrkHitTypeBit::ONE_DIMENSIONAL ] ) {

                      double strip_half_length = 50. ; //FIXME: need to get the proper strip length (from the surface !?) 
                      gear::Vector3D v( strip_half_length , h->getV()[1] ,  h->getV()[0] , gear::Vector3D::spherical ) ;
                      
                      gear::Vector3D x0 = p - v ;
                      gear::Vector3D x1 = p + v ;
                      
                      ced_line_ID( x0[0], x0[1], x0[2], x1[0], x1[1], x1[2],
                                   layer , 1. , color, h->id() );
                    }
                }
            }
            
            
        } else if( col->getTypeName() == LCIO::CALORIMETERHIT ){
            
            int color = 0xee0000 ;
            
            layer = ( layer > -1 ? layer : CALORIMETERHIT_LAYER ) ;
            drawParameters[np].Layer = layer ;
            
            //ced_describe_layer( colName.c_str() ,layer);
            MarlinCED::add_layer_description(colName, layer);
            
            
            LCTypedVector<CalorimeterHit> v( col ) ;
            MarlinCED::drawObjectsWithPosition( v.begin(), v.end() , marker, size , color, layer ) ;
            
            
        } else if( col->getTypeName() == LCIO::RECONSTRUCTEDPARTICLE ){ //hauke
            
            
            layer = ( layer > -1 ? layer : RECOPARTICLE_LAYER ) ;
            drawParameters[np].Layer = layer ;
            
            MarlinCED::add_layer_description(colName, layer);
            
            int nelem = col->getNumberOfElements();
            
            float TotEn = 0.0;
            float TotPX = 0.0;
            float TotPY = 0.0;
            float TotPZ = 0.0;
            
            for (int ip(0); ip < nelem; ++ip) {
                
                int color = _colors[ ip % _colors.size() ] ;
                
                ReconstructedParticle * part = dynamic_cast<ReconstructedParticle*>(col->getElementAt(ip));
                TrackVec trackVec = part->getTracks();
                unsigned nTracks =  (unsigned)trackVec.size();
                ClusterVec clusterVec = part->getClusters();
                unsigned nClusters = (unsigned)clusterVec.size();
                
                float ene = part->getEnergy();
                float px  = (float)part->getMomentum()[0];
                float py  = (float)part->getMomentum()[1];
                float pz  = (float)part->getMomentum()[2];
                int type = (int)part->getType();
                
                TotEn += ene;
                TotPX += px;
                TotPY += py;
                TotPZ += pz;
                
                streamlog_out( DEBUG2)  << "Particle : " << ip
                << " type : " << type
                << " PX = " << px
                << " PY = " << py
                << " PZ = " << pz
                << " E  = " << ene << std::endl;
                
                if (nClusters > 0 ) {
                    
                    for (unsigned icluster=0; icluster<nClusters; ++icluster) {
                        Cluster * cluster = clusterVec[icluster];
                        CalorimeterHitVec hitvec = cluster->getCalorimeterHits();
                        int nHits = (int)hitvec.size();
                        for (int iHit = 0; iHit < nHits; ++iHit) {
                            CalorimeterHit * hit = hitvec[iHit];
                            float x = hit->getPosition()[0];
                            float y = hit->getPosition()[1];
                            float z = hit->getPosition()[2];
                            ced_hit_ID(x,y,z,marker, layer ,size,color,part->id()); 
                        }
                        
                    }
                }
                
                if (nTracks > 0 ) {
                    
                    for (unsigned itrack=0; itrack<nTracks; ++itrack) {
                        
                        Track * trk = trackVec[itrack];
                        
                        // -- collect hits from all track segments
                        TrackerHitVec tHV ;
                        if( !trk->getTracks().empty() ){
                            
                            streamlog_out( DEBUG ) << " -- track has "<< trk->getTracks().size() << " subtracks - will use these for displaying hits "
                            << std::endl ;
                            
                            std::copy( trk->getTrackerHits().begin() , trk->getTrackerHits().end() , std::back_inserter(  tHV ) ) ;
                            
                            for( unsigned j=0 ,N = trk->getTracks().size() ; j<N ; ++j ){
                                
                                const Track* t = trk->getTracks()[j] ; 
                                
                                streamlog_out( DEBUG ) << " -- track j= "<< j << " has " <<  t->getTrackerHits().size() << " hits  "
                                << std::endl ;
                                
                                std::copy( t->getTrackerHits().begin() , t->getTrackerHits().end() , std::back_inserter(  tHV ) ) ;
                            }
                            
                        } 
                        const TrackerHitVec& hitvec = (  tHV.empty() ?  trk->getTrackerHits()  :  tHV  ) ; 
                        
                        int nHits = (int)hitvec.size();
                        
                        if(nHits > 0){
                            for (int iHit = 0; iHit < nHits; ++iHit) {
                                TrackerHit * hit = hitvec[iHit];
                                float x = (float)hit->getPosition()[0];
                                float y = (float)hit->getPosition()[1];
                                float z = (float)hit->getPosition()[2];
                                ced_hit_ID(x,y,z,marker, layer,size,color,part->id());
                            }
                        }
                        if((nHits==0 || _drawHelixForPFOs==1) && std::fabs(part->getCharge())>0.001){
                            
                            TrackState *ts=0;
                            
                            if (trk->getTrackStates().size() > 0) {
                                ts = trk->getTrackStates().at(0) ;
                            }
                            if(ts!=0){
                                
                                double bField = Global::GEAR->getBField().at(  gear::Vector3D(0,0,0)  ).z() ; 
                                
                                double pt;
                                if (bField != 0.0 && std::abs(ts->getOmega()) > 0.00001 ){
                                    pt = bField * 3e-4 / std::abs( ts->getOmega() ) ;
                                }else{
                                    pt = 1.e10;
                                }
                                double charge = ( ts->getOmega() > 0. ?  1. : -1. ) ;
                                
                                double px = pt * std::cos(  ts->getPhi() ) ;
                                double py = pt * std::sin(  ts->getPhi() ) ;
                                double pz = pt * ts->getTanLambda() ;
                                
                                
                                double xs = ts->getReferencePoint()[0] -  ts->getD0() * sin( ts->getPhi() ) ;
                                double ys = ts->getReferencePoint()[1] +  ts->getD0() * cos( ts->getPhi() ) ;
                                double zs = ts->getReferencePoint()[2] +  ts->getZ0() ;
                                
                                int helixColor = ( _useColorForHelixTracks ? color : 0xdddddd ) ;

                                //helix
                                MarlinCED::drawHelix(bField, charge, xs, ys, zs, px, py, pz, marker|(layer<<CED_LAYER_SHIFT), size/2, 
                                                     helixColor, 0.0, _helix_max_r, _helix_max_z, part->id() ); //hauke: add id
                            }
                        }
                    }
                }
            }
            
            streamlog_out( DEBUG2) << std::endl
            << "Total Energy and Momentum Balance of Event" << std::endl
            << "Energy = " << TotEn
            << " PX = " << TotPX
            << " PY = " << TotPY
            << " PZ = " << TotPZ 
            << std::endl 
            << std::endl;
            
            
            //LCTypedVector<CalorimeterHit> v( col ) ;
            //MarlinCED::drawObjectsWithPosition( v.begin(), v.end() , marker, size , color, layer ) ;
            
            
            //std::cout << std::endl << "#############        print RECONSTRUCTEDPARTICLE #################"  << std::endl;
        }
        
    } // while
    
    streamlog_out( MESSAGE ) << " ++++++++ collections shown on layer [ evt: " << evt->getEventNumber() 
    <<  " run: " << evt->getRunNumber() << " ] :   +++++++++++++ " << std::endl ;
    
    for(unsigned np=0 ; np < nCols ; ++np){
        
        const std::string & colName = drawParameters[np].ColName ;
        //     int size =   drawParameters[np].Size ;
        //     int marker = drawParameters[np].Marker ;
        int layer =  drawParameters[np].Layer ;
        
        streamlog_out( MESSAGE )  << "    +++++  " << colName <<  "\t  on layer: " << layer << std::endl ;
    }
    streamlog_out( MESSAGE ) << " ++++++++ use shift-[LN] for LN>10  +++++++++++++ " << std::endl ;
    
    
    
    //++++++++++++++++++++++++++++++++++++


    MarlinCED::draw(this, _waitForKeyboard );


    //++++++++++++++++++++++++++++++++++++
    
    _nEvt ++ ;
}



void CEDViewer::check( LCEvent * evt ) { 
    // nothing to check here - could be used to fill checkplots in reconstruction processor
}

void CEDViewer::printParticle(int id, LCEvent * evt){
    streamlog_out( MESSAGE )  << "CEDViewer::printParticle id: " << id << std::endl;
} 


void CEDViewer::end(){ 
    
    streamlog_out(DEBUG2 ) << "end() :" << " processed " << _nEvt 
    << " events in " << _nRun << " runs "
    << std::endl ;
    
}

