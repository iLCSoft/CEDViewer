/* -*- Mode: C++; tab-width: 2; indent-tabs-mode: nil; c-basic-offset: 2 -*- */
#include "DDCEDViewer.h"

#include <iostream>

#include <EVENT/LCCollection.h>
#include <EVENT/Cluster.h>
#include <EVENT/Track.h>
#include <EVENT/MCParticle.h>
#include <EVENT/ReconstructedParticle.h>
#include <EVENT/SimCalorimeterHit.h>
#include <EVENT/SimTrackerHit.h>

#include <EVENT/ReconstructedParticle.h>

#include <UTIL/LCTypedVector.h>
#include <UTIL/CellIDDecoder.h>
#include <UTIL/ILDConf.h>

#include "DDMarlinCED.h"
#include <cmath>
#include "DD4hep/LCDD.h"
#include "DD4hep/DD4hepUnits.h" 

#include "TVector3.h"

using namespace lcio ;
using namespace marlin ;
using namespace DD4hep::Geometry ;
using namespace DD4hep;
using namespace DD4hep::DDRec ;

#define SIMTRACKERHIT_LAYER 1
#define TRACKERHIT_LAYER 2
#define TRACK_LAYER 3
#define TRACKHELIX_LAYER 3

#define SIMCALORIMETERHIT_LAYER 4
#define CALORIMETERHIT_LAYER 5
#define CLUSTER_LAYER 6

#define MCPARTICLE_LAYER 7
#define RECOPARTICLE_LAYER 8


DDCEDViewer aDDCEDViewer ;

DDCEDViewer::DDCEDViewer() : Processor("DDCEDViewer") {
    
    // modify processor description
    _description = "DDCEDViewer: event display of LCIO objects  - based on CED by A.Zhelezov." ;    //todo
    
    StringVec cluExample ;
    cluExample.push_back("CalorimeterHits" ) ;
    cluExample.push_back("0" ) ;
    cluExample.push_back("3" ) ;
    
    registerOptionalParameter( "DrawCollection" ,
                              "collection to be displayed ( ColName, marker type[0-2], size)"  ,
                              _drawCollections ,
                              cluExample ,
                              cluExample.size() ) ;
    
    StringVec TPClayerExample ;
    TPClayerExample.push_back("TPCTracks" ) ;
    TPClayerExample.push_back("0" ) ;
    TPClayerExample.push_back("3" ) ;
    TPClayerExample.push_back("4" ) ;
    
    /*****options for LCIO drawing*****/
    registerOptionalParameter( "DrawInLayer" ,
                              "collection to be displayed ( ColName, marker type[0-2], size, layer)",
                              _drawCollectionsLayer ,
                              TPClayerExample ,
                              TPClayerExample.size() ) ;
    
    
    registerProcessorParameter( "DrawHelixForTrack" ,
                               "draw a helix for Track objects: -1: none, 0: default, 1: atIP, 2: atFirstHit, 3: atLastHit, 4: atCalorimeter",
                               _drawHelixForTracks ,
                               0  ) ;
        
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
    
    registerProcessorParameter( "UseTrackerExtentForLimitsOfHelix" ,
                               "Use the gear parameters to define the max extent of drawing a helix",
                               //_useTPCForLimitsOfHelix ,
                               _useTrackerForLimitsOfHelix ,
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

    /*****options for detector drawing*****/
    registerProcessorParameter( "DrawDetector" ,
                  "Call the drawDetector function else call draw function",
                  _begin ,
                  (bool)0  ) ;
    //example StringVec to indicate the function how the list is meant to be read out.
    StringVec DetailledLayerExample ; DetailledLayerExample.push_back("NameToBeDrawnDetailly" ) ;
    registerOptionalParameter( "DetailledDrawing" ,
            "List of detector names that are printed in more detail.",
            _detailled,
            DetailledLayerExample,
            1  ) ;  
    //  
    registerOptionalParameter( "DrawSurfaces" ,
            "Draw the geometry as a set of individual surfaces (if available) instead of simplified structures.",
            _surfaces ,
            (bool)0  ) ;
 
}

void DDCEDViewer::init() {
    
    // usually a good idea to
    streamlog_out(DEBUG) << "   init called  " << std::endl ;
    printParameters() ;
    DDMarlinCED::init(this) ;
    
    
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

void DDCEDViewer::processRunHeader( LCRunHeader* run) {
    _nRun++ ;
}

void DDCEDViewer::processEvent( LCEvent * evt ) {
    //-----------------------------------------------------------------------
    // Reset drawing buffer and START drawing collection

    DDMarlinCED::newEvent(this) ;
    
    //added by Thorben Quast for removing GEAR dependencies
    DD4hep::Geometry::LCDD& lcdd = DD4hep::Geometry::LCDD::getInstance();

    DDMarlinCED::drawDD4hepDetector(lcdd, _surfaces, _detailled);
    DDCEDViewer::drawDD4LCIO(evt, lcdd);

    DDMarlinCED::draw(this, _waitForKeyboard );

    //++++++++++++++++++++++++++++++++++++

    _nEvt ++ ;
}



void DDCEDViewer::check( LCEvent * evt ) { 
    // nothing to check here - could be used to fill checkplots in reconstruction processor
}

void DDCEDViewer::printParticle(int id, LCEvent * evt){
    streamlog_out( MESSAGE )  << "CEDViewer::printParticle id: " << id << std::endl;
} 


void DDCEDViewer::end(){ 
    streamlog_out(DEBUG2 ) << "end() :" << " processed " << _nEvt 
    << " events in " << _nRun << " runs "
    << std::endl ;
}

void DDCEDViewer::drawDD4LCIO(LCEvent * evt, DD4hep::Geometry::LCDD& lcdd){
DDCEDPickingHandler &pHandler=DDCEDPickingHandler::getInstance();
    
    pHandler.update(evt);
    
    //-----------------------------------------------------------------------
    if( _useTrackerForLimitsOfHelix ){
        //only draw the helix in the TPC within the tracker volume
        _helix_max_r = getTrackerExtent(lcdd)[0];
        _helix_max_z = getTrackerExtent(lcdd)[1];
    }
    //add DrawInLayer to
    //std::vector< DrawParameters > drawParameters ;
    this->drawParameters.reserve( _drawCollections.size() + _drawCollectionsLayer.size()  ) ;
    
    if( parameterSet( "DrawCollection" ) ) {
        
        unsigned index = 0 ;
        while( index < _drawCollections.size() ){
            
            const std::string & colName = _drawCollections[ index++ ] ;
            int marker = std::atoi( _drawCollections[ index++ ].c_str() ) ;
            int size = std::atoi( _drawCollections[ index++ ].c_str() ) ;
            //std::cout << "#######################      SIZE: " << size << " #######################################" << std::endl;
            int layer = -1 ;
            
            this->drawParameters.push_back(DrawParameters( colName,size,marker,layer ) );
        }
    }
    if( parameterSet( "DrawInLayer" ) ) {
        
        unsigned index = 0 ;
        while( index < _drawCollectionsLayer.size() ){
            
            const std::string & colName = _drawCollectionsLayer[ index++ ] ;
            int marker = std::atoi( _drawCollectionsLayer[ index++ ].c_str() ) ;
            int size   = std::atoi( _drawCollectionsLayer[ index++ ].c_str() ) ;
            int layer  = std::atoi( _drawCollectionsLayer[ index++ ].c_str() ) ;
            
            this->drawParameters.push_back(DrawParameters( colName,size,marker,layer ) );
            //std::cout << "layer: " << layer << " description: " << colName << std::endl; //hauke
        }
    }
    
    unsigned nCols = this->drawParameters.size() ;
    //drawing with 'standard' routines
    for(unsigned np=0 ; np < nCols ; ++np){
        
        const std::string & colName = this->drawParameters[np].ColName ;
        int size =   this->drawParameters[np].Size ;
        int marker = this->drawParameters[np].Marker ;
        int layer =  this->drawParameters[np].Layer ;
        
        LCCollection* col = 0 ;
        try{
            
            col = evt->getCollection( colName ) ;
            
        }catch(DataNotAvailableException &e){
            
            streamlog_out( WARNING ) << " collection " << colName <<  " not found ! "   << std::endl ;
            continue ;
        }
        
        //specific drawing of individual collections
        if (colName == "JetOut")
            DDCEDViewer::drawJets(lcdd, layer, np, colName, col);
        else{
            //draw the individual collections
            if( col->getTypeName() == LCIO::CLUSTER ){
                DDCEDViewer::drawCluster(lcdd, layer, np, colName, marker, col, size);
            } else if( col->getTypeName() == LCIO::TRACK ){
                DDCEDViewer::drawTrack(lcdd, layer, np, colName, marker, col, size);
            } else if( col->getTypeName() == LCIO::MCPARTICLE ){
                DDCEDViewer::drawMCParticle(lcdd, layer, np, colName, marker, col, size);
            } else if( col->getTypeName() == LCIO::SIMTRACKERHIT ){
                DDCEDViewer::drawSIMTrackerHit(layer, np, colName, marker, col, _colors, size);
            } else if( col->getTypeName() == LCIO::SIMCALORIMETERHIT ){
                DDCEDViewer::drawSIMCalorimeterHit(layer, np, colName, marker, col, _colors, size);
            } else if(  col->getTypeName() == LCIO::TRACKERHIT       ||
                      col->getTypeName() == LCIO::TRACKERHITPLANE  ||
                      col->getTypeName() == LCIO::TRACKERHITZCYLINDER   ){
                DDCEDViewer::drawTrackerHit(layer, np, colName, marker, col, size);
            } else if( col->getTypeName() == LCIO::CALORIMETERHIT ){
                DDCEDViewer::drawCalorimeterHit(layer, np, colName, marker, col, size);
            } else if( col->getTypeName() == LCIO::RECONSTRUCTEDPARTICLE ){ 
                DDCEDViewer::drawReconstructedParticle(lcdd, layer, np, colName, marker, col, size);
            }    
        }
    }
}

void DDCEDViewer::drawJets(DD4hep::Geometry::LCDD& lcdd, int& layer, unsigned& np, std::string colName, LCCollection* col){
    streamlog_out( MESSAGE )  << " drawing jets from collection " << colName << std::endl ;

    //some default color: to be extended
    float RGBAcolor[4] = {0., 0., 1.0, 0.3};
    int color = int(RGBAcolor[2]*(15*16+15)) + int(RGBAcolor[1]*(15*16+15))*16*16+ int(RGBAcolor[0]*(15*16+15))*16*16*16*16;

    DDMarlinCED::add_layer_description(colName, layer);

    for (int j=0; j < col->getNumberOfElements(); ++j) { //number of elements in a jet
        ReconstructedParticle * jet = dynamic_cast<ReconstructedParticle*>( col->getElementAt(j) );
        streamlog_out( DEBUG )  <<   "     - jet energy " << jet->getEnergy() << std::endl ;
        streamlog_out( MESSAGE )  <<   "     - jet energy " << jet->getEnergy() << std::endl ;
        TVector3 v(jet->getMomentum()[0], jet->getMomentum()[1], jet->getMomentum()[2]); 
        const ReconstructedParticleVec & pv = jet->getParticles();
        float pt_norm = 0.0;
        for (unsigned int k = 0; k<pv.size(); ++k){
            const double * pm = pv[k]->getMomentum();
            TVector3 pp(pm[0], pm[1] , pm[2]);
            TVector3 ju = v.Unit();
            TVector3 pt = pp - (ju.Dot(pp))*ju;

            pt_norm += pt.Mag();

            int LineSize = 1;
            // start point
            float refx = 0.0;
            float refy = 0.0;
            float refz = 0.0;
            float momScale = 100;

            int layerIp = layer;
            ced_line_ID(refx, refy, refz, momScale*pm[0], momScale*pm[1], momScale*pm[2], layerIp, LineSize, color, pv[k]->id());
        }
                       
        double center_c[3] = {0., 0., 0. };
        double rotation_c[3] = { 0.,  v.Theta()*180./M_PI , v.Phi()*180./M_PI };
        
        double scale_pt = 20;
        double scale_mom = 25;
        double min_pt = 50;
      
        ced_cone_r_ID( min_pt + scale_pt*pt_norm , scale_mom*v.Mag() , center_c, rotation_c, layer, RGBAcolor,jet->id()); 
        
    }
}


void DDCEDViewer::drawCluster(DD4hep::Geometry::LCDD& lcdd, int& layer, unsigned& np, std::string colName, int& marker, LCCollection* col, int& size){
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
        this->drawParameters[np].Layer = layer ;

        int color =  _colors[ i % _colors.size() ] ;
        for( CalorimeterHitVec::const_iterator it = hits.begin();  it != hits.end() ; it++ ) {
            //hauke hoelbe: add id for picking!
            ced_hit_ID( (*it)->getPosition()[0],
                       (*it)->getPosition()[1],
                       (*it)->getPosition()[2],
                       marker, layer, size , color, clu->id() ) ;
        } 
        float x = clu->getPosition()[0] ;
        float y = clu->getPosition()[1] ;
        float z = clu->getPosition()[2] ;
        ced_hit_ID( x,y,z, marker, layer , size*3 , color, clu->id() ) ;
        LCVector3D v(x,y,z) ;
        LCVector3D d(1.,1.,1.) ;
        float length = 100 + 500 * ( clu->getEnergy() - emin )  / ( emax - emin )  ;
        d.setMag( length  ) ;
        d.setTheta( clu->getITheta() ) ;
        d.setPhi( clu->getIPhi() ) ;
        LCVector3D dp( v + d ) , dm( v - d )   ;
        ced_line_ID( dp.x() , dp.y() , dp.z(),
                    dm.x() , dm.y() , dm.z(),
                    layer  , 1 , color, clu->id() );
    }
}

void DDCEDViewer::drawTrack(DD4hep::Geometry::LCDD& lcdd, int& layer, unsigned& np, std::string colName, int& marker, LCCollection* col, int& size){
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
        this->drawParameters[np].Layer = layer ;
        int color = _colors[ i % _colors.size() ] ;
        DDMarlinCED::add_layer_description(colName, layer);
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
            double* bFieldVector = new double[3];
            lcdd.field().combinedMagnetic(Position(0,0,0), bFieldVector) ;
            double bField = bFieldVector[2] / dd4hep::tesla;
            delete bFieldVector;
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
            this->drawParameters[np].Layer = layer ;
            if( _drawHelixForTracks >= 0 && pt > 0.01 ) {
                const int ml = marker | ( layer << CED_LAYER_SHIFT ) ;
                int helixColor = ( _useColorForHelixTracks ? color : 0xdddddd ) ;
                DDMarlinCED::drawHelix( bField , charge, xs, ys, zs ,
                                     px, py, pz, ml ,  2 , helixColor  ,
                                     0.0, _helix_max_r,
                                     _helix_max_z, trk->id() ) ;
            }
        }
    }
}

void DDCEDViewer::drawMCParticle(DD4hep::Geometry::LCDD& lcdd, int& layer, unsigned& np, std::string colName, int& marker, LCCollection* col, int& size) {
    streamlog_out( DEBUG ) << "  drawing MCParticle collection " << std::endl ;
    //new line drawing implemented by Thorben Quast 07 August 2015
    CalorimeterDrawParams ECalBarrelParams = getCalorimeterParameters(lcdd, "ECalBarrel");
    CalorimeterDrawParams ECalEndcapParams = getCalorimeterParameters(lcdd, "ECalEndcap");
    CalorimeterDrawParams HCalBarrelParams = getCalorimeterParameters(lcdd, "HCalBarrel");
    CalorimeterDrawParams HCalEndcapParams = getCalorimeterParameters(lcdd, "HCalEndcap");
    for(int i=0; i<col->getNumberOfElements() ; i++){
        MCParticle* mcp = dynamic_cast<MCParticle*> ( col->getElementAt( i ) ) ;
        float charge = mcp->getCharge ();
        if( mcp-> getGeneratorStatus() != 1 && _usingParticleGun == false ) 
            continue ; // stable particles only
        if ( mcp->getEnergy() < _mcpECut ) 
            continue ;           // ECut ?
        
        streamlog_out( DEBUG ) << "  drawing MCParticle pdg "
        << mcp->getPDG()
        << " genstat: " << mcp->getGeneratorStatus()
        << std::endl ;
        layer = ( layer > -1 ? layer : MCPARTICLE_LAYER ) ;
        this->drawParameters[np].Layer = layer ;
        DDMarlinCED::add_layer_description(colName, layer);
        
        double px = mcp->getMomentum()[0];
        double py = mcp->getMomentum()[1];
        double pz = mcp->getMomentum()[2];
        double pt = sqrt(pow(px,2)+pow(py,2));
        double p = sqrt(pow(pt,2)+pow(pz,2));
        double x = mcp->getVertex()[0] ;
        double y = mcp->getVertex()[1] ;
        double z = mcp->getVertex()[2] ;
        
        if( std::fabs( charge ) > 0.0001  ) {
            double* bFieldVector = new double[3];
            lcdd.field().combinedMagnetic(Position(0,0,0), bFieldVector) ;
            double bField = bFieldVector[2] / dd4hep::tesla;
            delete bFieldVector;
            streamlog_out( DEBUG ) << "  drawing MCParticle helix for p_t "
            << sqrt(px*px+py*py)
            << std::endl ;
            const int ml = marker | ( layer << CED_LAYER_SHIFT ) ;
            //temporary variables
            double _hmr = _helix_max_r;
            double _hmz = _helix_max_z;
            DDMarlinCED::drawHelix( bField , charge, x, y, z,
                                 px, py, pz, ml , size , 0x7af774  ,
                                 0.0,  _hmr ,
                                 _hmz, mcp->id() ) ;
        } else { // neutral
            int color  ;
            double length, yokeR, yokeZ;
            switch(  std::abs(mcp->getPDG() )  ){
                //refactored length calculation (T. Quast 7 Aug 15) 
                case 22:
                    color = 0xf9f920;          // photon
                    length = calculateTrackLength(ECalBarrelParams, ECalEndcapParams, x, y, z, px, py, pz);
                    break ;
                case 12:  case 14: case 16: // neutrino
                    color =  0xdddddd  ;
                    yokeR = getYokeExtent(lcdd)[0];
                    yokeZ = getYokeExtent(lcdd)[1];
                    length = (fabs(pt/pz) > yokeR/yokeZ) ?
                            yokeR * sqrt(1. + pow(pz/pt,2)):
                            yokeZ * sqrt(1. + pow(pt/pz,2));
                    break ;
                default:
                    color = 0xb900de  ;        // neutral hadron
                    length = calculateTrackLength(HCalBarrelParams, HCalEndcapParams, x, y, z, px, py, pz);
            }
            //tracks with vertex outside the according calorimeter are not drawn, length is passed as 0
            ced_line_ID( x , y , z ,
                        x + length*px/p ,  y + length*py/p ,  z + length*pz/p ,
                        layer  , size, color, mcp->id() );
        }
    }
}

void DDCEDViewer::drawSIMTrackerHit(int& layer, unsigned& np, std::string colName, int& marker, LCCollection* col, std::vector<int>& _colors, int& size){
    layer = ( layer > -1 ? layer : SIMTRACKERHIT_LAYER ) ;
    this->drawParameters[np].Layer = layer ;
    DDMarlinCED::add_layer_description(colName, layer);

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
}

void DDCEDViewer::drawSIMCalorimeterHit(int& layer, unsigned& np, std::string colName, int& marker, LCCollection* col, std::vector<int>& _colors, int& size){
    layer = ( layer > -1 ? layer : SIMCALORIMETERHIT_LAYER ) ;
    this->drawParameters[np].Layer = layer ;
    
    DDMarlinCED::add_layer_description(colName, layer);
    
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
}

void DDCEDViewer::drawTrackerHit(int& layer, unsigned& np, std::string colName, int& marker, LCCollection* col, int& size){
    int color = 0xee0044 ;
    layer = ( layer > -1 ? layer : TRACKERHIT_LAYER ) ;
    this->drawParameters[np].Layer = layer ;
    //ced_describe_layer( colName.c_str() ,layer);
    DDMarlinCED::add_layer_description(colName, layer);

    if( col->getTypeName() != LCIO::TRACKERHITPLANE ){
        // draw a marker at hit position
        LCTypedVector<TrackerHit> v( col ) ;
        DDMarlinCED::drawObjectsWithPosition( v.begin(), v.end() , marker, size , color, layer) ;
    } else { // LCIO::TRACKERHITPLANE
        lcio::CellIDDecoder<TrackerHitPlane> dec( col ) ;
        LCTypedVector<TrackerHitPlane> hits( col ) ;
        for( unsigned i=0,N=hits.size() ; i<N ; ++i){
            TrackerHitPlane* h = hits[i] ;
            TVector3 p(h->getPosition()[0] ,  h->getPosition()[1] ,  h->getPosition()[2]);
            ced_hit_ID( p.X(), p.Y(), p.Z(), marker, layer , size , color, h->id() );
            // draw an additional line for strip hits 
            if(  BitSet32( h->getType() )[ ILDTrkHitTypeBit::ONE_DIMENSIONAL ] ) {
                double strip_half_length = 50. ; //FIXME: need to get the proper strip length (from the surface !?) 
                TVector3 v(1,1,1);//must be initialized as non zero to set spherical coordinates
                v.SetMag(strip_half_length);
                v.SetPhi(h->getV()[1]);
                v.SetTheta(h->getV()[0]);
                TVector3 x0 = p - v;
                TVector3 x1 = p + v;
                ced_line_ID( x0.X(), x0.Y(), x0.Z(), x1.X(), x1.Y(), x1.Z(),
                           layer , 1. , color, h->id() );
            }
        }
    }   
}

void DDCEDViewer::drawCalorimeterHit(int& layer, unsigned& np, std::string colName, int& marker, LCCollection* col, int& size){
    int color = 0xee0000 ;
    layer = ( layer > -1 ? layer : CALORIMETERHIT_LAYER ) ;
    this->drawParameters[np].Layer = layer ;
    DDMarlinCED::add_layer_description(colName, layer);
    LCTypedVector<CalorimeterHit> v( col ) ;
    DDMarlinCED::drawObjectsWithPosition( v.begin(), v.end() , marker, size , color, layer ) ;
}

void DDCEDViewer::drawReconstructedParticle(DD4hep::Geometry::LCDD& lcdd, int& layer, unsigned& np, std::string colName, int& marker, LCCollection* col, int& size){
    //hauke
    layer = ( layer > -1 ? layer : RECOPARTICLE_LAYER ) ;
    this->drawParameters[np].Layer = layer ;
    
    DDMarlinCED::add_layer_description(colName, layer);
    
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
                        double* bFieldVector = new double[3];
                        lcdd.field().combinedMagnetic(Position(0,0,0), bFieldVector) ;
                        double bField = bFieldVector[2] / dd4hep::tesla;
                        delete bFieldVector;
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
                        DDMarlinCED::drawHelix(bField, charge, xs, ys, zs, px, py, pz, marker|(layer<<CED_LAYER_SHIFT), size/2, 
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
}


/*
void DDCEDViewer::drawJets(DD4hep::Geometry::LCDD& lcdd, int& layer, unsigned& np, std::string colName, LCCollection* col){
    streamlog_out( MESSAGE )  << " drawing jets from collection " << colName << std::endl ;

    //some default color: to be extended
    float RGBAcolor[4] = {0., 0., 1.0, 0.3};
    int color = int(RGBAcolor[2]*(15*16+15)) + int(RGBAcolor[1]*(15*16+15))*16*16+ int(RGBAcolor[0]*(15*16+15))*16*16*16*16;

    DDMarlinCED::add_layer_description(colName, layer);

    for (int j=0; j < col->getNumberOfElements(); ++j) { //number of elements in a jet
        ReconstructedParticle * jet = dynamic_cast<ReconstructedParticle*>( col->getElementAt(j) );
        streamlog_out( DEBUG )  <<   "     - jet energy " << jet->getEnergy() << std::endl ;
        streamlog_out( MESSAGE )  <<   "     - jet energy " << jet->getEnergy() << std::endl ;
        TVector3 v(jet->getMomentum()[0], jet->getMomentum()[1], jet->getMomentum()[2]); 
        const ReconstructedParticleVec & pv = jet->getParticles();
        float pt_norm = 0.0;
        for (unsigned int k = 0; k<pv.size(); ++k){
            const double * pm = pv[k]->getMomentum();
            TVector3 pp(pm[0], pm[1] , pm[2]);
            TVector3 ju = v.Unit();
            TVector3 pt = pp - (ju.Dot(pp))*ju;

            pt_norm += pt.Mag();

            int LineSize = 1;
            // start point
            float refx = 0.0;
            float refy = 0.0;
            float refz = 0.0;
            float momScale = 100;

            int layerIp = layer;
            ced_line_ID(refx, refy, refz, momScale*pm[0], momScale*pm[1], momScale*pm[2], layerIp, LineSize, color, pv[k]->id());
        }
                       
        double center_c[3] = {0., 0., 0. };
        double rotation_c[3] = { 0.,  v.Theta()*180./M_PI , v.Phi()*180./M_PI };
        
        double scale_pt = 20;
        double scale_mom = 25;
        double min_pt = 50;
      
        ced_cone_r_ID( min_pt + scale_pt*pt_norm , scale_mom*v.Mag() , center_c, rotation_c, layer, RGBAcolor,jet->id()); 
        
    }
}
*/

/********************************************************************
//Helper functions, might be exported into a utility.cc at some point
********************************************************************/

//get the outer extents of the tracker
double* getTrackerExtent(DD4hep::Geometry::LCDD& lcdd){
    double* extent = new double[2];
    extent[0] =lcdd.constant<double>("tracker_region_rmax")/dd4hep::mm;
    extent[1] = lcdd.constant<double>("tracker_region_zmax")/dd4hep::mm;
    return extent;
}
//get the outer extents of the yoke
double* getYokeExtent(DD4hep::Geometry::LCDD& lcdd){
    double* extent = new double[2];
    const std::vector< DD4hep::Geometry::DetElement>& calorimeters     = lcdd.detectors( "calorimeter" ) ;
    DD4hep::Geometry::DetElement yoke;
    for( unsigned i=0,n=calorimeters.size() ; i<n ; ++i ){
        std::string detName = calorimeters[i].name();
        bool isYokeBarrel = (detName == "YokeBarrel") ;
        bool isYokeEndcap = (detName == "YokeEndcap") ;
        if (!isYokeBarrel && !isYokeEndcap)
            continue;
        yoke = calorimeters[i];
        DD4hep::DDRec::LayeredCalorimeterData* yokeGeo;
        try{ 
            yokeGeo = yoke.extension<DD4hep::DDRec::LayeredCalorimeterData>() ; 
        } catch(std::runtime_error& e){
            streamlog_out( MESSAGE ) <<  "MC Particles for " << detName << " cannot be drawn"<<std::endl;
            extent[0] = extent[1] = 0;
            return extent;
        }
        if (isYokeBarrel)
            extent[0] = yokeGeo->extent[1]/dd4hep::mm;
        if (isYokeEndcap)
            extent[1] = yokeGeo->extent[3]/dd4hep::mm; 
    }
    return extent;
}

//calculates and returns the relevant calorimeter parameters for track length calculations.
CalorimeterDrawParams getCalorimeterParameters(DD4hep::Geometry::LCDD& lcdd, std::string name, bool selfCall){
    CalorimeterDrawParams params;
    if (selfCall)
        name[1] = tolower(name[1]); 
    const std::vector< DD4hep::Geometry::DetElement>& calorimeters     = lcdd.detectors( "calorimeter" ) ;
    DD4hep::Geometry::DetElement calo;
    for( unsigned i=0,n=calorimeters.size() ; i<n ; ++i ){
        if ((std::string)calorimeters[i].name() == name){
            calo = calorimeters[i];
            break;
        }
    }
    DD4hep::DDRec::LayeredCalorimeterData* caloGeo;
    try{ 
        caloGeo = calo.extension<DD4hep::DDRec::LayeredCalorimeterData>() ; 
    } catch(std::runtime_error& e){
            if (!selfCall)
                return getCalorimeterParameters(lcdd, name, true); 
            else{
                streamlog_out( MESSAGE ) <<  "MC Particles for " << name << " cannot be drawn"<<std::endl;
                params.delta_z = -1;    //no spatial extension --> no drawing
                return params;
            }
    }
    if (caloGeo->layoutType == DD4hep::DDRec::LayeredCalorimeterData::BarrelLayout){
        params.r_inner = caloGeo->extent[0]/dd4hep::mm;
        params.delta_r = caloGeo->extent[1]/dd4hep::mm - params.r_inner;
        params.z_0 = 0.;
        params.delta_z = caloGeo->extent[3]/dd4hep::mm;
    }else{
        params.r_inner = caloGeo->extent[0]/dd4hep::mm;
        params.delta_r = caloGeo->extent[1]/dd4hep::mm - params.r_inner;
        params.z_0 = caloGeo->extent[2]/dd4hep::mm;
        params.delta_z = (caloGeo->extent[3]/dd4hep::mm - params.z_0);  //we are interested in the full length! 
                                                        //CEDGeoTube only requires half length as an argument       
    }

    return params;
}

//It suffices to perform the calculations in the first quadrant due to the detector's symmetry.
//The signs of the tracks' directions are ultimately determined by the momenta.
double calculateTrackLength(CalorimeterDrawParams barrel, CalorimeterDrawParams endcap, double x, double y, double z, double px, double py, double pz){
    if (barrel.delta_z == -1 || endcap.delta_z == -1) return 0;   //the case if the parameters could not be loaded properly
    
    double length;
    double rel_X0 = 0.5;    //mean interaction length traversing the material perpendicularly; must be <= 1 !
    double pt = sqrt(px*px + py*py);
    double pt_over_pz = pt/fabs(pz);

    double r = sqrt(x*x+y*y);
    if (r > barrel.r_inner || r > endcap.r_inner) return 0;
    double p = 2 * (px * x + py * y)/pt;
    //double sign_r = (x*px + y*py >= 0) ? 1. : -1.;
    double q = r*r - barrel.r_inner * barrel.r_inner;
    double distance_to_barrel_r = -p/2 + sqrt(p*p/4 - q);
    
    q = r*r - endcap.r_inner * endcap.r_inner;
    double distance_to_endcap_r = -p/2 + sqrt(p*p/4 - q);
    
    double sign_pz = (pz >= 0) ? 1. : -1.;
    double distance_to_barrel_z = barrel.delta_z - sign_pz*z;
    double distance_to_endcap_z = endcap.z_0 - sign_pz*z;
    
    //case 1: barrel only
    if (pt_over_pz > (distance_to_barrel_r + barrel.delta_r)/distance_to_barrel_z){ 
        length = rel_X0 * barrel.delta_r + distance_to_barrel_r * sqrt(1. + pow(1./pt_over_pz,2));
    }
    //case 2: barrel + endcap hit
    else if(pt_over_pz > distance_to_barrel_r/distance_to_barrel_z){
        //x == path traversed in barrel, rotation symmetry is still assumed at this point which is a valid approximation most of the times
        double x = (distance_to_barrel_z - distance_to_barrel_r/pt_over_pz)*sqrt(1+pow(pt_over_pz,2));
        //case 2a: traversed path in the barrel is larger than defined interaction path --> case 1
        if (x>=rel_X0*barrel.delta_r)
            length = rel_X0 * barrel.delta_r + distance_to_barrel_r * sqrt(1. + pow(1./pt_over_pz,2));
        //case 2b: particle is not absorbed in barrel but reaches the endcap --> length as in case 3 minus x
        else{
            length = (rel_X0 - x / barrel.delta_r) * endcap.delta_z + distance_to_endcap_z * sqrt(1. + pow(pt_over_pz,2));
        }
        //case 2c: distance from z-axis exceeds endcap extension (e.g. if particle travels through gab)
        double length_r = length * pt_over_pz/sqrt(1. + pow(pt_over_pz,2));
        if (length_r > (distance_to_endcap_r + endcap.delta_r)){
          length = length * (distance_to_endcap_r + endcap.delta_r)/length_r;
        }
    }
    //case 3: endcap only
    else if(pt_over_pz > distance_to_endcap_r/distance_to_endcap_z){ 
        length = rel_X0 * endcap.delta_z + distance_to_endcap_z * sqrt(1. + pow(pt_over_pz,2));
    }
    //case 4: part of endcap hit
    else if(pt_over_pz > distance_to_endcap_r/(distance_to_endcap_z + endcap.delta_z)){
        double x = (distance_to_endcap_r/pt_over_pz - distance_to_endcap_z)*sqrt(1+pow(pt_over_pz,2));
        //case 4a: traversed path in endcap is larger than defined interaction path --> case 3
        if (x>=rel_X0*endcap.delta_z)
            length = rel_X0 * endcap.delta_z + distance_to_endcap_z * sqrt(1. + pow(pt_over_pz,2));
        //case 4b: particle is not fully absorbed the endcap --> draw up to the yoke
        else
            length = (distance_to_endcap_z + endcap.delta_z) * sqrt(1. + pow(pt_over_pz,2));
    }
    //case 5: neither the endcap nor a barrel is hit
    else{  
        length = (distance_to_endcap_z + endcap.delta_z) * sqrt(1. + pow(pt_over_pz,2));
    }
    return fabs(length);
}