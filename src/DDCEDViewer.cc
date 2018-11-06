/* -*- Mode: C++; tab-width: 2; indent-tabs-mode: nil; c-basic-offset: 2 -*- */
#include "DDCEDViewer.h"

#include <iostream>
#include <set>

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
#include "DD4hep/DD4hepUnits.h" 

#include "ColorMap.h"
#include "TVector3.h"

using namespace lcio ;
using namespace marlin ;

using dd4hep::Position;

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

    registerProcessorParameter( "ColorByEnergy" ,
                               "color recunstructed particle by energy",
                               _colorEnergy,
                               bool(false) ) ;

    registerProcessorParameter( "ColorByEnergyMin" ,
                               "Minimal value for energy which will be represented as blue",
                               _colorEnergyMin,
                               double(0.0) ) ;

    registerProcessorParameter( "ColorByEnergyMax" ,
                               "Maximal value for energy which will be represented as red",
                               _colorEnergyMax,
                               double(35.0) ) ;

    registerProcessorParameter( "ColorByEnergySaturation" ,
                               "Hue value that will be used to determine the pallete",
                               _colorEnergySaturation,
                               double(0.8) ) ;

    registerProcessorParameter( "ColorByEnergyBrightness" ,
                               "Brigtness value that will be used to determine the pallete",
                               _colorEnergyValue,
                               double(0.8) ) ;

    registerProcessorParameter( "ColorByEnergyAutoColor" ,
                               "Automatically adjust event by event the blue to min energy and red to max energy of event",
                               _colorEnergyAuto,
                               bool(false) ) ;

    registerProcessorParameter( "ScaleLineThickness" ,
                               "Scale the line thickness of drawn helixes",
                               _scaleLineThickness,
                               double(1.0) ) ;

    registerProcessorParameter( "ScaleMarkerSize" ,
                               "Scale the size of the markes",
                               _scaleMarkerSize,
                               double(1.0) ) ;
    
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

    registerProcessorParameter( "DrawEllipsoidForPFOClusters" ,
                                "draw an ellipsoid for the Clusters of PFO objects: 0: disabled, 1: enabled",
                                _drawEllipsoidForPFOClusters,
                                0  ) ;
        

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

void DDCEDViewer::processRunHeader( LCRunHeader* /*run*/) {
    _nRun++ ;
}

void DDCEDViewer::processEvent( LCEvent * evt ) {
    //-----------------------------------------------------------------------
    // Reset drawing buffer and START drawing collection

    DDMarlinCED::newEvent(this) ;
    
    //added by Thorben Quast for removing GEAR dependencies
    dd4hep::Detector& theDetector = dd4hep::Detector::getInstance();

    DDMarlinCED::drawDD4hepDetector(theDetector, _surfaces, _detailled);
    DDCEDViewer::drawDD4LCIO(evt, theDetector);

    DDMarlinCED::draw(this, _waitForKeyboard );

    //++++++++++++++++++++++++++++++++++++

    _nEvt ++ ;
}



void DDCEDViewer::check( LCEvent * /*evt*/ ) {
    // nothing to check here - could be used to fill checkplots in reconstruction processor
}

void DDCEDViewer::printParticle(int id, LCEvent * /*evt*/){
    streamlog_out( MESSAGE )  << "CEDViewer::printParticle id: " << id << std::endl;
} 


void DDCEDViewer::end(){ 
    streamlog_out(DEBUG2 ) << "end() :" << " processed " << _nEvt 
    << " events in " << _nRun << " runs "
    << std::endl ;
}

void DDCEDViewer::drawDD4LCIO(LCEvent * evt, dd4hep::Detector& theDetector){
DDCEDPickingHandler &pHandler=DDCEDPickingHandler::getInstance();
    
    pHandler.update(evt);
    
    //-----------------------------------------------------------------------
    if( _useTrackerForLimitsOfHelix ){
        //only draw the helix of charged tracks (other than muons) in the TPC within the tracker volume
        _helix_max_r = getTrackerExtent(theDetector)[0];
        _helix_max_z = getTrackerExtent(theDetector)[1];
    }
    // //add DrawInLayer to
    // //std::vector< DrawParameters > drawParameters ;
    // this->drawParameters.reserve( _drawCollections.size() + _drawCollectionsLayer.size()  ) ;
    
    // if( parameterSet( "DrawCollection" ) ) {
        
    //     unsigned index = 0 ;
    //     while( index < _drawCollections.size() ){
            
    //         const std::string & colName = _drawCollections[ index++ ] ;
    //         int marker = std::atoi( _drawCollections[ index++ ].c_str() ) ;
    //         int size = std::atoi( _drawCollections[ index++ ].c_str() ) ;
    //         //std::cout << "#######################      SIZE: " << size << " #######################################" << std::endl;
    //         int layer = -1 ;
            
    //         this->drawParameters.push_back(DrawParameters( colName,size,marker,layer ) );
    //     }
    // }
    // if( parameterSet( "DrawInLayer" ) ) {
        
    //     unsigned index = 0 ;
    //     while( index < _drawCollectionsLayer.size() ){
            
    //         const std::string & colName = _drawCollectionsLayer[ index++ ] ;
    //         int marker = std::atoi( _drawCollectionsLayer[ index++ ].c_str() ) ;
    //         int size   = std::atoi( _drawCollectionsLayer[ index++ ].c_str() ) ;
    //         int layer  = std::atoi( _drawCollectionsLayer[ index++ ].c_str() ) ;
            
    //         this->drawParameters.push_back(DrawParameters( colName,size,marker,layer ) );
    //         //std::cout << "layer: " << layer << " description: " << colName << std::endl; //hauke
    //     }
    // }
    
    std::set< std::string > colsInEvent ;

    unsigned nCols = this->drawParameters.size() ;
    //draw the individual collections as indicated in DrawInLayer
    for(unsigned np=0 ; np < nCols ; ++np){
        
        const std::string & colName = this->drawParameters[np].ColName ;
        int size =   this->drawParameters[np].Size ;
        int marker = this->drawParameters[np].Marker ;
        int layer =  this->drawParameters[np].Layer ;
        
        LCCollection* col = 0 ;
        try{
            col = evt->getCollection( colName ) ;

            colsInEvent.insert( colName ) ;

        }catch(DataNotAvailableException &e){
            streamlog_out( DEBUG5 ) << " collection " << colName <<  " not found ! "   << std::endl ;
            continue ;
        }
        
        if( colName.find("Jet") != std::string::npos ){
            //Usually jet collections contain the substring "Jet", e.g. "JetOut", "Durham_XJets", ...
            DDCEDViewer::drawJets(theDetector, layer, colName, col);
        } else if( col->getTypeName() == LCIO::CLUSTER ){
            DDCEDViewer::drawCluster(theDetector, layer, np, colName, marker, col, size);
        } else if( col->getTypeName() == LCIO::TRACK ){
            DDCEDViewer::drawTrack(theDetector, layer, np, colName, marker, col, size);
        } else if( col->getTypeName() == LCIO::MCPARTICLE ){
            DDCEDViewer::drawMCParticle(theDetector, layer, np, colName, marker, col, size);
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
            DDCEDViewer::drawReconstructedParticle(theDetector, layer, np, colName, marker, col, size);
        }    
    }

    streamlog_out( MESSAGE ) << " ++++++++ collections shown on layers:  +++++++++++++ " << std::endl ;
    
    for(unsigned np=0 ; np < nCols ; ++np){
        
        const std::string & colName = drawParameters[np].ColName ;
        //     int size =   drawParameters[np].Size ;
        //     int marker = drawParameters[np].Marker ;
        int layer =  drawParameters[np].Layer ;
        
        if( colsInEvent.find( colName ) != colsInEvent.end() ) {
          streamlog_out( MESSAGE )  << "    +++++  " << colName <<  "\t  on layer: " << layer << std::endl ;
        }
    }

    streamlog_out( MESSAGE ) << std::endl
                             <<  " [ evt: " << evt->getEventNumber()
                             <<  "   run: " << evt->getRunNumber() << " ] " << std::endl ;

    streamlog_out( MESSAGE ) << " ++++++++ use shift-[LN] for LN>10  +++++++++++++ " << std::endl ;
    


}

void DDCEDViewer::drawCluster(dd4hep::Detector& /*theDetector*/, int& layer, unsigned& np, std::string /*colName*/, int& marker, LCCollection* col, int& size){
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
                       marker, layer, size*_scaleMarkerSize , color, clu->id() ) ;
        } 
        float x = clu->getPosition()[0] ;
        float y = clu->getPosition()[1] ;
        float z = clu->getPosition()[2] ;
        ced_hit_ID( x,y,z, marker, layer , size*3*_scaleMarkerSize , color, clu->id() ) ;
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

void DDCEDViewer::drawTrack(dd4hep::Detector& theDetector, int& layer, unsigned& np, std::string colName, int& marker, LCCollection* col, int& size){
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
                       marker, layer , size*_scaleMarkerSize , color , trk->id() ) ;
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
                           1 , layer , size*10*_scaleMarkerSize , color , trk->id() ) ;
                break ;
        }
        if( ts !=0 ){
            double* bFieldVector = new double[3];
            theDetector.field().combinedMagnetic(Position(0,0,0), bFieldVector) ;
            double bField = bFieldVector[2] / dd4hep::tesla;
            delete[] bFieldVector;
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
                                     px, py, pz, ml ,  2*_scaleLineThickness , helixColor  ,
                                     0.0, _helix_max_r,
                                     _helix_max_z, trk->id() ) ;
            }
        }
    }
}

void DDCEDViewer::drawMCParticle(dd4hep::Detector& theDetector, int& layer, unsigned& np, std::string colName, int& marker, LCCollection* col, int& size) {
    streamlog_out( DEBUG ) << "  drawing MCParticle collection " << std::endl ;
    //new line drawing implemented by Thorben Quast 07 August 2015
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
            theDetector.field().combinedMagnetic(Position(0,0,0), bFieldVector) ;
            double bField = bFieldVector[2] / dd4hep::tesla;
            delete[] bFieldVector;
            streamlog_out( DEBUG ) << "  drawing MCParticle helix for p_t "
            << sqrt(px*px+py*py)
            << std::endl ;
            const int ml = marker | ( layer << CED_LAYER_SHIFT ) ;
            //maximal extension of all charged tracks
            double _hmr, _hmz;
            switch(std::abs(mcp->getPDG() ) ){
                case 13:
                    _hmr = getCalorimeterParameters(theDetector, "HCalBarrel").r_inner + getCalorimeterParameters(theDetector, "HCalBarrel").delta_r; 
                    _hmz = getCalorimeterParameters(theDetector, "HCalEndcap").z_0 + getCalorimeterParameters(theDetector, "HCalEndcap").delta_z;
                    break;
                default:
                    _hmr = _helix_max_r; 
                    _hmz = _helix_max_z;
            }
            DDMarlinCED::drawHelix( bField , charge, x, y, z,
                                 px, py, pz, ml , size*_scaleLineThickness , 0x7af774  ,
                                 0.0,  _hmr ,
                                 _hmz, mcp->id() ) ; 
        } else { // neutral
            int color  ;
            double length, yokeR, yokeZ;
            switch(  std::abs(mcp->getPDG() )  ){
                //refactored length calculation (T. Quast 7 Aug 15) 
                case 22:
                    color = 0xf9f920;          // photon
                    length = calculateTrackLength("ecal", theDetector, x, y, z, px, py, pz);
                    break ;
                case 12:  case 14: case 16: // neutrino
                    color =  0xdddddd  ;
                    yokeR = getYokeExtent(theDetector)[0];
                    yokeZ = getYokeExtent(theDetector)[1];
                    length = (fabs(pt/pz) > yokeR/yokeZ) ?
                            yokeR * sqrt(1. + pow(pz/pt,2)):
                            yokeZ * sqrt(1. + pow(pt/pz,2));
                    break ;
                default:
                    color = 0xb900de  ;        // neutral hadron
                    length = calculateTrackLength("hcal", theDetector, x, y, z, px, py, pz);
            }
            //tracks with vertex outside the according calorimeter are not drawn, length is passed as 0
            ced_line_ID( x , y , z ,
                        x + length*px/p ,  y + length*py/p ,  z + length*pz/p ,
                        layer  , size, color, mcp->id() );
        }
    }
}

void DDCEDViewer::drawSIMTrackerHit(int& layer, unsigned& np, std::string colName, int& marker, LCCollection* col, std::vector<int>& _Colors, int& size){
    layer = ( layer > -1 ? layer : SIMTRACKERHIT_LAYER ) ;
    this->drawParameters[np].Layer = layer ;
    DDMarlinCED::add_layer_description(colName, layer);

    for( int i=0, n=col->getNumberOfElements(); i<n ; i++ ){
        
        SimTrackerHit* h = dynamic_cast<SimTrackerHit*>( col->getElementAt(i) ) ;
        
        // color code by MCParticle
        const int mci = ( h->getMCParticle() ? h->getMCParticle()->id() :  0 ) ;
        int color = _Colors[  mci % _Colors.size() ] ;
        
        int id =    h->id();
        ced_hit_ID( h->getPosition()[0],
                   h->getPosition()[1],
                   h->getPosition()[2],
                   marker,layer, size*_scaleMarkerSize , color, id ) ;
        
    }
}

void DDCEDViewer::drawSIMCalorimeterHit(int& layer, unsigned& np, std::string colName, int& marker, LCCollection* col, std::vector<int>& _Colors, int& size){
    layer = ( layer > -1 ? layer : SIMCALORIMETERHIT_LAYER ) ;
    this->drawParameters[np].Layer = layer ;
    
    DDMarlinCED::add_layer_description(colName, layer);
    
    for( int i=0, n=col->getNumberOfElements(); i<n ; i++ ){
        
        SimCalorimeterHit* h = dynamic_cast<SimCalorimeterHit*>( col->getElementAt(i) ) ;
        
        // color code by MCParticle
        //        const int mci = (  h->getNMCContributions() !=0  ?  h->getParticleCont(0)->id()  :  0 ) ;
        const int mci = (  h->getNMCContributions() !=0  && h->getParticleCont(0) ?  h->getParticleCont(0)->id()  :  0 ) ;
        int color = _Colors[  mci % _Colors.size() ] ;
        
        int id =    h->id();
        ced_hit_ID( h->getPosition()[0],
                   h->getPosition()[1],
                   h->getPosition()[2],
                   marker,layer, size*_scaleMarkerSize , color, id ) ;
        
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
            ced_hit_ID( p.X(), p.Y(), p.Z(), marker, layer , size*_scaleMarkerSize , color, h->id() );
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

void DDCEDViewer::drawReconstructedParticle(dd4hep::Detector& theDetector, int& layer, unsigned& np, std::string colName, int& marker, LCCollection* col, int& size){
    //hauke
    layer = ( layer > -1 ? layer : RECOPARTICLE_LAYER ) ;
    this->drawParameters[np].Layer = layer ;
    
    DDMarlinCED::add_layer_description(colName, layer);
    int nelem = col->getNumberOfElements();
    
    float TotEn = 0.0;
    float TotPX = 0.0;
    float TotPY = 0.0;
    float TotPZ = 0.0;
    
    //Determine the maximal and minimal cluster energy depositions in the event for color scaling (-->when drawing ellipsoids/cylinders).
    double Emin = 99999.; double Emax = 0;
    double pEmin = 99999.; double pEmax = 0;
    for (int ip(0); ip < nelem; ++ip) {
        ReconstructedParticle * part = dynamic_cast<ReconstructedParticle*>(col->getElementAt(ip));
        ClusterVec clusterVec = part->getClusters();
        float pene = part->getEnergy();
        pEmin = fmin(pEmin, pene);
        pEmax = fmax(pEmax, pene);
        unsigned nClusters = (unsigned)clusterVec.size();
        if (nClusters > 0 ) {
            for (unsigned int p=0; p<nClusters; p++) {
                double e = clusterVec[p]->getEnergy();
                Emin = fmin(Emin, e);
                Emax = fmax(Emax, e);
            }   
        }
    }
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
        
        if( _colorEnergy ){
          if( _colorEnergyAuto ){
            color = ColorMap::NumberToTemperature(ene,pEmin,pEmax,_colorEnergySaturation,_colorEnergyValue);
          }else{
            color = ColorMap::NumberToTemperature(ene,_colorEnergyMin,_colorEnergyMax,_colorEnergySaturation,_colorEnergyValue);
          }
        }

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
              ced_hit_ID(x,y,z,marker, layer ,size*_scaleMarkerSize,color,part->id());
            }
          }

          //fg: this code below needs some work: the ellipsoids are flat and they should use the cluster
          //    shower parameters really, also this should probably be drawn in a different layer ....
          if( _drawEllipsoidForPFOClusters ) {
            //refactored Cluster drawing as ellipsoids 
            //by Thorben Quast, CERN Summer Student 2015
            //18 August 2015
            for (unsigned int p=0; p<nClusters; p++) {
              //Energy clusters are drawn as ellipsoids.
              //For each cluster, it's (energy weighted) central position, the deposited energy and the intrinsic direction in terms of sperical angles are given.
              //The minimal and maximal deposited energies among all clusters in the displayed event have been determined previously and will be needed for coloring.
              Cluster * cluster = clusterVec[p];
              double cluster_center[] = {cluster->getPosition()[0], cluster->getPosition()[1], cluster->getPosition()[2]};
              double phi = cluster->getIPhi();
              double theta = cluster->getITheta();
              //Use the direction of the cluster center w.r.t. the origin if no intrinsic angles are given.
              if( phi ==0. && theta==0.){
                theta = atan( sqrt( cluster_center[0]*cluster_center[0] + cluster_center[1]*cluster_center[1] ) / cluster_center[2]  ) ;
                phi = atan2( cluster_center[1] , cluster_center[0] ) ;
              }
              //Energy weighted moments of inertia are calculated. Ultimately, the eigenvalues of the 3x3 matrix will be a measure of the ellipsoids' extensions.
              CalorimeterHitVec hitvec = cluster->getCalorimeterHits();
              int nHits = (int)hitvec.size();
              double Etot = 0;                
              double I[3][3]; for(int i=0; i<3; i++) for(int j=0; j<3; j++) I[i][j] = 0;
              //The angles theta and phi are used to transform the coordinates of each hit into a c.s. in which the x'-axis is parrallel to
              //the ellisoid's intrinsic direction. Note that this does not describe an unambigious system as any rotation along the x'-axis does not touch this constraint.
              //The transformation is achieved by a typical multiplication of rotation matrices: R(theta, phi) = R_y(Pi/2 - theta)*R_z(phi)
              double R[3][3]; 
              R[0][0] = cos(phi) * sin(theta); R[1][0] = -sin(phi); R[2][0] = -cos(phi)*cos(theta); R[0][1] = sin(phi)*sin(theta); 
              R[1][1] = cos(phi); R[2][1] = -sin(phi)*cos(theta); R[0][2] = cos(theta); R[1][2] = 0; R[2][2] = sin(theta);
              double tot_x =0;
              double tot_y =0;
              double tot_z =0;
              for (int q = 0; q < nHits; q++){
                CalorimeterHit * hit = hitvec[q];
                float x = hit->getPosition()[0];
                float y = hit->getPosition()[1];
                float z = hit->getPosition()[2];
                float e = hit->getEnergy();
                //ced_hit_ID(x,y,z,marker, layer ,size,color,part->id());   //this line draws the indivdual hits within a cluster
                //translation and rotation of the coordinates
                x -= cluster_center[0];
                y -= cluster_center[1];
                z -= cluster_center[2];
                tot_x += x*e;
                tot_y += y*e;
                tot_z += z*e;
                
                double new_x = x * R[0][0] + y * R[0][1] + z * R[0][2];
                double new_y = x * R[1][0] + y * R[1][1] + z * R[1][2];
                double new_z = x * R[2][0] + y * R[2][1] + z * R[2][2];
                x = new_x; y = new_y; z = new_z;
                //calculate moments of inertia
                I[0][0] += x*x*e; I[1][1] += y*y*e; I[2][2] += z*z*e; 
                I[0][1] = I[1][0] += x*y*e; I[0][2] = I[2][0] += x*z*e; I[1][2] = I[2][1] += y*z*e;
                Etot += e;
              }
              //The result of the rotation by the matrix R, the following coordinates correspond with each other:
              //  (component in coordinate system with x' || intrinsic direction)       (system in which ellipsoids are initially placed)
              //                          x'                                      <-->        z       
              //                          y'                                      <-->        y       
              //                          z'                                      <-->       -x   
              //These assignments are corrected for by a modified rotation of the ellipsoid along its y-axis (see declaration of double rotate[])
              
              //I is not diagonal yet as only one axis was fixed when applying the rotation R.
              //The remaining lengths are determined by the solution of the 2x2 Eigenvalues (p-q formula).
              double lambda1 = 0.5*(I[2][2]+I[1][1]) + sqrt( pow(0.5*(I[2][2]+I[1][1]),2 ) + pow(I[2][1],2)-I[2][2]*I[1][1]);
              double lambda2 = 0.5*(I[2][2]+I[1][1]) - sqrt( pow(0.5*(I[2][2]+I[1][1]),2 ) + pow(I[2][1],2)-I[2][2]*I[1][1]);
              double sizes[3];
              sizes[0] = I[0][0]; sizes[1] = lambda1; sizes[2] = lambda2;
              //Remaining: (more or less) Arbitrary rescaling and transformation to a physical length (sqrt + energy division)
              for (int i=0; i<3; i++)  sizes[i] = sqrt(17.727)*sqrt(sizes[i])/Etot;
              
              double alpha = 0.5*asin(2*I[1][2]/fabs(lambda1-lambda2)) * (I[1][1]-I[2][2])/fabs(I[1][1]-I[2][2]);
              //We want to rotate the ellipsoid w.r.t. to the y-axis by 90deg - theta. Due to the rotation by R that maps x' <--> z, it is now
              //upside down such that an additional sign is needed.
              double rotate[] = {alpha, -(90-theta*180/M_PI), phi*180/M_PI};
              
              //The colors (blue and red) are set according the deposited energy in the cluster by comparison to other clusters in the event
              int ellipsoid_color = returnRGBClusterColor(cluster->getEnergy(), Emin, Emax, 256, 'a', 3);
              
              //Draw the ellipsoids, uncommenting the line with cylinders works as well.
              ced_ellipsoid_r(sizes, cluster_center, rotate, layer, ellipsoid_color); 
              //ced_geocylinder_r(0.25*(sizes[0]+sizes[1]), sizes[2], cluster_center, rotate, 36, ellipsoid_color, layer); 
            }
          } 
        }

        //************
        if (nTracks > 0) {
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
                        ced_hit_ID(x,y,z,marker, layer,size*_scaleMarkerSize,color,part->id());
                    }
                }
                if((nHits==0 || _drawHelixForPFOs==1) && std::fabs(part->getCharge())>0.001){
                    TrackState *ts=0;
                    if (trk->getTrackStates().size() > 0) {
                        ts = trk->getTrackStates().at(0) ;
                    }
                    if(ts!=0){
                        double* bFieldVector = new double[3];
                        theDetector.field().combinedMagnetic(Position(0,0,0), bFieldVector) ;
                        double bField = bFieldVector[2] / dd4hep::tesla;
                        delete[] bFieldVector;
                        double pt;
                        if (bField != 0.0 && std::abs(ts->getOmega()) > 0.00001 ){
                            pt = bField * 3e-4 / std::abs( ts->getOmega() ) ;
                        }else{
                            pt = 1.e10;
                        }
                        double charge = ( ts->getOmega() > 0. ?  1. : -1. ) ;
                        double Px = pt * std::cos(  ts->getPhi() ) ;
                        double Py = pt * std::sin(  ts->getPhi() ) ;
                        double Pz = pt * ts->getTanLambda() ;
                        double Xs = ts->getReferencePoint()[0] -  ts->getD0() * sin( ts->getPhi() ) ;
                        double Ys = ts->getReferencePoint()[1] +  ts->getD0() * cos( ts->getPhi() ) ;
                        double Zs = ts->getReferencePoint()[2] +  ts->getZ0() ;
                        int helixColor = ( _useColorForHelixTracks ? color : 0xdddddd ) ;
                        //helix
                        DDMarlinCED::drawHelix(bField, charge, Xs, Ys, Zs, Px, Py, Pz, marker|(layer<<CED_LAYER_SHIFT), size/2*_scaleLineThickness,
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

void DDCEDViewer::drawJets(dd4hep::Detector& theDetector, int layer, std::string colName, LCCollection* col){
    //default color is orange
    float RGBAcolor[4] = {.9, .7, .0, 0.25};
    int color = int(RGBAcolor[2]*(15*16+15)) + int(RGBAcolor[1]*(15*16+15))*16*16+ int(RGBAcolor[0]*(15*16+15))*16*16*16*16;

    //only one registration for all jets in CEDViewer
    DDMarlinCED::add_layer_description(colName, layer);

    for (int j=0; j < col->getNumberOfElements(); ++j) { //number of elements in a jet
        ReconstructedParticle * jet = dynamic_cast<ReconstructedParticle*>( col->getElementAt(j) );
        streamlog_out( DEBUG )  <<   "     - jet energy " << jet->getEnergy() << std::endl ;

        //total momentum of the jet
        TVector3 v(jet->getMomentum()[0], jet->getMomentum()[1], jet->getMomentum()[2]); 
        
        //init relevant objects for looping over all particles in the jet
        const ReconstructedParticleVec & pv = jet->getParticles();
        int N_elements = pv.size();
        float pt_tot = 0.0; float E_max = 0.0; float mean_tan_angle = 0.0;
        std::vector<TVector3> pp; std::vector<float> pt; std::vector<float> E;
        pp.reserve(N_elements); pt.reserve(N_elements); E.reserve(N_elements);

        //calculate longitudinal, transverse momentum (w.r. to jet axis) for each particle
        //from that deduce a pt-weighted mean (tan-) angle and determine the highest pt contribution
        for (int k = 0; k<N_elements; ++k){
            TVector3 pp_k(pv[k]->getMomentum()[0], pv[k]->getMomentum()[1], pv[k]->getMomentum()[2]);
            pp.push_back(pp_k);
            TVector3 ju = v.Unit();
            TVector3 pt_k = pp_k - (ju.Dot(pp_k))*ju;  
            TVector3 pl_k = pp_k - pt_k;
            pt.push_back(pt_k.Mag());
            pt_tot += pt[k];
            E.push_back(pv[k]->getEnergy());
            E_max = (E[k] > E_max) ? E[k]: E_max;
            mean_tan_angle += pt[k]*(pt[k]/pl_k.Mag());

        }
        mean_tan_angle /= pt_tot;

        //draw the line of movement for each particle in the jet
        for (int k = 0; k<N_elements; ++k){
            float center_ref[3] = {0., 0., 0.};
            //100% * distance = length holds for the entry with highest pt, the others obtain only a respective fraction
            double momLength = (E[k]/E_max)*calculateTrackLength("", theDetector, center_ref[0], center_ref[1], center_ref[2], pp[k].X(), pp[k].Y(), pp[k].Z());                    //line size
            //approximation: all lines start in origin (TODO, if jet origin known)
            ced_line_ID(center_ref[0], center_ref[1], center_ref[2], momLength*pp[k].X()/pp[k].Mag(), momLength*pp[k].Y()/pp[k].Mag(), momLength*pp[k].Z()/pp[k].Mag(), layer, 1, color, pv[k]->id());     
        }

        //calculate the parameters of the jet cone
        //approximation: all lines start in origin (TODO, if jet origin known)            
        double center_c[3] = {0., 0., 0. };
        double rotation_c[3] = { 0.,  v.Theta()*180./M_PI , v.Phi()*180./M_PI };
        double coneHeight = calculateTrackLength("", theDetector, center_c[0], center_c[1], center_c[2], v.X(), v.Y(), v.Z());
        //1. baseline radius, 2. height, 3. origin doublet, 4. rotation triplet,...
        ced_cone_r_ID( mean_tan_angle * coneHeight , coneHeight , center_c, rotation_c, layer, RGBAcolor,jet->id()); 
        
    }
}
/********************************************************************
//Helper functions, might be exported into a utility.cc at some point
********************************************************************/

//get the outer extents of the tracker
double* getTrackerExtent(dd4hep::Detector& theDetector){
    double* extent = new double[2];
    extent[0] =theDetector.constant<double>("tracker_region_rmax")/dd4hep::mm;
    extent[1] = theDetector.constant<double>("tracker_region_zmax")/dd4hep::mm;
    return extent;
}
//get the outer extents of the yoke
double* getYokeExtent(dd4hep::Detector& theDetector){
    double* extent = new double[2];
    const std::vector< dd4hep::DetElement>& calorimeters     = theDetector.detectors( "calorimeter" ) ;
    dd4hep::DetElement yoke;
    for( unsigned i=0,n=calorimeters.size() ; i<n ; ++i ){
        std::string detName = calorimeters[i].name();
        bool isYokeBarrel = (detName == "YokeBarrel") ;
        bool isYokeEndcap = (detName == "YokeEndcap") ;
        if (!isYokeBarrel && !isYokeEndcap)
            continue;
        yoke = calorimeters[i];
        dd4hep::rec::LayeredCalorimeterData* yokeGeo;
        try{ 
            yokeGeo = yoke.extension<dd4hep::rec::LayeredCalorimeterData>() ; 
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
CalorimeterDrawParams getCalorimeterParameters(dd4hep::Detector& theDetector, std::string name, bool selfCall){
    CalorimeterDrawParams params;
    if (selfCall)
        name[1] = tolower(name[1]); 
    const std::vector< dd4hep::DetElement>& calorimeters     = theDetector.detectors( "calorimeter" ) ;
    dd4hep::DetElement calo;
    for( unsigned i=0,n=calorimeters.size() ; i<n ; ++i ){
        if ((std::string)calorimeters[i].name() == name){
            calo = calorimeters[i];
            break;
        }
    }
    dd4hep::rec::LayeredCalorimeterData* caloGeo;
    try{ 
        caloGeo = calo.extension<dd4hep::rec::LayeredCalorimeterData>() ; 
    } catch(std::runtime_error& e){
            if (!selfCall)
                return getCalorimeterParameters(theDetector, name, true); 
            else{
                streamlog_out( MESSAGE ) <<  "MC Particles for " << name << " cannot be drawn"<<std::endl;
                params.delta_z = -1;    //no spatial extension --> no drawing
                return params;
            }
    }
    if (caloGeo->layoutType == dd4hep::rec::LayeredCalorimeterData::BarrelLayout){
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
double calculateTrackLength(std::string type, dd4hep::Detector& theDetector, double x, double y, double z, double px, double py, double pz){
    double rel_X0;
    CalorimeterDrawParams barrel; CalorimeterDrawParams endcap;
    if (type == "ecal"){
        barrel = getCalorimeterParameters(theDetector, "ECalBarrel");
        endcap = getCalorimeterParameters(theDetector, "ECalEndcap");
        rel_X0 = 0.5;
    }else if(type == "hcal"){
        barrel = getCalorimeterParameters(theDetector, "HCalBarrel");
        endcap = getCalorimeterParameters(theDetector, "HCalEndcap");
        rel_X0 = 0.5;
    }else{
        barrel = getCalorimeterParameters(theDetector, "ECalBarrel");
        endcap = getCalorimeterParameters(theDetector, "ECalEndcap");
        rel_X0 = 0.;
    }

    if (barrel.delta_z == -1 || endcap.delta_z == -1) return 0;   //the case if the parameters could not be loaded properly
    
    double length;
    double pt = sqrt(px*px + py*py);
    double pt_over_pz = pt/fabs(pz);

    double r = sqrt(x*x+y*y);
    if (r > barrel.r_inner || r > endcap.r_inner) return 0;
    double p = 2 * (px * x + py * y)/pt;

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
        double X = (distance_to_barrel_z - distance_to_barrel_r/pt_over_pz)*sqrt(1+pow(pt_over_pz,2));
        //case 2a: traversed path in the barrel is larger than defined interaction path --> case 1
        if (X>=rel_X0*barrel.delta_r)
            length = rel_X0 * barrel.delta_r + distance_to_barrel_r * sqrt(1. + pow(1./pt_over_pz,2));
        //case 2b: particle is not absorbed in barrel but reaches the endcap --> length as in case 3 minus x
        else{
            length = (rel_X0 - X / barrel.delta_r) * endcap.delta_z + distance_to_endcap_z * sqrt(1. + pow(pt_over_pz,2));
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
        double X = (distance_to_endcap_r/pt_over_pz - distance_to_endcap_z)*sqrt(1+pow(pt_over_pz,2));
        //case 4a: traversed path in endcap is larger than defined interaction path --> case 3
        if (X>=rel_X0*endcap.delta_z)
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

int returnRGBClusterColor(float eneCluster, float cutoff_min, float cutoff_max, int color_steps, char scale, int colorMap){
    int color = 0x000000; //default colour: black
    int color_delta = 0; //colour step in the [0, color_steps] spectrum
    unsigned int rgb[] = {0, 0, 0}; //array of RGB to be returned as one 0x000000 HEX value

    /**
    * Check the input values for sanity */
    if (cutoff_min > cutoff_max) {
        std::cout << "Error in 'DSTViewer::returnRGBClusterColor': cutoff_min < cutoff_max" << std::endl;
    }
    if (eneCluster < 0.0) {
        std::cout << "Error in 'DSTViewer::returnRGBClusterColor': eneCluster is negative!" << std::endl;
    }
    if (cutoff_min < 0.0) {
        std::cout << "Error in 'DSTViewer::returnRGBClusterColor': eneCluster is negative!" << std::endl;
    }
    if (colorMap < 0 || colorMap > 6) {
        std::cout << "Error in 'DSTViewer::returnRGBClusterColor': wrong colorMap param!" << std::endl;
    }
    // Input values in log-scale
    float log_ene = std::log(eneCluster+1);
    float log_min = std::log(cutoff_min+1);
    float log_max = std::log(cutoff_max+1);
    float log_delta = log_max - log_min;
    float log_step = log_delta/(float)color_steps;

    switch(scale){
        case 'a': default: //log
            color_delta = (int) ((log_ene-log_min)/log_step); // which colour bin does the value go to? We have [0x00,0xFF] bins
            break;
        case 'b': //linear
            color_delta = (int)((eneCluster - cutoff_min)/(cutoff_max - cutoff_min)*color_steps);
            break;
    }


    if (color_delta >= color_steps){
        color_delta = color_steps;
    }
    if (color_delta < 0){
        color_delta = 0;
    }

    ColorMap::selectColorMap(colorMap)(rgb, color_delta, 0, color_steps);
    color = ColorMap::RGB2HEX(rgb[0],rgb[1],rgb[2]);

    return color;
}



