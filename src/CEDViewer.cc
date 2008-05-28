/* -*- Mode: C++; tab-width: 2; indent-tabs-mode: nil; c-basic-offset: 2 -*- */
#ifdef USE_CLHEP  // only if CLHEP is available !
#include "CEDViewer.h"
#include <iostream>

#include <EVENT/LCCollection.h>
#include <EVENT/Cluster.h>
#include <EVENT/Track.h>
#include <EVENT/MCParticle.h>
#include <EVENT/SimCalorimeterHit.h>
#include <EVENT/SimTrackerHit.h>

#include <UTIL/LCTypedVector.h>

// #include <ced_cli.h>
#include "MarlinCED.h"
#include "CLHEP/Vector/ThreeVector.h"

// fix for transition from CLHEP 1.8 to 1.9
namespace CLHEP{}
using namespace CLHEP ;


#include <math.h>
#include <cmath>
#include <cstdlib>

#include <marlin/Global.h>
#include <gear/GEAR.h>
#include <gear/BField.h>
#include <gear/TPCParameters.h>
#include <gear/PadRowLayout2D.h>

using namespace lcio ;
using namespace marlin ;


CEDViewer aCEDViewer ;


CEDViewer::CEDViewer() : Processor("CEDViewer") {
  
  // modify processor description
  _description = "CEDViewer: event display of LCIO objects  - based on CED by A.Zhelezov." ;
  
  StringVec cluExample ;
  cluExample.push_back("DrawCollection" ) ;
  cluExample.push_back("0" ) ;
  cluExample.push_back("3" ) ;
  
  registerOptionalParameter( "DrawCollection" , 
                             "collection to be displayed ( ColName, marker type[0-2], size)"  ,
                             _drawCollections ,
                             cluExample ,
                             cluExample.size() ) ;
  
  
}

void CEDViewer::init() { 

  // usually a good idea to
  printParameters() ;

  MarlinCED::init(this) ;
//   ced_client_init("localhost",7286);
//   ced_register_elements();
  
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


  _nRun = 0 ;
  _nEvt = 0 ;
  
}

void CEDViewer::processRunHeader( LCRunHeader* run) { 
  _nRun++ ;
} 


void CEDViewer::processEvent( LCEvent * evt ) { 

//-----------------------------------------------------------------------
// Reset drawing buffer and START drawing collection
  MarlinCED::newEvent(this) ;
//   ced_new_event();  
//-----------------------------------------------------------------------

  const gear::TPCParameters& gearTPC = Global::GEAR->getTPCParameters() ;
  const gear::PadRowLayout2D& padLayout = gearTPC.getPadLayout() ;

  if( parameterSet( "DrawCollection" ) ) {
    

    unsigned index = 0 ;
    while( index < _drawCollections.size() ){
      

      std::string colName(   _drawCollections[ index++ ] ) ; 
      int marker( std::atoi( _drawCollections[ index++ ].c_str() ) ) ; 
      int size  ( std::atoi( _drawCollections[ index++ ].c_str() ) ) ; 
      
      
      LCCollection* col = 0 ;
      try{

        col = evt->getCollection( colName ) ;

      }catch(DataNotAvailableException &e){
        // if collection doesn't exist go to next
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
	  
          int ml = marker | ( 7 << CED_LAYER_SHIFT ) ;
          int color =  _colors[ i % _colors.size() ] ;
          for( CalorimeterHitVec::const_iterator it = hits.begin();  it != hits.end() ; it++ ) {
	    
            ced_hit( (*it)->getPosition()[0],
                     (*it)->getPosition()[1],
                     (*it)->getPosition()[2],
                     ml, size , color ) ;
	    
          } // hits
	  
          float x = clu->getPosition()[0] ;
          float y = clu->getPosition()[1] ;
          float z = clu->getPosition()[2] ;

          ced_hit( x,y,z, ml, size*3 , color ) ;

          Hep3Vector v(x,y,z) ;

          Hep3Vector d(1.,1.,1.) ;


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
 
          Hep3Vector dp( v + d ) , dm( v - d )   ;

          ced_line( dp.x() , dp.y() , dp.z(),  
                    dm.x() , dm.y() , dm.z(),
                    ml , 1 , color );	 
	  

        } // cluster
	
      } else if( col->getTypeName() == LCIO::TRACK ){

        for( int i=0 ; i< col->getNumberOfElements() ; i++ ){
	  
          Track* trk = dynamic_cast<Track*>( col->getElementAt(i) ) ;
          const TrackerHitVec& hits = trk->getTrackerHits() ;
	  
          int ml = marker | ( 8 << CED_LAYER_SHIFT ) ;
	  
          for( TrackerHitVec::const_iterator it = hits.begin();  it != hits.end() ; it++ ) {
	    
            ced_hit( (*it)->getPosition()[0],
                     (*it)->getPosition()[1],
                     (*it)->getPosition()[2],
                     ml , size , _colors[ i % _colors.size() ] ) ;
	    
          } // hits
	  
          // draw the helix:
          double bField = Global::GEAR->getBField().at(  gear::Vector3D(0,0,0)  ).z() ; 
          //double bField = gearTPC.getDoubleVal("tpcBField") ;


          double pt = bField * 3e-4 / std::abs( trk->getOmega() ) ;
          double charge = ( trk->getOmega() > 0. ?  1. : -1. ) ;
	 
          double px = pt * std::cos(  trk->getPhi() ) ;
          double py = pt * std::sin(  trk->getPhi() ) ;
          double pz = pt * trk->getTanLambda() ;


          double rx = trk->getReferencePoint()[0] ;
          double ry = trk->getReferencePoint()[1] ;
          double rz = trk->getReferencePoint()[2] ;

          ml = marker | ( 9 << CED_LAYER_SHIFT ) ;

          MarlinCED::drawHelix( bField , charge, rx, ry, rz , 
                                px, py, pz, ml , size ,  0xffffff ,
                                0.0, padLayout.getPlaneExtent()[1]+100. , 
                                gearTPC.getMaxDriftLength()+100. ) ;


        } // track

      } else if( col->getTypeName() == LCIO::MCPARTICLE ){
	
        streamlog_out( DEBUG ) << "  drawing MCParticle collection " << std::endl ;

        for(int i=0; i<col->getNumberOfElements() ; i++){
	  
          MCParticle* mcp = dynamic_cast<MCParticle*> ( col->getElementAt( i ) ) ;
	  
          float charge = mcp->getCharge (); 
	  
          streamlog_out( DEBUG ) << "  drawing MCParticle pdg " 
                                 << mcp->getPDG() 
                                 << " genstat: " << mcp->getGeneratorStatus() 
                                 << std::endl ;


          //if( mcp-> getGeneratorStatus() != 1 ) continue ; // stable particles only   
          //	  if( mcp-> getSimulatorStatus() != 0 ) continue ; // stable particles only   
          if( mcp->getDaughters().size() > 0  ) continue ;    // stable particles only   
          // FIXME: need definition of stable particles (partons, decays in flight,...)
	  
          if ( mcp->getEnergy() < 0.001 ) continue ;           // ECut ?
	  
          double px = mcp->getMomentum()[0]; 
          double py = mcp->getMomentum()[1]; 
          double pz = mcp->getMomentum()[2];


          double x = mcp->getVertex()[0] ;
          double y = mcp->getVertex()[1] ;
          double z = mcp->getVertex()[2] ;	  

          if( std::fabs( charge ) > 0.0001  ) { 
	    

            double bField = Global::GEAR->getBField().at(  gear::Vector3D(0,0,0)  ).z() ; 
            //gearTPC.getDoubleVal("tpcBField") ;

            //            ml = marker | ( 7 << CED_LAYER_SHIFT ) ;

          streamlog_out( DEBUG ) << "  drawing MCParticle helix for p_t " 
                                 << sqrt(px*px+py*py)
                                 << std::endl ;

            MarlinCED::drawHelix( bField , charge, x, y, z, 
                                  px, py, pz, marker , size , 0x7af774  ,
                                  0.0,  padLayout.getPlaneExtent()[1]+100. ,
                                  gearTPC.getMaxDriftLength()+100. ) ;	    
	    
          } else { // neutral
	    
            int color  ;
            double r_min = 300 ;
            double z_max ;
            double r_max ;
	    
            if(std::abs( mcp->getPDG() )==22){
	      
              color = 0xf9f920;  // photon
              r_max = 1800 ;  // somewhere in the ecal
              z_max = 2930 ;  // somewhere in the ecal
	      
            } else {
		
              color = 0xb900de  ;  // neutral hadron
              r_max = 2500 ;  // somewhere in the hcal
              z_max = 3600 ;  // somewhere in the hcal
            } 
	    
            double pt = hypot(px,py);	 
            double p = hypot(pt,pz); 
	    
            double length = ( std::abs( pt/pz) > r_max/z_max ) ?  // hit barrel or endcap ? 
              r_max * p / pt  :  z_max * p / pz  ;
	    
            ced_line( r_min*px/p ,  r_min*py/p ,  r_min*pz/p , 
                      length*px/p ,  length*py/p ,  length*pz/p , 
                      marker , size, color );	 
	    
          }
        }
      } else if( col->getTypeName() == LCIO::SIMTRACKERHIT ){

        int color = 0xff00ff ;
        int layer = 2 ;
        LCTypedVector<SimTrackerHit> v( col ) ;
        MarlinCED::drawObjectsWithPosition( v.begin(), v.end() , marker, size , color, layer) ;

      } else if( col->getTypeName() == LCIO::SIMCALORIMETERHIT ){

        int color = 0xff0000 ;
        int layer = 3 ;
        LCTypedVector<SimCalorimeterHit> v( col ) ;
        MarlinCED::drawObjectsWithPosition( v.begin(), v.end() , marker, size , color, layer ) ;

      } else if( col->getTypeName() == LCIO::TRACKERHIT ){

        int color = 0xee0044 ;
        int layer = 4 ;
        LCTypedVector<TrackerHit> v( col ) ;
        MarlinCED::drawObjectsWithPosition( v.begin(), v.end() , marker, size , color, layer) ;

      } else if( col->getTypeName() == LCIO::CALORIMETERHIT ){

        int color = 0xee0000 ;
        int layer = 5 ;
        LCTypedVector<CalorimeterHit> v( col ) ;
        MarlinCED::drawObjectsWithPosition( v.begin(), v.end() , marker, size , color, layer ) ;

      }



    } // while

  }
  
  

//++++++++++++++++++++++++++++++++++++
  MarlinCED::draw(this) ;
//++++++++++++++++++++++++++++++++++++
    
  _nEvt ++ ;
}



void CEDViewer::check( LCEvent * evt ) { 
  // nothing to check here - could be used to fill checkplots in reconstruction processor
}


void CEDViewer::end(){ 
  
  streamlog_out(DEBUG) << "end() :" << " processed " << _nEvt 
                       << " events in " << _nRun << " runs "
                       << std::endl ;
  
}

#endif
