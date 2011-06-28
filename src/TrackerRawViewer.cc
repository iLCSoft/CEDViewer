#include "TrackerRawViewer.h"

#include "MarlinCED.h"
#include "ColorMap.h"

#include <EVENT/LCCollection.h>
#include <EVENT/TrackerRawData.h>

//--- Marlin
//#include "marlin/VerbosityLevels.h"
#include "marlin/ConditionsProcessor.h"
#include "marlin/Global.h"


#include "ChannelPosMap.h"


//---- GEAR ----
#include "gear/TPCModule.h"
#include "gear/PadRowLayout2D.h"
#include "gear/BField.h"

using namespace lcio ;
using namespace marlin ;


TrackerRawViewer aTrackerRawViewer ;

// --- fix for lccd/ConditionsMap.hh tp print key made from pair
std::ostream& operator<<(  std::ostream& os , const std::pair<int,int>& p){
  os  << "[" << p.first << "," << p.second << "]" ; 
  return os ;
}
  

// --- prin TrackerRawData - need to go replace the one in LCIO ...
std::ostream& operator<<(  std::ostream& os , const TrackerRawData& rd){
  
  using namespace std ;

  os << setw(30) << setfill(' ') << left << "CellID0"<< setfill(' ') << right << setw(40) << dec << rd.getCellID0() << endl;
  os << setw(30) << setfill(' ') << left << "CellID1"<< setfill(' ') << right << setw(40) << dec << rd.getCellID1() << endl;
  os << setw(30) << setfill(' ') << left << "Time"   << setfill(' ') << right << setw(40) << dec << rd.getTime() << endl;

  const ShortVec& adc = rd.getADCValues() ;
	  
  os << setfill(' ') << left << "ADC values: ["    ;
  for(unsigned int j=0 ; j< adc.size() ; ++j ) {
    os  << setfill(' ') << right << setw(4) << dec << adc[j] << ","  ;
  }
  os << "]" << std::endl ;

  return os ;
}

void printTRD(const LCObject* rd){
  streamlog_out( MESSAGE ) << *(TrackerRawData* ) rd ;
}

//===============================================================================================================================

TrackerRawViewer::TrackerRawViewer() : Processor("TrackerRawViewer") ,
				       _tpcParams(0),
				       _posMap(0),
				       _chMap(0),
				       _nRun(0),
				       _nEvt(0) {
  
  // modify processor description
  _description = "TrackerRawViewer: visualization of TPC prototype raw data." ;
  

  // register steering parameters: name, description, class-variable, default value

  registerInputCollection( LCIO::TRACKERRAWDATA,
			   "CollectionName" , 
			   "Name of the TrackerRawData collection"  ,
			   _colName ,
			   std::string("AltroRawData") ) ;
  
  registerInputCollection( LCIO::TRACKERHIT,
			   "HitCollectionName" , 
			   "Name of the TrackerHit collection - optional"  ,
			   _hitColName ,
			   std::string("TPCTrackerHits") ) ;

  registerInputCollection( LCIO::TRACK,
			   "TrackCollectionName" , 
			   "Name of the Tracker collection - optional"  ,
			   _trkColName ,
			   std::string("TPCTracks") ) ;

  registerOptionalParameter( "ChannelMappingCollection" , 
			     "Name of the LCCD collection with channel mapping"  ,
			     _chMapCollection ,
			     std::string("ADCChannelMapping") ) ;
  
  registerOptionalParameter( "ChannelPositionTextFile" , 
			     "Optionally use a text file for the hardware channel to position mapping - overwrites mapping from LCCD and GEAR" ,
			     _chPosTextFile ,
			     std::string("channel_positions.txt") ) ;


  registerProcessorParameter( "DriftVelocity" , 
			      "drift velocity in mu/ns", 
			      _driftVelocity ,
			      float( 80. ) ) ;


  registerProcessorParameter( "ColorScaleMaxADC" , 
			      "ADC value used for the maximum of the color scale" ,
			      _colorScaleMaxADC ,
			      int( 32 ) ) ;


  registerProcessorParameter( "ColorScheme" , 
			      "color scheme - 1: hot , 2 : cold " , 
			      _colorScheme , 
			      int( 1 ) ) ;


}


void TrackerRawViewer::init() { 
  
  streamlog_out(DEBUG) << "   init called  "  << std::endl ;
  
  // usually a good idea to
  printParameters() ;
  
  if( ! parameterSet( "ChannelMappingCollection" ) && ! parameterSet( "ChannelPositionTextFile" ) ) {
    
    throw Exception(" No channel mapping provided - specify either \"ChannelMappingCollection\" or \"ChannelPositionTextFile\" " ) ;
  } 
  
  
  CEDPickingHandler::getInstance().registerFunction( LCIO::TRACKERRAWDATA , &printTRD  ) ;


  MarlinCED::init(this) ;
  
  
  if( parameterSet( "ChannelPositionTextFile" ) ) {
    
    streamlog_out( WARNING ) << "===================================================================== " << std::endl 
			     << "    reading channel position map from text file :     "  << std::endl 
			     << "     " <<  _chPosTextFile << " - ignoring LCCD and GEAR !!!!  "  << std::endl 
			     << "===================================================================== " << std::endl 
			     << std::endl ;
    
    _posMap = new ChannelPosMap( _chPosTextFile ) ; 
    
  }
  else {   /// if( parameterSet( "ChannelMappingCollection" ) ) {
    
    _chMap = new ChannelMap( &tpcconddata::ADCChannelMapping::hardwareToSoftwareKey ) ;
    
    // register it for automatic update
    ConditionsProcessor::registerChangeListener( _chMap , _chMapCollection ) ;
  }
  

  // define some colors schemes - one of which is used for displaying tracks  
  static unsigned bojum_colors[200]={
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



  // static int Red         =   0 ;
  // static int Orange      =   1 ;
  // static int Plum        =   2 ;
  // static int Violet      =   3 ;
  // static int Blue        =   4 ;
  // static int LightBlue   =   5 ;
  // static int Aquamarine  =   6 ;
  // static int Green       =   7 ;
  // static int Olive       =   8 ;
  // static int Yellow      =   9 ;

  // fixme: make this a parameter
  static int scheme = LightBlue  ;

  for(int i = 0 ; i<20 ; _colors.push_back( bojum_colors[  ( scheme * 20  )  + i++ ] )  ) ;


  
  _nRun = 0 ;
  _nEvt = 0 ;
  
}

void TrackerRawViewer::processRunHeader( LCRunHeader* run) { 
  
  _nRun++ ;
} 

void TrackerRawViewer::processEvent( LCEvent * evt ) { 
  
  
  MarlinCED::newEvent( this , -1 ) ; // -1 suppresses drawing of detector (done here)
  
  //--------  enable picking ------
  CEDPickingHandler &pHandler=CEDPickingHandler::getInstance();
  pHandler.update( evt ); 
  
  
  // ---- get geometry from GEAR ----
  _tpcParams = &Global::GEAR->getTPCParameters() ;
  
  float zAnode = Global::GEAR->getTPCParameters().getMaxDriftLength();
  
  float tpcOuterR = _tpcParams->getDoubleVal("TPC_outer_radius") ;
  double offx     = _tpcParams->getDoubleVal("TPC_cylinder_centre_c0") ;
  double offy     = _tpcParams->getDoubleVal("TPC_cylinder_centre_c1") ;

  // get readout frequency from first module (can this be different for different modules ? )
  double readOutFrequency = _tpcParams->getModule(0).getReadoutFrequency() ;
  
  double mmPerReadoutTick = _driftVelocity * 1.e6  / readOutFrequency ;

  streamlog_out( DEBUG ) << " readOutFrequency = " << readOutFrequency 
			 << " driftVelocity " << _driftVelocity
			 << " mmPerReadoutTick " << mmPerReadoutTick
			 << std::endl ;
  

  gear::Vector3D offset( offx ,offy , 0. ) ; 
  

  // --- draw a cylinder for the TPC  ----

  const static int nCorners =       40 ;
  const static int tpcColor = 0xaaaaaa ;
  const static int tpcPhi0  =        0. ;
  
  static CED_GeoCylinder geoCylindersANY[] = {                      
    { tpcOuterR , nCorners , tpcPhi0,  zAnode ,   -zAnode , tpcColor  }, 
  };
  
  ced_geocylinders( sizeof(geoCylindersANY)/sizeof(CED_GeoCylinder), geoCylindersANY );

  //========================================================================
 

  if( isFirstEvent() ) {   }
    
  int layer = 1 ;
  int size = 1 ;
  int marker = 0 ;

  int ml = marker | ( layer << CED_LAYER_SHIFT ) ;
  
  LCCollection* col = evt->getCollection( _colName ) ;
  
  MarlinCED::add_layer_description( _colName, layer); 

  int nData = col->getNumberOfElements()  ;
  
  for(int i=0; i< nData ; i++){
    
    TrackerRawData* data = dynamic_cast<TrackerRawData*>( col->getElementAt( i ) ) ;
    
    double padCenter[2] ;

    bool padFound  = getPadCenter ( padCenter , data->getCellID0() , data->getCellID1() ) ;

    if( !padFound ) {

      continue ; // nothing to draw 
    }

    gear::Vector3D pos( padCenter[0] - offset[0] , padCenter[1] - offset[1]  ,  0.  ) ;


    //======= draw the ADC spectrum ===========================================

    const ShortVec& adc = data->getADCValues() ;
    
    for(unsigned int j=0 ; j< adc.size() ; ++j ) {
      
      if( adc[j] > 0 ) { 
	
	
	pos[2] =  zAnode  - ( ( data->getTime() + j )  *  mmPerReadoutTick ) ;
	
	
	// select color based on maxADC and color scheme
	unsigned int rgb[3] ;
	
	ColorMap::selectColorMap( _colorScheme )( rgb, adc[j] , 0 ,  _colorScaleMaxADC  ) ;
	
	int color = ColorMap::RGB2HEX( rgb[0], rgb[1], rgb[2] ); 
	
	ced_hit_ID( pos.x()   ,
		    pos.y()  ,
		    pos.z() ,
		    ml, size , color, data->id() ) ;
      }
    }
  }    

  // ===============================================================================
  //       optionally draw hits
  // ===============================================================================

  layer = 2 ;
  size = 9 ;
  marker = 2 ;
  
  try{ 
    
    LCCollection* hitCol = evt->getCollection( _hitColName ) ;
    
    MarlinCED::add_layer_description( _hitColName, layer); 

    int ml = marker | ( layer << CED_LAYER_SHIFT ) ;
    
    int nData = hitCol->getNumberOfElements()  ;
    
    
    for(int i=0; i< nData ; i++){
      
      TrackerHit* hit = dynamic_cast<TrackerHit*>( hitCol->getElementAt( i ) ) ;
      
      int color =  0xe2df04 ; //  _colors[ i % _colors.size() ] ;
      //yellows:    0x565501,0x727002,0x898702,0xa5a303,0xb7b403,0xd1cd04,0xe2df04,
      
      
      gear::Vector3D pos( hit->getPosition()[0]  - offset[0], 
			  hit->getPosition()[1]  - offset[1], 
			  hit->getPosition()[2] ) ;
      
      ced_hit_ID( pos.x()   ,
		  pos.y()  ,
		  pos.z() ,
		  ml, size , color, hit->id() ) ;
      
    }
    
  } catch( lcio::DataNotAvailableException& e) {
    
    streamlog_out( DEBUG ) << " no hit collection in event .... " << std::endl ;
  }
  
    

  // ===============================================================================
  //       optionally draw tracks
  // ===============================================================================

  layer = 3 ;
  size = 4 ;
  marker = 0 ;
  
  try{ 
    
    LCCollection* trkCol = evt->getCollection( _trkColName ) ;
    
    MarlinCED::add_layer_description( _trkColName, layer); 
      

    int ml = marker | ( layer << CED_LAYER_SHIFT ) ;
    
    int nData = trkCol->getNumberOfElements()  ;
    
    for(int i=0; i< nData ; i++){
      
      Track* trk = dynamic_cast<Track*>( trkCol->getElementAt( i ) ) ;
      
 
      const TrackerHitVec& hits = trk->getTrackerHits() ;
      
      for(unsigned int j=0; j< hits.size() ; j++){
	
	TrackerHit* hit = hits[j] ;
	
	int color = _colors[ i % _colors.size() ] ;
	
	gear::Vector3D pos( hit->getPosition()[0]  - offset[0], 
			    hit->getPosition()[1]  - offset[1], 
			    hit->getPosition()[2] ) ;
	
	ced_hit_ID( pos.x()   ,
		    pos.y()  ,
		    pos.z() ,
		    ml, size , color, hit->id() ) ;
	
      }
    
      
      // draw the helix:
      double bField = Global::GEAR->getBField().at(  gear::Vector3D(0,0,0)  ).z() ; 
      
      double pt = bField * 3e-4 / std::abs( trk->getOmega() ) ;

      if( bField == 0.0 ){
	pt = 1000. ; // this will result in a straight line be drawn by MarlinCED::drawHelix
      }
      double charge = ( trk->getOmega() > 0. ?  1. : -1. ) ;
      
      double px = pt * std::cos(  trk->getPhi() ) ;
      double py = pt * std::sin(  trk->getPhi() ) ;
      double pz = pt * trk->getTanLambda() ;
      
#define Correct_Track_Params
#ifdef Correct_Track_Params
      
      double xs = trk->getReferencePoint()[0] -  trk->getD0() * sin( trk->getPhi() )  - offset[0] ;
      double ys = trk->getReferencePoint()[1] +  trk->getD0() * cos( trk->getPhi() )  - offset[1] ;
      double zs = trk->getReferencePoint()[2] +  trk->getZ0() ;
      
#else  // assume track params have reference point at origin ...
      
      double xs = 0 -  trk->getD0() * sin( trk->getPhi() ) ;
      double ys = 0 +  trk->getD0() * cos( trk->getPhi() ) ;
      double zs = 0 +  trk->getZ0() ;
#endif        
      
            
      ml = marker | ( layer << CED_LAYER_SHIFT ) ;
      
      //      if( pt > 0.01 ) // sanity check
      MarlinCED::drawHelix( bField , charge, xs, ys, zs , 
			    px, py, pz, ml , 1 ,  0xffffff ,
			    0.0,  tpcOuterR,  //_tpcParams->getPlaneExtent()[1] , 
			    zAnode+100., trk->id() ) ;
      
      gear::Vector3D vs( xs, ys, zs ) ;
      gear::Vector3D vp( px, py, pz ) ;

      streamlog_out( DEBUG ) << " *** drawing helix with start: " << vs << std::endl 
			     << "       momentum : " <<   vp<< std::endl 
			     << "    pt: " << pt<< std::endl 
			     << "    trk->getOmega() : " << trk->getOmega()<< std::endl 
			     << " bField " << bField
			     << std::endl ;
 
          
    }
    
  } catch( lcio::DataNotAvailableException& e) {
    
    streamlog_out( DEBUG ) << " no hit collection in event .... " << std::endl ;
  }


  //=======================================================================================



  streamlog_out(MESSAGE) << "   processing event: " << evt->getEventNumber() 
			 << "   in run:  "          << evt->getRunNumber() 
			 <<  std::endl ;
    
    
    
  int waitForKeyboard =  true ; //false ;
  
  //++++++++++++++++++++++++++++++++++++
  MarlinCED::draw( this , waitForKeyboard ) ;
  //++++++++++++++++++++++++++++++++++++
  
  _nEvt ++ ;
}


bool TrackerRawViewer::getPadCenter( double *padCenter, int cellID0, int cellID1 ) {

  bool padFound  = false ;

  // if text file was given use this 
  if( _posMap != 0 ){

    ChannelPosMap::iterator it = _posMap->find(   std::make_pair( cellID0 , cellID1 )  ) ;

    if( it != _posMap->end() ) {

      // text file has prototype coordinates - not global GEAR ones !
      double offx     = _tpcParams->getDoubleVal("TPC_cylinder_centre_c0") ;
      double offy     = _tpcParams->getDoubleVal("TPC_cylinder_centre_c1") ;
      
      //----- it seems that the text file has x and y exchanged wrt to GEAR !!!???
      padCenter[0] = it->second[1] + offx ;
      padCenter[1] = it->second[0] + offy ;

      padFound = true ;

    } else{

      streamlog_out( WARNING ) << "  no position for cellID0=" << cellID0 << "cellID1 " << cellID1 << std::endl ;
    } 

  } else {
    
    tpcconddata::ADCChannelMapping* cm = 0 ;
    
    try {
      
      cm = &_chMap->find(  std::make_pair(  cellID0, cellID1 )  ) ;   
      
      const gear::TPCModule& tpcMod = _tpcParams->getModule( cm->getModuleID() ) ;
      
      int padIndex = cm->getPadID() ;
      
      gear::Vector2D padC = tpcMod.getPadCenter( padIndex ) ;
      
      padCenter[0] = padC[0] ;
      padCenter[1] = padC[1] ;
      
      padFound  = true ;
    
    } catch(lccd::DataNotAvailableException& ){}

  }


  return padFound ;
}

void TrackerRawViewer::check( LCEvent * evt ) { 
  // nothing to check here - could be used to fill checkplots in reconstruction processor
}


void TrackerRawViewer::end(){ 
  
  //   std::cout << "TrackerRawViewer::end()  " << name() 
  // 	    << " processed " << _nEvt << " events in " << _nRun << " runs "
  // 	    << std::endl ;

}

