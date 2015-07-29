//by Thorben Quast, Summer Student 2015
//19/07/2015

#include "DrawDetectorDD4hep.h"
#include "MarlinCED.h"
#include <iostream>

#include "DD4hep/LCDD.h"
#include "DDRec/DetectorData.h"
#include "DD4hep/DD4hepUnits.h" 
//new:
#include "DDRec/SurfaceManager.h"
#include "DDRec/Surface.h"

#include "TColor.h"
#include <cmath>



// ----- include for verbosity dependend logging ---------
#include "marlin/VerbosityLevels.h"

using namespace lcio ;
using namespace marlin ;
using namespace DD4hep ;
using namespace DD4hep::Geometry ;
using namespace DD4hep::DDRec ;

DrawDetectorDD4hep aDrawDetectorDD4hep ;


DrawDetectorDD4hep::DrawDetectorDD4hep() : Processor("DrawDetectorDD4hep") {
    
  // modify processor description
  _description = "DrawDetectorDD4hep draws a detector in CED using a DD4hep model" ;
  
  
  registerProcessorParameter( "DrawDetector" ,
			      "Call the drawDetector function else call draw function",
			      _begin ,
			      (bool)0  ) ;
  //example StringVec to indicate the function how the list is meant to be read out.
  StringVec layerExample ; layerExample.push_back("NameToBeDrawnDetailly" ) ;
  registerOptionalParameter( "DetailledDrawing" ,
            "List of detector names that are printed in more detail.",
            _detailled,
            layerExample,
            1  ) ;  
  //  
  registerOptionalParameter( "DrawSurfaces" ,
            "Draw the geometry as a set of individual surfaces (if available) instead of simplified structures.",
            _surfaces ,
            (bool)0  ) ;
}



void DrawDetectorDD4hep::init() {
    
  streamlog_out(DEBUG) << "   init called  " << std::endl ;
  
  // usually a good idea to
  printParameters() ;
  
  MarlinCED::init(this);
    
  _nRun = 0 ;
  _nEvt = 0 ;
  
}


void DrawDetectorDD4hep::processRunHeader( LCRunHeader* run) {
    
    _nRun++ ;
}



void DrawDetectorDD4hep::processEvent( LCEvent * evt ) {
    
    MarlinCED::newEvent(this , -1 );

    drawDD4hepDetector(DD4hep::Geometry::LCDD::getInstance()  );

    MarlinCED::draw(this, 1 );

    _nEvt ++ ;
}



void DrawDetectorDD4hep::check( LCEvent * evt ) {
    // nothing to check here - could be used to fill checkplots in reconstruction processor
} 


void DrawDetectorDD4hep::end(){
    
  streamlog_out( MESSAGE)  << "DrawDetectorDD4hep::end()  " << name()
			   << " processed " << _nEvt << " events in " << _nRun << " runs "
			   << std::endl ;
  
}
/****************/
/*Helpers*/
/****************/

void getVisAttributes(DD4hep::Geometry::DetElement det, unsigned &color, bool &visible) {
    DD4hep::Geometry::VisAttr thisVisAttribute = det.volume().visAttributes();  
    if (thisVisAttribute.isValid()){
      TColor* c = new TColor(thisVisAttribute.color(), 1, 1, 1);
      //convert the given TColor into a hexadecimal whereby R_i, G_i, B_i integers
      //color = 0x00|R_1 R_2|G_1 G_2|B_1 B_2|
      color = ((int(255*c->GetRed())<<16) |  //multiply with 2^16 to get the last two bits
        (int(255*c->GetGreen())<<8)|          //multiply with 2^8 to get the middle bits
        (int(255*c->GetBlue())<<0));          //multiply with 2^0 to get the first two bits
                                              //the | operator is a bitwise addition of the number
      delete c;
    }
    else{
      streamlog_out( MESSAGE )<<"color: pointer does not exist"<<std::endl;
        color  = 8947848; //== 0xff999999
        visible = true;   
    }
    //TODO: Hard coded to make every element visible for now
    visible = true;
}

//converts the parameters in LayeredCalorimeterData given by the appropriate drivers
//into those required by the CEDGeoTube
CEDGeoTubeParams CalorimeterParameterConversion (LayeredCalorimeterData *calo){
  //get all the information from the lcdd class
  double rMin = calo->extent[0]/dd4hep::mm ;
  double rMax = calo->extent[1]/dd4hep::mm ;
  //attention! the meaning of these two variables depends on their context 
  //(see e.g. ECalEndcap_o2_v01_geo.cpp lines 78/79 vs. ECalBarrel_o1_v01.cpp lines 74/75)
  double zMin = calo->extent[2]/dd4hep::mm ; 
  double zMax = calo->extent[3]/dd4hep::mm ;
  
  double inner_symmetry = calo->inner_symmetry;
  double outer_symmetry = calo->outer_symmetry;
  //Note that the phi0's are given in degrees rather than rad in this implementation
  double inner_phi0 = calo->inner_phi0*180/M_PI;
  double outer_phi0 = calo->outer_phi0*180/M_PI;
  //new convention to prevent weird overlaps
  int NCircle = 36;
  //correct for the case of circle approximation, the given is conventional
  if (outer_symmetry < 1)  outer_symmetry = NCircle;
  if (inner_symmetry < 1)  inner_symmetry = NCircle;

  CEDGeoTubeParams returnParams;
  //given: distance middle point - center of edge section
  //required: distance middle point - edge intersection (corner) to prevent weird overlaps
  returnParams.Rmax = rMax/cos(M_PI/outer_symmetry); 
  returnParams.Rmin = rMin/cos(M_PI/inner_symmetry);
  
  //nothing to convert
  returnParams.inner_symmetry = inner_symmetry;
  returnParams.outer_symmetry = outer_symmetry;
  
  //by default, CED draws tubes with the corner facing downward
  //what we want: phi0 = 0 <=> straight edge parallel to x-y plane
  //therefore, the tube must be rotated by 360 - (90 + 360./(2*number of sides)) degrees since 
  //respecting the rotation symmetry in 360/n, the minimal phi0 is calculated to prevent interference with phi cuts in the CEDViewer
  //phi0 == 0 implies a normal vector of the first cell parallel to the +x-axis 
  //but CED starts drawing symmetrically from the +y-axis
  returnParams.phi0 = outer_phi0 + 270. - 180./outer_symmetry;
  returnParams.phi0 = returnParams.phi0 - (360./outer_symmetry)*(int (returnParams.phi0/(360./outer_symmetry))); 
  //in ILD and CLIC, both inner and outer shapes agree. For a generic solution, another parameter in LayeredCalorimeterData (e.g. phi_inner) should be introduced
  //i.e. delta_phi that is the rotational angle of the inner angle with respect to the outer layer
  returnParams.delta_phi = -returnParams.phi0 + inner_phi0 + 270. - 180./inner_symmetry;
  returnParams.delta_phi = returnParams.delta_phi - (360./inner_symmetry)*(int (returnParams.delta_phi/(360./inner_symmetry)));

  

  //endcaps and barrels take the same parameters in CED but the interpretation of the z-coordinates differs:
  returnParams.isBarrel = calo->layoutType == LayeredCalorimeterData::BarrelLayout;
  
  //type specific conversions
  if (returnParams.isBarrel){
    //barrels are drawn centered at z=0, i.e. they start at -zMax (z0) and their half length (delta_z) is zMax
    returnParams.delta_z = zMax;
    returnParams.z0 = -zMax;
  }
  else {
    //In the drivers, zMin is given to be the start of one disk. (Symmetry at the z-axis is assumed during placement.)
    //zMax is passed as zMin + thickness such that half the distance is correctly obtained by the line below.
    returnParams.delta_z = 0.5*(zMax - zMin);
    returnParams.z0 = zMin;
  }
  return returnParams;
}
//converts the parameters in ZDiskPetalsData given by the appropriate drivers
//into those required by the CEDGeoTube
CEDGeoTubeParams PetalParameterConversion (std::vector<DDRec::ZDiskPetalsData::LayerLayout>::iterator thisLayer){

  double phi0 = thisLayer->phi0*180/M_PI;
  double distanceSensitive = thisLayer->distanceSensitive/dd4hep::mm;
  //Supposedly, this is the actual width of the ring. The expansion along the z-axis is hard coded to 0.3.
  double lengthSensitive = thisLayer->lengthSensitive/dd4hep::mm;
  double zPosition = thisLayer->zPosition/dd4hep::mm;
  double thicknessSupport = thisLayer->thicknessSupport/dd4hep::mm;
  int petalNumber = thisLayer->petalNumber;

  CEDGeoTubeParams returnParams;
  /*   
  //Old implementation with GEAR
  returnParams.Rmax = distanceSensitive;
  returnParams.Rmin = distanceSensitive-1.0*lengthSensitive;
  */
  returnParams.Rmax = distanceSensitive+1.0*lengthSensitive;
  returnParams.Rmin = distanceSensitive;

  //(see comment in line 151)
  returnParams.Rmax = returnParams.Rmax/cos(M_PI/petalNumber);
  returnParams.Rmin = returnParams.Rmin/cos(M_PI/petalNumber);
  //Number of edges = number of petals
  returnParams.inner_symmetry = petalNumber;
  returnParams.outer_symmetry = petalNumber;
  //(see comment in line 159)
  returnParams.phi0 = phi0 + 270. - 180./petalNumber - (360./petalNumber)*(int ((270.-180./petalNumber)/(360./petalNumber)));
  returnParams.delta_phi = 0.0;
  
  //thicknessSupport is negligibly small
  returnParams.delta_z = thicknessSupport; 
  //Again: z0 is the left handed starting point for drawing
  returnParams.z0 = zPosition-returnParams.delta_z;
  return returnParams;
}

//converts the parameters from a LayeredCalorimeterData layer given by the appropriate drivers
//into those required by the CEDGeoTube
CEDGeoTubeParams CalorimeterLayerParameterConversion(std::vector<DDRec::LayeredCalorimeterData::Layer>::iterator thisLayer){
  double distance = thisLayer->distance/dd4hep::mm;
  double thickness = thisLayer->thickness/dd4hep::mm;
  double cellSize0 = thisLayer->cellSize0/dd4hep::mm;
  double cellSize1 = thisLayer->cellSize1/dd4hep::mm;
  
  int NCircle = 36;   //hard coded number of edges to form a circle

  CEDGeoTubeParams returnParams;
  returnParams.Rmax = distance + 1.0*thickness;
  returnParams.Rmin = distance;
  //(see comment in line 151)
  returnParams.Rmax = returnParams.Rmax/cos(M_PI/NCircle);
  returnParams.Rmin = returnParams.Rmin/cos(M_PI/NCircle);
  //assume round edges
  returnParams.inner_symmetry = NCircle;
  returnParams.outer_symmetry = NCircle;
  returnParams.phi0 = 270. - 180./NCircle - (360./NCircle)*(int ((270.-180./NCircle)/(360./NCircle)));
  returnParams.delta_phi = 0.0;

  //cellSize1 is half the length along the z-axis (see e.g. Solenoid_o1_v01_gep.cpp line 58/59)
  returnParams.delta_z = cellSize1;
  //cellSize0 is defined to be the middle point of the geometry 
  returnParams.z0 = -cellSize1 + cellSize0;

  return returnParams;
}

//converts the parameters from a FixedPadSizeTPCData given by the appropriate drivers
//into those required by the CEDGeoTube
CEDGeoTubeParams TPCParameterConversion(FixedPadSizeTPCData *tpc){
  double zHalf = tpc->zHalf/dd4hep::mm;
  //these radii include the insensitive space!
  double rMin = tpc->rMin/dd4hep::mm;
  double rMax = tpc->rMax/dd4hep::mm;
  
  int NCircle = 36;   //hard coded number of edges to form a circle

  CEDGeoTubeParams returnParams;

  returnParams.Rmax = rMax;
  returnParams.Rmin = rMin;
  //(see comment in line 151)
  returnParams.Rmax = returnParams.Rmax/cos(M_PI/NCircle);
  returnParams.Rmin = returnParams.Rmin/cos(M_PI/NCircle);
  //assume round edges
  returnParams.inner_symmetry = NCircle;
  returnParams.outer_symmetry = NCircle;
  returnParams.phi0 = 270. - 180./NCircle - (360./NCircle)*(int ((270.-180./NCircle)/(360./NCircle)));
  returnParams.delta_phi = 0.0;

  
  returnParams.delta_z = zHalf;
  returnParams.z0 = -zHalf;

  return returnParams;
}

//converts the parameters from a ZPlanarData::LayerLayout layer given by the appropriate drivers
//into those required by the CEDGeoBox (for drawing of staves) or by CEDGeoTube (for approximation of the set of staves into tubes)
LayerGeometry TrackerLayerParameterConversion(std::vector<DDRec::ZPlanarData::LayerLayout>::iterator thisLayer){
  int nLadders = thisLayer->ladderNumber;
  double phi0 = thisLayer->phi0*180/M_PI;
  
  double distance_sensitive = thisLayer->distanceSensitive/dd4hep::mm ;
  double thickness_sensitive = thisLayer->thicknessSensitive/dd4hep::mm ;
  double width_sensitive = thisLayer->widthSensitive/dd4hep::mm ;
  double offset_sensitive = thisLayer->offsetSensitive/dd4hep::mm ;
  double zHalf_sensitive = thisLayer->zHalfSensitive/dd4hep::mm ;  
  
  LayerGeometry Geometry;
  //two possibilites to draw: 
  //1) Draw layer consisting of staves
  double currentPhi; double radius;  double deltaPhi = 360./nLadders;
  for (int i=0; i<nLadders; i++){
    //placement shall begin along -z axis like for all other geometries
    currentPhi = phi0 + i*deltaPhi;
    
    //distance_sensitive is passed to be the minimal radius by the driver;
    //for drawing, it should be converted to the radius of the middle line
    radius = distance_sensitive+0.5* thickness_sensitive;

    CEDGeoBox stave;
    //place the center of the box at the appropriate coordinates
    //offset_sensitive is an additonal shift in placement orthogonal to the original 2D radius vector
    stave.center[0] = (radius*cos(currentPhi*M_PI/180) - offset_sensitive*sin(currentPhi*M_PI/180));
    stave.center[1] = (radius*sin(currentPhi*M_PI/180) + offset_sensitive*cos(currentPhi*M_PI/180));
    stave.center[2] = 0.0;  //placed z=0
    //dimensions are straight forward in an xyz coordinate system
    stave.sizes[0]  = thickness_sensitive;
    stave.sizes[1]  = width_sensitive;
    stave.sizes[2]  = zHalf_sensitive * 2;
    //the individual staves are finally rotated in the coordinate system to their appropriate position
    //herby, the rotation is performed around the z-axis
    stave.rotate[0] = 0.0;
    stave.rotate[1] = 0.0;
    stave.rotate[2] = currentPhi;

    Geometry.staves.push_back(stave);
  }


  //2) Summarize the set of staves into a geotube
  //analogous conversions as for the CalorimeterLayerParameterConversion (see above)
  Geometry.tube.Rmax = (distance_sensitive+thickness_sensitive)/cos(M_PI/nLadders); 
  Geometry.tube.Rmin = distance_sensitive/cos(M_PI/nLadders);
  Geometry.tube.inner_symmetry = nLadders;
  Geometry.tube.outer_symmetry = nLadders;
  Geometry.tube.phi0 = phi0 + 270. - 180./nLadders - (360./nLadders)*(int ((270.-180./nLadders)/(360./nLadders)));
  Geometry.tube.delta_phi = 0.0;
  Geometry.tube.delta_z = zHalf_sensitive;
  Geometry.tube.z0 = - zHalf_sensitive;

  //return both possible sets of parameters describing the geometry. The choice for either is implemented in the main draw routine
  return Geometry;
}

bool DrawSurfaces(DD4hep::DDRec::SurfaceManager &surfMan, std::string detName, unsigned color, int layer){
  typedef DD4hep::DDRec::SurfaceMap SMap;
  const SMap* sMap = surfMan.map(detName);
  int lineCounter = 0;
  if(sMap) {
    for (SMap::const_iterator it = sMap->begin(); it != sMap->end(); ++it){
      DD4hep::DDRec::Surface* surf = dynamic_cast<DD4hep::DDRec::Surface*> (it->second);
      if (!surf) continue;
      if (!(surf->type().isVisible())) continue;
      const std::vector<std::pair<Vector3D,Vector3D> > lines = surf->getLines();
      if (lines.empty()){
        streamlog_out( MESSAGE )<<" **** drawSurfaces(): empty lines vector for surface "<< *surf <<std::endl;
        continue;
      }
      for(unsigned i = 0; i <lines.size(); i++){
        unsigned default_width = 2;
        unsigned default_type = layer;
        ced_line(lines[i].first.x()/dd4hep::mm ,lines[i].first.y()/dd4hep::mm ,lines[i].first.z()/dd4hep::mm ,
                lines[i].second.x()/dd4hep::mm ,lines[i].second.y()/dd4hep::mm ,lines[i].second.z()/dd4hep::mm ,
                default_type,default_width, color);
        ced_line_ID(lines[i].first.x()/dd4hep::mm ,lines[i].first.y()/dd4hep::mm ,lines[i].first.z()/dd4hep::mm ,
                lines[i].second.x()/dd4hep::mm ,lines[i].second.y()/dd4hep::mm ,lines[i].second.z()/dd4hep::mm ,
                default_type,default_width, color, layer);
        lineCounter++; 
      }
    }
  } 
  return lineCounter > 0; //at least one line must have been drawn, otherwise the surface is considered to be undrawn
}

bool DrawDetectorDD4hep::detailledDrawing(std::string detName){
  if(! parameterSet("DetailledDrawing")) 
    return false;
  unsigned index = 0;
  while(index < _detailled.size()){
    if (detName.compare(_detailled[index])==0)
      return true;
    index++;
  }
  return false;
}
//End of helpers

void DrawDetectorDD4hep::drawDD4hepDetector( DD4hep::Geometry::LCDD& lcdd ){
  typedef std::vector< DD4hep::Geometry::DetElement> DetVec ;
  // get DetElements for the main sub detectors from dd4hep 
  const DetVec& trackers     = lcdd.detectors( "tracker" ) ;
  const DetVec& calorimeters = lcdd.detectors( "calorimeter" ) ;
  const DetVec& passiveDets  = lcdd.detectors( "passive" ) ;
  //allocate reference to the surface manager
  DD4hep::DDRec::SurfaceManager& surfMan = *lcdd.extension<DD4hep::DDRec::SurfaceManager>();
  //some temporary parameters for visualization
  unsigned color; bool visible;
  //temporary objects
  DetElement det; std::string detName;

  std::vector<CEDGeoTube> gTV ; 
  int detLayer = NUMBER_DATA_LAYER; 


  //--- loop over all calorimeters
  for( unsigned i=0,n=calorimeters.size() ; i<n ; ++i ){
    det = calorimeters[i] ;
    detName = det.name() ;
    LayeredCalorimeterData* calo = 0 ;
    streamlog_out( MESSAGE ) << " ......processing " << detName << std::endl;  
    //try to get the appropriate extension
    try{ 
      calo = det.extension<LayeredCalorimeterData>() ; 
    } catch(std::runtime_error& e){
      streamlog_out( MESSAGE ) <<  detName 
             << " has no extension of type LayeredCalorimeterData - nothing to draw "
             <<   std::endl;  
    }

    //get the visAttributes of the detElement's volume or (if not existing) some default values
    getVisAttributes(det, color, visible);

    //draw if the object exists and its visibility is set true
    if( calo != 0 && visible) {                    
      streamlog_out( MESSAGE ) <<" DRAWING " << detName << std::endl;
      int layer = detLayer++;
      if (this->_surfaces && DrawSurfaces(surfMan, detName, color, layer)){
        streamlog_out( MESSAGE )<<"Surfaces are drawn"<<std::endl;
      }else{
        //get the required parameter set for drawing a CEDGeoTube
        CEDGeoTubeParams params;    
        params = CalorimeterParameterConversion(calo);
        //consistently allocated this helper for linking the element correctly to the GUI control via ID
        gTV.push_back( CEDGeoTube( params.Rmax, params.Rmin, params.outer_symmetry, params.inner_symmetry, params.phi0, params.delta_phi, params.delta_z,  params.z0, color , layer  ,1,1 ) ) ; 
        //gTV.push_back( CEDGeoTube( params.Rmax, params.Rmin, params.outer_symmetry, params.inner_symmetry, 11.5, 0.0, params.delta_z,  params.z0, color , layer  ,1,1 ) ) ; 
        MarlinCED::set_layer_description( detName , layer );
        if (!params.isBarrel){
          //place the second one symmetric to the z-axis. An additional shift by the width of the layer is needed since the appropriate argument represents the left handed start of the geometry
          gTV.push_back( CEDGeoTube( params.Rmax, params.Rmin, params.outer_symmetry, params.inner_symmetry, params.phi0, params.delta_phi, params.delta_z,  - (params.z0+2.0*params.delta_z), color , layer  ,1,1 ) ) ; 
        }
      }
    }
  }
  //--- draw trackers 
  //mostly repetition of the calorimeter 
  for( unsigned i=0,n=trackers.size() ; i<n ; ++i ){
    det = trackers[i] ;
    detName = trackers[i].name() ;
    ZPlanarData* trkPlanar = 0; ZDiskPetalsData* trkDisk = 0; FixedPadSizeTPCData* trkTPC = 0;
    streamlog_out( MESSAGE ) << " ......processing" <<  detName << std::endl; 
    
    bool drawflag = true;
    try{ 
      trkPlanar = det.extension<ZPlanarData>();
    } catch(std::runtime_error& e){
      try{
        trkDisk = det.extension<ZDiskPetalsData>();
      }catch(std::runtime_error& e){
        try{
          trkTPC = det.extension<FixedPadSizeTPCData>();
        } catch(std::runtime_error& e){
            drawflag = false; 
            streamlog_out( MESSAGE ) <<  detName 
             << " has no extension of type ZPlanarData/ZDiskPetalsData - nothing to draw "
             <<   std::endl;           
        }
      }
    }

    getVisAttributes(det, color, visible);
    //get the visAttributes of the detElement's volume or (if not existing) some default values
    if (drawflag && visible){
      streamlog_out( MESSAGE ) <<" DRAWING " << detName << std::endl;
      int layer = detLayer++;
      
      if (this->_surfaces && DrawSurfaces(surfMan, detName, color, layer)){
        streamlog_out( MESSAGE )<<"Surfaces are drawn"<<std::endl;
      }else{
        //the following if statements are exclusive, i.e. only one may apply
        if(trkPlanar){
          for (std::vector<DDRec::ZPlanarData::LayerLayout>::iterator thisLayer = trkPlanar->layers.begin(); thisLayer != trkPlanar->layers.end(); thisLayer++){
            LayerGeometry Geo;
            Geo = TrackerLayerParameterConversion(thisLayer);
            if (detailledDrawing(detName)){
              for( unsigned stave_i=0; stave_i<Geo.staves.size() ; ++stave_i ){
                ced_geobox_r_ID( Geo.staves[stave_i].sizes, Geo.staves[stave_i].center, Geo.staves[stave_i].rotate, color, layer,0);
                ced_geobox_r_solid( Geo.staves[stave_i].sizes, Geo.staves[stave_i].center, Geo.staves[stave_i].rotate, color, layer);
              }
            }
            else{ 
              gTV.push_back( CEDGeoTube( Geo.tube.Rmax, Geo.tube.Rmin, Geo.tube.inner_symmetry, Geo.tube.outer_symmetry, Geo.tube.phi0, Geo.tube.delta_phi, Geo.tube.delta_z,  Geo.tube.z0, color , layer  ,1,1 ) );
            }
          }
        }

        if(trkDisk){ 
          if (detailledDrawing(detName)){
            streamlog_out( MESSAGE )<<detName<<": Not drawn for now (appropriate geometry does not exist)"<<std::endl;
          }
          else{
            for (std::vector<DDRec::ZDiskPetalsData::LayerLayout>::iterator thisLayer = trkDisk->layers.begin(); thisLayer != trkDisk->layers.end(); thisLayer++){
              CEDGeoTubeParams params;    
              params = PetalParameterConversion(thisLayer);
              gTV.push_back( CEDGeoTube( params.Rmax, params.Rmin, params.outer_symmetry, params.inner_symmetry, params.phi0, params.delta_phi, params.delta_z,  params.z0, color , layer  ,1,1 ) ) ; 
              //place the second one symmetric to the z-axis. An additional shift by the width of the layer is needed since the appropriate argument represents the left handed start of the geometry
              gTV.push_back( CEDGeoTube( params.Rmax, params.Rmin, params.outer_symmetry, params.inner_symmetry, params.phi0, params.delta_phi, params.delta_z,  -(params.z0+2*params.delta_z), color , layer  ,1,1 ) ) ; 
              MarlinCED::set_layer_description( detName , layer );
            }
          }
        }

        if(trkTPC) {
          CEDGeoTubeParams params;    
          params = TPCParameterConversion(trkTPC);
          gTV.push_back( CEDGeoTube( params.Rmax, params.Rmin, params.outer_symmetry, params.inner_symmetry, params.phi0, params.delta_phi, params.delta_z,  params.z0, color , layer  ,1,1 ) ) ; 
        }
      }
      MarlinCED::set_layer_description( detName , layer );
    }
  }


  //--- draw passive 
  for( unsigned i=0,n=passiveDets.size() ; i<n ; ++i ){
    det = passiveDets[i] ;
    detName = passiveDets[i].name() ;
    ConicalSupportData* passiveConical = 0; LayeredCalorimeterData* passiveCalo = 0;
    streamlog_out( MESSAGE ) << " ......processing " <<  detName << std::endl;  
    bool drawflag = true;
    try{ 
        passiveConical = det.extension<ConicalSupportData>();
    } catch(std::runtime_error& e){
      try{
        passiveCalo = det.extension<LayeredCalorimeterData>();
      } catch(std::runtime_error& e){
        drawflag = false;
          streamlog_out( MESSAGE ) <<  detName 
             << " has no extension of type ConicalSupportData/LayeredCalorimeterData - nothing to draw "
             <<   std::endl;  
      }
    }
    //get the visAttributes of the detElement's volume or (if not existing) some default values
    getVisAttributes(det, color, visible);
    if (drawflag && visible){
      int layer = detLayer++;
      streamlog_out( MESSAGE ) <<" DRAWING " << detName << std::endl;
      if (this->_surfaces && DrawSurfaces(surfMan, detName, color, layer)){
        streamlog_out( MESSAGE )<<"Surfaces are drawn"<<std::endl;
      }else{
        if (passiveConical){
          streamlog_out( MESSAGE )<<detName<<" is not drawn for now (not needed)."<<std::endl;
          for (std::vector<DDRec::ConicalSupportData::Section>::iterator thisSection = passiveConical->sections.begin(); thisSection != passiveConical->sections.end(); thisSection++){
          }
        }
        if (passiveCalo){
          for (std::vector<DDRec::LayeredCalorimeterData::Layer>::iterator thisLayer = passiveCalo->layers.begin(); thisLayer != passiveCalo->layers.end(); thisLayer++){
            CEDGeoTubeParams params;
            params = CalorimeterLayerParameterConversion(thisLayer);
            gTV.push_back( CEDGeoTube( params.Rmax, params.Rmin, params.outer_symmetry, params.inner_symmetry,  params.phi0, params.delta_phi, params.delta_z,  params.z0, color , layer  ,1,1 ) ) ; 
          }
        }
      }
      MarlinCED::set_layer_description( detName , layer );
    }
  }
  // ========================================================================
  //Draw the tubes:
  ced_geotubes( gTV.size() ,  (CED_GeoTube*) &gTV[0] );

  // ========================================================================
  

   
  MarlinCED::write_layer_description();
}


