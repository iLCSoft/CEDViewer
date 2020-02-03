/** DSTViewer - event display for DST files
 *  @author: S.Daraszewicz
 *  @date: 15.08.08
 **/
#include "DSTViewer.h"
#include <EVENT/ReconstructedParticle.h>
#include <iostream>
#include "ClusterShapes.h"

#include "HelixClass.h"
#include <math.h>
#include "MarlinCED.h"

#include <gear/GEAR.h>
#include <gear/BField.h>
#include <gear/TPCParameters.h>
#include <gear/PadRowLayout2D.h>
#include <gearimpl/Vector3D.h>
#include "ColorMap.h"

#include "ced_cli.h"

using namespace lcio ;
using namespace marlin ;

/** Define the key layers in the simulation */
#define PION_LAYER		1
#define PHOTON_LAYER	        2
#define NEUTRON_LAYER	        3		
#define NHADR_LAYER  	        4
#define CHADR_LAYER  	        5
#define TPC_LAYER    	        6
#define ECAL_LAYER   	        7
#define HCAL_LAYER   	        8
#define CLUSTER_LAYER	        9
#define HIT_LAYER		0
/** Jet layers... */
#define JET_DEFAULT_LAYER      11 
#define JET2_LAYER		1
#define JET3_LAYER		1
#define JET4_LAYER		1
#define JET5_LAYER		1
#define JET6_LAYER		1
/** momentum at ip layer */
#define MOM_LAYER		1
/** Two alternative cluster representations */
#define BACKUP_LAYER	        18
#define BACKUP_LAYER2	        19

/** Two alternative cluster representations */
#define IP_JET2			20
#define IP_JET3			21
#define IP_JET4			22
#define IP_JET5			23
#define IP_JET6			24

DSTViewer aDSTViewer ;

void DSTViewer::writeLayerDescription(void){
  /*
    int i;
    for(i=0;i<25;i++){
    //ced_describe_layer("", i);    //delete all old descriptions
    MarlinCED::set_layer_description("", i);
    }
  */
  MarlinCED::add_layer_description("Pions", PION_LAYER	);    
  MarlinCED::add_layer_description("Photons", PHOTON_LAYER	);    
  MarlinCED::add_layer_description("Neutrons", NEUTRON_LAYER);	   	
  MarlinCED::add_layer_description("NHadr", NHADR_LAYER  );	
  MarlinCED::add_layer_description("CHadr", CHADR_LAYER  );	
  MarlinCED::add_layer_description("TPC", TPC_LAYER    );	
  MarlinCED::add_layer_description("ECAL", ECAL_LAYER   );	
  MarlinCED::add_layer_description("HCAL", HCAL_LAYER   );	
  MarlinCED::add_layer_description("Clusters", CLUSTER_LAYER);	
  MarlinCED::add_layer_description("Hits", HIT_LAYER	);	
  MarlinCED::add_layer_description("JET2", JET2_LAYER	);	
  MarlinCED::add_layer_description("JET3", JET3_LAYER	);	
  MarlinCED::add_layer_description("JET4", JET4_LAYER	);	
  MarlinCED::add_layer_description("JET5", JET5_LAYER	);	
  MarlinCED::add_layer_description("JET6", JET6_LAYER	);	
  MarlinCED::add_layer_description("MOM", MOM_LAYER	);	
  MarlinCED::add_layer_description("Backup1", BACKUP_LAYER	);    
  MarlinCED::add_layer_description("Backup2", BACKUP_LAYER2);	
  MarlinCED::add_layer_description("IP_JET2", IP_JET2		);	
  MarlinCED::add_layer_description("IP_JET3", IP_JET3		);	
  MarlinCED::add_layer_description("IP_JET4", IP_JET4		);	
  MarlinCED::add_layer_description("IP_JET5", IP_JET5		);	
  MarlinCED::add_layer_description("IP_JET6", IP_JET6		);	
}


/**
 * DSTViewer description
 **/	
DSTViewer::DSTViewer() : Processor("DSTViewer") {

  _description = "CED based event display for DST files";  

  registerInputCollection( LCIO::RECONSTRUCTEDPARTICLE,
			   "ParticleCollection" , 
			   "Particle Collection Name" , 
			   _particleCollection , 
			   std::string("RecoParticles"));

  StringVec  jetColNames ;// typedef std::vector< std::string > StringVec ;  
  // jetColNames.push_back( "Durham_2Jets" ) ;
  // jetColNames.push_back( "Durham_3Jets" ) ;
  // jetColNames.push_back( "Durham_4Jets" ) ;
  // jetColNames.push_back( "Durham_5Jets" ) ;
  // jetColNames.push_back( "Durham_6Jets" ) ;

  registerInputCollections( LCIO::RECONSTRUCTEDPARTICLE,
			    "JetCollections" , 
			    "Jet Collections Names'" , 
			    _jetCollections , 
			    jetColNames );

  registerProcessorParameter("MagneticField",
			     "Magnetic Field",
			     _bField,
			     (float)3.5);

  registerProcessorParameter("WaitForKeyboard",
			     "Wait for Keyboard before proceed",
			     _waitForKeyboard,
			     (int)1);

  registerProcessorParameter("LayerReco",
			     "Layer for Reco Particles",
			     _layerReco,
			     (int)9);

  registerProcessorParameter( "DrawDetectorID" , 
			      "draw detector from GEAR file with given ID (see MarlinCED::newEvent() ) : 0 ILD, -1 none",
			      _detModel ,
			      0 ) ;
}

/**
 * DSTViewer initialise */
void DSTViewer::init() {

  _nRun = -1;
  _nEvt = 0;
  
  MarlinCED::init(this);

  // usually a good idea to
  printParameters() ;
}


/**
 * Process run header (?) */
void DSTViewer::processRunHeader( LCRunHeader* /*run*/) { 
  _nRun++ ;
  _nEvt = 0;
} 


/**
 * Main function processEvent */
void DSTViewer::processEvent( LCEvent * evt ) { 
  CEDPickingHandler &pHandler=CEDPickingHandler::getInstance();
  pHandler.update(evt); 
  //DSTViewer::writeLayerDescription();

  // Gets the GEAR global geometry params
  const gear::TPCParameters& gearTPC = Global::GEAR->getTPCParameters() ;
  const gear::PadRowLayout2D& padLayout = gearTPC.getPadLayout() ;

  // this defines the number of layers (?)
  int marker = 9;
	
  streamlog_out(MESSAGE)  << std::endl
			  << " (+++ This is the DSTViewer +++) " << std::endl
			  << "  DSTViewer  : event number "      << _nEvt << std::endl
			  << std::endl;
	
  _mcpList.clear();


  /**
   * Drawing Geometry
   * Reset drawing buffer and START drawing collection
   */ 
  MarlinCED::newEvent( this, _detModel ) ;

  /** Drawing Reconstructed Particles */
  try {

    LCCollection * col = evt->getCollection(_particleCollection.c_str());
    int nelem = col->getNumberOfElements();

    /** Total E, p, q values for the event */
    float TotEn = 0.0;
    float TotPX = 0.0;
    float TotPY = 0.0;
    float TotPZ = 0.0;
    float TotCharge = 0.0;
    float bField = Global::GEAR->getBField().at(  gear::Vector3D(0,0,0)  ).z() ;

    float ene_min = 0.0;
    float ene_max = 100.0;

    for (int ip(0); ip < nelem; ++ip) {
      ReconstructedParticle * part = dynamic_cast<ReconstructedParticle*>(col->getElementAt(ip));

      //TrackVec trackVec = part->getTracks();
      //int nTracks =  (int)trackVec.size();
      // this is needed for Calorimeter Hits display
      ClusterVec clusterVec = part->getClusters();
      int nClusters = (int)clusterVec.size();

      /** get particle properties */
      float ene = (float) part->getEnergy();
      float px  = (float)(part->getMomentum()[0]);
      float py  = (float)(part->getMomentum()[1]);
      float pz  = (float)(part->getMomentum()[2]);
      const double * pm = part->getMomentum();
      gear::Vector3D pp( pm[0], pm[1] , pm[2] ) ;
      float ptot = pp.r();
      float charge = (float)part->getCharge();
      int type = (int)part->getType();

      //Vertex *startVert = part->getStartVertex();
      //float ver_x = (float) startVert->getPosition()[0];
      //std::cout << "start vertex x = " << startVertx << std::endl;

      // start point
      float refx = 0.0;
      float refy = 0.0;
      float refz = 0.0;

      TotEn += ene;
      TotPX += px;
      TotPY += py;
      TotPZ += pz;
      TotCharge += charge;

      streamlog_out( DEBUG3 )	<< "Particle : " << ip
				<< " Type : " << type 
				<< " PX = " << px
				<< " PY = " << py
				<< " PZ = " << pz
				<< " P_tot = " << ptot
				<< " Charge = " << charge
				<< " E  = " << ene << std::endl;
				
      int ml = marker | ( TPC_LAYER<<CED_LAYER_SHIFT) ; // this defines the default layer 1
      if (type == 22){
	ml = marker | (ECAL_LAYER<<CED_LAYER_SHIFT); // photons on layer 2
      }
      if (type ==2112){
	ml = marker | (NEUTRON_LAYER<<CED_LAYER_SHIFT); //neutrons on layer 3
      }

      int color = returnTrackColor(type);
      int size = returnTrackSize(type);

      /** Draw the helix and straight lines */
      //MarlinCED::drawHelix(bField, charge, refx, refy, refz, px, py, pz, ml, size, color, 0.0, padLayout.getPlaneExtent()[1], 
      //				gearTPC.getMaxDriftLength());
							
      MarlinCED::drawHelix(bField, charge, refx, refy, refz, px, py, pz, ml, size, color, 0.0, padLayout.getPlaneExtent()[1], 
			   gearTPC.getMaxDriftLength(), part->id() ); //hauke: add id

      //				/** Draw momentum lines from the ip */
      char Mscale = 'b'; // 'b': linear, 'a': log
      int McolorMap = 2; //hot: 3
      int McolorSteps = 256;
      //float p = sqrt(px*px+py*py+pz*pz);
      float ptot_min = 0.0;
      float ptot_max = 25.0;
      int Mcolor = returnRGBClusterColor(ptot, ptot_min, ptot_max, McolorSteps, Mscale, McolorMap);
      int LineSize = 1;
      float momScale = 100;
      //ced_line(refx, refy, refz, momScale*px, momScale*py, momScale*pz, MOM_LAYER, LineSize, Mcolor);

      //hauke
      ced_line_ID(refx, refy, refz, momScale*px, momScale*py, momScale*pz, MOM_LAYER, LineSize, Mcolor, part->id()); //the right id?




      if (nClusters > 0 ) {
	// std::cout 	<< "nCluster > 0" << std::endl;
	// std::cout 	<<  clusterVec.size() <<std::endl;

	Cluster * cluster = clusterVec[0];

	double * center  = new double[3];
	float * center_r  = new float[3];
	center[0] = cluster->getPosition()[0];
	center[1] = cluster->getPosition()[1];
	center[2] = cluster->getPosition()[2];

	center_r[0] = (float)cluster->getPosition()[0];
	center_r[1] = (float)cluster->getPosition()[1];
	center_r[2] = (float)cluster->getPosition()[2];

	float phi = cluster->getIPhi();
	float theta = cluster->getITheta();

	if( phi ==0. && theta==0.){
	  // use the cluster postion:
	  theta = atan( sqrt( center[0]*center[0] + center[1]*center[1] ) / center[2]  ) ;
	  phi = atan2( center[1] , center[0] ) ;
	}

	float eneCluster = cluster->getEnergy();

	double rotate[] = {0.0, 0.0, 0.0};
	// for particles centered at the origin
	if (refx==0 && refy==0 && refz==0){
	  rotate[0] = 0.0;
	  rotate[1] = theta*180/M_PI;
	  rotate[2] = phi*180/M_PI;
	}
	else {
	  streamlog_out( DEBUG3 )	<< "This is not yet implemented" << std::endl;
	  return;
	} 


	double sizes[] = {100.0, 100.0, 100.0};
	int scale_z = 4;
	sizes[0] = returnClusterSize(eneCluster, ene_min, ene_max);
	sizes[1] = sizes[0];
	sizes[2] = scale_z*returnClusterSize(eneCluster, ene_min, ene_max);

	char scale = 'a'; // 'b': linear, 'a': log
	int colorMap = 3; //jet: 3
	int colorSteps = 256;
	int color2 = returnRGBClusterColor(eneCluster, ene_min, ene_max, colorSteps, scale, colorMap);

	//	int hit_type = 1 | HIT_LAYER;
					
	int cylinder_sides = 30;
	//problem with fisheye view
	ced_geocylinder_r(sizes[0]/2, sizes[2], center, rotate, cylinder_sides, color2, CLUSTER_LAYER); 



	//ced_hit(center[0],center[1],center[2], hit_type, (int)(sqrt(2)*sizes[0]/4), color2);
	ced_hit_ID(center[0],center[1],center[2], 1, HIT_LAYER, (int)(sqrt(2)*sizes[0]/4), color2, cluster->id()); //hauke

					
	int transparency = 0x66;
	int rgba = addAlphaChannelToColor(color2, transparency);
					
	//ced_cluellipse_r((float)sizes[0], (float)sizes[2], center_r, rotate, BACKUP_LAYER, rgba);
	ced_cluellipse_r_ID((float)sizes[0], (float)sizes[2], center_r, rotate, BACKUP_LAYER, rgba, cluster->id()); //hauke

					
	transparency = 0xCC;
	rgba = addAlphaChannelToColor(color2, transparency);
					
	//ced_ellipsoid_r(sizes, center, rotate, BACKUP_LAYER2, rgba); 
	ced_ellipsoid_r_ID(sizes, center, rotate, BACKUP_LAYER2, rgba, cluster->id()); //hauke

	/*
	 * End points for the arrow
	 */

	int sizeLine = 2;
	//int ml_line = 9;
	float radius = (float)0.5*sizes[2];
	float arrow = 0.2*radius;
	int color_arrow = color2 + 0x009900; //colour of the needle

	float xEnd = center[0] + (radius - arrow) * std::sin(theta)*std::cos(phi);
	float yEnd = center[1] + (radius - arrow) * std::sin(theta)*std::sin(phi);
	float zEnd = center[2] + (radius - arrow) * std::cos(theta);

	float xArrowEnd = center[0] + (radius) * std::sin(theta)*std::cos(phi);
	float yArrowEnd = center[1] + (radius) * std::sin(theta)*std::sin(phi);
	float zArrowEnd = center[2] + (radius) * std::cos(theta);

	// this is the direction arrow
	//ced_line(center[0],center[1],center[2], xEnd, yEnd, zEnd, CLUSTER_LAYER, sizeLine, color2);					
	ced_line_ID(center[0],center[1],center[2], xEnd, yEnd, zEnd, CLUSTER_LAYER, sizeLine, color2,cluster->id()); //hauke

	//ced_line(xEnd, yEnd, zEnd, xArrowEnd, yArrowEnd, zArrowEnd, CLUSTER_LAYER, sizeLine, color_arrow);
	ced_line_ID(xEnd, yEnd, zEnd, xArrowEnd, yArrowEnd, zArrowEnd, CLUSTER_LAYER, sizeLine, color_arrow, cluster->id());//hauke

      }
    }

    int color = 0x00000000;
    // ---- now we draw the jets (if any ) --------------------
    for(unsigned int i=0;i<_jetCollections.size();++i){ //number of jets

      streamlog_out( DEBUG )  << " drawing jets from collection " << _jetCollections[i] << std::endl ;

      LCCollection * col2 = evt->getCollection( _jetCollections[i] );
						
      int nelem2 = col2->getNumberOfElements();
						
      //int color = 0x000000;
						
						
      for (int j=0; j < nelem2; ++j) { //number of elements in a jet
	ReconstructedParticle * jet = dynamic_cast<ReconstructedParticle*>( col2->getElementAt(j) );
	streamlog_out( DEBUG )  <<   "     - jet energy " << jet->getEnergy() << std::endl ;
						
	const double * mom = jet->getMomentum() ;
	//const double jet_ene = jet->getEnergy();
	gear::Vector3D v( mom[0], mom[1] , mom[2] ) ;
							
	int layer = 0;
	const ReconstructedParticleVec & pv = jet->getParticles();
	float pt_norm = 0.0;
	for (unsigned int k = 0; k<pv.size(); ++k){
	  const double * pm = pv[k]->getMomentum();
	  gear::Vector3D pp( pm[0], pm[1] , pm[2] ) ;
	  gear::Vector3D ju = v.unit();
	  gear::Vector3D pt = pp - (ju.dot(pp))*ju;
	  pt_norm += pt.r();
								
	  layer = returnJetLayer(_jetCollections[i]);
	  color = returnJetColor(_jetCollections[i], j);
	  //								
	  //								MarlinCED::drawHelix(bField, pv[k]->getCharge(), 0.0, 0.0, 0.0, pm[0], pm[1], pm[2], layer, 1, color, 0.0, padLayout.getPlaneExtent()[1], 
	  //								gearTPC.getMaxDriftLength());
								
	  /** Draw momentum lines from the ip */
	  //char Mscale = 'b'; // 'b': linear, 'a': log
	  //int McolorMap = 2; //hot: 3
	  //int McolorSteps = 256;
	  //float ptot_min = 0.0;
	  //float ptot_max = 25.0;
	  //int Mcolor = returnRGBClusterColor(ptot, ptot_min, ptot_max, McolorSteps, Mscale, McolorMap);
	  //layer = returnJetLayer(_jetCollections[i]);
	  int LineSize = 1;
	  // start point
	  float refx = 0.0;
	  float refy = 0.0;
	  float refz = 0.0;
	  float momScale = 100;
	  int layerIp = returnIpLayer(_jetCollections[i]);
	  //ced_line(refx, refy, refz, momScale*pm[0], momScale*pm[1], momScale*pm[2], layerIp, LineSize, color);
	  ced_line_ID(refx, refy, refz, momScale*pm[0], momScale*pm[1], momScale*pm[2], layerIp, LineSize, color, pv[k]->id()); //hauke
	}
						   					   
	double center_c[3] = {0., 0., 0. };
	double rotation_c[3] = { 0.,  v.theta()*180./M_PI , v.phi()*180./M_PI };
	float RGBAcolor[4] = { float(0.2+0.2*i), float(0.2+0.2*i), float(1.0-0.25*i), float(0.3)};
	if (i==0){
	  RGBAcolor[0] = 1.0;
	  RGBAcolor[1] = 0.3;
	  RGBAcolor[2] = 0.1;
	  RGBAcolor[3] = 0.3;
	}
	double scale_pt = 20;
	double scale_mom = 25;
	double min_pt = 50;
	//ced_cone_r( min_pt + scale_pt*pt_norm , scale_mom*v.r() , center_c, rotation_c, layer, RGBAcolor);

	//std::cout<<"CONE: from 0,0,0 to: " << pm[0] << " " << pm[1] << " " << pm[2] << " ID: " << jet->id() <<  std::endl;
                            

	const double *pm=jet->getMomentum(); 
	int color2 = int(RGBAcolor[2]*(15*16+15)) + int(RGBAcolor[1]*(15*16+15))*16*16+ int(RGBAcolor[0]*(15*16+15))*16*16*16*16;

	//printf("RGBAcolor: %f %f %f color: %x\n", RGBAcolor[0], RGBAcolor[1], RGBAcolor[2], color2);
	//ced_line_ID(0.0, 0.0, 0.0, 33*pm[0], 33*pm[1], 33*pm[2], layer, 1, color2, jet->id()); //hauke
	int i2;
	for(i2=1;i2<104;i2++){
	  ced_line_ID((i2-1)/4*pm[0], (i2-1)/4*pm[1], (i2-1)/4*pm[2], i2/4*pm[0], i2/4*pm[1], i2/4*pm[2], layer, 1, color2, jet->id()); //hauke
	  //ced_hit_ID(i2/4*pm[0], i2/4*pm[1], i2/4*pm[2], layer, 2, color, jet->id()); //hauke

	}
	ced_hit_ID((i2-1)/4*pm[0], (i2-1)/4*pm[1], (i2-1)/4*pm[2], 0, layer, 2, color2, jet->id()); //hauke

                            


	ced_cone_r_ID( min_pt + scale_pt*pt_norm , scale_mom*v.r() , center_c, rotation_c, layer, RGBAcolor,jet->id()); //hauke
      }
    }



    char scale = 'a'; // 'b': linear, 'a': log
    int colorMap = 3; //jet: 3
    unsigned int colorSteps = 256;
    unsigned int ticks = 6; //middle 
    DSTViewer::showLegendSpectrum(colorSteps, scale, colorMap, ene_min, ene_max, ticks);

    streamlog_out( MESSAGE ) << std::endl
			     << "Total Energy and Momentum Balance of Event" << std::endl
			     << "Energy = " << TotEn
			     << " PX = " << TotPX
			     << " PY = " << TotPY
			     << " PZ = " << TotPZ
			     << " Charge = " << TotCharge << std::endl
			     << std::endl ;

    streamlog_out( DEBUG2 ) << "Setup properties" << std::endl
			    << "B_Field = " << bField << " T" << std::endl ;
				
  }
  catch( DataNotAvailableException &e){

    streamlog_out( DEBUG4 ) << "  data not available : " << e.what() << std::endl ;
  }

  DSTViewer::writeLayerDescription();
  /*
   * This refreshes the view? ...
   */
  MarlinCED::draw(this, _waitForKeyboard );
  _nEvt++;
}

void DSTViewer::check( LCEvent * /*evt*/ ) { }

void DSTViewer::end(){ } 

/**
 * PRELIM */
int DSTViewer::returnTrackColor(int type) {

  int kcol = 0x999999; //default black - unknown

  if (type==211) 			kcol = 0x990000; //red - pi^+
  else if (type==-211) 	kcol = 0x996600; //orange - pi^-
  else if (type==22) 		kcol = 0x00ff00; //yellow - photon
  else if (type==11) 		kcol = 0x660066; //dark blue - e^-
  else if (type==-11) 	kcol = 0x660099; //violet - e^+
  else if (type==2112)	kcol = 0x99FFFF; //white - n0
  else if (type==-2112)	kcol = 0x99FFFF; //white n bar
  else {
    streamlog_out( DEBUG ) << "Unassigned type of colour: default" << std::endl;
  }
  return kcol;
}


int DSTViewer::returnClusterColor(float eneCluster, float cutoff_min, float cutoff_max){

  /**
   * Check the input values for sanity */
  if (cutoff_min > cutoff_max) {
    std::cout << "Error in 'DSTViewer::returnClusterColor': cutoff_min < cutoff_max" << std::endl;
  }
  if (eneCluster < 0.0) {
    std::cout << "Error in 'DSTViewer::returnClusterColor': eneCluster is negative!" << std::endl;
  }
  if (cutoff_min < 0.0) {
    std::cout << "Error in 'DSTViewer::returnClusterColor': eneCluster is negative!" << std::endl;
  }

  int color = 0x000000; //default colour: black

  // Input values in log-scale
  float log_ene = std::log(eneCluster+1);
  float log_min = std::log(cutoff_min+1);
  float log_max = std::log(cutoff_max+1);

  int color_steps = 16*16; // 16*16 different steps from 0x600000 to 0x6000FF
  float log_delta = log_max - log_min;
  float log_step = log_delta/color_steps;

  int color_delta = (int) ((log_ene-log_min)/log_step); // which colour bin does the value go to? We have [0,16*16] bins

  //int color_delta_2 = (int)color_delta*16*16/2;
  //int color_delta_3 = (int)color_delta_2*16*16;

  if (color_delta >= 16*16){
    color = 0x990000;
  }
  else if (color_delta < 0){
    color = 0x9900FF;
  }
  else {
    color = 0x9900FF - color_delta;// + color_delta_2 - color_delta_3;
  }	

  //		// debug
  //		std::cout << "DEBUG: DSTViewer::returnClusterColor()" << std::endl;
  //		std::cout << "log_ene = " << log_ene << std::endl;
  //		std::cout << "log_min = " << log_min << std::endl;
  //		std::cout << "log_max = " << log_max << std::endl;
  //		std::cout << "log_step = " << log_step << std::endl;
  //		std::cout << "color_delta = " << color_delta << std::endl;
  //		using std::hex;
  //      	std::cout << hex << color << std::endl;

  return color;
}


void DSTViewer::showLegendSpectrum(const unsigned int color_steps, char scale, int colorMap, float ene_min, float ene_max, unsigned int ticks){

  const unsigned int numberOfColours = 3;
  unsigned int** rgb_matrix = new unsigned int*[color_steps]; //legend colour matrix;

  for (unsigned int i=0; i<color_steps; ++i){
    rgb_matrix[i] = new unsigned int[numberOfColours];
    ColorMap::selectColorMap(colorMap)(rgb_matrix[i], (float)i, 0.0, (float)color_steps);
    //			std::cout << "----------------------" << std::endl;
    //			std::cout << "i = " << i << std::endl;
    //			std::cout << "red = " << rgb_matrix[i][0] << std::endl;
    //			std::cout << "green = " << rgb_matrix[i][1] << std::endl;
    //			std::cout << "blue = " << rgb_matrix[i][2] << std::endl;
  }
  ced_legend(ene_min, ene_max, color_steps, rgb_matrix, ticks, scale); 
}

int DSTViewer::returnRGBClusterColor(float eneCluster, float cutoff_min, float cutoff_max, int color_steps, char scale, int colorMap){

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

  /**
   * Different color scales */
  switch(scale){
  case 'a': default: //log

    // log scale color_delta increment
    color_delta = (int) ((log_ene-log_min)/log_step); // which colour bin does the value go to? We have [0x00,0xFF] bins

    //std::cout << "Log color scale choosen" << std::endl;
    break;

  case 'b': //linear

    color_delta = (int)((eneCluster - cutoff_min)/(cutoff_max - cutoff_min)*color_steps);

    //std::cout << "linear color scale choosen" << std::endl;;
    break;
  }

  /*
   * Bins outside of the colour range */
  if (color_delta >= color_steps){
    color_delta = color_steps;
    //std::cout << "color_delta > color_steps" << std::endl;
  }
  if (color_delta < 0){
    color_delta = 0;
    //std::cout << "color_delta < 0" << std::endl;
  }

  ColorMap::selectColorMap(colorMap)(rgb, color_delta, 0, color_steps);
  //		// DEBUG
  //		std::cout << "color_delta = " << color_delta << std::endl;
  //		std::cout << "color_delta = " << color_delta << std::endl;
  //		std::cout << "red = " << rgb[0] << std::endl;
  //		std::cout << "green = " << rgb[1] << std::endl;
  //		std::cout << "blue = " << rgb[2] << std::endl;

  color = ColorMap::RGB2HEX(rgb[0],rgb[1],rgb[2]);

  return color;
}

/**
 * PRELIM */
int DSTViewer::returnTrackSize(float type){

  int size = 1; //default size

  if (type==2112)			size = 1; //positive hadron
  else if (type==-2112)	size = 1; //negative hadron

  return size;
}

int DSTViewer::returnClusterSize(float eneCluster, float cutoff_min, float cutoff_max){ 

  /**
   * Check the input values for sanity */
  if (cutoff_min > cutoff_max) {
    std::cout << "Error in 'DSTViewer::returnClusterSize': cutoff_min < cutoff_max" << std::endl;
  }
  if (eneCluster < 0.0) {
    std::cout << "Error in 'DSTViewer::returnClusterSize': eneCluster is negative!" << std::endl;
  }
  if (cutoff_min < 0.0) {
    std::cout << "Error in 'DSTViewer::returnClusterSize': eneCluster is negative!" << std::endl;
  }

  int size = 0; //default size: zero

		// sizes
  int size_min = 10;
  int size_max = 120;

  // Input values in log-scale
  float log_ene = std::log(eneCluster+1);
  float log_min = std::log(cutoff_min+1);
  float log_max = std::log(cutoff_max+1);

  int size_steps = size_max - size_min; // default: 90 step sizes
  float log_delta = log_max - log_min;
  float log_step = log_delta/size_steps;

  int size_delta = (int) ((log_ene-log_min)/log_step); // which size bin does the value go to?

  if (size_delta >= size_steps){
    size = size_max;
  }
  else if (size_delta < size_min){
    size = size_min;
  }
  else {
    size = size_min + size_delta;
  }	

  /**
   * Check the output */
  if (size <=0){
    std::cout << "Error in 'DSTViewer::returnClusterSize': return size is negative!" << std::endl;
  }

  //		std::cout << "DEBUG: DSTViewer::returnClusterSize()" << std::endl;
  //		std::cout << "log_ene = " << log_ene << std::endl;
  //		std::cout << "log_min = " << log_min << std::endl;
  //		std::cout << "log_max = " << log_max << std::endl;
  //		std::cout << "log_step = " << log_step << std::endl;
  //		std::cout << "size_delta = " << size_delta << std::endl;
  //		std::cout << size << std::endl;

  return size;
}
	
int DSTViewer::returnIpLayer(std::string jetColName){
		
  int layer = BACKUP_LAYER;

  if (jetColName == "Durham_2Jets"){
    layer = IP_JET2;
  }
  else if (jetColName == "Durham_3Jets"){
    layer = IP_JET3;
  }
  else if (jetColName == "Durham_4Jets"){
    layer = IP_JET4;
  }
  else if (jetColName == "Durham_5Jets"){
    layer = IP_JET5;
  }
  else if (jetColName == "Durham_6Jets"){
    layer = IP_JET6;
  }
  else {
    layer = BACKUP_LAYER;
  }
		
  return layer;
		
}
	
int DSTViewer::returnJetLayer(std::string jetColName){
		
  int layer = JET_DEFAULT_LAYER;

  if (jetColName == "Durham_2Jets"){
    layer = JET2_LAYER;
  }
  else if (jetColName == "Durham_3Jets"){
    layer = JET3_LAYER;
  }
  else if (jetColName == "Durham_4Jets"){
    layer = JET4_LAYER;
  }
  else if (jetColName == "Durham_5Jets"){
    layer = JET5_LAYER;
  }
  else if (jetColName == "Durham_6Jets"){
    layer = JET6_LAYER;
  }
	
  return layer;	
		
}
	
//	float * DSTViewer::returnConeColor(std::string jetColName){
//		float RGBAcolor[4] = {0.4+0.1*i, 0.2+0.2*i, 1.0-0.25*i, 0.3};
//		
//		float opacity = 0.3;
//
//		if (jetColName == "Durham_2Jets"){
//			RGBAcolor[4] = {1, 0.1, 0.0, opacity};
//		}
//		else if (jetColName == "Durham_3Jets"){
//			RGBAcolor[4] = {1, 0.1, 0.0, opacity};
//		}
//		else if (jetColName == "Durham_4Jets"){
//			RGBAcolor[4] = {1, 0.1, 0.0, opacity};	
//		}
//		else if (jetColName == "Durham_5Jets"){
//			RGBAcolor[4] = {1, 0.1, 0.0, opacity};
//		}
//		else if (jetColName == "Durham_6Jets"){
//			RGBAcolor[4] = {1, 0.1, 0.0, opacity};
//		}
//		else {
//			RGBAcolor[4] = {1, 0.1, 0.0, opacity};		
//		}
//		
//		return RGBAcolor;
//	}
	
int DSTViewer::addAlphaChannelToColor(int color, int alphaChannel){
  int rgba = 0xEE000000;
		
  rgba = (alphaChannel<<24) + color;
		
  //std::cout << "color = " << rgba << std::endl;
		
  return rgba;
}
	
	
int DSTViewer::returnJetColor(std::string /*jetColName*/, int colNumber){
  int color = 0x00000000;
		
  int white = 0x55555555;
  int black = 0x55000000;
  int red = 0x55550000;
  int green = 0x55008000;
  int yellow = 0x55555500;
  int fuchsia = 0x55550055;
  //int blue = 0x55000055;
  //int orange = 0x5555a500;
  //int violet = 0x55ee82ee;
  //int purple = 0x55800080;
  //int silver = 0x55c0c0c0;
  //int gold = 0x5555d700;
  //int gray = 0x55808080;
  //int aqua = 0x55005555;
  //int skyblue = 0x5587ceeb;
  //int lightblue = 0x55a558e6;
  //int khaki = 0x55f0e68c;

  // assign the color based on the index in the parameter JetCollections in the steering file

  switch (colNumber){
  case 0:
    color = yellow;
    break ;
  case 1:
    color = red;
    break ;
  case 2:
    color = white;
    break ;
  case 3:
    color = fuchsia;
    break ;
  case 4:
    color = black;
    break ;
  case 5:
    color = green;
    break ;

  } 
	
  return color;	
}

