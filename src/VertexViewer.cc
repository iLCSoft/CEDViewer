#include "VertexViewer.h"
#include <EVENT/LCObject.h>
#include <EVENT/LCCollection.h>
#include <EVENT/SimCalorimeterHit.h>
#include <EVENT/CalorimeterHit.h>
#include <EVENT/SimTrackerHit.h>
#include <EVENT/TrackerHit.h>
#include <EVENT/Cluster.h>
#include <EVENT/Track.h>
#include <EVENT/ReconstructedParticle.h>
#include <EVENT/MCParticle.h>
#include <iostream>
#include "ClusterShapes.h"
#include <iomanip>
#include <UTIL/LCTOOLS.h>
#include <UTIL/LCRelationNavigator.h>
#include "HelixClass.h"
#include <math.h>
#include <ced_cli.h>

//#include "MarlinCED.h"

#include "marlin/Global.h"

//Gear
#include <gear/GEAR.h>
#include <gear/TPCParameters.h>
#include <gear/CalorimeterParameters.h>
#include <gear/LayerLayout.h>
#include <gear/VXDParameters.h>
#include <gear/VXDLayerLayout.h>
#include <gear/BField.h>

using namespace lcio ;
using namespace marlin ;

extern "C" {
  void tfithl_(int & NPT, double * XF,double * YF, float * RF, float *PF, double * WF,float * ZF,
	       float * WZF, int & IOPT, float * VV0,
	       float * EE0, float & CH2PH, float & CH2Z);
}

extern "C" {
  extern struct {
    int nplmat;          //total number of ladders
    float xplmat[200];   //X coordinate of the center of the ladder
    float yplmat[200];   //Y coordinate of the center of the ladder
    float zplmat[200];   //Z coordinate of the center of the ladder
    float widplmat[200]; //width of the ladder
    float lenplmat[200]; //length of the ladder (ladder parallel to z axis)
    float phiplmat[200]; //angle between normal vector of ladder and x axis
    float xrlpl[200];    //radiation length
    float xelospl[200];  //energy loss
  } fkddes2_; 
}


VertexViewer aVertexViewer ;


VertexViewer::VertexViewer() : Processor("VertexViewer") {

  _description = "Vertex Drawing Utility" ;
  
  registerProcessorParameter("TrueTrackCollection" , 
			     "True Track Collection Name" , 
			     _trueTracksCollection , 
			     std::string("TrueTracks"));

  registerProcessorParameter("TrueTrackMCPRelCollection" , 
			     "True Track MCP Relation Collection Name" , 
			     _trueTracksMCPCollection , 
			     std::string("TrueTracksMCP"));

  std::vector<std::string> collections;
  collections.push_back(std::string("VTXTrackerHits"));
  collections.push_back(std::string("FTDTrackerHits"));
  collections.push_back(std::string("SITTrackerHits"));
  collections.push_back(std::string("TPCTrackerHits"));

  std::vector<std::string> simCollections;
  simCollections.push_back(std::string("VXDCollection"));
  
  registerProcessorParameter("TrackerHitCollection" , 
			     "TrackerHit Collection Name" , 
			     _trackerHitCollection , 
			     collections);

  registerProcessorParameter("SimTrakerHitCollection",
			     "SimTrakerHit Collection Name",
			     _simTrackerHitCollection,
			     simCollections);
  
  registerProcessorParameter("TrackCollection" , 
			     "Track Collection Name" , 
			     _tracksCollection , 
			     std::string("SiTracks"));

  registerProcessorParameter("TrackMCPRelCollection" , 
			     "Track MCP Rel Collection Name" , 
			     _tracksMCPCollection , 
			     std::string("SiTracksMCP"));

  registerProcessorParameter("LayerTrueTracks", 
			     "Layer for True Tracks",
			     _layerTrueTracks, 
			     (int)-1);

  registerProcessorParameter("LayerTrackerHits", 
			     "Layer for TrackerHits",
			     _layerTrackerHits, 
			     (int)-1);

  registerProcessorParameter("LayerTracks", 
			     "Layer for Tracks",
			     _layerTracks, 
			     (int)-1);

  registerProcessorParameter("LayerSimTrackerHits", 
			     "Layer for SimTrackerHits",
			     _layerSimTrackerHits, 
			     (int)-1);

  registerProcessorParameter("CutOnD0",
			     "Cut On D0",
			     _cutOnD0,
			     float(1e+20));
  
  registerProcessorParameter("CutOnZ0",
			     "Cut On Z0",
			     _cutOnZ0,
			     float(1e+20));

  registerProcessorParameter("DetectorModel",
           "Detector Model",
           _detModel,
           (int)0);

}

void VertexViewer::init() {

    _nRun = -1;
    _nEvt = 0;

    //    MarlinCED::init(this) ;
    ced_client_init("localhost",7286);
    ced_register_elements();

}


void VertexViewer::processRunHeader( LCRunHeader* /*run*/) { 
  _nRun++ ;
  _nEvt = 0;
} 

void VertexViewer::processEvent( LCEvent * evt ) { 

    std::cout << std::endl;
    std::cout << " +++++++++++++++++++++++++++" << std::endl;
    std::cout << " Vertex Viewer : event number " << _nEvt << std::endl;
    std::cout << std::endl;

    //    _bField = Global::parameters->getFloatVal("BField");
    _bField = Global::GEAR->getBField().at( gear::Vector3D( 0., 0., 0.) ).z() ;

    
    ced_new_event();  

    //Draw Vertex detector
    if (_detModel == 1)
    {
    /*
      const gear::TPCParameters& gearTPC = Global::GEAR->getTPCParameters();
      const gear::CalorimeterParameters& gearEcalBarrel = Global::GEAR->getEcalBarrelParameters();
      const gear::CalorimeterParameters& gearEcalEndcap = Global::GEAR->getEcalEndcapParameters();
      const gear::LayerLayout& gearEcalBarrelLayerLayout = gearEcalBarrel.getLayerLayout();
      const gear::LayerLayout& gearEcalEndcapLayerLayout = gearEcalEndcap.getLayerLayout();
      const gear::CalorimeterParameters& gearHcalBarrel = Global::GEAR->getHcalBarrelParameters();
      const gear::CalorimeterParameters& gearHcalEndcap = Global::GEAR->getHcalEndcapParameters();
      const gear::LayerLayout& gearHcalBarrelLayerLayout = gearHcalBarrel.getLayerLayout();
      const gear::LayerLayout& gearHcalEndcapLayerLayout = gearHcalEndcap.getLayerLayout();
    */
  
      const gear::VXDParameters& vxdDetector = Global::GEAR->getVXDParameters();
      const gear::VXDLayerLayout& vxdLayerLayout = vxdDetector.getVXDLayerLayout();
  
      //const double deg2rad = Pi / 180.;
      //const double rad2deg = 180. / Pi;
  
  
      //Vertex detector
      int nLayersVTX = vxdLayerLayout.getNLayers();
      int nLadders = 0;
      float ladder_phi0, ladder_distance, ladder_thickness, ladder_width, ladder_length, ladder_offset;
      float sensitive_distance, sensitive_thickness;
      float x0,y0,x1,y1;
  
      for (int i=0; i<nLayersVTX; i++) {
        nLadders = vxdLayerLayout.getNLadders(i);
  
        ladder_phi0 = float(vxdLayerLayout.getPhi0(i));
        ladder_distance = float(vxdLayerLayout.getLadderDistance(i));
        ladder_thickness = float(vxdLayerLayout.getLadderThickness(i));
        ladder_width = float(vxdLayerLayout.getLadderWidth(i));
        ladder_length = float (vxdLayerLayout.getLadderLength(i));
        ladder_offset = float (vxdLayerLayout.getLadderOffset(i));
  
        sensitive_distance = float(vxdLayerLayout.getSensitiveDistance(i));
        sensitive_thickness = float(vxdLayerLayout.getSensitiveThickness(i));
  
        float ft_xplmat,ft_yplmat;
  
        ladder_distance += 0.5* ladder_thickness;
        sensitive_distance += 0.5* sensitive_thickness;

        ladder_length = 10.0*ladder_length;
  
        if (nLadders >= 3)
        {
          //Realistic ladder structure
          float currPhi;
          float angleLadders = 2*M_PI/ nLadders;
          float cosphi, sinphi;
  
          for (int j=0; j<nLadders; ++j) {
  
	    currPhi = ladder_phi0 + (angleLadders * j);
            cosphi = cos(currPhi);
            sinphi = sin(currPhi);
  
            ft_xplmat = (sensitive_distance*cosphi - ladder_offset*sinphi);
            ft_yplmat = (sensitive_distance*sinphi + ladder_offset*cosphi);
  
            x0 = 10.0*(ft_xplmat + (sinphi*0.5*ladder_width));
            y0 = 10.0*(ft_yplmat - (cosphi*0.5*ladder_width));
            x1 = 10.0*(ft_xplmat - (sinphi*0.5*ladder_width));
            y1 = 10.0*(ft_yplmat + (cosphi*0.5*ladder_width));
  
            ced_line(x0,y0,0.0,x1,y1,0.0,0,0,0xffffff); 
            ced_line(x0,y0,-ladder_length,x1,y1,-ladder_length,0,0,0xffffff);
            ced_line(x0,y0,ladder_length,x1,y1,ladder_length,0,0,0xffffff);
            ced_line(x0,y0,-ladder_length,x0,y0,ladder_length,0,0,0xffffff);
            ced_line(x1,y1,-ladder_length,x1,y1,ladder_length,0,0,0xffffff);
          }
        }
        else
        {
          ced_geocylinder(10.0*ladder_distance,24,0.0,ladder_length,-ladder_length,0xffffff);
        }
      }
    }

//---------------------------------------------------------------------------

    LCCollection * MCPartCol = evt->getCollection("MCParticle");
    int nPart = MCPartCol->getNumberOfElements();
    //     std::cout << " PDG     Px        Py        Pz       Energy" << std::endl;
    //     //           "  11   126.278    57.085   102.390     0.000"
    for (int iP=0;iP<nPart;++iP) {
      MCParticle * part = dynamic_cast<MCParticle*>(MCPartCol->getElementAt(iP));
      //       if (part->getEnergy()>5.0&&!(part->isBackscatter())) {
      // 	printf("%5i  %8.3f  %8.3f  %8.3f  %8.3f\n",
      // 	       part->getPDG(),
      // 	       part->getMomentum()[0],
      // 	       part->getMomentum()[1],
      // 	       part->getMomentum()[2],
      // 	       part->getEnergy());
      int pdg = part->getPDG();
      if (fabs(pdg)>=400 && fabs(pdg)<600) {
	std::cout << " PDG = " << pdg 
		  << "(x,y,z) = " 
		  << part->getEndpoint()[0] << " "
		  << part->getEndpoint()[1] << " "
		  << part->getEndpoint()[2] << std::endl;
	//float x = part->getEndpoint()[0];
	//float y = part->getEndpoint()[1];
	//float z = part->getEndpoint()[2];	 
	//	ced_hit(20.*x,20.*y,20.*z,_layerTrueTracks<<CED_LAYER_SHIFT,10,0xffffff);
      }
    }
    
    std::cout << std::endl;
    std::map<LCObject*,Track*> ObjTrack;
    if (_layerTracks >=0 ) {
      try {
	LCCollection * col = evt->getCollection(_tracksCollection.c_str());
	LCCollection * relCol = evt->getCollection(_tracksMCPCollection.c_str());
	LCRelationNavigator navigator(relCol);
	int nelem = col->getNumberOfElements();
	for (int iclust(0); iclust < nelem; ++iclust) {
	  Track * track = 
	    dynamic_cast<Track*>(col->getElementAt(iclust));
	  LCObjectVec objVec = navigator.getRelatedToObjects(track);
	  FloatVec weights = navigator.getRelatedToWeights(track);
	  LCObject * obj = NULL;
	  float wght = 0;
	  int nObj = objVec.size();
	  for (int iObj=0;iObj<nObj;++iObj) {
	    if (weights[iObj]>wght) {
	      wght = weights[iObj];
	      obj = objVec[iObj];
	    }
	  }
	  ObjTrack[obj] = track;
	}		
      }
      catch(DataNotAvailableException &e) {};
    }

// Drawing Tracker Hits 
    if (_layerTrackerHits >= 0 ) {
      int nCol = int(_trackerHitCollection.size());
      for (int iCol=0;iCol<nCol;++iCol) {
	try {
	  LCCollection * col = 
	    evt->getCollection(_trackerHitCollection[iCol].c_str()) ;
	  int nelem = col->getNumberOfElements();
	  std::cout << "Collection " << _trackerHitCollection[iCol].c_str() 
		    << " ; # of hits = " << nelem << std::endl; 
	  for (int ielem = 0; ielem < nelem; ++ielem) {
	    TrackerHit * hit = 
	      dynamic_cast<TrackerHit*>( col->getElementAt(ielem));
	    float x = (float)hit->getPosition()[0];
	    float y = (float)hit->getPosition()[1];
	    float z = (float)hit->getPosition()[2];
      //std::cout << x << " " << y << " " << z << std::endl;
	    ced_hit(10.*x,10.*y,10.*z, 
		    _layerTrackerHits<<CED_LAYER_SHIFT,6,0x7cf774);
	  }	
	}
	catch(DataNotAvailableException &e) {}
      }
      ced_send_event();
    }


// Drawing True Tracks
    std::map<LCObject*,int> McpInt; 
    if (_layerTrueTracks >= 0) {
	try {

	    LCCollection * col = 
	      evt->getCollection(_trueTracksCollection.c_str());
	    LCCollection * relCol = 
	      evt->getCollection(_trueTracksMCPCollection.c_str());
	    LCRelationNavigator navigator(relCol);
	    int nelem = col->getNumberOfElements();
	    std::cout << std::endl;
	    std::cout << "++++++++++++++++++++++++++++++++++++" << std::endl;
	    std::cout << "      DISPLAYING TRUE TRACKS " << std::endl;
	    std::cout << "++++++++++++++++++++++++++++++++++++" << std::endl;
	    std::cout << "# of True Tracks : " << nelem << std::endl;
	    printf("Track   VTX   FTD   SIT   TPC       D0        Z0        Px        Py        Pz       PDG  \n");
	    printf("--------------------------------------------------------------------------------------------------------------\n");
	    float totPx=0.0;
	    float totPy=0.0;
	    float totPz=0.0;
	    float totEn=0.0;
	    float totPxN=0.0;
	    float totPyN=0.0;
	    float totPzN=0.0;
	    float totEnN=0.0;	    
	    for (int iclust(0); iclust < nelem; ++iclust) {
	      Track * track = 
		dynamic_cast<Track*>(col->getElementAt(iclust));
	      LCObjectVec objVec = navigator.getRelatedToObjects(track);
	      FloatVec weights = navigator.getRelatedToWeights(track);
	      LCObject * obj = NULL;
	      float wght = 0;
	      int nObj = objVec.size();
	      for (int iObj=0;iObj<nObj;++iObj) {
		if (weights[iObj]>wght) {
		  wght = weights[iObj];
		  obj = objVec[iObj];
		}
	      }
	      McpInt[obj] = iclust;
	      int color = returnColor(iclust);
	      TrackerHitVec hitvec = track->getTrackerHits();
	      int nhits = (int)hitvec.size();
	      float * ah = new float[nhits];
	      float * xh = new float[nhits];
	      float * yh = new float[nhits];
	      float * zh = new float[nhits];
	      float zmin = 1.0E+10;
	      float zmax = -1.0E+10;
	      for (int ihit(0); ihit < nhits; ++ihit) {
		TrackerHit * hit = hitvec[ihit];
		float x = (float)hit->getPosition()[0];
		float y = (float)hit->getPosition()[1];
		float z = (float)hit->getPosition()[2];
		ced_hit(10.*x,10.*y,10.*z,_layerTrueTracks<<CED_LAYER_SHIFT,6,color);
		ah[ihit] = 1.0;
		xh[ihit] = x;
		yh[ihit] = y;
		zh[ihit] = z;
		if (z < zmin)
		  zmin = z;
		if (z > zmax)
		  zmax = z;
	      } 	
	      ClusterShapes * shapes = new ClusterShapes(nhits,ah,xh,yh,zh);
	      //float zBegin, zEnd;
	      if (fabs(zmin)<fabs(zmax)) {
		//zBegin = zmin;
		//zEnd   = zmax;
	      }
	      else {
		//zBegin = zmax;
		//zEnd   = zmin;
	      }
	      //	      float signPz = zEnd - zBegin;		
	      //	      float dz = (zmax - zmin) / 100.;
	      float par[5];
	      float dpar[5];
	      float chi2;
	      float distmax;
	      shapes->FitHelix(500, 0, 1, par, dpar, chi2, distmax);
	      float x0 = par[0];
	      float y0 = par[1];
	      float r0 = par[2];
	      //	      float bz = par[3];
	      float phi0 = par[4];
	      HelixClass * helix = new HelixClass();
	      MCParticle * part = dynamic_cast<MCParticle*>(obj);	      
	      float mom[3];
	      float ver[3];
	      for (int iC=0;iC<3;++iC) {
		mom[iC] = float(part->getMomentum()[iC]);
		ver[iC] = float(part->getVertex()[iC]);
	      }
	      //	      float vv = sqrt(ver[0]*ver[0]+ver[1]*ver[1]);
		
	      float charge = part->getCharge();
	      helix->Initialize_VP(ver,mom,charge,_bField);
	      //	      helix->Initialize_Canonical(track->getPhi(),track->getD0(),track->getZ0(),
	      //					  track->getOmega(),track->getTanLambda(),_bField);
	      zmin = float(part->getVertex()[2]);
	      zmax = float(part->getEndpoint()[2]);
	      float tmax = (zmax - zmin)/mom[2];
	      float dt = tmax/1000;
	      x0 = float(part->getVertex()[0]);
	      y0 = float(part->getVertex()[1]);
	      r0 = helix->getRadius();
	      float pxy = sqrt(mom[0]*mom[0]+mom[1]*mom[1]);
	      float xc = helix->getXC();
	      float yc = helix->getYC();
	      phi0 = atan2(y0-yc,x0-xc);
	      if (fabs(track->getD0())<_cutOnD0&&fabs(track->getZ0())<_cutOnZ0) {
		int nHitsVTX = track->getSubdetectorHitNumbers()[0];
		int nHitsFTD = track->getSubdetectorHitNumbers()[1];
		int nHitsSIT = track->getSubdetectorHitNumbers()[2];
		int nHitsTPC = track->getSubdetectorHitNumbers()[3];	    		
 		printf(" %2i     %3i   %3i   %3i  %4i  %9.3f %9.3f   %7.2f   %7.2f   %7.2f   %5i\n",
 		       iclust, nHitsVTX, nHitsFTD, nHitsSIT, nHitsTPC, track->getD0(), track->getZ0(), 
 		       part->getMomentum()[0],part->getMomentum()[1],part->getMomentum()[2],part->getPDG());
		// 		std::cout << "     " << part ;
		if (ObjTrack[obj] == NULL)
		  std::cout << "   *** " ;
		std::cout << std::endl;
		totPx +=  helix->getMomentum()[0];
		totPy +=  helix->getMomentum()[1];
		totPz +=  helix->getMomentum()[2];	
		float Energy = sqrt(helix->getMomentum()[0]*helix->getMomentum()[0]+
				    helix->getMomentum()[1]*helix->getMomentum()[1]+
				    helix->getMomentum()[2]*helix->getMomentum()[2]);
		totEn += Energy;
		if (ObjTrack[obj]==NULL) {
		  totPxN +=  helix->getMomentum()[0];
		  totPyN +=  helix->getMomentum()[1];
		  totPzN +=  helix->getMomentum()[2];	
		  totEnN += Energy;
		}
	      }      

	      //	      if (chi2 > 0. && chi2 < 10000.) {
		for (int it=0; it < 1000; ++it) {
		  float t1 = it*dt;
		  float t2 = t1 + dt;
		  float z1 = zmin + t1*mom[2];
		  float z2 = zmin + t2*mom[2];
		  float x1 = xc + r0*cos(-charge*pxy*t1/r0+phi0);
		  float y1 = yc + r0*sin(-charge*pxy*t1/r0+phi0);
		  float x2 = xc + r0*cos(-charge*pxy*t2/r0+phi0);
		  float y2 = yc + r0*sin(-charge*pxy*t2/r0+phi0);			
		  ced_line(10.*x1,10.*y1,10.*z1,10.*x2,10.*y2,10.*z2,
			   _layerTrueTracks<<CED_LAYER_SHIFT,2,color);
		  
		}		
		//	      }
	      delete shapes;
	      delete helix;
	      delete[] xh;
	      delete[] yh;
	      delete[] zh;
	      delete[] ah;
	      ced_send_event();	      
	    }
	    printf("===========================================================================================================\n");
	    std::cout << "Statistics (True Tracks) ---> " << std::endl;
	    printf("Total  4P --> Px = %7.2f  Py = %7.2f  Pz = %7.2f  E = %7.2f\n",
		   totPx,totPy,totPz,totEn);
	    printf("Missed 4P --> Px = %7.2f  Py = %7.2f  Pz = %7.2f  E = %7.2f\n",
		   totPxN,totPyN,totPzN,totEnN);




	}
	catch(DataNotAvailableException &e){}	
    }


// Drawing Reconstructed Tracks
    if (_layerTracks >= 0) {
	try {
	  LCCollection * col = evt->getCollection(_tracksCollection.c_str());
	  LCCollection * relCol = evt->getCollection(_tracksMCPCollection.c_str());
	  LCRelationNavigator navigator(relCol);
	  int nelem = col->getNumberOfElements();
	  std::cout << std::endl;
	  std::cout << "++++++++++++++++++++++++++++++++++++" << std::endl;
	  std::cout << "  DISPLAYING RECONSTRUCTED TRACKS   " << std::endl;
	  std::cout << "++++++++++++++++++++++++++++++++++++" << std::endl;
	  std::cout << "# of Reconstructed Tracks = " << nelem << std::endl;
	  printf(" Trk VTX  FTD  SIT     TPC         D0        Z0      Px      Py      Pz      L    W     chi2/ndf\n");
	  //     "  0  5(5) 0(0) 2(2)   96(458)    -0.014    -0.008    0.36    0.47    0.06    0 1.000   1.198     0xd6bf5f0"
	  printf("-----------------------------------------------------------------------------------------------------------------------\n");
	  float totPx=0.0;
	  float totPy=0.0;
	  float totPz=0.0;
	  float totEn=0.0;
	  for (int iclust(0); iclust < nelem; ++iclust) {
	    Track * track = 
	      dynamic_cast<Track*>(col->getElementAt(iclust));
	    LCObjectVec objVec = navigator.getRelatedToObjects(track);
	    FloatVec weights = navigator.getRelatedToWeights(track);
	    LCObject * obj = NULL;
	    float wght = 0;
	    int nObj = objVec.size();
	    for (int iObj=0;iObj<nObj;++iObj) {
	      if (weights[iObj]>wght) {
		wght = weights[iObj];
		obj = objVec[iObj];
	      }
	    }
	    TrackerHitVec hitvec = track->getTrackerHits();
	    int color = returnColor(iclust+5);
	    int nhits = (int)hitvec.size();
	    float * ah = new float[nhits];
	    float * xh = new float[nhits];
	    float * yh = new float[nhits];
	    float * zh = new float[nhits];
	    float zmin = 1.0E+10;
	    float zmax = -1.0E+10;
  	    int nHitsVTX = track->getSubdetectorHitNumbers()[0];
 	    int nHitsFTD = track->getSubdetectorHitNumbers()[1];
 	    int nHitsSIT = track->getSubdetectorHitNumbers()[2];
 	    int nHitsTPC = track->getSubdetectorHitNumbers()[3];
	    //int nHitsSET = track->getSubdetectorHitNumbers()[4];
	    //int nHitsETD = track->getSubdetectorHitNumbers()[5];
 	    int nHitsVTXInFit = track->getSubdetectorHitNumbers()[6];
 	    int nHitsFTDInFit = track->getSubdetectorHitNumbers()[7];
 	    int nHitsSITInFit = track->getSubdetectorHitNumbers()[8];
 	    int nHitsTPCInFit = track->getSubdetectorHitNumbers()[9];
 	    //int nHitsSETInFit = track->getSubdetectorHitNumbers()[10];
 	    //int nHitsETDInFit = track->getSubdetectorHitNumbers()[11];

	    float xprev = 0.;
	    float yprev = 0.;
	    float zprev = 0.;
	    HelixClass * helix = new HelixClass();
	    helix->Initialize_Canonical(track->getPhi(),track->getD0(),track->getZ0(),
					track->getOmega(),track->getTanLambda(),_bField);
	    float PX = helix->getMomentum()[0];
	    float PY = helix->getMomentum()[1];
	    float PT = sqrt(PX*PX+PY*PY);
	    float rmax = 0;
	    for (int ihit(0); ihit < nhits; ++ihit) {
	      TrackerHit * hit = hitvec[ihit];
	      float x = (float)hit->getPosition()[0];
	      float y = (float)hit->getPosition()[1];
	      float z = (float)hit->getPosition()[2];
	      ced_hit(10.*x,10.*y,10.*z,_layerTracks<<CED_LAYER_SHIFT,6,color);
	      ced_line(10.*x,10.*y,10.*z,10.*xprev,10.*yprev,10.*zprev,_layerTracks<<CED_LAYER_SHIFT,2,color);
	      if (PT<1.2&&nHitsTPC>180) {
		ced_hit(10.*x,10.*y,10.*z,_layerTracks<<CED_LAYER_SHIFT,6,color);
	      }
	      float rH = sqrt(x*x+y*y);
	      if (rH>rmax)
		rmax = rH;
	      xprev = x;
	      yprev = y;
	      zprev = z;
	      ah[ihit] = 1.0;
	      xh[ihit] = x;
	      yh[ihit] = y;
	      zh[ihit] = z;
	      if (z < zmin)
		zmin = z;
	      if (z > zmax)
		zmax = z;
	    } 	
	    ClusterShapes * shapes = new ClusterShapes(nhits,ah,xh,yh,zh);
	    //float zBegin, zEnd;
	    if (fabs(zmin)<fabs(zmax)) {
	      //zBegin = zmin;
	      //zEnd   = zmax;
	    }
	    else {
	      //zBegin = zmax;
	      //zEnd   = zmin;
	    }
	    //	    float signPz = zEnd - zBegin;		
	    //float dz = (zmax - zmin) / 100.;
	    float par[5];
	    float dpar[5];
	    float chi2;
	    float distmax;
	    shapes->FitHelix(500, 0, 1, par, dpar, chi2, distmax);
	    //float x0 = par[0];
	    //float y0 = par[1];
	    //float r0 = par[2];
	    //float bz = par[3];
	    //float phi0 = par[4];
	    if (fabs(track->getD0())<_cutOnD0&&fabs(track->getZ0())<_cutOnZ0) {
	      float ndf = float(track->getNdf());
	      if (ndf<1.0)
		ndf = 1.0;
	      printf(" %2i  %1i(%1i) %1i(%1i) %1i(%1i) %3i(%4i) %9.3f %9.3f %7.2f %7.2f %7.2f  %3i   %4.2f  %7.3f\n",
		     iclust, 
		     nHitsVTX, nHitsVTXInFit,nHitsFTD,nHitsFTDInFit,
		     nHitsSIT, nHitsSITInFit,nHitsTPC,nHitsTPCInFit,
		     track->getD0(), track->getZ0(), 
		     helix->getMomentum()[0],helix->getMomentum()[1],helix->getMomentum()[2],
		     McpInt[obj],wght,track->getChi2()/ndf);	      
	      //	      std::cout << "     " << obj << std::endl;
	      //	      if (PT<1.2&&nHitsTPC>180)
	      //		std::cout << " * rmax = " << rmax << std::endl;
	      //	      else 
	      //		std::cout << std::endl;
	      totPx +=  helix->getMomentum()[0];
	      totPy +=  helix->getMomentum()[1];
	      totPz +=  helix->getMomentum()[2];	
	      float Energy = sqrt(helix->getMomentum()[0]*helix->getMomentum()[0]+
				  helix->getMomentum()[1]*helix->getMomentum()[1]+
				  helix->getMomentum()[2]*helix->getMomentum()[2]);
	      totEn += Energy;	      
	    }      
 	      if (nHitsVTX>5 || nHitsFTD>7 || nHitsSIT > 2) {
 		for (int ihit=0;ihit<nhits;++ihit) {
		  int det = hitvec[ihit]->getType()/100;
		  if (det <=4 )
		    std::cout << ihit << " "  << hitvec[ihit]->getType() << " " 
			      << hitvec[ihit]->getPosition()[0] << " " 
			      << hitvec[ihit]->getPosition()[1] << " " 
			      << hitvec[ihit]->getPosition()[2] << " "  << std::endl;
 		}
 	      }
          
	      if (chi2 > 0. && chi2 < 10000.) {
		//for (int iz(0); iz < 100; ++iz) {
		//  float z1 = zmin + iz*dz;
		//  float z2 = z1 + dz;
		//  float x1 = x0 + r0*cos(bz*z1+phi0);
		//  float y1 = y0 + r0*sin(bz*z1+phi0);
		//  float x2 = x0 + r0*cos(bz*z2+phi0);
		//  float y2 = y0 + r0*sin(bz*z2+phi0);			
		//  //		  ced_line(10.*x1,10.*y1,10.*z1,10.*x2,10.*y2,10.*z2,_layerTracks<<CED_LAYER_SHIFT,1,color);
		//}
		delete shapes;
		delete helix;
		delete[] xh;
		delete[] yh;
		delete[] zh;
		delete[] ah;
	      }
	      ced_send_event();
	      //	      getchar();
	  }
	  printf("======================================================================================================================\n");
	  std::cout << "Statistics (Reco Tracks) ---> " << std::endl;
	  printf("Total 4P --> Px = %7.2f  Py = %7.2f  Pz = %7.2f  E = %7.2f\n",
		 totPx,totPy,totPz,totEn);
	}
	catch(DataNotAvailableException &e){}	
    }

    ced_send_event();
    std::cout << std::endl;

    // Drawing Sim Tracker Hits ---->
    if (_layerSimTrackerHits >= 0 ) {
      int nCol = int(_simTrackerHitCollection.size());
      for (int iCol=0;iCol<nCol;++iCol) {
	try {
	  LCCollection * col = 
	    evt->getCollection(_simTrackerHitCollection[iCol].c_str()) ;
	  int nelem = col->getNumberOfElements();
	  std::cout << "SimTrackerHit Collection " << _simTrackerHitCollection[iCol].c_str() 
		    << " ; # of hits = " << nelem << std::endl; 
	  for (int ielem = 0; ielem < nelem; ++ielem) {
	    SimTrackerHit * hit = 
	      dynamic_cast<SimTrackerHit*>( col->getElementAt(ielem));
	    float x = (float)hit->getPosition()[0];
	    float y = (float)hit->getPosition()[1];
	    float z = (float)hit->getPosition()[2];
	    ced_hit(10.*x,10.*y,10.*z, 
		    _layerSimTrackerHits<<CED_LAYER_SHIFT,6,0xffbd22);
	  }	
	}
	catch(DataNotAvailableException &e) {}
      }
      ced_send_event();
    }

    std::cout << std::endl;
    std::cout << "================================================" << std::endl;
    std::cout << "== Vertex Viewer : Press any key to continue --"  << std::endl;
    std::cout << "================================================" << std::endl;    
    getchar();
    std::cout << std::endl;

//++++++++++++++++++++++++++++++++++++
//  MarlinCED::draw(this) ;
//++++++++++++++++++++++++++++++++++++

    _nEvt++;

}


void VertexViewer::check( LCEvent * /*evt*/ ) { }
  
void VertexViewer::end(){ } 

int VertexViewer::returnColor(int counter) {

	int icol =  counter % 32;
	int kcol =  0x000000;
	if (icol==0)  kcol = 0x00ff00;
	if (icol==1)  kcol = 0xAA00ff;
	if (icol==2)  kcol = 0xff0000;
	if (icol==3)  kcol = 0x00ffff;
	if (icol==4)  kcol = 0xffff00;
	if (icol==5)  kcol = 0xff00ff;
	if (icol==6)  kcol = 0xffffff;
	if (icol==7)  kcol = 0xf0f0f0;
	if (icol==8)  kcol = 0x000088;
	if (icol==9)  kcol = 0x008800;
	if (icol==10) kcol = 0x880000;
	if (icol==11) kcol = 0x008888;
	if (icol==12) kcol = 0x888800;
	if (icol==13) kcol = 0x880088;
	if (icol==14) kcol = 0x888888;
	if (icol==15) kcol = 0x00A888;
	if (icol==16) kcol = 0xffefd5;
	if (icol==17) kcol = 0xffdab9;
	if (icol==18) kcol = 0xffe4e1;
	if (icol==19) kcol = 0x6495ed;
	if (icol==20) kcol = 0x00bfff;
	if (icol==21) kcol = 0x40e0d0;
	if (icol==22) kcol = 0x7cfc00;
	if (icol==23) kcol = 0xffefd5;	
	if (icol==24) kcol = 0xcd5c5c;
	if (icol==25) kcol = 0xff1493;
	if (icol==26) kcol = 0xc0ff3e;
	if (icol==27) kcol = 0xff8247;
	if (icol==28) kcol = 0xcd6090;
	if (icol==29) kcol = 0xff34b3;
	if (icol==30) kcol = 0xcdad00;
	if (icol==31) kcol = 0xee6363;

	return kcol;

}
