#include "HEP_RecordProcessor.h"

#include "Phys_Geom_Database.h"
#include "layers.h"

#include <iostream>
#include <cmath>

#include "lcio.h"

#include "EVENT/LCCollection.h"
#include "EVENT/MCParticle.h"

// #include <ced_cli.h>
#include "MarlinCED.h"
#include <gear/GEAR.h>
#include <gear/TPCParameters.h>
#include <gearimpl/TPCParametersImpl.h>
#include <gearimpl/FixedPadSizeDiskLayout.h>


namespace marlin{


  void MC_Particles_Balance(LCEvent * evt);
  void print_MC_Particles(LCEvent * evt, float ecut);

  void draw_MC_Photons(LCEvent* evt, float ecut);
  void draw_MC_Neutral_Hadrons(LCEvent* evt, float ecut);
  void draw_MC_Charged_Hadrons(LCEvent * evt, float ecut);

  double My_abs(double);
  double My_abs(double a){return a < 0.0 ? -a: a;}

  int    My_iabs(int);
  int    My_iabs(int a){return a < 0.0 ? -a: a;}

// ?????   7 should be a global constant  ????
//  PGdb HEP_PGDB(7); 
  PGdb HEP_PGDB; 
    
  HEP_RecordProcessor aHEP_RecordProcessor ;
  HEP_RecordProcessor::HEP_RecordProcessor() : Processor("HEP_RecordProcessor") {
    _description = "HEP record balance and drawing" ;
  }

  void HEP_RecordProcessor::init() { 

  //    Initialize event display   
  MarlinCED::init(this) ;
//fg     ced_client_init("localhost",7286);
//fg     ced_register_elements();
    //    std::cout << "HEP_RecordProcessor::init()  " << name() << std::endl 
    //	      << "  parameters: " << std::endl << *parameters()  ;
 
    _nRun = 0 ;
    _nEvt = 0 ;
}

  void HEP_RecordProcessor::processRunHeader( LCRunHeader* run) { 
    std::cout << "BoojumProcessor::processRun()  " << name() 
	      << " in run " << run->getRunNumber() << std::endl ;
    _nRun++ ;
  } 
//-----------------------------------------------------------------------
//-----------------------------------------------------------------------
void HEP_RecordProcessor::processEvent( LCEvent * evt ) { 
//-----------------------------------------------------------------------
//-----------------------------------------------------------------------
//-----------------------------------------------------------------------
// Reset drawing buffer and START drawing collection
  MarlinCED::newEvent(this) ;
//   ced_new_event();  
//-----------------------------------------------------------------------
//fg   ced_new_event();  // Reset drawing buffer and START drawing collection

// fg - done in MarlinCED::newEvent(this) :
//   //          The Siplest geometry 
//   static CED_GeoCylinder geoCylinders[] = {       // for TESLA Detector Geometry
//     {    50.0,  3,  0.0, 5658.5, -5658.5, 0xff     }, // beam tube
//     {   380.0,  4,  0.0, 2658.5, -2658.5, 0xff     }, // inner TPC
//     {  1840.0,  8, 22.5, 2700.0, -2700.0, 0x7f7f1f }, // inner ECAL
//     {  2045.7,  8, 22.5, 2700.0, -2700.0, 0x7f7f1f }, // outer ECAL
//     {  2045.7,  8, 22.5, 101.00,  2820.0, 0x7f7f1f }, // endcap ECAL
//     {  2045.7,  8, 22.5, 101.00, -3022.0, 0x7f7f1f }, // endcap ECAL
//     {  3000.0, 16,  0.0, 2658.5, -2658.5, 0xcf00   }, // outer HCAL
//     {  3000.0,  8, 22.5, 702.25,  2826.0, 0xcf00   }, // endcap HCAL
//     {  3000.0,  8, 22.5, 702.25, -4230.5, 0xcf00   }, // endcap HCAL
//   }; 
//   /*
//   static CED_GeoCylinder geoCylinders[] = {    // for Prototype
//     {    180.0,  4,  45.0, 110.0, 0.0, 0xff }, // beam tube
//     {    500.0,  4,  45.0, 250.0, 220., 0xff } // inner TPC
//   };
//   */
//   ced_geocylinders(sizeof(geoCylinders)/sizeof(CED_GeoCylinder),geoCylinders);
 
//-----------------------------------------------------------------------

//fg   ced_send_event();    // Does not clean up screen

//  std::cout << std::endl
//	    << "        ====  HEP record for event "<<evt->getEventNumber()<<" ==============" << std::endl ;
  print_MC_Particles(evt, 0.02);
  MC_Particles_Balance(evt);
  
  draw_MC_Photons        (evt,0.2);  //  Draw Photons from HEP record
  draw_MC_Neutral_Hadrons(evt,0.2);  //  Draw Neutral Hadrons from HEP record
  draw_MC_Charged_Hadrons(evt,0.2);  //  Draw Charged Hadrons from HEP record
  
//fg   getchar();

//++++++++++++++++++++++++++++++++++++
  MarlinCED::draw(this) ;
  //++++++++++++++++++++++++++++++++++++
//   ced_draw_event();    // Make clean up screen before drawing

  _nEvt ++ ;
  
  return ;  // from HEP_RecordProcessor::processEvent( LCEvent * evt ) 
}
  
void HEP_RecordProcessor::check( LCEvent * evt ) { 
  std::cout << "HEP_RecordProcessor::check()  " << name() 
	    << " in event " << evt->getEventNumber() << " (run " 
	    << evt->getRunNumber() << ") "
	    << std::endl ;
}

void HEP_RecordProcessor::end(){ 
  std::cout << "HEP_RecordProcessor::end()  " << name() 
	    << " processed " << _nEvt << " events in " << _nRun << " runs "
	    << std::endl ;
}

//-----------------------------------------------------------------------
void draw_MC_Photons(LCEvent * evt, float ecut){
//-----------------------------------------------------------------------

//   double TPC_z           = HEP_PGDB.p_tpc_z_size();
//   double TPC_out         = HEP_PGDB.p_tpc_out_rad();
  
  const gear::TPCParameters  &pTPC      = Global::GEAR->getTPCParameters();
  const gear::PadRowLayout2D &padLayout = pTPC.getPadLayout();
  const gear::DoubleVec      &planeExt  = padLayout.getPlaneExtent();
  double TPC_z           = pTPC.getMaxDriftLength() ;
  double TPC_out         = planeExt[1]; 

  double add =  300.;          // Length of the photon "tracks"
  double tz = TPC_z   + 250.;  // to start close to ECAL surface
  double to = TPC_out + 100.;  // to start close to ECAL surface
  
  double x1,x2,y1,y2,z1,z2;
  double px,py,pz,pt,pr, cx,cy,cz, rsq,d;
  int idpdg;
  const double* mom;
  float enr;
  
  LCCollection* mcpCol = evt->getCollection("MCParticle" ) ;  
  for(int i=0; i<mcpCol->getNumberOfElements() ; i++){
    MCParticle* imc = dynamic_cast<MCParticle*> ( mcpCol->getElementAt( i ) ) ;
    idpdg = imc-> getPDG (); 
    mom = imc-> getMomentum (); 
    enr = imc-> getEnergy (); 
    //           Draw Pi0 s -- even if they are not the stable, (was decayed by simulator)
    if(idpdg == 111) { //  Pi_0s    
	if (enr > ecut) {//  All Pi0s with E > 20 MeV
	  add = 100. + 2000.*log(1.+enr);
	  /*
	     std::cout <<" PI_0 # "<< i <<'\t'<<" Energy = "<<enr<<'\t'
	     <<",   Mom = "<<mom[0]<<",  "<<mom[1]<<",  "<<mom[2]
	     << std::endl;
	  */
	  px = mom[0]; 
	  py = mom[1]; 
	  pz = mom[2];
	  pt = hypot(px,py);	 
	  pr = hypot(pt,pz);	 
	  //	  pr = sqrt(px*px + py*py + pz*pz);	 
	  // Direct Cosines
	  cx =  px/pr;
	  cy =  py/pr;
	  cz =  pz/pr;
	  x1 = 10.*cx;
	  y1 = 10.*cy;
	  z1 = 10.*cz;
	  x2 = x1 + add*cx;
	  y2 = y1 + add*cy;
	  z2 = z1 + add*cz;
	  //           yellow = 0xf9f920,   pink = 0xb900de
	  ced_line(x1,y1,z1,x2,y2,z2, NHEP_LAYER , 1, 0xf9f920);	 
	}
    }
    if( imc-> getGeneratorStatus() == 1 ) { // stable particles only   
    //           Draw Pi0 s   if they were not decayed by simulator
      if(idpdg == 111) { //  Pi_0s    
	if (enr > ecut) {//  All Pi0s with E > 20 MeV
	  add = 100. + 2000.*log(1.+enr);
	  /*
	     std::cout <<" PI_0 # "<< i <<'\t'<<" Energy = "<<enr<<'\t'
	     <<",   Mom = "<<mom[0]<<",  "<<mom[1]<<",  "<<mom[2]
	     << std::endl;
	  */
	  px = mom[0]; 
	  py = mom[1]; 
	  pz = mom[2];
	  pt = hypot(px,py);	 
	  pr = hypot(pt,pz);	 
	  //	  pr = sqrt(px*px + py*py + pz*pz);	 
	  // Direct Cosines
	  cx =  px/pr;
	  cy =  py/pr;
	  cz =  pz/pr;
	  x1 = 10.*cx;
	  y1 = 10.*cy;
	  z1 = 10.*cz;
	  x2 = x1 + add*cx;
	  y2 = y1 + add*cy;
	  z2 = z1 + add*cz;
	  //           yellow = 0xf9f920,   pink = 0xb900de
	  ced_line(x1,y1,z1,x2,y2,z2, NHEP_LAYER , 1, 0xf9f920);	 
	}
    }
      if(idpdg == 22) { 
	if (enr > ecut) {//  All gammas with E > 20 MeV
	  add = 10. + 300.*log(1.+enr);
	  /*
	     std::cout <<" Photon # "<< i <<'\t'<<" Energy = "<<enr<<'\t'
	     <<",   Mom = "<<mom[0]<<",  "<<mom[1]<<",  "<<mom[2]
	     << std::endl;
	  */
	  px = mom[0]; 
	  py = mom[1]; 
	  pz = mom[2];
	  pt = hypot(px,py);	 
	  pr = hypot(pt,pz);	 
	  //	  pr = sqrt(px*px + py*py + pz*pz);	 
	  // Direct Cosines
	  cx =  px/pr;
	  cy =  py/pr;
	  cz =  pz/pr;
	  //   Full length
	  d = tz/My_abs(cz);
	  //     2-D radius 
	  rsq = d*hypot(cx,cy);
	  if (rsq > to)       // start point at the cylinder of TPC
	    d = to/hypot(cx,cy);
	  x1 = d*cx;
	  y1 = d*cy;
	  z1 = d*cz;
	  x2 = x1 + add*cx;
	  y2 = y1 + add*cy;
	  z2 = z1 + add*cz;
	  //           yellow = 0xf9f920,   pink = 0xb900de
	  ced_line(x1,y1,z1,x2,y2,z2, NHEP_LAYER , 2, 0xf9f920);	 
	}
      }
    }
  }
//fg   ced_send_event();    // Does not clean up screen
  
  return;
}   //  end draw_MC_Photons
  
//-----------------------------------------------------------------------
void draw_MC_Neutral_Hadrons(LCEvent * evt, float ecut){
//-----------------------------------------------------------------------
//   double TPC_z           = HEP_PGDB.p_tpc_z_size();
//   double TPC_out         = HEP_PGDB.p_tpc_out_rad();
  
  const gear::TPCParameters  &pTPC      = Global::GEAR->getTPCParameters();
  const gear::PadRowLayout2D &padLayout = pTPC.getPadLayout();
  const gear::DoubleVec      &planeExt  = padLayout.getPlaneExtent();
  double TPC_z           = pTPC.getMaxDriftLength() ;
  double TPC_out         = planeExt[1]; 

  
  double add =  300.;          // Length of the hadron "tracks"
  double tz = TPC_z   + 350.;  // to start close to ECAL surface
  double to = TPC_out + 170.;  // to start close to ECAL surface
  
  double x1,x2,y1,y2,z1,z2;
  double px,py,pz,pt,pr, cx,cy,cz, rsq,d;
  
  int idpdg;
  const double* mom;
  float enr;
  LCCollection* mcpCol = evt->getCollection("MCParticle" ) ;  
  for(int i=0; i<mcpCol->getNumberOfElements() ; i++){
    MCParticle* imc = dynamic_cast<MCParticle*> ( mcpCol->getElementAt( i ) ) ;
    idpdg = imc-> getPDG (); 
    mom = imc-> getMomentum (); 
    enr = imc-> getEnergy (); 
    if( imc-> getGeneratorStatus() == 1 ) { // stable particles only   
      if(   // long lived neutral hadrons
	 (My_iabs(idpdg)==2112)|| // neutron
	 (My_iabs(idpdg)== 130)   // KoL
	 ) {
	if (enr > ecut) {//  neut. hadr. with E > 20 MeV
	  add = 50. + 300.*log(1.+enr);
	  /*
	     std::cout <<" Neut. hadron # "<< idpdg <<'\t'<<" Energy = "<<enr<<'\t'
	     <<",   Mom = "<<mom[0]<<",  "<<mom[1]<<",  "<<mom[2]
	     << std::endl;
	  */
	  px = mom[0]; 
	  py = mom[1]; 
	  pz = mom[2];
	  pt = hypot(px,py);	 
	  pr = hypot(pt,pz);	 
	  //  pr = sqrt(px*px + py*py + pz*pz);	 
	  // Direct Cosines
	  cx =  px/pr;
	  cy =  py/pr;
	  cz =  pz/pr;
	  //   Full length
	  d = tz/My_abs(cz);
	  //     2-D radius 
	  rsq = d*hypot(cx,cy);
	  if (rsq > to)       // start point at the cylinder of TPC
	    d = to/hypot(cx,cy);
	  x1 = d*cx;
	  y1 = d*cy;
	  z1 = d*cz;
	  x2 = x1 + add*cx;
	  y2 = y1 + add*cy;
	  z2 = z1 + add*cz;
	  //           yellow = 0xf9f920,   pink = 0xb900de
	  ced_line(x1,y1,z1,x2,y2,z2, NHEP_LAYER, 5, 0xf6a8f7);	 
	  x2 = x1;
	  y2 = y1;
	  z2 = z1;
	  x1 = 100.*cx;
	  y1 = 100.*cy;
	  z1 = 100.*cz;
	  ced_line(x1,y1,z1,x2,y2,z2, NHEP_LAYER , 1, 0xb900de);	 
	}
      }
      if(// short lived neutral hadrons
	 (My_iabs(idpdg)== 310)|| // KoS
	 (My_iabs(idpdg)==3122)|| // Lambda
	 (My_iabs(idpdg)==3212)|| // Sigma
	 (My_iabs(idpdg)==3322)   // Xi
	 ) {                      
	if (enr > ecut) {//  neut. hadr. with E > 20 MeV
	  add = 50. + 300.*log(1.+enr);
	  px = mom[0]; 
	  py = mom[1]; 
	  pz = mom[2];
	  pt = hypot(px,py);	 
	  pr = hypot(pt,pz);	 
	  // Direct Cosines
	  cx =  px/pr;
	  cy =  py/pr;
	  cz =  pz/pr;
	  x1 = 100.*cx;
	  y1 = 100.*cy;
	  z1 = 100.*cz;
	  x2 = x1 + add*cx;
	  y2 = y1 + add*cy;
	  z2 = z1 + add*cz;
	  ced_line(x1,y1,z1,x2,y2,z2,  NHEP_LAYER, 5, 0xb900de);
	}
      }
    }
  }
//fg   ced_send_event();    // Does not clean up screen
  
  return;
}   //  end draw_MC_Neutral_Hadrons
  
//-----------------------------------------------------------------------
void draw_MC_Charged_Hadrons(LCEvent * evt, float ecut){
//-----------------------------------------------------------------------
  unsigned kcol;
  int lwid;
  double x1,x2,y1,y2,z1,z2;
  double px,py,pz,pt,pr, cx,cy,cz;
  float rsq;
  float rad3;

  int idpdg;
  const double* mom;
  float enr;

  float EC = -2.9979251E-4;
  float BFIELD = 40.0; // !!!!!!!!!!!!!! ATTENTION !!!!!!!   MAGNETIC FIELD
  float charge;
  float step;
  float tet;
  float sint;
  float sintt;
  float tsint;
  float cos1t;
  float f1;
  float f2;
  float f3;
  int   npt1 = 0;
  float xp = 0.0;
  float yp = 0.0;
  float zp = 0.0;

  LCCollection* mcpCol = evt->getCollection("MCParticle" ) ;  
  for(int i=0; i<mcpCol->getNumberOfElements() ; i++){
    MCParticle* imc = dynamic_cast<MCParticle*> ( mcpCol->getElementAt( i ) ) ;
    idpdg = imc-> getPDG (); 
    mom   = imc-> getMomentum (); 
    enr   = imc-> getEnergy (); 
    //    float      charg = imc-> getCharge (); // this is equal zero for any particle
    if( imc-> getGeneratorStatus() == 1 ) { // stable particles only   
      if(!((My_iabs(idpdg)==2112)||(My_iabs(idpdg)== 311)||
	   (My_iabs(idpdg)== 130)||(My_iabs(idpdg)== 310)||  // neutral hadrons
	   (My_iabs(idpdg)==3122)||(My_iabs(idpdg)==3212))
	 && !((idpdg == 22)||(idpdg == 111))  // photons and Pi0s
	 && !((My_iabs(idpdg)==12)||(My_iabs(idpdg)==14)||(My_iabs(idpdg)==16)) // neutrino
	 ) { 
	//  Charge is not easy to calculate correctly from HEP record
	charge = 1.0;
	if(idpdg < 0) charge = -1.0;
	if((My_iabs(idpdg)==13)||(My_iabs(idpdg)==11)) charge = -charge;
	if (enr > ecut) {//  neut. hadr. with E > 20 MeV
	  px = mom[0]; 
	  py = mom[1]; 
	  pz = mom[2];
	  pt = hypot(px,py);
	  pr = hypot(pt,pz);
	  if(pt > 0.2) {
	    // Direct Cosines
	    cx =  px/pr;
	    cy =  py/pr;
	    cz =  pz/pr;
	    //     SUBROUTINE helix_geant_draw (charge,xp,yp,zp,cx,cy,cz,p,icol)
	    // ------------------------------------------------------------------
	    //       units are kgauss,centimeters,gev/c
	    // ------------------------------------------------------------------
	    int npt = 60;
	    float st = 500.0/npt;
	    if (pr < 0.7) {
	      npt = 300;
	      st = 500.0/npt;
	    }
	    npt1 = 0;
	    xp = 0.0;
	    yp = 0.0;
	    zp = 0.0;
	    x1 = 0.0;
	    y1 = 0.0;
	    z1 = 0.0;
	    for (int i2 = 0; i2 < npt+1; i2++) {
	      step = st*i2;                     // helix length
	      tet = charge*BFIELD*EC*step/pr;
	      sint = sin(tet);
	      sintt = (sint/tet);
	      tsint = (tet-sint)/tet;
	      cos1t = 2.*(sin(0.5*tet))*(sin(0.5*tet))/tet;
	      f1 = step * sintt;
	      f2 = step * cos1t;
	      f3 = step * tsint * cz;
	      x2 = 10.*(xp + (f1*cx - f2*cy)); // factor 10 because [mm] 
	      y2 = 10.*(yp + (f1*cy + f2*cx)); // factor 10 because [mm] 
	      z2 = 10.*(zp + (f1*cz + f3));    // factor 10 because [mm] 
// end     SUBROUTINE helix_geant_draw (charge,xp,yp,zp,cx,cy,cz,p,icol)
	      rsq  = hypot(x2,y2);
	      rad3 = hypot(rsq,z2);
	      npt1++ ;
	      if((npt1>0)&&(rad3>300.)&&(My_abs(z2)<6000.)&&(rsq<3000.)) { //do not draw outside volume
		kcol = 0x7af774;
		lwid = 1;
		if(My_iabs(idpdg)==13){ kcol = 0x0be228; lwid = 2;}  // muon
		if(My_iabs(idpdg)==11){ kcol = 0xe851a1; lwid = 2;}  // electron
		ced_line(x1,y1,z1,x2,y2,z2, CHEP_LAYER , lwid, kcol);	 
	      }
	      x1 = x2;
	      y1 = y2;
	      z1 = z2;
	    }
	  }
	}
      }
    }
  }
//fg   ced_send_event();    // Does not clean up screen
  
  return;
}  //  end draw_MC_Charged_Hadrons

//-----------------------------------------------------------------------
void MC_Particles_Balance(LCEvent * evt){
//-----------------------------------------------------------------------
  int idpdg;
  const double* mom;
  float enr;
  //double mass;
   LCCollection* mcpCol = evt->getCollection("MCParticle" ) ;  
//-----------------------------------------------------------------------
// Calculate balance at IP taking into account everything
//-----------------------------------------------------------------------
   double px,py,pz,pt,ttet;

   double e_to_tube  = 0.;
   double e_to_tubex = 0.;
   double e_to_tubey = 0.;
   double e_to_tubez = 0.;
   int    n_to_tube = 0; 

   double e_neutr = 0.;   
   double e_neutrx= 0.;   
   double e_neutry= 0.;   
   double e_neutrz= 0.;   
   int    n_neutr= 0;    
   
   double e_muon = 0.;  
   double e_muonx= 0.;  
   double e_muony= 0.;  
   double e_muonz= 0.;  
   int    n_muon= 0;   
   
   double e_elect = 0.; 
   double e_electx= 0.; 
   double e_electy= 0.; 
   double e_electz= 0.; 
   int    n_elect= 0;  
   
   double e_photon = 0.;
   double e_photonx= 0.;
   double e_photony= 0.;
   double e_photonz= 0.;
   int    n_photon= 0; 
   
   double e_pi0 = 0.;
   double e_pi0x= 0.;
   double e_pi0y= 0.;
   double e_pi0z= 0.;
   int    n_pi0= 0; 
   
   double e_llhadr = 0.;
   double e_llhadrx= 0.;
   double e_llhadry= 0.;
   double e_llhadrz= 0.;
   int    n_llhadr= 0; 
   
   double e_slhadr = 0.;
   double e_slhadrx= 0.;
   double e_slhadry= 0.;
   double e_slhadrz= 0.;
   int    n_slhadr= 0; 
   
   double e_chadr = 0.; 
   double e_chadrx= 0.; 
   double e_chadry= 0.; 
   double e_chadrz= 0.; 
   int    n_chadr= 0;  

   for(int i=0; i<mcpCol->getNumberOfElements() ; i++){
     MCParticle* imc = dynamic_cast<MCParticle*> ( mcpCol->getElementAt( i ) ) ;
     idpdg = imc-> getPDG (); 
     mom = imc-> getMomentum (); 
     enr = imc-> getEnergy (); 
     //mass = imc-> getMass (); 
     //    float cch = imc-> getCharge (); 
     if( imc-> getGeneratorStatus() == 1 ) { // stable particles only   
       px = mom[0]; 
       py = mom[1]; 
       pz = mom[2];
       pt = hypot(px,py);
       ttet = atan2(pt,pz);
       if ((fabs(ttet) < 0.1) || (fabs(M_PI-ttet) < 0.1)) {
	 e_to_tube  += enr;
	 e_to_tubex += px;
	 e_to_tubey += py;
	 e_to_tubez += pz;
	 n_to_tube ++;
       } 
       if((My_iabs(idpdg)==12)||(My_iabs(idpdg)==14)||(My_iabs(idpdg)==16)) {
	 e_neutr  += enr;
	 e_neutrx += px;
	 e_neutry += py;
	 e_neutrz += pz;
	 n_neutr ++;
       } 
       if(My_iabs(idpdg)==13) { // mu+ mu- 
	 e_muon  += enr;
	 e_muonx += px;
	 e_muony += py;
	 e_muonz += pz;
	 n_muon ++;
       } 
       if(My_iabs(idpdg)==11) { //  e+ e-
	 e_elect  += enr;
	 e_electx += px;
	 e_electy += py;
	 e_electz += pz;
	 n_elect ++;
       } 
       if(idpdg == 111) { // Pi0 as stable 
	 e_pi0  += enr;
	 e_pi0x += px;
	 e_pi0y += py;
	 e_pi0z += pz;
	 n_pi0 ++;
       } 
       if(idpdg == 22) { // photon
	 e_photon  += enr;
	 e_photonx += px;
	 e_photony += py;
	 e_photonz += pz;
	 n_photon ++;
       } 
       if(    // long lived neutral hadrons
	  (My_iabs(idpdg)==2112)|| // neutron
	  (My_iabs(idpdg)== 130)   // KoL
	  ) {                     
	 e_llhadr  += enr;
	 e_llhadrx += px;
	 e_llhadry += py;
	 e_llhadrz += pz;
	 n_llhadr ++;
       }
       if(  // short lived neutral hadrons
	  (My_iabs(idpdg)== 310)|| // KoS
	  (My_iabs(idpdg)==3122)|| // Lambda
	  (My_iabs(idpdg)==3212)|| // Sigma
	  (My_iabs(idpdg)==3322)   // Xi
	  ) {
	 e_slhadr  += enr;
	 e_slhadrx += px;
	 e_slhadry += py;
	 e_slhadrz += pz;
	 n_slhadr ++;
       }
       if(!(My_iabs(idpdg)==12) && !(My_iabs(idpdg)==14) && !(My_iabs(idpdg)==16) &&  // neutrinos   
	  !(My_iabs(idpdg)==13) && // mu+ mu- 
	  !(My_iabs(idpdg)==11) && //  e+ e-
	  !(idpdg == 111) && // Pi0
	  !(idpdg == 22) &&  // photon
	  !(My_iabs(idpdg)==2112) && !(My_iabs(idpdg)== 311) && // neutral hadrons
	  !(My_iabs(idpdg)== 130) && !(My_iabs(idpdg)== 310) && // neutral hadrons
	  !(My_iabs(idpdg)==3122) && !(My_iabs(idpdg)==3212)) {  // charged hadrons
	 e_chadr  += enr;
	 e_chadrx += px;
	 e_chadry += py;
	 e_chadrz += pz;
	 n_chadr ++;
       }
     }
   }
  
  double e_sum = e_elect+e_muon+e_chadr+e_pi0+e_photon+e_llhadr+e_slhadr+e_neutr;
  double evt_energy = e_sum;
  //   Minus muon loosed energy = 1.3 GeV in average
  double e_mu_lost = e_muon  - n_muon*1.3;
  double e_lost = e_neutr + e_mu_lost + e_to_tube;
  double e_real = evt_energy - e_lost;
  if(e_real< 0.0) e_real = 0.000001;
  
  std::cout <<" =============================================================="<< std::endl;
  std::cout << " ======== HEP Record Balance for event: " << evt->getEventNumber()<<" ======="<< std::endl ;
  std::cout <<" =============================================================="<< std::endl;
  std::cout <<" ==============  Possible lost  ==================="<< std::endl;
  std::cout <<"  Neutrino energy      = "<<e_neutr<<",  in "<< n_neutr<<" neutrinos"<< std::endl;
  std::cout <<"  Energy to beam tube  = "<<e_to_tube<<",  in "<<n_to_tube<<" particles"<< std::endl;
  std::cout <<"  Muons energy lost    = "<<e_mu_lost<<"  in "<<n_muon<<" muons"<< std::endl;
  std::cout <<"  --------------------------------------------------"<< std::endl;
  std::cout <<"  Total Event energy at IP = "<<evt_energy<<" [GeV]"<< std::endl;
  std::cout <<"  --------------------------------------------------"<< std::endl;
  std::cout <<"  Whole lost Energy        = "<<e_lost<< std::endl;
  std::cout <<"  Available Energy in calo.= "<<e_real<< std::endl;
  std::cout <<" =============================================================="<< std::endl;
  std::cout <<"  Muon energy               = "<<e_muon  <<'\t'<<"  in "<< n_muon  <<" muons"<< std::endl;
  std::cout <<"  Electron energy           = "<<e_elect <<'\t'<<"  in "<< n_elect <<" electrons"<< std::endl;
  std::cout <<"  Charged hadron energy     = "<<e_chadr <<'\t'<<"  in "<< n_chadr <<" hadrons"<< std::endl;
  std::cout <<"  -------------------------------------------------------------"<< std::endl;
  std::cout <<"  Pi0 energy (if stable)    = "<<e_pi0   <<'\t'<<"  in "<< n_pi0   <<" Pi zeros"<< std::endl;
  std::cout <<"  Photon energy             = "<<e_photon<<'\t'<<"  in "<< n_photon<<" photons"<< std::endl;
  std::cout <<"  -------------------------------------------------------------"<< std::endl;
  std::cout <<"  Long lived hadron energy  = "<<e_llhadr<<'\t'<<"  in "<< n_llhadr<<" hadrons"<< std::endl;
  std::cout <<"  Short lived hadron energy = "<<e_slhadr<<'\t'<<"  in "<< n_slhadr<<" hadrons"<< std::endl;
  std::cout <<" =============================================================="<< std::endl;
}

//-----------------------------------------------------------------------
void print_MC_Particles(LCEvent * evt, float ecut){
  //-----------------------------------------------------------------------
  int idpdg;
  const double* mom;
  float enr;
  //float mass;
  LCCollection* mcpCol = evt->getCollection("MCParticle" ) ;  
//-----------------------------------------------------------------------
//            Print  Pi zeros
//-----------------------------------------------------------------------
/*
  for(int i=0; i<mcpCol->getNumberOfElements() ; i++){
    MCParticle* imc = dynamic_cast<MCParticle*> ( mcpCol->getElementAt( i ) ) ;
    idpdg = imc-> getPDG (); 
    mom = imc-> getMomentum (); 
    enr = imc-> getEnergy (); 
    mass = imc-> getMass (); 
    if(idpdg == 111) { //  Pi_0s    
      std::cout <<" Pi_0 # "<< idpdg <<'\t'<<" Energy = "<<enr<<'\t'
		<<",   Mom = "<<mom[0]<<",  "<<'\t'<<mom[1]<<",  "<<'\t'<<mom[2]
		<< std::endl;
    }
  }
  std::cout <<"      All Pi_0s are decayed into gammas  !"<< std::endl<<std::endl;
*/
//-----------------------------------------------------------------------
//              Print Neutrinos
//-----------------------------------------------------------------------
   int n_nu = 0;
 for(int i=0; i<mcpCol->getNumberOfElements() ; i++){
    MCParticle* imc = dynamic_cast<MCParticle*> ( mcpCol->getElementAt( i ) ) ;
    idpdg = imc-> getPDG (); 
    mom = imc-> getMomentum (); 
    enr = imc-> getEnergy (); 
    //mass = imc-> getMass (); 
    if( imc-> getGeneratorStatus() == 1 ) { // stable particles only   
      if((My_iabs(idpdg)==12)||(My_iabs(idpdg)==14)||(My_iabs(idpdg)==16)) { // neutrinos   
	  n_nu ++;
	std::cout <<" Neutrino # "<< idpdg <<'\t'<<" Energy = "<<enr<<'\t'
		  <<",   Mom = "<<mom[0]<<",  "<<'\t'<<mom[1]<<",  "<<'\t'<<mom[2]
		  << std::endl;
      }
    }
  }
  std::cout <<"        Number of neutrinos, E > "<<ecut<< " GeV = "<< n_nu <<std::endl<<std::endl; 
//-----------------------------------------------------------------------
//              Print  mu+ mu- 
//-----------------------------------------------------------------------
  int n_mu = 0;
  for(int i=0; i<mcpCol->getNumberOfElements() ; i++){
    MCParticle* imc = dynamic_cast<MCParticle*> ( mcpCol->getElementAt( i ) ) ;
    idpdg = imc-> getPDG (); 
    mom = imc-> getMomentum (); 
    enr = imc-> getEnergy (); 
    //mass = imc-> getMass (); 
    if( imc-> getGeneratorStatus() == 1 ) { // stable particles only   
      if(My_iabs(idpdg)==13) { // mu+ mu- 
	  n_mu ++;
	std::cout <<" Muon # "<< idpdg <<'\t'<<" Energy = "<<enr<<'\t'
		  <<",   Mom = "<<mom[0]<<",  "<<'\t'<<mom[1]<<",  "<<'\t'<<mom[2]
		  << std::endl;
      }
    }
  }
  std::cout <<"        Number of muons, E > "<<ecut<< " GeV = "<< n_mu <<std::endl<<std::endl; 
//-----------------------------------------------------------------------
//              Print  e+ e- 
//-----------------------------------------------------------------------
  int n_elect = 0;
 for(int i=0; i<mcpCol->getNumberOfElements() ; i++){
    MCParticle* imc = dynamic_cast<MCParticle*> ( mcpCol->getElementAt( i ) ) ;
    idpdg = imc-> getPDG (); 
    mom = imc-> getMomentum (); 
    enr = imc-> getEnergy (); 
    //mass = imc-> getMass (); 
    if( imc-> getGeneratorStatus() == 1 ) { // stable particles only   
      if(My_iabs(idpdg)==11) { //  e+ e-
	  n_elect ++;
	std::cout <<" Electron # "<< idpdg <<'\t'<<" Energy = "<<enr<<'\t'
		  <<",   Mom = "<<mom[0]<<",  "<<'\t'<<mom[1]<<",  "<<'\t'<<mom[2]
		  << std::endl;
      }
    }
  }
  std::cout <<"        Number of electrons, E > "<<ecut<< " GeV = "<< n_elect <<std::endl<<std::endl; 
//-----------------------------------------------------------------------
//              Print photons
//-----------------------------------------------------------------------
  int n_photons = 0;
  for(int i=0; i<mcpCol->getNumberOfElements() ; i++){
    MCParticle* imc = dynamic_cast<MCParticle*> ( mcpCol->getElementAt( i ) ) ;
    idpdg = imc-> getPDG (); 
    mom = imc-> getMomentum (); 
    enr = imc-> getEnergy (); 
    //mass = imc-> getMass (); 
    if( imc-> getGeneratorStatus() == 1 ) { // stable particles only   
      if(idpdg == 22) { // photon
	if (enr > ecut) {//  All gammas with E > 20 MeV
	  n_photons ++;
	  std::cout <<" Photon # "<< idpdg <<'\t'<<" Energy = "<<enr<<'\t'
		    <<",   Mom = "<<mom[0]<<",  "<<'\t'<<mom[1]<<",  "<<'\t'<<mom[2]
		    << std::endl;
	}
      }
    }
  }
  std::cout <<"       Number of photons, E > "<<ecut<< " GeV = "<< n_photons <<std::endl<<std::endl; 
//-----------------------------------------------------------------------
//              Print neutral hadrons
//-----------------------------------------------------------------------
  int n_hadr = 0;
  for(int i=0; i<mcpCol->getNumberOfElements() ; i++){
    MCParticle* imc = dynamic_cast<MCParticle*> ( mcpCol->getElementAt( i ) ) ;
    idpdg = imc-> getPDG (); 
    mom = imc-> getMomentum (); 
    enr = imc-> getEnergy (); 
    //mass = imc-> getMass (); 
    if( imc-> getGeneratorStatus() == 1 ) { // stable particles only   
      if((My_iabs(idpdg)==2112)||(My_iabs(idpdg)== 311)||
	 (My_iabs(idpdg)== 130)||(My_iabs(idpdg)== 310)||
	 (My_iabs(idpdg)==3122)||(My_iabs(idpdg)==3212)) { // neutron hadrons
	if (enr > ecut) {// neutral hadrons with E > 20 MeV
	  n_hadr ++;
	  std::cout <<" Neutral hadron # "<< idpdg <<'\t'<<" Energy = "<<enr<<'\t'
		    <<",   Mom = "<<mom[0]<<",  "<<'\t'<<mom[1]<<",  "<<'\t'<<mom[2]
		    << std::endl;
	}
      }
    }
  }
  std::cout <<"        Number of neutral hadrons, E > "<<ecut<< " GeV = "<< n_hadr <<std::endl<<std::endl; 
//-----------------------------------------------------------------------
//              Print charged hadrons
//-----------------------------------------------------------------------
  int c_hadr = 0;
  for(int i=0; i<mcpCol->getNumberOfElements() ; i++){
    MCParticle* imc = dynamic_cast<MCParticle*> ( mcpCol->getElementAt( i ) ) ;
    idpdg = imc-> getPDG (); 
    mom = imc-> getMomentum (); 
    enr = imc-> getEnergy (); 
    //mass = imc-> getMass (); 
    if( imc-> getGeneratorStatus() == 1 ) { // stable particles only   
      if(!(My_iabs(idpdg)==12) && !(My_iabs(idpdg)==14) && !(My_iabs(idpdg)==16) &&  // neutrinos   
	 !(My_iabs(idpdg)==13) && // mu+ mu- 
	 !(My_iabs(idpdg)==11) && //  e+ e-
	 !(idpdg == 22) && // photon
	 !(My_iabs(idpdg)==2112) && !(My_iabs(idpdg)== 311) && // neutral hadrons
	 !(My_iabs(idpdg)== 130) && !(My_iabs(idpdg)== 310) && // neutral hadrons
	 !(My_iabs(idpdg)==3122) && !(My_iabs(idpdg)==3212)) { // neutral hadrons
	if (enr > ecut) {// charged hadrons with E > 20 MeV
	  c_hadr ++;
	  std::cout <<" Charged hadron # "<< idpdg <<'\t'<<" Energy = "<<enr<<'\t'
		    <<",   Mom = "<<mom[0]<<",  "<<'\t'<<mom[1]<<",  "<<'\t'<<mom[2]
		    << std::endl;
	}
      }
    }
  }
  std::cout <<"        Number of charged hadrons, E > "<<ecut<< " GeV = "<< c_hadr <<std::endl<<std::endl; 
  
  return;
}   //  end print_MC_Particles

}// namespace marlin
