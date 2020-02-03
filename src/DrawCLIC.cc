#include "DrawCLIC.h"
#include "MarlinCED.h"
#include <iostream>

#include <gear/CalorimeterParameters.h>
#include <gear/LayerLayout.h>
//#include <gear/SiPlanesParameters.h>
//#include <gear/SiPlanesLayerLayout.h>
#include <gear/VXDLayerLayout.h>
#include <gear/VXDParameters.h>
#include <gear/FTDParameters.h>
#include <gear/FTDLayerLayout.h>
#include <gear/ZPlanarParameters.h>
#include <gear/ZPlanarLayerLayout.h>


#include <signal.h>
#include <unistd.h>
#include <cstdlib>
#include <cstdio>
#include <ctime>
#include <time.h>
#include <termios.h>
#include <poll.h>
#include <sys/select.h>
#include <termios.h>


// ----- include for verbosity dependend logging ---------
#include "marlin/VerbosityLevels.h"

using namespace lcio ;
using namespace marlin ;


DrawCLIC aDrawCLIC ;


DrawCLIC::DrawCLIC() : Processor("DrawCLIC") {
    
    // modify processor description
    _description = "DrawCLIC draws the CLIC detector in CED" ;
    
    
    registerProcessorParameter( "DrawDetector" ,
                               "Call the drawDetector function else call draw function",
                               _begin ,
                               (bool)0  ) ;
    
}



void DrawCLIC::init() {
    
    streamlog_out(DEBUG) << "   init called  " << std::endl ;
    
    // usually a good idea to
    printParameters() ;
    
    MarlinCED::init(this);
    
    _nRun = 0 ;
    _nEvt = 0 ;
    
}


void DrawCLIC::processRunHeader( LCRunHeader* /*run*/) {
    
    _nRun++ ;
}



void DrawCLIC::processEvent( LCEvent * /*evt*/ ) {
    
    MarlinCED::newEvent(this , -1 );
    drawGEARDetector();
    MarlinCED::draw(this, 1 );
    _nEvt ++ ;
}



void DrawCLIC::check( LCEvent * /*evt*/ ) {
    // nothing to check here - could be used to fill checkplots in reconstruction processor
}


void DrawCLIC::end(){
    
    //   std::cout << "DrawCLIC::end()  " << name()
    // 	    << " processed " << _nEvt << " events in " << _nRun << " runs "
    // 	    << std::endl ;
    
}


void DrawCLIC::drawGEARDetector(){
    
    drawDetectorFromGear( Global::GEAR ) ;
}



void DrawCLIC::drawDetectorFromGear( gear::GearMgr* gearMgr ){
    //
    // based on original code from V.Morgunov, MPI
    //
    
    //############ TPC #########################
  bool showTPC=false ;//true;
  
  float r_min_tpc = 0;
  float r_max_tpc = 0;
  float z_max_tpc = 0;
  
  // try{
  //       const gear::TPCParameters&  pTPC      = gearMgr->getTPCParameters();
        
  //       // Multi-module support
  //       const gear::DoubleVec&      planeExt  = pTPC.getPlaneExtent();
  //       r_min_tpc = planeExt[0];
  //       r_max_tpc = planeExt[1];
  //       z_max_tpc = pTPC.getMaxDriftLength();
  //   }catch(gear::UnknownParameterException& e){
  //       showTPC=false;
  //   }
    
    
    // ########## F #####################
    bool showECAL = true;
    bool showECALEndcap = true;
    float r_min_ecal_bar = 0;
    float r_max_ecal_bar = 0;
    //float z_min_ecal_bar = 0;
    float z_max_ecal_bar = 0;
    float r_max_ecal_ecap = 0;
    float z_min_ecal_ecap = 0;
    float z_max_ecal_ecap = 0;
    float r_min_ecal_ecap = 0;
    
    try{
        const gear::CalorimeterParameters& pECAL_B =
        gearMgr->getEcalBarrelParameters();
        r_min_ecal_bar = pECAL_B.getExtent()[0];
        r_max_ecal_bar = pECAL_B.getExtent()[1];
        // float z_min_ecal_bar = pECAL_B.getExtent()[2];
        z_max_ecal_bar = pECAL_B.getExtent()[3];
        const gear::CalorimeterParameters& pECAL_E =
        gearMgr->getEcalEndcapParameters();
        r_min_ecal_ecap = pECAL_E.getExtent()[0];
        r_max_ecal_ecap = pECAL_E.getExtent()[1];
        z_min_ecal_ecap = pECAL_E.getExtent()[2];
        z_max_ecal_ecap = pECAL_E.getExtent()[3];
    }catch(gear::UnknownParameterException& e){
        showECAL=false;
        showECALEndcap=false;
    }
    
    
    
    
    
    //############# HCAL ##########################
    bool showHCAL=true;
    bool showHCALRing=true;
    bool showHCALEndcap=true;
    float r_min_hcal_bar =0;
    float r_max_hcal_bar =0;
    //    float z_min_hcal_bar =0;
    float z_max_hcal_bar =0;
    float r_min_hcal_ring=0;
    float r_max_hcal_ring=0;
    float z_min_hcal_ring=0;
    float z_max_hcal_ring=0;
    float r_min_hcal_ecap=0;
    float r_max_hcal_ecap=0;
    float z_min_hcal_ecap=0;
    float z_max_hcal_ecap=0;
    
    try{
        const gear::CalorimeterParameters& pHCAL_B =
        gearMgr->getHcalBarrelParameters();
        //  _innerHcalRadius = float(pHcalBarrel.getExtent()[0]);
        r_min_hcal_bar = pHCAL_B.getExtent()[0];
        r_max_hcal_bar = pHCAL_B.getExtent()[1];
        //float z_min_hcal_bar = pHCAL_B.getExtent()[2];
        z_max_hcal_bar = pHCAL_B.getExtent()[3];
        const gear::CalorimeterParameters& pHCAL_R =
        gearMgr->getHcalRingParameters();
        r_min_hcal_ring = pHCAL_R.getExtent()[0];
        r_max_hcal_ring = pHCAL_R.getExtent()[1];
        z_min_hcal_ring = pHCAL_R.getExtent()[2];
        z_max_hcal_ring = pHCAL_R.getExtent()[3];
        const gear::CalorimeterParameters& pHCAL_E =
        gearMgr->getHcalEndcapParameters();
        r_min_hcal_ecap = pHCAL_E.getExtent()[0];
        r_max_hcal_ecap = pHCAL_E.getExtent()[1];
        z_min_hcal_ecap = pHCAL_E.getExtent()[2];
        z_max_hcal_ecap = pHCAL_E.getExtent()[3];
    }catch(gear::UnknownParameterException& e){
        showHCAL=false;
        showHCALRing=false;
        showHCALEndcap=false;
    }
    
    
    
    float r_min_lhcal = 0.0 ;
    float r_max_lhcal = 0.0 ;
    float z_min_lhcal = 0.0 ;
    
    float z_max_lhcal = 0.0 ;
    
    bool showLHcal = false ;
    // make this optional = as CLIC does not have an LHcal
    try{
        
        const gear::CalorimeterParameters& pLHCal =
        gearMgr->getLHcalParameters();
        
        r_min_lhcal = pLHCal.getExtent()[0];
        r_max_lhcal = pLHCal.getExtent()[1];
        z_min_lhcal = pLHCal.getExtent()[2];
        z_max_lhcal = pLHCal.getExtent()[3];
        
        showLHcal = true ;
    }
    catch( gear::UnknownParameterException& e){
    }
    
    bool showLCal = false ;
    
    
    float r_min_lcal=0.0;
    float r_max_lcal=0.0;
    float z_min_lcal=0.0;
    float z_max_lcal=0.0;
    
    try{
        const gear::CalorimeterParameters& pLCal =
        gearMgr->getLcalParameters();
        r_min_lcal = pLCal.getExtent()[0];
        r_max_lcal = pLCal.getExtent()[1];
        z_min_lcal = pLCal.getExtent()[2];
        z_max_lcal = pLCal.getExtent()[3];
        showLCal = true ;
    }
    catch( gear::UnknownParameterException& e){
    }
    
    
    //######## Beamcal ##############################
    bool showBeamcal=true;
    float r_min_beamcal=0;
    float r_max_beamcal=0;
    float z_min_beamcal=0;
    float z_max_beamcal=0;
    
    try{
        const gear::CalorimeterParameters& pBeamcal =
        gearMgr->getBeamCalParameters();
        r_min_beamcal = pBeamcal.getExtent()[0];
        r_max_beamcal = pBeamcal.getExtent()[1];
        z_min_beamcal = pBeamcal.getExtent()[2];
        z_max_beamcal = pBeamcal.getExtent()[3];
    }catch( gear::UnknownParameterException& e){
        showBeamcal=false;
    }
    
    
    //############## Yoke ############################
    bool showYoke = true;
    
    float r_min_yoke_bar=0;
    float r_max_yoke_bar=0;
    float z_max_yoke_bar=0;
    float r_min_yoke_plug=0;
    float r_max_yoke_plug=0;
    float z_min_yoke_plug=0;
    float z_max_yoke_plug=0;
    float r_min_yoke_ecap=0;
    float r_max_yoke_ecap=0;
    float z_min_yoke_ecap=0;
    float z_max_yoke_ecap=0;
    
    try{
        const gear::CalorimeterParameters& pYOKE_B =
        gearMgr->getYokeBarrelParameters();
        //  _innerYokeRadius = float(pYokeBarrel.getExtent()[0]);
        r_min_yoke_bar = pYOKE_B.getExtent()[0];
        r_max_yoke_bar = pYOKE_B.getExtent()[1];
        //float z_min_yoke_bar = pYOKE_B.getExtent()[2];
        z_max_yoke_bar = pYOKE_B.getExtent()[3];
        const gear::CalorimeterParameters& pYOKE_R =
        gearMgr->getYokePlugParameters();
        r_min_yoke_plug = pYOKE_R.getExtent()[0];
        r_max_yoke_plug = pYOKE_R.getExtent()[1];
        z_min_yoke_plug = pYOKE_R.getExtent()[2];
        z_max_yoke_plug = pYOKE_R.getExtent()[3];
        const gear::CalorimeterParameters& pYOKE_E =
        gearMgr->getYokeEndcapParameters();
        r_min_yoke_ecap = pYOKE_E.getExtent()[0];
        r_max_yoke_ecap = pYOKE_E.getExtent()[1];
        z_min_yoke_ecap = pYOKE_E.getExtent()[2];
        z_max_yoke_ecap = pYOKE_E.getExtent()[3];
    }catch( gear::UnknownParameterException& e){
        showYoke=false;
    }
    
    
    
    
    // ------- coil parameters have changed in ILD_01
    bool showCoil = true;
    
    float coil_half_z        =  0 ;
    float coil_inner_radius  =  0 ;
    float coil_outer_radius  =  0 ;
    
    
    
    try{
        const gear::GearParameters&  pCoil      = gearMgr->getGearParameters("CoilParameters");
        try {
            
            coil_half_z         =  pCoil.getDoubleVal("Coil_cryostat_half_z" ) ;
            coil_inner_radius  =   pCoil.getDoubleVal("Coil_cryostat_inner_radius" ) ;
            coil_outer_radius  =   pCoil.getDoubleVal("Coil_cryostat_outer_radius" ) ;
            
        }  catch( gear::UnknownParameterException& e){
            // the parameters named _inner_cyl_ seem to be the ones that define the envelope (strangely enough)....
            coil_half_z         =  pCoil.getDoubleVal("Coil_cryostat_inner_cyl_half_z" ) ;
            coil_inner_radius  =   pCoil.getDoubleVal("Coil_cryostat_inner_cyl_inner_radius" ) ;
            coil_outer_radius  =   pCoil.getDoubleVal("Coil_cryostat_inner_cyl_outer_radius" ) ;
        }
    }catch( gear::UnknownParameterException& e){
        showCoil=false;
    }
    
    
    //============================================================================================================
    // here we might either have default GearParameters for the FTD or
    // starting from ILD_01 proper FTDParameters ...
    // we fill the layer (disk) parameters into four DoubleVecs that are then used for drawing
    // the detctor:
    
    DoubleVec ftd_d  ;  // thickness
    DoubleVec ftd_ri ;  // inner r
    DoubleVec ftd_ro ;  // outer r
    DoubleVec ftd_z  ;  // z position
    
    try{
        
        //     const gear::FTDParameters&  pFTD = gearMgr->getFTDParameters();
        const gear::FTDLayerLayout&  pFTD = gearMgr->getFTDParameters().getFTDLayerLayout()  ;
        
        streamlog_out( DEBUG2 ) << " filling FTD parameters from gear::FTDParameters - n layers : " <<  pFTD.getNLayers() << std::endl ;
        
        for( unsigned i=0, N = pFTD.getNLayers() ; i<N ; ++i ){
            
            // this only really works for the staggered design
            //create the even numbered petall
            if( pFTD.getAlpha(i) != 0  ) {
                streamlog_out( ERROR ) << "drawCLICdetector: Cannot draw design for tilt angle (alpha) != 0.0 " << pFTD.getAlpha(i)  << " exit(1) called" << std::endl ;
                exit(1);
            }
            
            int nsensors = pFTD.getNSensors(i) ;
            
            // create a disk to represent even number petals front side
            ftd_d .push_back( pFTD.getSensitiveThickness(i) ) ;
            ftd_ri.push_back( pFTD.getSensitiveRinner(i) ) ;
            ftd_ro.push_back( pFTD.getMaxRadius(i) ) ;
            ftd_z .push_back( pFTD.getSensitiveZposition(i, 0, 1) ) ;
            
            // create a disk to represent odd number petals front side
            ftd_d .push_back( pFTD.getSensitiveThickness(i) ) ;
            ftd_ri.push_back( pFTD.getSensitiveRinner(i) ) ;
            ftd_ro.push_back( pFTD.getMaxRadius(i) ) ;
            ftd_z .push_back( pFTD.getSensitiveZposition(i, 1, 1) ) ;
            
            
            if (pFTD.isDoubleSided(i)) {
                
                // create a disk to represent even number petals rear side
                ftd_d .push_back( pFTD.getSensitiveThickness(i) ) ;
                ftd_ri.push_back( pFTD.getSensitiveRinner(i) ) ;
                ftd_ro.push_back( pFTD.getMaxRadius(i) ) ;
                ftd_z .push_back( pFTD.getSensitiveZposition(i, 0, (nsensors/2))+1 ) ;
                
                
                // create a disk to represent odd number petals rear side
                ftd_d .push_back( pFTD.getSensitiveThickness(i) ) ;
                ftd_ri.push_back( pFTD.getSensitiveRinner(i) ) ;
                ftd_ro.push_back( pFTD.getMaxRadius(i) ) ;
                ftd_z .push_back( pFTD.getSensitiveZposition(i, 1, (nsensors/2))+1 ) ;
                
            }
            
            
            
        }
        
        
    }
    catch( gear::UnknownParameterException& e){}
    
    
    try{
        
        const gear::GearParameters& pFTD = gearMgr->getGearParameters("FTD");
        
        streamlog_out( DEBUG2 ) << " filling FTD parameters from old gear::GearParameters " << std::endl ;
        
        const DoubleVec& FTD_d   =  pFTD.getDoubleVals("FTDDiskSupportThickness" )  ;
        const DoubleVec& FTD_ri  =  pFTD.getDoubleVals("FTDInnerRadius" )  ;
        const DoubleVec& FTD_ro  =  pFTD.getDoubleVals("FTDOuterRadius" )  ;
        const DoubleVec& FTD_z   =  pFTD.getDoubleVals("FTDZCoordinate" )  ;
        
        
        std::copy( FTD_d.begin() , FTD_d.end() , std::back_inserter(  ftd_d  )  ) ;
        std::copy( FTD_ri.begin(), FTD_ri.end(), std::back_inserter(  ftd_ri )  ) ;
        std::copy( FTD_ro.begin(), FTD_ro.end(), std::back_inserter(  ftd_ro )  ) ;
        std::copy( FTD_z.begin() , FTD_z.end() , std::back_inserter(  ftd_z  )  ) ;
    }
    catch( gear::UnknownParameterException& e){
    }
    
    //note: if both try blocks fail, the ftd vectors simply will be empty and no disks will be drawn
    //============================================================================================================
    
    
    //-- VXD Parameters--
    int nLayersVTX = 0 ;
    const gear::VXDParameters* pVXDDetMain = 0;
    const gear::VXDLayerLayout* pVXDLayerLayout = 0;
    
    try{
        pVXDDetMain = &gearMgr->getVXDParameters();
        pVXDLayerLayout = &(pVXDDetMain->getVXDLayerLayout());
        nLayersVTX = pVXDLayerLayout->getNLayers();
    }
    catch( gear::UnknownParameterException& e){
    }
    
    
    //-- SET Parameters--
    int nLayersSET = 0 ;
    const gear::ZPlanarParameters* pSETDetMain = 0;
    const gear::ZPlanarLayerLayout* pSETLayerLayout =0;
    
    try{
        pSETDetMain = &gearMgr->getSETParameters();
        pSETLayerLayout = &(pSETDetMain->getZPlanarLayerLayout());
        nLayersSET = pSETLayerLayout->getNLayers();
    }
    catch( gear::UnknownParameterException& e){
    }
    
    //-- SIT Parameters--
    int nLayersSIT = 0 ;
    const gear::ZPlanarParameters* pSITDetMain = 0;
    const gear::ZPlanarLayerLayout* pSITLayerLayout = 0;
    
    try{
        pSITDetMain = &gearMgr->getSITParameters();
        pSITLayerLayout = &(pSITDetMain->getZPlanarLayerLayout());
        nLayersSIT = pSITLayerLayout->getNLayers();
    }
    catch( gear::UnknownParameterException& e){
    }
    
    
    
    float rad2deg = 180.0 / M_PI;
    DoubleVec rSIT ;
    DoubleVec lSIT ;
    
    //old SIT using cylinders
    try{
        const gear::GearParameters& pSITDet = gearMgr->getGearParameters("SIT");
        
        const DoubleVec& rSIT_temp  = pSITDet.getDoubleVals("SITLayerRadius")  ;
        const DoubleVec& lSIT_temp  = pSITDet.getDoubleVals("SITLayerHalfLength") ;
        // only in ILD_01
        //   const DoubleVec& thSIT = pSITDet.getDoubleVals("SITLayerThickness") ; // SITSupportLayerThickness ?
        std::copy( rSIT_temp.begin() , rSIT_temp.end() , std::back_inserter(  rSIT )  ) ;
        std::copy( lSIT_temp.begin() , lSIT_temp.end() , std::back_inserter(  lSIT )  ) ;
    }
    catch( gear::UnknownParameterException& e){
    }
    
    
    // ======= ================================================================
    //To convert inner radius of polygone to its outer radius
    //   float Cos4  = cos(M_PI/4.0);
    
    float Cos8  = cos(M_PI/8.0);
    float Cos12  = cos(M_PI/12.0);
    //float Cos16 = cos(M_PI/16.);
    // convertion of  inner radius of polygone to its outer radius
    float r_inn_ecal_bar     = r_min_ecal_bar/Cos12;
    float r_out_ecal_bar     = (r_max_ecal_bar)/Cos12 ;
    float r_inn_ecal_ecap    = r_min_ecal_ecap;///Cos12;
    float r_out_ecal_ecap    = r_max_ecal_ecap/Cos12;
    float thick_ecal_ecap    = 0.5*(z_max_ecal_ecap - z_min_ecal_ecap);
    float shift_ecal_z_plus  = z_min_ecal_ecap;
    float shift_ecal_z_minus = z_min_ecal_ecap + 2.0*thick_ecal_ecap;
    
    float r_inn_hcal_bar     = r_min_hcal_bar/Cos12;
    float r_out_hcal_bar     = r_max_hcal_bar/Cos12;
    
    float r_inn_hcal_ring    = r_min_hcal_ring/Cos12; //fg: was cos4
    float r_out_hcal_ring    = r_max_hcal_ring/Cos12;
    float thick_hcal_ring    = 0.5*(z_max_hcal_ring -
                                    z_min_hcal_ring + 20.0); // +20 by hand to see hits inside
    float shift_hcalr_z_plus  = z_min_hcal_ring;
    float shift_hcalr_z_minus = z_min_hcal_ring + 2.0*thick_hcal_ring;
    //   float r_inn_hcal_ecap    = r_min_hcal_ecap/Cos4;
    
    float r_out_hcal_ecap    = r_max_hcal_ecap/Cos12;

    //float r_out_hcal_ecap    = r_max_hcal_ecap  ; //fg: the encap driver writes out the outer radius...
    
    float thick_hcal_ecap    = 0.5*(z_max_hcal_ecap -
                                    z_min_hcal_ecap + 20.0); // +20 by hand to see hits inside
    float shift_hcal_z_plus  = z_min_hcal_ecap;
    float shift_hcal_z_minus = z_min_hcal_ecap + 2.0*thick_hcal_ecap;
    
    float thick_lhcal         = 0.5 * ( z_max_lhcal -  z_min_lhcal ) ;
    float shift_lhcal_z_plus  = z_min_lhcal;
    float shift_lhcal_z_minus = z_min_lhcal +  2.0 * thick_lhcal ;
    
    float thick_lcal         = 0.5 * ( z_max_lcal -  z_min_lcal ) ;
    float shift_lcal_z_plus  = z_min_lcal;
    float shift_lcal_z_minus = z_min_lcal +  2.0 * thick_lcal ;
    
    float thick_beamcal         = 0.5 * ( z_max_beamcal -  z_min_beamcal ) ;
    float shift_beamcal_z_plus  = z_min_beamcal;
    float shift_beamcal_z_minus = z_min_beamcal +  2.0 * thick_beamcal ;
    
    
    float r_inn_yoke_bar     = r_min_yoke_bar/Cos8;
    float r_out_yoke_bar     = r_max_yoke_bar/Cos8;
    
    float r_inn_yoke_plug    = r_min_yoke_plug/Cos8; //fg: was cos4
    float r_out_yoke_plug    = r_max_yoke_plug/Cos12;
    float thick_yoke_plug    = 0.5*(z_max_yoke_plug -
                                    z_min_yoke_plug + 20.0); // +20 by hand to see hits inside
    float shift_yoker_z_plus  = z_min_yoke_plug;
    float shift_yoker_z_minus = z_min_yoke_plug + 2.0*thick_yoke_plug;
    float r_inn_yoke_ecap    = r_min_yoke_ecap/Cos8;
    float r_out_yoke_ecap    = r_max_yoke_ecap/Cos8;
    float thick_yoke_ecap    = 0.5*(z_max_yoke_ecap -
                                    z_min_yoke_ecap + 20.0); // +20 by hand to see hits inside
    float shift_yoke_z_plus  = z_min_yoke_ecap;
    float shift_yoke_z_minus = z_min_yoke_ecap + 2.0*thick_yoke_ecap;
    
    // ========================================================================
    
    
    // colors used in Mokka ILD_00
    static const unsigned vxdCol  = 0xafafaf ;
    static const unsigned sitCol  = 0x8c6e0a ; // light grey
    static const unsigned setCol  = 0xdddddd ;
    static const unsigned tpcCol  = 0xf5f300 ;
    static const unsigned ecalCol = 0x339900 ;
    static const unsigned hcalCol = 0x141e96 ;
    static const unsigned yokeCol = 0x990000 ;
    static const unsigned coilCol = 0x333333 ;
    static const unsigned ftdCol  = 0x8c6e0a ;
    static const unsigned fcalCol = 0x5a0078 ;
    static const unsigned bcalCol = 0x5a0078 ;
    static const unsigned lcalCol = 0x006666 ;

    // define layers for sub detectors
    static const int fDL = NUMBER_DATA_LAYER; //  first non data layer
    static const unsigned  vxdLayer = fDL + 0 ;
    static const unsigned  sitLayer = fDL + 1 ;
    static const unsigned  ftdLayer = fDL + 2 ;
    static const unsigned  tpcLayer = fDL + 3 ;
    static const unsigned ecalLayer = fDL + 4 ;
    static const unsigned ecalEndcapLayer = fDL + 5 ;
    static const unsigned hcalLayer = fDL + 6 ;
    static const unsigned hcalRingLayer = fDL + 7 ;
    static const unsigned hcalEndcapLayer = fDL + 8 ;
    static const unsigned yokeLayer = fDL + 9 ;
    static const unsigned coilLayer = fDL + 10 ;
    static const unsigned fcalLayer = fDL + 11 ;
    static const unsigned  setLayer = fDL + 12 ;
    
    
    //------------------ draw VXD first -------------------------
    
    
    for (int i=0; i<nLayersVTX; ++i) {
        
        int nLadders = pVXDLayerLayout->getNLadders(i);
        
        float _ladder_phi0 = float(pVXDLayerLayout->getPhi0(i));
        
        float _sensitive_distance = float(pVXDLayerLayout->getSensitiveDistance(i));
        float _sensitive_thickness = float(pVXDLayerLayout->getSensitiveThickness(i));
        float _sensitive_width = float(pVXDLayerLayout->getSensitiveWidth(i));
        
        float _sensitive_length = float(pVXDLayerLayout->getSensitiveLength(i)  * 2.  ); // lenght is half length really !!!
        
        float _sensitive_offset = float (pVXDLayerLayout->getSensitiveOffset(i));
        
        float currPhi;
        float angleLadders = 2*M_PI / nLadders;
        float cosphi, sinphi;
        
        
        //        _sensitive_distance +=0.5* _sensitive_thickness;
        //        _sensitive_distance -=0.065;
        
        for (int j=0; j<nLadders; ++j) {
            float vxd_radius=0;
            
            currPhi = _ladder_phi0 + (angleLadders * j);
            cosphi = cos(currPhi);
            sinphi = sin(currPhi);
            
            double r_offset = 0;
            
            
            switch (nLadders) {
                    
                case 18 :
                    r_offset = -1.14/2.0;
                    break;
                case 24 :
                    r_offset = -0.89/2.0;
                    break;
                case 30 :
                    r_offset = 0.81/2.0;
                    break;
                case 36 :
                    r_offset = 0.77/2.0;
                    break;
            }
            
            if(j%2==0){
                r_offset=-1.0*r_offset;
            }
            
            vxd_radius = _sensitive_distance+r_offset;
            
            double  sizes[3] ;
            double  center[3] ;
            unsigned int color = vxdCol;//0xFFFFFF;

            center[0] = (vxd_radius*cosphi - _sensitive_offset*sinphi);
            center[1] = (vxd_radius*sinphi + _sensitive_offset*cosphi);
            center[2] = 0.0;
            sizes[0]  = _sensitive_thickness;
            sizes[1]  = _sensitive_width;
            sizes[2]  = _sensitive_length ;
            
            double rotate[3];
            rotate[0] = 0.0;
            rotate[1] = 0.0;
            rotate[2] = currPhi*rad2deg;
            
            //ced_geobox_r( sizes, center, rotate, color, vxdLayer);
            
            ced_geobox_r_ID( sizes, center, rotate, color, vxdLayer,0);
            ced_geobox_r_solid( sizes, center, rotate, color, vxdLayer);
            
        }
    }
    
    
    //------------------ draw SIT Planar -------------------------
    
    
    // for (int i=0; i<nLayersSIT; ++i) {
    
    //   int nLadders = pSITLayerLayout->getNLadders(i);
    
    //   float _ladder_phi0 = float(pSITLayerLayout->getPhi0(i));
    
    //   float _sensitive_distance = float(pSITLayerLayout->getSensitiveDistance(i));
    //   float _sensitive_thickness = float(pSITLayerLayout->getSensitiveThickness(i));
    //   float _sensitive_width = float(pSITLayerLayout->getSensitiveWidth(i));
    
    //   float _sensitive_length = float(pSITLayerLayout->getSensitiveLength(i)  * 2.  ); // lenght is half length really !!!
    
    //   float _sensitive_offset = float (pSITLayerLayout->getSensitiveOffset(i));
    
    //   float currPhi;
    //   float angleLadders = 2*M_PI / nLadders;
    //   float cosphi, sinphi;
    
    //   _sensitive_distance +=0.5* _sensitive_thickness;
    
    //   for (int j=0; j<nLadders; ++j) {
    
    //     currPhi = _ladder_phi0 + (angleLadders * j);
    //     cosphi = cos(currPhi);
    //     sinphi = sin(currPhi);
    
    //     double  sizes[3] ;
    //     double  center[3] ;
    //     unsigned int color = 0xFFFFFF;
    
    //     center[0] = (_sensitive_distance*cosphi - _sensitive_offset*sinphi);
    //     center[1] = (_sensitive_distance*sinphi + _sensitive_offset*cosphi);
    //     center[2] = 0.0;
    //     sizes[0]  = _sensitive_thickness;
    //     sizes[1]  = _sensitive_width;
    //     sizes[2]  = _sensitive_length ;
    
    //     double rotate[3];
    //     rotate[0] = 0.0;
    //     rotate[1] = 0.0;
    //     rotate[2] = currPhi*rad2deg;
    
    //     ced_geobox_r( sizes, center, rotate, color, sitLayer);
    //     //      ced_geobox_r_solid( sizes, center, rotate, color, sitLayer);
    
    //   }
    // }
    
    
    
    //-----------------------------------------------------------
    
    std::vector<CEDGeoTube> gTV ;
    
    for( unsigned i=0, N = ftd_z.size(); i<N ; ++i) {
        gTV.push_back( CEDGeoTube( ftd_ri[i],          ftd_ro[i],  40,  40,   0.0   , 0, ftd_d[i]  ,    ftd_z[i] ,   ftdCol, ftdLayer, 0,0 ) ) ;  //  FTD
        gTV.push_back( CEDGeoTube( ftd_ri[i],          ftd_ro[i],  40,  40,   0.0   , 0, ftd_d[i]  ,  - ftd_z[i] ,   ftdCol, ftdLayer, 0,0 ) ) ;  //  FTD
    }
    
    // new sit
    for (int i=0; i<nLayersSIT; i++   ) {
        
        int nl_sit = pSITLayerLayout->getNLadders( i );
        float phi0_sit = float( pSITLayerLayout->getPhi0( i ) );
        float r_inn_sit = float( pSITLayerLayout->getSensitiveDistance( i  )  ) / cos( M_PI / nl_sit )  ;
        float r_out_sit = r_inn_sit ; //float( pSITLayerLayout->getSensitiveDistance( i + 1   ) ) / cos( M_PI / nl_sit )  ;
        float z_sit = float( pSITLayerLayout->getSensitiveLength( i ) ) ;
        
        gTV.push_back( CEDGeoTube( r_out_sit,     r_inn_sit,    nl_sit , nl_sit , phi0_sit , phi0_sit,  z_sit,   -z_sit,          sitCol, sitLayer ,0,1) ) ;  //  SIT
        
    }
    //old SIT
    for(unsigned i=0,N= rSIT.size() ; i<N ; ++i){
        gTV.push_back( CEDGeoTube( rSIT[i],          rSIT[i]-0.1 ,                 40, 40,  0.0, 0, lSIT[i],        -lSIT[i],            sitCol, sitLayer ,0,1) ) ;  //  SIT
    }
    
    
    // new set
    for (int i=0; i<nLayersSET; i++   ) {
        
        int nl_set = pSETLayerLayout->getNLadders( i );
        float phi0_set = float( pSETLayerLayout->getPhi0( i ) );
        float r_inn_set = float( pSETLayerLayout->getSensitiveDistance( i  )  ) / cos( M_PI / nl_set )  ;
        float r_out_set = r_inn_set ; //float( pSETLayerLayout->getSensitiveDistance( i+1 ) ) / cos( M_PI / nl_set )  ;
        float z_set = float( pSETLayerLayout->getSensitiveLength( i ) ) ;
        
        gTV.push_back( CEDGeoTube( r_out_set,     r_inn_set,    nl_set , nl_set , phi0_set , phi0_set,  z_set,   -z_set,          setCol, setLayer ,0,1) ) ;  //  SET
        
    }
    
    if(showBeamcal){
        gTV.push_back( CEDGeoTube( r_max_beamcal,      r_min_beamcal,              40, 40,    0., 0, thick_beamcal,  shift_beamcal_z_plus,   bcalCol, fcalLayer ,0,0) ) ; //    BEAMCAL +Z
        gTV.push_back( CEDGeoTube( r_max_beamcal,      r_min_beamcal,              40, 40, 0.,    0, thick_beamcal, -shift_beamcal_z_minus,  bcalCol, fcalLayer ,0,0) ) ;  //   BEAMCAL -Z
    }
    
    
    if(showTPC){
        gTV.push_back( CEDGeoTube( r_max_tpc,          r_min_tpc,                  40, 40,  0.0, 0, z_max_tpc,        -z_max_tpc,            tpcCol, tpcLayer ,1,1) ) ; //  TPC
    }
    
    if(showECAL){
        gTV.push_back( CEDGeoTube( r_out_ecal_bar,     r_inn_ecal_bar,              12,  12, 15, 0,  z_max_ecal_bar,   -z_max_ecal_bar,       ecalCol, ecalLayer ,0,1) ) ; //  ECAL Barrel
    }
    if(showECALEndcap) {
        gTV.push_back( CEDGeoTube( r_out_ecal_ecap,    r_inn_ecal_ecap,            12, 40, 15, 0,  thick_ecal_ecap,   shift_ecal_z_plus,    ecalCol, ecalEndcapLayer ,0,0) ) ; //  endcap ECAL +Z
        gTV.push_back( CEDGeoTube( r_out_ecal_ecap,    r_inn_ecal_ecap,            12, 40, 15, 0,  thick_ecal_ecap,  -shift_ecal_z_minus,   ecalCol, ecalEndcapLayer ,0,0) ) ; //  endcap ECAL -Z
    }
    
    if(showHCAL){
        gTV.push_back( CEDGeoTube( r_out_hcal_bar,     r_inn_hcal_bar,             12,  12, 15, 0, z_max_hcal_bar,  -z_max_hcal_bar,      hcalCol, hcalLayer ,0,1) ) ; //  HCAL Barrel
    }
    if(false) {
        gTV.push_back( CEDGeoTube( r_out_hcal_ring,    r_inn_hcal_ring,             40,  40,  15,     0, thick_hcal_ring,  shift_hcalr_z_plus,  hcalCol, hcalRingLayer ,0,1) ) ; //  ring HCAL +Z
        gTV.push_back( CEDGeoTube( r_out_hcal_ring,    r_inn_hcal_ring,             40,  40,  15,     0, thick_hcal_ring, -shift_hcalr_z_minus, hcalCol, hcalRingLayer ,0,1) ) ; //  ring HCAL -Z
    }
    if(showHCALEndcap) {
        gTV.push_back( CEDGeoTube( r_out_hcal_ecap,    r_min_hcal_ecap,             12, 40,  15,     0, thick_hcal_ecap,  shift_hcal_z_plus,   hcalCol, hcalEndcapLayer ,0,1) ) ; //  endcap HCAL +Z
        gTV.push_back( CEDGeoTube( r_out_hcal_ecap,    r_min_hcal_ecap,             12, 40,  15,     0, thick_hcal_ecap, -shift_hcal_z_minus,  hcalCol, hcalEndcapLayer ,0,1) ) ;  //  endcap HCAL -Z
    }
    
    if(showCoil){
        gTV.push_back( CEDGeoTube( coil_outer_radius,  coil_inner_radius,          40, 40,        0.0,   0, coil_half_z, -coil_half_z,              coilCol, coilLayer ,0,0) ) ;  //  coil
    }
    
    if(showYoke){
        gTV.push_back( CEDGeoTube( r_out_yoke_plug,    r_inn_yoke_plug,            12, 8,        15.0,   7.5, thick_yoke_plug,  shift_yoker_z_plus,  yokeCol, yokeLayer ,0,0) ) ; //  plug YOKE +Z
        gTV.push_back( CEDGeoTube( r_out_yoke_plug,    r_inn_yoke_plug,            12, 8,        15.0,   7.5, thick_yoke_plug, -shift_yoker_z_minus, yokeCol, yokeLayer ,0,0) ) ; //  plug YOKE -Z
    }
    
    
    
    if( showLHcal ){
        gTV.push_back( CEDGeoTube( r_max_lhcal,        r_min_lhcal,                40, 40,    0., 0, thick_lhcal,  shift_lhcal_z_plus,   fcalCol , fcalLayer ,0,0) ) ; //    LHCAL +Z
        gTV.push_back( CEDGeoTube( r_max_lhcal,        r_min_lhcal,                40, 40,    0., 0, thick_lhcal, -shift_lhcal_z_minus,  fcalCol , fcalLayer ,0,0) ) ;  //   LHCAL -Z
    }
    
    if ( showLCal ){
        gTV.push_back( CEDGeoTube( r_max_lcal,         r_min_lcal,                 40, 40,    0., 0, thick_lcal,  shift_lcal_z_plus,   lcalCol, fcalLayer ,0,0) ) ; //    LCAL +Z
        gTV.push_back( CEDGeoTube( r_max_lcal,         r_min_lcal,                 40, 40,    0., 0, thick_lcal, -shift_lcal_z_minus,  lcalCol, fcalLayer ,0,0) ) ;  //   LCAL -Z
    }
    
    
    if(showYoke){
        gTV.push_back( CEDGeoTube( r_out_yoke_bar,     r_inn_yoke_bar,             8, 8,        22.5,   0, z_max_yoke_bar,  - z_max_yoke_bar,     yokeCol,  yokeLayer, 0, 0) ) ; //  YOKE Barrel
        gTV.push_back( CEDGeoTube( r_out_yoke_ecap,    r_inn_yoke_ecap,            8, 8,        22.5,   0, thick_yoke_ecap,  shift_yoke_z_plus,   yokeCol,  yokeLayer, 0, 0) ) ; //  endcap YOKE +Z
        gTV.push_back( CEDGeoTube( r_out_yoke_ecap,    r_inn_yoke_ecap,            8, 8,        22.5,   0, thick_yoke_ecap, -shift_yoke_z_minus,  yokeCol,  yokeLayer, 0, 0) ) ;  //  endcap YOKE -Z
    }
    
    // ========================================================================
    
    ced_geotubes( gTV.size() ,  (CED_GeoTube*) &gTV[0] );
    
    // ========================================================================
    
    MarlinCED::set_layer_description("FTD", ftdLayer );
    MarlinCED::set_layer_description("VXD", vxdLayer );
    MarlinCED::set_layer_description("SIT", sitLayer );
    MarlinCED::set_layer_description("SET", setLayer );
    
    if(showTPC){
        MarlinCED::set_layer_description("TPC", tpcLayer );
    }
    if(showECAL){
        MarlinCED::set_layer_description("ECAL", ecalLayer );
    }
    if(showECALEndcap){
        MarlinCED::set_layer_description("ECALEndcap", ecalEndcapLayer );
    }
    if(showHCAL){
        MarlinCED::set_layer_description("HCAL", hcalLayer );
    }
    if(showHCALRing){
        MarlinCED::set_layer_description("HCALRing", hcalRingLayer );
    }
    if(showHCALEndcap){
        MarlinCED::set_layer_description("HCALEndcap", hcalEndcapLayer );
    }
    if(showCoil){
        MarlinCED::set_layer_description("Coil", coilLayer );
    }
    
    if(showYoke){
        MarlinCED::set_layer_description("Yoke", yokeLayer );
    }
    
    if( showLHcal  && showBeamcal){
        MarlinCED::set_layer_description("LCAL, Beamcal, LHcal", fcalLayer );
    }else if(showLHcal){
        MarlinCED::set_layer_description("LCAL, LHcal",fcalLayer );
    }else if(showBeamcal){
        MarlinCED::set_layer_description("LCAL, Beamcal",fcalLayer );
    }else if(showLCal){
        MarlinCED::set_layer_description("LCAL",fcalLayer );
    }
    
    MarlinCED::write_layer_description();
    
    
} // drawGEARDetector


