////////////////////////////////////////////////////////////////////////
// This Alg takes in the good hits from WCHitFinderAlg and creates    //
// a track with various kinematic and geometric values.  In effect,   //
// this alg is the track building part of WCTrackBuilderAlg_new.cxx       //
// in functions run after "FinalizeGoodHits".                         //
// Other algs can be created in lieu of this one if a better track    //
// is developed.  This is merely the version we used when hit finding //
// and track builder were run together in WCTrackBuilderAlg_new.cxx       //
// Author: Greg Pulliam gkpullia@syr.edu                              //
//////////////////////////////////////////////////////////////////////// 

//Framework includes
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "CLHEP/Units/SystemOfUnits.h"
#include "larcore/Geometry/AuxDetGeo.h"
#include "art/Framework/Services/Optional/TFileService.h"

// LArIAT includes
#include "LArIATRecoAlg/WCTrackBuilderAlg_new.h"


#include <iostream>
#include <cmath>
#include <cstdlib>
#include <TH1F.h>
#include <string>
#include <TH2F.h>



//--------------------------------------------------------------
//Constructor  TRACK
WCTrackBuilderAlg_new::WCTrackBuilderAlg_new( fhicl::ParameterSet const& pset )
{
  this->reconfigure(pset);

  //Testing the AuxDetGeo capabilitites
  std::vector<geo::AuxDetGeo*> const & theAuxDetGeoVect = fGeo->AuxDetGeoVec();
  double centerOfDet[3] = {0,0,0};
  for( size_t iDet = 0; iDet < fGeo->NAuxDets() ; ++iDet ){
    geo::AuxDetGeo* anAuxDetGeo = theAuxDetGeoVect.at(iDet);
    anAuxDetGeo->GetCenter(centerOfDet);

    /*
    std::cout << "************** AuxDet: " << iDet << "******************"  << std::endl;
    std::cout << "Length: " << anAuxDetGeo->Length() << std::endl;
    std::cout << "Half Width 1: " << anAuxDetGeo->HalfWidth1() << std::endl;
    std::cout << "Half Width 2: " << anAuxDetGeo->HalfWidth2() << std::endl;
    std::cout << "Test print" << std::endl;
    std::cout << "Center (x,y,z): (" 
	      << centerOfDet[0] << "," 
	      << centerOfDet[1] << ","
	      << centerOfDet[2] 
	      << ")" << std::endl;
    */

    //Setting the TGeo world locations of the MWPCs in mm
    if(iDet == 1){ //WC1
      fX_cntr_1 = centerOfDet[0] * CLHEP::cm;
      fY_cntr_1 = centerOfDet[1] * CLHEP::cm;
      fZ_cntr_1 = centerOfDet[2] * CLHEP::cm;
    }
   if(iDet == 2){ //WC2
      fX_cntr_2 = centerOfDet[0] * CLHEP::cm; //correcting where the gdml has an 18mm offset for WC2X
      fY_cntr_2 = centerOfDet[1] * CLHEP::cm;
      fZ_cntr_2 = centerOfDet[2] * CLHEP::cm;
    }
   if(iDet == 3){ //WC3
      fX_cntr_3 = centerOfDet[0] * CLHEP::cm;
      fY_cntr_3 = centerOfDet[1] * CLHEP::cm;
      fZ_cntr_3 = centerOfDet[2] * CLHEP::cm;
    }
   if(iDet == 4){ //WC4
      fX_cntr_4 = centerOfDet[0] * CLHEP::cm;
      fY_cntr_4 = centerOfDet[1] * CLHEP::cm;
      fZ_cntr_4 = centerOfDet[2] * CLHEP::cm;
    }
  }
  auto tpcGeo = fGeo->begin_TPC_id().get();
  double tpcLocalCenter[3] = {0.};
  tpcGeo->LocalToWorld(tpcLocalCenter, fCenter_of_tpc);
   // std::cout << "fCenter of TPC: "
  	//    << fCenter_of_tpc[0] << ","
  	  //  << fCenter_of_tpc[1] << ","
  	    //<< fCenter_of_tpc[2] << std::endl;
//std::cout<<CLHEP::cm<<std::endl;
  // put the center of the TPC in world coordinates into mm
  for(int i = 0; i < 3; ++i) {fCenter_of_tpc[i] *= CLHEP::cm;}
      std::cout << "fCenter of TPC: "
  	    << fCenter_of_tpc[0] << ","
  	    << fCenter_of_tpc[1] << ","
  	    << fCenter_of_tpc[2] << std::endl;
  float dxus= fX_cntr_2-fX_cntr_1;
  float dyus= fY_cntr_2-fY_cntr_1;
  float dzus= fZ_cntr_2-fZ_cntr_1;
  fDelta_z_us=pow(dxus*dxus+dyus*dyus+dzus*dzus,.5);
  float dxds= fX_cntr_4-fX_cntr_3;
  float dyds= fY_cntr_4-fY_cntr_3;
  float dzds= fZ_cntr_4-fZ_cntr_3;
  fDelta_z_ds=pow(dxds*dxds+dyds*dyds+dzds*dzds,.5);
  
  //std::cout<<"delta_us: "<<fDelta_z_us<<" delta_ds: "<<fDelta_z_ds<<std::endl;
  // get the active half width and length of the TPC in mm, geometry
  // returns the values in cm, so have to multiply by CLHEP::cm
  fHalf_z_length_of_tpc = 0.5*tpcGeo->ActiveLength() * CLHEP::cm;
  fHalf_x_length_of_tpc = tpcGeo->ActiveHalfWidth()  * CLHEP::cm;

}

//--------------------------------------------------------------  
//Destructor //BOTH
WCTrackBuilderAlg_new::~WCTrackBuilderAlg_new()
{

}

//--------------------------------------------------------------
void WCTrackBuilderAlg_new::reconfigure( fhicl::ParameterSet const& pset )  //BOTH
{

  fB_field_tesla        = pset.get<float >("BFieldInTesla",      0.       );


  fCentralYKink         = pset.get<float >("CentralYKink",        -0.01    ); //These four are parameters from histos I produced from picky-good tracks
  fSigmaYKink           = pset.get<float >("SigmaYKink",          0.03      );
  fCentralYDist         = pset.get<float >("CentralYDist",        0.69      );
  fSigmaYDist           = pset.get<float >("SigmaYDist",          18.0      );

  fPrintDisambiguation = false;
  fPickyTracks          = pset.get<bool  >("PickyTracks",         false     );
  fHighYield            = pset.get<bool  >("HighYield",           false     );
  
  //Survey constants
  fDelta_z_us           = pset.get<float >("DeltaZus",            1551.15   );  //this will recalculated using geometry instead of hardcoding. 
  fDelta_z_ds 		= pset.get<float >("DeltaZds",   	  1570.06   );   
  fL_eff        	= pset.get<float >("LEffective", 	  1145.34706);	
  fmm_to_m    		= pset.get<float >("MMtoM",        	  0.001     );	
  fGeV_to_MeV 		= pset.get<float >("GeVToMeV",    	  1000.0    );   
  fMP_X                 = pset.get<float> ("MidplaneInterceptFactor", 1);
  fMP_M                 = pset.get<float> ("MidplaneSlopeFactor", 1);

  
  return;
}
//--------------------------------------------------------------
//This is the function that is called to load correct row of the lariat_xml_database table for a run. This must be called within the beginSubRun() method of your analysis module
void WCTrackBuilderAlg_new::loadXMLDatabaseTableForBField( int run, int subrun ) //TRACK
{
  fRun = run;
  fSubRun = subrun;
  current=std::stod(fDatabaseUtility->GetIFBeamValue("mid_f_mc7an",fRun,fSubRun));
  //std::cout<<"Current: "<<current<<std::endl;
  if(fabs(current)>90){fB_field_tesla= .003375*current;}
  if(fabs(current)<90 && fabs(current)>70){fB_field_tesla= .0034875*current;}
  if(fabs(current)<70 && fabs(current)>50){fB_field_tesla= .003525*current;}
  if(fabs(current)<50 && fabs(current)>30){fB_field_tesla= .003525*current;}  
  if(fabs(current)<30){fB_field_tesla= .0035375*current;}  
 // fB_field_tesla = 0.0035*std::stod(fDatabaseUtility->GetIFBeamValue("mid_f_mc7an",fRun,fSubRun));
  std::cout << "Run: " << fRun << ", Subrun: " << fSubRun << ", B-field: " << fB_field_tesla << std::endl;
}

//--------------------------------------------------------------
//Main function called for each trigger
void WCTrackBuilderAlg_new::reconstructTracks(std::vector<double> & reco_pz_list,               
					   std::vector<double> & y_kink_list,
					   std::vector<double> & x_dist_list,
					   std::vector<double> & y_dist_list,
					   std::vector<double> & z_dist_list,
					   std::vector<double> & x_face_list,
					   std::vector<double> & y_face_list,
					   std::vector<double> & incoming_theta_list,
					   std::vector<double> & incoming_phi_list,
					   std::vector<WCHitList> & trigger_final_tracks,
					   std::vector<std::vector<WCHitList> > & good_hits,
					   bool verbose,
					   bool pickytracks,
					   bool highyield,
					   int & track_count,/* ,
					   std::vector<TH1F*> & WCHitsGoodTracks,
					   std::vector<TH2F*> & WCMult,
					   std::vector<TH1F*> & WireHits_TheTrack,
					   std::vector<TH1F*> & BadTrackHits,
					   TH2F* & TargetXY,
					   TH2F* & PickyTracksTargetXY,
					   TH1F* & ResSquare,
					   TH1F* & Reco4pt,
					   TH2F* & Reco4ptdiff,
					   std::vector<TH2F*> & TimingXY,
					   std::vector<TH2F*> & RegressionPlots,
					   std::vector<TH1F*> & RegressionPlots1D,
					   */std::vector<TH2F*> & Recoplots/*,
					   TH1F* & Bfield */)
{					   
  fVerbose = verbose;
  fPickyTracks = pickytracks;
  fHighYield = highyield;
  for(int i=0; i<3; ++i){
    BestYStats.push_back(999);
  }					 	
  //Determine if one should skip this trigger based on whether there is exactly one good hit in each wire chamber and axis
  //If there isn't, continue after adding a new empty vector to the reco_pz_array contianer.
  //If there is, then move on with the function.
  //This can be modified to permit more than one good hit in each wire chamber axis - see comments in function
  bool skip = shouldSkipTrigger(good_hits);
  if( skip == true ) return;
/*   Bfield->Fill(fB_field_tesla) */;
  
  //At this point, we should have a list of good hits with at least one good hit in X,Y for each WC.
  //Now find all possible combinations of these hits that could form a track, sans consideration
  //of the kinks or end displacements (for now). For the "exactly one" condition set in the above
  //step, this won't matter, but if you want to set the condition to "at least one hit in each WC axis,
  //this will give many combinations.
  bool lonely_hit_bool = buildTracksFromHits(good_hits,
					     reco_pz_list,
					     y_kink_list,
					     x_dist_list,
					     y_dist_list,
					     z_dist_list,
					     track_count,
					     x_face_list,
					     y_face_list,
					     incoming_theta_list,
					     incoming_phi_list,
					     trigger_final_tracks,/* ,
					     WCHitsGoodTracks,
					     WCMult,
					     BadTrackHits,
					     TargetXY,
					     PickyTracksTargetXY,
					     ResSquare,
					     Reco4pt
					     Reco4ptdiff,
					     TimingXY,
					     RegressionPlots,
					     RegressionPlots1D,*/
					     Recoplots);
					     
  //Need to use the cut information to whittle down track candidates
  if( !fPickyTracks ){
    disambiguateTracks( reco_pz_list,
			y_kink_list,
			x_dist_list,
			y_dist_list,
			z_dist_list,
			track_count,
			x_face_list,
			y_face_list,
			incoming_theta_list,
			incoming_phi_list,
			trigger_final_tracks,
			lonely_hit_bool);			
  }
/*  if(trigger_final_tracks.size()==1 && trigger_final_tracks[0].hits.size()==8){ //if we have one track, with 8 hits, fill the wire hits of the good track
   for(size_t iter=0; iter<8; ++iter){
     WireHits_TheTrack[iter]->Fill(trigger_final_tracks[0].hits[iter].wire);
     BadTrackHits[iter]->Fill(trigger_final_tracks[0].hits[iter].wire,-1);}
   } */
}

//=====================================================================
//Find the ykink, x/y/z_dist variables, and reco_pz TRACK
void WCTrackBuilderAlg_new::getTrackMom_Kink_End(WCHitList track,
					     float & reco_pz,
					     float & y_kink,
					     float (&dist_array)[3],/* ,
					     TH2F* & TargetXY,
					     TH1F* & ResSquare,
					     TH1F* & Reco4pt,
					     TH2F* & Reco4ptdiff,
					     std::vector<TH2F*> & RegressionPlots,
					     std::vector<TH1F*> & RegressionPlots1D,*/
					     std::vector<TH2F*> & Recoplots )
{
  std::vector<float> wire_x;
  std::vector<float> time_x;
  std::vector<float> wire_y;
  std::vector<float> time_y;
  
  //Loop through the hits in the track and fill the wire/time x/y
  for( size_t iHit = 0; iHit < track.hits.size() ; ++iHit ){
    int var = iHit % 2;
    if( var == 0 ){
	wire_x.push_back(track.hits.at(iHit).wire);
	time_x.push_back(track.hits.at(iHit).time);
	std::cout<<"for iHit: "<<iHit<<"x wire :"<<track.hits.at(iHit).wire<<std::endl;
      }
    if( var == 1 ){
	wire_y.push_back(track.hits.at(iHit).wire);
	time_y.push_back(track.hits.at(iHit).time);
	std::cout<<"for iHit: "<<iHit<<"y wire :"<<track.hits.at(iHit).wire<<std::endl;
      }   
  }
  float sin13=sin((3.141592654/180)*13.0);
  float sin3=sin((3.141592654/180)*3.0);
  float cos3=cos((3.141592654/180)*3.0);
  float cos13=cos((3.141592654/180)*13.0);

//For 4 Hit Tracks, we just need the angles of the upstream and downstream ends to calculate the momentum 
  // correct x and z for rotation of MWPCs about the vertical
    float x[4] = {fX_cntr_1 + wire_x[0] * cos13,
		  fX_cntr_2 + wire_x[1] * cos13,
		  fX_cntr_3 + wire_x[2] * cos3,
		  fX_cntr_4 + wire_x[3] * cos3 };
    //float y[4] = { fY_cntr_1 + wire_y[0],
//		 fY_cntr_2 + wire_y[1],
//		 fY_cntr_3 + wire_y[2],
//		 fY_cntr_4 + wire_y[3] };
    float z[4] = { fZ_cntr_1 + wire_x[0] * sin13,
		 fZ_cntr_2 + wire_x[1] * sin13,
		 fZ_cntr_3 + wire_x[2] * sin3,
 		 fZ_cntr_4 + wire_x[3] * sin3};
    float dx_us=x[1]-x[0];
    float dx_ds=x[3]-x[2];
   // float dy_us=y[1]-y[0];
   // float dy_ds=y[3]-y[2];
    float dz_us=z[1]-z[0];
    float dz_ds=z[3]-z[2];
    float theta_x_us= atan(dx_us/dz_us);
    float theta_x_ds= atan(dx_ds/dz_ds);
    //float theta_y_us= atan(dy_us/dz_us);
    //float theta_y_ds= atan(dy_ds/dz_ds);
  //Calculate angles to be used in momentum reco
  float delta_x_us = wire_x[1] - wire_x[0];
  float delta_y_us = wire_y[1] - wire_y[0];

  float delta_x_ds = wire_x[3] - wire_x[2];
  float delta_y_ds = wire_y[3] - wire_y[2];

  float atan_x_us = atan(delta_x_us / fDelta_z_us);
  float atan_x_ds = atan(delta_x_ds / fDelta_z_ds);
  //std::cout<<"us: "<< fDelta_z_us<<"ds: "<<fDelta_z_ds<<std::endl;
  float atan_y_us = atan(delta_y_us / fDelta_z_us);
  float atan_y_ds = atan(delta_y_ds / fDelta_z_ds);
     std::cout<<"tan us:"<<theta_x_us*180/3.1415<<"tan2 us: "<<atan2(dx_us,dz_us)*180/3.1415<<std::endl;
     std::cout<<"tan ds:"<<theta_x_ds*180/3.1415<<"tan2 ds: "<<atan2(dx_ds,dz_ds)*180/3.1415<<std::endl;  
  //Calculate momentum and y_kink
  //reco_pz = (fabs(fB_field_tesla) * fL_eff * fmm_to_m * fGeV_to_MeV ) / (float(3.3 * ((10.0*3.141592654/180.0) + (sin(atan_x_ds) - sin(atan_x_us))))) ;// cos(atan_y_ds);
  
  if(NHits==4){reco_pz = (fabs(fB_field_tesla) * fL_eff * fmm_to_m * fGeV_to_MeV ) / float(3.3 * ((10.0*3.141592654/180.0) + (atan_x_ds - atan_x_us)))/cos(atan_y_ds);} //use the normal 4 point momentum formula
   
    
  
  
  y_kink = atan_y_us - atan_y_ds;

  //Calculate the X/Y/Y Track End Distances
  float pos_us[3] = {0.0,0.0,0.0};
  float pos_ds[3] = {0.0,0.0,0.0};
  midPlaneExtrapolation(wire_x,wire_y,pos_us,pos_ds,reco_pz,WCMissed,Recoplots);
  for( int iDist = 0; iDist < 3 ; ++iDist ){
    dist_array[iDist] = pos_ds[iDist]-pos_us[iDist];
  }       
}


//=====================================================================
//More geometry //TRACK
void WCTrackBuilderAlg_new::midPlaneExtrapolation(std::vector<float> x_wires,
					      std::vector<float> y_wires,
					      float (&pos_us)[3],
					      float (&pos_ds)[3],
					      double  reco_pz,
					      int WCMissed,/*,
 					      TH2F* & TargetXY,
					      TH1F* & ResSquare,
					      TH1F* & Reco4pt,
					      TH2F* & Reco4ptdiff,
					      std::vector<TH2F*> & RegressionPlots,
					      std::vector<TH1F*> & RegressionPlots1D,*/
					      std::vector<TH2F*> & Recoplots)
{
  float sin13=sin((3.141592654/180)*13.0);
  float sin3=sin((3.141592654/180)*3.0);
  float cos3=cos((3.141592654/180)*3.0);
  float cos13=cos((3.141592654/180)*13.0);
  // correct x and z for rotation of MWPCs about the vertical
  float x[4] = { fX_cntr_1 + x_wires[0] * cos13,
		 fX_cntr_2 + x_wires[1] * cos13,
		 fX_cntr_3 + x_wires[2] * cos3,
		 fX_cntr_4 + x_wires[3] * cos3 };
  float y[4] = { fY_cntr_1 + y_wires[0],
		 fY_cntr_2 + y_wires[1],
		 fY_cntr_3 + y_wires[2],
		 fY_cntr_4 + y_wires[3] };
  float z[4] = { fZ_cntr_1 + x_wires[0] * sin13,
		 fZ_cntr_2 + x_wires[1] * sin13,
		 fZ_cntr_3 + x_wires[2] * sin3,
 		 fZ_cntr_4 + x_wires[3] * sin3};
  		 
  for(int iWC=0; iWC<4; ++iWC){
  std::cout<<"For WC: "<<iWC+1<<"The position is :"<<x[iWC]<<", "<<y[iWC]<<", "<<z[iWC]<<std::endl;
  }
  
//upstream leg of WC's tracked forward to midplane, using coordinate system with TPC as the origin  
  float slope_xz_us = (z[1] - z[0])/(x[1] - x[0]); 
  float z_int_xz_us_1 = z[0] - slope_xz_us * x[0];
  float z_int_xz_us_2 = z[1] - slope_xz_us * x[1];
  float z_int_xz_us = (z_int_xz_us_1 + z_int_xz_us_2) / 2.0;
  float x_us = (z_int_xz_us - fMid_plane_z_int_xz) / (fMid_plane_slope_xz - slope_xz_us);
  float z_us = fMid_plane_z_int_xz + fMid_plane_slope_xz * x_us; 
  float slope_yz_us = (y[1] - y[0]) / (z[1] - z[0]);
  float y_int_yz_us_1 = y[0] - slope_yz_us * z[0];
  float y_int_yz_us_2 = y[1] - slope_yz_us * z[1];
  float y_int_yz_us = (y_int_yz_us_1 + y_int_yz_us_2) / 2.0;
  float y_us = y_int_yz_us + slope_yz_us * z_us;    
//downstream leg of WCs tracked back to midplane, using coodinate system with TPC as the origin
  float slope_xz_ds = (z[3] - z[2]) / (x[3] - x[2]);
  float z_int_xz_ds_1 = z[2] - slope_xz_ds * x[2];
  float z_int_xz_ds_2 = z[3] - slope_xz_ds * x[3];
  float z_int_xz_ds = (z_int_xz_ds_1 + z_int_xz_ds_2) / 2.0;
    
  float x_ds = (z_int_xz_ds - fMid_plane_z_int_xz) / (fMid_plane_slope_xz - slope_xz_ds);
  float z_ds = fMid_plane_z_int_xz + fMid_plane_slope_xz * x_ds;

  float slope_yz_ds = (y[3] - y[2])/(z[3] - z[2]);
  float z_int_yz_ds_1 = y[2] - slope_yz_ds*z[2];
  float z_int_yz_ds_2 = y[3] - slope_yz_ds*z[3];
  float z_int_yz_ds = (z_int_yz_ds_1 + z_int_yz_ds_2)/2.0;

  float y_ds = z_int_yz_ds + slope_yz_ds * z_ds;
 

  pos_us[0] = x_us; pos_us[1] = y_us; pos_us[2] = z_us;
  pos_ds[0] = x_ds; pos_ds[1] = y_ds; pos_ds[2] = z_ds;  
  
 


//Hardcoding survey data from Run 1, Docdb 1371-final.  These are the coordinates of each wire chamber from target  
/*   float fX_cntr_target[4]={-394.081,-740.6132,-1064.9204,-1192.6316};
  float fY_cntr_target[4]={-.381, .0254, -3.1496, -20.4216};
  float fZ_cntr_target[4]={1698.244, 3196.0058, 5163.8484, 7602.0168};
   */
//TPC origin system calculation
//Better (i.e correct) version of calculating the midplane.  All lines (us, ds, mp) are done in the form X=mZ+b, so m must be some form of dx/dz
  //finding m, b for us leg
  float center_slope_us=(fX_cntr_2- fX_cntr_1)/(fZ_cntr_2 - fZ_cntr_1);
  float center_intercept_us=fX_cntr_2 - center_slope_us * fZ_cntr_2;
  //finding m,b for ds leg
  float center_slope_ds=(fX_cntr_4 - fX_cntr_3)/(fZ_cntr_4 - fZ_cntr_3);
  float center_intercept_ds=fX_cntr_4 - center_slope_ds * fZ_cntr_4;
  //At the midplane, these must intersect, setting the equations equal, we solve for Zmp, then plug back into either us, ds equation to get Xmp
  float midplane_z= (center_intercept_ds - center_intercept_us)/(center_slope_us - center_slope_ds);  
  float midplane_x= center_slope_us * midplane_z + center_intercept_us;
  //We expect the midplane to be oriented at 8 degrees, which means Mmp==cot(8), and we get Bmp from that, Xmp, Zmp
  float tan8=tan(8.0*3.141592654/180);
  std::cout<<"Midplane intercept factor: "<<fMP_X<<" Midplane slope factor: "<<fMP_M<<std::endl;
  float midplane_slope=tan8 * fMP_M;
  float midplane_intercept=(midplane_x-midplane_z*tan8) *fMP_X;
  std::cout<<"Midplane intercept:"<<midplane_intercept<<std::endl;
// 
// //TARGET origin system calculation
// 
// 
//     
//    //Better (i.e correct) version of calculating the midplane.  All lines (us, ds, mp) are done in the form X=mZ+b, so m must be some form of dx/dz
//   //finding m, b for us leg
//   float center_slope_us_tgt=(fX_cntr_target[1] - fX_cntr_target[0])/(fZ_cntr_target[1] - fZ_cntr_target[0]);
//   float center_intercept_us_tgt=fX_cntr_target[1] - center_slope_us_tgt * fZ_cntr_target[1];
//   //finding m,b for ds leg
//   float center_slope_ds_tgt=(fX_cntr_target[3] - fX_cntr_target[2])/(fZ_cntr_target[3] - fZ_cntr_target[2]);
//   float center_intercept_ds_tgt=fX_cntr_target[3] - center_slope_ds_tgt * fZ_cntr_target[3];
//   //At the midplane, these must intersect, setting the equations equal, we solve for Zmp, then plug back into either us, ds equation to get Xmp
//   float midplane_z_tgt= (center_intercept_ds_tgt - center_intercept_us_tgt)/(center_slope_us_tgt - center_slope_ds_tgt);  
//   float midplane_x_tgt= center_slope_us_tgt * midplane_z_tgt + center_intercept_us_tgt;
//   //We expect the midplane to be oriented at 8 degrees, which means Mmp==cot(8), and we get Bmp from that, Xmp, Zmp
//   float midplane_intercept_tgt=(midplane_x_tgt-midplane_z_tgt*tan8*fMP_M)*fMP_X; //fMP_M is the slope shift to calibrate where the midplane should be from where we currently have it
//   std::cout<<"TGT Midplane Intercept: "<<midplane_intercept_tgt<<std::endl;
//   
// //Find the position of the hit 		 
//   float X[4]={fX_cntr_target[0] + x_wires[0] * cos13 ,fX_cntr_target[1] + x_wires[1] * cos13 ,fX_cntr_target[2] + x_wires[2] * cos3 ,fX_cntr_target[3] + x_wires[3] * cos3};
//   float Y[4]={fY_cntr_target[0] + y_wires[0], fY_cntr_target[1] + y_wires[1], fY_cntr_target[2]+ y_wires[2], fY_cntr_target[3]+ y_wires[3]};
//   std::cout<<"Y[0]"<<Y[0]<<std::endl;
//   float Z[4]={fZ_cntr_target[0] + x_wires[0] * sin13, fZ_cntr_target[1] + x_wires[1] * sin13, fZ_cntr_target[2] + x_wires[2] * sin3, fZ_cntr_target[3] + x_wires[3] * sin3};
// //If we have a 4 point hit, we can calculate the angle between the us and ds legs and use that for the momentum.
//  
// 
//   float delta_y_ds = y_wires[3] - y_wires[2];
//   float atan_y_ds = atan(delta_y_ds / fDelta_z_ds); 
//   
  //To see how picky tracks works, I'm going to change Nhits and WCMissed so both versions of momentum is calculated to compare with the standard 4pt.  This is for checks only and should be deleted after verification. 
  NHits=3;
  WCMissed=2;
    float dy_us=y[1]-y[0];
    float dy_ds=y[3]-y[2];
    float dz_us=z[1]-z[0];
    float dz_ds=z[3]-z[2];
  float theta_y_us= atan(dy_us/dz_us);
  float theta_y_ds= atan(dy_ds/dz_ds);
  float S2_reco_pz;
  float S3_reco_pz;
  if(NHits==3){
  //float zmid=4205;
    if(WCMissed==2){
      float ds_dz=z[3]-z[2];
      float ds_dx=x[3]-x[2];
      float ds_slope=ds_dx/ds_dz;
      float ds_int_x= x[3]-ds_slope*z[3];
      float z_ds=(midplane_intercept-ds_int_x)/(ds_slope-midplane_slope); //Solving tan8*Z+Bmp == Mds*Z+Bds
      float x_ds=ds_slope*z_ds+ ds_int_x;
      std::cout<<"x hits "<<x[2]<<", "<<x[3]<<std::endl;
      std::cout<<"z hits "<<z[2]<<", "<<z[3]<<std::endl;
      std::cout<<"X centers "<<fX_cntr_3<<", "<<fX_cntr_4<<std::endl;
      std::cout<<"Z centers "<<fZ_cntr_3<<", "<<fZ_cntr_4<<std::endl;
      std::cout<<"ds_dx, ds_dz, slope, int"<<ds_dx<<", "<<ds_dz<<", "<<ds_slope<<", "<<ds_int_x<<std::endl;
      std::cout<<"midplane x,z"<<x_ds<<" "<<z_ds<<std::endl;
      float us_dz=(z_ds-z[0]);
      float us_dx=(x_ds-x[0]);       
      float us_theta=asin(us_dx/pow(us_dx*us_dx+us_dz*us_dz,.5));
      float ds_theta=asin(ds_dx/pow(ds_dx*ds_dx+ds_dz*ds_dz,.5));
      std::cout<<"s2 thetas :"<<us_theta<<", "<<ds_theta<<std::endl;
       S2_reco_pz=(fabs(fB_field_tesla) * fL_eff * fmm_to_m * fGeV_to_MeV ) / float(3.3 * (ds_theta - us_theta))/cos(theta_y_ds);
       std::cout<<"S2:"<<S2_reco_pz<<std::endl;
/*        Recoplots[0]->Fill(reco_pz,(reco_pz-S2_reco_pz)/S2_reco_pz);
       Recoplots[6]->Fill(reco_pz,S2_reco_pz); */
        if(current<70 && current> 50){
       
	 //Recoplots[2]->Fill(reco_pz,(reco_pz-S2_reco_pz)/S2_reco_pz);*/
       Recoplots[8]->Fill(reco_pz,S2_reco_pz);
       } 
        if(current>90){
	 //Recoplots[4]->Fill(reco_pz,(reco_pz-S2_reco_pz)/S2_reco_pz);
	 Recoplots[10]->Fill(reco_pz,S2_reco_pz);
       } 
              
    }
    WCMissed=3;
    if(WCMissed==3){
    
  
  
      float us_dz=z[1]-z[0];
      float us_dx=x[1]-x[0];
      float us_slope=us_dx/us_dz;
      float us_int_x= x[1]-us_slope*z[1];
      float z_us=(midplane_intercept-us_int_x)/(us_slope-midplane_slope); //Solving tan8*Z+Bmp == Mus*Z+Bus
      float x_us=us_slope*z_us + us_int_x;
      //std::cout<<"x_us_midplane: "<<x_us_target<<std::endl;
      //std::cout<<"z_us_midplane: "<<z_us_target<<std::endl;
      float ds_dz=-(z_us-z[3]);
      float ds_dx=-(x_us-x[3]);    
      float us_theta=asin(us_dx/pow(us_dx*us_dx+us_dz*us_dz,.5));
      float ds_theta=asin(ds_dx/pow(ds_dx*ds_dx+ds_dz*ds_dz,.5)); 
      //std::cout<<"us_theta: "<<us_theta_target<<std::endl;
      //std::cout<<"ds_theta: "<<ds_theta_target<<std::endl;
      S3_reco_pz=(fabs(fB_field_tesla) * fL_eff * fmm_to_m * fGeV_to_MeV ) / float(3.3 * (ds_theta - us_theta))/cos(theta_y_us);
      //std::cout<<"S3:"<<S3_reco_pz<<std::endl;
      // Recoplots[1]->Fill(reco_pz,(reco_pz-S3_reco_pz)/S3_reco_pz);
      //Recoplots[7]->Fill(reco_pz,S3_reco_pz);
       if(current<70 && current> 50){
	// Recoplots[3]->Fill(reco_pz,(reco_pz-S3_reco_pz)/S3_reco_pz);
         Recoplots[9]->Fill(reco_pz,S3_reco_pz);
       } 
        if(current>90){
	 //Recoplots[5]->Fill(reco_pz,(reco_pz-S3_reco_pz)/S3_reco_pz);
         Recoplots[11]->Fill(reco_pz,S3_reco_pz);	 
       } 
        
    }
  NHits=4;  
  WCMissed=-1;
  //Recoplots[0]->Fill(S2_reco_pz,reco_pz);
  //Recoplots[1]->Fill(S3_reco_pz,reco_pz);
} 

//Lets track back to the plane that intersects where the target should be, and find the y point where the track intersects.

/*   float slope_yz_us_tgt=(Y[1]-Y[0])/ (Z[1]-Z[0]);
  float y_int_yz_us_1_tgt = Y[1] - slope_yz_us_tgt * Z[1];



  float Target_intercept =X[1]-((X[1]-X[0])/(Z[1]-Z[0]))*Z[1];
  
  float Target_x=Target_intercept/(1-(sin13/cos13)*((X[1]-X[0])/(Z[1]-Z[0])));
  float Target_z=Target_intercept/(cos13/sin13-((X[1]-X[0])/(Z[1]-Z[0])));
 
  float Target_apparent=pow(pow(Target_x,2)+pow(Target_z,2),.5);
  
if the X position was negative, we need the distance to be negative as well.
  if(Target_x<0){Target_apparent=-Target_apparent;} */
/*   TargetXY->Fill(Target_apparent, y_int_yz_us_1_tgt); */

  

//Calculate the Residual skipping wc2 and wc3 to see which is best.  
/*   std::vector<float> residualsquare_skipsecond = Regression(y,z, ResSquare, 2);//REMEMBER TO DELETE THIS LATER!!!!!!!
  std::vector<float> residualsquare_skipthird = Regression(y,z, ResSquare, 3);//REMEMBER TO DELETE THIS LATER!!!!!!!
  std::vector<float> residualsquare_allfour = Regression(y,z, ResSquare, 0); //Using WC=0 so none get skipped
  RegressionPlots[0]->Fill(residualsquare_skipsecond[2], residualsquare_allfour[2]);
  RegressionPlots[1]->Fill(residualsquare_skipthird[2], residualsquare_allfour[2]);
  RegressionPlots[2]->Fill(residualsquare_skipsecond[2],residualsquare_skipthird[2]);
  RegressionPlots1D[0]->Fill(residualsquare_allfour[2]-residualsquare_skipsecond[2]);
  RegressionPlots1D[1]->Fill(residualsquare_allfour[2]-residualsquare_skipthird[2]);
  RegressionPlots1D[2]->Fill(residualsquare_skipthird[2]-residualsquare_skipsecond[2]); */
  
//This is the code to keep that will take the best Y points available, using 3 or 4 (if available) point regression
/*   if(fHighYield){
    if(NHits==4){
      std::vector<float> FourRegression=Regression(y,z, ResSquare, 0);
      std::vector<float> ThreeRegression_SkipTwo=Regression(y,z,ResSquare,2);
      std::vector<float> ThreeRegression_SkipThree=Regression(y,z,ResSquare,3);
      if(FourRegression[2]<ThreeRegression_SkipTwo[2] && FourRegression[2]<ThreeRegression_SkipThree[2] && FourRegression[2]<BestYStats[2]){
        BestYStats=FourRegression;
	BestYHits.clear();
	for(int i=0; i<4; ++i){
	BestYHits.push_back(y[i]);
	}
      }
      if(ThreeRegression_SkipTwo[2]<ThreeRegression_SkipThree[2] && ThreeRegression_SkipTwo[2]<FourRegression[2] && ThreeRegression_SkipTwo[2]<BestYStats[2]){
        BestYStats=ThreeRegression_SkipTwo;
	BestYHits.clear();
	for(int i=0; i<4; ++i){
	  if(i!=1){
	    BestYHits.push_back(y[i]);
	  }
	}
      }
      if(ThreeRegression_SkipThree[2]<ThreeRegression_SkipTwo[2] && ThreeRegression_SkipThree[2]<FourRegression[2] && ThreeRegression_SkipThree[2]<BestYStats[2]){
        BestYStats=ThreeRegression_SkipThree;
	BestYHits.clear();
	for(int i=0; i<4; ++i){
	  if(i!=2){
	    BestYHits.push_back(y[i]);
	  }
	}               	       
      }
    }
    if(NHits==3){
      std::vector<float> ThreeRegression=Regression(y,z,ResSquare,WCMissed);
      if(ThreeRegression[2]<BestYStats[2]){
        BestYStats=ThreeRegression;
	BestYHits.clear();
	for(int i=0; i<4; ++i){
	  if(i != WCMissed-1){
	    BestYHits.push_back(y[i]);
	  }
	}   
      }
    }  
  } */
}
//===================================================			   
/* //Testing my ability to do a regression and make a plot.
std::vector<float> WCTrackBuilderAlg_new::Regression(float (&y)[4],
						 float (&z)[4],
						 TH1F* & ResSquare,
						 int skippedWC)
{
 std::vector<float> RegressionValues;
  int Npoints=0;
  float sum_z=0;
  float sum_zz=0;
  float sum_y=0;
  float sum_yz=0;
  float intercept=0;
  float slope=0;
  float residualsquare=0;
  float residual=0;
    for(int i=0; i<4; ++i){
      if(i != skippedWC-1){ //turning WC# to array index
        sum_y  += y[i];
        sum_zz += z[i]*z[i];
        sum_z  += z[i];
        sum_yz += y[i]*z[i];  
	++Npoints; //Residual is a function of number of points used.
      }
    }
  slope=(Npoints*sum_yz-sum_y*sum_z)/(Npoints*sum_zz-sum_z*sum_z);
  intercept=(sum_y-slope*sum_z)/Npoints;
  for(int i=0; i<4;++i){
    if(i != skippedWC-1){
    residual= (y[i]-slope*z[i]-intercept)/std::sqrt(1+slope*slope);
    residualsquare += residual*residual/Npoints; 
    }
  }
  if(Npoints==4){
    ResSquare->Fill(residualsquare);
}
  RegressionValues.push_back(slope);
  RegressionValues.push_back(intercept);
  RegressionValues.push_back(residualsquare);
  return RegressionValues;
} */

//=====================================================================
//Taking the set of good hits and finding all combinations of possible tracks. These may not be physically
//reasonable, but could just be anything with a hit on each wire plane axis.
//TRACK!!!!!!!!!!!!!!!!!!!!!
bool WCTrackBuilderAlg_new::buildTracksFromHits(std::vector<std::vector<WCHitList> > & good_hits,
					    std::vector<double> & reco_pz_list,
					    std::vector<double> & y_kink_list,
					    std::vector<double> & x_dist_list,
					    std::vector<double> & y_dist_list,
					    std::vector<double> & z_dist_list,
					    int & track_count,
					    std::vector<double> & x_on_tpc_face_list,
					    std::vector<double> & y_on_tpc_face_list,
					    std::vector<double> & incoming_theta_list,
					    std::vector<double> & incoming_phi_list,
					    std::vector<WCHitList> & track_list,/* ,
					    std::vector<TH1F*> & GoodTrackWireHits,
					    std::vector<TH2F*> & WCMult,
					    std::vector<TH1F*> & BadTrackHits,
					    TH2F* & TargetXY,
					    TH2F* & PickyTracksTargetXY,
					    TH1F* & ResSquare,
					    TH1F* & Reco4pt,
					    TH2F* & Reco4ptdiff,
					    std::vector<TH2F*> & TimingXY,
					    std::vector<TH2F*> & RegressionPlots,
					    std::vector<TH1F*> & RegressionPlots1D,*/
					    std::vector<TH2F*> & Recoplots)
					    
					    
					   
{
  //Reconstructed momentum buffer for storing pz for all combinations in this trigger
  std::vector<double> reco_pz_buffer;
  int track_count_this_trigger = 0;

  //Clear semi-persistent vectors
  fGoodTrackCandidateErrors.clear();
  fGoodTrackCandidateHitLists.clear();

  bool lonely_hit_bool = false;   //Single hit on at most 7 WCAxes, at least 1
  for( size_t iWC = 0; iWC < 4 ; ++iWC ){
    for( size_t iAx = 0; iAx < 2 ; ++iAx ){
      if( good_hits.at(iWC).at(iAx).hits.size() == 1 ) lonely_hit_bool = true;
    }
  }

  //Loop through all combinations of tracks
  //int total_possible_tracks=0;
  //int track_counter=0;
  //total_possible_tracks=good_hits[0][0].hits.size()*good_hits[1][0].hits.size()*good_hits[2][0].hits.size()*good_hits[3][0].hits.size()*good_hits[0][1].hits.size()*good_hits[2][1].hits.size()*good_hits[1][1].hits.size()*good_hits[3][1].hits.size();
  //std::cout<<"Total number of track combos: "<<total_possible_tracks<<std::endl;
  for( size_t iHit0 = 0; iHit0 < good_hits[0][0].hits.size(); ++iHit0 ){
    for( size_t iHit1 = 0; iHit1 < good_hits[0][1].hits.size(); ++iHit1 ){
      for( size_t iHit2 = 0; iHit2 < good_hits[1][0].hits.size(); ++iHit2 ){
	for( size_t iHit3 = 0; iHit3 < good_hits[1][1].hits.size(); ++iHit3 ){
	  for( size_t iHit4 = 0; iHit4 < good_hits[2][0].hits.size(); ++iHit4 ){
	    for( size_t iHit5 = 0; iHit5 < good_hits[2][1].hits.size(); ++iHit5 ){
	      for( size_t iHit6 = 0; iHit6 < good_hits[3][0].hits.size(); ++iHit6 ){
		for( size_t iHit7 = 0; iHit7 < good_hits[3][1].hits.size(); ++iHit7 ){

		  //Push back a track
		  WCHitList track;
		  track.hits.push_back(good_hits[0][0].hits[iHit0]);
		  track.hits.push_back(good_hits[0][1].hits[iHit1]);
		  track.hits.push_back(good_hits[1][0].hits[iHit2]);
		  track.hits.push_back(good_hits[1][1].hits[iHit3]);
		  track.hits.push_back(good_hits[2][0].hits[iHit4]);
		  track.hits.push_back(good_hits[2][1].hits[iHit5]);
		  track.hits.push_back(good_hits[3][0].hits[iHit6]);
		  track.hits.push_back(good_hits[3][1].hits[iHit7]);
/*                   if(fPickyTracks){
		  TimingXY[0]->Fill(track.hits[0].time,track.hits[1].time);
		  TimingXY[1]->Fill(track.hits[2].time,track.hits[3].time);
		  TimingXY[2]->Fill(track.hits[4].time,track.hits[5].time);
		  TimingXY[3]->Fill(track.hits[6].time,track.hits[7].time);
		  } */

		  //Track reco info
		  float reco_pz = 0;
		  float y_kink = 0;
		  float dist_array[3] = {0,0,0};
		  float x_on_tpc_face = 0;
		  float y_on_tpc_face = 0;
		  float incoming_theta = 0;
		  float incoming_phi = 0;
		  		  
		  //Previously track_p_extrap_dists
		  getTrackMom_Kink_End(track,reco_pz,y_kink,dist_array,/* , TargetXY, ResSquare, Reco4pt, Reco4ptdiff, RegressionPlots, RegressionPlots1D, */Recoplots );

		  //Get things like x/y on the tpc face, theta, phi for track
		  findTrackOnTPCInfo(track,x_on_tpc_face,y_on_tpc_face,incoming_theta,incoming_phi);		

		  ////////////////////////////////////////////////////////////////////////////////////////////////
		  //                                                                                            //
		  // If you want to cut on quality of track, put a conditional here and use                     //
		  // the values of y_kink, x/y/z distance, etc. to do the cut. Using this, you                  //
		  // can decide whether to push back the track list (holds all of the hits for that             //
		  // track) and the lists of interesting track quantities (momenta, x/y face, theta/phi, etc.). //
		  //                                                                                            //
		  // The lists are updated to include all triggers, but don't worry - the module using the      //
		  // reconstructTracks function has a way of identifying which tracks come from each trigger.   //
		  // I can't remember what this ^ means, but it's obsolete now that we have the slicer. -Ryan   //
		  //                                                                                            //
		  ////////////////////////////////////////////////////////////////////////////////////////////////

		  //Seems nonintuitive, but when there are picky tracks, we want to get all of them,
		  //and when we cut, some of them get removed.
		  bool is_good_track = false;
		  if( !fPickyTracks )
		    is_good_track = cutOnGoodTracks(track,y_kink,dist_array,track_list.size());
		 		  		  
		  //Convert x on tpc face to convention
		  //		  x_on_tpc_face = x_on_tpc_face+fHalf_x_length_of_tpc;

		  //Add the track to the track list
		  
		  if( fPickyTracks || (!fPickyTracks && is_good_track) ){
		    track_list.push_back(track);
		    //Hit wires for WC1X, WC1Y, WC2X, WC2Y etc.....
		    //Also Fill the Histograms for Tracks we won't use.  Will subract off the wires of the "best" track from the good tracks
		    //leaving only the hits we discard when we disambiguateTracks.
/* 		    GoodTrackWireHits[0]->Fill(good_hits[0][0].hits[iHit0].wire);
		    GoodTrackWireHits[1]->Fill(good_hits[0][1].hits[iHit1].wire);
		    GoodTrackWireHits[2]->Fill(good_hits[1][0].hits[iHit2].wire);
		    GoodTrackWireHits[3]->Fill(good_hits[1][1].hits[iHit3].wire);
		    GoodTrackWireHits[4]->Fill(good_hits[2][0].hits[iHit4].wire);
		    GoodTrackWireHits[5]->Fill(good_hits[2][1].hits[iHit5].wire);
		    GoodTrackWireHits[6]->Fill(good_hits[3][0].hits[iHit6].wire);
		    GoodTrackWireHits[7]->Fill(good_hits[3][1].hits[iHit7].wire);
		    BadTrackHits[0]->Fill(good_hits[0][0].hits[iHit0].wire);
		    BadTrackHits[1]->Fill(good_hits[0][1].hits[iHit1].wire);
		    BadTrackHits[2]->Fill(good_hits[1][0].hits[iHit2].wire);
		    BadTrackHits[3]->Fill(good_hits[1][1].hits[iHit3].wire);
		    BadTrackHits[4]->Fill(good_hits[2][0].hits[iHit4].wire);
		    BadTrackHits[5]->Fill(good_hits[2][1].hits[iHit5].wire);
		    BadTrackHits[6]->Fill(good_hits[3][0].hits[iHit6].wire);
		    BadTrackHits[7]->Fill(good_hits[3][1].hits[iHit7].wire); */

		    
		    //Storing the momentum in the buffer that will be
		    //pushed back into the final reco_pz_array for this trigger
		    //POSSIBLY IRRELEVANT NOW - MAY NEED TO TRIM THIS IF HAVE TIME, BUT IT DOESN'T INTERFERE
		    reco_pz_buffer.push_back(reco_pz);
		    
		    //Filling full info lists
		    reco_pz_list.push_back(reco_pz);
		    y_kink_list.push_back(y_kink);
		    x_dist_list.push_back(dist_array[0]);
		    y_dist_list.push_back(dist_array[1]);
		    z_dist_list.push_back(dist_array[2]);
		    x_on_tpc_face_list.push_back(x_on_tpc_face);
		    y_on_tpc_face_list.push_back(y_on_tpc_face);
		    incoming_theta_list.push_back(incoming_theta);
		    incoming_phi_list.push_back(incoming_phi);
		    track_count_this_trigger++;
		    track_count++;
		  }
		}
	      }
	    }
	  }
	}
      }  
    }
  }
/*   int mult=1;
  std::vector<int> visitedwires;
  //std::cout<<"Track Count: "<<track_count<<std::endl;
  for(size_t planeIter=0; planeIter<8; ++planeIter){  
    for(size_t i=0; i<track_list.size(); ++i){
      bool DidChange = false;
      bool WireCounted = false;
      int the_wire=track_list[i].hits[planeIter].wire;
      //visitedwires.push_back(the_wire);
      track_list[i].hits[planeIter].isVisited=true;

      for(size_t j=track_list.size()-1; j>i; --j){
      
        if(track_list[j].hits[planeIter].wire == the_wire && track_list[j].hits[planeIter].isVisited == false){++mult;  DidChange = true; track_list[j].hits[planeIter].isVisited=true;
	
	}
      }
      if(DidChange)
         {
         WCMult[planeIter]->Fill(the_wire,mult);
	 }
      if(!DidChange){
        for(size_t wireiter=0; wireiter<visitedwires.size(); ++wireiter){
	  if(visitedwires[wireiter]==the_wire){WireCounted = true;}
	}
	if(WireCounted==false){WCMult[planeIter]->Fill(the_wire,1);}
	}
	   
      mult=1;
      visitedwires.push_back(the_wire);
    }
  }  */  

  
  //Clear the hit lists for each WC/axis
  for( size_t iWC = 0; iWC < good_hits.size() ; ++iWC ){
    for( size_t iAx = 0; iAx < good_hits.at(iWC).size() ; ++iAx ){
      good_hits.at(iWC).at(iAx).hits.clear();
    }
  }
  if( lonely_hit_bool ) return true;
  else{ return false; }
	     
}

//=====================================================================
//See if trigger has a good enough hit set to continue //TRACK!!!
bool WCTrackBuilderAlg_new::shouldSkipTrigger(std::vector<std::vector<WCHitList> > & good_hits)
{
  //Now determine if we want to skip
  bool skip = false;
  NHits=0;
  WCMissed=-1;
  for( size_t iWC = 0; iWC < 4 ; ++iWC ){
  if(good_hits[iWC][0].hits.size()>0 && good_hits[iWC][1].hits.size()>0){++NHits;}
  else{WCMissed=iWC+1;}
  } 
  if(fPickyTracks && fHighYield){
    for(size_t iWC=0; iWC<4; ++iWC){
      for(size_t iAX=0; iAX<2; ++iAX){
        if(good_hits[iWC][iAX].hits.size() > 1 || NHits<3 || WCMissed==0 || WCMissed==3){  
	//skip events that have more than 1 hit in an axis, have less than 3 X or Y hits, or the missed X is the first or last WC
          skip = true;
	  break;
	}
      }
    }   
  }
  if(!fPickyTracks && fHighYield){
    if(NHits<3 || WCMissed==0 || WCMissed==3){ //skip events with less than 3 X/Y hits or is missing the first or last WC
      skip = true;
    }  
  }
  for( size_t iWC = 0; iWC < good_hits.size() ; ++iWC ){
    for( size_t iAx = 0; iAx < good_hits.at(iWC).size() ; ++iAx ){
      if( fPickyTracks  && !fHighYield){
	if( good_hits.at(iWC).at(iAx).hits.size() != 1 ){ //require 1, unabigious 4 (x,y) track
	  skip = true;
	  break;
	}
      }
      if( !fPickyTracks && !fHighYield ){
	if( good_hits.at(iWC).at(iAx).hits.size() < 1 ){ //allow many 4 (x,y) point tracks and cut down later
	  skip = true;
	  break;
	}
      }
    }
  }
  
  if( skip == true ){
    if( fVerbose ){
      std::cout << "skipping this trigger." << std::endl;
    }
    //Clear the hit lists for each WC/axis
    for( size_t iWC = 0; iWC < good_hits.size() ; ++iWC ){
      for( size_t iAx = 0; iAx < good_hits.at(iWC).size() ; ++iAx ){
	good_hits.at(iWC).at(iAx).hits.clear();
      }
    }
    return true;
  }
  else return false;

}

//=====================================================================
//TRACKS!!!
void WCTrackBuilderAlg_new::findTrackOnTPCInfo(WCHitList track, float &x, float &y, float &theta, float &phi )
{
  
  //Get position vectors of the points on WC3 and WC4 
  float WC3_point[3] = { fX_cntr_3 + track.hits.at(4).wire*float(cos(3.141592654/180*(3.0))),
			   fY_cntr_3 + track.hits.at(5).wire,
			   fZ_cntr_3 + track.hits.at(4).wire*float(sin(3.141592654/180*(3.0))) };
  float WC4_point[3] = { fX_cntr_4 + track.hits.at(6).wire*float(cos(3.141592654/180*(3.0))),
			   fY_cntr_4 + track.hits.at(7).wire,
			   fZ_cntr_4 + track.hits.at(6).wire*float(sin(3.141592654/180*(3.0))) };


  //  transformWCHits(WC3_point,WC4_point);


  /*
  //OLD WAY OF FINDING X, Y, THETA, PHI AT FACE
  //Now have hit vectors in the frame of the TPC. Now we recreate the second track and find its
  //intersection with the upstream plane of the TPC. In this new frame, the upstream plane is just
  //Z = -450 mm. So we parametrize the track with t and find at which t Z = -450. We then use that
  //to get X and Y intercepts.  
  float parameter_t = (-1*fHalf_z_length_of_tpc-WC3_point[2])/(WC4_point[2]-WC3_point[2]);
  float x_at_US_plane = (WC3_point[0])+parameter_t*(WC4_point[0]-WC3_point[0]);
  float y_at_US_plane = (WC3_point[1])+parameter_t*(WC4_point[1]-WC3_point[1]);  
  x = x_at_US_plane;
  y = y_at_US_plane;
  float r = pow(pow(x-WC4_point[0],2)+pow(y-WC4_point[1],2),0.5);
  theta = atan(r/(-1*fHalf_z_length_of_tpc-WC4_point[2]));
  */

  //Now we have hit vectors in the frame of the TPC. Now we recreate the second track and find its
  //intersection with the upstream plane of the TPC. In this new frame, the upstream plane is just
  //Z = 0 mm, as it is in the TPC coordinate system
  float parameter_t = (WC3_point[2])/(WC3_point[2]-WC4_point[2]);
  float x_at_US_plane = (WC3_point[0])+parameter_t*(WC4_point[0]-WC3_point[0]);
  float y_at_US_plane = (WC3_point[1])+parameter_t*(WC4_point[1]-WC3_point[1]);  
  x = x_at_US_plane;
  y = y_at_US_plane;
  float r = pow(pow(x-WC4_point[0],2)+pow(y-WC4_point[1],2),0.5);
  theta = atan(-r/(WC4_point[2]));

  //Calculating phi (degeneracy elimination for the atan function)
  float dY = WC4_point[1]-WC3_point[1];
  float dX = WC4_point[0]-WC3_point[0];
  if( dY > 0 && dX > 0 ){ phi = atan(dY/dX); }
  else if( dY > 0 && dX < 0 ){ phi = atan(dY/dX)+3.141592654; }
  else if( dY < 0 && dX < 0 ){ phi = atan(dY/dX)+3.141592654; }
  else if( dY < 0 && dX > 0 ){ phi = atan(dY/dX)+6.28318; }
  else if( dY == 0 && dX == 0 ){ phi = 0; }//defined by convention
  else if( dY == 0 ){
    if( dX > 0 ){ phi = 0; }
    else{ phi = 3.141592654; }
  }
  else if( dX == 0 ){
    if( dY > 0 ){ phi = 3.141592654/2; }
    else{ phi = 3.141592654*3/2; }
  }
  else{ std::cout << "One should never reach here." << std::endl; }
  

}
  
//=====================================================================
//Transform these into the coordinate system of the TPC
//TRACKS!!!!!
void WCTrackBuilderAlg_new::transformWCHits( float (&WC3_point)[3],
		      float (&WC4_point)[3])
{
  //First transformation: a translation by the location of the TPC
  for( int iDim = 0; iDim < 3; ++iDim ){
    WC3_point[iDim] = WC3_point[iDim] - fCenter_of_tpc[iDim];
    WC4_point[iDim] = WC4_point[iDim] - fCenter_of_tpc[iDim];
  }
   
}

//=====================================================================
//Take a look at the track produced by the combinations of hits and determine
//if it's good enough to be considered a track.
//TRACKS!!!!
bool WCTrackBuilderAlg_new::cutOnGoodTracks( WCHitList track,
					 float & y_kink,
					 float (&dist_array)[3],
					 size_t track_index)
{
  //Hard cut on ykink
  if( y_kink > (fCentralYKink + fSigmaYKink) || y_kink < (fCentralYKink - fSigmaYKink) )
    return false;

  //Hard cut on ydist
  if( dist_array[1] > (fCentralYDist + fSigmaYDist ) || dist_array[1] < (fCentralYDist - fSigmaYDist) )
    return false;

  std::pair<float,float> the_error_pair(y_kink,dist_array[1]);
  std::pair<WCHitList,size_t> hitlist_pair(track,track_index);
  fGoodTrackCandidateErrors.emplace(track_index,the_error_pair);  
  fGoodTrackCandidateHitLists.push_back(hitlist_pair);
  return true;

}
					 
//=====================================================================
//For all tracks passing the hard cuts, narrow down on those with identical hits
//in any of the WCAxes
//TRACKS!!!!
void WCTrackBuilderAlg_new::disambiguateTracks( std::vector<double> & reco_pz_list,
					    std::vector<double> & y_kink_list,
					    std::vector<double> & x_dist_list,
					    std::vector<double> & y_dist_list,
					    std::vector<double> & z_dist_list,
					    int & track_count,
					    std::vector<double> & x_on_tpc_face_list,
					    std::vector<double> & y_on_tpc_face_list,
					    std::vector<double> & incoming_theta_list,
					    std::vector<double> & incoming_phi_list,
					    std::vector<WCHitList> & track_list,
					    bool lonely_hit_bool)
{
  if( fGoodTrackCandidateHitLists.size() > 1 ){
    if( lonely_hit_bool ){ //For now, only allow single tracks to be made
      if( fPrintDisambiguation ){
	std::cout << "Candidate Hit Lists Vector size: " << fGoodTrackCandidateHitLists.size() << std::endl;      
	for( size_t iTrack = 0; iTrack < fGoodTrackCandidateHitLists.size() ; ++iTrack ){
	  std::cout << "Candidate Hit List Index: " << fGoodTrackCandidateHitLists.at(iTrack).second << std::endl;
	}
      }
      
      //Loop through all possible good tracks and find the one with the lowest error
      float theSmallestError = 99999;
      float theSmallestErrorIndex = 99999;
      float scalingFactor = 0.001818;          //Used to give same weight to deviations in y_kink and y_dist
      for( size_t iHitList1 = 0; iHitList1 < fGoodTrackCandidateHitLists.size() ; ++iHitList1 ){
	float trackError = pow((pow(fGoodTrackCandidateErrors.at(fGoodTrackCandidateHitLists.at(iHitList1).second).first,2) +
				pow(fGoodTrackCandidateErrors.at(fGoodTrackCandidateHitLists.at(iHitList1).second).second*scalingFactor,2)),2);
	if( trackError < theSmallestError ){
	  theSmallestError = trackError;
	  theSmallestErrorIndex = fGoodTrackCandidateHitLists.at(iHitList1).second;
	}
      }

      //Keep only the good-index track from the initially good tracks (remove everything else)
      bool continue_bool = false;
      size_t tSize = track_list.size();
      for( size_t iTrack = 0; iTrack < tSize ; ++iTrack ){
	if( iTrack == theSmallestErrorIndex ){
	  continue_bool = true;
	  continue;
	}
	if( !continue_bool ){
	  reco_pz_list.erase(reco_pz_list.begin());
	  y_kink_list.erase(y_kink_list.begin());
	  x_dist_list.erase(x_dist_list.begin());
	  y_dist_list.erase(y_dist_list.begin());
	  z_dist_list.erase(z_dist_list.begin());
	  track_count--;
	  x_on_tpc_face_list.erase(x_on_tpc_face_list.begin());
	  y_on_tpc_face_list.erase(y_on_tpc_face_list.begin());
	  incoming_theta_list.erase(incoming_theta_list.begin());
	  incoming_phi_list.erase(incoming_phi_list.begin());
	  track_list.erase(track_list.begin());
	}
	else{
 	  reco_pz_list.pop_back();
	  y_kink_list.pop_back();
	  x_dist_list.pop_back();
	  y_dist_list.pop_back();
	  z_dist_list.pop_back();
	  track_count--;
	  x_on_tpc_face_list.pop_back();
	  y_on_tpc_face_list.pop_back();
	  incoming_theta_list.pop_back();
	  incoming_phi_list.pop_back();
	  track_list.pop_back();
	  
	}	
      }
    }
    else{ //Otherwise, kill every track - too complex to disambiguate at the moment.
      reco_pz_list.clear();
      y_kink_list.clear();
      x_dist_list.clear();
      y_dist_list.clear();
      z_dist_list.clear();
      track_count = 0;
      x_on_tpc_face_list.clear();
      y_on_tpc_face_list.clear();
      incoming_theta_list.clear();
      incoming_phi_list.clear();
      track_list.clear();
    }      
  }

  if( fPrintDisambiguation ){
    std::cout << "Number of tracks left in tracklist: " << track_list.size() << std::endl;
    std::cout << "Other list sizes: Pz:" << reco_pz_list.size() << ", Yk: " << y_kink_list.size() << ", y_dist: " << y_dist_list.size() << ", phi: " << incoming_phi_list.size() << std::endl;
  }
  
  

    /*
    //Loop through all possible ambiguous pairs of hitlists
    for( size_t iHitList1 = 0; iHitList1 < fGoodTrackCandidateHitLists.size() ; ++iHitList1 ){
      for( size_t iHitList2 = 0; iHitList2 < fGoodTrackCandidateHitLists.size() ; ++iHitList2 ){
	if( iHitList1 >= iHitList2 ) continue;
	//Loop through all of the hits in the two hit lists.
	//If any are the same, then run a comparison using the errors.
	for( size_t iHit1 = 0; iHit1 < fGoodTrackCandidateHitLists.at(iHitList1).first.hits.size(); ++iHit1 ){
	  bool break_bool = false;
	  for( size_t iHit2 = 0; iHit2 < fGoodTrackCandidateHitLists.at(iHitList2).first.hits.size(); ++iHit2 ){
	    if( (fGoodTrackCandidateHitLists.at(iHitList1).first.hits.at(iHit1).wire == 
		 fGoodTrackCandidateHitLists.at(iHitList2).first.hits.at(iHit2).wire) &&
		(fGoodTrackCandidateHitLists.at(iHitList1).first.hits.at(iHit1).time ==
		 fGoodTrackCandidateHitLists.at(iHitList2).first.hits.at(iHit2).time) ){
	      //If the hits are equal, identify which of the two tracks is better and push the bad one's index
	      //into a vector for later elimination
	      float t1FullError = pow((pow(fGoodTrackCandidateErrors.at(fGoodTrackCandidateHitLists.at(iHitList1).second).first,2) +
				       pow(fGoodTrackCandidateErrors.at(fGoodTrackCandidateHitLists.at(iHitList1).second).second,2)),2);
	      float t2FullError = pow((pow(fGoodTrackCandidateErrors.at(fGoodTrackCandidateHitLists.at(iHitList2).second).first,2) +
				       pow(fGoodTrackCandidateErrors.at(fGoodTrackCandidateHitLists.at(iHitList2).second).second,2)),2);
	      if( t1FullError <= t2FullError ){
		if( fPrintDisambiguation )
		  std::cout << "To-be-cut TrackIndex: " << fGoodTrackCandidateHitLists.at(iHitList2).second << std::endl; 
		bad_track_indices.push_back(fGoodTrackCandidateHitLists.at(iHitList2).second);
	      }
	      if( t1FullError > t2FullError ){
		if( fPrintDisambiguation )
		  std::cout << "To-be-cut TrackIndex: " << fGoodTrackCandidateHitLists.at(iHitList1).second << std::endl; 
		bad_track_indices.push_back(fGoodTrackCandidateHitLists.at(iHitList1).second);
	      }
	      break_bool = true;
	      break;
	    }
	    if( break_bool ) break;
	  }
	}
      }
    }
    
    //Debug printing
    if( fPrintDisambiguation ){
      for( size_t iBad = 0; iBad < bad_track_indices.size(); ++iBad ){
	std::cout << "Bad track index: " << bad_track_indices.at(iBad) << std::endl;
      }
    }


    //Go through the bad track indices and determine a unique list of track indices (since there may be repeats)
    std::vector<size_t> unique_bad_indices;
    for( size_t iInd = 0; iInd < bad_track_indices.size() ; ++iInd ){
      bool uniqueBool = true;
      for( size_t iUInd = 0; iUInd < unique_bad_indices.size() ; ++iInd ){
	if( bad_track_indices.at(iInd) == unique_bad_indices.at(iUInd) )
	  uniqueBool = false;
      }
      if( uniqueBool ) unique_bad_indices.push_back(bad_track_indices.at(iInd));
    }

    //Debug printing
    if( fPrintDisambiguation ){
      for( size_t iBad = 0; iBad < unique_bad_indices.size(); ++iBad ){
	std::cout << "Unique Bad track index: " << unique_bad_indices.at(iBad) << std::endl;
      }
    }
    	  
    //Go through the bad track indices and remove those tracks from the existing lists
    for( size_t iBad = 0; iBad < unique_bad_indices.size(); ++iBad ){
      reco_pz_list.erase(reco_pz_list.begin()+unique_bad_indices.at(iBad));
      y_kink_list.erase(y_kink_list.begin()+unique_bad_indices.at(iBad));
      x_dist_list.erase(x_dist_list.begin()+unique_bad_indices.at(iBad));
      y_dist_list.erase(y_dist_list.begin()+unique_bad_indices.at(iBad));
      z_dist_list.erase(z_dist_list.begin()+unique_bad_indices.at(iBad));
      track_count--;
      x_on_tpc_face_list.erase(x_on_tpc_face_list.begin()+unique_bad_indices.at(iBad));
      y_on_tpc_face_list.erase(y_on_tpc_face_list.begin()+unique_bad_indices.at(iBad));
      incoming_theta_list.erase(incoming_theta_list.begin()+unique_bad_indices.at(iBad));
      incoming_phi_list.erase(incoming_phi_list.begin()+unique_bad_indices.at(iBad));
      track_list.erase(track_list.begin()+unique_bad_indices.at(iBad));
      //For all additional unique bad indices that are greater than that corresponding
      //to iBad's, lower them by one, since the vectors have been shortened by one.
      for( size_t iUB = 0; iUB < unique_bad_indices.size(); ++iUB ){
	if( unique_bad_indices.at(iUB) > unique_bad_indices.at(iBad) ){
	  unique_bad_indices.at(iUB)--;
	}
      }
    }
      
  }
    */
}	
					     
