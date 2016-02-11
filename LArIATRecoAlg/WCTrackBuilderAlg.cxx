////////////////////////////////////////////////////////////////
//                                                            //
// This is a class definition for the wire chamber track      //
// builder algorithm, used to reconstruct momentum and other  //
// geometrical properties of test-beam particles passing      //
// through LArIAT's four wire chambers                        //
//                                                            //
// Authors: Ryan Linehan, rlinehan@stanford.edu               //                           
//          Johnny Ho, johnnyho@uchicago.edu                  //
//          Jason St. John, stjohn@fnal.gov                   //
//                                                            //
////////////////////////////////////////////////////////////////

//Framework includes
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "CLHEP/Units/SystemOfUnits.h"
#include "larcore/Geometry/AuxDetGeo.h"
#include "art/Framework/Services/Optional/TFileService.h"

// LArIAT includes
#include "LArIATRecoAlg/WCTrackBuilderAlg.h"


#include <iostream>
#include <cmath>
#include <cstdlib>
#include <TH1F.h>
#include <string>

//--------------------------------------------------------------
//Constructor
WCTrackBuilderAlg::WCTrackBuilderAlg( fhicl::ParameterSet const& pset )
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
      fX_cntr_2 = centerOfDet[0] * CLHEP::cm;
      fY_cntr_2 = centerOfDet[1] * CLHEP::cm;
      fZ_cntr_2 = centerOfDet[2] * CLHEP::cm;
    }
   if(iDet == 3){ //WC1
      fX_cntr_3 = centerOfDet[0] * CLHEP::cm;
      fY_cntr_3 = centerOfDet[1] * CLHEP::cm;
      fZ_cntr_3 = centerOfDet[2] * CLHEP::cm;
    }
   if(iDet == 4){ //WC2
      fX_cntr_4 = centerOfDet[0] * CLHEP::cm;
      fY_cntr_4 = centerOfDet[1] * CLHEP::cm;
      fZ_cntr_4 = centerOfDet[2] * CLHEP::cm;
    }
  }


    
  // rough guess beginning 2015.02: moved WC4 33" along the 3deg beam.
  // x = -1146.64 mm = -1102.77 mm - sin(3 deg)*33"
  // z =  7587.99 mm =  6750.94 mm + cos(3 deg)*33"
  // Positions from 2014 runs: (fX_cntr_4, fY_cntr_4, fZ_cntr_4) = (-1102.77, 0.0, 6750.94)
  
  // derive some useful geometric/al constants
  
  // upstream (us) leg center line projected onto xz-plane
  fCntr_slope_xz_us = (fZ_cntr_2 - fZ_cntr_1)/(fX_cntr_2 - fX_cntr_1);
  fCntr_z_int_xz_us_1 = fZ_cntr_1 - (fCntr_slope_xz_us * fX_cntr_1);
  fCntr_z_int_xz_us_2 = fZ_cntr_2 - (fCntr_slope_xz_us * fX_cntr_2);
  fCntr_z_int_xz_us = (fCntr_z_int_xz_us_1 + fCntr_z_int_xz_us_2) / 2.0;
  
  // downstream (ds) leg center line projected onto xz-plane
  fCntr_slope_xz_ds = (fZ_cntr_4 - fZ_cntr_3)/(fX_cntr_4 - fX_cntr_3);
  fCntr_z_int_xz_ds_3 = fZ_cntr_3 - (fCntr_slope_xz_ds * fX_cntr_3);
  fCntr_z_int_xz_ds_4 = fZ_cntr_4 - (fCntr_slope_xz_ds * fX_cntr_4);
  fCntr_z_int_xz_ds = (fCntr_z_int_xz_ds_3 + fCntr_z_int_xz_ds_4) / 2.0;
  
  // mid-plane
  // includes intersection point of upstream and downstream
  // center lines above and taken to be at 8 y-degrees from the
  // xy-plane
  fMid_plane_x = -1.0 * (fCntr_z_int_xz_ds - fCntr_z_int_xz_us) / (fCntr_slope_xz_ds - fCntr_slope_xz_us);
  
  fMid_plane_z_us = fCntr_z_int_xz_us + fCntr_slope_xz_us * fMid_plane_x;
  fMid_plane_z_ds = fCntr_z_int_xz_ds + fCntr_slope_xz_ds * fMid_plane_x;
  fMid_plane_z = (fMid_plane_z_us + fMid_plane_z_ds) / 2.0;
  fMid_plane_slope_xz = tan(8.0*3.141592654/180);
  fMid_plane_z_int_xz = fMid_plane_z - fMid_plane_slope_xz * fMid_plane_x;

  // Use the Geometry service to get positions and lengths
  // of the TPC
  auto tpcGeo = fGeo->begin_TPC_id().get();
  double tpcLocalCenter[3] = {0.};
  tpcGeo->LocalToWorld(tpcLocalCenter, fCenter_of_tpc);
  //  std::cout << "fCenter of TPC: "
  //	    << fCenter_of_tpc[0] << ","
  //	    << fCenter_of_tpc[1] << ","
  //	    << fCenter_of_tpc[2] << std::endl;

  // put the center of the TPC in world coordinates into mm
  for(int i = 0; i < 3; ++i) fCenter_of_tpc[i] *= CLHEP::cm;

  // get the active half width and length of the TPC in mm, geometry
  // returns the values in cm, so have to multiply by CLHEP::cm
  fHalf_z_length_of_tpc = 0.5*tpcGeo->ActiveLength() * CLHEP::cm;
  fHalf_x_length_of_tpc = tpcGeo->ActiveHalfWidth()  * CLHEP::cm;

}

//--------------------------------------------------------------  
//Destructor
WCTrackBuilderAlg::~WCTrackBuilderAlg()
{

}

//--------------------------------------------------------------
void WCTrackBuilderAlg::reconfigure( fhicl::ParameterSet const& pset )
{
  fNumber_tdcs          = pset.get<int   >("NumberTDCs",         16         );
  fNumber_wire_chambers = pset.get<int   >("NumberWireChambers", 4          );
  fB_field_tesla        = pset.get<float >("BFieldInTesla",      0.       );

  fTime_bin_scaling     = pset.get<double>("TimeBinScaling",      1.0/1280.0);
  fWire_scaling         = pset.get<double>("WireScaling",         1.0/64.0  ); 
  fGoodHitAveragingEps 	= pset.get<double>("GotHitAveragingEps",  2.0       );

  fDBSCANEpsilon        = pset.get<float >("DBSCANEpsilon",       1.0/16.0  );
  fDBSCANMinHits        = pset.get<int   >("DBSCANMinHits",       1         );

  fCentralYKink         = pset.get<float >("CentralYKink",        -0.01    ); //These four are parameters from histos I produced from picky-good tracks
  fSigmaYKink           = pset.get<float >("SigmaYKink",          0.03      );
  fCentralYDist         = pset.get<float >("CentralYDist",        0.69      );
  fSigmaYDist           = pset.get<float >("SigmaYDist",          18.0      );

  fPrintDisambiguation = false;
  fPickyTracks          = pset.get<bool  >("PickyTracks",         false     );
  
  //Survey constants
  fDelta_z_us           = pset.get<float >("DeltaZus",            1551.15   );   
  fDelta_z_ds 		= pset.get<float >("DeltaZds",   	  1570.06   );   
  fL_eff        	= pset.get<float >("LEffective", 	  1145.34706);	
  fmm_to_m    		= pset.get<float >("MMtoM",        	  0.001     );	
  fGeV_to_MeV 		= pset.get<float >("GeVToMeV",    	  1000.0    );   

  fTrack_Type = 9999;
  
  return;
}

//--------------------------------------------------------------
int WCTrackBuilderAlg::getTrackType()
{
  return fTrack_Type;
}

//--------------------------------------------------------------
//This is the function that is called to load correct row of the lariat_xml_database table for a run. This must be called within the beginSubRun() method of your analysis module
void WCTrackBuilderAlg::loadXMLDatabaseTableForBField( int run, int subrun )
{
  fRun = run;
  fSubRun = subrun;
  fB_field_tesla = 0.0035*std::stod(fDatabaseUtility->GetIFBeamValue("mid_f_mc7an",fRun,fSubRun));
  std::cout << "Run: " << fRun << ", Subrun: " << fSubRun << ", B-field: " << fB_field_tesla << std::endl;
}

//--------------------------------------------------------------
//Main function called for each trigger
void WCTrackBuilderAlg::reconstructTracks( std::vector<int> tdc_number_vect,
					   std::vector<float> hit_channel_vect,
					   std::vector<float> hit_time_bin_vect,
					   std::vector<std::vector<double> > & reco_pz_array,
					   std::vector<double> & reco_pz_list,               
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
					   int & good_trigger_counter,
					   int trigger_number,
					   int & track_count)
{
  fVerbose = verbose;
  
  //Initialize hit/cluster time/channel buffers that are only defined for one trigger value
  std::vector<std::vector<float> > hit_time_buffer;
  std::vector<std::vector<float> > hit_wire_buffer;
  std::vector<std::vector<float> > cluster_time_buffer;
  std::vector<std::vector<float> > cluster_wire_buffer;

  //Create a set of buffers with 1st-Dim length of 8 (for 8 wire chamber axes: 1X, 1Y, 2X, 2Y, etc.)
  initializeBuffers( hit_time_buffer, hit_wire_buffer, cluster_time_buffer, cluster_wire_buffer );
  
  //Fill the buffers with the hit times and wires for this trigger value
  fillTimeAndWireBuffers( tdc_number_vect, hit_time_buffer, hit_wire_buffer, hit_time_bin_vect, hit_channel_vect );
  
  //Create a vector of clusters with DBSCAN, where each entry into cluster_time_buffer/cluster_wire_buffer is a different cluster
  createClusters( trigger_number, hit_time_buffer, hit_wire_buffer, cluster_time_buffer, cluster_wire_buffer );
  
  //Finding the hits to be used in momentum reconstruction
  //Note here that only one hit may be accepted from a cluster. The idea is that a particle passing through the MWPC causes noise
  //that is spatially and temporally clustered around the initial hit of the particle. 
  //In addition, if sufficiently close together, two or more good hits can be averaged  
  findGoodHits(cluster_time_buffer,cluster_wire_buffer,good_hits);
  
  //Determine if one should skip this trigger based on whether there is at least one good hit in each wire chamber and axis
  //If there isn't, continue after adding a new empty vector to the reco_pz_array contianer
  //If there is, then move on with the function 

  /*  
  //Sanity check
  std::cout << "Good hit check: ";
  for( int iWC = 0; iWC < 4 ; ++iWC ){
    for (int iAx = 0; iAx < 2 ; ++iAx ){
      if( good_hits.at(iWC).at(iAx).hits.size() > 0 ) std::cout << good_hits.at(iWC).at(iAx).hits.size();
      else std::cout << "0";
    }
  }
  std::cout << std::endl;
  */  

  
  //Determine if one should skip this trigger based on whether there is exactly one good hit in each wire chamber and axis
  //If there isn't, continue after adding a new empty vector to the reco_pz_array contianer.
  //If there is, then move on with the function.
  //This can be modified to permit more than one good hit in each wire chamber axis - see comments in function
  bool skip = shouldSkipTrigger(good_hits,reco_pz_array);
  if( fVerbose ){
    if( skip == false ){ 
      std::cout << "Trigger: " << trigger_number << " is a good trigger." << std::endl;
      good_trigger_counter++;
    }
  }
  if( skip == true ) return;
  
  
  //At this point, we should have a list of good hits with at least one good hit in X,Y for each WC.
  //Now find all possible combinations of these hits that could form a track, sans consideration
  //of the kinks or end displacements (for now). For the "exactly one" condition set in the above
  //step, this won't matter, but if you want to set the condition to "at least one hit in each WC axis,
  //this will give many combinations.
  bool lonely_hit_bool = buildTracksFromHits(good_hits,
					     reco_pz_array,
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
					     trigger_final_tracks);
  
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
}


//=====================================================================
//Find the ykink, x/y/z_dist variables, and reco_pz
void WCTrackBuilderAlg::getTrackMom_Kink_End(WCHitList track,
					     float & reco_pz,
					     float & y_kink,
					     float (&dist_array)[3])
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
      }
    if( var == 1 ){
	wire_y.push_back(track.hits.at(iHit).wire);
	time_y.push_back(track.hits.at(iHit).time);
      }   
  }
  
  //Calculate angles to be used in momentum reco
  float delta_x_us = wire_x[1] - wire_x[0];
  float delta_y_us = wire_y[1] - wire_y[0];

  float delta_x_ds = wire_x[3] - wire_x[2];
  float delta_y_ds = wire_y[3] - wire_y[2];

  float atan_x_us = atan(delta_x_us / fDelta_z_us);
  float atan_x_ds = atan(delta_x_ds / fDelta_z_ds);

  float atan_y_us = atan(delta_y_us / fDelta_z_us);
  float atan_y_ds = atan(delta_y_ds / fDelta_z_ds);

  //Calculate momentum and y_kink
  reco_pz = (fabs(fB_field_tesla) * fL_eff * fmm_to_m * fGeV_to_MeV ) / float(3.3 * ((10.0*3.141592654/180.0) + (atan_x_ds - atan_x_us))) / cos(atan_y_ds);


  y_kink = atan_y_us - atan_y_ds;

  //Calculate the X/Y/Y Track End Distances
  float pos_us[3] = {0.0,0.0,0.0};
  float pos_ds[3] = {0.0,0.0,0.0};
  midPlaneExtrapolation(wire_x,wire_y,pos_us,pos_ds);
  for( int iDist = 0; iDist < 3 ; ++iDist ){
    dist_array[iDist] = pos_ds[iDist]-pos_us[iDist];
  }       
}


//=====================================================================
//More geometry
void WCTrackBuilderAlg::midPlaneExtrapolation(std::vector<float> x_wires,
					      std::vector<float> y_wires,
					      float (&pos_us)[3],
					      float (&pos_ds)[3])
{
  // correct x and z for rotation of MWPCs about the vertical
  float x[4] = { fX_cntr_1 + x_wires[0] * float(cos(3.141592654/180*(13.0))),
		 fX_cntr_2 + x_wires[1] * float(cos(3.141592654/180*(13.0))),
		 fX_cntr_3 + x_wires[2] * float(cos(3.141592654/180*(3.0))),
		 fX_cntr_4 + x_wires[3] * float(cos(3.141592654/180*(3.0))) };
  float y[4] = { fY_cntr_1 + y_wires[0],
		 fY_cntr_2 + y_wires[1],
		 fY_cntr_3 + y_wires[2],
		 fY_cntr_4 + y_wires[3] };
  float z[4] = { fZ_cntr_1 + x_wires[0] * float(sin(3.141592654/180*(13.0))),
		 fZ_cntr_2 + x_wires[1] * float(sin(3.141592654/180*(13.0))),
		 fZ_cntr_3 + x_wires[2] * float(sin(3.141592654/180*(3.0))),
		 fZ_cntr_4 + x_wires[3] * float(sin(3.141592654/180*(3.0))) };
  
  // upstream (us) leg, extrapolating forward to mid-plane
  float slope_xz_us = (z[1] - z[0])/(x[1] - x[0]);
  float z_int_xz_us_1 = z[0] - slope_xz_us * x[0];
  float z_int_xz_us_2 = z[1] - slope_xz_us * x[1];
  float z_int_xz_us = (z_int_xz_us_1 + z_int_xz_us_2) / 2.0;

  float x_us = (z_int_xz_us - fMid_plane_z_int_xz) / (fMid_plane_slope_xz - slope_xz_us);
  float z_us = fMid_plane_z_int_xz + fMid_plane_slope_xz * x_us;

  float slope_yz_us = (y[1] - y[0]) / (z[1] - z[0]);
  float z_int_yz_us_1 = y[0] - slope_yz_us * z[0];
  float z_int_yz_us_2 = y[1] - slope_yz_us * z[1];
  float z_int_yz_us = (z_int_yz_us_1 + z_int_yz_us_2) / 2.0;

  float y_us = z_int_yz_us + slope_yz_us * z_us;

  //downstream (ds) leg, extrapolating back to mid-plane
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

}
			   

//=====================================================================
//Taking the set of good hits and finding all combinations of possible tracks. These may not be physically
//reasonable, but could just be anything with a hit on each wire plane axis.
bool WCTrackBuilderAlg::buildTracksFromHits(std::vector<std::vector<WCHitList> > & good_hits,
					    std::vector<std::vector<double> > & reco_pz_array,
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
					    std::vector<WCHitList> & track_list)
					    
					    
					   
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


		  //Track reco info
		  float reco_pz = 0;
		  float y_kink = 0;
		  float dist_array[3] = {0,0,0};
		  float x_on_tpc_face = 0;
		  float y_on_tpc_face = 0;
		  float incoming_theta = 0;
		  float incoming_phi = 0;
		  		  
		  //Previously track_p_extrap_dists
		  getTrackMom_Kink_End(track,reco_pz,y_kink,dist_array);

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

  //Fill the array with the found tracks' reconstructed momenta //AGAIN, OBSOLETE, BUT IT DOESN'T HURT ANYTHING TECHNICALLY
  reco_pz_array.push_back(reco_pz_buffer);
  
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
//See if trigger has a good enough hit set to continue
bool WCTrackBuilderAlg::shouldSkipTrigger(std::vector<std::vector<WCHitList> > & good_hits,
					    std::vector<std::vector<double> > & reco_pz_array)
{
  //Check to see if there is only a single hit on one of the WCAxes
  bool lonely_hit_bool = false;   //Single hit on at most 7 WCAxes, at least 1
  bool unique_hit_bool = true;   //Single hit on all WCAxes
  bool missing_hit_bool = false;  //Missing hit on 1 or more WCAxes
  for( size_t iWC = 0; iWC < 4 ; ++iWC ){
    for( size_t iAx = 0; iAx < 2 ; ++iAx ){
      if( good_hits.at(iWC).at(iAx).hits.size() == 0 ) missing_hit_bool = true;
      if( good_hits.at(iWC).at(iAx).hits.size() == 1 ) lonely_hit_bool = true;
      if( good_hits.at(iWC).at(iAx).hits.size() > 1 ) unique_hit_bool = false;
    }
  }
  if( missing_hit_bool ) fTrack_Type = 0;
  else if( lonely_hit_bool && unique_hit_bool ) fTrack_Type = 1;
  else if( lonely_hit_bool && !unique_hit_bool ) fTrack_Type = 2;
  else if( !lonely_hit_bool && !unique_hit_bool ) fTrack_Type = 3;
  else{ std::cout << "Unknown track condition. Check me." << std::endl;
  }

  //Now determine if we want to skip
  bool skip = false;
  for( size_t iWC = 0; iWC < good_hits.size() ; ++iWC ){
    for( size_t iAx = 0; iAx < good_hits.at(iWC).size() ; ++iAx ){
      if( fPickyTracks ){
	if( good_hits.at(iWC).at(iAx).hits.size() != 1 ){ //<-----------------------------------------TO CHANGE THE ONE-HIT-PER-PLANE TRACK RESTRICTION!!!!
	  skip = true;
	  break;
	}
      }
      if( !fPickyTracks ){
	if( good_hits.at(iWC).at(iAx).hits.size() < 1 ){ //<-----------------------------------------TO CHANGE THE ONE-HIT-PER-PLANE TRACK RESTRICTION!!!!
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
    //Push back a placeholder vect into the reco_pz_array
    std::vector<double> emptyRecoPVect;
    reco_pz_array.push_back(emptyRecoPVect);
    return true;
  }
  else return false;

}

//=====================================================================
//NOTE: BE CAREFUL ABOUT REFERENCING STRUCTS - MIGHT NOT REFERENCE STUFF INSIDE? CHECK
void WCTrackBuilderAlg::findGoodHits( std::vector<std::vector<float> > cluster_time_buffer,
				      std::vector<std::vector<float> > cluster_wire_buffer,
				      std::vector<std::vector<WCHitList> > & good_hits)
{
  //Loop through wire chamber axes (remember, 0-7)
  for( int iWCAx = 0; iWCAx < fNumber_wire_chambers*2 ; ++iWCAx ){
    
    //Sanity check
    if(fVerbose){
      if( cluster_time_buffer.at(iWCAx).size() != cluster_wire_buffer.at(iWCAx).size() ){
	std::cout << "Cluster wire/time buffer size mismatch! Error!" << std::endl;
	return;
      }
    }    
    size_t number_clusters = cluster_time_buffer.at(iWCAx).size();   
    

    //  std::cout << "number clusters WCAx: " << iWCAx << " is " << number_clusters << std::endl;
 
    //Loop through clusters
    for( size_t iClust = 0; iClust < number_clusters; ++iClust ){
      float time = float(cluster_time_buffer.at(iWCAx).at(iClust));
      float wire = float(cluster_wire_buffer.at(iWCAx).at(iClust));
      //Convert back to wire chamber and axis for good hits. Lack of foresight on my part (as everywhere here)
      int wire_chamber = int(iWCAx/2);
      int axis = iWCAx % 2;
      finalizeGoodHits(wire,time,good_hits.at(wire_chamber).at(axis));
    }
  }
  
  if( fVerbose ){
    std::cout << "Number of good hits in each wire plane: "; 
    for( int iWC = 0; iWC < 4 ; ++iWC ){
      for (int iAx = 0; iAx < 2 ; ++iAx ){
	std::cout << good_hits.at(iWC).at(iAx).hits.size() << ", ";
      }
    }
    std::cout << std::endl;
  }
}

//=====================================================================
//Finalize the hits and place them into the final good hit list
void WCTrackBuilderAlg::finalizeGoodHits(float wire,
					 float time,
					 WCHitList & finalGoodHitList)
{
  //Loop through the existing good hit list
  for( size_t iHit = 0; iHit < finalGoodHitList.hits.size(); ++iHit ){
    WCHit hit = finalGoodHitList.hits.at(iHit);
    //If there are good hits close enough to each other, average them and get rid of the old hit
    if( fabs(hit.wire-wire) < fGoodHitAveragingEps && fabs(hit.time-time) < fGoodHitAveragingEps ){
      float average_wire = (hit.wire+wire)/2;
      float average_time = (hit.time+time)/2;
      finalGoodHitList.hits.erase(finalGoodHitList.hits.begin()+iHit);
      wire = average_wire;
      time = average_time;
    }
  }
  //Now that we have averaged, push back the final good hit list with the (possibly) averaged hit
  WCHit finalGoodHit;
  finalGoodHit.wire = wire;
  finalGoodHit.time = time;
  finalGoodHitList.hits.push_back(finalGoodHit);
}




//=====================================================================
//Take in tdc and channel and convert to wire number and axis
//UPDATED
void WCTrackBuilderAlg::convertToWireNumber(int channel,
			 int TDC_index,
			 float & wire)
{
  wire = ((TDC_index+1) % 2) *64 - channel;                  //Want wire number to be zero in middle
}

//=====================================================================
//Use DBSCAN to create clusters of hits that are spatially and temporally close
//UPDATED
void WCTrackBuilderAlg::createClusters( int trigger_number,
					std::vector<std::vector<float> > hit_time_buffer,
					std::vector<std::vector<float> > hit_wire_buffer,
					std::vector<std::vector<float> > & cluster_time_buffer,
					std::vector<std::vector<float> > & cluster_wire_buffer)
{

  //Parameter for clusters being too large
  size_t max_hits = 10;

  //Loop through WC Axes
  for( int iWCAx = 0; iWCAx < fNumber_wire_chambers*2; ++iWCAx ){

    //Create hit and scaled hit vectors for use in DBSCAN's clustering
    WCHitList hits;
    WCHitList scaled_hits;
    createHitAndScaledHitVectors( iWCAx, hit_time_buffer, hit_wire_buffer, scaled_hits ); //(fill the above vectors with scaled/hits)
   
    std::vector<WCHitList> cluster_list;
    if( scaled_hits.hits.size() != 0 )
      run_DBSCAN( trigger_number, iWCAx, scaled_hits,cluster_list);
    
    //Loop through clusters and see if they are too big (pancake-like cross-talk)
    for( size_t iClust = 0; iClust < cluster_list.size(); ++iClust ){
      if( cluster_list.at(iClust).hits.size() > max_hits ){ continue; }
      float wire = 9999;
      float time = 9998;
      findLowestTimeHitInCluster( cluster_list.at(iClust), wire, time );
      cluster_time_buffer.at(iWCAx).push_back(time/fTime_bin_scaling);
      cluster_wire_buffer.at(iWCAx).push_back(wire/fWire_scaling);
    }
  }
}

//=====================================================================
//Finding the one hit with the first time to represent the true
//passing point of the particle
void WCTrackBuilderAlg::findLowestTimeHitInCluster( WCHitList cluster,
				 float & wire,
				 float & time )
{
  float lowest_time = 9997;
  float lowest_time_index = 9996;
  for( size_t iHit = 0; iHit < cluster.hits.size() ; ++iHit ){
    if( cluster.hits.at(iHit).time < lowest_time ){
      lowest_time = cluster.hits.at(iHit).time;
      lowest_time_index = iHit;
    }
  }
  wire = cluster.hits.at(lowest_time_index).wire;
  time = cluster.hits.at(lowest_time_index).time;
}


//=====================================================================
//function that runs dbscan algorithm on the hit wire/time set
//Input is ^
//Output is vector of same length as # of hits with labels
void WCTrackBuilderAlg::run_DBSCAN( int trigger_number,
				    int WCAx_number,
				    WCHitList scaled_hits, 
				    std::vector<WCHitList> & cluster_list )
{
  //Parameters for algorithim
  float epsilon = fDBSCANEpsilon;
  size_t min_hits = fDBSCANMinHits;

  //Create matrix of neighborhoods of hits, given epsilon
  std::vector<WCHitList> neighborhood_matrix = createNeighborhoodMatrix( scaled_hits, epsilon );
  
  //Cluster counter
  int cluster_counter = 0;
  for( size_t iSH = 0; iSH < scaled_hits.hits.size(); ++iSH ){
    //If hit has been visited, ignore it
    if( scaled_hits.hits.at(iSH).isVisited == true ){ continue; }
    scaled_hits.hits.at(iSH).isVisited = true;
    WCHitList neighbor_hits = regionQuery( neighborhood_matrix, scaled_hits.hits.at(iSH), scaled_hits, epsilon );
    
    //If there aren't enough nearest neighbors, count this hit as noise for now
    if( neighbor_hits.hits.size() < min_hits ){ scaled_hits.hits.at(iSH).cluster_index = -1; }
    
    //Otherwise, create a new cluster from it
    else{ 
      cluster_counter++;
      scaled_hits.hits.at(iSH).cluster_index = cluster_counter-1;
      expandCluster( neighborhood_matrix, scaled_hits.hits.at(iSH), scaled_hits, neighbor_hits, epsilon, min_hits, cluster_counter-1 );
      
    }
  }
  
 
   //Extract cluster info from scaled_hits and put into clusters
   for( int iClust = 0; iClust < cluster_counter; ++iClust ){
     WCHitList cluster;
     cluster_list.push_back(cluster);
   }
   for( size_t iSH = 0; iSH < scaled_hits.hits.size(); ++iSH ){
     if( scaled_hits.hits.at(iSH).cluster_index == -1 ){ 
       continue;
     }   
     cluster_list.at(scaled_hits.hits.at(iSH).cluster_index).hits.push_back(scaled_hits.hits.at(iSH));
   }
}

//=====================================================================
//Creating a matrix of neighbors for use in efficient DBSCAN operation
std::vector<WCHitList> WCTrackBuilderAlg::createNeighborhoodMatrix( WCHitList scaled_hits,
						 float epsilon )
{
  //Final vector
  std::vector<WCHitList> neighborhoods_vector;

  //Loop through all hits
  for( size_t iHit = 0; iHit < scaled_hits.hits.size() ; ++iHit ){
    //Create a hit list representing all hits (including itself) that are within epsilon of that hit
    WCHitList neighborhood_hits;
    for( size_t iSubHit = 0; iSubHit < scaled_hits.hits.size(); ++iSubHit ){
      float distance = pow(pow(scaled_hits.hits.at(iHit).wire-scaled_hits.hits.at(iSubHit).wire,2) + pow(scaled_hits.hits.at(iHit).time-scaled_hits.hits.at(iSubHit).time,2),0.5);
      if( distance < epsilon ){ neighborhood_hits.hits.push_back(scaled_hits.hits.at(iSubHit)); }
    }
    neighborhoods_vector.push_back(neighborhood_hits);
  }
  return neighborhoods_vector;
}

//=====================================================================
//DBSCAN function - once a cluster is seeded, reach out from it and find
//other hits that are also in this cluster
void WCTrackBuilderAlg::expandCluster( std::vector<WCHitList> neighborhood_matrix,
		    WCHit the_hit,
		    WCHitList & scaled_hits,
		    WCHitList neighbor_hits,
		    float epsilon,
		    size_t min_hits,
		    int cluster_index)
{

   //Loop through neighbors
  for( size_t iNB = 0; iNB < neighbor_hits.hits.size() ; ++iNB ){
    if( neighbor_hits.hits.at(iNB).isVisited == false ){
      //Need to do each setting for the neighbor hit and the scaled hits cluster (scaled hits retains all info)
      neighbor_hits.hits.at(iNB).isVisited = true;
      scaled_hits.hits.at(neighbor_hits.hits.at(iNB).hit_index).isVisited = true;
      WCHitList next_neighbors = regionQuery( neighborhood_matrix, neighbor_hits.hits.at(iNB), scaled_hits, epsilon );
      //If there are enough next-neighbors for this neighbor, append the neighbor hits list
      if( next_neighbors.hits.size() >= min_hits ){
	for( size_t iNN = 0; iNN < next_neighbors.hits.size() ; ++iNN ){
	  //If the hit already exists in the neighbors vector, then continue
	  bool isAlreadyFound = false;
	  for( size_t iHit2 = 0; iHit2 < neighbor_hits.hits.size() ; ++iHit2 ){
	    if( next_neighbors.hits.at(iNN).hit_index == neighbor_hits.hits.at(iHit2).hit_index ){
	      isAlreadyFound = true;
	    }
	  }
	  if( isAlreadyFound == false ){
	    // if( scaled_hits.hits.at(next_neighbors.hits.at(iNN).hit_index).isVisited == false ){
	    neighbor_hits.hits.push_back(next_neighbors.hits.at(iNN));
	  }
	}
      }
    }
    
    //    std::cout << "Scaled hits size: " << scaled_hits.hits.size() << ", searched hit index: " << neighbor_hits.hits.at(iNB).hit_index << std::endl;
    
    //If this neighbor hit is not yet part of a cluster, add it
    if( neighbor_hits.hits.at(iNB).cluster_index == -1 ){ 
      neighbor_hits.hits.at(iNB).cluster_index = cluster_index;
      scaled_hits.hits.at(neighbor_hits.hits.at(iNB).hit_index).cluster_index = cluster_index;
    }
    
  }
}


 
//=====================================================================
//DBSCAN function - finds hits within epsilon of the central hit
WCHitList WCTrackBuilderAlg::regionQuery( std::vector<WCHitList> neighborhood_matrix,
					  WCHit the_hit,
					  WCHitList & scaled_hits,
					  float epsilon )
{
  //Create a final hit list to return
  WCHitList neighbor_hit_list;

  for( size_t iNH = 0; iNH < neighborhood_matrix.at(the_hit.hit_index).hits.size(); ++iNH ){
    //    if( scaled_hits.hits.at(neighborhood_matrix.at(the_hit.hit_index).hits.at(iNH).hit_index).isVisited == true ){ continue; }
    neighbor_hit_list.hits.push_back(scaled_hits.hits.at(neighborhood_matrix.at(the_hit.hit_index).hits.at(iNH).hit_index));
  }
  return neighbor_hit_list;
}


//=====================================================================
//Create hit vectors for convenient use in DBSCAN's clustering
void WCTrackBuilderAlg::createHitAndScaledHitVectors( int WCAx_number,
						      const std::vector<std::vector<float> > hit_time_buffer,
						      const std::vector<std::vector<float> > hit_wire_buffer,
						      WCHitList & scaled_hit_list)
{
  //Sanity Check
  if( fVerbose ){
    if( hit_time_buffer.size() != hit_wire_buffer.size() ){ std::cout << "Error: vector size mismatch." << std::endl; }
    if( hit_time_buffer.at(WCAx_number).size() != hit_wire_buffer.at(WCAx_number).size() ){ std::cout << "Error: sub-vector size mismatch." << std::endl; }
  }

  //For each element in the hit time buffer (for each hit)
  for( size_t iHit = 0; iHit < hit_time_buffer.at(WCAx_number).size(); ++iHit ){
    //Create a scaled hit and fill it with time, wire number, and whether it has been visited
    //Also with a hit index and cluster index (-1 means noise cluster)
    WCHit scaled_hit;
    scaled_hit.time = hit_time_buffer.at(WCAx_number).at(iHit)*fTime_bin_scaling;
    scaled_hit.wire = hit_wire_buffer.at(WCAx_number).at(iHit)*fWire_scaling;
    scaled_hit.isVisited = false;
    scaled_hit.hit_index = iHit;
    scaled_hit.cluster_index = -1;
    scaled_hit_list.hits.push_back(scaled_hit);
  }
}

//=====================================================================
//Fill buffers that are indexed by wire chamber # (0-7)
//Here, 0 is WC1X, 1 is WC1Y, 2 is WC2X, and so on...
//This standard is used throughout this file
//UPDATED
void WCTrackBuilderAlg::fillTimeAndWireBuffers( const std::vector<int> & tdc_number_vect,
						std::vector<std::vector<float> > & hit_time_buffer,
						std::vector<std::vector<float> > & hit_wire_buffer,
						const std::vector<float> & hit_time_vect,
						const std::vector<float> & hit_channel_vect)
{
  //Sanity check
  //if( fVerbose ){
  //  std::cout << "*************** Buffer Filling Info *****************" << std::endl;
  // }

  //Loop over the wire plane axes
  for( int iWCAx = 0; iWCAx < fNumber_wire_chambers*2 ; ++iWCAx ){
    //Loop over the tdc labels (possibly repeated) within this trigger and find those with the WCAxis = iWCAx
    for( size_t iTDC = 0; iTDC < tdc_number_vect.size(); ++iTDC ){

      int hit_wire_chamber_axis = int((tdc_number_vect.at(iTDC)-1)/2); //-1 is for tdc index
      //      std::cout << "TEST: tdc: " << tdc_number_vect.at(iTDC)-1 << ", WCAx: " << hit_wire_chamber_axis << ", iWCAx: " << iWCAx << std::endl;
      if( hit_wire_chamber_axis == iWCAx ){
	float wire = 0;
	convertToWireNumber( hit_channel_vect.at(iTDC), tdc_number_vect.at(iTDC)-1, wire );
	hit_wire_buffer.at(iWCAx).push_back(wire);
	hit_time_buffer.at(iWCAx).push_back(float(hit_time_vect.at(iTDC)));

	//Sanity Check
	//	if( fVerbose ){
	//  std::cout << "(iTDC,Channel): (" << tdc_number_vect.at(iTDC)-1 << "," << hit_channel_vect.at(trigger_number).at(iTDC) << "), (WCAx,Wire): (" << iWCAx << "," << wire << ")" << ", time: " << hit_time_vect.at(trigger_number).at(iTDC) << std::endl;
	//	}
      }
    }
  }
}	       

//=====================================================================
//Set the length of the time and wire buffers to hold 8 arrays (the number of the wire chamber axes)
//UPDATED
void WCTrackBuilderAlg::initializeBuffers( std::vector<std::vector<float> > & hit_time_buffer,
					   std::vector<std::vector<float> > & hit_wire_buffer,
					   std::vector<std::vector<float> > & cluster_time_buffer,
					   std::vector<std::vector<float> > & cluster_wire_buffer ) 
{
  std::vector<float> temp_buffer;
  for( int iWCAx = 0; iWCAx < fNumber_wire_chambers*2 ; ++iWCAx ){
    hit_time_buffer.push_back(temp_buffer);
    hit_wire_buffer.push_back(temp_buffer);
    cluster_time_buffer.push_back(temp_buffer);
    cluster_wire_buffer.push_back(temp_buffer);
  }
  
  //Sanity check
  //  if( fVerbose == true ){
  //  std::cout << "Length of the hit time buffer: " << hit_time_buffer.size() << std::endl;
  //}

}




//=====================================================================
void WCTrackBuilderAlg::findTrackOnTPCInfo(WCHitList track, float &x, float &y, float &theta, float &phi )
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
void WCTrackBuilderAlg::transformWCHits( float (&WC3_point)[3],
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
bool WCTrackBuilderAlg::cutOnGoodTracks( WCHitList track,
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
void WCTrackBuilderAlg::disambiguateTracks( std::vector<double> & reco_pz_list,
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
