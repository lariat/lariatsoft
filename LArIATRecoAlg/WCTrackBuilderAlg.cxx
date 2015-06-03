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

#include "WCTrackBuilderAlg.h"

#include <iostream>
#include <cmath>
#include <cstdlib>





//--------------------------------------------------------------
//Constructor
WCTrackBuilderAlg::WCTrackBuilderAlg( fhicl::ParameterSet const& pset )
{
  this->reconfigure(pset);
    
  fNumber_tdcs = 16;
  fNumber_wire_chambers = 4;
  fB_field_tesla = 0.35;

  fTime_bin_scaling = 1.0/1280.0;
  fWire_scaling = 1.0/64.0;  
  fGoodHitAveragingEps = 2.0;

  //Survey constants
  delta_z_us = 1551.15;  // millimeters
  delta_z_ds = 1570.06;  // millimeters
  l_eff = 1145.34706;  // the effective magnetic length of the
  mm_to_m = 0.001;
  GeV_to_MeV = 1000.0;
  
  // center (cntr) of multi-wire proportional chambers
  x_cntr_1 = -392.33;
  y_cntr_1 = 0.0;
  z_cntr_1 = 1683.59;
  x_cntr_2 = -738.35;
  y_cntr_2 = 0.0;
  z_cntr_2 = 3195.65;
  x_cntr_3 = -1019.89;
  y_cntr_3 = 0.0;
  z_cntr_3 = 5183.07;
  x_cntr_4 = -1146.64;
  y_cntr_4 = 0.0;
  z_cntr_4 = 7587.99;
  
  // rough guess beginning 2015.02: moved WC4 33" along the 3deg beam.
  // x = -1146.64 mm = -1102.77 mm - sin(3 deg)*33"
  // z =  7587.99 mm =  6750.94 mm + cos(3 deg)*33"
  // Positions from 2014 runs: (x_cntr_4, y_cntr_4, z_cntr_4) = (-1102.77, 0.0, 6750.94)
  
  // derive some useful geometric/al constants
  
  // upstream (us) leg center line projected onto xz-plane
  cntr_slope_xz_us = (z_cntr_2 - z_cntr_1)/(x_cntr_2 - x_cntr_1);
  cntr_z_int_xz_us_1 = z_cntr_1 - (cntr_slope_xz_us * x_cntr_1);
  cntr_z_int_xz_us_2 = z_cntr_2 - (cntr_slope_xz_us * x_cntr_2);
  cntr_z_int_xz_us = (cntr_z_int_xz_us_1 + cntr_z_int_xz_us_2) / 2.0;
  
  // downstream (ds) leg center line projected onto xz-plane
  cntr_slope_xz_ds = (z_cntr_4 - z_cntr_3)/(x_cntr_4 - x_cntr_3);
  cntr_z_int_xz_ds_3 = z_cntr_3 - (cntr_slope_xz_ds * x_cntr_3);
  cntr_z_int_xz_ds_4 = z_cntr_4 - (cntr_slope_xz_ds * x_cntr_4);
  cntr_z_int_xz_ds = (cntr_z_int_xz_ds_3 + cntr_z_int_xz_ds_4) / 2.0;
  
  // mid-plane
  // includes intersection point of upstream and downstream
  // center lines above and taken to be at 8 y-degrees from the
  // xy-plane
  mid_plane_x = -1.0 * (cntr_z_int_xz_ds - cntr_z_int_xz_us) / (cntr_slope_xz_ds - cntr_slope_xz_us);
  
  mid_plane_z_us = cntr_z_int_xz_us + cntr_slope_xz_us * mid_plane_x;
  mid_plane_z_ds = cntr_z_int_xz_ds + cntr_slope_xz_ds * mid_plane_x;
  mid_plane_z = (mid_plane_z_us + mid_plane_z_ds) / 2.0;
  mid_plane_slope_xz = tan(8.0*3.141592654/180);
  mid_plane_z_int_xz = mid_plane_z - mid_plane_slope_xz * mid_plane_x;

  center_of_tpc[0] = -1200;  //First appx
  center_of_tpc[1] = 0;
  center_of_tpc[2] = 8500;
  half_length_of_tpc = 450; //mm	      
  euler_phi = -90.0*3.141592654/180; //ra
  euler_theta = -3*3.141592654/180; //rad
  euler_psi = 90.0*3.141592654/180; //rad
}

//--------------------------------------------------------------  
//Destructor
WCTrackBuilderAlg::~WCTrackBuilderAlg()
{

}

//--------------------------------------------------------------
void WCTrackBuilderAlg::reconfigure( fhicl::ParameterSet const& pset )
{
  
}

//--------------------------------------------------------------
void WCTrackBuilderAlg::firstFunction()
{

  std::cout << "First function called." << std::endl;

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
					   std::vector<std::vector<WCHitList> > & good_hits,
					   bool verbose,
					   int & good_trigger_counter,
					   int trigger_number,
					   int track_count)
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
  
  //GOOD TO HERE
  
  //Create a vector of clusters with DBSCAN, where each entry into cluster_time_buffer/cluster_wire_buffer is a different cluster
  createClusters( trigger_number, hit_time_buffer, hit_wire_buffer, cluster_time_buffer, cluster_wire_buffer );
  
  //GOOD TO HERE, MOST LIKELY
  
  //Finding the hits to be used in momentum reconstruction
  //Note here that only one hit may be accepted from a cluster. The idea is that a particle passing through the MWPC causes noise
  //that is spatially and temporally clustered around the initial hit of the particle. 
   //In addition, if sufficiently close together, two or more good hits can be averaged
  
  findGoodHits(cluster_time_buffer,cluster_wire_buffer,good_hits);
  
  //Determine if one should skip this trigger based on whether there is at least one good hit in each wire chamber and axis
  //If there isn't, continue after adding a new empty vector to the reco_pz_array contianer
  //If there is, then move on with the function 
  
  //Sanity check
  std::cout << "Good hit check: ";
  for( int iWC = 0; iWC < 4 ; ++iWC ){
    for (int iAx = 0; iAx < 2 ; ++iAx ){
      if( good_hits.at(iWC).at(iAx).hits.size() > 0 ) std::cout << good_hits.at(iWC).at(iAx).hits.size();
      else std::cout << "0";
    }
  }
  std::cout << std::endl;
  
  
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
  //of the kinks or end displacements (for now)



  buildTracksFromHits(good_hits,
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
		      incoming_phi_list);

  

  
  


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
  
  //print out tracks
  //  for( size_t iHit = 0; iHit < track.hits.size() ; ++iHit ){
  //  std::cout << "Hit: " << iHit << ", Hit wire: " << track.hits.at(iHit).wire << std::endl;
  // }
    
  for( size_t iPlane = 0; iPlane < wire_x.size() ; ++iPlane ){
    //    std::cout << "Wire_x: " << wire_x.at(iPlane) << std::endl;
  }

  float delta_x_us = wire_x[1] - wire_x[0];
  float delta_y_us = wire_y[1] - wire_y[0];
  float delta_x_ds = wire_x[3] - wire_x[2];
  float delta_y_ds = wire_y[3] - wire_y[2];

  float atan_x_us = atan(delta_x_us / delta_z_us);
  float atan_x_ds = atan(delta_x_ds / delta_z_ds);

  float atan_y_us = atan(delta_y_us / delta_z_us);
  float atan_y_ds = atan(delta_y_ds / delta_z_ds);

  

  //  std::cout << "Atan DS/US: " << atan_x_ds << "," << atan_x_us << std::endl;

  
  reco_pz = (fabs(fB_field_tesla) * l_eff * mm_to_m * GeV_to_MeV ) / float(3.3 * ((10.0*3.141592654/180.0) + (atan_x_ds - atan_x_us)));
  // std::cout << "reco_pz inner: " << reco_pz << std::endl;
  y_kink = atan_y_us - atan_y_ds;

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
  float x[4] = { x_cntr_1 + x_wires[0] * float(cos(3.141592654/180*(13.0))),
		   x_cntr_2 + x_wires[1] * float(cos(3.141592654/180*(13.0))),
		   x_cntr_3 + x_wires[2] * float(cos(3.141592654/180*(3.0))),
		   x_cntr_4 + x_wires[3] * float(cos(3.141592654/180*(3.0))) };
  float y[4] = { y_cntr_1 + y_wires[0],
		   y_cntr_2 + y_wires[1],
		   y_cntr_3 + y_wires[2],
		   y_cntr_4 + y_wires[3] };
  float z[4] = { z_cntr_1 + x_wires[0] * float(sin(3.141592654/180*(13.0))),
		   z_cntr_2 + x_wires[1] * float(sin(3.141592654/180*(13.0))),
		   z_cntr_3 + x_wires[2] * float(sin(3.141592654/180*(3.0))),
		   z_cntr_4 + x_wires[3] * float(sin(3.141592654/180*(3.0))) };
  
  // upstream (us) leg, extrapolating forward to mid-plane
  float slope_xz_us = (z[1] - z[0])/(x[1] - x[0]);
  float z_int_xz_us_1 = z[0] - slope_xz_us * x[0];
  float z_int_xz_us_2 = z[1] - slope_xz_us * x[1];
  float z_int_xz_us = (z_int_xz_us_1 + z_int_xz_us_2) / 2.0;

  float x_us = (z_int_xz_us - mid_plane_z_int_xz) / (mid_plane_slope_xz - slope_xz_us);
  float z_us = mid_plane_z_int_xz + mid_plane_slope_xz * x_us;

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
    
  float x_ds = (z_int_xz_ds - mid_plane_z_int_xz) / (mid_plane_slope_xz - slope_xz_ds);
  float z_ds = mid_plane_z_int_xz + mid_plane_slope_xz * x_ds;

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
void WCTrackBuilderAlg::buildTracksFromHits(std::vector<std::vector<WCHitList> > & good_hits,
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
					    std::vector<double> & incoming_phi_list)
					    
					   
{
  //Reconstructed momentum buffer for storing pz for all combinations in this trigger
  std::vector<double> reco_pz_buffer;
  int track_count_this_trigger = 0;

  std::vector<WCHitList> track_list;

  //Loop through all combinations of tracks
  for( size_t iHit0 = 0; iHit0 < good_hits.at(0).at(0).hits.size(); ++iHit0 ){
    for( size_t iHit1 = 0; iHit1 < good_hits.at(0).at(1).hits.size(); ++iHit1 ){
      for( size_t iHit2 = 0; iHit2 < good_hits.at(1).at(0).hits.size(); ++iHit2 ){
	for( size_t iHit3 = 0; iHit3 < good_hits.at(1).at(1).hits.size(); ++iHit3 ){
	  for( size_t iHit4 = 0; iHit4 < good_hits.at(2).at(0).hits.size(); ++iHit4 ){
	    for( size_t iHit5 = 0; iHit5 < good_hits.at(2).at(1).hits.size(); ++iHit5 ){
	      for( size_t iHit6 = 0; iHit6 < good_hits.at(3).at(0).hits.size(); ++iHit6 ){
		for( size_t iHit7 = 0; iHit7 < good_hits.at(3).at(1).hits.size(); ++iHit7 ){

		  
		  WCHitList track;
		  track.hits.push_back(good_hits.at(0).at(0).hits.at(iHit0));
		  track.hits.push_back(good_hits.at(0).at(1).hits.at(iHit1));
		  track.hits.push_back(good_hits.at(1).at(0).hits.at(iHit2));
		  track.hits.push_back(good_hits.at(1).at(1).hits.at(iHit3));
		  track.hits.push_back(good_hits.at(2).at(0).hits.at(iHit4));
		  track.hits.push_back(good_hits.at(2).at(1).hits.at(iHit5));
		  track.hits.push_back(good_hits.at(3).at(0).hits.at(iHit6));
		  track.hits.push_back(good_hits.at(3).at(1).hits.at(iHit7));
		  track_list.push_back(track);

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
		  
		  //Storing the momentum in the buffer that will be
		  //pushed back into the final reco_pz_array for this trigger
		  reco_pz_buffer.push_back(reco_pz);
		  
		  //Filling full info lists
		  reco_pz_list.push_back(reco_pz);
		  y_kink_list.push_back(y_kink*180.0/3.1415926);
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

  //Fill the array with the found tracks' reconstructed momenta
  reco_pz_array.push_back(reco_pz_buffer);
  
  //Clear the hit lists for each WC/axis
  for( size_t iWC = 0; iWC < good_hits.size() ; ++iWC ){
    for( size_t iAx = 0; iAx < good_hits.at(iWC).size() ; ++iAx ){
      good_hits.at(iWC).at(iAx).hits.clear();
    }
  }	     
}

//=====================================================================
//See if trigger has a good enough hit set to continue
bool WCTrackBuilderAlg::shouldSkipTrigger(std::vector<std::vector<WCHitList> > & good_hits,
					    std::vector<std::vector<double> > & reco_pz_array)
{
  bool skip = false;
  for( size_t iWC = 0; iWC < good_hits.size() ; ++iWC ){
    for( size_t iAx = 0; iAx < good_hits.at(iWC).size() ; ++iAx ){
      if( good_hits.at(iWC).at(iAx).hits.size() != 1 ){
	skip = true;
	break;
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
  }
  for( int iWC = 0; iWC < 4 ; ++iWC ){
    for (int iAx = 0; iAx < 2 ; ++iAx ){
      std::cout << good_hits.at(iWC).at(iAx).hits.size() << ", ";
    }
  }
  std::cout << std::endl;
  
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
    
    //SHOULD BE GOOD UP TO THIS POINT
    
 
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
  float epsilon = 1.0/16.0;
  size_t min_hits = 1;

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
       //   std::cout << "Hit belongs to noise cluster." << std::endl;
       continue;
     }   
     cluster_list.at(scaled_hits.hits.at(iSH).cluster_index).hits.push_back(scaled_hits.hits.at(iSH));
   }
  //  std::cout << "number of clusters: " << cluster_list.size() << ", number of hits: " << scaled_hits.hits.size() << std::endl;

  //Plot hits and clusters in root (happens for each trigger)
   //plotHitsAndClustersInRootFile(trigger_number,WCAx_number,scaled_hits,cluster_list);
  
  

				 
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


/*
//=====================================================================
//Plotting the hits and clusters formed in a root file
void WCTrackBuilderAlg::plotHitsAndClustersInRootFile( int trigger_number,
				    int WCAx_number,
				    WCHitList scaled_hits,
				    std::vector<WCHitList> cluster_list )
{
  //Open a new root file
  TFile * hit_cluster_file = new TFile("/lariat/app/users/linehan3/Summer_2015/DQM_Transfer/hit_cluster_file.root","UPDATE");
  char name[40];
  char name2[40];
  sprintf(name,"hit_hist_Trig%d_WCAx%d",trigger_number,WCAx_number);
  TH2F * hit_2D_hist = new TH2F(name,name,128,-1,1,1280,0,1);
  sprintf(name2,"clust_hist_Trig%d_WCAx%d",trigger_number,WCAx_number);
  TH2F * cluster_2D_hist = new TH2F(name2,name2,128,-1,1,1280,0,1);
  
  //Loop through hits and fill hit histo
  for( size_t iHit = 0; iHit < scaled_hits.hits.size() ; ++iHit ){
    //    std::cout << "Scaled Wire: " << scaled_hits.hits.at(iHit).wire << ", Scaled Time: " << scaled_hits.hits.at(iHit).time << std::endl;
    hit_2D_hist->SetBinContent((scaled_hits.hits.at(iHit).wire*64+64),(scaled_hits.hits.at(iHit).time*1280),1);
  }

  //Loop through clusters and fill cluster histo with different color for each cluster
  for( size_t iClust = 0; iClust < cluster_list.size() ; ++iClust ){
    for( size_t iHit = 0; iHit < cluster_list.at(iClust).hits.size() ; ++iHit ){
      cluster_2D_hist->SetBinContent((cluster_list.at(iClust).hits.at(iHit).wire*64+64),(cluster_list.at(iClust).hits.at(iHit).time*1280),iClust+1);
    }
  }

  hit_2D_hist->Write();
  cluster_2D_hist->Write();
  hit_cluster_file->Write();
  hit_cluster_file->Close();



}
*/

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
  float WC3_point[3] = { x_cntr_3 + track.hits.at(4).wire*float(cos(3.141592654/180*(3.0))),
			   y_cntr_3 + track.hits.at(5).wire,
			   z_cntr_3 + track.hits.at(4).wire*float(sin(3.141592654/180*(3.0))) };
  float WC4_point[3] = { x_cntr_4 + track.hits.at(6).wire*float(cos(3.141592654/180*(3.0))),
			   y_cntr_4 + track.hits.at(7).wire,
			   z_cntr_4 + track.hits.at(6).wire*float(sin(3.141592654/180*(3.0))) };

  std::cout << "Pre: WC3: (" << WC3_point[0] << "," << WC3_point[1] << "," << WC3_point[2] << ")" << " ----> WC4: " <<
    "(" << WC4_point[0] << "," << WC4_point[1] << "," << WC4_point[2] << ")" << std::endl;


  transformWCHits(WC3_point,WC4_point);

  std::cout << "Post: WC3: (" << WC3_point[0] << "," << WC3_point[1] << "," << WC3_point[2] << ")" << " ----> WC4: " <<
    "(" << WC4_point[0] << "," << WC4_point[1] << "," << WC4_point[2] << ")" << std::endl;

  //Now have hit vectors in the frame of the TPC. Now we recreate the second track and find its
  //intersection with the upstream plane of the TPC. In this new frame, the upstream plane is just
  //Z = -450 mm. So we parametrize the track with t and find at which t Z = -450. We then use that
  //to get X and Y intercepts.
  
  float parameter_t = (-1*half_length_of_tpc-WC3_point[2])/(WC4_point[2]-WC3_point[2]);
  float x_at_US_plane = (WC3_point[0])+parameter_t*(WC4_point[0]-WC3_point[0]);
  float y_at_US_plane = (WC3_point[1])+parameter_t*(WC4_point[1]-WC3_point[1]);
  
  x = x_at_US_plane;
  y = y_at_US_plane;
  float r = pow(pow(x-WC4_point[0],2)+pow(y-WC4_point[1],2),0.5);
  std::cout << "r: " << r << ", denom: " << -1*half_length_of_tpc-WC4_point[2] << ", x_face: " << x << ", WC4_point[2]: " << WC4_point[0] << ", WC3Point[0]: " << WC3_point[0] << std::endl;
  theta = atan(r/(-1*half_length_of_tpc-WC4_point[2]));

  //Calculating phi
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
    WC3_point[iDim] = WC3_point[iDim] - center_of_tpc[iDim];
    WC4_point[iDim] = WC4_point[iDim] - center_of_tpc[iDim];
  }
   
}

