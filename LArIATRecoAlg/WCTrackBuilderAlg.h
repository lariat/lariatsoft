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


#ifndef WCTRACKBUILDERALG_H
#define WCTRACKBUILDERALG_H

//C++ includes
#include <vector>
#include <cmath>
#include <iostream>

//Framework includes
#include "fhiclcpp/ParameterSet.h"
#include "art/Persistency/Common/Ptr.h"
#include "art/Persistency/Common/PtrVector.h"

//Structs for organizational purposes
struct WCHit
{
  float wire;
  float time;
  bool isVisited;
  int hit_index;
  int cluster_index;
};

struct WCHitList
{
  std::vector<WCHit> hits;
};




//--------------------------------------------
class WCTrackBuilderAlg{
 public:
  
  //Constructor/destructor
  WCTrackBuilderAlg( fhicl::ParameterSet const& pset );
  ~WCTrackBuilderAlg();
  
  
  void reconfigure( fhicl::ParameterSet const& pset );
  
  void firstFunction();
  
  
  void reconstructTracks( std::vector<int> tdc_number_vect,
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
					     int track_count);


  void getTrackMom_Kink_End(WCHitList track,
					       float & reco_pz,
					       float & y_kink,
					       float (&dist_array)[3]);
  
  void midPlaneExtrapolation(std::vector<float> x_wires,
						std::vector<float> y_wires,
						float (&pos_us)[3],
						float (&pos_ds)[3]);
  
void buildTracksFromHits(std::vector<std::vector<WCHitList> > & good_hits,
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
					    std::vector<double> & incoming_phi_list);

  bool shouldSkipTrigger(std::vector<std::vector<WCHitList> > & good_hits,
					      std::vector<std::vector<double> > & reco_pz_array);

  void findGoodHits( std::vector<std::vector<float> > cluster_time_buffer,
					std::vector<std::vector<float> > cluster_wire_buffer,
					std::vector<std::vector<WCHitList> > & good_hits);
  
  void finalizeGoodHits(float wire,
					   float time,
					   WCHitList & finalGoodHitList);

  void convertToWireNumber(int channel,
					      int TDC_index,
					      float & wire);

  void createClusters( int trigger_number,
					  std::vector<std::vector<float> > hit_time_buffer,
					  std::vector<std::vector<float> > hit_wire_buffer,
					  std::vector<std::vector<float> > & cluster_time_buffer,
					  std::vector<std::vector<float> > & cluster_wire_buffer);

  void findLowestTimeHitInCluster( WCHitList cluster,
						      float & wire,
						      float & time );

  void run_DBSCAN( int trigger_number,
				      int WCAx_number,
				      WCHitList scaled_hits, 
				      std::vector<WCHitList> & cluster_list );

  std::vector<WCHitList> createNeighborhoodMatrix( WCHitList scaled_hits,
								      float epsilon );
  /*
  void plotHitsAndClustersInRootFile( int trigger_number,
							 int WCAx_number,
							 WCHitList scaled_hits,
							 std::vector<WCHitList> cluster_list );
  */
  void expandCluster( std::vector<WCHitList> neighborhood_matrix,
					 WCHit the_hit,
					 WCHitList & scaled_hits,
					 WCHitList neighbor_hits,
					 float epsilon,
					 size_t min_hits,
					 int cluster_index);

WCHitList regionQuery( std::vector<WCHitList> neighborhood_matrix,
					  WCHit the_hit,
					  WCHitList & scaled_hits,
					  float epsilon );

void createHitAndScaledHitVectors( int WCAx_number,
						      const std::vector<std::vector<float> > hit_time_buffer,
						      const std::vector<std::vector<float> > hit_wire_buffer,
						      WCHitList & scaled_hit_list);

void fillTimeAndWireBuffers( const std::vector<int> & tdc_number_vect,
			     std::vector<std::vector<float> > & hit_time_buffer,
			     std::vector<std::vector<float> > & hit_wire_buffer,
			     const std::vector<float> & hit_time_vect,
			     const std::vector<float> & hit_channel_vect);

void initializeBuffers( std::vector<std::vector<float> > & hit_time_buffer,
					   std::vector<std::vector<float> > & hit_wire_buffer,
					   std::vector<std::vector<float> > & cluster_time_buffer,
					   std::vector<std::vector<float> > & cluster_wire_buffer );

 void findTrackOnTPCInfo(WCHitList track,
			 float &x,
			 float &y,
			 float &theta,
			 float &phi);
 
 void transformWCHits( float (&WC3_point)[3],
		       float (&WC4_point)[3]);
 
  



  
  
 private:
  //Hardware constants
  int fNumber_tdcs;
  int fNumber_wire_chambers;

  //Parameters that can be varied
  double fTime_bin_scaling;
  double fWire_scaling;
  double fGoodHitAveragingEps;
  
 

  /////////////////////////////////
  // CONSTANTS FROM SURVEY       //
  /////////////////////////////////

  // constants and survey data for position reconstruction
  float fDelta_z_us;
  float fDelta_z_ds;
  float fL_eff;
  float fmm_to_m;
  float fGeV_to_MeV;

  // center (cntr) of multi-wire proportional chambers
  float fX_cntr_1;
  float fY_cntr_1;
  float fZ_cntr_1;
  float fX_cntr_2;
  float fY_cntr_2;
  float fZ_cntr_2;
  float fX_cntr_3;
  float fY_cntr_3;
  float fZ_cntr_3;
  float fX_cntr_4;
  float fY_cntr_4;
  float fZ_cntr_4;

  // derive some useful geometric/al constants
  // upstream (us) leg center line projected onto xz-plane
  float fCntr_slope_xz_us;
  float fCntr_z_int_xz_us_1;
  float fCntr_z_int_xz_us_2;
  float fCntr_z_int_xz_us;

  // downstream (ds) leg center line projected onto xz-plane
  float fCntr_slope_xz_ds;
  float fCntr_z_int_xz_ds_3;
  float fCntr_z_int_xz_ds_4;
  float fCntr_z_int_xz_ds;

  // mid-plane
  // includes intersection point of upstream and downstream
  // center lines above and taken to be at 8 y-degrees from the
  // xy-plane

  
  float fMid_plane_x;
  float fMid_plane_z_us;
  float fMid_plane_z_ds;
  float fMid_plane_z;
  float fMid_plane_slope_xz;
  float fMid_plane_z_int_xz;

  float fCenter_of_tpc[3];
  float fHalf_z_length_of_tpc;


  //Misc
  float fB_field_tesla;
  bool fVerbose;
  
};


#endif
