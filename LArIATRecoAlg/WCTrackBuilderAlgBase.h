
#ifndef WCTRACKBUILDERALGBASE_H
#define WCTRACKBUILDERALGBASE_H

//C++ includes
#include <vector>
#include <cmath>
#include <iostream>

//Framework includes
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"


// LArSoft includes
#include "Geometry/Geometry.h"

//LArIAT include
#include "Utilities/DatabaseUtilityT1034.h"

//ROOT includes
#include <TH1F.h>

//Structs for organizational purposes pulling from the HitFinderAlg
#include "LArIATRecoAlg/WCHitFinderAlg.h"



//--------------------------------------------
class WCTrackBuilderAlgBase{
 public:
  
  //Constructor/destructor
  WCTrackBuilderAlgBase( fhicl::ParameterSet const& pset );
  ~WCTrackBuilderAlgBase();
  
  
  void reconfigure( fhicl::ParameterSet const& pset );
  

  
  void reconstructTracks( std::vector<double> & reco_pz_list,               
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
			  int & track_count);


  void getTrackMom_Kink_End(WCHitList track,
			    float & reco_pz,
			    float & y_kink,
			    float (&dist_array)[3]);
  
  void midPlaneExtrapolation(std::vector<float> x_wires,
			     std::vector<float> y_wires,
			     float (&pos_us)[3],
			     float (&pos_ds)[3]);
  
  bool buildTracksFromHits(std::vector<std::vector<WCHitList> > & good_hits,
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
			   std::vector<WCHitList> & track_list);

  bool shouldSkipTrigger(std::vector<std::vector<WCHitList> > & good_hits);

  
  void findTrackOnTPCInfo(WCHitList track,
			  float &x,
			  float &y,
			  float &theta,
			  float &phi);
  
  void transformWCHits( float (&WC3_point)[3],
			float (&WC4_point)[3]);

  bool cutOnGoodTracks( WCHitList track,
			float & y_kink,
			float (&dist_array)[3],
			size_t track_index);

  void disambiguateTracks( std::vector<double> & reco_pz_list,
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
			   bool lonely_hit_bool);

  void loadXMLDatabaseTableForBField( int run, int subrun );

 private:
  float  fCentralYKink;
  float  fSigmaYKink;
  float  fCentralYDist;
  float  fSigmaYDist;
  float  fPrintDisambiguation;
  bool   fPickyTracks;


  art::ServiceHandle<util::DatabaseUtilityT1034> fDatabaseUtility;


  //Semi-Persistent vectors
  std::map<size_t,std::pair<float,float> > fGoodTrackCandidateErrors;
  std::vector<std::pair<WCHitList,size_t> > fGoodTrackCandidateHitLists;

  /////////////////////////////////
  // CONSTANTS FROM SURVEY       //
  /////////////////////////////////

  // constants and survey data for position reconstruction
  float fDelta_z_us; /// millimeters
  float fDelta_z_ds; /// millimeters
  float fL_eff;      /// the effective magnetic length of the
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

  double fCenter_of_tpc[3];  //<------------------- CENTER OF TPC HERE !!!!!!!!
  float fHalf_z_length_of_tpc;
  float fHalf_x_length_of_tpc;


  //Misc
  float fB_field_tesla;
  bool fVerbose;
  int fRun;
  int fSubRun;

  art::ServiceHandle<geo::Geometry> fGeo;  /// Geometry Service Handle
  
};


#endif
