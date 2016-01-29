
#ifndef WCTRACKBUILDERALG_NEW_H
#define WCTRACKBUILDERALG_NEW_H

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
#include <TH2F.h>

//Structs for organizational purposes pulling from the HitFinderAlg
#include "LArIATRecoAlg/WCHitFinderAlg.h"



//--------------------------------------------
class WCTrackBuilderAlg_new{
 public:
  
  //Constructor/destructor
  WCTrackBuilderAlg_new( fhicl::ParameterSet const& pset );
  ~WCTrackBuilderAlg_new();
  
  
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
			  bool pickytracks,
			  bool highyield,
			  int & track_count,/*, 
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
			  std::vector<TH1F*> & RegressionPlots1D,*/
			  std::vector<TH2F*> & Recoplots/*
			  TH1F* & Bfield */);


  void getTrackMom_Kink_End(WCHitList track,
			    float & reco_pz,
			    float & y_kink,
			    float (&dist_array)[3],/* ,
			    TH2F* & TargetXY,
			    TH1F* & ResSquare,
			    TH1F* & Reco4pt,
			    TH2F* & Reco4ptdiff,
			    std::vector<TH2F*> & RegressionPlots,
			    std::vector<TH1F*> & RegressionPlots1D,*/
			    std::vector<TH2F*> & Recoplots );

  void midPlaneExtrapolation(std::vector<float> x_wires,
			     std::vector<float> y_wires,
			     float (&pos_us)[3],
			     float (&pos_ds)[3],
			     double reco_pz,
			     int WCMissed,/*,
			     TH2F* & TargetXY,
			     TH1F* & ResSquare,
 			     TH1F* & Reco4pt,
			     TH2F* & Reco4ptdiff,
			     std::vector<TH2F*> & RegressionPlots,
 			     std::vector<TH1F*> & RegressionPlots1D,*/
 			     std::vector<TH2F*> & Recoplots);
			     
  std::vector<float> Regression(float (&y)[4],
                            float (&z)[4],/*
 			    TH1F* & ResSquare */
			    int skippedWC);
  
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
			   std::vector<WCHitList> & track_list,/* ,
			   std::vector<TH1F*> & WCHitsGoodTracks,
			   std::vector<TH2F*> & WCMult,
			   std::vector<TH1F*> & BadTrackHits,
			   TH2F* & TargetXY,
			   TH2F* & PickyTracksXY,
			   TH1F* & ResSquare,
			   TH1F* & Reco4pt,
			   TH2F* & Reco4ptdiff,
			   std::vector<TH2F*> & TimingXY,
			   std::vector<TH2F*> & RegressionPlots,
			   std::vector<TH1F*> & RegressionPlots1D,*/
			   std::vector<TH2F*> & Recoplots);

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
  bool   fHighYield;
  int NHits;
  int WCMissed;
  std::vector<float> BestYStats;
  std::vector<float> BestYHits;




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
  float current;
  float fMP_X;
  float fMP_M;
  art::ServiceHandle<geo::Geometry> fGeo;  /// Geometry Service Handle
  
};


#endif
