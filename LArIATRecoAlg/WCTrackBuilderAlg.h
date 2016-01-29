
#ifndef WCTRACKBUILDERALG_H
#define WCTRACKBUILDERALG_H

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
class WCTrackBuilderAlg{
 public:
  
  //Constructor/destructor
  WCTrackBuilderAlg( fhicl::ParameterSet const& pset );
  ~WCTrackBuilderAlg();
  
   void reconfigure( fhicl::ParameterSet const& pset );
   
   void loadXMLDatabaseTableForBField( int run, int subrun );
   
   int reconstructTracks(std::vector<double> & reco_pz_list,               
					     std::vector<double> & x_face_list,
					     std::vector<double> & y_face_list,
					     std::vector<double> & incoming_theta_list,
					     std::vector<double> & incoming_phi_list,
					     std::vector<WCHitList> & event_final_tracks,
					     std::vector<std::vector<WCHitList> > & good_hits,
					     bool pickytracks,
					     bool highyield,
					     std::vector<double> & y_kink_list,
					     std::vector<double> & x_dist_list,
					     std::vector<double> & y_dist_list,
					     std::vector<double> & z_dist_list);
		
   bool shouldSkipTrigger(std::vector<std::vector<WCHitList> > & good_hits);
   
   void buildFourPointTracks(std::vector<std::vector<WCHitList> > & good_hits,
	                      std::vector<double> & reco_pz_list,
			      std::vector<double> & x_face_list,
		              std::vector<double> & y_face_list,
			      std::vector<double> & incoming_theta_list,
			      std::vector<double> & incoming_phi_list,
			      std::vector<WCHitList> & event_final_tracks);
			      
   void findTheHitPositions(WCHitList & track,
		            float (&x)[4],
         	            float (&y)[4],
		            float (&z)[4]);
		       
   std::vector<float> Regression(float (&y)[4],
			         float (&z)[4]);
				 
   void calculateTheMomentum(WCHitList & best_track,
		             float (&x)[4],
			     float (&y)[4],
			     float (&z)[4],
			     float & reco_pz,
			     std::vector<float> & BestTrackStats);				 
				 
   void projectToTPC(WCHitList & best_track,
		    float (&x)[4],
		    float (&y)[4],
		    float (&z)[4],
		    std::vector<float> & bestRegressionStats,
		    std::vector<double> & x_face_list,
		    std::vector<double> & y_face_list,
		    std::vector<double> & incoming_theta_list,
		    std::vector<double> & incoming_phi_list);
		    
		    				 
   void buildThreePointTracks(std::vector<std::vector<WCHitList> > & good_hits,
	                      std::vector<double> & reco_pz_list,
			      std::vector<double> & x_face_list,
		              std::vector<double> & y_face_list,
			      std::vector<double> & incoming_theta_list,
			      std::vector<double> & incoming_phi_list,
			      std::vector<WCHitList> & event_final_tracks);
			      
   void calculateTheThreePointMomentum(WCHitList & best_track,
				       float(&x)[4],
				       float(&y)[4],
				       float(&z)[4],
				       float & reco_pz,
				       std::vector<float> & BestTrackStats);
				       				 
				 
  void extrapolateTheMissedPoint(WCHitList & best_track,
			         float(&x)[4],
			         float(&y)[4],
			         float(&z)[4],
				 float & reco_pz,
				 std::vector<float> & BestTrackStats,
				 std::vector<float> & missed_wires);
		
  private:
  
  bool   fPickyTracks;
  bool   fHighYield;
  int NHits;
  int WCMissed;
  				 
  art::ServiceHandle<util::DatabaseUtilityT1034> fDatabaseUtility;
 
 
  /////////////////////////////////
  // CONSTANTS FROM SURVEY       //
  /////////////////////////////////
				 

  float fL_eff;      /// the effective magnetic length of the magnet
  float fmm_to_m;
  float fGeV_to_MeV;				 

  //center of the WC's
  float fX_cntr[4];
  float fY_cntr[4];
  float fZ_cntr[4];
  
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
  float fMidplane_intercept;
  float fMidplane_slope;
  int initialconst;
  art::ServiceHandle<geo::Geometry> fGeo;  /// Geometry Service Handle				 
				 
};

#endif				 
				 
				 
				 
				 
				 
				 
				 
				 
				 
				 
				 
				 
				 
				 
				 
