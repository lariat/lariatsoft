#include <vector>
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "CLHEP/Units/SystemOfUnits.h"
#include "Geometry/AuxDetGeo.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "Utilities/DatabaseUtilityT1034.h"
class WCTrackAlgBase
{
   public:

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
   int fRun;
   int fSubRun;
   float fB_field_tesla;
   art::ServiceHandle<util::DatabaseUtilityT1034> fDatabaseUtility;
   art::ServiceHandle<geo::Geometry> fGeo; 
   
          virtual void Hello() = 0;
	  void reconfigure( fhicl::ParameterSet const& pset )
	  {  fB_field_tesla        = pset.get<float >("BFieldInTesla",      0.       );
	  }
	  //virtual void JumpThroughHoopFactory()=0;
	   void loadXMLDatabaseTableForBField( int run, int subrun )
{
  fRun = run;
  fSubRun = subrun;
  fB_field_tesla = 0.0035*std::stod(fDatabaseUtility->GetIFBeamValue("mid_f_mc7an",fRun,fSubRun));
  std::cout << "Run: " << fRun << ", Subrun: " << fSubRun << ", B-field: " << fB_field_tesla << std::endl;
}

 
/*    virtual void reconstructTracks( std::vector<double> & reco_pz_list,               
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
	  	  */

   void InitializeGeometry()
   {
   
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
  std::cout<<fHalf_x_length_of_tpc<<std::endl;
  }
 };
typedef WCTrackAlgBase* (*CreateWCTrackBuilderAlgFn)(void);
