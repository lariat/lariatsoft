////////////////////////////////////////////////////////////////////////
// This Alg takes in the good hits from WCHitFinderAlg and creates    //
// a track with various kinematic and geometric values.  In effect,   //
// this alg is the track building part of WCTrackBuilderAlg.cxx       //
// in functions run after "FinalizeGoodHits".                         //
// Other algs can be created in lieu of this one if a better track    //
// is developed.  This is merely the version we used when hit finding //
// and track builder were run together in WCTrackBuilderAlg.cxx       //
// *** All distances in here are in mm! ***                           //
// Author: Greg Pulliam gkpullia@syr.edu                              //
//////////////////////////////////////////////////////////////////////// 

//Framework includes
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "CLHEP/Units/SystemOfUnits.h"
#include "larcorealg/Geometry/AuxDetGeo.h"
#include "art/Framework/Services/Optional/TFileService.h"

// LArIAT includes
#include "LArIATRecoAlg/WCTrackBuilderAlg.h"


#include <iostream>
#include <cmath>
#include <cstdlib>
#include <TH1F.h>
#include <string>
#include <TH2F.h>
#include <TVector3.h>
#include <TMath.h>

WCTrackBuilderAlg::WCTrackBuilderAlg( fhicl::ParameterSet const& pset )
{
  this->reconfigure(pset);
  fB_field_tesla=fMCMagneticField;
  //Testing the AuxDetGeo capabilitites
  double centerOfDet[3] = {0,0,0};
  for( size_t iDet = 0; iDet < fGeo->NAuxDets() ; ++iDet ){
    geo::AuxDetGeo const& anAuxDetGeo = fGeo->AuxDet(iDet);
    std::string detName = anAuxDetGeo.Name();
    size_t wcnum = 999;
    if( detName == "volAuxDetSensitiveWC1") wcnum = 1;
    if( detName == "volAuxDetSensitiveWC2") wcnum = 2;
    if( detName == "volAuxDetSensitiveWC3") wcnum = 3;
    if( detName == "volAuxDetSensitiveWC4") wcnum = 4;
    if( wcnum != 999 ){
      anAuxDetGeo.GetCenter(centerOfDet);
      fX_cntr[wcnum-1] = centerOfDet[0] * CLHEP::cm;
      fY_cntr[wcnum-1] = centerOfDet[1] * CLHEP::cm;
      fZ_cntr[wcnum-1] = centerOfDet[2] * CLHEP::cm;
    }
  }
  auto tpcGeo = fGeo->begin_TPC_id().get();
  double tpcLocalCenter[3] = {0.};
  tpcGeo->LocalToWorld(tpcLocalCenter, fCenter_of_tpc);
  // put the center of the TPC in world coordinates into mm
  for(int i = 0; i < 3; ++i) {fCenter_of_tpc[i] *= CLHEP::cm;}
  // get the active half width and length of the TPC in mm, geometry
  // returns the values in cm, so have to multiply by CLHEP::cm
  fHalf_z_length_of_tpc = 0.5*tpcGeo->ActiveLength() * CLHEP::cm;
  fHalf_x_length_of_tpc = tpcGeo->ActiveHalfWidth()  * CLHEP::cm;

  // Conversion constants
  fmm_to_m      = 0.001;	
  fGeV_to_MeV 	= 1000.0;
  fDeg_to_Rad   = TMath::Pi()/180.; 
}
//------------------------------------------------------------------------------
WCTrackBuilderAlg::~WCTrackBuilderAlg()
{

}
//------------------------------------------------------------------------------
void WCTrackBuilderAlg::reconfigure( fhicl::ParameterSet const& pset )
{

  fB_field_tesla        = pset.get<float >("BFieldInTesla",      0.       );
  fMCMagneticField      = pset.get<float >("MCMagneticFieldTesla", 0.0);


  //fCentralYKink         = pset.get<float >("CentralYKink",        -0.01    ); //These four are parameters from histos I produced from picky-good tracks
  //fSigmaYKink           = pset.get<float >("SigmaYKink",          0.03      );
  //fCentralYDist         = pset.get<float >("CentralYDist",        0.69      );
  //fSigmaYDist           = pset.get<float >("SigmaYDist",          18.0      );

 // fPrintDisambiguation = false;
  fPickyTracks          = pset.get<bool  >("PickyTracks",         false     );
  fDiagnostics          = pset.get<bool  >("Diagnostics",         false     );
  
  //Survey constants
  //fDelta_z_us           = pset.get<float >("DeltaZus",            1551.15   );  //this will recalculated using geometry instead of hardcoding. 
  //fDelta_z_ds 		= pset.get<float >("DeltaZds",   	  1570.06   );   
  //fL_eff        	= pset.get<float >("LEffective", 	  1145.34706);	
  fL_eff        	= pset.get<float >("LEffective", 	  1204);	
  fMP_X                 = pset.get<float> ("MidplaneInterceptFactor", 1);
  fMP_M                 = pset.get<float> ("MidplaneSlopeFactor", 1);
  fMidplane_intercept   = pset.get<float> ("MidplaneIntercept", 31067.4); 
  fMidplane_slope       = pset.get<float> ("MidplaneSlope", 8.0);

// Where the midplane was before allowing 3 point tracks.  For 3 pt tracks, the intercept and slope will be multiplied by a factor, depending on current setting and WC missed. 
// Don't change these values unless you plan to recalibrate for all current (A) runs. -G Pulliam.
  
  
  
  return;
}
//-------------------------------------------------------------------------------
void WCTrackBuilderAlg::loadXMLDatabaseTableForBField( int run, int subrun )
{
  fRun = run;
  fSubRun = subrun;
  current=fabs(std::stod(fDatabaseUtility->GetIFBeamValue("mid_f_mc7an",fRun,fSubRun)));
//if(fabs(current)>90){fB_field_tesla= .003375*current;}
//if(fabs(current)<90 && fabs(current)>70){fB_field_tesla= .0034875*current;}	
//if(fabs(current)<70 && fabs(current)>50){fB_field_tesla= .003525*current;}	
//if(fabs(current)<50 && fabs(current)>30){fB_field_tesla= .003525*current;}  	
//if(fabs(current)<30){fB_field_tesla= .0035375*current;}  
  fB_field_tesla= (-.1538*pow(10,-4)*pow(current,3)+.2245*pow(10,-2)*pow(current,2)-.1012*current+36.59)*current/10000; // Doug Jensen's cubic equation for magnetic field as a function of current.
  std::cout << "Run: " << fRun << ", Subrun: " << fSubRun << ", B-field: " << fB_field_tesla << std::endl;
}
//--------------------------------------------------------------
//Main function called for each trigger
void WCTrackBuilderAlg::reconstructTracks(std::vector<double> & reco_pz_list,
                                             std::vector<double> & reco_pz2M_list,               
					     std::vector<double> & x_face_list,
					     std::vector<double> & y_face_list,
					     std::vector<double> & incoming_theta_list,
					     std::vector<double> & incoming_phi_list,
					     std::vector<WCHitList> & event_final_tracks,
					     std::vector<std::vector<WCHitList> > & good_hits,
					     bool pickytracks,
					     bool diagnostics, 
					     std::vector<double> & y_kink_list,
					     std::vector<double> & x_dist_list,
					     std::vector<double> & y_dist_list,
					     std::vector<double> & z_dist_list,
					     int & WCMissed,
					     std::vector<TH2F*>  & Recodiff,
					     TH1F* & WCdistribution,
					     float & residual,
					     float (&hit_position_vect)[4][3],
                                             float offset)
{					   
  fPickyTracks = pickytracks;
  fDiagnostics= diagnostics;
//  for(int i=0; i<4; ++i){
//    for(int j=0; j<3; ++j){
//      hit_position_vect_alg[i][j]=hit_position_vect[i][j];
//    }
//  }
  
  //std::cout<<"Time to make some tracks!"<<std::endl;
  //std::cout<<"PickyTracks : "<<fPickyTracks<<"High Yield : "<<fHighYield<<"Diagnostics : "<<fDiagnostics<<std::endl;
  
  initialconst=-999;  //Just a number to use to initialize things before they get filled correctly.
  WCMissed=initialconst;
  
  //Determine if one should skip this trigger based on whether there is exactly one good hit in each wire chamber and axis
  //If there isn't, continue after adding a new empty vector to the reco_pz_array contianer.
  //If there is, then move on with the function.
  //This can be modified to permit more than one good hit in each wire chamber axis - see comments in function
  bool skip = shouldSkipTrigger(good_hits,WCMissed,WCdistribution);
  if( skip == true ) return;
  
  // Assign value to class member variable
  fWCMissed = WCMissed;

  //Depending on if an event has a hit in all 4 WC or whether it missed WC2 or WC3 (but not both), we reconstruct the momentum differently. This code doesn't change from before we allowed 3 point tracks.
  if(fNHits==4){
  
  //At this point, we should have a list of good hits with at least one good hit in X,Y for each WC.
  //Now find all possible combinations of these hits that could form a track, sans consideration
  //of the kinks or end displacements (for now). For the "exactly one" condition set in the above
  //step, this won't matter, but if you want to set the condition to "at least one hit in each WC axis,
  //this will give many combinations.
     trackres=buildFourPointTracks(good_hits,
	                 reco_pz_list,
	                 reco_pz2M_list,
			 x_face_list,
		         y_face_list,
			 incoming_theta_list,
			 incoming_phi_list,
			 event_final_tracks,
			 y_kink_list,
			 x_dist_list,
			 y_dist_list,
			 z_dist_list,
			 WCMissed,
			 hit_position_vect,
                         offset);
					   
   //std::cout<<"Built four point track"<<std::endl;					     
  //Need to use the cut information to whittle down track candidates
//     if( !fPickyTracks ){
//       disambiguateTracks( reco_pz_list,
// 			  y_kink_list,
// 			  x_dist_list,
// 			  y_dist_list,
// 			  z_dist_list,
// 			  track_count,
// 			  x_face_list,
// 			  y_face_list,
// 			  incoming_theta_list,
// 			  incoming_phi_list,
// 			  event_final_tracks,
// 			  lonely_hit_bool);			
//     }
  }
  //If we have a three point track, we do some geometry and then scale the momentum based on some tuning on data.
  if(fNHits==3)
  {
     trackres=buildThreePointTracks(good_hits,
    		          reco_pz_list,
			  x_face_list,
			  y_face_list,
			  incoming_theta_list,
			  incoming_phi_list,
			  event_final_tracks,
			  y_kink_list,
			  x_dist_list,
			  y_dist_list,
			  z_dist_list,
			  WCMissed,
			  hit_position_vect,
			  offset);		
  //std::cout<<"Build three point track"<<std::endl;
  }
  residual=trackres;
  //To compare four point tracks to three point tracks
  if(fDiagnostics==true && fNHits==4){
    MakeDiagnosticPlots( good_hits,
    		         Recodiff,	                
                         reco_pz_list,
			 x_face_list,
		         y_face_list,
			 incoming_theta_list,
			 incoming_phi_list,
			 y_kink_list,
			 x_dist_list,
			 y_dist_list,
			 z_dist_list,
			 WCMissed);
  } 
}
//=====================================================================
bool WCTrackBuilderAlg::shouldSkipTrigger(std::vector<std::vector<WCHitList> > & good_hits,
					  int & WCMissed,
					  TH1F* & WCDist)
{
  //Now determine if we want to skip
  bool skip = false;
  fNHits=0;
  for( size_t iWC = 0; iWC < 4 ; ++iWC ){
    if(good_hits[iWC][0].hits.size()>0 && good_hits[iWC][1].hits.size()>0){++fNHits;}
    else{WCMissed=iWC+1;}
  }
  
  WCDist->Fill(0);
  if(fNHits>2){WCDist->Fill(1);}
  if(fNHits>2 && (WCMissed==1)){WCDist->Fill(2);}
  if(fNHits>2 && (WCMissed==2)){WCDist->Fill(3);}
  if(fNHits>2 && (WCMissed==3)){WCDist->Fill(4);}
  if(fNHits>2 && (WCMissed==4)){WCDist->Fill(5);}
  if(fNHits==4){WCDist->Fill(6);}
  
  //If we don't have 3 or 4 hits, skip.
  if(fNHits<3){
    skip = true;
  } 
  if(fPickyTracks){
    for(size_t iWC=0; iWC<4; ++iWC){
      for(size_t iAX=0; iAX<2; ++iAX){
        if(good_hits[iWC][iAX].hits.size() != 1){  
	//Only allow events with 1 and only 1 hit on each axis
          skip = true;
	  break;
	}
      }
    }   
  }
  if(!fPickyTracks){
    if(WCMissed==1 || WCMissed==4){ //skip events with less than 3 X/Y hits or is missing the first or last WC
      skip = true;
    }  
  }
  //if(WCMissed !=2 && WCMissed !=3){
    //skip=true;
    //}
  if( skip == true ){
    if( fVerbose ){
      std::cout << "skipping this event." << std::endl;
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

//==================================================================================
float WCTrackBuilderAlg::calculateRecoPz(float theta_x_us, float theta_x_ds, float bestTrackSlope )
{
  float num   = (fabs(fB_field_tesla) * fL_eff * fmm_to_m * fGeV_to_MeV );
  float denom = (3.3*(sin(theta_x_ds) - sin(theta_x_us)))*cos(atan(bestTrackSlope));
  return num / denom;
}

//===================================================================================
float WCTrackBuilderAlg::buildFourPointTracks(std::vector<std::vector<WCHitList> > & good_hits,
	                      std::vector<double> & reco_pz_list,
			      std::vector<double> & reco_pz2M_list,
			      std::vector<double> & x_face_list,
		              std::vector<double> & y_face_list,
			      std::vector<double> & incoming_theta_list,
			      std::vector<double> & incoming_phi_list,
			      std::vector<WCHitList> & event_final_tracks,
			      std::vector<double> & y_kink_list,
			      std::vector<double> & x_dist_list,
			      std::vector<double> & y_dist_list,
			      std::vector<double> & z_dist_list,
			      int & WCMissed,
			      float (&hit_position_vect)[4][3],
                              float offset)
{
  float x[4]{0,0,0,0};
  float y[4]{0,0,0,0};
  float z[4]{0,0,0,0};
  float reco_pz=0;
  float reco_pz2M = 0;
  WCHitList best_track;
  std::vector<float> track_stats;
  std::vector<float> bestRegressionStats;
  float bestResSq=initialconst;
    for(int i=0;i<3;i++){
    bestRegressionStats.push_back(initialconst);
    }
  //Loop over all combinations of hits, and find the positions
  for( size_t iHit0 = 0; iHit0 < good_hits[0][0].hits.size(); ++iHit0 ){
    for( size_t iHit1 = 0; iHit1 < good_hits[0][1].hits.size(); ++iHit1 ){
      for( size_t iHit2 = 0; iHit2 < good_hits[1][0].hits.size(); ++iHit2 ){
	for( size_t iHit3 = 0; iHit3 < good_hits[1][1].hits.size(); ++iHit3 ){
	  for( size_t iHit4 = 0; iHit4 < good_hits[2][0].hits.size(); ++iHit4 ){
	    for( size_t iHit5 = 0; iHit5 < good_hits[2][1].hits.size(); ++iHit5 ){
	      for( size_t iHit6 = 0; iHit6 < good_hits[3][0].hits.size(); ++iHit6 ){
		for( size_t iHit7 = 0; iHit7 < good_hits[3][1].hits.size(); ++iHit7 ){
		  WCHitList track;
		  track.hits.push_back(good_hits[0][0].hits[iHit0]);
		  track.hits.push_back(good_hits[0][1].hits[iHit1]);
		  track.hits.push_back(good_hits[1][0].hits[iHit2]);
		  track.hits.push_back(good_hits[1][1].hits[iHit3]);
		  track.hits.push_back(good_hits[2][0].hits[iHit4]);
		  track.hits.push_back(good_hits[2][1].hits[iHit5]);
		  track.hits.push_back(good_hits[3][0].hits[iHit6]);
		  track.hits.push_back(good_hits[3][1].hits[iHit7]);
		  //std::cout<<"track to test made"<<std::endl;
		  findTheHitPositions(track,x,y,z,WCMissed);
		  //std::cout<<"Found Hit position"<<std::endl;
		  //Do regression on the four points to find the straightest track in Y
		  track_stats=Regression(y,z,WCMissed);
		  //std::cout<<"Regressed"<<std::endl;
		  //std::cout<<"Regression value : "<<track_stats[2]<<std::endl;
		  if(track_stats[2]<fabs(bestResSq)){
		    best_track=track;
		    bestRegressionStats=track_stats;
		    bestResSq=track_stats[2];
		  }
		  //std::cout<<"Tested best track"<<std::endl;
		}
	      }	
	    }	
	  }	
	}	
      }		
    }		
  }
  //std::cout<<"All tracks Completed"<<std::endl;
  //std::cout<<bestResSq<<std::endl;
  if(bestResSq<12){
//Now we should have the straightest track in Y, which will be the track that goes to the event.  Now we get the momentum and projections onto the TPC  
    calculateTheMomentum(best_track,x,y,z,reco_pz, reco_pz2M, bestRegressionStats, offset);
  //std::cout<<"Setting 4 point position arrays"<<std::endl;
  for(size_t i=0; i<4; ++i){
    hit_position_vect[i][0]=x[i];
    hit_position_vect[i][1]=y[i];
    hit_position_vect[i][2]=z[i];
  }
  //std::cout<<"4 point position array set!"<<std::endl;
  //float mom_error= CalculateTheMomentumError(x,y,z);
  //std::cout<<"Momentum Calculated"<<std::endl;
  reco_pz_list.push_back(reco_pz);
  reco_pz2M_list.push_back(reco_pz2M);

//We should also have the x,y,z points of the best_track, so now find where it hits the TPC
  projectToTPC(best_track,x,y,z,bestRegressionStats,x_face_list,y_face_list,incoming_theta_list,incoming_phi_list);
  //std::cout<<"TPC Projected"<<std::endl;
  calculateTrackKink_Dists(x,y,z,bestRegressionStats,y_kink_list,x_dist_list,y_dist_list,z_dist_list);
  //std::cout<<"Track Kink/Dist calculated"<<std::endl;
  event_final_tracks.push_back(best_track);
} 
 return bestResSq; 
}
//==================================================================================
void WCTrackBuilderAlg::findTheHitPositions(WCHitList & track,
					    float (&x)[4],
					    float (&y)[4],
					    float (&z)[4],
					    int & WCMissed)
{
  std::vector<float> x_wires;
  std::vector<float> y_wires;
  for(int i=0; i<4; ++i){
     x_wires.push_back(initialconst);
     y_wires.push_back(initialconst);
  }
  //std::cout<<"Track Size: "<<track.hits.size()<<std::endl;
  for(size_t iHit = 0; iHit < track.hits.size() ; ++iHit ){
    int var = iHit % 2;
    int jHit= iHit;
    if( var == 0  && int (jHit != 2*(WCMissed-1)+var) ){ // Skip iHit for the missed WC 
      x_wires[(iHit-var)/2]=(track.hits.at(iHit).wire);
      //std::cout<<"For iHit: "<<iHit<<"x wire: "<<track.hits.at(iHit).wire<<std::endl;
    }
    if( var == 1 && int(jHit != 2*(WCMissed-1)+var)){
      y_wires[(iHit-var)/2]=(track.hits.at(iHit).wire);
      //std::cout<<"For iHit: "<<iHit<<"y wire: "<<track.hits.at(iHit).wire<<std::endl;
    }   
  }
  //std::cout<<"Wires set"<<std::endl;
  float sin13 = sin( fDeg_to_Rad * 13.0 );
  float sin3  = sin( fDeg_to_Rad * 3.0 );
  float cos13 = cos( fDeg_to_Rad * 13.0 );
  float cos3  = cos( fDeg_to_Rad * 3.0 );
  float cosangle[4]={cos13,cos13,cos3,cos3};
  float sinangle[4]={sin13,sin13,sin3,sin3};
  for(int iWC=0; iWC<4; ++iWC){
    if(iWC != WCMissed-1){
  // Find the position of the hit, correcting for the rotation of the WC of 13 degrees for US or 3 degrees for DS
      x[iWC] =  fX_cntr[iWC] + x_wires[iWC] * cosangle[iWC];
      y[iWC] =  fY_cntr[iWC] + y_wires[iWC];
      z[iWC] =  fZ_cntr[iWC] + x_wires[iWC] * sinangle[iWC];
      //std::cout<<"For WC: "<<iWC+1<<"The position is :"<<x[iWC]<<", "<<y[iWC]<<", "<<z[iWC]<<std::endl;
    }
    if(iWC==WCMissed-1){ //Just making sure a value is set, should never be used.
    //std::cout<<"Faking position for missed hit now"<<std::endl;
      x[iWC]=initialconst;
      y[iWC]=initialconst; 
      z[iWC]=initialconst;
        
    }
  }  		 		 	 
}
//=================================================================================
std::vector<float> WCTrackBuilderAlg::Regression(float (&y)[4],
						float (&z)[4],
						int & WCMissed)
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
      if(i != WCMissed-1){ //turning WC# to array index
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
    if(i != WCMissed-1){
    residual= (y[i]-slope*z[i]-intercept)/std::sqrt(1+slope*slope);
    residualsquare += (residual)*(residual); 
    }
  }
  float avgresidual= std::sqrt(residualsquare)/(Npoints-2);
  RegressionValues.push_back(slope);
  RegressionValues.push_back(intercept);
  RegressionValues.push_back(avgresidual);
  return RegressionValues;
} 
//================================================================================    
void WCTrackBuilderAlg::calculateTheMomentum(WCHitList & best_track,
					     float (&x)[4],
					     float (&y)[4],
					     float (&z)[4],
					     float & reco_pz,
					     float & reco_pz2M,
					     std::vector<float> & BestTrackStats,
                                             float offset)
{
  // We need the x,y,z again, which would be for the last track combination tried. 
  // So we need to find the positions again
  
  // Use fWCMissed (class member) since no local version of the variable was
  // passed to this function. 
  findTheHitPositions(best_track,x,y,z,fWCMissed);

  //std::cout<<"Best track hit position found"<<std::endl;
  //Calculate the angle of the track, in the x,z and y,z planes, in upstream(us) and downstream(ds) ends of the WC
  float dx_us=x[1]-x[0];
  float dx_ds=x[3]-x[2];
  //float dy_us=y[1]-y[0];
  //float dy_ds=y[3]-y[2];
  float dz_us=z[1]-z[0];
  float dz_ds=z[3]-z[2];
  float theta_x_us= atan(dx_us/dz_us);
  float theta_x_ds= atan(dx_ds/dz_ds);
  //float theta_y_us= atan(dy_us/dz_us);
  //float theta_y_ds= atan(dy_ds/dz_ds);
  //reco_pz = (fabs(fB_field_tesla) * fL_eff * fmm_to_m * fGeV_to_MeV ) / (3.3*(sin(theta_x_ds) - sin(theta_x_us)))/cos(atan(BestTrackStats[0]));
  reco_pz = calculateRecoPz(theta_x_us,theta_x_ds,BestTrackStats[0]);
  
  //std::cout<<"B: "<<fB_field_tesla<<" momentum: "<<reco_pz<<std::endl;
  
  //Calculating the event using the 2 magnets aproximation
  float theta_central = (theta_x_us*(1+offset)+theta_x_ds*(1-offset))/2.;
  //std::cout<<" Theta Central: "<<theta_central<<std::endl;
  float thetaM1_1 = theta_x_us + 13.*fDeg_to_Rad;
  //std::cout<<" Theta M1_1: "<<thetaM1_1<<std::endl;
  float thetaM1_2 = theta_central + 13.*fDeg_to_Rad;   
  //std::cout<<" Theta M1_2: "<<thetaM1_2<<std::endl;
  double pM1 = double(fabs(fB_field_tesla*(1+offset)) * fL_eff/2. * fmm_to_m * fGeV_to_MeV ) / double(3.3*(sin(thetaM1_2) - sin(thetaM1_1)))/cos(atan(BestTrackStats[0]));
  float thetaM2_1 = theta_central + 3.*fDeg_to_Rad;
  float thetaM2_2 = theta_x_ds + 3.*fDeg_to_Rad;   
  double pM2 = double(fabs(fB_field_tesla*(1-offset)) * fL_eff/2. * fmm_to_m * fGeV_to_MeV ) / double (3.3*(sin(thetaM2_2) - sin(thetaM2_1)))/cos(atan(BestTrackStats[0]));
  reco_pz2M = double(pM2+pM1)/2.;
  //std::cout<<"Reco pz 1 magnet "<<reco_pz<<" Reco Pz 2 magnets "<<reco_pz2M<<std::endl;
  //std::cout<<"Dispersion of the momentum: "<<fabs(reco_pz - reco_pz2M)/reco_pz<<std::endl;
  
}
//==================================================================================
void WCTrackBuilderAlg::projectToTPC(WCHitList & best_track,
					     float (&x)[4],
					     float (&y)[4],
					     float (&z)[4],
					     std::vector<float> & bestRegressionStats,
					     std::vector<double> & x_face_list,
					     std::vector<double> & y_face_list,
					     std::vector<double> & incoming_theta_list,
					     std::vector<double> & incoming_phi_list)
{					         
  // The intercept of the regression is where the track hits at z=0, which is the TPC face, 
  // so we get the y position on the TPC for free
  y_face_list.push_back(bestRegressionStats[1]);

  // Though the regression probably gives the best extrapolation to the TPC, we can only use 
  // it for y projections. If we need x values, we have to use the positions ds, which could 
  // introduce error if we have either WC3,4 wrong.  For the x position at the TPC, we need 
  // to use the ds WC to extrapolate a line to the TPC (z==0)
  float dx_ds = x[3] - x[2];  
  float dy_ds = y[3] - y[2];
  float dz_ds = z[3] - z[2];
  float ds_slopex = dx_ds / dz_ds;
  x_face_list.push_back( x[3] - ds_slopex*z[3] );
  
  // Phi is measured CCW around a circle starting at (x,y) = (1,0), in radians. 
  // Theta is the opening angle relative to the z-axis.
  // Make TVector3 out of the dx,dy and extract phi/theta from it.
  TVector3 vec(dx_ds, dy_ds, dz_ds);
  incoming_phi_list.push_back(vec.Phi());
  incoming_theta_list.push_back(vec.Theta());
}
//====================================================================================
void WCTrackBuilderAlg::calculateTrackKink_Dists(float (&x)[4],
						 float (&y)[4],
						 float (&z)[4],
						 std::vector<float> & track_stats,
						 std::vector<double> & y_kink_list,
						 std::vector<double> & x_dist_list,
						 std::vector<double> & y_dist_list,
						 std::vector<double> & z_dist_list)
{
  float dx_us=x[1]-x[0];
  float dy_us=y[1]-y[0];
  float dz_us=z[1]-z[0];
  float dx_ds=x[3]-x[2];
  float dy_ds=y[3]-y[2];
  float dz_ds=z[3]-z[2];
  float x_us_slope=dx_us/dz_us;
  float y_us_slope=dy_us/dz_us;
  float x_ds_slope=dx_ds/dz_ds;
  float y_ds_slope=dy_ds/dz_ds;
  float x_us_int=x[1]-x_us_slope*z[1];
  float y_us_int=y[1]-y_us_slope*z[1];
  float x_ds_int=x[3]-x_ds_slope*z[3];
  float y_ds_int=y[3]-y_ds_slope*z[3]; 
//Now we have the equations of the lines, x=mz+b and y=mz+b, for the US and DS legs.  First Y_kink, the angle difference between the us and ds legs.
  y_kink_list.push_back(atan(y_ds_slope)-atan(y_us_slope));
//Because we need a way to compare tracks, regardless of current setting, we have to have a standard for the midplane. We will use the normal midplane (fMP_X=fMP_M=1)
   float z_mp_us=(fMidplane_intercept-x_us_int)/(x_us_slope-1/(tan(8.0*fDeg_to_Rad)));  //X,Y,Z where US interesects Midplane.
   float x_mp_us=x_us_slope*z_mp_us+x_us_int;  
   float y_mp_us=y_us_slope*z_mp_us+y_us_int;
   float z_mp_ds=(fMidplane_intercept-x_ds_int)/(x_ds_slope-1/(tan(8.0*fDeg_to_Rad)));  //X,Y,Z where DS interesects Midplane.
   float x_mp_ds=x_ds_slope*z_mp_ds+x_ds_int;
   float y_mp_ds=y_ds_slope*z_mp_ds+y_ds_int;
   x_dist_list.push_back(x_mp_ds-x_mp_us);
   y_dist_list.push_back(y_mp_ds-y_mp_us); 
   z_dist_list.push_back(z_mp_ds-z_mp_us);  
}
//====================================================================================
float WCTrackBuilderAlg::buildThreePointTracks(std::vector<std::vector<WCHitList> > & good_hits,
	                     		      std::vector<double> & reco_pz_list,
			     	              std::vector<double> & x_face_list,
		                              std::vector<double> & y_face_list,
			                      std::vector<double> & incoming_theta_list,
			                      std::vector<double> & incoming_phi_list,
			                      std::vector<WCHitList> & event_final_tracks,
					      std::vector<double> & y_kink_list,
					      std::vector<double> & x_dist_list,
					      std::vector<double> & y_dist_list,
					      std::vector<double> & z_dist_list,
					      int & WCMissed,
					       float (&hit_position_vect)[4][3],
					       float offset)
{
//Code here is similar to the buildFourPointTrack version, but we have to allow for missed WC, and require some additional geometry to find the momentum
  float x[4]{0,0,0,0};
  float y[4]{0,0,0,0};
  float z[4]{0,0,0,0};
  std::vector<float> missed_wire_hits;
  std::vector<float> bestRegressionStats;
  WCHitList best_track;
  float bestResSq=initialconst;
  std::vector<float> track_stats;
    for(int i=0;i<3;i++){
      bestRegressionStats.push_back(initialconst);
    }
    for(int i=0; i<2; ++i){
      missed_wire_hits.push_back(initialconst);
    }
  //std::cout<<"WCMissed"<<WCMissed<<std::endl;
  //Loop over all combinations of hits, and find the positions
  //Because I cannot figure a better way of doing this, we will do two loops, one if WC2 is missed, the other if WC3 is missed. The code will be the same, just filling an empty hit list in place of the missed WC. -GP
  WCHitList NullList;
  WCHit FakeHit;
  FakeHit.wire=initialconst;
  FakeHit.time=initialconst;
  FakeHit.hit_index=initialconst;
  FakeHit.cluster_index=initialconst;
  FakeHit.isVisited=false;
  NullList.hits.push_back(FakeHit);

  if(WCMissed==2){
    for( size_t iHit0 = 0; iHit0 < good_hits[0][0].hits.size(); ++iHit0 ){
      for( size_t iHit1 = 0; iHit1 < good_hits[0][1].hits.size(); ++iHit1 ){
        //for( size_t iHit2 = 0; iHit2 < good_hits[1][0].hits.size(); ++iHit2 ){
	  //for( size_t iHit3 = 0; iHit3 < good_hits[1][1].hits.size(); ++iHit3 ){
	    for( size_t iHit4 = 0; iHit4 < good_hits[2][0].hits.size(); ++iHit4 ){
	      for( size_t iHit5 = 0; iHit5 < good_hits[2][1].hits.size(); ++iHit5 ){
	        for( size_t iHit6 = 0; iHit6 < good_hits[3][0].hits.size(); ++iHit6 ){
		  for( size_t iHit7 = 0; iHit7 < good_hits[3][1].hits.size(); ++iHit7 ){
		    WCHitList track;
		    track.hits.push_back(good_hits[0][0].hits[iHit0]);
		    track.hits.push_back(good_hits[0][1].hits[iHit1]);
		    //std::cout<<"faking the hit"<<std::endl;
		    track.hits.push_back(NullList.hits[0]);
		    track.hits.push_back(NullList.hits[0]);
		    //std::cout<<"Faked Hit"<<std::endl;
		    //track.hits.push_back(good_hits[1][0].hits[iHit2]);
		    //track.hits.push_back(good_hits[1][1].hits[iHit3]);
		    track.hits.push_back(good_hits[2][0].hits[iHit4]);
		    track.hits.push_back(good_hits[2][1].hits[iHit5]);
		    track.hits.push_back(good_hits[3][0].hits[iHit6]);
		    track.hits.push_back(good_hits[3][1].hits[iHit7]);
		    findTheHitPositions(track,x,y,z,WCMissed);
		    //std::cout<<"Three Point Hit position found"<<std::endl;
		    //Do regression on the four points to find the straightest track in Y
		    track_stats=Regression(y,z,WCMissed);
		    //std::cout<<"Track Res Sq :"<<track_stats[2]<<std::endl;
		    //std::cout<<"Three Regressed"<<std::endl;
		    if(track_stats[2]<fabs(bestResSq)){
		      best_track=track;
		      bestRegressionStats=track_stats;
		      bestResSq=track_stats[2];
		    }
		    //std::cout<<"Three track checked for best track"<<std::endl;
		  }
	        }	
	      }	
	    }	
	  //}	
        //}		
      }		
    }
  }
  if(WCMissed==3){
    for( size_t iHit0 = 0; iHit0 < good_hits[0][0].hits.size(); ++iHit0 ){
      for( size_t iHit1 = 0; iHit1 < good_hits[0][1].hits.size(); ++iHit1 ){
        for( size_t iHit2 = 0; iHit2 < good_hits[1][0].hits.size(); ++iHit2 ){
	  for( size_t iHit3 = 0; iHit3 < good_hits[1][1].hits.size(); ++iHit3 ){
	    //for( size_t iHit4 = 0; iHit4 < good_hits[2][0].hits.size(); ++iHit4 ){
	      //for( size_t iHit5 = 0; iHit5 < good_hits[2][1].hits.size(); ++iHit5 ){
	        for( size_t iHit6 = 0; iHit6 < good_hits[3][0].hits.size(); ++iHit6 ){
		  for( size_t iHit7 = 0; iHit7 < good_hits[3][1].hits.size(); ++iHit7 ){
		    WCHitList track;
		    track.hits.push_back(good_hits[0][0].hits[iHit0]);
		    track.hits.push_back(good_hits[0][1].hits[iHit1]);
		    track.hits.push_back(good_hits[1][0].hits[iHit2]);
		    track.hits.push_back(good_hits[1][1].hits[iHit3]);
		    track.hits.push_back(NullList.hits[0]);
		    track.hits.push_back(NullList.hits[0]);
		    //track.hits.push_back(good_hits[2][0].hits[iHit4]);
		    //track.hits.push_back(good_hits[2][1].hits[iHit5]);
		    track.hits.push_back(good_hits[3][0].hits[iHit6]);
		    track.hits.push_back(good_hits[3][1].hits[iHit7]);
		    findTheHitPositions(track,x,y,z,WCMissed);
		    //Do regression on the four points to find the straightest track in Y
		    track_stats=Regression(y,z,WCMissed);
		     //std::cout<<"Track Res Sq :"<<track_stats[2]<<std::endl;
		    if(track_stats[2]<fabs(bestResSq)){
		      best_track=track;
		      bestRegressionStats=track_stats;
		      bestResSq=track_stats[2];
		      //std::cout<<track_stats[2]<<std::endl;
		    }
		  }
	        }	
	      //}	
	    //}	
	  }	
        }		
      }		
    }
  }
    if(WCMissed==4){
    for( size_t iHit0 = 0; iHit0 < good_hits[0][0].hits.size(); ++iHit0 ){
      for( size_t iHit1 = 0; iHit1 < good_hits[0][1].hits.size(); ++iHit1 ){
        for( size_t iHit2 = 0; iHit2 < good_hits[1][0].hits.size(); ++iHit2 ){
	  for( size_t iHit3 = 0; iHit3 < good_hits[1][1].hits.size(); ++iHit3 ){
	    for( size_t iHit4 = 0; iHit4 < good_hits[2][0].hits.size(); ++iHit4 ){
	      for( size_t iHit5 = 0; iHit5 < good_hits[2][1].hits.size(); ++iHit5 ){
	        //for( size_t iHit6 = 0; iHit6 < good_hits[3][0].hits.size(); ++iHit6 ){
		  //for( size_t iHit7 = 0; iHit7 < good_hits[3][1].hits.size(); ++iHit7 ){
		    WCHitList track;
		    track.hits.push_back(good_hits[0][0].hits[iHit0]);
		    track.hits.push_back(good_hits[0][1].hits[iHit1]);
		    track.hits.push_back(good_hits[1][0].hits[iHit2]);
		    track.hits.push_back(good_hits[1][1].hits[iHit3]);
		    //track.hits.push_back(NullList.hits[0]);
		    //track.hits.push_back(NullList.hits[0]);
		    track.hits.push_back(good_hits[2][0].hits[iHit4]);
		    track.hits.push_back(good_hits[2][1].hits[iHit5]);
		    track.hits.push_back(NullList.hits[0]);
		    track.hits.push_back(NullList.hits[0]);
		    //track.hits.push_back(good_hits[3][0].hits[iHit6]);
		    //track.hits.push_back(good_hits[3][1].hits[iHit7]);
		    findTheHitPositions(track,x,y,z,WCMissed);
		    //Do regression on the four points to find the straightest track in Y
		     track_stats=Regression(y,z,WCMissed);
		     //std::cout<<"Track Res Sq :"<<track_stats[2]<<std::endl;
		    if(track_stats[2]<fabs(bestResSq)){
		      best_track=track;
		      bestRegressionStats=track_stats;
		      bestResSq=track_stats[2];
		      //std::cout<<track_stats[2]<<std::endl;
		    }
		  //}
	        //}	
	      }	
	    }	
	  }	
        }		
      }		
    }
  }
  if(bestResSq<12){
  //std::cout<<"Best track residual"<<bestRegressionStats[2]<<std::endl;
  float reco_pz_three=0;  
//Now we should have the straightest track in Y, which will be the track that goes to the event.  Now we get the momentum and projections onto the TPC  
  calculateTheThreePointMomentum(best_track,x,y,z,reco_pz_three,bestRegressionStats,WCMissed);
 // std::cout<<"Three Momentum calculated"<<std::endl;
  reco_pz_list.push_back(reco_pz_three);
//Now that the momentum is scaled back to what we would get with four points, try to find where the missed point would be to reconstruct this scaled momentum
  extrapolateTheMissedPoint(best_track,x,y,z,reco_pz_three,bestRegressionStats, missed_wire_hits,WCMissed);
  //std::cout<<"Missed Point Extrapolated"<<std::endl;

//Put the extrapolated hit onto the track
  if(WCMissed==2){
    best_track.hits[2].wire=missed_wire_hits[0];
    best_track.hits[3].wire=missed_wire_hits[1];
  }
  if(WCMissed==3){
    best_track.hits[4].wire=missed_wire_hits[0];
    best_track.hits[5].wire=missed_wire_hits[1];
  } 
  if(WCMissed==4){
    best_track.hits[6].wire=missed_wire_hits[0];
    best_track.hits[7].wire=missed_wire_hits[1];
  } 
  findTheHitPositions(best_track,x,y,z,initialconst); //Find the hit positions again, with the now complete track, with WCMissed=initialconst to avoid skipping the hit we extrapolated.
 // std::cout<<"Setting the 3 point position vector"<<std::endl;
  for(int i=0; i<4; ++i){
    hit_position_vect[i][0]=x[i];
    hit_position_vect[i][1]=y[i];
    hit_position_vect[i][2]=z[i];
  }
  //std::cout<<"3 point position set!"<<std::endl;
//We should also have the x,y,z points of the best_track, so now find where it hits the TPC
   calculateTrackKink_Dists(x,y,z,bestRegressionStats,y_kink_list,x_dist_list,y_dist_list,z_dist_list);
  projectToTPC(best_track,x,y,z,bestRegressionStats,x_face_list,y_face_list,incoming_theta_list,incoming_phi_list);
event_final_tracks.push_back(best_track);
}
return bestResSq;
}
//======================================================================================
void WCTrackBuilderAlg::calculateTheThreePointMomentum(WCHitList & best_track,
						       float(&x)[4],
						       float(&y)[4],
						       float(&z)[4],
						       float & reco_pz,
						       std::vector<float> & BestTrackStats,
						       int & WCMissed)
{
  findTheHitPositions(best_track,x,y,z,WCMissed);
  //std::cout<<"3 momentum position found"<<std::endl;
  //std::cout<<current<<std::endl;
  //Now depending on which WC was missed and which current the magnets ran, we use different calibration constants for the line we force both the US and DS legs to cross. This was tuned to a subset of data in Run 1.
  //Only have data to tune for 60A runs and 100A runs
  //WC 2 missed calibration
  if(WCMissed==2){
    if(current>50 && current< 70){
      fMP_M=.98;
      fMP_X=1.05;
    }
    if(current>90){
      fMP_M=1.01;
      fMP_X=1;
    }
    fMP_M=1;
    fMP_X=1;
    //std::cout<<"Midplane set"<<std::endl;
    //float midplane_slope=1/(tan((3.141592654/180))*8.0)*fMP_M;
    float midplane_slope=1./tan(fDeg_to_Rad*8.0)*fMP_M;
    float midplane_intercept=fMidplane_intercept*fMP_X;
    //std::cout<<"fMidplane int "<<fMidplane_intercept<<std::endl;
    //std::cout<<"mp int "<<midplane_intercept<<std::endl;
    float ds_dz=z[3]-z[2];
    float ds_dx=x[3]-x[2];
    float ds_slope=ds_dx/ds_dz;
    float ds_int_x= x[2]-ds_slope*z[2];
    float z_ds=(midplane_intercept-ds_int_x)/(ds_slope-midplane_slope); //Solving tan8*Z+Bmp == Mds*Z+Bds
    float x_ds=ds_slope*z_ds+ ds_int_x;
    //std::cout<<"x hits "<<x[2]<<", "<<x[3]<<std::endl;
    //std::cout<<"z hits "<<z[2]<<", "<<z[3]<<std::endl;
    //std::cout<<"X centers "<<fX_cntr[2]<<", "<<fX_cntr[3]<<std::endl;
    //std::cout<<"Z centers "<<fZ_cntr[2]<<", "<<fZ_cntr[3]<<std::endl;
    //std::cout<<"ds_dx, ds_dz, slope, int"<<ds_dx<<", "<<ds_dz<<", "<<ds_slope<<", "<<ds_int_x<<std::endl;
    //std::cout<<"midplane x,z"<<x_ds<<" "<<z_ds<<std::endl;
    float us_dz=(z_ds-z[0]);
    float us_dx=(x_ds-x[0]);       
    float theta_x_us=atan(us_dx/us_dz);
    float theta_x_ds=atan(ds_dx/ds_dz); 
    //std::cout<<"WC 2 Theta US: "<<theta_x_us<<std::endl;
    //std::cout<<"WC 2 Theta DS: "<<theta_x_ds<<std::endl;
    //float theta_y_ds=atan(((y[3]-y[2])/(z[3]-z[2])));
    //reco_pz = (fabs(fB_field_tesla) * fL_eff * fmm_to_m * fGeV_to_MeV ) / (float (3.3*(theta_x_ds - theta_x_us)*cos(theta_y_ds)));
    //std::cout<<"S2 mom: "<<reco_pz<<std::endl;
    //reco_pz = (fabs(fB_field_tesla) * fL_eff * fmm_to_m * fGeV_to_MeV ) / (float (3.3*(sin(theta_x_ds) - sin(theta_x_us))*cos(atan(BestTrackStats[0]))));
    reco_pz = calculateRecoPz(theta_x_us,theta_x_ds,BestTrackStats[0]);
    //std::cout<<"momentum set"<<std::endl;
    reco_pz=(reco_pz+13)/1.15;


    
    //std::cout<<"momentum scaled"<<std::endl;
  }
  //WC 3 missed calibration
  if(WCMissed==3){
    if(current>50 && current< 70){
      fMP_M=1;
      fMP_X=1.02;
    }
    if(current>90){
      fMP_M=1.05;
      fMP_X=1.05;
    }  
      fMP_M=1;
      fMP_X=1;
      float midplane_slope=1/tan(fDeg_to_Rad*8.0)*fMP_M;
      float midplane_intercept=fMidplane_intercept*fMP_X;    
      float us_dz=z[1]-z[0];
      float us_dx=x[1]-x[0];
      float us_slope=us_dx/us_dz;
      float us_int_x= x[1]-us_slope*z[1];
      float z_us=(midplane_intercept-us_int_x)/(us_slope-midplane_slope); //Solving tan8*Z+Bmp == Mus*Z+Bus
      float x_us=us_slope*z_us + us_int_x;
      float ds_dz=-(z_us-z[3]);
      float ds_dx=-(x_us-x[3]);    
      float theta_x_us=asin(us_dx/pow(us_dx*us_dx+us_dz*us_dz,.5));
      float theta_x_ds=asin(ds_dx/pow(ds_dx*ds_dx+ds_dz*ds_dz,.5)); 
      //std::cout<<"WC 3 Theta US: "<<theta_x_us<<std::endl;
      //std::cout<<"WC 3 Theta DS: "<<theta_x_ds<<std::endl;
     // float theta_y_us=atan(((y[1]-y[0])/(z[1]-z[0])));
      //reco_pz =(fabs(fB_field_tesla) * fL_eff * fmm_to_m * fGeV_to_MeV ) /( float(3.3 * (theta_x_ds-theta_x_us)*cos(theta_y_us))); 
      //reco_pz =(fabs(fB_field_tesla) * fL_eff * fmm_to_m * fGeV_to_MeV ) /( float(3.3 * (sin(theta_x_ds)-sin(theta_x_us))*cos(atan(BestTrackStats[0])))); 
      reco_pz = calculateRecoPz(theta_x_us,theta_x_ds,BestTrackStats[0]);
      //Caibrate depending on WCMissed
      reco_pz=(reco_pz-7)/.91;
    
  }
   if(WCMissed==4){
    if(current>50 && current< 70){
      fMP_M=1;
      fMP_X=1.02;
    }
    if(current>90){
      fMP_M=1.05;
      fMP_X=1.05;
    }  
      fMP_M=1;
      fMP_X=1;
      float midplane_slope=1./tan(fDeg_to_Rad*8.0)*fMP_M;
      float midplane_intercept=fMidplane_intercept*fMP_X;    
      float us_dz=z[1]-z[0];
      float us_dx=x[1]-x[0];
      float us_slope=us_dx/us_dz;
      float us_int_x= x[1]-us_slope*z[1];
      float z_us=(midplane_intercept-us_int_x)/(us_slope-midplane_slope); //Solving tan8*Z+Bmp == Mus*Z+Bus
      float x_us=us_slope*z_us + us_int_x;
      float ds_dz=-(z_us-z[2]);
      float ds_dx=-(x_us-x[2]);    
      float theta_x_us=asin(us_dx/pow(us_dx*us_dx+us_dz*us_dz,.5));
      float theta_x_ds=asin(ds_dx/pow(ds_dx*ds_dx+ds_dz*ds_dz,.5)); 
     // std::cout<<"WC 4 Theta US: "<<theta_x_us<<std::endl;
      //std::cout<<"WC 4 Theta DS: "<<theta_x_ds<<std::endl;
     // float theta_y_us=atan(((y[1]-y[0])/(z[1]-z[0])));
      //reco_pz =(fabs(fB_field_tesla) * fL_eff * fmm_to_m * fGeV_to_MeV ) /( float(3.3 * (theta_x_ds-theta_x_us)*cos(theta_y_us))); 
      //reco_pz =(fabs(fB_field_tesla) * fL_eff * fmm_to_m * fGeV_to_MeV ) /( float(3.3 * (sin(theta_x_ds)-sin(theta_x_us))*cos(atan(BestTrackStats[0])))); 
      reco_pz = calculateRecoPz(theta_x_us,theta_x_ds,BestTrackStats[0]);
   //Calibrate depending on current
      reco_pz=(reco_pz-37.7883)/.701274;
    
  }   
}
//======================================================================================
void WCTrackBuilderAlg::extrapolateTheMissedPoint(WCHitList & best_track,
					          float(&x)[4],
					          float(&y)[4],
					          float(&z)[4],
					          float & reco_pz,
					          std::vector<float> & BestTrackStats,
						  std::vector<float> & missed_wires,
						  int & WCMissed)
{
  if(WCMissed==2){
    //float dx_us=x[1]-x[0];
    float dx_ds=x[3]-x[2];
    //float dy_us=y[1]-y[0];
    //float dy_ds=y[3]-y[2];
    //float dz_us=z[1]-z[0];
    float dz_ds=z[3]-z[2];
    //float theta_x_us= atan(dx_us/dz_us);
    float theta_x_ds= atan(dx_ds/dz_ds);
    //float theta_y_us= atan(dy_us/dz_us);
    //float theta_y_ds= atan(dy_ds/dz_ds);
    //reco_pz = (fabs(fB_field_tesla) * fL_eff * fmm_to_m * fGeV_to_MeV ) / (3.3*(theta_x_ds - theta_x_us))/cos(theta_y_ds);
    float theta_x_us_reco=asin(sin(theta_x_ds)-(fabs(fB_field_tesla)*fL_eff*fmm_to_m*fGeV_to_MeV/(3.3*reco_pz*cos(atan(BestTrackStats[0])))));
    //float theta_x_us_reco=-(fabs(fB_field_tesla) * fL_eff * fmm_to_m * fGeV_to_MeV )/(3.3*reco_pz*cos(atan(BestTrackStats[0]))) + theta_x_ds;
    //We still need X2,Z2, both of which are functions of the X wire of the hit.  We solve for the x wire in WC2 and then find X2,Z2.
    float missedXwire=-(fX_cntr[1]-x[0]-(fZ_cntr[1]-z[0])*tan(theta_x_us_reco))/(cos(fDeg_to_Rad*13.0)-sin(fDeg_to_Rad*13.0)*tan(theta_x_us_reco));
    x[1]=fX_cntr[1]+cos(fDeg_to_Rad*13.0)*missedXwire;
    z[1]=fZ_cntr[1]+sin(fDeg_to_Rad*13.0)*missedXwire;
    //Then use the regression to get y from z.
    y[1]=BestTrackStats[0]*z[1]+BestTrackStats[1];
    float missedYwire=y[1]-fY_cntr[1];
    //std::cout<<"Missed X wire: "<<missedXwire<<" Missed Y Wire : "<<missedYwire<<std::endl;
    //std::cout<<"extrapolated x: "<<x[1]<<" extrapolated y: "<<y[1]<<std::endl;
    missed_wires[0]=missedXwire;
    missed_wires[1]=missedYwire;
  }
  if(WCMissed==3){
    float dx_us=x[1]-x[0];
    //float dx_ds=x[3]-x[2];
    //float dy_us=y[1]-y[0];
    //float dy_ds=y[3]-y[2];
    float dz_us=z[1]-z[0];
    //float dz_ds=z[3]-z[2];
    float theta_x_us= atan(dx_us/dz_us);
    //float theta_x_ds= atan(dx_ds/dz_ds);
    //float theta_y_us= atan(dy_us/dz_us);
    //float theta_y_ds= atan(dy_ds/dz_ds);
    float theta_x_ds_reco=asin(sin(theta_x_us)+(fabs(fB_field_tesla)*fL_eff*fmm_to_m*fGeV_to_MeV/(3.3*reco_pz*cos(atan(BestTrackStats[0])))));
    float missedXwire=(x[3]-fX_cntr[2]-(z[3]-fZ_cntr[2])*tan(theta_x_ds_reco))/(cos(fDeg_to_Rad*3.0)-(sin(fDeg_to_Rad*3.0))*tan(theta_x_ds_reco));
    x[2]=fX_cntr[2]+cos(fDeg_to_Rad*3.0)*missedXwire;
    z[2]=fZ_cntr[2]+sin(fDeg_to_Rad*3.0)*missedXwire;
    y[2]=BestTrackStats[0]*z[2]+BestTrackStats[1];
    float missedYwire=y[2]-fY_cntr[2];
    //std::cout<<"Missed X wire: "<<missedXwire<<" Missed Y Wire : "<<missedYwire<<std::endl;
    //std::cout<<"extrapolated x: "<<x[2]<<" extrapolated y: "<<y[2]<<std::endl;
    missed_wires[0]=missedXwire;
    missed_wires[1]=missedYwire;
  }
    if(WCMissed==4){
    float dx_us=x[1]-x[0];
    //float dx_ds=x[3]-x[2];
    //float dy_us=y[1]-y[0];
    //float dy_ds=y[3]-y[2];
    float dz_us=z[1]-z[0];
    //float dz_ds=z[3]-z[2];
    float theta_x_us= atan(dx_us/dz_us);
    //float theta_x_ds= atan(dx_ds/dz_ds);
    //float theta_y_us= atan(dy_us/dz_us);
    //float theta_y_ds= atan(dy_ds/dz_ds);
    float theta_x_ds_reco=asin(sin(theta_x_us)+(fabs(fB_field_tesla)*fL_eff*fmm_to_m*fGeV_to_MeV/(3.3*reco_pz*cos(atan(BestTrackStats[0])))));
    float missedXwire=(x[2]-fX_cntr[3]-(z[2]-fZ_cntr[3])*tan(theta_x_ds_reco))/(cos(fDeg_to_Rad*3.0)-(sin(fDeg_to_Rad*3.0))*tan(theta_x_ds_reco));
    x[3]=fX_cntr[3]+cos(fDeg_to_Rad*3.0)*missedXwire;
    z[3]=fZ_cntr[3]+sin(fDeg_to_Rad*3.0)*missedXwire;
    y[3]=BestTrackStats[0]*z[3]+BestTrackStats[1];
    float missedYwire=y[3]-fY_cntr[3];
    //std::cout<<"Missed X wire: "<<missedXwire<<" Missed Y Wire : "<<missedYwire<<std::endl;
    //std::cout<<"extrapolated x: "<<x[2]<<" extrapolated y: "<<y[2]<<std::endl;
    missed_wires[0]=missedXwire;
    missed_wires[1]=missedYwire;
  }   
}
//======================================================================================
void WCTrackBuilderAlg::MakeDiagnosticPlots(std::vector<std::vector<WCHitList> > & good_hits,
					    std::vector<TH2F*> & Recodiff,	                      
					    std::vector<double> & reco_pz_list,
			     		    std::vector<double> & x_face_list,
		              		    std::vector<double> & y_face_list,
			      		    std::vector<double> & incoming_theta_list,
			      		    std::vector<double> & incoming_phi_list,
			      		    std::vector<double> & y_kink_list,
			      		    std::vector<double> & x_dist_list,
			      	            std::vector<double> & y_dist_list,
			      	            std::vector<double> & z_dist_list,
			      	            int & WCMissed)
{
//We recalculate all the variables with four points, then do both the WCMissed=2 and WCMissed=3 functions to compare each to the four point versions.
  float x[4]{0,0,0,0};
  float y[4]{0,0,0,0};
  float z[4]{0,0,0,0};
  float reco_pz=0;
  float reco_pz2M = 0;
  WCHitList best_track;
  std::vector<float> bestRegressionStats;
  std::vector<float> missed_wire_hits;
  float bestResSq=initialconst;
    for(int i=0;i<3;i++){
    bestRegressionStats.push_back(initialconst);
    }
    for(int i=0; i<2; ++i){
    missed_wire_hits.push_back(initialconst);
    }
  //Loop over all combinations of hits, and find the positions
  for( size_t iHit0 = 0; iHit0 < good_hits[0][0].hits.size(); ++iHit0 ){
    for( size_t iHit1 = 0; iHit1 < good_hits[0][1].hits.size(); ++iHit1 ){
      for( size_t iHit2 = 0; iHit2 < good_hits[1][0].hits.size(); ++iHit2 ){
	for( size_t iHit3 = 0; iHit3 < good_hits[1][1].hits.size(); ++iHit3 ){
	  for( size_t iHit4 = 0; iHit4 < good_hits[2][0].hits.size(); ++iHit4 ){
	    for( size_t iHit5 = 0; iHit5 < good_hits[2][1].hits.size(); ++iHit5 ){
	      for( size_t iHit6 = 0; iHit6 < good_hits[3][0].hits.size(); ++iHit6 ){
		for( size_t iHit7 = 0; iHit7 < good_hits[3][1].hits.size(); ++iHit7 ){
		//std::cout<<iHit7<<std::endl;
		  WCHitList track;
		  //std::cout<<"At the top of the loop"<<std::endl;
		  //std::cout<<iHit0<<iHit1<<iHit2<<iHit3<<iHit4<<iHit5<<iHit6<<iHit7<<std::endl;
		  track.hits.push_back(good_hits[0][0].hits[iHit0]);
		  track.hits.push_back(good_hits[0][1].hits[iHit1]);
		  track.hits.push_back(good_hits[1][0].hits[iHit2]);
		  track.hits.push_back(good_hits[1][1].hits[iHit3]);
		  track.hits.push_back(good_hits[2][0].hits[iHit4]);
		  track.hits.push_back(good_hits[2][1].hits[iHit5]);
		  track.hits.push_back(good_hits[3][0].hits[iHit6]);
		  track.hits.push_back(good_hits[3][1].hits[iHit7]);
		  //std::cout<<"track to test made"<<std::endl;
		  findTheHitPositions(track,x,y,z,WCMissed);
		  //std::cout<<"Found Hit position"<<std::endl;
		  //Do regression on the four points to find the straightest track in Y
		  std::vector<float> track_stats=Regression(y,z,WCMissed);
		  //std::cout<<"Regressed"<<std::endl;
		  //std::cout<<"Regression value : "<<track_stats[2]<<std::endl;
		  if(track_stats[2]<fabs(bestResSq)){
		    best_track=track;
		    bestRegressionStats=track_stats;
		    bestResSq=track_stats[2];
		  }
		  //std::cout<<"Tested best track"<<std::endl;
		}
	      }	
	    }	
	  }	
	}	
      }		
    }		
  }
  //std::cout<<"All tracks Completed"<<std::endl;
  //std::cout<<bestResSq<<std::endl;
  //if(bestResSq<10){
  Recodiff[58]->Fill(bestResSq,bestResSq);
  float fourres=bestResSq;
//Now we should have the straightest track in Y, which will be the track that goes to the event.  Now we get the momentum and projections onto the TPC  
  float offset = 0.;
  calculateTheMomentum(best_track,x,y,z,reco_pz, reco_pz2M, bestRegressionStats, offset);
  float fourmom=reco_pz;
  std::cout<<"Reco 2M: "<<reco_pz2M<<std::endl;
  float fourx[4]={x[0],x[1],x[2],x[3]};
  float foury[4]={y[0],y[1],y[2],y[3]};
  float fourz[4]={z[0],z[1],z[2],z[3]};
  //float xmidplane_calculated=0;
  //float ymidplane_calculated=0;
  //float zmidplane_calculated=0;
  float closest_distance_fake=0;
  //Recodiff[97]->Fill(zmidplane_calculated,xmidplane_calculated);
  //Recodiff[98]->Fill(closest_distance,closest_distance);
  TVector3 MidplaneVect=PlotTheMidplane(fourx,foury,fourz,closest_distance_fake);
  float closest_distance=PlotTheMidplane(fourx,foury,fourz);
  Recodiff[97]->Fill(MidplaneVect.Z(),MidplaneVect.X());
  Recodiff[98]->Fill(closest_distance,closest_distance);
  float x2minus[4]= {x[0],x[1]-3,x[2],x[3]};
  calculateTheMomentumGiven(best_track,x2minus,y,z,reco_pz,bestRegressionStats);
  float mom2minus=reco_pz;
  float x2plus[4]= {x[0],x[1]+3,x[2],x[3]};
  calculateTheMomentumGiven(best_track,x2plus,y,z,reco_pz,bestRegressionStats);
  float mom2plus=reco_pz;
  float x3minus[4]= {x[0],x[1],x[2]-3,x[3]};
  calculateTheMomentumGiven(best_track,x3minus,y,z,reco_pz,bestRegressionStats);
  float mom3minus=reco_pz;
  float x3plus[4]= {x[0],x[1],x[2]+3,x[3]};
  calculateTheMomentumGiven(best_track,x3plus,y,z,reco_pz,bestRegressionStats);
  float mom3plus=reco_pz;
  
  float xplusplus[4]= {x[0],x[1]+3,x[2]+3,x[3]};
  calculateTheMomentumGiven(best_track,xplusplus,y,z,reco_pz,bestRegressionStats);
  float momplusplus=reco_pz;
    float xplusminus[4]= {x[0],x[1]+3,x[2]-3,x[3]};
  calculateTheMomentumGiven(best_track,xplusminus,y,z,reco_pz,bestRegressionStats);
  float momplusminus=reco_pz;
    float xminusplus[4]= {x[0],x[1]-3,x[2]+3,x[3]};
  calculateTheMomentumGiven(best_track,xminusplus,y,z,reco_pz,bestRegressionStats);
  float momminusplus=reco_pz;
    float xminusminus[4]= {x[0],x[1]-3,x[2]-3,x[3]};
  calculateTheMomentumGiven(best_track,xminusminus,y,z,reco_pz,bestRegressionStats);
  float momminusminus=reco_pz;
  
  calculateTheMomentum(best_track,x,y,z,reco_pz, reco_pz2M, bestRegressionStats, offset);
  Recodiff[89]->Fill(fourmom,-(fourmom-mom2minus)/fourmom);
  Recodiff[90]->Fill(fourmom,-(fourmom-mom2plus)/fourmom);
  Recodiff[91]->Fill(fourmom,-(fourmom-mom3minus)/fourmom);
  Recodiff[92]->Fill(fourmom,-(fourmom-mom3plus)/fourmom);
  Recodiff[93]->Fill(fourmom,-(fourmom-momplusplus)/fourmom);
  Recodiff[94]->Fill(fourmom,-(fourmom-momplusminus)/fourmom);
  Recodiff[95]->Fill(fourmom,-(fourmom-momminusplus)/fourmom);
  Recodiff[96]->Fill(fourmom,-(fourmom-momminusminus)/fourmom);
  float mom_error= CalculateTheMomentumError(fourx,foury,fourz,reco_pz);
  Recodiff[88]->Fill(reco_pz,mom_error/reco_pz);
  //Checking Doug Jensen's method of finding residual of WC2 and WC3 point to line through WC1 and WC4
  float slope_doug=(y[3]-y[0])/(z[3]-z[0]);
  float intercept_doug=y[0]-slope_doug*z[0];
  float restwo_doug=(y[1]-slope_doug*z[1]-intercept_doug)/(std::sqrt(1+slope_doug*slope_doug));
  float resthree_doug=(y[2]-slope_doug*z[2]-intercept_doug)/(std::sqrt(1+slope_doug*slope_doug));
  Recodiff[86]->Fill(reco_pz,restwo_doug);
  Recodiff[87]->Fill(reco_pz,resthree_doug);
  //Ending Doug's method.
  //float fourmom=reco_pz;
  float fourwires[8]={0,0,0,0,0,0,0,0};
  float fourtimes[8]={0,0,0,0,0,0,0,0};
  for(int i=0; i<8; ++i){
    fourwires[i]=best_track.hits[i].wire;
    fourtimes[i]=best_track.hits[i].time;
  }
  //std::cout<<"Momentum Calculated"<<std::endl;
  //reco_pz_list.push_back(reco_pz);
//We should also have the x,y,z points of the best_track, so now find where it hits the TPC
  projectToTPC(best_track,x,y,z,bestRegressionStats,x_face_list,y_face_list,incoming_theta_list,incoming_phi_list);
  double fourface[2]={x_face_list[0],y_face_list[0]};
  double fourangles[2]={incoming_phi_list[0],incoming_theta_list[0]};
  //std::cout<<"TPC Projected"<<std::endl;
  calculateTrackKink_Dists(x,y,z,bestRegressionStats,y_kink_list,x_dist_list,y_dist_list,z_dist_list);
  double fourdistlist[4]={x_dist_list[0],y_dist_list[0],z_dist_list[0],y_kink_list[0]};
  
  
  
  
  fNHits=3;
  WCMissed=2;
  
  for(int i=0; i<4; ++i){
    x[i]=0;
    y[i]=0;
    z[i]=0;
  }
  bestRegressionStats.clear();
  missed_wire_hits.clear();
 
  bestResSq=initialconst;
    for(int i=0;i<3;i++){
      bestRegressionStats.push_back(initialconst);
    }
    for(int i=0; i<2; ++i){
    missed_wire_hits.push_back(initialconst);
    }
  //std::cout<<"WCMissed"<<WCMissed<<std::endl;
  //Loop over all combinations of hits, and find the positions
  //Because I cannot figure a better way of doing this, we will do two loops, one if WC2 is missed, the other if WC3 is missed. The code will be the same, just filling an empty hit list in place of the missed WC. -GP
  WCHitList NullList;
  WCHit FakeHit;
  FakeHit.wire=initialconst;
  FakeHit.time=initialconst;
  FakeHit.hit_index=initialconst;
  FakeHit.cluster_index=initialconst;
  FakeHit.isVisited=false;
  NullList.hits.push_back(FakeHit);

  if(WCMissed==2){
    for( size_t iHit0 = 0; iHit0 < good_hits[0][0].hits.size(); ++iHit0 ){
      for( size_t iHit1 = 0; iHit1 < good_hits[0][1].hits.size(); ++iHit1 ){
        //for( size_t iHit2 = 0; iHit2 < good_hits[1][0].hits.size(); ++iHit2 ){
	  //for( size_t iHit3 = 0; iHit3 < good_hits[1][1].hits.size(); ++iHit3 ){
	    for( size_t iHit4 = 0; iHit4 < good_hits[2][0].hits.size(); ++iHit4 ){
	      for( size_t iHit5 = 0; iHit5 < good_hits[2][1].hits.size(); ++iHit5 ){
	        for( size_t iHit6 = 0; iHit6 < good_hits[3][0].hits.size(); ++iHit6 ){
		  for( size_t iHit7 = 0; iHit7 < good_hits[3][1].hits.size(); ++iHit7 ){
		    WCHitList track;
		    track.hits.push_back(good_hits[0][0].hits[iHit0]);
		    track.hits.push_back(good_hits[0][1].hits[iHit1]);
		    //std::cout<<"faking the hit"<<std::endl;
		    track.hits.push_back(NullList.hits[0]);
		    track.hits.push_back(NullList.hits[0]);
		    //std::cout<<"Faked Hit"<<std::endl;
		    //track.hits.push_back(good_hits[1][0].hits[iHit2]);
		    //track.hits.push_back(good_hits[1][1].hits[iHit3]);
		    track.hits.push_back(good_hits[2][0].hits[iHit4]);
		    track.hits.push_back(good_hits[2][1].hits[iHit5]);
		    track.hits.push_back(good_hits[3][0].hits[iHit6]);
		    track.hits.push_back(good_hits[3][1].hits[iHit7]);
		    findTheHitPositions(track,x,y,z,WCMissed);
		    //std::cout<<"Three Point Hit position found"<<std::endl;
		    //Do regression on the four points to find the straightest track in Y
		    std::vector<float> track_stats=Regression(y,z,WCMissed);
		    //std::cout<<"Track Res Sq :"<<track_stats[2]<<std::endl;
		    //std::cout<<"Three Regressed"<<std::endl;
		    if(track_stats[2]<fabs(bestResSq)){
		      best_track=track;
		      bestRegressionStats=track_stats;
		    }
		    //std::cout<<"Three track checked for best track"<<std::endl;
		  }
	        }	
	      }	
	    }	
	  //}	
        //}		
      }		
    }
  }
  //std::cout<<"Best track residual"<<bestRegressionStats[2]<<std::endl;
  float twores=bestRegressionStats[2];
  Recodiff[59]->Fill(bestRegressionStats[2],bestRegressionStats[2]);
  float reco_pz_three=0;  
//Now we should have the straightest track in Y, which will be the track that goes to the event.  Now we get the momentum and projections onto the TPC  
  calculateTheThreePointMomentum(best_track,x,y,z,reco_pz_three,bestRegressionStats,WCMissed);
  float twomom=reco_pz_three;
 // std::cout<<"Three Momentum calculated"<<std::endl;
  //reco_pz_list.push_back(reco_pz_three);
  //Recodiff[0]->Fill(reco_four,reco_pz_three);//Now that the momentum is scaled back to what we would get with four points, try to find where the missed point would be to reconstruct this scaled momentum
  extrapolateTheMissedPoint(best_track,x,y,z,reco_pz_three,bestRegressionStats, missed_wire_hits,WCMissed);
  //std::cout<<"Missed Point Extrapolated"<<std::endl;
  //std::cout<<"Missed X wire: "<<missed_wire_hits[0]<<" Missed Y Wire : "<<missed_wire_hits[1]<<std::endl;
  //MissedWireHits[0]->Fill(missed_wire_hits[0]);
  //MissedWireHits[1]->Fill(missed_wire_hits[1]);
//Put the extrapolated hit onto the track
  if(WCMissed==2){
    best_track.hits[2].wire=missed_wire_hits[0];
    best_track.hits[3].wire=missed_wire_hits[1];
  }
  if(WCMissed==3){
    best_track.hits[4].wire=missed_wire_hits[0];
    best_track.hits[5].wire=missed_wire_hits[1];
  }
  if(WCMissed==4){
    best_track.hits[6].wire=missed_wire_hits[0];
    best_track.hits[7].wire=missed_wire_hits[1];
  } 
  x_face_list.clear();
  y_face_list.clear();
  incoming_phi_list.clear();
  incoming_theta_list.clear();
  x_dist_list.clear();
  y_dist_list.clear();
  z_dist_list.clear();
  y_kink_list.clear();
  calculateTrackKink_Dists(x,y,z,bestRegressionStats,y_kink_list,x_dist_list,y_dist_list,z_dist_list);
  projectToTPC(best_track,x,y,z,bestRegressionStats,x_face_list,y_face_list,incoming_theta_list,incoming_phi_list);
  float twox[4]={x[0],x[1],x[2],x[3]};
  float twoy[4]={y[0],y[1],y[2],y[3]};
  float twoz[4]={z[0],z[1],z[2],z[3]};
  float twowires[8]={0,0,0,0,0,0,0,0};
  for(int i=0; i<8; ++i){
    twowires[i]=best_track.hits[i].wire;
  }
  double twoface[2]={x_face_list[0],y_face_list[0]};
  double twoangles[2]={incoming_phi_list[0],incoming_theta_list[0]}; 
  double twodistlist[4]={x_dist_list[0],y_dist_list[0],z_dist_list[0],y_kink_list[0]};
  
  for(int i=0; i<4; ++i){
    x[i]=0;
    y[i]=0;
    z[i]=0;
  }
   bestRegressionStats.clear();
   missed_wire_hits.clear();
  bestResSq=initialconst;
    for(int i=0;i<3;i++){
      bestRegressionStats.push_back(initialconst);
    }
    for(int i=0; i<2; ++i){
      missed_wire_hits.push_back(initialconst);
    }
  //std::cout<<"WCMissed"<<WCMissed<<std::endl;
  //Loop over all combinations of hits, and find the positions
  //Because I cannot figure a better way of doing this, we will do two loops, one if WC2 is missed, the other if WC3 is missed. The code will be the same, just filling an empty hit list in place of the missed WC. -GP

  
    
WCMissed=3;  
  if(WCMissed==3){
    for( size_t iHit0 = 0; iHit0 < good_hits[0][0].hits.size(); ++iHit0 ){
      for( size_t iHit1 = 0; iHit1 < good_hits[0][1].hits.size(); ++iHit1 ){
        for( size_t iHit2 = 0; iHit2 < good_hits[1][0].hits.size(); ++iHit2 ){
	  for( size_t iHit3 = 0; iHit3 < good_hits[1][1].hits.size(); ++iHit3 ){
	    //for( size_t iHit4 = 0; iHit4 < good_hits[2][0].hits.size(); ++iHit4 ){
	      //for( size_t iHit5 = 0; iHit5 < good_hits[2][1].hits.size(); ++iHit5 ){
	        for( size_t iHit6 = 0; iHit6 < good_hits[3][0].hits.size(); ++iHit6 ){
		  for( size_t iHit7 = 0; iHit7 < good_hits[3][1].hits.size(); ++iHit7 ){
		    WCHitList track;
		    track.hits.push_back(good_hits[0][0].hits[iHit0]);
		    track.hits.push_back(good_hits[0][1].hits[iHit1]);
		    track.hits.push_back(good_hits[1][0].hits[iHit2]);
		    track.hits.push_back(good_hits[1][1].hits[iHit3]);
		    track.hits.push_back(NullList.hits[0]);
		    track.hits.push_back(NullList.hits[0]);
		    //track.hits.push_back(good_hits[2][0].hits[iHit4]);
		    //track.hits.push_back(good_hits[2][1].hits[iHit5]);
		    track.hits.push_back(good_hits[3][0].hits[iHit6]);
		    track.hits.push_back(good_hits[3][1].hits[iHit7]);
		    findTheHitPositions(track,x,y,z,WCMissed);
		    //Do regression on the four points to find the straightest track in Y
		    std::vector<float> track_stats=Regression(y,z,WCMissed);
		     //std::cout<<"Track Res Sq :"<<track_stats[2]<<std::endl;
		    if(track_stats[2]<fabs(bestResSq)){
		      best_track=track;
		      bestRegressionStats=track_stats;
		      //std::cout<<track_stats[2]<<std::endl;
		    }
		  }
	        }	
	      //}	
	    //}	
	  }	
        }		
      }		
    }
  }
 Recodiff[60]->Fill(bestRegressionStats[2],bestRegressionStats[2]);
  float threeres=bestRegressionStats[2];
 reco_pz_three=0;  
//Now we should have the straightest track in Y, which will be the track that goes to the event.  Now we get the momentum and projections onto the TPC  
  calculateTheThreePointMomentum(best_track,x,y,z,reco_pz_three,bestRegressionStats,WCMissed);
  float threemom=reco_pz_three;
 // std::cout<<"Three Momentum calculated"<<std::endl;
  //reco_pz_list.push_back(reco_pz_three);
  //Recodiff[0]->Fill(reco_four,reco_pz_three);
//Now that the momentum is scaled back to what we would get with four points, try to find where the missed point would be to reconstruct this scaled momentum
  extrapolateTheMissedPoint(best_track,x,y,z,reco_pz_three,bestRegressionStats, missed_wire_hits,WCMissed);
  //std::cout<<"Missed Point Extrapolated"<<std::endl;
  //std::cout<<"Missed X wire: "<<missed_wire_hits[0]<<" Missed Y Wire : "<<missed_wire_hits[1]<<std::endl;
//Put the extrapolated hit onto the track
  if(WCMissed==2){
    best_track.hits[2].wire=missed_wire_hits[0];
    best_track.hits[3].wire=missed_wire_hits[1];
  }
  if(WCMissed==3){
    best_track.hits[4].wire=missed_wire_hits[0];
    best_track.hits[5].wire=missed_wire_hits[1];
  }
  if(WCMissed==4){
    best_track.hits[6].wire=missed_wire_hits[0];
    best_track.hits[7].wire=missed_wire_hits[1];
  }   
//We should also have the x,y,z points of the best_track, so now find where it hits the TPC
  x_face_list.clear();
  y_face_list.clear();
  incoming_phi_list.clear();
  incoming_theta_list.clear();
  x_dist_list.clear();
  y_dist_list.clear();
  z_dist_list.clear();
  y_kink_list.clear();
   calculateTrackKink_Dists(x,y,z,bestRegressionStats,y_kink_list,x_dist_list,y_dist_list,z_dist_list);
  projectToTPC(best_track,x,y,z,bestRegressionStats,x_face_list,y_face_list,incoming_theta_list,incoming_phi_list);  
  float threex[4]={x[0],x[1],x[2],x[3]};
  float threey[4]={y[0],y[1],y[2],y[3]};
  float threez[4]={z[0],z[1],z[2],z[3]};
  float threewires[8]={0,0,0,0,0,0,0,0};
  for(int i=0; i<8; ++i){
    threewires[i]=best_track.hits[i].wire;
  }
  double threeface[2]={x_face_list[0],y_face_list[0]};
  double threeangles[2]={incoming_phi_list[0],incoming_theta_list[0]}; 
  double threedistlist[4]={x_dist_list[0],y_dist_list[0],z_dist_list[0],y_kink_list[0]};
  
for(int i=0; i<4; ++i){
    x[i]=0;
    y[i]=0;
    z[i]=0;
  }
   bestRegressionStats.clear();
   missed_wire_hits.clear();
  bestResSq=initialconst;
    for(int i=0;i<3;i++){
      bestRegressionStats.push_back(initialconst);
    }
    for(int i=0; i<2; ++i){
      missed_wire_hits.push_back(initialconst);
    }
  //std::cout<<"WCMissed"<<WCMissed<<std::endl;  
 WCMissed=4;
if(WCMissed==4){
    for( size_t iHit0 = 0; iHit0 < good_hits[0][0].hits.size(); ++iHit0 ){
      for( size_t iHit1 = 0; iHit1 < good_hits[0][1].hits.size(); ++iHit1 ){
        for( size_t iHit2 = 0; iHit2 < good_hits[1][0].hits.size(); ++iHit2 ){
	  for( size_t iHit3 = 0; iHit3 < good_hits[1][1].hits.size(); ++iHit3 ){
	    for( size_t iHit4 = 0; iHit4 < good_hits[2][0].hits.size(); ++iHit4 ){
	      for( size_t iHit5 = 0; iHit5 < good_hits[2][1].hits.size(); ++iHit5 ){
	        //for( size_t iHit6 = 0; iHit6 < good_hits[3][0].hits.size(); ++iHit6 ){
		  //for( size_t iHit7 = 0; iHit7 < good_hits[3][1].hits.size(); ++iHit7 ){
		    WCHitList track;
		    track.hits.push_back(good_hits[0][0].hits[iHit0]);
		    track.hits.push_back(good_hits[0][1].hits[iHit1]);
		    track.hits.push_back(good_hits[1][0].hits[iHit2]);
		    track.hits.push_back(good_hits[1][1].hits[iHit3]);
		    //track.hits.push_back(NullList.hits[0]);
		    //track.hits.push_back(NullList.hits[0]);
		    track.hits.push_back(good_hits[2][0].hits[iHit4]);
		    track.hits.push_back(good_hits[2][1].hits[iHit5]);
		    track.hits.push_back(NullList.hits[0]);
		    track.hits.push_back(NullList.hits[0]);
		    //track.hits.push_back(good_hits[3][0].hits[iHit6]);
		    //track.hits.push_back(good_hits[3][1].hits[iHit7]);
		    findTheHitPositions(track,x,y,z,WCMissed);
		    //Do regression on the four points to find the straightest track in Y
		    std::vector<float> track_stats=Regression(y,z,WCMissed);
		    //std::cout<<"Track Res Sq :"<<track_stats[2]<<std::endl;
		    if(track_stats[2]<fabs(bestResSq)){
		      best_track=track;
		      bestRegressionStats=track_stats;
		      bestResSq=track_stats[2];
		      //std::cout<<track_stats[2]<<std::endl;
		    }
		  //}
	        //}	
	      }	
	    }	
	  }	
        }		
      }		
    }
  } 
 Recodiff[61]->Fill(bestRegressionStats[2],bestRegressionStats[2]);
  float fourmres=bestRegressionStats[2];
 reco_pz_three=0;  
//Now we should have the straightest track in Y, which will be the track that goes to the event.  Now we get the momentum and projections onto the TPC  
  calculateTheThreePointMomentum(best_track,x,y,z,reco_pz_three,bestRegressionStats,WCMissed);
  float fourmmom=reco_pz_three;
 // std::cout<<"Three Momentum calculated"<<std::endl;
  //reco_pz_list.push_back(reco_pz_three);
  //Recodiff[0]->Fill(reco_four,reco_pz_three);
//Now that the momentum is scaled back to what we would get with four points, try to find where the missed point would be to reconstruct this scaled momentum
  extrapolateTheMissedPoint(best_track,x,y,z,reco_pz_three,bestRegressionStats, missed_wire_hits,WCMissed);
  //std::cout<<"Missed Point Extrapolated"<<std::endl;
  //std::cout<<"Missed X wire: "<<missed_wire_hits[0]<<" Missed Y Wire : "<<missed_wire_hits[1]<<std::endl;
//Put the extrapolated hit onto the track
  if(WCMissed==2){
    best_track.hits[2].wire=missed_wire_hits[0];
    best_track.hits[3].wire=missed_wire_hits[1];
  }
  if(WCMissed==3){
    best_track.hits[4].wire=missed_wire_hits[0];
    best_track.hits[5].wire=missed_wire_hits[1];
  }
  if(WCMissed==4){
    best_track.hits[6].wire=missed_wire_hits[0];
    best_track.hits[7].wire=missed_wire_hits[1];
  }   
//We should also have the x,y,z points of the best_track, so now find where it hits the TPC
  x_face_list.clear();
  y_face_list.clear();
  incoming_phi_list.clear();
  incoming_theta_list.clear();
  x_dist_list.clear();
  y_dist_list.clear();
  z_dist_list.clear();
  y_kink_list.clear();
   calculateTrackKink_Dists(x,y,z,bestRegressionStats,y_kink_list,x_dist_list,y_dist_list,z_dist_list);
  projectToTPC(best_track,x,y,z,bestRegressionStats,x_face_list,y_face_list,incoming_theta_list,incoming_phi_list);  
  float fourmx[4]={x[0],x[1],x[2],x[3]};
  float fourmy[4]={y[0],y[1],y[2],y[3]};
  float fourmz[4]={z[0],z[1],z[2],z[3]};
  float fourmwires[8]={0,0,0,0,0,0,0,0};
  for(int i=0; i<8; ++i){
    fourmwires[i]=best_track.hits[i].wire;
  }
  double fourmface[2]={x_face_list[0],y_face_list[0]};
  double fourmangles[2]={incoming_phi_list[0],incoming_theta_list[0]}; 
  double fourmdistlist[4]={x_dist_list[0],y_dist_list[0],z_dist_list[0],y_kink_list[0]};
  for(int i=0; i<8; ++i){
    Recodiff[i]->Fill(fourwires[i],twowires[i]);
    Recodiff[i+8]->Fill(fourwires[i],threewires[i]);
  
  }
  for(int i=0; i<4; ++i){
    Recodiff[i+16]->Fill(fourx[i],twox[i]);
    Recodiff[i+20]->Fill(fourx[i],threex[i]);
    Recodiff[i+24]->Fill(foury[i],twoy[i]);
    Recodiff[i+28]->Fill(foury[i],threey[i]);
    Recodiff[i+32]->Fill(fourz[i],twoz[i]);
    Recodiff[i+36]->Fill(fourz[i],threez[i]);    
    Recodiff[i+48]->Fill(fourdistlist[i],twodistlist[i]);
    Recodiff[i+52]->Fill(fourdistlist[i],threedistlist[i]);
    Recodiff[i+66]->Fill(fourdistlist[i],fourmdistlist[i]);

    
  }
  for(int i=0; i<2; ++i){
    Recodiff[i+40]->Fill(fourface[i],twoface[i]);
    Recodiff[i+42]->Fill(fourface[i],threeface[i]);
    Recodiff[i+44]->Fill(fourangles[i],twoangles[i]);
    Recodiff[i+46]->Fill(fourangles[i],threeangles[i]);
    Recodiff[i+70]->Fill(fourangles[i],fourmangles[i]);
    Recodiff[i+72]->Fill(fourwires[i+6],fourmwires[i+6]);
    Recodiff[i+74]->Fill(fourface[0]-fourmface[0],fourface[1]-fourmface[1]);
  }
  Recodiff[56]->Fill(fourmom,(twomom-fourmom)/fourmom);
  Recodiff[57]->Fill(fourmom,(threemom-fourmom)/fourmom);
  Recodiff[62]->Fill(fourmom,(fourmmom-fourmom)/fourmom);
  Recodiff[63]->Fill(fourx[3],fourmx[3]);
  Recodiff[64]->Fill(foury[3],fourmy[3]);
  Recodiff[65]->Fill(fourz[3],fourmz[3]);
  Recodiff[76]->Fill(fourres,twores);
  Recodiff[77]->Fill(fourres,threeres);
  Recodiff[78]->Fill(fourres,fourmres);
  Recodiff[79]->Fill(fourwires[4]-threewires[4],fourwires[4]-threewires[4]);
  Recodiff[80]->Fill(fourmom,fourmom);
  Recodiff[81]->Fill(fourtimes[1]-fourtimes[0],fourtimes[1]-fourtimes[0]);
  Recodiff[82]->Fill(fourtimes[3]-fourtimes[2],fourtimes[3]-fourtimes[2]);
  Recodiff[83]->Fill(fourtimes[5]-fourtimes[4],fourtimes[5]-fourtimes[4]);
  Recodiff[84]->Fill(fourtimes[7]-fourtimes[6],fourtimes[7]-fourtimes[6]);
  Recodiff[85]->Fill(fourmom,fourres);
//}
}
//=====================================================================================
float WCTrackBuilderAlg::CalculateTheMomentumError(float (&x)[4],
						   float (&y)[4],
						   float (&z)[4],
						   float & reco_pz) 
{
  float dB=fB_field_tesla*.02;
  float dwire=10;
  float dxalign=25.4;
  float dzalign=25.4;
  float alpha=1/3.3;
  float theta_ds=atan((x[3]-x[2])/(z[3]-z[2]));
  float theta_us=atan((x[1]-x[0])/(z[1]-z[0]));
  float xerror_us=dxalign+dwire*cos(theta_us);
  float xerror_ds=dxalign+dwire*cos(theta_ds);
  float zerror_us=dzalign+dwire*sin(theta_us);
  float zerror_ds=dzalign+dwire*sin(theta_ds);
  float dz_ds=z[3]-z[2];
  float dz_us=z[1]-z[0];
  float Berror=pow(alpha*fL_eff*dB/(sin(theta_ds)-sin(theta_us)),2);
  float gamma_us=pow(alpha*fB_field_tesla*fL_eff/(sin(theta_ds)-sin(theta_us))/cos(theta_us)/dz_us,2);
  float gamma_ds=pow(alpha*fB_field_tesla*fL_eff/(sin(theta_ds)-sin(theta_us))/cos(theta_ds)/dz_ds,2);
  float error_us=gamma_us*(2*pow(xerror_us,2)+pow(theta_us,2)*(2*pow(zerror_us,2)));
  float error_ds=gamma_ds*(2*pow(xerror_ds,2)+pow(theta_ds,2)*(2*pow(zerror_ds,2)));
  float error_mom=pow(Berror+error_us+error_ds,.5);
  std::cout<<theta_us<<" "<<theta_ds<<" "<<xerror_us<<" "<<xerror_ds<<" "<<zerror_us<<" "<<zerror_ds<<" "<<Berror<<" "<<error_us<<" "<<error_ds<<" "<<error_mom<<" "<<reco_pz<<" "<<error_mom/reco_pz<<std::endl;
  return error_mom;
}
//========================================================================================
void WCTrackBuilderAlg::calculateTheMomentumGiven(WCHitList & best_track,
					     float (&x)[4],
					     float (&y)[4],
					     float (&z)[4],
					     float & reco_pz,
					     std::vector<float> & BestTrackStats)
{
//We need the x,y,z again, which would be for the last track combination tried. So we need to find the positions again
  //findTheHitPositions(best_track,x,y,z,WCMissed);
  //std::cout<<"Best track hit position found"<<std::endl;
//Calculate the angle of the track, in the x,z and y,z planes, in upstream(us) and downstream(ds) ends of the WC
  float dx_us=x[1]-x[0];
  float dx_ds=x[3]-x[2];
  //float dy_us=y[1]-y[0];
  //float dy_ds=y[3]-y[2];
  float dz_us=z[1]-z[0];
  float dz_ds=z[3]-z[2];
  float theta_x_us= atan(dx_us/dz_us);
  float theta_x_ds= atan(dx_ds/dz_ds);
  //float theta_y_us= atan(dy_us/dz_us);
  //float theta_y_ds= atan(dy_ds/dz_ds);
  //reco_pz = (fabs(fB_field_tesla) * fL_eff * fmm_to_m * fGeV_to_MeV ) / (3.3*(sin(theta_x_ds) - sin(theta_x_us))) / cos(atan(BestTrackStats[0]));
  reco_pz = calculateRecoPz(theta_x_us,theta_x_ds,BestTrackStats[0]);
  
}  
//=======================================================================================
TVector3 WCTrackBuilderAlg::PlotTheMidplane(float (&x)[4],
					float (&y)[4],
					float (&z)[4],
					float dist)
{

TVector3 v0(x[0],y[0],z[0]);
TVector3 v1(x[1],y[1],z[1]);
TVector3 v2(x[2],y[2],z[2]);
TVector3 v3(x[3],y[3],z[3]);
TVector3 v10=v1-v0;
TVector3 v23=v2-v3;
TVector3 v03=v0-v3;
TVector3 v10_unit=v10.Unit();
TVector3 v23_unit=v23.Unit();
double b=v10_unit.Dot(v23_unit);
double d=v10_unit.Dot(v03);
double e=v23_unit.Dot(v03);
double s_us=(b*e-d)/(1-b*b);
double s_ds=(e-b*d)/(1-b*b);
TVector3 closest_dist_vect=v03+((b*e-d)*v10_unit-(e-b*d)*v23_unit)*(1/(1-b*b));
float dist_temp=float(closest_dist_vect.Mag());
dist=dist_temp;
TVector3 v23_close=v3+s_ds*v23_unit; //Q(t_c)
TVector3 closest_vect=v03+s_us*v10_unit-s_ds*v23_unit; //w_c
TVector3 midpointvect=v23_close + 0.5*closest_vect;
return midpointvect;
//std::cout<<"midpoint: ["<<x_mid<<", "<<y_mid<<", "<<z_mid<<"]"<<std::endl;
//std::cout<<"distance of closest approach: "<<dist<<std::endl;
}
//-===========================================================
float WCTrackBuilderAlg::PlotTheMidplane(float (&x)[4],
					float (&y)[4],
					float (&z)[4])
{

TVector3 v0(x[0],y[0],z[0]);
TVector3 v1(x[1],y[1],z[1]);
TVector3 v2(x[2],y[2],z[2]);
TVector3 v3(x[3],y[3],z[3]);
TVector3 v10=v1-v0;
TVector3 v23=v2-v3;
TVector3 v03=v0-v3;
TVector3 v10_unit=v10.Unit();
TVector3 v23_unit=v23.Unit();
double b=v10_unit.Dot(v23_unit);
double d=v10_unit.Dot(v03);
double e=v23_unit.Dot(v03);
//double s_us=(b*e-d)/(1-b*b);
//double s_ds=(e-b*d)/(1-b*b);
TVector3 closest_dist_vect=v03+((b*e-d)*v10_unit-(e-b*d)*v23_unit)*(1/(1-b*b));
float dist_temp=float(closest_dist_vect.Mag());
float dist=dist_temp;
//TVector3 v23_close=v3+s_ds*v23_unit;
//TVector3 closest_vect=v03+s_us*v10_unit-s_ds*v23_unit;
//TVector3 midpointvect=v3+s_ds*v23_unit+.5*closest_vect;
return dist;
//std::cout<<"midpoint: ["<<x_mid<<", "<<y_mid<<", "<<z_mid<<"]"<<std::endl;
//std::cout<<"distance of closest approach: "<<dist<<std::endl;
}					
