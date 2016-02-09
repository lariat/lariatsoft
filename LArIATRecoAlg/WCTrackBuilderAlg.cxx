////////////////////////////////////////////////////////////////////////
// This Alg takes in the good hits from WCHitFinderAlg and creates    //
// a track with various kinematic and geometric values.  In effect,   //
// this alg is the track building part of WCTrackBuilderAlg.cxx       //
// in functions run after "FinalizeGoodHits".                         //
// Other algs can be created in lieu of this one if a better track    //
// is developed.  This is merely the version we used when hit finding //
// and track builder were run together in WCTrackBuilderAlg.cxx       //
// Author: Greg Pulliam gkpullia@syr.edu                              //
//////////////////////////////////////////////////////////////////////// 

//Framework includes
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "CLHEP/Units/SystemOfUnits.h"
#include "Geometry/AuxDetGeo.h"
#include "art/Framework/Services/Optional/TFileService.h"

// LArIAT includes
#include "LArIATRecoAlg/WCTrackBuilderAlg.h"


#include <iostream>
#include <cmath>
#include <cstdlib>
#include <TH1F.h>
#include <string>
#include <TH2F.h>

WCTrackBuilderAlg::WCTrackBuilderAlg( fhicl::ParameterSet const& pset )
{
  this->reconfigure(pset);

  //Testing the AuxDetGeo capabilitites
  std::vector<geo::AuxDetGeo*> const & theAuxDetGeoVect = fGeo->AuxDetGeoVec();
  double centerOfDet[3] = {0,0,0};
  for( size_t iDet = 0; iDet < fGeo->NAuxDets() ; ++iDet ){
    geo::AuxDetGeo* anAuxDetGeo = theAuxDetGeoVect.at(iDet);
    anAuxDetGeo->GetCenter(centerOfDet);

    //Setting the TGeo world locations of the MWPCs in mm
    if(iDet == 1){ //WC1
      fX_cntr[0] = centerOfDet[0] * CLHEP::cm;
      fY_cntr[0] = centerOfDet[1] * CLHEP::cm;
      fZ_cntr[0] = centerOfDet[2] * CLHEP::cm;
    }
   if(iDet == 2){ //WC2
      fX_cntr[1] = centerOfDet[0] * CLHEP::cm;
      fY_cntr[1] = centerOfDet[1] * CLHEP::cm;
      fZ_cntr[1] = centerOfDet[2] * CLHEP::cm;
    }
   if(iDet == 3){ //WC3
      fX_cntr[2] = centerOfDet[0] * CLHEP::cm;
      fY_cntr[2] = centerOfDet[1] * CLHEP::cm;
      fZ_cntr[2] = centerOfDet[2] * CLHEP::cm;
    }
   if(iDet == 4){ //WC4
      fX_cntr[3] = centerOfDet[0] * CLHEP::cm;
      fY_cntr[3] = centerOfDet[1] * CLHEP::cm;
      fZ_cntr[3] = centerOfDet[2] * CLHEP::cm;
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

}
//------------------------------------------------------------------------------
WCTrackBuilderAlg::~WCTrackBuilderAlg()
{

}
//------------------------------------------------------------------------------
void WCTrackBuilderAlg::reconfigure( fhicl::ParameterSet const& pset )
{

  fB_field_tesla        = pset.get<float >("BFieldInTesla",      0.       );


  //fCentralYKink         = pset.get<float >("CentralYKink",        -0.01    ); //These four are parameters from histos I produced from picky-good tracks
  //fSigmaYKink           = pset.get<float >("SigmaYKink",          0.03      );
  //fCentralYDist         = pset.get<float >("CentralYDist",        0.69      );
  //fSigmaYDist           = pset.get<float >("SigmaYDist",          18.0      );

 // fPrintDisambiguation = false;
  fPickyTracks          = pset.get<bool  >("PickyTracks",         false     );
  fHighYield            = pset.get<bool  >("HighYield",           false     );
  fDiagnostics          = pset.get<bool  >("Diagnostics",         false     );
  
  //Survey constants
  //fDelta_z_us           = pset.get<float >("DeltaZus",            1551.15   );  //this will recalculated using geometry instead of hardcoding. 
  //fDelta_z_ds 		= pset.get<float >("DeltaZds",   	  1570.06   );   
  fL_eff        	= pset.get<float >("LEffective", 	  1145.34706);	
  fmm_to_m    		= pset.get<float >("MMtoM",        	  0.001     );	
  fGeV_to_MeV 		= pset.get<float >("GeVToMeV",    	  1000.0    );   
  fMP_X                 = pset.get<float> ("MidplaneInterceptFactor", 1);
  fMP_M                 = pset.get<float> ("MidplaneSlopeFactor", 1);
  fMidplane_intercept   = pset.get<float> ("MidplaneIntercept", 1027.65); 
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
  current=std::stod(fDatabaseUtility->GetIFBeamValue("mid_f_mc7an",fRun,fSubRun));
  if(fabs(current)>90){fB_field_tesla= .003375*current;}
  if(fabs(current)<90 && fabs(current)>70){fB_field_tesla= .0034875*current;}
  if(fabs(current)<70 && fabs(current)>50){fB_field_tesla= .003525*current;}
  if(fabs(current)<50 && fabs(current)>30){fB_field_tesla= .003525*current;}  
  if(fabs(current)<30){fB_field_tesla= .0035375*current;}  
  std::cout << "Run: " << fRun << ", Subrun: " << fSubRun << ", B-field: " << fB_field_tesla << std::endl;
}
//--------------------------------------------------------------
//Main function called for each trigger
void WCTrackBuilderAlg::reconstructTracks(std::vector<double> & reco_pz_list,               
					     std::vector<double> & x_face_list,
					     std::vector<double> & y_face_list,
					     std::vector<double> & incoming_theta_list,
					     std::vector<double> & incoming_phi_list,
					     std::vector<WCHitList> & event_final_tracks,
					     std::vector<std::vector<WCHitList> > & good_hits,
					     bool pickytracks,
					     bool highyield,
					     bool diagnostics, 
					     std::vector<double> & y_kink_list,
					     std::vector<double> & x_dist_list,
					     std::vector<double> & y_dist_list,
					     std::vector<double> & z_dist_list,
					     int & WCMissed,
					     std::vector<TH2F*>  & Recodiff)
{					   
  fPickyTracks = pickytracks;
  fHighYield = highyield;
  fDiagnostics= diagnostics;
  std::cout<<"PickyTracks : "<<fPickyTracks<<"High Yield : "<<fHighYield<<std::endl;
  initialconst=-999999999;  //Just a number to use to initialize things before they get filled correctly.
  WCMissed=initialconst;  					 	
  //Determine if one should skip this trigger based on whether there is exactly one good hit in each wire chamber and axis
  //If there isn't, continue after adding a new empty vector to the reco_pz_array contianer.
  //If there is, then move on with the function.
  //This can be modified to permit more than one good hit in each wire chamber axis - see comments in function
  bool skip = shouldSkipTrigger(good_hits,WCMissed);
  //std::cout<<"should skip trigger done"<<std::endl;
  std::cout<<"Hits : "<<NHits<<"WC Missed : "<<WCMissed<<std::endl;
  if( skip == true ) return;
  
  //Depending on if an event has a hit in all 4 WC or whether it missed WC2 or WC3 (but not both), we reconstruct the momentum differently. This code doesn't change from before we allowed 3 point tracks.
  if(NHits==4){
  
  //At this point, we should have a list of good hits with at least one good hit in X,Y for each WC.
  //Now find all possible combinations of these hits that could form a track, sans consideration
  //of the kinks or end displacements (for now). For the "exactly one" condition set in the above
  //step, this won't matter, but if you want to set the condition to "at least one hit in each WC axis,
  //this will give many combinations.
    buildFourPointTracks(good_hits,
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
			 WCMissed);
					   
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
  if(NHits==3)
  {
    buildThreePointTracks(good_hits,
    		          reco_pz_list,
			  x_face_list,
			  y_face_list,
			  incoming_phi_list,
			  incoming_theta_list,
			  event_final_tracks,
			  y_kink_list,
			  x_dist_list,
			  y_dist_list,
			  z_dist_list,
			  WCMissed);
  //std::cout<<"Build three point track"<<std::endl;
  }
  //To compare four point tracks to three point tracks
  if(fDiagnostics==true){
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
					  int & WCMissed)
{
  //Now determine if we want to skip
  bool skip = false;
  NHits=0;
  for( size_t iWC = 0; iWC < 4 ; ++iWC ){
  if(good_hits[iWC][0].hits.size()>0 && good_hits[iWC][1].hits.size()>0){++NHits;}
  else{WCMissed=iWC+1;}
  }
  //If we don't have 3 or 4 hits, skip.
  if(NHits<3){
    skip = true;
  } 
  if(fPickyTracks && fHighYield){
    for(size_t iWC=0; iWC<4; ++iWC){
      for(size_t iAX=0; iAX<2; ++iAX){
        if(good_hits[iWC][iAX].hits.size() > 1 || WCMissed==1 || WCMissed==4){  
	//skip events that have more than 1 hit in an axis, or the missed WC is the first or last WC
          skip = true;
	  break;
	}
      }
    }   
  }
  if(!fPickyTracks && fHighYield){
    if(WCMissed==1 || WCMissed==4){ //skip events with less than 3 X/Y hits or is missing the first or last WC
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
//===================================================================================
void WCTrackBuilderAlg::buildFourPointTracks(std::vector<std::vector<WCHitList> > & good_hits,
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
			      int & WCMissed)
{
  float x[4]{0,0,0,0};
  float y[4]{0,0,0,0};
  float z[4]{0,0,0,0};
  float reco_pz=0;
  WCHitList best_track;
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
		  std::cout<<"Regression value : "<<track_stats[2]<<std::endl;
		  if(track_stats[2]<fabs(bestResSq)){
		    best_track=track;
		    bestRegressionStats=track_stats;
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
  std::cout<<bestResSq<<std::endl;
//Now we should have the straightest track in Y, which will be the track that goes to the event.  Now we get the momentum and projections onto the TPC  
  calculateTheMomentum(best_track,x,y,z,reco_pz,bestRegressionStats);
  //std::cout<<"Momentum Calculated"<<std::endl;
  reco_pz_list.push_back(reco_pz);
//We should also have the x,y,z points of the best_track, so now find where it hits the TPC
  projectToTPC(best_track,x,y,z,bestRegressionStats,x_face_list,y_face_list,incoming_theta_list,incoming_phi_list);
  //std::cout<<"TPC Projected"<<std::endl;
  calculateTrackKink_Dists(x,y,z,bestRegressionStats,y_kink_list,x_dist_list,y_dist_list,z_dist_list);
  //std::cout<<"Track Kink/Dist calculated"<<std::endl;
  event_final_tracks.push_back(best_track);
  
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
  std::cout<<"Track Size: "<<track.hits.size()<<std::endl;
  for(size_t iHit = 0; iHit < track.hits.size() ; ++iHit ){
    int var = iHit % 2;
    int jHit= iHit;
    if( var == 0  && int (jHit != 2*(WCMissed-1)+var) ){ // Skip iHit for the missed WC 
      x_wires.push_back(track.hits.at(iHit).wire);
    }
    if( var == 1 && int(jHit != 2*(WCMissed-1)+var)){
      y_wires.push_back(track.hits.at(iHit).wire);
    }   
  }
  //std::cout<<"Wires set"<<std::endl;
  float sin13=sin((3.141592654/180)*13.0);
  float sin3=sin((3.141592654/180)*3.0);
  float cos3=cos((3.141592654/180)*3.0);
  float cos13=cos((3.141592654/180)*13.0);
  float cosangle[4]={cos13,cos13,cos3,cos3};
  float sinangle[4]={sin13,sin13,sin3,sin3};
  for(int iWC=0; iWC<4; ++iWC){
    if(iWC != WCMissed-1){
  // Find the position of the hit, correcting for the rotation of the WC of 13 degrees for US or 3 degrees for DS
      x[iWC] =  fX_cntr[iWC] + x_wires[iWC] * cosangle[iWC];
      y[iWC] =  fY_cntr[iWC] + y_wires[iWC];
      z[iWC] =  fZ_cntr[iWC] + x_wires[iWC] * sinangle[iWC];
      std::cout<<"For WC: "<<iWC+1<<"The position is set"<<std::endl;
    }
    if(iWC==WCMissed-1){ //Just making sure a value is set, should never be used.
    std::cout<<"Faking position for missed hit now"<<std::endl;
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
    residualsquare += residual*residual/Npoints; 
    }
  }
  RegressionValues.push_back(slope);
  RegressionValues.push_back(intercept);
  RegressionValues.push_back(residualsquare);
  return RegressionValues;
} 
//================================================================================    
void WCTrackBuilderAlg::calculateTheMomentum(WCHitList & best_track,
					     float (&x)[4],
					     float (&y)[4],
					     float (&z)[4],
					     float & reco_pz,
					     std::vector<float> & BestTrackStats)
{
//We need the x,y,z again, which would be for the last track combination tried. So we need to find the positions again
  findTheHitPositions(best_track,x,y,z,WCMissed);
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
  reco_pz = (fabs(fB_field_tesla) * fL_eff * fmm_to_m * fGeV_to_MeV ) / (3.3*(theta_x_ds - theta_x_us))/cos(atan(BestTrackStats[0]));
  
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
//The intercept of the regression is where the track hits at z=0, which is the TPC face, so we get the y position on the TPC for free
  double y_face=bestRegressionStats[1];
  y_face_list.push_back(y_face);
//Though the regression probably gives the best extrapolation to the TPC, we can only use it for y projections. If we need x values, we have to use the positions ds, which could introduce error if we have either WC3,4 wrong
//For the x position at the TPC, we need to use the ds WC to extrapolate a line to the TPC (z==0)
  float dx_ds=x[3]-x[2];  
  float dz_ds=z[3]-z[2];
  float dy_ds=y[3]-y[2];
  float ds_slope=dx_ds/dz_ds;
  float ds_x_intercept=x[3]-ds_slope*z[3];
  double x_face=ds_x_intercept;
  x_face_list.push_back(x_face);			      
//We can use the atan2 to find phi (the angle in the XY plane) of the track, which will go CCW around the circle starting at (x,y)=(1,0). This will be in radians
  incoming_phi_list.push_back(atan2(dy_ds,dx_ds));
//The vector pointing from WC4(x4,y4,z4) to the X,Y face hit (Xface,Yface,0) can be rotated around the Z axis, which would trace out a cone.  The opening angle of that cone would be 2*theta.  So we can use the equation of
//a cone to find theta 

  double r=pow(pow(x_face-x[3],2)+pow(y_face-y[3],2),0.5);
//Unforunately, this angle will always be positive, even though we expect tracks to enter around -3 degrees to the z axis.  However, the previous iteration of WCTrackBuilder left theta positive, so I will too.  The negative
//accounts for the fact that I use z4 as the height, which is a negative number, which we need to make positive.
  double theta=atan(-r/z[3]);
  incoming_theta_list.push_back(theta);
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
   float z_mp_us=(fMidplane_intercept-x_us_int)/(x_us_slope-1/(tan(8.0*3.141592654/180)));  //X,Y,Z where US interesects Midplane.
   float x_mp_us=x_us_slope*z_mp_us+x_us_int;  
   float y_mp_us=y_us_slope*z_mp_us+y_us_int;
   float z_mp_ds=(fMidplane_intercept-x_ds_int)/(x_ds_slope-1/(tan(8.0*3.141592654/180)));  //X,Y,Z where DS interesects Midplane.
   float x_mp_ds=x_ds_slope*z_mp_ds+x_ds_int;
   float y_mp_ds=y_ds_slope*z_mp_ds+y_ds_int;
   x_dist_list.push_back(x_mp_ds-x_mp_us);
   y_dist_list.push_back(y_mp_ds-y_mp_us); 
   z_dist_list.push_back(z_mp_ds-z_mp_us);  
}
//====================================================================================
void WCTrackBuilderAlg::buildThreePointTracks(std::vector<std::vector<WCHitList> > & good_hits,
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
					      int & WCMissed)
{
//Code here is similar to the buildFourPointTrack version, but we have to allow for missed WC, and require some additional geometry to find the momentum
  float x[4]{0,0,0,0};
  float y[4]{0,0,0,0};
  float z[4]{0,0,0,0};
  std::vector<float> missed_wire_hits;
  std::vector<float> bestRegressionStats;
  WCHitList best_track;
  float bestResSq=initialconst;
    for(int i=0;i<3;i++){
      bestRegressionStats.push_back(initialconst);
    }
    for(int i=0; i<2; ++i){
      missed_wire_hits.push_back(initialconst);
    }
  std::cout<<"WCMissed"<<WCMissed<<std::endl;
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
		    std::cout<<"Track Res Sq :"<<track_stats[2]<<std::endl;
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
		     std::cout<<"Track Res Sq :"<<track_stats[2]<<std::endl;
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
  std::cout<<"Best track residual"<<bestRegressionStats[2]<<std::endl;
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
//We should also have the x,y,z points of the best_track, so now find where it hits the TPC
   calculateTrackKink_Dists(x,y,z,bestRegressionStats,y_kink_list,x_dist_list,y_dist_list,z_dist_list);
  projectToTPC(best_track,x,y,z,bestRegressionStats,x_face_list,y_face_list,incoming_theta_list,incoming_phi_list);
event_final_tracks.push_back(best_track);
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
    //std::cout<<"Midplane set"<<std::endl;
    float midplane_slope=tan((3.141592654/180)*8.0)*fMP_M;
    float midplane_intercept=fMidplane_intercept*fMP_X;
    float ds_dz=z[3]-z[2];
    float ds_dx=x[3]-x[2];
    float ds_slope=ds_dx/ds_dz;
    float ds_int_x= x[3]-ds_slope*z[3];
    float z_ds=(midplane_intercept-ds_int_x)/(ds_slope-midplane_slope); //Solving tan8*Z+Bmp == Mds*Z+Bds
    float x_ds=ds_slope*z_ds+ ds_int_x;
    float us_dz=(z_ds-z[0]);
    float us_dx=(x_ds-x[0]);       
    float theta_x_us=asin(us_dx/pow(us_dx*us_dx+us_dz*us_dz,.5));
    float theta_x_ds=asin(ds_dx/pow(ds_dx*ds_dx+ds_dz*ds_dz,.5));
    reco_pz = (fabs(fB_field_tesla) * fL_eff * fmm_to_m * fGeV_to_MeV ) / float (3.3*(theta_x_ds - theta_x_us)*cos(atan(BestTrackStats[0])));
    //std::cout<<"momentum set"<<std::endl;
    //Scale depending on current setting
    if(current>50 && current<70){
      reco_pz=(reco_pz+147.3)/1.358;
    }
    if(current>90){
    reco_pz=(reco_pz+214.7)/1.468;
    }
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
      float midplane_slope=tan((3.141592654/180)*8.0)*fMP_M;
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
      reco_pz =(fabs(fB_field_tesla) * fL_eff * fmm_to_m * fGeV_to_MeV ) / float(3.3 * (theta_x_ds-theta_x_us)*cos(atan(BestTrackStats[0]))); 
   //Calibrate depending on current
    if(current>50 && current<70){
      reco_pz=(reco_pz-79.91)/.7646;
    }
    if(current>90){
      reco_pz=(reco_pz-121.35)/77.12;
    }
    
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
    float theta_x_us_reco=-(fabs(fB_field_tesla) * fL_eff * fmm_to_m * fGeV_to_MeV )/(3.3*reco_pz*cos(atan(BestTrackStats[0]))) + theta_x_ds;
    //We still need X2,Z2, both of which are functions of the X wire of the hit.  We solve for the x wire in WC2 and then find X2,Z2.
    float missedXwire=(fX_cntr[1]-x[0]-(fZ_cntr[1]-z[0])*tan(theta_x_us_reco))/(cos((3.141592654/180)*13.0)-sin((3.141592654/180)*13.0)*tan(theta_x_us_reco));
    x[1]=fX_cntr[1]+cos((3.141592654/180)*13.0)*missedXwire;
    z[1]=fZ_cntr[1]+sin((3.141592654/180)*13.0)*missedXwire;
    //Then use the regression to get y from z.
    y[1]=BestTrackStats[0]*z[1]+BestTrackStats[1];
    float missedYwire=y[1]-fY_cntr[1];
    std::cout<<"Missed X wire: "<<missedXwire<<" Missed Y Wire : "<<missedYwire<<std::endl;
    std::cout<<"extrapolated x: "<<x[1]<<" extrapolated y: "<<y[1]<<std::endl;
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
    float theta_x_ds_reco=(fabs(fB_field_tesla) * fL_eff * fmm_to_m * fGeV_to_MeV )/(3.3*reco_pz*cos(atan(BestTrackStats[0])))+theta_x_us; 
    float missedXwire=(x[3]-fX_cntr[2]-(z[3]-fZ_cntr[2])*tan(theta_x_ds_reco))/(cos((3.141592654/180)*3.0)-(sin((3.141592654/180)*3.0))*tan(theta_x_ds_reco));
    x[2]=fX_cntr[2]+cos((3.141592654/180)*3.0)*missedXwire;
    z[2]=fZ_cntr[2]+sin((3.141592654/180)*3.0)*missedXwire;
    y[2]=BestTrackStats[0]*z[2]+BestTrackStats[1];
    float missedYwire=y[2]-fY_cntr[2];
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
		  std::cout<<"Regression value : "<<track_stats[2]<<std::endl;
		  if(track_stats[2]<fabs(bestResSq)){
		    best_track=track;
		    bestRegressionStats=track_stats;
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
  std::cout<<bestResSq<<std::endl;
//Now we should have the straightest track in Y, which will be the track that goes to the event.  Now we get the momentum and projections onto the TPC  
  calculateTheMomentum(best_track,x,y,z,reco_pz,bestRegressionStats);
  float fourx[4]={x[0],x[1],x[2],x[3]};
  float foury[4]={y[0],y[1],y[2],y[3]};
  float fourz[4]={z[0],z[1],z[2],z[3]};
  float fourmom=reco_pz;
  float fourwires[8]={0,0,0,0,0,0,0,0};
  for(int i=0; i<8; ++i){
    fourwires[i]=best_track.hits[i].wire;
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
  
  
  
  
  NHits=3;
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
  std::cout<<"WCMissed"<<WCMissed<<std::endl;
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
		    std::cout<<"Track Res Sq :"<<track_stats[2]<<std::endl;
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
  std::cout<<"Best track residual"<<bestRegressionStats[2]<<std::endl;
  float reco_pz_three=0;  
//Now we should have the straightest track in Y, which will be the track that goes to the event.  Now we get the momentum and projections onto the TPC  
  calculateTheThreePointMomentum(best_track,x,y,z,reco_pz_three,bestRegressionStats,WCMissed);
  float twomom=reco_pz_three;
 // std::cout<<"Three Momentum calculated"<<std::endl;
  //reco_pz_list.push_back(reco_pz_three);
  //Recodiff[0]->Fill(reco_four,reco_pz_three);//Now that the momentum is scaled back to what we would get with four points, try to find where the missed point would be to reconstruct this scaled momentum
  extrapolateTheMissedPoint(best_track,x,y,z,reco_pz_three,bestRegressionStats, missed_wire_hits,WCMissed);
  //std::cout<<"Missed Point Extrapolated"<<std::endl;
  std::cout<<"Missed X wire: "<<missed_wire_hits[0]<<" Missed Y Wire : "<<missed_wire_hits[1]<<std::endl;
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
  std::cout<<"WCMissed"<<WCMissed<<std::endl;
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
		     std::cout<<"Track Res Sq :"<<track_stats[2]<<std::endl;
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
  std::cout<<"Missed X wire: "<<missed_wire_hits[0]<<" Missed Y Wire : "<<missed_wire_hits[1]<<std::endl;
//Put the extrapolated hit onto the track
  if(WCMissed==2){
    best_track.hits[2].wire=missed_wire_hits[0];
    best_track.hits[3].wire=missed_wire_hits[1];
  }
  if(WCMissed==3){
    best_track.hits[4].wire=missed_wire_hits[0];
    best_track.hits[5].wire=missed_wire_hits[1];
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
    
  }
  for(int i=0; i<2; ++i){
    Recodiff[i+40]->Fill(fourface[i],twoface[i]);
    Recodiff[i+42]->Fill(fourface[i],threeface[i]);
    Recodiff[i+44]->Fill(fourangles[i],twoangles[i]);
    Recodiff[i+46]->Fill(fourangles[i],threeangles[i]);
  }
  Recodiff[56]->Fill(fourmom,twomom);
  Recodiff[57]->Fill(fourmom,threemom);
}
