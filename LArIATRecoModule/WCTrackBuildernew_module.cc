////////////////////////////////////////////////////////////////////////
// Class:       WCTrackBuildernew
// Module Type: producer
// File:        WCTrackBuildernew_module.cc
//
// Generated at Fri Oct 16 14:58:18 2015 by Greg Pulliam using artmod
// from cetpkgsupport v1_08_06.
////////////////////////////////////////////////////////////////////////

#ifndef WCTRACKBUILDERNEW_H
#define WCTRACKBUILDERNEW_H


#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include <iostream>

#include <vector>
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "art/Framework/Services/Optional/TFileService.h"


//ROOT Things
#include <TH1F.h>
#include <TH2F.h>
#include <TTree.h>

//LArIAT Things
#include "RawDataUtilities/TriggerDigitUtility.h"
#include "LArIATRecoAlg/WCTrackBuilderAlg.h"
#include "LArIATRecoAlg/WCHitFinderAlg.h"
#include "LArIATDataProducts/WCTrack.h"
#include "Utilities/DatabaseUtilityT1034.h"

#include <memory>
#include <utility>
#include <string>
#include <fstream>
namespace wct{
class WCTrackBuildernew;

class WCTrackBuildernew : public art::EDProducer {
public:
  explicit WCTrackBuildernew(fhicl::ParameterSet const & p);
  // The destructor generated by the compiler is fine for classes
  // without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  WCTrackBuildernew(WCTrackBuildernew const &) = delete;
  WCTrackBuildernew(WCTrackBuildernew &&) = delete;
  WCTrackBuildernew & operator = (WCTrackBuildernew const &) = delete;
  WCTrackBuildernew & operator = (WCTrackBuildernew &&) = delete;

  // Required functions.
  void produce(art::Event & e) override;

  // Selected optional functions.
  void beginJob() override;
  //void beginJob(fhicl::ParameterSet const & p);
  void beginRun(art::Run & r) override;
  void beginSubRun(art::SubRun & sr) override;
  void endJob() override;
//  void endRun(art::Run & r) override;
//  void endSubRun(art::SubRun & sr) override;
  void reconfigure(fhicl::ParameterSet const & p) override;
//  void respondToCloseInputFile(art::FileBlock const & fb) override;
//  void respondToCloseOutputFiles(art::FileBlock const & fb) override;
//  void respondToOpenInputFile(art::FileBlock const & fb) override;
//  void respondToOpenOutputFiles(art::FileBlock const & fb) override;
  void convertDigitsToVectors( std::vector<raw::AuxDetDigit> the_digits_1,
			       std::vector<raw::AuxDetDigit> the_digits_2,
			       std::vector<raw::AuxDetDigit> the_digits_3,
			       std::vector<raw::AuxDetDigit> the_digits_4,
			       std::vector<int> & tdc_number_vect,
			       std::vector<float> & hit_channel_vect,
			       std::vector<float> & hit_time_bin_vect );

  void createAuxDetStyleVectorsFromHitLists(WCHitList final_track,
					    std::vector<int> & WC_vect,
					    std::vector<float> & hit_wire_vect);
					    
					    
  void plotTheTrackInformation( std::vector<double> reco_pz_list,
				std::vector<double> x_face_list,
				std::vector<double> y_face_list,
				std::vector<double> theta_list,
				std::vector<double> phi_list,
				std::vector<double> y_kink_list,
				std::vector<double> x_dist_list,
				std::vector<double> y_dist_list,
				std::vector<double> z_dist_list);
  void ResetTree();
				
				
  void MakeSomePlotsFromHits(std::vector<std::vector<WCHitList> > good_hits);
  void PlotHitsBeforeClustering(std::vector<int> tdcs, std::vector<float> channels, std::vector<float> times);
private:

  // Declare member data here.

    //Offset ont he B field
    float offset;

    int subRun = -999;
    int run = -999;

    //Algorithm object for track building
    WCTrackBuilderAlg fWCTrackBuilderAlg;
    std::string       fSlicerSourceLabel;
    WCHitFinderAlg    fWCHitFinderAlg;

    //Hardware constants
    int fNumber_wire_chambers;
    int fNumber_wires_per_tdc;
    
    //int evtcounter = 0;
    //Histograms for plotting
    TH1F* fReco_Pz;
    TH1F* fY_Kink;
    TH1F* fX_Dist;
    TH1F* fY_Dist;
    TH1F* fZ_Dist;
    TH1F* fX_Face_Dist;
    TH1F* fY_Face_Dist;
    TH1F* fTheta_Dist;
    TH1F* fPhi_Dist;
    TH1F* fTrack_Type;
    std::vector<TH2F*> fRecodiff;
    TH1F* fWCDist;
    TTree* fTree; 
    
    int nevt=0;     
    //Misc
    bool fVerbose;
    bool fPickyTracks;
    bool fCheckTracks;
    int WCx1[200][100],WCx2[200][100],WCx3[200][100],WCx4[200][100],WCy1[200][100],WCy2[200][100],WCy3[200][100],WCy4[200][100]; //Arrays for the wires hit in each WC

    //Hits and Times before clustering. I want to see the noise hits.
    float WC1Xhit[200][100],WC2Xhit[200][100],WC3Xhit[200][100],WC4Xhit[200][100];
    float WC1Yhit[200][100],WC2Yhit[200][100],WC3Yhit[200][100],WC4Yhit[200][100];
    float WC1Xtime[200][100],WC2Xtime[200][100],WC3Xtime[200][100],WC4Xtime[200][100];
    float WC1Ytime[200][100],WC2Ytime[200][100],WC3Ytime[200][100],WC4Ytime[200][100];   
    float WC1XMult,WC2XMult,WC3XMult,WC4XMult, WC1YMult,WC2YMult,WC3YMult,WC4YMult;
    bool PickyTrackCheck;
};


WCTrackBuildernew::WCTrackBuildernew(fhicl::ParameterSet const & p)
 : fWCTrackBuilderAlg(p.get< fhicl::ParameterSet > ("WCTrackBuilderAlg")) // these should be initialized
 , fWCHitFinderAlg(p.get< fhicl::ParameterSet >("WCHitFinderAlg"))            // here instead of reconfigure()
{
  // Call appropriate produces<>() functions here.
      this->reconfigure(p);

    // Call appropriate produces<>() functions here.  
    produces<std::vector<ldp::WCTrack> >();
}

void WCTrackBuildernew::produce(art::Event & e)
{
    // If this is a new subrun in *real* data, then load B-field from DB.
    // (For MC, this is taken from the fhicl parameter "MCMagneticFieldTesla")
    if( e.isRealData() ){
      bool newEvent = (run != (int)e.run() && subRun != (int)e.subRun() );
      if( newEvent ) {
        fWCTrackBuilderAlg.loadXMLDatabaseTableForBField( e.run(), e.subRun() );
        subRun  = e.subRun();
        run     = e.run();
      }
    }
  
    //std::cout<<"Event start"<<std::endl;
    ResetTree();
    //Creating the WCTrack Collection
    std::unique_ptr<std::vector<ldp::WCTrack> > WCTrackCol(new std::vector<ldp::WCTrack> );  

    //Retrieving the digits from the sliced event
    art::Handle< std::vector<raw::AuxDetDigit> > AuxDetDigitHandle;
    e.getByLabel(fSlicerSourceLabel,AuxDetDigitHandle);
    
    //Loop through the auxdetdigits and collect those that are from the WCs
    std::vector<raw::AuxDetDigit> WC1Digits;
    std::vector<raw::AuxDetDigit> WC2Digits;
    std::vector<raw::AuxDetDigit> WC3Digits;
    std::vector<raw::AuxDetDigit> WC4Digits;
    for( size_t iDig = 0; iDig < AuxDetDigitHandle->size(); ++iDig ){
      if( AuxDetDigitHandle->at(iDig).AuxDetName() == "MWPC1" )
	WC1Digits.push_back(AuxDetDigitHandle->at(iDig));
      if( AuxDetDigitHandle->at(iDig).AuxDetName() == "MWPC2" )
	WC2Digits.push_back(AuxDetDigitHandle->at(iDig));
      if( AuxDetDigitHandle->at(iDig).AuxDetName() == "MWPC3" )
	WC3Digits.push_back(AuxDetDigitHandle->at(iDig));
      if( AuxDetDigitHandle->at(iDig).AuxDetName() == "MWPC4" )
	WC4Digits.push_back(AuxDetDigitHandle->at(iDig));
    }  
    std::vector<int> tdc_number_vect;
    std::vector<float> hit_channel_vect;
    std::vector<float> hit_time_bin_vect;
    convertDigitsToVectors( WC1Digits,
			    WC2Digits,
			    WC3Digits,
			    WC4Digits,
			    tdc_number_vect,
			    hit_channel_vect,
			    hit_time_bin_vect ); 
			      
   		    		    
    std::vector<double> reco_pz_list;                  //Final reco pz result for full_track_info = true, not indexed by trigger
    std::vector<double> reco_pz2M_list;
    std::vector<double> unscaled_reco_pz_list;
    std::vector<double> x_face_list;
    std::vector<double> y_face_list;
    std::vector<double> theta_list;
    std::vector<double> phi_list;
    std::vector<double> y_kink_list;
    std::vector<double> x_dist_list;
    std::vector<double> y_dist_list;
    std::vector<double> z_dist_list;
    std::vector<WCHitList> final_tracks;
    float hit_position_vect[4][3];
    float residual;  
    std::vector<std::vector<WCHitList> > good_hits; //Two vectors: WC#, axis. - Will be cleared for each trigger 
    int WCMissed; //The WC missed for the event, if there is one.
    //Initializing the good hit arrays to a default state - these clear for every trigger
    //Have 2-dimensional array of hitlists:
    //1st Dim: WC
    //2nd Dim: Axis    
    WCHitList hitList;
    std::vector<WCHitList> hitListAxis;
    for( int iAx = 0; iAx < 2; ++iAx ){ hitListAxis.push_back(hitList); }
    for( int iWC = 0; iWC < fNumber_wire_chambers; ++iWC ){ good_hits.push_back(hitListAxis); }
    
    //initialize the position array for the hits in the track put on the event
    for(int i=0; i<4; ++i){
      for(int j=0; j<3; ++j){
       hit_position_vect[i][j]=99999;
       }
     }
    
    //int good_trigger_counter = 0;
    fWCHitFinderAlg.createHits(tdc_number_vect,
    			       hit_channel_vect,
			       hit_time_bin_vect,
			       good_hits,
			       fVerbose);
			       
   		
  WC1XMult=(float)(good_hits[0][0].hits.size());
  WC2XMult=(float)(good_hits[1][0].hits.size());
  WC3XMult=(float)(good_hits[2][0].hits.size());
  WC4XMult=(float)(good_hits[3][0].hits.size());
  
  WC1YMult=(float)(good_hits[0][1].hits.size());
  WC2YMult=(float)(good_hits[1][1].hits.size());
  WC3YMult=(float)(good_hits[2][1].hits.size());
  WC4YMult=(float)(good_hits[3][1].hits.size());        
  fTrack_Type->Fill(fWCHitFinderAlg.getTrackType(good_hits));
  double ScalingFactor=fWCTrackBuilderAlg.GetScalingFactor();
//std::cout<<"Hit Finding done, going to Track Building"<<std::endl;
  fWCTrackBuilderAlg.reconstructTracks(  reco_pz_list,
                                         reco_pz2M_list,
					 x_face_list,
					 y_face_list,
					 theta_list,
					 phi_list,
					 final_tracks,
					 good_hits,
					 fPickyTracks,
					 fCheckTracks,
					 y_kink_list,
					 x_dist_list,
					 y_dist_list,
					 z_dist_list,
					 WCMissed,
					 fRecodiff,
					 fWCDist,
					 residual,
					 hit_position_vect,
                                         offset);			       
  // Convert quantities from mm to cm, since WCTrackBuilderalg uses mm, 
  // but LArSoft uses cm for everything
  for (size_t i = 0; i < x_face_list.size(); i++) x_face_list[i] *= 0.1;
  for (size_t i = 0; i < y_face_list.size(); i++) y_face_list[i] *= 0.1;
  for (size_t i = 0; i < x_dist_list.size(); i++) x_dist_list[i] *= 0.1;
  for (size_t i = 0; i < y_dist_list.size(); i++) y_dist_list[i] *= 0.1;
  for (size_t i = 0; i < z_dist_list.size(); i++) z_dist_list[i] *= 0.1;
  for (size_t i = 0; i < 4; i++)
  {
    for (size_t j = 0; j < 3; j++)
    {
        hit_position_vect[i][j] *= 0.1;
    }
  }
  //std::cout<<" reco_pz : "<<reco_pz_list.size()<<std::endl;
  //if(reco_pz2M_list.size())  std::cout<<"Checking the filling up of reco_pz2m: "<<reco_pz2M_list.size()<<" "<<reco_pz2M_list[0]<<" "<<reco_pz_list[0]<<std::endl;
     //Pick out the tracks created under this current trigger and fill WCTrack objects with info.
    //(This must be done because the track/etc. lists encompass all triggers
    int tracknumber=final_tracks.size();
    for( int iNewTrack = 0; iNewTrack<tracknumber; ++iNewTrack ){
      std::vector<int> WC_vect;
      std::vector<float> hit_wire_vect;
      
      WCHitList final_track = final_tracks[iNewTrack];
      unscaled_reco_pz_list.push_back(reco_pz_list[iNewTrack]/(1/(1+ScalingFactor)));
      
      
      //Filling as done above, but formats the WC and hit wire vectors in the WCAuxDetDigit style
      createAuxDetStyleVectorsFromHitLists(final_track,
					   WC_vect,
					   hit_wire_vect);
 					   
  if(WC1XMult==1 && WC2XMult==1 && WC3XMult==1 && WC4XMult==1 && WC1YMult==1 && WC2YMult==1 && WC3YMult==1 && WC4YMult==1) { PickyTrackCheck=true;} //We probably already know from .fcl params, but we're gonna check if this is a picky track to add to the WCTrack Object
 else{PickyTrackCheck=false;}      

      //WCTrack object creation and association with trigger created
if(reco_pz2M_list.size() > 0){      ldp::WCTrack the_track(reco_pz_list[iNewTrack],
	                     reco_pz2M_list[iNewTrack],
                             y_kink_list[iNewTrack],
			     x_dist_list[iNewTrack],
			     y_dist_list[iNewTrack],
			     z_dist_list[iNewTrack],
			     x_face_list[iNewTrack],
			     y_face_list[iNewTrack],
			     theta_list[iNewTrack],
			     phi_list[iNewTrack],
			     WC_vect,
			     hit_wire_vect,
			     hit_position_vect,
			     WCMissed,
			     residual,			     
			     WC1XMult,
			     WC2XMult,
			     WC3XMult,
			     WC4XMult,
			     WC1YMult,
			     WC2YMult,
			     WC3YMult,
			     WC4YMult,			     
			     PickyTrackCheck,
			     unscaled_reco_pz_list[iNewTrack]);
      (*WCTrackCol).push_back( the_track );}   
else{
ldp::WCTrack the_track(reco_pz_list[iNewTrack],
                             y_kink_list[iNewTrack],
			     x_dist_list[iNewTrack],
			     y_dist_list[iNewTrack],
			     z_dist_list[iNewTrack],
			     x_face_list[iNewTrack],
			     y_face_list[iNewTrack],
			     theta_list[iNewTrack],
			     phi_list[iNewTrack],
			     WC_vect,
			     hit_wire_vect,
			     hit_position_vect,
			     WCMissed,
			     residual,
			     WC1XMult,
			     WC2XMult,
			     WC3XMult,
			     WC4XMult,
			     WC1YMult,
			     WC2YMult,
			     WC3YMult,
			     WC4YMult,	
			     PickyTrackCheck,
			     unscaled_reco_pz_list[iNewTrack]);
      (*WCTrackCol).push_back( the_track );

}
    }

    //Plot the reconstructed momentum, y_kink, and delta X, Y, Z in histos
    plotTheTrackInformation(reco_pz_list,
			    x_face_list,
			    y_face_list,
			    theta_list,
			    phi_list,
			    y_kink_list,
			    x_dist_list,
			    y_dist_list,
			    z_dist_list);
    
    
    //Put objects into event (root file)
    e.put(std::move(WCTrackCol)); 
    //std::cout<<"Event Number "<<nevt<<std::endl;
    ++nevt;
   
    fTree->Fill();
   // std::cout<<"Event Done"<<std::endl;
    //hit_position_vect.clear();//clear the position vector for the next event.				
}
//==================================================================================================
void WCTrackBuildernew::ResetTree()
{
  for(int i=0; i<100; ++i)
    for(int j=0; j<200;j++)
  {
    WCx1[j][i]=-99999;
    WCx2[j][i]=-99999;
    WCx3[j][i]=-99999;
    WCx4[j][i]=-99999;
    WCy1[j][i]=-99999;
    WCy2[j][i]=-99999;
    WCy3[j][i]=-99999;
    WCy4[j][i]=-99999;
    WC1Xhit[j][i]=-99999;
    WC2Xhit[j][i]=-99999;
    WC3Xhit[j][i]=-99999;
    WC4Xhit[j][i]=-99999;
    WC1Yhit[j][i]=-99999;
    WC2Yhit[j][i]=-99999;
    WC3Yhit[j][i]=-99999;
    WC4Yhit[j][i]=-99999;
    WC1Xtime[j][i]=-99999;
    WC2Xtime[j][i]=-99999;
    WC3Xtime[j][i]=-99999;
    WC4Xtime[j][i]=-99999;
    WC1Ytime[j][i]=-99999;
    WC2Ytime[j][i]=-99999;
    WC3Ytime[j][i]=-99999;
    WC4Ytime[j][i]=-99999;
    WC1XMult=-99999;
    WC2XMult=-99999;
    WC3XMult=-99999;
    WC4XMult=-99999;
    WC1YMult=-99999;
    WC2YMult=-99999;
    WC3YMult=-99999;
    WC4YMult=-99999;  
    
    
  }
PickyTrackCheck=false;  
}
//==================================================================================================
void WCTrackBuildernew::MakeSomePlotsFromHits(std::vector<std::vector<WCHitList> > good_hits)
{  
  for(size_t iHit=0; iHit<good_hits[0][0].hits.size(); ++iHit)
  {
    WCx1[nevt][iHit]=good_hits[0][0].hits[iHit].wire;
  }
    for(size_t iHit=0; iHit<good_hits[1][0].hits.size(); ++iHit)
  {
    WCx2[nevt][iHit]=good_hits[1][0].hits[iHit].wire;
  }
    for(size_t iHit=0; iHit<good_hits[2][0].hits.size(); ++iHit)
  {
    WCx3[nevt][iHit]=good_hits[2][0].hits[iHit].wire;
  }
    for(size_t iHit=0; iHit<good_hits[3][0].hits.size(); ++iHit)
  {
    WCx4[nevt][iHit]=good_hits[3][0].hits[iHit].wire;
  }
    for(size_t iHit=0; iHit<good_hits[0][1].hits.size(); ++iHit)
  {
    WCy1[nevt][iHit]=good_hits[0][1].hits[iHit].wire;
  }
    for(size_t iHit=0; iHit<good_hits[1][1].hits.size(); ++iHit)
  {
    WCy2[nevt][iHit]=good_hits[1][1].hits[iHit].wire;
  }
    for(size_t iHit=0; iHit<good_hits[2][1].hits.size(); ++iHit)
  {
    WCy3[nevt][iHit]=good_hits[2][1].hits[iHit].wire;
  }
    for(size_t iHit=0; iHit<good_hits[3][1].hits.size(); ++iHit)
  {
    WCy4[nevt][iHit]=good_hits[3][1].hits[iHit].wire;
  }

}
//==================================================================================================
void WCTrackBuildernew::createAuxDetStyleVectorsFromHitLists(WCHitList final_track,
								   std::vector<int> & WC_vect,
								   std::vector<float> & hit_wire_vect)
  {
    for( size_t iHit = 0; iHit < final_track.hits.size() ; ++iHit ){
      WC_vect.push_back(int(iHit/2)+1);          //Look at how hits are pushed into the tracks in buildTracksFromHits (alg)

      float the_wire = (final_track.hits.at(iHit).wire*-1)+64+(128*(iHit%2));
      if (fVerbose) { std::cout << "Old WCAxis/Wire: " << iHit << "/" << final_track.hits.at(iHit).wire << ", New WC/Wire: " << int(iHit/2)+1 << "/" << the_wire << std::endl; }
      hit_wire_vect.push_back(the_wire);    
    }
  }
  //===================================================================================
  void WCTrackBuildernew::plotTheTrackInformation( std::vector<double> reco_pz_list,
							 std::vector<double> x_face_list,
							 std::vector<double> y_face_list,
							 std::vector<double> theta_list,
							 std::vector<double> phi_list,
							 std::vector<double> y_kink_list,
							 std::vector<double> x_dist_list,
							 std::vector<double> y_dist_list,
							 std::vector<double> z_dist_list)
  {
    //Loop through the tracks and fill
    for( size_t iTrack = 0; iTrack < reco_pz_list.size(); ++iTrack ){
      fReco_Pz->Fill(reco_pz_list.at(iTrack));
      fY_Kink->Fill(y_kink_list.at(iTrack));
      fX_Dist->Fill(x_dist_list.at(iTrack));
      fY_Dist->Fill(y_dist_list.at(iTrack));
      fZ_Dist->Fill(z_dist_list.at(iTrack));
      fX_Face_Dist->Fill(x_face_list.at(iTrack));
      fY_Face_Dist->Fill(y_face_list.at(iTrack));
      fTheta_Dist->Fill(theta_list.at(iTrack));
      fPhi_Dist->Fill(phi_list.at(iTrack));
    }
    
  }
//=========================================================================================
void WCTrackBuildernew::beginJob()//fhicl::ParameterSet const & p)
{



  // Implementation of optional member function here.
      // Implementation of optional member function here.
art::ServiceHandle<art::TFileService> tfs;
fWCDist= tfs->make<TH1F>("WCCond","WC Conditions",7,0,7); 
fTree=tfs->make<TTree>("WCTree","WCTree");
fTree->Branch("WC1XMult",&WC1XMult,"WC1XMult/F");
fTree->Branch("WC2XMult",&WC2XMult,"WC2XMult/F");
fTree->Branch("WC3XMult",&WC3XMult,"WC3XMult/F");
fTree->Branch("WC4XMult",&WC4XMult,"WC4XMult/F");
fTree->Branch("WC1YMult",&WC1YMult,"WC1YMult/F");
fTree->Branch("WC2YMult",&WC2YMult,"WC2YMult/F");
fTree->Branch("WC3YMult",&WC3YMult,"WC3YMult/F");
fTree->Branch("WC4YMult",&WC4YMult,"WC4YMult/F");
//Hists that should be used for diagnostics and deleted before production
if(fCheckTracks){
  for(int i=0; i<99; ++i){
    fRecodiff.push_back(tfs->make<TH2F>());
  }
fRecodiff[0]= tfs->make<TH2F>("WC1XWire4v2","WC1XWire4v2",500,-250,250,500,-250,250);
fRecodiff[1]= tfs->make<TH2F>("WC1YWire4v2","WC1YWire4v2",500,-250,250,500,-250,250);
fRecodiff[2]= tfs->make<TH2F>("WC2XWire4v2","WC2XWire4v2",500,-250,250,500,-250,250);
fRecodiff[3]= tfs->make<TH2F>("WC2YWire4v2","WC2YWire4v2",500,-250,250,500,-250,250);
fRecodiff[4]= tfs->make<TH2F>("WC3XWire4v2","WC3XWire4v2",500,-250,250,500,-250,250);
fRecodiff[5]= tfs->make<TH2F>("WC3YWire4v2","WC3YWire4v2",500,-250,250,500,-250,250);
fRecodiff[6]= tfs->make<TH2F>("WC4XWire4v2","WC4XWire4v2",500,-250,250,500,-250,250);
fRecodiff[7]= tfs->make<TH2F>("WC4YWire4v2","WC4YWire4v2",500,-250,250,500,-250,250);


fRecodiff[8]= tfs->make<TH2F>("WC1XWire4v3","WC1XWire4v3",500,-250,250,500,-250,250);
fRecodiff[9]= tfs->make<TH2F>("WC1YWire4v3","WC1YWire4v3",500,-250,250,500,-250,250);
fRecodiff[10]= tfs->make<TH2F>("WC2XWire4v3","WC2XWire4v3",500,-250,250,500,-250,250);
fRecodiff[11]= tfs->make<TH2F>("WC2YWire4v3","WC2YWire4v3",500,-250,250,500,-250,250);
fRecodiff[12]= tfs->make<TH2F>("WC3XWire4v3","WC3XWire4v3",500,-250,250,500,-250,250);
fRecodiff[13]= tfs->make<TH2F>("WC3YWire4v3","WC3YWire4v3",500,-250,250,500,-250,250);
fRecodiff[14]= tfs->make<TH2F>("WC4XWire4v3","WC4XWire4v3",500,-250,250,500,-250,250);
fRecodiff[15]= tfs->make<TH2F>("WC4YWire4v3","WC4YWire4v3",500,-250,250,500,-250,250);

fRecodiff[16]= tfs->make<TH2F>("X14v2","X14v2",4000,0,4000,4000,0,4000);
fRecodiff[17]= tfs->make<TH2F>("X24v2","X24v2",300,550,850,300,550,850);
fRecodiff[18]= tfs->make<TH2F>("X34v2","X34v2",4000,0,4000,4000,0,4000);
fRecodiff[19]= tfs->make<TH2F>("X44v2","X44v2",4000,0,4000,4000,0,4000);

fRecodiff[20]= tfs->make<TH2F>("X14v3","X14v3",4000,0,4000,4000,0,4000);
fRecodiff[21]= tfs->make<TH2F>("X24v3","X24v3",4000,0,4000,4000,0,4000);
fRecodiff[22]= tfs->make<TH2F>("X34v3","X34v3",4000,0,4000,4000,0,4000);
fRecodiff[23]= tfs->make<TH2F>("X44v3","X44v3",4000,0,4000,4000,0,4000);

fRecodiff[24]= tfs->make<TH2F>("Y14v2","Y14v2",200,-100,100,200,-100,100);
fRecodiff[25]= tfs->make<TH2F>("Y24v2","Y24v2",200,-100,100,200,-100,100);
fRecodiff[26]= tfs->make<TH2F>("Y34v2","Y34v2",200,-100,100,200,-100,100);
fRecodiff[27]= tfs->make<TH2F>("Y44v2","Y44v2",200,-100,100,200,-100,100);

fRecodiff[28]= tfs->make<TH2F>("Y14v3","Y14v3",200,-100,100,200,-100,100);
fRecodiff[29]= tfs->make<TH2F>("Y24v3","Y24v3",200,-100,100,200,-100,100);
fRecodiff[30]= tfs->make<TH2F>("Y34v3","Y34v3",200,-100,100,200,-100,100);
fRecodiff[31]= tfs->make<TH2F>("Y44v3","Y44v3",200,-100,100,200,-100,100);

fRecodiff[32]= tfs->make<TH2F>("Z14v2","Z14v2",10000,-10000,0,10000,-10000,0);
fRecodiff[33]= tfs->make<TH2F>("Z24v2","Z24v2",100,-5400,-5300,100,-5400,-5300);
fRecodiff[34]= tfs->make<TH2F>("Z34v2","Z34v2",10000,-10000,0,10000,-10000,0);
fRecodiff[35]= tfs->make<TH2F>("Z44v2","Z44v2",10000,-10000,0,10000,-10000,0);

fRecodiff[36]= tfs->make<TH2F>("Z14v3","Z14v3",10000,-10000,0,10000,-10000,0);
fRecodiff[37]= tfs->make<TH2F>("Z24v3","Z24v3",10000,-10000,0,10000,-10000,0);
fRecodiff[38]= tfs->make<TH2F>("Z34v3","Z34v3",10000,-10000,0,10000,-10000,0);
fRecodiff[39]= tfs->make<TH2F>("Z44v3","Z44v3",10000,-10000,0,10000,-10000,0);

fRecodiff[40]=tfs->make<TH2F>("TPCx4v2","TPCx4v2",800,-400,400,800,-400,400);
fRecodiff[41]=tfs->make<TH2F>("TPCy4v2","TPCy4v2",400,-200,200,400,-200,200);
fRecodiff[42]=tfs->make<TH2F>("TPCx4v3","TPCx4v3",300,100,400,300,100,400);
fRecodiff[43]=tfs->make<TH2F>("TPCy4v3","TPCy4v3",400,-200,200,400,-200,200);

fRecodiff[44]=tfs->make<TH2F>("TPCPhi4v2","TPCPhi4v2",80,-4,4,80,-4,4);
fRecodiff[45]=tfs->make<TH2F>("TPCTheta4v2","TPCTheta4v2",100,-.5,.5,100,-.5,.5);

fRecodiff[46]=tfs->make<TH2F>("TPCPhi4v3","TPCPhi4v3",80,-4,4,80,-4,4);
fRecodiff[47]=tfs->make<TH2F>("TPCTheta4v3","TPCTheta4v3",100,-.5,.5,100,-.5,.5);

fRecodiff[48]=tfs->make<TH2F>("XDist4v2","XDist4v2",500,-250,250,500,-250,250);
fRecodiff[49]=tfs->make<TH2F>("YDist4v2","YDist4v2",100,-50,50,100,-50,50);
fRecodiff[50]=tfs->make<TH2F>("ZDist4v2","ZDist4v2",100,-50,50,100,-50,50);
fRecodiff[51]=tfs->make<TH2F>("YKink4v2","YKink4v2",200,-.1,.1,200,-.1,.1);

fRecodiff[52]=tfs->make<TH2F>("XDist4v3","XDist4v3",500,-250,250,500,-250,250);
fRecodiff[53]=tfs->make<TH2F>("YDist4v3","YDist4v3",100,-50,50,100,-50,50);
fRecodiff[54]=tfs->make<TH2F>("ZDist4v3","ZDist4v3",100,-50,50,100,-50,50);
fRecodiff[55]=tfs->make<TH2F>("YKink4v3","YKink4v3",200,-.1,.1,200,-.1,.1);

fRecodiff[56]=tfs->make<TH2F>("mom4v2","relative error of S2 mom",1000,0,2000,200,-1,1);
fRecodiff[57]=tfs->make<TH2F>("mom4v3","relative error of S3 mom",1000,0,2000,200,-1,1);
fRecodiff[62]=tfs->make<TH2F>("mom4v4","relative error of S4 mom",1000,0,2000,200,-1,1);

fRecodiff[58]=tfs->make<TH2F>("Best residual all four","best residual all four",3000,0,300,3000,0,300);
fRecodiff[59]=tfs->make<TH2F>("Best residual Skip 2","best residual Skip 2",3000,0,300,3000,0,300);
fRecodiff[60]=tfs->make<TH2F>("Best residual Skip 3","best residual Skip 3",3000,0,300,3000,0,300);
fRecodiff[61]=tfs->make<TH2F>("Best residual Skip 4","best residual Skip 4",3000,0,300,3000,0,300);

fRecodiff[63]= tfs->make<TH2F>("X44v4","X44v4",4000,0,4000,4000,0,4000);
fRecodiff[64]= tfs->make<TH2F>("Y44v4","Y44v4",200,-100,100,200,-100,100);
fRecodiff[65]= tfs->make<TH2F>("Z44v4","Z44v4",10000,-10000,0,10000,-10000,0);

fRecodiff[66]=tfs->make<TH2F>("XDist4v4","XDist4v4",500,-250,250,500,-250,250);
fRecodiff[67]=tfs->make<TH2F>("YDist4v4","YDist4v4",100,-50,50,100,-50,50);
fRecodiff[68]=tfs->make<TH2F>("ZDist4v4","ZDist4v4",100,-50,50,100,-50,50);
fRecodiff[69]=tfs->make<TH2F>("YKink4v4","YKink4v4",200,-.1,.1,200,-.1,.1);

fRecodiff[70]=tfs->make<TH2F>("TPCPhi4v4","TPCPhi4v4",80,-4,4,80,-4,4);
fRecodiff[71]=tfs->make<TH2F>("TPCTheta4v4","TPCTheta4v4",100,-.5,.5,100,-.5,.5);

fRecodiff[72]= tfs->make<TH2F>("WC4XWire4v4","WC4XWire4v4",500,-250,250,500,-250,250);
fRecodiff[73]= tfs->make<TH2F>("WC4YWire4v4","WC4YWire4v4",500,-250,250,500,-250,250);

fRecodiff[74]=tfs->make<TH2F>("TPCx4v4","TPCx4v4",800,-400,400,800,-400,400);
fRecodiff[75]=tfs->make<TH2F>("TPCy4v4","TPCy4v4",400,-200,200,400,-200,200);

fRecodiff[76]=tfs->make<TH2F>("Best residual 4v2","best residual 4v2",3000,0,300,3000,0,300);
fRecodiff[77]=tfs->make<TH2F>("Best residual 4v3","best residual 4v3",3000,0,300,3000,0,300);
fRecodiff[78]=tfs->make<TH2F>("Best residual 4v4","best residual 4v4",3000,0,300,3000,0,300);

fRecodiff[79]=tfs->make<TH2F>("XwirediffS3","Difference in Xwire hit in WC3 (Actual-recreated)",600,-300,300,600,-300,300);
fRecodiff[80]=tfs->make<TH2F>("momentum","Reconstructed momentum",2000,0,2000,2000,0,2000);
fRecodiff[81]=tfs->make<TH2F>("WC1Timediff","Difference in time tick of X and Y hit in WC1",300,-150,150,300,-150,150);
fRecodiff[82]=tfs->make<TH2F>("WC2Timediff","Difference in time tick of X and Y hit in WC2",300,-150,150,300,-150,150);
fRecodiff[83]=tfs->make<TH2F>("WC3Timediff","Difference in time tick of X and Y hit in WC3",300,-150,150,300,-150,150);
fRecodiff[84]=tfs->make<TH2F>("WC4Timediff","Difference in time tick of X and Y hit in WC4",300,-150,150,300,-150,150);
fRecodiff[85]=tfs->make<TH2F>("4momvsres"," Residual versus momentum for 4 point tracks",180,0,1800,120,0,12);
fRecodiff[86]=tfs->make<TH2F>("dougsresidualtwo", "Distance WC2 Misses line through WC1 and WC4 (y,z) vs Momentum", 180,0,1800,1000,-500,500);
fRecodiff[87]=tfs->make<TH2F>("dougsresidualthree", "Distance WC3 Misses line through WC1 and WC4 (y,z) vs Momentum", 180,0,1800,1000,-500,500);
fRecodiff[88]=tfs->make<TH2F>("MomentumError","Fractional error of momentum",180,0,1800,200,0,.2);
fRecodiff[89]=tfs->make<TH2F>("mom2minus","Fractional Change in Momentum with -3mm shift in WC2X vs Orginal Momentum",1800,0,1800,500,-.05,.05);
fRecodiff[90]=tfs->make<TH2F>("mom2plus","Fractional Change in Momentum with +3mm shift in WC2X vs Orginal Momentum",1800,0,1800,500,-.05,.05);
fRecodiff[91]=tfs->make<TH2F>("mom3minus","Fractional Change in Momentum with -3mm shift in WC3X vs Orginal Momentum",1800,0,1800,500,-.05,.05);
fRecodiff[92]=tfs->make<TH2F>("mom3plus","Fractional Change in Momentum with +3mm shift in WC3X vs Orginal Momentum",1800,0,1800,500,-.05,.05);
fRecodiff[93]=tfs->make<TH2F>("Momplusplus","Fractional Change with WC 2 +3mm amd WC 3 + 3mm",1800,0,1800,1000,-.1,.1);
fRecodiff[94]=tfs->make<TH2F>("Momplusminus","Fractional Change with WC 2 +3mm amd WC 3 - 3mm",1800,0,1800,1000,-.1,.1);
fRecodiff[95]=tfs->make<TH2F>("Momminusplus","Fractional Change with WC 2 -3mm amd WC 3 + 3mm",1800,0,1800,1000,-.1,.1);
fRecodiff[96]=tfs->make<TH2F>("Momminusminus","Fractional Change with WC 2 -3mm amd WC 3 - 3mm",1800,0,1800,1000,-.1,.1);
fRecodiff[97]=tfs->make<TH2F>("XZMidplane", "XZ point halfway between line of closest approach", 1000,-5000,-4000,1000,0,1000);
fRecodiff[98]=tfs->make<TH2F>("dist","Distance of Closest Approach", 100,0,100,100,0,100);

}

//Hists that should probably stay for the production run.    
    fReco_Pz = tfs->make<TH1F>("Reco_Pz","Reconstructed momentum in XZ plane", 180, 0, 1800);
    fY_Kink = tfs->make<TH1F>("Y_Kink","Angle between US/DS tracks in Y direction (degrees)",200,-10*3.1415926/180,10*3.141592654/180);
    fX_Dist = tfs->make<TH1F>("X_Dist","X distance between US/DS tracks at midplane (mm)",1200,-60,1260);
    fY_Dist = tfs->make<TH1F>("Y_Dist","Y distance between US/DS tracks at midplane (mm)",1200,-600,600);
    fZ_Dist = tfs->make<TH1F>("Z_Dist","Z distance between US/DS tracks at midplane (mm)",1200,-60,1260);
    fX_Face_Dist = tfs->make<TH1F>("X_Face","X Location of Track's TPC Entry (mm)",1600,-200,1400);
    fY_Face_Dist = tfs->make<TH1F>("Y_Face","Y Location of Track's TPC Entry (mm)",800,-400,400);
    fTheta_Dist = tfs->make<TH1F>("Theta","Track Theta (w.r.t. TPC Z axis), (radians),",400,-.4,0.4);
    fPhi_Dist = tfs->make<TH1F>("Phi","Track Phi (w.r.t. TPC X axis), (radians)",2000,-6.28318,6.28318);                   
    fReco_Pz->GetXaxis()->SetTitle("Reconstructed momentum (MeV/c)");
    fReco_Pz->GetYaxis()->SetTitle("Tracks per 10 MeV/c");
    fY_Kink->GetXaxis()->SetTitle("Reconstructed y_kink (radians)");
    fY_Kink->GetYaxis()->SetTitle("Tracks per 0.000872 radians");
    fX_Dist->GetXaxis()->SetTitle("X distance between US and DS track ends");
    fX_Dist->GetYaxis()->SetTitle("Tracks per 1 mm");
    fY_Dist->GetXaxis()->SetTitle("Y distance between US and DS track ends");
    fY_Dist->GetYaxis()->SetTitle("Tracks per 1 mm");
    fZ_Dist->GetXaxis()->SetTitle("Z distance between US and DS track ends");
    fZ_Dist->GetYaxis()->SetTitle("Tracks per 1 mm");
    fX_Face_Dist->GetXaxis()->SetTitle("X (mm)");
    fX_Face_Dist->GetYaxis()->SetTitle("Tracks per 1 mm");
    fY_Face_Dist->GetXaxis()->SetTitle("Y (mm)");
    fY_Face_Dist->GetYaxis()->SetTitle("Tracks per 1 mm");
    fTheta_Dist->GetXaxis()->SetTitle("Theta (radians)");
    fTheta_Dist->GetYaxis()->SetTitle("Tracks per .002 radians");
    fPhi_Dist->GetXaxis()->SetTitle("Phi (radians)");
    fPhi_Dist->GetYaxis()->SetTitle("Tracks per 0.00628 radians");
    
    fTrack_Type = tfs->make<TH1F>("TrackType","WCTrack conditions: 1=missHit,2=uniqueHits,3=lonelyHit,4=socialHits",4,0,4);
    fTrack_Type->GetYaxis()->SetTitle("# Events");
    fTrack_Type->GetXaxis()->SetTitle("Track Conditions");
}
//===============================================================================================
  void WCTrackBuildernew::convertDigitsToVectors( std::vector<raw::AuxDetDigit> the_digits_1,
						     std::vector<raw::AuxDetDigit> the_digits_2,
						     std::vector<raw::AuxDetDigit> the_digits_3,
						     std::vector<raw::AuxDetDigit> the_digits_4,
						     std::vector<int> & tdc_number_vect,
						     std::vector<float> & hit_channel_vect,
						     std::vector<float> & hit_time_bin_vect )
  {  
    //std::ofstream outfile;
    //outfile.open("rawhits100A.txt", fstream::app);
    //outfile<<"Event number: "<<evtcounter<<std::endl;	
    //int hitcounter=0;
    if (fVerbose) {  std::cout << "Digits' sizes, 1:2:3:4: " << the_digits_1.size() << ":" << the_digits_2.size() << ":"  << the_digits_3.size() << ":" << the_digits_4.size() << std::endl;}
    //Loop through digits for WC1
    for( size_t iDigit = 0; iDigit < the_digits_1.size() ; ++iDigit ){
      raw::AuxDetDigit a_digit = (the_digits_1.at(iDigit));
      for( size_t iHit = 0; iHit < a_digit.NADC() ; ++iHit ){
        if (fVerbose) { std::cout << "(TDC,channel,time): (" << int(a_digit.Channel()/fNumber_wires_per_tdc)+1 << "," << a_digit.Channel() % 64 << "," << a_digit.ADC(iHit)<< ")" << std::endl; }
	hit_time_bin_vect.push_back(a_digit.ADC(iHit));
	tdc_number_vect.push_back(int(a_digit.Channel()/fNumber_wires_per_tdc)+1);
	hit_channel_vect.push_back(a_digit.Channel() % 64);
	//outfile<<"TDC: "<<int(a_digit.Channel()/fNumber_wires_per_tdc)+1<<", Wire: "<<a_digit.Channel() % 64<<", Time: "<<a_digit.ADC(iHit)<<std::endl;
	//hitcounter++;
      }
    }

    //Loop through digits for WC2
    for( size_t iDigit = 0; iDigit < the_digits_2.size() ; ++iDigit ){
      raw::AuxDetDigit a_digit = (the_digits_2.at(iDigit));
      for( size_t iHit = 0; iHit < a_digit.NADC() ; ++iHit ){
	if( a_digit.ADC(iHit) != 0 ){
	  if (fVerbose){ std::cout << "(TDC,channel,time): (" << a_digit.Channel()/fNumber_wires_per_tdc + 5 << "," << a_digit.Channel() % 64 << "," << a_digit.ADC(iHit)<< "), --> a_digit.Channel(): " << a_digit.Channel() << ", fNumber_wires...: " << fNumber_wires_per_tdc << std::endl;}
	  hit_time_bin_vect.push_back(a_digit.ADC(iHit));
	  tdc_number_vect.push_back(int(a_digit.Channel()/fNumber_wires_per_tdc)+5);
	  hit_channel_vect.push_back(a_digit.Channel() % 64);
	  //outfile<<"TDC: "<<int(a_digit.Channel()/fNumber_wires_per_tdc)+5<<", Wire: "<<a_digit.Channel() % 64<<", Time: "<<a_digit.ADC(iHit)<<std::endl;
	  //hitcounter++;
	}
      }
    }

    //Loop through digits for WC3
    for( size_t iDigit = 0; iDigit < the_digits_3.size() ; ++iDigit ){
      raw::AuxDetDigit a_digit = (the_digits_3.at(iDigit));
      for( size_t iHit = 0; iHit < a_digit.NADC() ; ++iHit ){
	if( a_digit.ADC(iHit) != 0 ){
	  if (fVerbose){std::cout << "(TDC,channel,time): (" << int(a_digit.Channel()/fNumber_wires_per_tdc)+9 << "," << a_digit.Channel() % 64 << "," << a_digit.ADC(iHit)<< ")" << std::endl;}
	  hit_time_bin_vect.push_back(a_digit.ADC(iHit));
	  tdc_number_vect.push_back(int(a_digit.Channel()/fNumber_wires_per_tdc)+9);
	  hit_channel_vect.push_back(a_digit.Channel() % 64);
	  //outfile<<"TDC: "<<int(a_digit.Channel()/fNumber_wires_per_tdc)+9<<", Wire: "<<a_digit.Channel() % 64<<", Time: "<<a_digit.ADC(iHit)<<std::endl;
	  //hitcounter++;
	}
      }
    }

    //Loop through digits for WC4
    for( size_t iDigit = 0; iDigit < the_digits_4.size() ; ++iDigit ){
      raw::AuxDetDigit a_digit = (the_digits_4.at(iDigit));
      for( size_t iHit = 0; iHit < a_digit.NADC() ; ++iHit ){
	if( a_digit.ADC(iHit) != 0 ){
	  if (fVerbose){std::cout << "(TDC,channel,time): (" << int(a_digit.Channel()/fNumber_wires_per_tdc)+13 << "," << a_digit.Channel() % 64 << "," << a_digit.ADC(iHit)<< ")" << std::endl;}
	  hit_time_bin_vect.push_back(a_digit.ADC(iHit));
	  tdc_number_vect.push_back(int(a_digit.Channel()/fNumber_wires_per_tdc)+13);
	  hit_channel_vect.push_back(a_digit.Channel() % 64);
	  //outfile<<"TDC: "<<int(a_digit.Channel()/fNumber_wires_per_tdc)+13<<", Wire: "<<a_digit.Channel() % 64<<", Time: "<<a_digit.ADC(iHit)<<std::endl;
	  //hitcounter++;
	}
      }
    }
   //outfile<<"Number of Hits for event: "<<hitcounter<<std::endl;
   //outfile.close();
  }
void WCTrackBuildernew::PlotHitsBeforeClustering(std::vector<int> tdcs, std::vector<float> channels, std::vector<float> times)
{
  int size=tdcs.size();
  for(int itdc=0; itdc<size;++itdc)
  {
    int WCax=int((tdcs.at(itdc)-1)/2); //TDCs run from 1-16, two per WC axis. Shifts to 0-7. WC1x, WC1y, WC2x, etc.
    float wire=(((tdcs.at(itdc))%2)*64-channels.at(itdc));  // Wires will iterate from -64 to 64.
    if(WCax==0)
    {
      WC1Xhit[nevt][itdc]=wire;
      WC1Xtime[nevt][itdc]=times.at(itdc);
    }    
    if(WCax==1)
    {
      WC1Yhit[nevt][itdc]=wire;
      WC1Ytime[nevt][itdc]=times.at(itdc);
    }
    if(WCax==2)
    {
      WC2Xhit[nevt][itdc]=wire;
      WC2Xtime[nevt][itdc]=times.at(itdc);
    }
    if(WCax==3)
    {
      WC2Yhit[nevt][itdc]=wire;
      WC2Ytime[nevt][itdc]=times.at(itdc);
    }
    if(WCax==4)
    {
      WC3Xhit[nevt][itdc]=wire;
      WC3Xtime[nevt][itdc]=times.at(itdc);
    }
    if(WCax==5)
    {
      WC3Yhit[nevt][itdc]=wire;
      WC3Ytime[nevt][itdc]=times.at(itdc);
    }
    if(WCax==6)
    {
      WC4Xhit[nevt][itdc]=wire;
      WC4Xtime[nevt][itdc]=times.at(itdc);
    }
    if(WCax==7)
    {
      WC4Yhit[nevt][itdc]=wire;
      WC4Ytime[nevt][itdc]=times.at(itdc);
    }                        
  }
}
void WCTrackBuildernew::beginRun(art::Run & r)
{
  // Implementation of optional member function here.
}

void WCTrackBuildernew::beginSubRun(art::SubRun & sr)
{
    // If the field override is not set, get the actual magnetic field value for the alg
    // Check if MC -- if not, load the magnetic field from a DB:
//    if( !sr.isRealData() ) fWCTrackBuilderAlg.loadXMLDatabaseTableForBField( sr.run(), sr.subRun() );}
    
    //if( fWCTrackBuilderAlg.fMCMagneticField == 0 ){ 
    //	fWCTrackBuilderAlg.loadXMLDatabaseTableForBField( sr.run(), sr.subRun() );}
}

void WCTrackBuildernew::endJob()
{
  // Implementation of optional member function here.
}

// void WCTrackBuildernew::endRun(art::Run & r)
// {
//   // Implementation of optional member function here.
// }
// 
// void WCTrackBuildernew::endSubRun(art::SubRun & sr)
// {
//   // Implementation of optional member function here.
// }

void WCTrackBuildernew::reconfigure(fhicl::ParameterSet const & p)
{
  // Implementation of optional member function here.
    fNumber_wire_chambers = p.get<int>("NWC"); //4;  
    fNumber_wires_per_tdc = p.get<int>("NWperTDC"); //64;
    fVerbose = p.get<bool>("Verbose", false);
    fSlicerSourceLabel = p.get<std::string>("SourceLabel");
    std::cout<<"Label WC: "<<fSlicerSourceLabel<<std::endl;
    fPickyTracks=p.get<bool>("PickyTracks");
    fCheckTracks=p.get<bool>("CheckTracks");
    offset = p.get<float>("BFieldOffset");


}
// 
// void WCTrackBuildernew::respondToCloseInputFile(art::FileBlock const & fb)
// {
//   // Implementation of optional member function here.
// }
// 
// void WCTrackBuildernew::respondToCloseOutputFiles(art::FileBlock const & fb)
// {
//   // Implementation of optional member function here.
// }
// 
// void WCTrackBuildernew::respondToOpenInputFile(art::FileBlock const & fb)
// {
//   // Implementation of optional member function here.
// }
// 
// void WCTrackBuilder::respondToOpenOutputFiles(art::FileBlock const & fb)
// {
//   // Implementation of optional member function here.
// }

DEFINE_ART_MODULE(WCTrackBuildernew)
}//end namespace
#endif //WCTrackBuildernew_H
