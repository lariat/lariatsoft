////////////////////////////////////////////////////////////////////////
// Class:       WCTrackBuilder
// Module Type: producer
// File:        WCTrackBulderSlicing_module.cc
//
// Generated at Tue Jun 23 23:05:11 2015 by Matthew Smylie using artmod
// from cetpkgsupport v1_08_06.
//
// Original code by Ryan Linehan.
//
////////////////////////////////////////////////////////////////////////

#ifndef WCTRACKBUILDERSLICING_H
#define WCTRACKBUILDERSLICING_H

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include <iostream>
#include <fstream>
#include <vector>
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "Utilities/AssociationUtil.h"
#include "art/Framework/Services/Optional/TFileService.h"

//ROOT Things
#include <TH1F.h>

//LArIAT Things
#include "RawDataUtilities/TriggerDigitUtility.h"
#include "LArIATRecoAlg/WCTrackBuilderAlg.h"
#include "LArIATDataProducts/WCTrack.h"
#include "Utilities/DatabaseUtilityT1034.h"

#include <memory>
#include <utility>
#include <string>

namespace wct {

  class WCTrackBuilderSlicing;

  class WCTrackBuilderSlicing : public art::EDProducer {
  public:
    explicit WCTrackBuilderSlicing(fhicl::ParameterSet const & p);
    // The destructor generated by the compiler is fine for classes
    // without bare pointers or other resource use.

    // Plugins should not be copied or assigned.
    WCTrackBuilderSlicing(WCTrackBuilderSlicing const &) = delete;
    WCTrackBuilderSlicing(WCTrackBuilderSlicing &&) = delete;
    WCTrackBuilderSlicing & operator = (WCTrackBuilderSlicing const &) = delete;
    WCTrackBuilderSlicing & operator = (WCTrackBuilderSlicing &&) = delete;

    void beginRun(art::Run & r) override;
    void beginSubRun(art::SubRun & sr) override;

    // Required functions.
    void produce(art::Event & e) override;

    void reconfigure(fhicl::ParameterSet const & p) override;
    void beginJob();
    void convertDigitsToVectors( std::vector<raw::AuxDetDigit> the_digits_1,
			       std::vector<raw::AuxDetDigit> the_digits_2,
			       std::vector<raw::AuxDetDigit> the_digits_3,
			       std::vector<raw::AuxDetDigit> the_digits_4,
			       std::vector<int> & tdc_number_vect,
			       std::vector<float> & hit_channel_vect,
			       std::vector<float> & hit_time_bin_vect );

    void createAuxDetStyleVectorsFromHitLists(WCHitList final_track,
					      std::vector<int> & WC_vect,
					      std::vector<float> & hit_wire_vect,
					      std::vector<float> & hit_time_vect);
  void plotTheTrackInformation( std::vector<double> reco_pz_list,
				std::vector<double> y_kink_list,
				std::vector<double> x_dist_list,
				std::vector<double> y_dist_list,
				std::vector<double> z_dist_list,
				std::vector<double> x_face_list,
				std::vector<double> y_face_list,
				std::vector<double> theta_list,
				std::vector<double> phi_list );


  private:
    // Declare member data here.

    //Algorithm object for track building
    WCTrackBuilderAlg fWCTrackBuilderAlg;
    std::string       fSlicerSourceLabel;

    //Hardware constants
    int fNumber_wire_chambers;
    int fNumber_wires_per_tdc;

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

 
    //Misc
    bool fVerbose;
  };

  //===============================================================================================
  WCTrackBuilderSlicing::WCTrackBuilderSlicing(fhicl::ParameterSet const & p) :fWCTrackBuilderAlg(p.get< fhicl::ParameterSet > ("WCTrackBuilderAlg"))
  // :
 
  // Initialize member data here.
  {      
    this->reconfigure(p);

    // Call appropriate produces<>() functions here.  
    produces<std::vector<ldp::WCTrack> >();
  }

  //===================================================================================
  void WCTrackBuilderSlicing::beginJob()
  {
    // Implementation of optional member function here.
    art::ServiceHandle<art::TFileService> tfs;
    fReco_Pz = tfs->make<TH1F>("Reco_Pz","Reconstructed momentum in XZ plane", 180, 0, 1800);
    fY_Kink = tfs->make<TH1F>("Y_Kink","Angle between US/DS tracks in Y direction (degrees)",100,-5*3.1415926/180,5*3.141592654/180);
    fX_Dist = tfs->make<TH1F>("X_Dist","X distance between US/DS tracks at midplane (mm)",120,-60,60);
    fY_Dist = tfs->make<TH1F>("Y_Dist","Y distance between US/DS tracks at midplane (mm)",120,-60,60);
    fZ_Dist = tfs->make<TH1F>("Z_Dist","Z distance between US/DS tracks at midplane (mm)",120,-60,60);
    fX_Face_Dist = tfs->make<TH1F>("X_Face","X Location of Track's TPC Entry (mm)",800,-200,600);
    fY_Face_Dist = tfs->make<TH1F>("Y_Face","Y Location of Track's TPC Entry (mm)",800,-400,400);
    fTheta_Dist = tfs->make<TH1F>("Theta","Track Theta (w.r.t. TPC Z axis), (radians),",100,0,0.2);
    fPhi_Dist = tfs->make<TH1F>("Phi","Track Phi (w.r.t. TPC X axis), (radians)",100,0,6.28318);

    
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
    fPhi_Dist->GetYaxis()->SetTitle("Tracks per 0.0628 radians");

    fTrack_Type = tfs->make<TH1F>("TrackType","WCTrack conditions: 1=missHit,2=uniqueHits,3=lonelyHit,4=socialHits",4,0,4);
    fTrack_Type->GetYaxis()->SetTitle("# Events");
    fTrack_Type->GetXaxis()->SetTitle("Track Conditions");
 
   
    
  }



  //===============================================================================================
  void WCTrackBuilderSlicing::produce(art::Event & e)
  { 
    
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

      
    //Track information variables
    int track_count = 0;
    std::vector<std::vector<double> > reco_pz_array;   // IGNORE!!! OBSOLETE AT THE MOMENT, BUT A BIT INVOLVED TO REMOVE, AND IT DOESN'T HURT, SO I'LL FIX LATER.
    std::vector<double> reco_pz_list;                  //Final reco pz result for full_track_info = true, not indexed by trigger
    std::vector<double> y_kink_list;
    std::vector<double> x_dist_list;
    std::vector<double> y_dist_list;
    std::vector<double> z_dist_list;
    std::vector<double> x_face_list;
    std::vector<double> y_face_list;
    std::vector<double> theta_list;
    std::vector<double> phi_list;
    std::vector<WCHitList> final_tracks;  
    std::vector<std::vector<WCHitList> > good_hits; //Two vectors: WC#, axis. - Will be cleared for each trigger

    //Initializing the good hit arrays to a default state - these clear for every trigger
    //Have 2-dimensional array of hitlists:
    //1st Dim: WC
    //2nd Dim: Axis
    WCHitList hitList;
    std::vector<WCHitList> hitListAxis;
    for( int iAx = 0; iAx < 2; ++iAx ){ hitListAxis.push_back(hitList); }
    for( int iWC = 0; iWC < fNumber_wire_chambers; ++iWC ){ good_hits.push_back(hitListAxis); }

    int good_trigger_counter = 0;
    
    //Debug printing
    //    if(fVerbose ){ std::cout << std::endl; std::cout << "OOOOOOOOOOOOOOOOOOO TRIGGER " << iTrig << " READOUT OOOOOOOOOOOOOOOOOOOOOOO" << std::endl;
    //}    
    
    //Take the AuxDetDigit information and put it into vectors that can
    //be parsed by the track reco algorithm.
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
    
    /*    
    //DO NOT ERASE - GOOD FOR SANITY CHECKS IN CASE AUXDETDIGITS AREN'T CONFIGURED CORRECTLY
    //Getting the dqm data for testing the module - temporarily read in from file
    int tdc_num = 0;
    float channel = 0;
    float time_bin = 0;
    size_t trig_counter = 0;
    ifstream myfile;
    myfile.open("run_005103_spill_0014.txt");
    while(1){
    myfile >> tdc_num >> channel >> time_bin;
    if(!myfile.good()){ break; }
    if( tdc_num == 9999 ){ 
    trig_counter++;
    continue;
    }
    if( trig_counter == iTrig ){
    tdc_number_vect.push_back(tdc_num);
    hit_channel_vect.push_back(channel);
    hit_time_bin_vect.push_back(time_bin);
    }
    }
    myfile.close();
    */
    
    //Do the track reconstruction lists are passed as reference,
      //and get filled here. By the end of this alg call, the lists
      //will have been appended with the track info for this trigger.
    int track_count_pre = track_count;
    fWCTrackBuilderAlg.reconstructTracks(tdc_number_vect,
					 hit_channel_vect,
					 hit_time_bin_vect,
					 reco_pz_array,
					 reco_pz_list,
					 y_kink_list,
					 x_dist_list,
					 y_dist_list,
					 z_dist_list,
					 x_face_list,
					 y_face_list,
					 theta_list,
					 phi_list,
					 final_tracks,
					 good_hits,
					 fVerbose,
					 good_trigger_counter,
					 0, //hardcoded in from iTrig, but never used (deprecated)
					 track_count);

    fTrack_Type->Fill(fWCTrackBuilderAlg.getTrackType());
    
    //Pick out the tracks created under this current trigger and fill WCTrack objects with info.
    //(This must be done because the track/etc. lists encompass all triggers
    for( int iNewTrack = 0; iNewTrack < track_count-track_count_pre; ++iNewTrack ){
      std::vector<int> WC_vect;
      std::vector<float> hit_wire_vect;
      std::vector<float> hit_time_vect;
      
      WCHitList final_track = final_tracks.at(final_tracks.size()-1-iNewTrack);
      
      
      //Filling as done above, but formats the WC and hit wire vectors in the WCAuxDetDigit style
      createAuxDetStyleVectorsFromHitLists(final_track,
					   WC_vect,
					   hit_wire_vect,
					   hit_time_vect);
      
      //WCTrack object creation and association with trigger created
      ldp::WCTrack the_track(reco_pz_list.at(reco_pz_list.size()-1-iNewTrack),
			     y_kink_list.at(y_kink_list.size()-1-iNewTrack),
			     x_dist_list.at(x_dist_list.size()-1-iNewTrack),
			     y_dist_list.at(y_dist_list.size()-1-iNewTrack),
			     z_dist_list.at(z_dist_list.size()-1-iNewTrack),
			     x_face_list.at(x_face_list.size()-1-iNewTrack),
			     y_face_list.at(y_face_list.size()-1-iNewTrack),
			     theta_list.at(theta_list.size()-1-iNewTrack),
			     phi_list.at(phi_list.size()-1-iNewTrack),
			     WC_vect,
			     hit_wire_vect,
			     hit_time_vect);
      (*WCTrackCol).push_back( the_track );
    }

    //Plot the reconstructed momentum, y_kink, and delta X, Y, Z in histos
    plotTheTrackInformation(reco_pz_list,
			    y_kink_list,
			    x_dist_list,
			    y_dist_list,
			    z_dist_list,
			    x_face_list,
			    y_face_list,
			    theta_list,
			    phi_list );
    
    
    //Put objects into event (root file)
    e.put(std::move(WCTrackCol));  
  }

  //=================================================================================================
  void WCTrackBuilderSlicing::reconfigure(fhicl::ParameterSet const & p)
  {
    fNumber_wire_chambers = p.get<int>("NWC"); //4;  
    fNumber_wires_per_tdc = p.get<int>("NWperTDC"); //64;
    fVerbose = p.get<bool>("Verbose", false);
    fSlicerSourceLabel = p.get<std::string>("SourceLabel");
  }

  //==================================================================================================
  void WCTrackBuilderSlicing::createAuxDetStyleVectorsFromHitLists(WCHitList final_track,
								   std::vector<int> & WC_vect,
								   std::vector<float> & hit_wire_vect,
								   std::vector<float> & hit_time_vect)
  {
    for( size_t iHit = 0; iHit < final_track.hits.size() ; ++iHit ){
      WC_vect.push_back(int(iHit/2)+1);          //Look at how hits are pushed into the tracks in buildTracksFromHits (alg)

      float the_wire = (final_track.hits.at(iHit).wire*-1)+64+(128*(iHit%2));
      if (fVerbose) { std::cout << "Old WCAxis/Wire: " << iHit << "/" << final_track.hits.at(iHit).wire << ", New WC/Wire: " << int(iHit/2)+1 << "/" << the_wire << std::endl; }
      hit_wire_vect.push_back(the_wire);
      hit_time_vect.push_back(final_track.hits.at(iHit).time);    
    }
  }

  //===================================================================================
  void WCTrackBuilderSlicing::plotTheTrackInformation( std::vector<double> reco_pz_list,
							 std::vector<double> y_kink_list,
							 std::vector<double> x_dist_list,
							 std::vector<double> y_dist_list,
							 std::vector<double> z_dist_list,
							 std::vector<double> x_face_list,
							 std::vector<double> y_face_list,
							 std::vector<double> theta_list,
							 std::vector<double> phi_list )
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

  //===================================================================================
  void WCTrackBuilderSlicing::beginRun(art::Run & r)
  {
    // Implementation of optional member function here.
  }
 
  //===================================================================================
  void WCTrackBuilderSlicing::beginSubRun(art::SubRun & sr)
  {
    // Implementation of optional member function here.
    fWCTrackBuilderAlg.loadXMLDatabaseTableForBField( sr.run(), sr.subRun() );
  }


  //===============================================================================================
  void WCTrackBuilderSlicing::convertDigitsToVectors( std::vector<raw::AuxDetDigit> the_digits_1,
						     std::vector<raw::AuxDetDigit> the_digits_2,
						     std::vector<raw::AuxDetDigit> the_digits_3,
						     std::vector<raw::AuxDetDigit> the_digits_4,
						     std::vector<int> & tdc_number_vect,
						     std::vector<float> & hit_channel_vect,
						     std::vector<float> & hit_time_bin_vect )
  {  
    if (fVerbose) {  std::cout << "Digits' sizes, 1:2:3:4: " << the_digits_1.size() << ":" << the_digits_2.size() << ":"  << the_digits_3.size() << ":" << the_digits_4.size() << std::endl;}
    
    //Loop through digits for WC1
    for( size_t iDigit = 0; iDigit < the_digits_1.size() ; ++iDigit ){
      raw::AuxDetDigit a_digit = (the_digits_1.at(iDigit));
      for( size_t iHit = 0; iHit < a_digit.NADC() ; ++iHit ){
        if (fVerbose) { std::cout << "(TDC,channel,time): (" << int(a_digit.Channel()/fNumber_wires_per_tdc)+1 << "," << a_digit.Channel() % 64 << "," << a_digit.ADC(iHit)<< ")" << std::endl; }
	hit_time_bin_vect.push_back(a_digit.ADC(iHit));
	tdc_number_vect.push_back(int(a_digit.Channel()/fNumber_wires_per_tdc)+1);
	hit_channel_vect.push_back(a_digit.Channel() % 64);
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
	}
      }
    }
  }

  DEFINE_ART_MODULE(WCTrackBuilderSlicing)

} //end namespace

#endif //WCTRACKBUILDER_H
