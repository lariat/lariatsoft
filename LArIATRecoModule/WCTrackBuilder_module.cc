////////////////////////////////////////////////////////////////////////
// Class:       WCTrackBuilder
// Module Type: producer
// File:        WCTrackBuilder_module.cc
//
// Generated at Tue Jun 23 23:05:11 2015 by Matthew Smylie using artmod
// from cetpkgsupport v1_08_06.
//
// Original code by Ryan Linehan.
//
////////////////////////////////////////////////////////////////////////

#ifndef WCTRACKBUILDER_H
#define WCTRACKBUILDER_H

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

//LArIAT Things
#include "RawDataUtilities/TriggerDigitUtility.h"
#include "LArIATRecoAlg/WCTrackBuilderAlg.h"
#include "LArIATDataProducts/WCTrack.h"

#include <memory>
#include <utility>

namespace wct {

  class WCTrackBuilder;

  class WCTrackBuilder : public art::EDProducer {
  public:
    explicit WCTrackBuilder(fhicl::ParameterSet const & p);
    // The destructor generated by the compiler is fine for classes
    // without bare pointers or other resource use.

    // Plugins should not be copied or assigned.
    WCTrackBuilder(WCTrackBuilder const &) = delete;
    WCTrackBuilder(WCTrackBuilder &&) = delete;
    WCTrackBuilder & operator = (WCTrackBuilder const &) = delete;
    WCTrackBuilder & operator = (WCTrackBuilder &&) = delete;

    // Required functions.
    void produce(art::Event & e) override;

    void reconfigure(fhicl::ParameterSet const & p) override;

    void convertDigitsToVectors( art::PtrVector<raw::AuxDetDigit> the_digits_1,
			       art::PtrVector<raw::AuxDetDigit> the_digits_2,
			       art::PtrVector<raw::AuxDetDigit> the_digits_3,
			       art::PtrVector<raw::AuxDetDigit> the_digits_4,
			       std::vector<int> & tdc_number_vect,
			       std::vector<float> & hit_channel_vect,
			       std::vector<float> & hit_time_bin_vect );

    void createAuxDetStyleVectorsFromHitLists(WCHitList final_track,
					      std::vector<int> & WC_vect,
					      std::vector<float> & hit_wire_vect,
					      std::vector<float> & hit_time_vect);


  private:
    // Declare member data here.

    //Trigger Utility
    std::string fTriggerUtility;

    //Algorithm object for track building
    WCTrackBuilderAlg fWCTrackBuilderAlg;

    //Hardware constants
    int fNumber_wire_chambers;
    int fNumber_wires_per_tdc;
 
    //Misc
    bool fVerbose;
  };

  //===============================================================================================
  WCTrackBuilder::WCTrackBuilder(fhicl::ParameterSet const & p) :fWCTrackBuilderAlg(p.get< fhicl::ParameterSet > ("WCTrackBuilderAlg"))
  // :
  // Initialize member data here.
  {      
    this->reconfigure(p);

    // Call appropriate produces<>() functions here.  
    produces<std::vector<ldp::WCTrack> >();
    produces<art::Assns<raw::Trigger, ldp::WCTrack> >();
  }

  //===============================================================================================
  void WCTrackBuilder::produce(art::Event & e)
  {    
    //Creating an association between the WireChamberTrack collection and the trigger
    std::unique_ptr<art::Assns<raw::Trigger, ldp::WCTrack> > TriggerWCTrackAssn(new art::Assns<raw::Trigger, ldp::WCTrack>);
  
    //Creating the WCTrack Collection
    std::unique_ptr<std::vector<ldp::WCTrack> > WCTrackCol(new std::vector<ldp::WCTrack> );  


    // ###########################################
    // ### Grab the trigger data utility (tdu) ###
    // ###########################################
   
    rdu::TriggerDigitUtility tdu(e, fTriggerUtility);
      
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

   
    // ##############################
    // ### Loop over the triggers ###
    // ##############################

    //Getting the trigger objects
    art::PtrVector<raw::Trigger> const& EventTriggersPtr = tdu.EventTriggersPtr();
    int good_trigger_counter = 0;
    for( size_t iTrig = 0; iTrig < tdu.NTriggers(); ++iTrig ){
    
      //Getting a new trigger
      art::Ptr<raw::Trigger> theTrigger = (EventTriggersPtr.at(iTrig));

      //Getting the wire chamber information
      art::PtrVector<raw::AuxDetDigit> WireChamber1Digits = tdu.TriggerMWPC1DigitsPtr(iTrig);
      art::PtrVector<raw::AuxDetDigit> WireChamber2Digits = tdu.TriggerMWPC2DigitsPtr(iTrig);
      art::PtrVector<raw::AuxDetDigit> WireChamber3Digits = tdu.TriggerMWPC3DigitsPtr(iTrig);
      art::PtrVector<raw::AuxDetDigit> WireChamber4Digits = tdu.TriggerMWPC4DigitsPtr(iTrig);

      //Debug printing
      if(fVerbose ){ std::cout << std::endl; std::cout << "OOOOOOOOOOOOOOOOOOO TRIGGER " << iTrig << " READOUT OOOOOOOOOOOOOOOOOOOOOOO" << std::endl;
      }    

      //Take the AuxDetDigit information and put it into vectors that can
      //be parsed by the track reco algorithm.
      std::vector<int> tdc_number_vect;
      std::vector<float> hit_channel_vect;
      std::vector<float> hit_time_bin_vect;
      convertDigitsToVectors( WireChamber1Digits,
			    WireChamber2Digits,
			    WireChamber3Digits,
			    WireChamber4Digits,
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
					   iTrig,
					   track_count);

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
	util::CreateAssn(*this, e, *(WCTrackCol.get()), theTrigger, *(TriggerWCTrackAssn.get()));
      }
    }
    
    //Put objects and associations into event (root file)
    e.put(std::move(TriggerWCTrackAssn));
    e.put(std::move(WCTrackCol));  
  }

  //=================================================================================================
  void WCTrackBuilder::reconfigure(fhicl::ParameterSet const & p)
  {
    fTriggerUtility = p.get< std::string >("TriggerUtility");
    fNumber_wire_chambers = p.get<int>("NWC"); //4;  
    fNumber_wires_per_tdc = p.get<int>("NWperTDC"); //64;
    fVerbose = p.get<bool>("Verbose");
  }

  //==================================================================================================
  void WCTrackBuilder::createAuxDetStyleVectorsFromHitLists(WCHitList final_track,
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

  //===============================================================================================
  void WCTrackBuilder::convertDigitsToVectors( art::PtrVector<raw::AuxDetDigit> the_digits_1,
						      art::PtrVector<raw::AuxDetDigit> the_digits_2,
						      art::PtrVector<raw::AuxDetDigit> the_digits_3,
						      art::PtrVector<raw::AuxDetDigit> the_digits_4,
						      std::vector<int> & tdc_number_vect,
						      std::vector<float> & hit_channel_vect,
						      std::vector<float> & hit_time_bin_vect )
  {  
    if (fVerbose) {  std::cout << "Digits' sizes, 1:2:3:4: " << the_digits_1.size() << ":" << the_digits_2.size() << ":"  << the_digits_3.size() << ":" << the_digits_4.size() << std::endl;}

    //Loop through digits for WC1
    for( size_t iDigit = 0; iDigit < the_digits_1.size() ; ++iDigit ){
      raw::AuxDetDigit a_digit = *(the_digits_1.at(iDigit));
      for( size_t iHit = 0; iHit < a_digit.NADC() ; ++iHit ){
        if (fVerbose) { std::cout << "(TDC,channel,time): (" << int(a_digit.Channel()/fNumber_wires_per_tdc)+1 << "," << a_digit.Channel() % 64 << "," << a_digit.ADC(iHit)<< ")" << std::endl; }
	hit_time_bin_vect.push_back(a_digit.ADC(iHit));
	tdc_number_vect.push_back(int(a_digit.Channel()/fNumber_wires_per_tdc)+1);
	hit_channel_vect.push_back(a_digit.Channel() % 64);
      }
    }

    //Loop through digits for WC2
    for( size_t iDigit = 0; iDigit < the_digits_2.size() ; ++iDigit ){
      raw::AuxDetDigit a_digit = *(the_digits_2.at(iDigit));
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
      raw::AuxDetDigit a_digit = *(the_digits_3.at(iDigit));
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
      raw::AuxDetDigit a_digit = *(the_digits_4.at(iDigit));
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

  DEFINE_ART_MODULE(WCTrackBuilder)

} //end namespace

#endif //WCTRACKBUILDER_H
