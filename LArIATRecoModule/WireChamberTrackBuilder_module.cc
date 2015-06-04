////////////////////////////////////////////////////////////////////////
// Class:       WireChamberTrackBuilder
// Module Type: producer
// File:        WireChamberTrackBuilder_module.cc
//
// Generated at Sat May 30 03:42:56 2015 by Ryan Linehan using artmod
// from cetpkgsupport v1_08_05.
////////////////////////////////////////////////////////////////////////

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
#include <TH1F.h>
#include <vector>
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "Utilities/AssociationUtil.h"

//LArIAT Things
#include "RawDataUtilities/TriggerDigitUtility.h"
#include "LArIATRecoAlg/WCTrackBuilderAlg.h"
#include "LArIATDataProducts/WCTrack.h"

#include <memory>
#include <utility>

class WireChamberTrackBuilder;

class WireChamberTrackBuilder : public art::EDProducer {
public:
  explicit WireChamberTrackBuilder(fhicl::ParameterSet const & pset);
  // The destructor generated by the compiler is fine for classes
  // without bare pointers or other resource use.

   // Plugins should not be copied or assigned.
  WireChamberTrackBuilder(WireChamberTrackBuilder const &) = delete;
  WireChamberTrackBuilder(WireChamberTrackBuilder &&) = delete;
  WireChamberTrackBuilder & operator = (WireChamberTrackBuilder const &) = delete;
  WireChamberTrackBuilder & operator = (WireChamberTrackBuilder &&) = delete;

  // Required functions.
  void produce(art::Event & e) override;

  // Selected optional functions.
  void beginJob() override;
  void beginRun(art::Run & r) override;
  void beginSubRun(art::SubRun & sr) override;
  void endJob() override;
  void endRun(art::Run & r) override;
  void endSubRun(art::SubRun & sr) override;
  void reconfigure(fhicl::ParameterSet const & p) override;
  void respondToCloseInputFile(art::FileBlock const & fb) override;
  void respondToCloseOutputFiles(art::FileBlock const & fb) override;
  void respondToOpenInputFile(art::FileBlock const & fb) override;
  void respondToOpenOutputFiles(art::FileBlock const & fb) override;
  void plotTheTrackInformation( std::vector<double> reco_pz_list,
				std::vector<double> y_kink_list,
				std::vector<double> x_dist_list,
				std::vector<double> y_dist_list,
				std::vector<double> z_dist_list,
				std::vector<double> x_face_list,
				std::vector<double> y_face_list,
				std::vector<double> theta_list,
				std::vector<double> phi_list );

  void convertDigitsToVectors( art::PtrVector<raw::AuxDetDigit> the_digits_1,
			       art::PtrVector<raw::AuxDetDigit> the_digits_2,
			       art::PtrVector<raw::AuxDetDigit> the_digits_3,
			       art::PtrVector<raw::AuxDetDigit> the_digits_4,
			       std::vector<int> & tdc_number_vect,
			       std::vector<float> & hit_channel_vect,
			       std::vector<float> & hit_time_bin_vect );

  void createVectorsFromHitLists(WCHitList final_track,
				 std::vector<int> & WC_axis_vect,
				 std::vector<float> & hit_wire_vect,
				 std::vector<float> & hit_time_vect);
    

private:
  // Declare member data here.

  std::string fTriggerUtility;

  //Algorithm object for track building
  WCTrackBuilderAlg fWCTrackBuilderAlg;

  //Hardware constants
  int fNumber_wire_chambers;
  int fNumber_wires_per_tdc;
 
  //Misc
  bool fVerbose;

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


};


WireChamberTrackBuilder::WireChamberTrackBuilder(fhicl::ParameterSet const & pset)
 :fWCTrackBuilderAlg(pset.get< fhicl::ParameterSet > ("WCTrackBuilderAlg"))
// :
// Initialize member data here.
{  
  fNumber_wire_chambers = 4;  
  fNumber_wires_per_tdc = 64;
  fVerbose = false;
  // Call appropriate produces<>() functions here.
  
  produces<std::vector<ldp::WCTrack> >();
  produces<art::Assns<raw::Trigger, ldp::WCTrack> >();
  

}

void WireChamberTrackBuilder::produce(art::Event & e)
{
  // Implementation of required member function here.

   //Creating an association between the WireChamberTrack collection and the trigger
  std::unique_ptr<art::Assns<raw::Trigger, ldp::WCTrack> > TriggerWCTrackAssn(new art::Assns<raw::Trigger, ldp::WCTrack>);
  
  //Creating the WCTrack Collection
  std::unique_ptr<std::vector<ldp::WCTrack> > WCTrackCol(new std::vector<ldp::WCTrack> );  
  
  //Creating the trigger Collection
  //  std::unique_ptr<std::vector<raw::Trigger> > TriggerCol(new std::vector<raw::Trigger> );

  // ###########################################
  // ### Grab the trigger data utility (tdu) ###
  // ###########################################

  //Bad way to do this...
  fTriggerUtility = "FragmentToDigit";
  rdu::TriggerDigitUtility tdu(e, fTriggerUtility);
  
  //Track information variables
  int track_count = 0;
  std::vector<std::vector<double> > reco_pz_array;         //Final reco pz result for full_track_info = false, indexed by trigger
  std::vector<double> reco_pz_list;                     //Final reco pz result for full_track_info = true, not indexed by trigger
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
    
    //Getting the trigger object
    art::Ptr<raw::Trigger> theTrigger = (EventTriggersPtr.at(iTrig));
    //    (*TriggerCol).push_back(theTrigger);

    //Creating the Track Vector object for this trigger
    //    art::PtrVector<ldp::WCTrack> WCTrackColTrigger;
    
    //Getting the wire chamber information
    art::PtrVector<raw::AuxDetDigit> WireChamber1Digits = tdu.TriggerMWPC1DigitsPtr(iTrig);
    art::PtrVector<raw::AuxDetDigit> WireChamber2Digits = tdu.TriggerMWPC2DigitsPtr(iTrig);
    art::PtrVector<raw::AuxDetDigit> WireChamber3Digits = tdu.TriggerMWPC3DigitsPtr(iTrig);
    art::PtrVector<raw::AuxDetDigit> WireChamber4Digits = tdu.TriggerMWPC4DigitsPtr(iTrig);
    
    //Debug printing
    if( fVerbose ){ std::cout << std::endl; std::cout << "OOOOOOOOOOOOOOOOOOO TRIGGER " << iTrig << " READOUT OOOOOOOOOOOOOOOOOOOOOOO" << std::endl;
    }
    std::cout << std::endl; std::cout << "OOOOOOOOOOOOOOOOOOO TRIGGER " << iTrig << " READOUT OOOOOOOOOOOOOOOOOOOOOOO" << std::endl;
    

    //Convert the current form of the WC info to the desired form
    std::vector<int> tdc_number_vect;
    std::vector<float> hit_channel_vect;
    std::vector<float> hit_time_bin_vect;
    /*
    convertDigitsToVectors( WireChamber1Digits,
			    WireChamber2Digits,
			    WireChamber3Digits,
			    WireChamber4Digits,
			    tdc_number_vect,
			    hit_channel_vect,
			    hit_time_bin_vect );
    */



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


    //Do the track reconstruction
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

    std::cout << "Number of tracks created: " << track_count-track_count_pre << std::endl;

    //Pick out the tracks created under this current trigger and fill WCTrack objects with info
    for( int iNewTrack = 0; iNewTrack < track_count-track_count_pre; ++iNewTrack ){
      std::vector<int> WC_axis_vect;
      std::vector<float> hit_wire_vect;
      std::vector<float> hit_time_vect;
      createVectorsFromHitLists(final_tracks.at(final_tracks.size()-1-iNewTrack),
				WC_axis_vect,
				hit_wire_vect,
				hit_time_vect);
      ldp::WCTrack the_track(reco_pz_list.at(reco_pz_list.size()-1-iNewTrack),
			     y_kink_list.at(y_kink_list.size()-1-iNewTrack),
			     x_dist_list.at(x_dist_list.size()-1-iNewTrack),
			     y_dist_list.at(y_dist_list.size()-1-iNewTrack),
			     z_dist_list.at(z_dist_list.size()-1-iNewTrack),
			     x_face_list.at(x_face_list.size()-1-iNewTrack),
			     y_face_list.at(y_face_list.size()-1-iNewTrack),
			     theta_list.at(theta_list.size()-1-iNewTrack),
			     phi_list.at(phi_list.size()-1-iNewTrack),
			     WC_axis_vect,
			     hit_wire_vect,
			     hit_time_vect);
      (*WCTrackCol).push_back( the_track );
      //      WCTrackColTrigger.push_back( the_track );
      //      util::CreateAssn(*this, e, *(TriggerCol.get()), WCTrackColTrigger, *(TriggerWCTrackAssn.get()));
      util::CreateAssn(*this, e, *(WCTrackCol.get()), theTrigger, *(TriggerWCTrackAssn.get()));
    }
  
  }
  
  //Plot the reconstructed momentum, y_kink, and delta X, Y, Z
  plotTheTrackInformation(reco_pz_list,
			  y_kink_list,
			  x_dist_list,
			  y_dist_list,
			  z_dist_list,
			  x_face_list,
			  y_face_list,
			  theta_list,
			  phi_list );
  
  //Put objects and associations into event (root file)
  e.put(std::move(TriggerWCTrackAssn));
  e.put(std::move(WCTrackCol));  



}

void WireChamberTrackBuilder::createVectorsFromHitLists(WCHitList final_track,
							std::vector<int> & WC_axis_vect,
							std::vector<float> & hit_wire_vect,
							std::vector<float> & hit_time_vect)
{
  for( size_t iHit = 0; iHit < final_track.hits.size() ; ++iHit ){
    WC_axis_vect.push_back(iHit);                                      //Look at how hits are pushed into the tracks in buildTracksFromHits (alg)
    hit_wire_vect.push_back(final_track.hits.at(iHit).wire);
    hit_time_vect.push_back(final_track.hits.at(iHit).time);
  }
}


void WireChamberTrackBuilder::beginJob()
{
  // Implementation of optional member function here.
  art::ServiceHandle<art::TFileService> tfs;
  fReco_Pz = tfs->make<TH1F>("Reco_Pz","Reconstructed momentum in XZ plane", 180, 0, 1800);
  fY_Kink = tfs->make<TH1F>("Y_Kink","Angle between US/DS tracks in Y direction (degrees)",100,-5,5);
  fX_Dist = tfs->make<TH1F>("X_Dist","X distance between US/DS tracks at midplane (mm)",120,-60,60);
  fY_Dist = tfs->make<TH1F>("Y_Dist","Y distance between US/DS tracks at midplane (mm)",120,-60,60);
  fZ_Dist = tfs->make<TH1F>("Z_Dist","Z distance between US/DS tracks at midplane (mm)",120,-60,60);
  fX_Face_Dist = tfs->make<TH1F>("X_Face_Dist","X Location of Track's TPC Entry (mm)",200,-400,400);
  fY_Face_Dist = tfs->make<TH1F>("Y_Face_Dist","Y Location of Track's TPC Entry (mm)",200,-400,400);
  fTheta_Dist = tfs->make<TH1F>("Theta_Dist","Track Theta (w.r.t. TPC Z axis), (radians),",100,0,0.2);
  fPhi_Dist = tfs->make<TH1F>("Phi_Dist","Track Phi (w.r.t. TPC X axis), (radians)",100,0,6.28318);

  fReco_Pz->GetXaxis()->SetTitle("Reconstructed momentum (MeV/c)");
  fReco_Pz->GetYaxis()->SetTitle("Tracks per 10 MeV");

}

void WireChamberTrackBuilder::beginRun(art::Run & r)
{
  // Implementation of optional member function here.
}

void WireChamberTrackBuilder::beginSubRun(art::SubRun & sr)
{
  // Implementation of optional member function here.
}

void WireChamberTrackBuilder::endJob()
{
  // Implementation of optional member function here.
}

void WireChamberTrackBuilder::endRun(art::Run & r)
{
  // Implementation of optional member function here.
}

void WireChamberTrackBuilder::endSubRun(art::SubRun & sr)
{
  // Implementation of optional member function here.
}

void WireChamberTrackBuilder::reconfigure(fhicl::ParameterSet const & p)
{
  // Implementation of optional member function here.
  fTriggerUtility = p.get< std::string >("TriggerUtility","FragmentToDigit");
}

void WireChamberTrackBuilder::respondToCloseInputFile(art::FileBlock const & fb)
{
  // Implementation of optional member function here.
}

void WireChamberTrackBuilder::respondToCloseOutputFiles(art::FileBlock const & fb)
{
  // Implementation of optional member function here.
}

void WireChamberTrackBuilder::respondToOpenInputFile(art::FileBlock const & fb)
{
  // Implementation of optional member function here.
}

void WireChamberTrackBuilder::respondToOpenOutputFiles(art::FileBlock const & fb)
{
  // Implementation of optional member function here.
}

void WireChamberTrackBuilder::plotTheTrackInformation( std::vector<double> reco_pz_list,
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

//Temporary - meant to deal with inconvenient data format of digits
void WireChamberTrackBuilder::convertDigitsToVectors( art::PtrVector<raw::AuxDetDigit> the_digits_1,
						      art::PtrVector<raw::AuxDetDigit> the_digits_2,
						      art::PtrVector<raw::AuxDetDigit> the_digits_3,
						      art::PtrVector<raw::AuxDetDigit> the_digits_4,
						      std::vector<int> & tdc_number_vect,
						      std::vector<float> & hit_channel_vect,
						      std::vector<float> & hit_time_bin_vect )
{
  
  std::cout << "Digits' sizes, 1:2:3:4: " << the_digits_1.size() << ":" << the_digits_2.size() << ":"  << the_digits_3.size() << ":" << the_digits_4.size() << std::endl;

  //Loop through digits for WC1
  for( size_t iDigit = 0; iDigit < the_digits_1.size() ; ++iDigit ){
    raw::AuxDetDigit a_digit = *(the_digits_1.at(iDigit));
    for( size_t iHit = 0; iHit < a_digit.NADC() ; ++iHit ){
      if( a_digit.ADC(iHit) != 0 ){
	std::cout << "(TDC,channel,time): (" << int(a_digit.Channel()/fNumber_wires_per_tdc)+1 << "," << a_digit.Channel() % 64 << "," << a_digit.ADC(iHit)<< ")" << std::endl;
	hit_time_bin_vect.push_back(a_digit.ADC(iHit));
	tdc_number_vect.push_back(int(a_digit.Channel()/fNumber_wires_per_tdc)+1);
	hit_channel_vect.push_back(a_digit.Channel() % 64);
      }
    }
  }

  //Loop through digits for WC2
  for( size_t iDigit = 0; iDigit < the_digits_2.size() ; ++iDigit ){
    raw::AuxDetDigit a_digit = *(the_digits_2.at(iDigit));
    for( size_t iHit = 0; iHit < a_digit.NADC() ; ++iHit ){
      if( a_digit.ADC(iHit) != 0 ){
	std::cout << "(TDC,channel,time): (" << a_digit.Channel()/fNumber_wires_per_tdc + 5 << "," << a_digit.Channel() % 64 << "," << a_digit.ADC(iHit)<< "), --> a_digit.Channel(): " << a_digit.Channel() << ", fNumber_wires...: " << fNumber_wires_per_tdc << std::endl;
	hit_time_bin_vect.push_back(a_digit.ADC(iHit));
	tdc_number_vect.push_back(int(a_digit.Channel()/fNumber_wires_per_tdc)+5);
	hit_channel_vect.push_back(a_digit.Channel() % 64);
      }
    }
  }

  //Loop through digits
  for( size_t iDigit = 0; iDigit < the_digits_3.size() ; ++iDigit ){
    raw::AuxDetDigit a_digit = *(the_digits_3.at(iDigit));
    for( size_t iHit = 0; iHit < a_digit.NADC() ; ++iHit ){
      if( a_digit.ADC(iHit) != 0 ){
	std::cout << "(TDC,channel,time): (" << int(a_digit.Channel()/fNumber_wires_per_tdc)+9 << "," << a_digit.Channel() % 64 << "," << a_digit.ADC(iHit)<< ")" << std::endl;
	hit_time_bin_vect.push_back(a_digit.ADC(iHit));
	tdc_number_vect.push_back(int(a_digit.Channel()/fNumber_wires_per_tdc)+9);
	hit_channel_vect.push_back(a_digit.Channel() % 64);
      }
    }
  }

  //Loop through digits
  for( size_t iDigit = 0; iDigit < the_digits_4.size() ; ++iDigit ){
    raw::AuxDetDigit a_digit = *(the_digits_4.at(iDigit));
    for( size_t iHit = 0; iHit < a_digit.NADC() ; ++iHit ){
      if( a_digit.ADC(iHit) != 0 ){
	std::cout << "(TDC,channel,time): (" << int(a_digit.Channel()/fNumber_wires_per_tdc)+13 << "," << a_digit.Channel() % 64 << "," << a_digit.ADC(iHit)<< ")" << std::endl;
	hit_time_bin_vect.push_back(a_digit.ADC(iHit));
	tdc_number_vect.push_back(int(a_digit.Channel()/fNumber_wires_per_tdc)+13);
	hit_channel_vect.push_back(a_digit.Channel() % 64);
      }
    }
  }


}
	










DEFINE_ART_MODULE(WireChamberTrackBuilder)
