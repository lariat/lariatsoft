////////////////////////////////////////////////////////////////////////
// Class:       TimeOfFlight
// Module Type: producer
// File:        TimeOfFlight_module.cc
// This is just to change the time stamp of this file
// Generated at Fri May 29 10:13:48 2015 by Daniel Smith using artmod
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
#include <TH1F.h>
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
// ### LArIAT Things ###
#include "RawDataUtilities/TriggerDigitUtility.h"

#include <memory>

namespace lrm {
  class TimeOfFlight;
}

class lrm::TimeOfFlight : public art::EDProducer {
public:
  explicit TimeOfFlight(fhicl::ParameterSet const & p);
  // The destructor generated by the compiler is fine for classes
  // without bare pointers or other resource use.

  // Plugins should not be copied or assigned.

  TimeOfFlight(TimeOfFlight const &) = delete;
  TimeOfFlight(TimeOfFlight &&) = delete;
  TimeOfFlight & operator = (TimeOfFlight const &) = delete;
  TimeOfFlight & operator = (TimeOfFlight &&) = delete;

  // Required functions.
  void produce(art::Event & e) override;

  std::vector<short> find_hits(std::vector<short> wv);
  std::vector<short> match_hits(std::vector<short> hits1, std::vector<short> hits2);

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

private:
  
  std::string fTriggerUtility; //<---Label for the module producing the triggers

  // Declare member data here.

};


lrm::TimeOfFlight::TimeOfFlight(fhicl::ParameterSet const & p)
// :
// Initialize member data here.
{


  // Call appropriate produces<>() functions here.
}

void lrm::TimeOfFlight::produce(art::Event & e)
{
  // Gets the trigger data 
  // ### This is a crap way to pass a string.....FIX ME!!!! ####
  fTriggerUtility = "FragmentToDigit";
  rdu::TriggerDigitUtility tdu(e, fTriggerUtility);    

  std::vector<short> tof;

  // Loops over the triggers
  for(size_t trig = 0; trig < tdu.NTriggers(); ++trig) {
 
    // Gets the data for the paddles
    std::vector<const raw::AuxDetDigit*> ust_wv = tdu.TriggerUpStreamTOFDigits(trig);
    std::vector<const raw::AuxDetDigit*> dst_wv = tdu.TriggerDownStreamTOFDigits(trig);

    if(ust_wv.size() == 2) {
    // Converting the information into vectors ... might not be a good thing
      std::vector<short> ust_v0;
      std::vector<short> ust_v1;
      std::vector<short> dst_v0;
      std::vector<short> dst_v1;
      for(size_t i = 1; i < ust_wv.at(0)->NADC(); ++i) { 
      	ust_v0.insert(ust_v0.end(), ust_wv.at(0)->ADC(i));
      	ust_v1.insert(ust_v1.end(), ust_wv.at(1)->ADC(i));
      	dst_v0.insert(dst_v0.end(), dst_wv.at(0)->ADC(i));
      	dst_v1.insert(dst_v1.end(), dst_wv.at(1)->ADC(i));
      }
    
      std::vector<short> ustof1_hits = find_hits(ust_v0);
      std::vector<short> ustof2_hits = find_hits(ust_v1); // Faults here

      std::vector<short> ustof_hits = match_hits(ustof1_hits, ustof2_hits);

      std::vector<short> dstof1_hits = find_hits(ust_v0);
      std::vector<short> dstof2_hits = find_hits(ust_v1);

      std::vector<short> dstof_hits = match_hits(dstof1_hits, dstof2_hits);

      std::vector<short> tof;
   
      for(size_t dst_hit = 0; dst_hit < dstof_hits.size(); ++dst_hit) {
	for(size_t ust_hit = 0; ust_hit < ustof_hits.size(); ++ust_hit) {
	  short differ = dstof_hits[dst_hit] - ustof_hits[ust_hit];
	  if(differ > 20 and differ < 100) {
	    tof.insert(tof.end(), dstof_hits[dst_hit]-ustof_hits[ust_hit]);
	  }
	}
      }
      
      for(size_t x = 0; x < tof.size(); ++x) { 
	std::cout << tof[x] << " ";
      }

      // At this point, we have the hits for a given trigger stored in tof
    }
  }

}

std::vector<short> lrm::TimeOfFlight::find_hits(std::vector<short> wv) {
  float threshold = -40;
  std::vector<float> gradient;
  std::vector<short> hits;
  
  bool rising_edge = false;

  for(unsigned short i = 2; i < wv.size(); ++i) {
    gradient.insert(gradient.end(),float(wv.at(i-2)-wv.at(i))/2);
    
    if(gradient.back() < threshold && rising_edge == false) {

      hits.insert(hits.end(),i);

      rising_edge = true;
    }

    if(gradient.back() > abs(threshold) && rising_edge == true) {

      rising_edge = false;
    }

  }
  
  if(hits.size() == 0) { hits.insert(hits.end(), 0); }

  /*
  for(size_t h = 0; h < hits.size(); ++h) {
    std::cout << hits[h] << ", ";
    } */
  
  return hits;

}

std::vector<short> lrm::TimeOfFlight::match_hits(std::vector<short> hits1, std::vector<short> hits2) {

  short time_threshold = 10; // Nanoseconds

  const size_t len_of_hits1 = hits1.size();
  const size_t len_of_hits2 = hits2.size();

  std::vector<std::vector<short>> diff_array (len_of_hits1, std::vector<short> (len_of_hits2, 0));

  std::vector<short> matched_hits;

  // Creates a table of differences for hits1 on the x and hits 2 on the y
  for(size_t row = 0; row < len_of_hits1; row++) {
    for(size_t col = 0; col < len_of_hits2; col++) {
      diff_array[row][col] = abs(int(hits1[row]-hits2[col])); //wary of int
    }
  }

  for(size_t row = 0; row < len_of_hits1; row++) {

    // Finds the index of the lowest element in 'diff_array' for a given 'row'
    // Index is stored into 'lowest'

    size_t lowest = 0;
    for(size_t col = 0; col < len_of_hits2; col++) {
      if(diff_array[row][col] < diff_array[row][lowest]) {
	lowest = col;
      }
    }

    // Saves the hit as matched if it is in the time_threshold
    if(diff_array[row][lowest] < time_threshold && hits1[row] != 0 && hits2[lowest] != 0) {
      matched_hits.insert(matched_hits.end(),std::min(hits1[row],hits2[lowest]));
    }
  }
  
  for(size_t x = 0; x < matched_hits.size(); ++x) {
    std::cout << matched_hits[x] << ", ";
  }

  // Reversed because the python script reverses it
  std::reverse(matched_hits.begin(),matched_hits.end());
  return matched_hits;
}

void lrm::TimeOfFlight::beginJob()
{
  // Implementation of optional member function here.
  art::ServiceHandle<art::TFileService> tfs;
 
}

void lrm::TimeOfFlight::beginRun(art::Run & r)
{
  // Implementation of optional member function here.
}

void lrm::TimeOfFlight::beginSubRun(art::SubRun & sr)
{

}

void lrm::TimeOfFlight::endJob()
{
  // Implementation of optional member function here.
}

void lrm::TimeOfFlight::endRun(art::Run & r)
{
  // Implementation of optional member function here.
}

void lrm::TimeOfFlight::endSubRun(art::SubRun & sr)
{
  // Implementation of optional member function here.
}

void lrm::TimeOfFlight::reconfigure(fhicl::ParameterSet const & p)
{

   /// Here lies some pseudo-code for making objects and their associations 
   
   //produces< std::vector<recob::WireChamberTrack>>(); (DON'T USE ME EXPLICITLY....FOR ILLUSTRATION ONLY)
   //produces< art::Assns<raw::Trigger, recob::WireChamberTrack>>(); (DON'T USE ME EXPLICITLY....FOR ILLUSTRATION ONLY)
   
   // implementing the ability to pass the name of the TriggerUtililty 
   fTriggerUtility   = p.get< std::string >("TriggerUtility", "FragmentToDigit");

  // Implementation of optional member function here.
}

void lrm::TimeOfFlight::respondToCloseInputFile(art::FileBlock const & fb)
{
  // Implementation of optional member function here.
}

void lrm::TimeOfFlight::respondToCloseOutputFiles(art::FileBlock const & fb)
{
  // Implementation of optional member function here.
}

void lrm::TimeOfFlight::respondToOpenInputFile(art::FileBlock const & fb)
{
  // Implementation of optional member function here.
}

void lrm::TimeOfFlight::respondToOpenOutputFiles(art::FileBlock const & fb)
{
  // Implementation of optional member function here.
}

DEFINE_ART_MODULE(lrm::TimeOfFlight)
