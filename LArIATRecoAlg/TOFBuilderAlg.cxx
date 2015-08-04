////////////////////////////////////////////////////////////////
//                                                            //
// This is a class definition for the wire chamber track      //
// builder algorithm, used to reconstruct momentum and other  //
// geometrical properties of test-beam particles passing      //
// through LArIAT's four wire chambers                        //
//                                                            //
// Authors: Ryan Linehan, rlinehan@stanford.edu               //                           
//          Johnny Ho, johnnyho@uchicago.edu                  //
//          Jason St. John, stjohn@fnal.gov                   //
//                                                            //
////////////////////////////////////////////////////////////////

//Framework includes
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "CLHEP/Units/SystemOfUnits.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"


// LArIAT includes
#include "LArIATRecoAlg/TOFBuilderAlg.h"

#include <iostream>
#include <cmath>
#include <cstdlib>

//--------------------------------------------------------------
//Constructor
TOFBuilderAlg::TOFBuilderAlg( fhicl::ParameterSet const& pset )
{
  this->reconfigure(pset);
  art::ServiceHandle<art::TFileService> tfs;
  tof_counts  = tfs->make<TH1F>("tof_counts" , "tof_counts" ,   100, 0.,   100.);
  timestamp_histo = tfs->make<TH1F>("timestamp_histo", "timestamp_histo", 10000, 0., 9000000000.);
  ustof_histo = tfs->make<TH1F>("ustof_histo", "ustof_histo", 1000, 0., 0.);
  tof_counts->GetXaxis()->SetTitle("ToF (ns)");
  tof_counts->GetYaxis()->SetTitle("N counts");

}

//--------------------------------------------------------------  
//Destructor
TOFBuilderAlg::~TOFBuilderAlg()
{

}

//--------------------------------------------------------------
void TOFBuilderAlg::reconfigure( fhicl::ParameterSet const& pset )
{

}

//--------------------------------------------------------------
std::vector<short> TOFBuilderAlg::find_hits(std::vector<short> wv) {
  // Hit finder for an inputted waveform

  // The threshold for what qualifies for a hit
  // Negative because we are looping for negative pulses
  float threshold = -40;

  // Vector that will be storing all of the found hits
  std::vector<short> hits;
  
  bool rising_edge = false;

  // Takes the gradient of the whole waveform
  // Uses the Centered Difference Quotient
  // Starts from 2 because the gradient can't be taken from the edges
  for(unsigned short i = 2; i < wv.size(); ++i) {
    
    float gradient = float(wv[i]-wv[i-2])*0.5;
    
    // Uses the rising_edge variable to test if a hit has already been found or not
    if(gradient < threshold && rising_edge == false) {
      hits.insert(hits.end(),i);
      rising_edge = true;
    }
    
    // Uses the negative threshold to indicate that we have left a hit
    if(gradient > std::abs(threshold) && rising_edge == true) {
      rising_edge = false;
    }

  }
  
  // Gives a hit of 0 if nothing else was found to prevent seg faults
  if(hits.size() == 0) { hits.insert(hits.end(), 0); }
 
  return hits;

}

//--------------------------------------------------------------

std::vector<short> TOFBuilderAlg::match_hits(std::vector<short> hits1, std::vector<short> hits2) {
  // Matches the hits found on a single TOF paddle

  // The hits have to be within the threshold, measured in nanoseconds
  short time_threshold = 10; 

  const size_t len_of_hits1 = hits1.size();
  const size_t len_of_hits2 = hits2.size();

  std::vector<std::vector<short>> diff_array (len_of_hits1, std::vector<short> (len_of_hits2, 0));

  std::vector<short> matched_hits;

  // Creates a table of differences, with  hits1 on the x and hits2 on the y
  for(size_t row = 0; row < len_of_hits1; row++) {
    for(size_t col = 0; col < len_of_hits2; col++) {
      diff_array[row][col] = abs(int(hits1[row]-hits2[col])); 
    }
  }

  // Using the difference table, we find all of lowest time differences between 
  //  hits found and, if that is in a given threshold, return that hit's TDC time

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
  
  return matched_hits;
}


std::pair <std::vector<short>, std::vector<long> > TOFBuilderAlg::get_TOF_and_TimeStamp(std::vector<const raw::AuxDetDigit*> ust_wv, std::vector<const raw::AuxDetDigit*> dst_wv){


  std::pair<std::vector<short>, std::vector<long> > p;
  // Tests if the event has 2 PMTs, the amount needed for analysis
  if(ust_wv.size() == 2) { 
    
    // Converts the digits into vectors to pass into the functions
    //   might be the wrong approach for this
    std::vector<short> ust_v0, ust_v1, dst_v0, dst_v1;
    for(size_t i = 1; i < ust_wv[0]->NADC(); ++i) { 
      ust_v0.insert(ust_v0.end(), ust_wv[0]->ADC(i));
      ust_v1.insert(ust_v1.end(), ust_wv[1]->ADC(i));
      dst_v0.insert(dst_v0.end(), dst_wv[0]->ADC(i));
      dst_v1.insert(dst_v1.end(), dst_wv[1]->ADC(i));
    }
  
  // Calls the hit finders for each waveform and then matches the hits, using functions
    std::vector<short> ustof1_hits = find_hits(ust_v0);
    std::vector<short> ustof2_hits = find_hits(ust_v1);
    
    std::vector<short> ustof_hits = match_hits(ustof1_hits, ustof2_hits);
    
    std::vector<short> dstof1_hits = find_hits(dst_v0);
    std::vector<short> dstof2_hits = find_hits(dst_v1);
    
    std::vector<short> dstof_hits = match_hits(dstof1_hits, dstof2_hits);
    
    // Loops over each of the found matched hits
    // This compares each found ust hit with each found dst hit
    for(size_t dst_hit = 0; dst_hit < dstof_hits.size(); ++dst_hit) {
      for(size_t ust_hit = 0; ust_hit < ustof_hits.size(); ++ust_hit) {
	
	// Actual calculation for the time of flight
	short differ = dstof_hits[dst_hit] - ustof_hits[ust_hit];
	
	// Continues only if the TOF is in an expected range
	if(differ > 20 and differ < 100) {
	  
	  // Adds the calculated TOF to the vector tof
	  tof.insert(tof.end(),differ);
	  
	  
	  // Adds the timestamp for each downstream hit to the vector timeStampDst
	  // dst_wv.at(0)->TimeStamp() gives the TTT (Trigger Time Tag) since a spill
	  //    each tick of the TTT is 8 ns
	  // Then add the downstream hit to that number to get out final timetamp


	  double time = (dst_wv[0]->TimeStamp()*8)+dstof_hits[dst_hit];
	  timeStampDst.insert(timeStampDst.end(), time);
	  
	  // Fills debug histos with tof, timestamp, and hit time for the ust hits
	  tof_counts->Fill(differ);
	  timestamp_histo->Fill(time);
	  ustof_histo->Fill(ustof_hits[ust_hit]); // only ust_hits that matched with dst hits   
	  
	}
      }
    }
  }

  p = std::make_pair( tof, timeStampDst);  
  return p;
}  

std::vector<short> TOFBuilderAlg::get_tof(std::vector<const raw::AuxDetDigit*> ust_wv,std::vector<const raw::AuxDetDigit*> dst_wv)
{
  tof = (get_TOF_and_TimeStamp(ust_wv, dst_wv)).first;
  return tof;
}

std::vector<long> TOFBuilderAlg::get_timeStampDst(std::vector<const raw::AuxDetDigit*> ust_wv,std::vector<const raw::AuxDetDigit*> dst_wv)
{
  timeStampDst = (get_TOF_and_TimeStamp(ust_wv, dst_wv)).second;
  return timeStampDst;
}

void TOFBuilderAlg::clear_tof_and_timeStampDst()
{
  tof.clear();
  timeStampDst.clear();
}
