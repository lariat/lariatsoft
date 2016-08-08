/////////////////////////////////////////////////////////////////
//                                                             //
// This is a class definition for the time of flight algorithm // 
//                                                             //
// Authors: Daniel Smith, dansmith@bu.edu                      //                           
//                                                             //
/////////////////////////////////////////////////////////////////

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

  fderHit = tfs->make<TH1F>("fderhit","fderHit",600,-100.,100.);
  fdeltaHit = tfs->make<TH1F>("fdeltaHit","fdeltaHit",50,-5.,5.);
  fLenHit = tfs->make<TH1F>("fLenHit","fLenHit",100,0.,100.);

  fDerUSA = tfs->make<TH1F>("fDerUSA","fDerUSA",300, -15.,15.0);
  fDerUSB = tfs->make<TH1F>("fDerUSB","fDerUSB",300, -15.,15.0);
  fDerDSA = tfs->make<TH1F>("fDerDSA","fDerDSA",300, -15.,15.0);
  fDerDSB = tfs->make<TH1F>("fDerDSB","fDerDSB",300, -15.,15.0);


}

//--------------------------------------------------------------  
//Destructor
TOFBuilderAlg::~TOFBuilderAlg()
{

}

//--------------------------------------------------------------
void TOFBuilderAlg::reconfigure( fhicl::ParameterSet const& pset )
{
  // TOF Adjustment parameters
  fLinear = pset.get<float>("Linear",-10); // ns
  fMultiple = pset.get<float>("Multiple",1.0);  

  // Hit parameters
  fHitMatchThreshold = pset.get<double>("HitMatchThreshold", 2.0); // ns
  fHitThreshold = pset.get<double>("HitThreshold", -2.5); // -2.5 for run II and -9.0 for run I
  fHitWait = pset.get<double>("HitWait", 10);

}

std::pair <std::vector<float>, std::vector<long> > TOFBuilderAlg::get_TOF_and_TimeStamp(std::vector<const raw::AuxDetDigit*> ust_wv, std::vector<const raw::AuxDetDigit*> dst_wv) { 

  // Tests if the event has 2 PMTs, the amount needed for analysis
  std::pair<std::vector<float>, std::vector<long> > p;

  std::vector<float> tof_vector;
  std::vector<long> time_vector;

  if(ust_wv.size() >= 2 and dst_wv.size() >= 2) {
  
    // Converts the digits into vectors to pass into the functions

    std::vector<float> ust_v0, ust_v1, dst_v0, dst_v1;

    for(size_t iADC = 1; iADC < ust_wv[0]->NADC(); ++iADC) { 
      ust_v0.push_back(ust_wv[0]->ADC(iADC));
      ust_v1.push_back(ust_wv[1]->ADC(iADC));
      dst_v0.push_back(dst_wv[0]->ADC(iADC));
      dst_v1.push_back(dst_wv[1]->ADC(iADC));
    }

    // Histogram fun time
    for(unsigned short i = 2; i < (ust_v0.size()-2); ++i) {    
      fDerUSA->Fill(float(-ust_v0[i+2]+8.*ust_v0[i+1]-8*ust_v0[i-1]+ust_v0[i-2])/12.);
    }
    for(unsigned short i = 2; i < (ust_v1.size()-2); ++i) {    
      fDerUSB->Fill(float(-ust_v1[i+2]+8.*ust_v1[i+1]-8*ust_v1[i-1]+ust_v1[i-2])/12.);
    }
    for(unsigned short i = 2; i < (dst_v0.size()-2); ++i) {    
      fDerDSA->Fill(float(-dst_v0[i+2]+8.*dst_v0[i+1]-8*dst_v0[i-1]+dst_v0[i-2])/12.);
    }
    for(unsigned short i = 2; i < (dst_v1.size()-2); ++i) {    
      fDerDSB->Fill(float(-dst_v1[i+2]+8.*dst_v1[i+1]-8*dst_v1[i-1]+dst_v1[i-2])/12.);
    }
    
  
    // Calls the hit finders for each waveform and then matches the hits, using functions
    std::vector<float> ustof_hits0 = find_hits(ust_v0);
    std::vector<float> ustof_hits1 = find_hits(ust_v1);
    std::vector<float> dstof_hits0 = find_hits(dst_v0);
    std::vector<float> dstof_hits1 = find_hits(dst_v1);

    std::vector<float> ust_hits = match_hits(ustof_hits0, ustof_hits1);
    std::vector<float> dst_hits = match_hits(dstof_hits0, dstof_hits1);

    // Loops over each of the found matched hits
    // This compares each found ust hit with each found dst hit
    if(ust_hits.size() != 0 and dst_hits.size() != 0) {
      for(size_t idst_hit = 0; idst_hit < dst_hits.size(); ++idst_hit) {
	for(size_t iust_hit = 0; iust_hit < ust_hits.size(); ++iust_hit) {
	
	  // Actual calculation for the time of flight
	  float Corrected_TOF = fMultiple*(dst_hits[idst_hit] - ust_hits[iust_hit]) + fLinear;
	
	  // Continues only if the TOF is in an expected range
	  if(Corrected_TOF > 10 and Corrected_TOF < 100) {
	  
	    // Adds the calculated TOF to the vector tof
	    TOF.push_back(Corrected_TOF);
	    tof_vector.push_back(Corrected_TOF);
	    // Adds the timestamp for each downstream hit to the vector timeStampDst
	    // dst_wv.at(0)->TimeStamp() gives the TTT (Trigger Time Tag) since a spill
	    //    each tick of the TTT is 8 ns
	    // Then add the downstream hit to that number to get out final timetamp

	    Dst_Timestamp.push_back((dst_wv.at(0)->TimeStamp()*8) + dst_hits.at(idst_hit));
	    time_vector.push_back((dst_wv.at(0)->TimeStamp()*8) + dst_hits.at(idst_hit));	    

	  }
	}
      }
    }
  }
 
  p = std::make_pair(tof_vector, time_vector);
  return p;

}  

std::vector<float> TOFBuilderAlg::match_hits(std::vector<float> hits1, std::vector<float> hits2) {
  // Matches the hits found on a single TOF paddle

  // The hits have to be within the threshold, measured in nanoseconds
  std::vector<float> matched_hits;

  if(hits1.size() == 0 or hits2.size() == 0) {
    return matched_hits;
  }

  std::vector< std::vector<float> > diff_array (hits1.size(), std::vector<float> (hits2.size(), 0));

  // Creates a table of differences, with  hits1 on the x and hits2 on the y
  std::vector<float> column;
  for(size_t row = 0; row < hits1.size(); row++) {
    for(size_t col = 0; col < hits2.size(); col++) {
      diff_array.at(row).at(col) = (std::abs(hits1[row]-hits2[col])); 
    }
  }

  // Using the difference table, we find all of lowest time differences between 
  //  hits found and, if that is in a given threshold, return that hit's TDC time

  for(size_t row = 0; row < diff_array.size(); row++) {

    // Finds the index of the lowest element in 'diff_array' for a given 'row'
    // Index is stored into 'lowest'
    /*
    size_t lowest = 0;
    for(size_t col = 0; col < diff_array.at(row).size(); col++) {
      if(diff_array.at(row).at(col) < diff_array.at(row).at(lowest)) {
	lowest = col;
      }
    }
    */

    // Saves the hit as matched if it is in the time_threshold
    for(size_t col = 0; col < diff_array.at(row).size(); col++) {
      if(hits1.at(row) != 0 && hits2.at(col) != 0) {
	fdeltaHit->Fill(hits1.at(row)-hits2.at(col));
	if(diff_array.at(row).at(col) < 1.4 && diff_array.at(row).at(col) > -0.4) { 
	  matched_hits.push_back((hits1.at(row)+hits2.at(col))/2.);
	}
      }
    } // End of loop over col
  } // End of loop over row
  
  return matched_hits;

}

//--------------------------------------------------------------
std::vector<float> TOFBuilderAlg::find_hits(std::vector<float> wv) { 
  // Hit finder for an inputted waveform

  // Vector that will be storing all of the found hits
  std::vector<float> hits;

  bool rising_edge = false;
  int width = 0;
  int length_of_hit = 0;

  bool start_clock = false;

  for(unsigned short i = 2; i < (wv.size()-2); ++i) {
    
    float gradient = float(-wv[i+2]+8.*wv[i+1]-8*wv[i-1]+wv[i-2])/12.;

    fderHit->Fill(gradient);

    // Try to slightly update this 
    if(gradient < fHitThreshold and rising_edge == false and start_clock == false) {

      // y = gradients, x = time
      // y = mx + b
      // m = gradient - grad_before
      // gradient = m * i + b
      // b = gradient - (gradient - grad_before)*i
      // grad_before = (gradient - grad_before) * (new_time) + gradient - (gradient - grad_before)*i
      // (new_time) = (grad_before - gradient + (gradient - grad_before)*i) / (gradient - grad_before)

      float grad_before = float(-wv[i-1+2]+8.*wv[i-1+1]-8*wv[i-1-1]+wv[i-1-2])/12.;
      float m = float(gradient) - float(grad_before);
      float time = i;
      float b = gradient - m * time;
      float new_time = (fHitThreshold - b) / m;

      //printf(" old %i - m*x+b = %f*x+%f - new %lf \n", i, m, b, new_time);
      //std::cout << "Old time : " << i << " new time : " << new_time << " \n";

      rising_edge = true;
      hits.push_back(new_time);
    }

    if(rising_edge == true) {

      length_of_hit++;

      if(gradient > fHitThreshold) { 
	if(start_clock == false) {
	  fLenHit->Fill(length_of_hit); 
	}
	start_clock = true; 
      }

      if(start_clock) { width++; }

      if(width > fHitWait) {
	rising_edge = false;
	start_clock = false;
	width = 0;
	length_of_hit = 0;
      }

    }

  }
  
  // Gives a hit of 0 if nothing else was found to prevent seg faults
  if(hits.size() == 0) { hits.push_back(0); }
 
  return hits;

}

std::vector<float> TOFBuilderAlg::get_tof(std::vector<const raw::AuxDetDigit*> ust_wv,std::vector<const raw::AuxDetDigit*> dst_wv) {
  return (get_TOF_and_TimeStamp(ust_wv, dst_wv)).first;
}

std::vector<long> TOFBuilderAlg::get_timeStampDst(std::vector<const raw::AuxDetDigit*> ust_wv,std::vector<const raw::AuxDetDigit*> dst_wv) {
  return (get_TOF_and_TimeStamp(ust_wv, dst_wv)).second;
}

void TOFBuilderAlg::clear_tof_and_timeStampDst()
{
  TOF.clear();
  Dst_Timestamp.clear();
}
