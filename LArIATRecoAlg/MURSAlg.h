////////////////////////////////////////////////////////////////
//                                                            //
// This is a class definition for the muon range stack        //
//         algorithm,                                         //
//                                                            //
// Authors: Pawel Kryczynski pkryczyn@fnal.gov, based on code //
// by Elena Gramellini elena.gramellini@yale.edu  and         //
// Greg Pulliam                                               //
////////////////////////////////////////////////////////////////


#ifndef MURSALG_H
#define MURSALG_H

//C++ includes
#include <vector>
#include <cmath>
#include <iostream>
#include <utility>

//Framework includes
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"

// LArSoft includes
#include "larcore/Geometry/Geometry.h"
#include "lardata/RawData/AuxDetDigit.h"

//ROOT
#include <TH1F.h>



//--------------------------------------------
class MURSAlg{
 public:
  
  //Constructor/destructor
  MURSAlg( fhicl::ParameterSet const& pset );
  ~MURSAlg();

  void reconfigure( fhicl::ParameterSet const& pset );
  void makeTheMuRSTracks( std::map<int, std::vector<int> > MuonRangeStackMap,
						   std::vector<MuRSTrack> & finalMuRSTrackVect,
						   std::vector<size_t> const punchHits );
  void disambiguateTracks( std::vector<MuRSTrack> & theMuRSTrackVect,
						    std::vector<MuRSTrack> & finalMuRSTrackVect);

  void trackArchitect( std::vector<std::vector<int> > ptHits,
						std::vector<std::vector<int> > p1Hits,
						std::vector<std::vector<int> > p2Hits,
						std::vector<std::vector<int> > p3Hits,
						std::vector<std::vector<int> > p4Hits,
						std::vector<MuRSTrack> & theMuRSTrackVect );

  void comparePlanes( std::vector<std::vector<int> > & thePlaneVector,
					       std::vector<MuRSTrack> & aNewMuRSTrackVect);


 private:
  
  std::string fSlicerSourceLabel;
  int fThreshold;
  bool fVerbose;
  size_t fNPaddles;
  size_t fNPlanes;
  size_t fNumberEventsToPlotWFs;
  int fEventCounter;
  int fEpsilonTime;
  
};


#endif
