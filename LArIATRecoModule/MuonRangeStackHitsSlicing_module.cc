////////////////////////////////////////////////////////////////////////
// Class:       MuonRangeStackHitsSlicing
// Module Type: producer
// File:        MuonRangeStackHitsSlicing_module.cc
//
// Generated at Thu Jun  4 10:04:24 2015 by Greg Pulliam using artmod
// from cetpkgsupport v1_08_05.
////////////////////////////////////////////////////////////////////////
#ifndef MUONRANGESTACKHITSSLICING_H
#define MUONRANGESTACKHITSSLICING_H


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
#include <TH2F.h>
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "RawDataUtilities/TriggerDigitUtility.h"
#include "LArIATDataProducts/MuonRangeStackHits.h"
#include <map>
#include <vector>
#include <memory>
#include <sstream>
namespace lrm{
	class MuonRangeStackHitsSlicing;
}

class MuonRangeStackHitsSlicing : public art::EDProducer {
public:
  explicit MuonRangeStackHitsSlicing(fhicl::ParameterSet const & p);
  // The destructor generated by the compiler is fine for classes
  // without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  MuonRangeStackHitsSlicing(MuonRangeStackHitsSlicing const &) = delete;
  MuonRangeStackHitsSlicing(MuonRangeStackHitsSlicing &&) = delete;
  MuonRangeStackHitsSlicing & operator = (MuonRangeStackHitsSlicing const &) = delete;
  MuonRangeStackHitsSlicing & operator = (MuonRangeStackHitsSlicing &&) = delete;

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
  void makeTheMuRSTracks( std::map<int, std::vector<int> > fMuonRangeStackMap,
                         std::vector<ldp::MuRSTrack> & MuRSTrackVector,
                         std::vector<size_t> const punchHits );
  void disambiguateTracks( std::vector<ldp::MuRSTrack> & theMuRSTrackVect,
                          std::vector<ldp::MuRSTrack> & finalMuRSTrackVect );
  void trackArchitect( std::vector<std::vector<int> > ptHits,
                      std::vector<std::vector<int> > p1Hits,
                      std::vector<std::vector<int> > p2Hits,
                      std::vector<std::vector<int> > p3Hits,
                      std::vector<std::vector<int> > p4Hits,
                      std::vector<ldp::MuRSTrack> & theMuRSTrackVect );
  void comparePlanes( std::vector<std::vector<int> > & thePlaneVector,
                     std::vector<ldp::MuRSTrack> & aNewMuRSTrackVect );



private:
  std::string fTriggerUtility;
  std::string fSlicerSourceLabel;
  TH1F* fMuRSHitTiming;
  TH1F* fPaddleHits;
  TH1F* fTotalPaddleHits;
  TH2F* fAmpVsPaddle;
  TH1F* fNumPlanesPenetrated;
  TH1F* fMuRS_PT_Logic_Events;

  std::vector<TH1F*> fPlaneOccupancyPerEvent;
  std::vector<TH2F*> fPaddleHitLocationsVsTime;

  int fThreshold;
  size_t fNPaddles;
  size_t fNPlanes;

  std::vector<int> fMuRSPaddleHits;
  std::map<int, std::vector<int> > fMuonRangeStackMap;
  
  size_t fNumberEventsToPlotWFs;
  size_t fEventCounter;
  int fEpsilonTime;
  // Declare member data here.

};


MuonRangeStackHitsSlicing::MuonRangeStackHitsSlicing(fhicl::ParameterSet const & p)
// :
// Initialize member data here.
{
  this->reconfigure(p);
  
  // Call appropriate produces<>() functions here.
  produces<std::vector<ldp::MuonRangeStackHits> >();
}

void MuonRangeStackHitsSlicing::produce(art::Event & e)
{
  // Implementation of required member function here.
  std::unique_ptr<std::vector<ldp::MuonRangeStackHits> > MuonRangeStackCol(new std::vector<ldp::MuonRangeStackHits>  );
  
  //Retrieving AuxDetDigits from event record
  art::Handle<std::vector<raw::AuxDetDigit>> AuxDetDigitHandle;
  e.getByLabel(fSlicerSourceLabel,AuxDetDigitHandle);
  
  //Finding the MuRS and punchthrough digits from the AuxDetDigits
  //Note that the punchthrough is just the trigger logic, but is still useful
  std::vector<raw::AuxDetDigit> MuRSDigits;
  std::vector<raw::AuxDetDigit> PunchthroughDigits;
  for(size_t iDig=0; iDig<AuxDetDigitHandle->size(); ++iDig){
    if(AuxDetDigitHandle->at(iDig).AuxDetName()=="MuonRangeStack"){
      MuRSDigits.push_back(AuxDetDigitHandle->at(iDig));
    }
    if(AuxDetDigitHandle->at(iDig).AuxDetName()=="PUNCH"){
      PunchthroughDigits.push_back(AuxDetDigitHandle->at(iDig));
    }
  }

  //Fill some efficiency histos and kill this event if either the punchthrough or the MuRS digits are nonexistent
  fMuRS_PT_Logic_Events->Fill(0);
  if( PunchthroughDigits.size() == 1 && MuRSDigits.size() == 0 ) fMuRS_PT_Logic_Events->Fill(1);
  else if( PunchthroughDigits.size() == 0 && MuRSDigits.size() > 0 ) fMuRS_PT_Logic_Events->Fill(2);
  else if( PunchthroughDigits.size() == 1 && MuRSDigits.size() > 0 ) fMuRS_PT_Logic_Events->Fill(3);
  //Actually, don't kill. We're testing something right meow to see if providing empty MuRS objects increases efficiency
  if( PunchthroughDigits.size() == 0 || MuRSDigits.size() == 0 ){
    std::vector<ldp::MuRSTrack> mursTrackVector;
    std::map<int,std::vector<int> > theMap;
    ldp::MuonRangeStackHits the_MuRS(theMap,mursTrackVector);
    (*MuonRangeStackCol).push_back(the_MuRS);		
    e.put(std::move(MuonRangeStackCol));	   
    return;
  }

  ////////////////// PUNCHTHROUGH INFO GETTING ///////////////////
  //Sanity check: should only be one digit for the Punchthrough
  if( PunchthroughDigits.size() > 1 ) std::cout << "WARNING: MORE THAN ONE PUNCHTHROUGH DIGIT." << std::endl;
  
  //Loop through the punchthrough and get a vector of its hits
  std::vector<size_t> punchHits;
  
  for( size_t iADC = 0; iADC < PunchthroughDigits.at(0).NADC() ; ++iADC ){
    if( PunchthroughDigits.at(0).ADC(iADC) < fThreshold ){
      punchHits.push_back(iADC);
    }
  }

  ////////////////// MURS INFO GETTING ///////////////////  
  int size=MuRSDigits.size();
  std::cout << "Size of MuRSDigits: " << size << std::endl;
  /*  
  //Debug sanity check - only temporary
  bool isThereAHit = false;
  if( size == 16 ){
    int TrigMult=size/16;
    for (int TrigIter=0; TrigIter<TrigMult; ++TrigIter){
      for (int nPaddle=TrigIter*16; nPaddle<(TrigIter+1)*16; ++nPaddle){
	auto PaddleDigit=MuRSDigits[nPaddle];
	for (size_t i=0; i<PaddleDigit.NADC(); ++i){
	  if(PaddleDigit.ADC(i)<fThreshold){
	    isThereAHit = true;
	  }
	}
      }
    }
  }
  if( isThereAHit == false ) return;
  std::cout << "Event number: " << e.event() << ", histo: " << fEventCounter << std::endl;
  */

  //Now doing legitimate filling
  int TrigMult=size/16;
  for (int TrigIter=0; TrigIter<TrigMult; ++TrigIter){
    for (int nPaddle=TrigIter*16; nPaddle<(TrigIter+1)*16; ++nPaddle){
      auto PaddleDigit=MuRSDigits[nPaddle];
      for (size_t i=0; i<PaddleDigit.NADC(); ++i){
        if(PaddleDigit.ADC(i)<fThreshold){
          fMuRSPaddleHits.push_back(i);
          fMuRSHitTiming->Fill(i);
          fAmpVsPaddle->Fill(nPaddle-TrigIter*16,fThreshold-PaddleDigit.ADC(i));
          fTotalPaddleHits->Fill(nPaddle-TrigIter*16);
          fPaddleHits->Fill(nPaddle-TrigIter*16);
        }
      }
      fMuonRangeStackMap.emplace(nPaddle-TrigIter*16,fMuRSPaddleHits);
      fMuRSPaddleHits.clear();
    }
    std::vector<ldp::MuRSTrack> mursTrackVector;
    makeTheMuRSTracks(fMuonRangeStackMap,mursTrackVector,punchHits);
    
    ldp::MuonRangeStackHits the_MuRS(fMuonRangeStackMap,mursTrackVector);
    (*MuonRangeStackCol).push_back(the_MuRS);
    fMuonRangeStackMap.clear();
  }

  e.put(std::move(MuonRangeStackCol));	   
  fEventCounter++;
}//produce()

//Sort hits in the MuRS by time and make tracks
void MuonRangeStackHitsSlicing::makeTheMuRSTracks( std::map<int, std::vector<int> > MuonRangeStackMap,
						   std::vector<ldp::MuRSTrack> & finalMuRSTrackVect,
						   std::vector<size_t> const punchHits )
{
  LOG_DEBUG("MuonRangeStackHitsSlicing")
  << "makeTheMuRSTracks called.";

  //Filling some histos about plane multiplicity per event
  std::vector<int> planeOccVect;
  for( size_t iPlane = 0; iPlane < fNPlanes ; ++iPlane )
    planeOccVect.push_back(0);
  for( size_t iPaddle = 0; iPaddle < fNPaddles; ++iPaddle ){
    int plane = floor( iPaddle/4 );
    std::vector<int> theMuRSHitTimes = MuonRangeStackMap.at(iPaddle);
    planeOccVect.at(plane) += theMuRSHitTimes.size();
  }
  for( size_t iPlane = 0; iPlane < fNPlanes ; ++iPlane )
    fPlaneOccupancyPerEvent.at(iPlane)->Fill(planeOccVect.at(iPlane));

  //Filling histos to show some events' hits' spatial vs. temporal locations
  if( fEventCounter < fNumberEventsToPlotWFs ){
    LOG_DEBUG("MuonRangeStackHitsSlicing")
    << "fEventCounter is small enough to fill histos.";
    for( size_t iPaddle = 0; iPaddle < fNPaddles; ++iPaddle ){
      for( size_t iHit = 0; iHit < MuonRangeStackMap.at(iPaddle).size(); ++iHit ){
        LOG_VERBATIM("MuonRangeStackHitsSlicing")
        << "Filling with paddle: "
        << iPaddle
        << ", time: "
        << MuonRangeStackMap.at(iPaddle)[iHit];
        fPaddleHitLocationsVsTime.at(fEventCounter)->Fill(MuonRangeStackMap.at(iPaddle)[iHit],iPaddle);
      }
    }
    //Fill the punchthrough info as the the -1th paddle.
    for( size_t iPTHit = 0; iPTHit < punchHits.size(); ++iPTHit ){
      fPaddleHitLocationsVsTime.at(fEventCounter)->Fill(punchHits.at(iPTHit),-1);
    }
  }

  //Pseudocode for track building
  //For each plane
  // Make vector of vectors: 1.) Plane, 2.) paddle, 3.) time

  //Creating vectors of hits
  std::vector<std::vector<int> > punchThroughHits;
  std::vector<std::vector<int> > plane1Hits;
  std::vector<std::vector<int> > plane2Hits;
  std::vector<std::vector<int> > plane3Hits;
  std::vector<std::vector<int> > plane4Hits;
  
  //Looping through the track vectors to sort the hits into planes
  for( size_t iPaddle = 0; iPaddle < fNPaddles ; ++iPaddle ){
    int plane = floor( iPaddle/4 );
    for( size_t iHit = 0; iHit < MuonRangeStackMap.at(iPaddle).size() ; ++iHit ){
      //Creating the hit
      std::vector<int> theHit;
      theHit.push_back(plane);
      theHit.push_back(iPaddle);
      theHit.push_back(MuonRangeStackMap.at(iPaddle)[iHit]);
      
      //Now pushing it into the right vector
      if( plane == 0 ) plane1Hits.push_back(theHit);
      if( plane == 1 ) plane2Hits.push_back(theHit);
      if( plane == 2 ) plane3Hits.push_back(theHit);
      if( plane == 3 ) plane4Hits.push_back(theHit);
    }
  }
  for( size_t iPTHit = 0; iPTHit < punchHits.size(); ++iPTHit ){
    std::vector<int> theHit;
    theHit.push_back(-1); //By convention, -1 for punch through
    theHit.push_back(-1); //By convention, 0 for punch through triggering (we don't have paddle info for the PT)
    theHit.push_back(punchHits.at(iPTHit));
    punchThroughHits.push_back(theHit);
  }

  //Now we have five vectors of hits represenging the 5 planes of scintillators.
  std::vector<ldp::MuRSTrack> theMuRSTrackVect;
  trackArchitect( punchThroughHits,
		  plane1Hits,
		  plane2Hits,
		  plane3Hits,
		  plane4Hits,
		  theMuRSTrackVect );

  //Now remove the murs tracks that have only one hit
  for( size_t iTrack = 0; iTrack < theMuRSTrackVect.size() ; ++iTrack ){
    if( theMuRSTrackVect[iTrack].HitVect.size() <= 1 ){
      theMuRSTrackVect.erase(theMuRSTrackVect.begin()+iTrack);
      iTrack--;
    }
  }

  //Now track disambiguation
  disambiguateTracks(theMuRSTrackVect,finalMuRSTrackVect);

  //Print out the rest
  for( size_t iTrack = 0; iTrack < finalMuRSTrackVect.size() ; ++iTrack ){
    LOG_VERBATIM("MuonRangeStackSlicing")
    << "************** MuRS TRACK "
    << iTrack
    << " ****************";
    for( size_t iHit = 0; iHit < finalMuRSTrackVect[iTrack].HitVect.size() ; ++iHit ){
      LOG_VERBATIM("MuonRangeStackSlicing")
      << "Plane: " << finalMuRSTrackVect[iTrack].HitVect[iHit].at(0)
      << ", Paddle: " << finalMuRSTrackVect[iTrack].HitVect[iHit].at(1)
      << ", Time: " << finalMuRSTrackVect[iTrack].HitVect[iHit].at(2);
    }
  }


}


//================================================================================================
//Some tracks will have multiple hits in the same plane. Make a track for each possible combination of these
void MuonRangeStackHitsSlicing::disambiguateTracks( std::vector<ldp::MuRSTrack> & theMuRSTrackVect,
                                                   std::vector<ldp::MuRSTrack> & finalMuRSTrackVect )
{
  //First we have to merge multiple hits in the same paddle that fall
  //within the allotted time window for trackbuilding. "Merge" means take the
  //lowest-time one.
  for( size_t iTrack = 0 ; iTrack < theMuRSTrackVect.size() ; ++iTrack ){
    int paddleCounts[17] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};    
    for( size_t iHit = 0 ; iHit < theMuRSTrackVect[iTrack].HitVect.size() ; ++iHit ){
      if( paddleCounts[theMuRSTrackVect[iTrack].HitVect[iHit].at(1)+1] == 1 ){
        theMuRSTrackVect[iTrack].HitVect.erase(theMuRSTrackVect[iTrack].HitVect.begin()+iHit);
        --iHit;
        continue;
      }
      paddleCounts[theMuRSTrackVect[iTrack].HitVect[iHit].at(1)+1]++;
    }
  }
  /*
  //Debug printing
  for( size_t iTrack = 0; iTrack < theMuRSTrackVect.size() ; ++iTrack ){
    std::cout << "************** MuRS TRACK " << iTrack << " ****************" << std::endl;
    for( size_t iHit = 0; iHit < theMuRSTrackVect[iTrack].HitVect.size() ; ++iHit ){
      std::cout << "Plane: " << theMuRSTrackVect[iTrack].HitVect[iHit].at(0)
		<< ", Paddle: " << theMuRSTrackVect[iTrack].HitVect[iHit].at(1) 
		<< ", Time: " << theMuRSTrackVect[iTrack].HitVect[iHit].at(2) << std::endl;
    }
  }
  */

    
  //Now do the disambiguation
  for( size_t iTrack = 0 ; iTrack < theMuRSTrackVect.size() ; ++iTrack ){
    std::vector<std::vector<int> > ptHits;
    std::vector<std::vector<int> > p1Hits;
    std::vector<std::vector<int> > p2Hits;
    std::vector<std::vector<int> > p3Hits;
    std::vector<std::vector<int> > p4Hits;
    for( size_t iHit = 0; iHit < theMuRSTrackVect[iTrack].HitVect.size(); ++iHit ){
      if( theMuRSTrackVect[iTrack].HitVect[iHit].at(0) == -1 ) ptHits.push_back(theMuRSTrackVect[iTrack].HitVect[iHit]);
      if( theMuRSTrackVect[iTrack].HitVect[iHit].at(0) == 0 ) p1Hits.push_back(theMuRSTrackVect[iTrack].HitVect[iHit]);
      if( theMuRSTrackVect[iTrack].HitVect[iHit].at(0) == 1 ) p2Hits.push_back(theMuRSTrackVect[iTrack].HitVect[iHit]);
      if( theMuRSTrackVect[iTrack].HitVect[iHit].at(0) == 2 ) p3Hits.push_back(theMuRSTrackVect[iTrack].HitVect[iHit]);
      if( theMuRSTrackVect[iTrack].HitVect[iHit].at(0) == 3 ) p4Hits.push_back(theMuRSTrackVect[iTrack].HitVect[iHit]);
    }
    //Now fill empties with a fake hit for the permutations. We'll delete them later.
    std::vector<int> fakeHit;
    if( ptHits.size() == 0 ) ptHits.push_back(fakeHit);
    if( p1Hits.size() == 0 ) p1Hits.push_back(fakeHit);
    if( p2Hits.size() == 0 ) p2Hits.push_back(fakeHit);
    if( p3Hits.size() == 0 ) p3Hits.push_back(fakeHit);
    if( p4Hits.size() == 0 ) p4Hits.push_back(fakeHit);
    
    //Now create new tracks for all possible combos of hits
    std::vector<ldp::MuRSTrack> tempMuRSTrackVect;
    for( size_t iPTH = 0; iPTH < ptHits.size() ; ++iPTH ){
      for( size_t iP1H = 0; iP1H < p1Hits.size() ; ++iP1H ){
        for( size_t iP2H = 0; iP2H < p2Hits.size() ; ++iP2H ){
          for( size_t iP3H = 0; iP3H < p3Hits.size() ; ++iP3H ){
            for( size_t iP4H = 0; iP4H < p4Hits.size() ; ++iP4H ){
              ldp::MuRSTrack theNewTrack;
              theNewTrack.HitVect.push_back(ptHits.at(iPTH));
              theNewTrack.HitVect.push_back(p1Hits.at(iP1H));
              theNewTrack.HitVect.push_back(p2Hits.at(iP2H));
              theNewTrack.HitVect.push_back(p3Hits.at(iP3H));
              theNewTrack.HitVect.push_back(p4Hits.at(iP4H));
              tempMuRSTrackVect.push_back(theNewTrack);
            }
          }
        }
      }
    }
    
    //Kill the fake hits
    for( size_t iTrack = 0; iTrack < tempMuRSTrackVect.size(); ++iTrack ){
      for( size_t iHit = 0; iHit < tempMuRSTrackVect[iTrack].HitVect.size() ; ++iHit ){
        if( tempMuRSTrackVect[iTrack].HitVect[iHit].size() == 0 ){
          tempMuRSTrackVect[iTrack].HitVect.erase(tempMuRSTrackVect[iTrack].HitVect.begin()+iHit);
          iHit--;
        }
      }
    }

    //Finally (and I do mean finally), push the remaining tracks into the final track vector.
    for( size_t iTrack = 0; iTrack < tempMuRSTrackVect.size(); ++iTrack )
      finalMuRSTrackVect.push_back(tempMuRSTrackVect[iTrack]);
   
    //(Whoops, I guess I didn't mean finally...) Clear the temp vector
    tempMuRSTrackVect.clear();
  }
}


//================================================================================================
//This does the legwork for matching hits that are close in time to form tracks
void MuonRangeStackHitsSlicing::trackArchitect( std::vector<std::vector<int> > ptHits,
						std::vector<std::vector<int> > p1Hits,
						std::vector<std::vector<int> > p2Hits,
						std::vector<std::vector<int> > p3Hits,
						std::vector<std::vector<int> > p4Hits,
						std::vector<ldp::MuRSTrack> & theMuRSTrackVect )
{
  //PSEUDOCODE!!!!
  //Do the below process, but for punch through first
  //For all hits in first plane
  // start a track vector with each hit as the beginning
  //For all track vectors, 
  // For all hits in second plane
  //  if there is a hit close enough to any hit in the newly made vector, push that vector back with the hit
  //  if not, then make a new track vector
  //For all (updated) track vectors
  // For all hits in the third plane,
  //  if there is a hit close enough to any hit in the newly made vector, push that vector back with the hit
  //  if not, then make a new track vector
  //For all (updated) track vectors
  // For all hits in the fourth plane
  //  if there is a hit close enough to any hit in the newly made vector, push that vector back with the hit
  //  if not, then make a new track vector
  
  std::vector<ldp::MuRSTrack> aNewMuRSTrackVect;
  //For all hits in punchthrough plane
  for( size_t iPTHit = 0; iPTHit < ptHits.size() ; ++iPTHit ){
    ldp::MuRSTrack theTrack;
    std::vector<int> theHit;
    theHit.push_back(ptHits.at(iPTHit).at(0));
    theHit.push_back(ptHits.at(iPTHit).at(1));
    theHit.push_back(ptHits.at(iPTHit).at(2));		  
    //make a new track vector, with each hit as the beginning entry
    theTrack.HitVect.push_back(theHit);
    aNewMuRSTrackVect.push_back(theTrack);
  }
  
  //Now do the comparisons to look for hit agreement
  comparePlanes( p1Hits, aNewMuRSTrackVect );
  comparePlanes( p2Hits, aNewMuRSTrackVect );
  comparePlanes( p3Hits, aNewMuRSTrackVect );
  comparePlanes( p4Hits, aNewMuRSTrackVect );
  theMuRSTrackVect = aNewMuRSTrackVect;
}						

//================================================================================================
void MuonRangeStackHitsSlicing::comparePlanes( std::vector<std::vector<int> > & thePlaneVector,
					       std::vector<ldp::MuRSTrack> & aNewMuRSTrackVect )
{
  //For all hits in the first plane, and for all tracks, compare hit times to check for agreement
  for( size_t iHit = 0; iHit < thePlaneVector.size() ; ++iHit ){
    bool isAbsorbed = false;
    for( size_t iTrack = 0; iTrack < aNewMuRSTrackVect.size() ; ++iTrack ){
      if( fabs(thePlaneVector[iHit].at(2) - aNewMuRSTrackVect[iTrack].HitVect.at(0).at(2)) < fEpsilonTime ){
        aNewMuRSTrackVect[iTrack].HitVect.push_back(thePlaneVector[iHit]);
        thePlaneVector.erase(thePlaneVector.begin()+iHit);
        iHit--; //Do so we don't skip any hits when the vector gets shrunk by the above line.
        isAbsorbed = true;
        break;
      }
    }
      //If a hit is not absorbed by an existing track, check to see if it is absorbed
      //by an existing track in this plane. If not, a new track is made
    if( !isAbsorbed ){
      ldp::MuRSTrack theTrack;
      theTrack.HitVect.push_back(thePlaneVector[iHit]);
      aNewMuRSTrackVect.push_back(theTrack);
      thePlaneVector.erase(thePlaneVector.begin()+iHit);
      --iHit;
    }
  }
}
/*
//================================================================================================
//Check the hit against all other hits in this plane to see if it belongs to a track
bool MuonRangeStackHitsSlicing::lastCheck( size_t theHit,
					   std::vector<std::vector<int> > & thePlaneVector,
					   std::vector<ldp::MuRSTrack> & aNewMuRSTrackVect)
{
  
  for( size_t iHit = 0; iHit < thePlaneVector.size(); ++iHit ){
    if( iHit == theHit ) continue;
    if( fabs(thePlaneVector[iHit].at(2) - thePlaneVector.at(theHit).at(2)) < fEpsilonTime ){
      aNewMuRSTrackVect[iTrack].HitVect.push_back(thePlaneVector[iHit]);
      thePlaneVector.erase(thePlaneVector.begin()+iHit);
      iHit--; //Do so we don't skip any hits when the vector gets shrunk by the above line.
      



}
*/

//================================================================================================
void MuonRangeStackHitsSlicing::beginJob()
{
  art::ServiceHandle<art::TFileService> tfs;
  fMuRSHitTiming = tfs->make<TH1F>("MuRSHitTicks","MuRSHitTicks", 3073, 0, 3073);
  fMuRS_PT_Logic_Events = tfs->make<TH1F>("MuRS_PT_Logic_Events","0: All events, 1: Only PT, 2: Only MuRS, 3: Both",4,0,4);
  fPaddleHits=tfs->make<TH1F>("PaddleHits", "PaddleHits",16, 0, 16);
  fTotalPaddleHits=tfs->make<TH1F>("TotalPaddleHits","TotalPaddleHits",16, 0, 16);
  fAmpVsPaddle=tfs->make<TH2F>("AmplitudeVsPaddle", "AmplitudeVsPaddle",16,0,16,2021,0,2021);
  fNumPlanesPenetrated=tfs->make<TH1F>("NumPlanesPenetrated","Number of Planes penetrated in event",4,0,4);

  for( size_t iPlane = 0; iPlane < fNPlanes; ++iPlane ){
    char name[40];
    sprintf(name,"Plane %d Hit Occupancy Per Event",int(iPlane));
    fPlaneOccupancyPerEvent.push_back(tfs->make<TH1F>(name,name,30,0,30));
  }

  for( size_t iEvt = 0; iEvt < fNumberEventsToPlotWFs ; ++iEvt ){
    char name[40];
    sprintf(name,"Event %d Hit Time vs. Paddle",int(iEvt));
    fPaddleHitLocationsVsTime.push_back(tfs->make<TH2F>(name,name,3072,0,3072,17,-1,16));
  }

  // Implementation of optional member function here.
}

void MuonRangeStackHitsSlicing::beginRun(art::Run & r)
{
  // Implementation of optional member function here.
}

void MuonRangeStackHitsSlicing::beginSubRun(art::SubRun & sr)
{

  // Implementation of optional member function here.
}

void MuonRangeStackHitsSlicing::endJob()
{
  // Implementation of optional member function here.
}

void MuonRangeStackHitsSlicing::endRun(art::Run & r)
{
  // Implementation of optional member function here.
}

void MuonRangeStackHitsSlicing::endSubRun(art::SubRun & sr)
{
  // Implementation of optional member function here.
}

void MuonRangeStackHitsSlicing::reconfigure(fhicl::ParameterSet const & p)
{
  fSlicerSourceLabel=p.get<std::string>("SourceLabel");
  fThreshold=p.get<int>("Threshold");
  fNPaddles = p.get<size_t>("NumPaddles",16);
  fNPlanes = p.get<size_t>("NumPlanes",4);
  fNumberEventsToPlotWFs = p.get<size_t>("NumberEventsToPlotWFs",50);
  fEventCounter = 0;
  fEpsilonTime = p.get<int>("EpsilonTime",3);
}

void MuonRangeStackHitsSlicing::respondToCloseInputFile(art::FileBlock const & fb)
{
  // Implementation of optional member function here.
}

void MuonRangeStackHitsSlicing::respondToCloseOutputFiles(art::FileBlock const & fb)
{
  // Implementation of optional member function here.
}

void MuonRangeStackHitsSlicing::respondToOpenInputFile(art::FileBlock const & fb)
{
  // Implementation of optional member function here.
}

void MuonRangeStackHitsSlicing::respondToOpenOutputFiles(art::FileBlock const & fb)
{
  // Implementation of optional member function here.
}

DEFINE_ART_MODULE(MuonRangeStackHitsSlicing)

#endif
