/////////////////////////////////////////////////////////////////
//                                                             //
// This is a class definition for the muon range stack hit     //
// finding algorithm                                           //  
//                                                             //
// Authors: Pawel Kryczynski pkryczyn@fnal.gov, but he's just  //                           
// adapted the code by Greg Pulliam                            //
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
#include "LArIATDataProducts/MuonRangeStackHits.h"
#include "LArIATRecoAlg/MURSAlg.h"
#include <iostream>
#include <cmath>
#include <cstdlib>

//--------------------------------------------------------------
//Constructor
MURSAlg::MURSAlg( fhicl::ParameterSet const& pset )
{

  this->reconfigure(pset);
  art::ServiceHandle<art::TFileService> tfs;

}

//--------------------------------------------------------------  
//Destructor
MURSAlg::~MURSAlg()
{

}

//--------------------------------------------------------------
void MURSAlg::reconfigure( fhicl::ParameterSet const& p )
{
  fSlicerSourceLabel=p.get<std::string>("SourceLabel");
  fThreshold=p.get<int>("Threshold");
  fVerbose = p.get<bool>("Verbose");
  fNPaddles = p.get<size_t>("NumPaddles",16);
  fNPlanes = p.get<size_t>("NumPlanes",4);
  fNumberEventsToPlotWFs = p.get<size_t>("NumberEventsToPlotWFs",50);
  fEventCounter = 0;
  fEpsilonTime = p.get<int>("EpsilonTime",3);

}






//--------------------------------------------------------------
//Sort hits in the MuRS by time and make tracks
void MURSAlg::makeTheMuRSTracks( std::map<int, std::vector<int> > MuonRangeStackMap,
						   std::vector<MuRSTrack> & finalMuRSTrackVect,
						   std::vector<size_t> const punchHits )
{
  if( fVerbose ) std::cout << "makeTheMuRSTracks called." << std::endl;

  //Filling some histos about plane multiplicity per event
  std::vector<int> planeOccVect;
  for( size_t iPlane = 0; iPlane < fNPlanes ; ++iPlane )
    planeOccVect.push_back(0);
  for( size_t iPaddle = 0; iPaddle < fNPaddles; ++iPaddle ){
    int plane = floor( iPaddle/4 );
    std::vector<int> theMuRSHitTimes = MuonRangeStackMap.at(iPaddle);
    planeOccVect.at(plane) += theMuRSHitTimes.size();
  }
  //for( size_t iPlane = 0; iPlane < fNPlanes ; ++iPlane )
    //fPlaneOccupancyPerEvent.at(iPlane)->Fill(planeOccVect.at(iPlane));

  //Filling histos to show some events' hits' spatial vs. temporal locations
 /* if( fEventCounter < fNumberEventsToPlotWFs ){
    if( fVerbose ) std::cout << "fEventCounter is small enough to fill histos." << std::endl;
    for( size_t iPaddle = 0; iPaddle < fNPaddles; ++iPaddle ){
      for( size_t iHit = 0; iHit < MuonRangeStackMap.at(iPaddle).size(); ++iHit ){
	//std::cout << "Filling with paddle: " << iPaddle << ", time: " << MuonRangeStackMap.at(iPaddle).at(iHit) << std::endl;
	//fPaddleHitLocationsVsTime.at(fEventCounter)->Fill(MuonRangeStackMap.at(iPaddle).at(iHit),iPaddle);
      }
    }
    //Fill the punchthrough info as the the -1th paddle.
    for( size_t iPTHit = 0; iPTHit < punchHits.size(); ++iPTHit ){
      fPaddleHitLocationsVsTime.at(fEventCounter)->Fill(punchHits.at(iPTHit),-1);
    }
  } */

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
      theHit.push_back(MuonRangeStackMap.at(iPaddle).at(iHit));
      
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
  std::vector<MuRSTrack> theMuRSTrackVect;
  trackArchitect( punchThroughHits,
		  plane1Hits,
		  plane2Hits,
		  plane3Hits,
		  plane4Hits,
		  theMuRSTrackVect );

  //Now remove the murs tracks that have only one hit
  for( size_t iTrack = 0; iTrack < theMuRSTrackVect.size() ; ++iTrack ){
    if( theMuRSTrackVect.at(iTrack).HitVect.size() <= 1 ){
      theMuRSTrackVect.erase(theMuRSTrackVect.begin()+iTrack);
      iTrack--;
    }
  }

  //Now track disambiguation
  disambiguateTracks(theMuRSTrackVect,finalMuRSTrackVect);

  //Print out the rest
  for( size_t iTrack = 0; iTrack < finalMuRSTrackVect.size() ; ++iTrack ){
    std::cout << "************** MuRS TRACK " << iTrack << " ****************" << std::endl;
    for( size_t iHit = 0; iHit < finalMuRSTrackVect.at(iTrack).HitVect.size() ; ++iHit ){
      std::cout << "Plane: " << finalMuRSTrackVect.at(iTrack).HitVect.at(iHit).at(0)
		<< ", Paddle: " << finalMuRSTrackVect.at(iTrack).HitVect.at(iHit).at(1) 
		<< ", Time: " << finalMuRSTrackVect.at(iTrack).HitVect.at(iHit).at(2) << std::endl;
    }
  }


}


//================================================================================================
//Some tracks will have multiple hits in the same plane. Make a track for each possible combination of these
void MURSAlg::disambiguateTracks( std::vector<MuRSTrack> & theMuRSTrackVect,
						    std::vector<MuRSTrack> & finalMuRSTrackVect )
{
  //First we have to merge multiple hits in the same paddle that fall
  //within the allotted time window for trackbuilding. "Merge" means take the
  //lowest-time one.
  for( size_t iTrack = 0 ; iTrack < theMuRSTrackVect.size() ; ++iTrack ){
    int paddleCounts[17] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};    
    for( size_t iHit = 0 ; iHit < theMuRSTrackVect.at(iTrack).HitVect.size() ; ++iHit ){
      if( paddleCounts[theMuRSTrackVect.at(iTrack).HitVect.at(iHit).at(1)+1] == 1 ){
	theMuRSTrackVect.at(iTrack).HitVect.erase(theMuRSTrackVect.at(iTrack).HitVect.begin()+iHit);
	--iHit;
	continue;
      }
      paddleCounts[theMuRSTrackVect.at(iTrack).HitVect.at(iHit).at(1)+1]++;
    }
  }
  /*
  //Debug printing
  for( size_t iTrack = 0; iTrack < theMuRSTrackVect.size() ; ++iTrack ){
    std::cout << "************** MuRS TRACK " << iTrack << " ****************" << std::endl;
    for( size_t iHit = 0; iHit < theMuRSTrackVect.at(iTrack).HitVect.size() ; ++iHit ){
      std::cout << "Plane: " << theMuRSTrackVect.at(iTrack).HitVect.at(iHit).at(0)
		<< ", Paddle: " << theMuRSTrackVect.at(iTrack).HitVect.at(iHit).at(1) 
		<< ", Time: " << theMuRSTrackVect.at(iTrack).HitVect.at(iHit).at(2) << std::endl;
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
    for( size_t iHit = 0; iHit < theMuRSTrackVect.at(iTrack).HitVect.size(); ++iHit ){
      if( theMuRSTrackVect.at(iTrack).HitVect.at(iHit).at(0) == -1 ) ptHits.push_back(theMuRSTrackVect.at(iTrack).HitVect.at(iHit));
      if( theMuRSTrackVect.at(iTrack).HitVect.at(iHit).at(0) == 0 ) p1Hits.push_back(theMuRSTrackVect.at(iTrack).HitVect.at(iHit));
      if( theMuRSTrackVect.at(iTrack).HitVect.at(iHit).at(0) == 1 ) p2Hits.push_back(theMuRSTrackVect.at(iTrack).HitVect.at(iHit));
      if( theMuRSTrackVect.at(iTrack).HitVect.at(iHit).at(0) == 2 ) p3Hits.push_back(theMuRSTrackVect.at(iTrack).HitVect.at(iHit));
      if( theMuRSTrackVect.at(iTrack).HitVect.at(iHit).at(0) == 3 ) p4Hits.push_back(theMuRSTrackVect.at(iTrack).HitVect.at(iHit));
    }
    //Now fill empties with a fake hit for the permutations. We'll delete them later.
    std::vector<int> fakeHit;
    if( ptHits.size() == 0 ) ptHits.push_back(fakeHit);
    if( p1Hits.size() == 0 ) p1Hits.push_back(fakeHit);
    if( p2Hits.size() == 0 ) p2Hits.push_back(fakeHit);
    if( p3Hits.size() == 0 ) p3Hits.push_back(fakeHit);
    if( p4Hits.size() == 0 ) p4Hits.push_back(fakeHit);
    
    //Now create new tracks for all possible combos of hits
    std::vector<MuRSTrack> tempMuRSTrackVect;
    for( size_t iPTH = 0; iPTH < ptHits.size() ; ++iPTH ){
      for( size_t iP1H = 0; iP1H < p1Hits.size() ; ++iP1H ){
	for( size_t iP2H = 0; iP2H < p2Hits.size() ; ++iP2H ){
	  for( size_t iP3H = 0; iP3H < p3Hits.size() ; ++iP3H ){
	    for( size_t iP4H = 0; iP4H < p4Hits.size() ; ++iP4H ){
	      MuRSTrack theNewTrack;
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
      for( size_t iHit = 0; iHit < tempMuRSTrackVect.at(iTrack).HitVect.size() ; ++iHit ){
	if( tempMuRSTrackVect.at(iTrack).HitVect.at(iHit).size() == 0 ){
	  tempMuRSTrackVect.at(iTrack).HitVect.erase(tempMuRSTrackVect.at(iTrack).HitVect.begin()+iHit);
	  iHit--;
	}
      }
    }

    //Finally (and I do mean finally), push the remaining tracks into the final track vector.
    for( size_t iTrack = 0; iTrack < tempMuRSTrackVect.size(); ++iTrack )
      finalMuRSTrackVect.push_back(tempMuRSTrackVect.at(iTrack));
   
    //(Whoops, I guess I didn't mean finally...) Clear the temp vector
    tempMuRSTrackVect.clear();
  }
}


//================================================================================================
//This does the legwork for matching hits that are close in time to form tracks
void MURSAlg::trackArchitect( std::vector<std::vector<int> > ptHits,
						std::vector<std::vector<int> > p1Hits,
						std::vector<std::vector<int> > p2Hits,
						std::vector<std::vector<int> > p3Hits,
						std::vector<std::vector<int> > p4Hits,
						std::vector<MuRSTrack> & theMuRSTrackVect )
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
  
  std::vector<MuRSTrack> aNewMuRSTrackVect;
  //For all hits in punchthrough plane
  for( size_t iPTHit = 0; iPTHit < ptHits.size() ; ++iPTHit ){
    MuRSTrack theTrack;
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
void MURSAlg::comparePlanes( std::vector<std::vector<int> > & thePlaneVector,
					       std::vector<MuRSTrack> & aNewMuRSTrackVect )
{
  //For all hits in the first plane, and for all tracks, compare hit times to check for agreement
  for( size_t iHit = 0; iHit < thePlaneVector.size() ; ++iHit ){
    bool isAbsorbed = false;
    for( size_t iTrack = 0; iTrack < aNewMuRSTrackVect.size() ; ++iTrack ){
      if( fabs(thePlaneVector.at(iHit).at(2) - aNewMuRSTrackVect.at(iTrack).HitVect.at(0).at(2)) < fEpsilonTime ){
	aNewMuRSTrackVect.at(iTrack).HitVect.push_back(thePlaneVector.at(iHit));
	thePlaneVector.erase(thePlaneVector.begin()+iHit);
	iHit--; //Do so we don't skip any hits when the vector gets shrunk by the above line.
	isAbsorbed = true;
	break;
      }
    }
    //If a hit is not absorbed by an existing track, check to see if it is absorbed
    //by an existing track in this plane. If not, a new track is made
    if( !isAbsorbed ){
      MuRSTrack theTrack;
      theTrack.HitVect.push_back(thePlaneVector.at(iHit));
      aNewMuRSTrackVect.push_back(theTrack);
      thePlaneVector.erase(thePlaneVector.begin()+iHit);
      iHit--;
    }
  }
}
