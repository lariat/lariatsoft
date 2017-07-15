////////////////////////////////////////////////////////////////////////
// This alg takes in TDC information for the LArIAT Wire Chambers     //
// and returns a vector of "good_hits" that can be used in            //
// finding good tracks later, using a track building alg of your      //
// choice.  In effect this is just the hit finding parts of the       //
// code in WCTrackBuildingSlicing_module.cc and WCTrackBuilderAlg.cxx //
// Author: Greg Pulliam, gkpullia@syr.edu                             //
////////////////////////////////////////////////////////////////////////

//Framework includes
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "CLHEP/Units/SystemOfUnits.h"
#include "larcorealg/Geometry/AuxDetGeo.h"
#include "art/Framework/Services/Optional/TFileService.h"

// LArIAT includes
#include "LArIATRecoAlg/WCHitFinderAlg.h"


#include <iostream>
#include <cmath>
#include <cstdlib>
#include <TH1F.h>
#include <string>

//-------------------------------------------------------------------
//Constructor
WCHitFinderAlg::WCHitFinderAlg( fhicl::ParameterSet const& pset )
{
  this->reconfigure(pset);
}
//-------------------------------------------------------------------
//Destructor
WCHitFinderAlg::~WCHitFinderAlg()
{

}
//-------------------------------------------------------------------
void WCHitFinderAlg::reconfigure( fhicl::ParameterSet const& pset )
{
  fNumber_tdcs          = pset.get<int   >("NumberTDCs",         16         );
  fNumber_wire_chambers = pset.get<int   >("NumberWireChambers", 4          );
  fTime_bin_scaling     = pset.get<double>("TimeBinScaling",      1.0/1280.0);
  fWire_scaling         = pset.get<double>("WireScaling",         1.0/64.0  ); 
  fGoodHitAveragingEps 	= pset.get<double>("GotHitAveragingEps",  2.0       );

  fDBSCANEpsilon        = pset.get<float >("DBSCANEpsilon",       1.0/16.0  );
  fDBSCANMinHits        = pset.get<int   >("DBSCANMinHits",       1         );

  fPrintDisambiguation = false;
  
  fTrack_Type = 9999;
  
  return;
}
//----------------------------------------------------------------------------------
int WCHitFinderAlg::getTrackType(std::vector<std::vector<WCHitList> > & good_hits) //HIT but need to pull parts of shouldSkipTrigger
{
  int fTrack_Type = -1;
  //Check to see if there is only a single hit on one of the WCAxes
  bool lonely_hit_bool = false;   //Single hit on at most 7 WCAxes, at least 1
  bool unique_hit_bool = true;   //Single hit on all WCAxes
  bool missing_hit_bool = false;  //Missing hit on 1 or more WCAxes
   for( size_t iWC = 0; iWC < 4 ; ++iWC ){
    for( size_t iAx = 0; iAx < 2 ; ++iAx ){
    //std::cout<<"IAx= "<<iAx<<" iWC= "<<iWC<<" size= "<<good_hits.at(iWC).at(iAx).hits.size()<<std::endl;
      if( good_hits.at(iWC).at(iAx).hits.size() == 0 ) missing_hit_bool = true;
      if( good_hits.at(iWC).at(iAx).hits.size() == 1 ) lonely_hit_bool = true;
      if( good_hits.at(iWC).at(iAx).hits.size() > 1 ) unique_hit_bool = false;
    }
  } 
  if( missing_hit_bool ) fTrack_Type = 0;
  else if( lonely_hit_bool && unique_hit_bool ) fTrack_Type = 1;
  else if( lonely_hit_bool && !unique_hit_bool ) fTrack_Type = 2;
  else if( !lonely_hit_bool && !unique_hit_bool ) fTrack_Type = 3;
  else{ std::cout << "Unknown track condition. Check me." << std::endl;
  }
  
  return fTrack_Type;
}

//---------------------------------------------------------------------------------
//Main function called for each trigger
void WCHitFinderAlg::createHits( std::vector<int> tdc_number_vect,
					   std::vector<float> hit_channel_vect,
					   std::vector<float> hit_time_bin_vect,
					   std::vector<std::vector<WCHitList> > & good_hits,
					   bool verbose)
{
  fVerbose = verbose;
  
  //Initialize hit/cluster time/channel buffers that are only defined for one trigger value
  std::vector<std::vector<float> > hit_time_buffer;
  std::vector<std::vector<float> > hit_wire_buffer;
  std::vector<std::vector<float> > cluster_time_buffer;
  std::vector<std::vector<float> > cluster_wire_buffer;

  //Create a set of buffers with 1st-Dim length of 8 (for 8 wire chamber axes: 1X, 1Y, 2X, 2Y, etc.)
  initializeBuffers( hit_time_buffer, hit_wire_buffer, cluster_time_buffer, cluster_wire_buffer );
  
  //Fill the buffers with the hit times and wires for this trigger value
  fillTimeAndWireBuffers( tdc_number_vect, hit_time_buffer, hit_wire_buffer, hit_time_bin_vect, hit_channel_vect );
  
  //Create a vector of clusters with DBSCAN, where each entry into cluster_time_buffer/cluster_wire_buffer is a different cluster
  createClusters(hit_time_buffer, hit_wire_buffer, cluster_time_buffer, cluster_wire_buffer );
  
  //Finding the hits to be used in momentum reconstruction
  //Note here that only one hit may be accepted from a cluster. The idea is that a particle passing through the MWPC causes noise
  //that is spatially and temporally clustered around the initial hit of the particle. 
  //In addition, if sufficiently close together, two or more good hits can be averaged  
  findGoodHits(cluster_time_buffer,cluster_wire_buffer,good_hits);
  getTrackType(good_hits);
}
//------------------------------------------------------------------------------------------
//Set the length of the time and wire buffers to hold 8 arrays (the number of the wire chamber axes)
void WCHitFinderAlg::initializeBuffers( std::vector<std::vector<float> > & hit_time_buffer,
					   std::vector<std::vector<float> > & hit_wire_buffer,
					   std::vector<std::vector<float> > & cluster_time_buffer,
					   std::vector<std::vector<float> > & cluster_wire_buffer ) 
{
  std::vector<float> temp_buffer;
  for( int iWCAx = 0; iWCAx < fNumber_wire_chambers*2 ; ++iWCAx ){
    hit_time_buffer.push_back(temp_buffer);
    hit_wire_buffer.push_back(temp_buffer);
    cluster_time_buffer.push_back(temp_buffer);
    cluster_wire_buffer.push_back(temp_buffer);
  }
  
//Sanity check
 if( fVerbose == true ){
 std::cout << "Length of the hit time buffer: " << hit_time_buffer.size() << std::endl;
}

}
//=====================================================================
//Fill buffers that are indexed by wire chamber # (0-7)
//Here, 0 is WC1X, 1 is WC1Y, 2 is WC2X, and so on...
//This standard is used throughout this file
//UPDATED
//HITS!!!!
void WCHitFinderAlg::fillTimeAndWireBuffers( const std::vector<int> & tdc_number_vect,
						std::vector<std::vector<float> > & hit_time_buffer,
						std::vector<std::vector<float> > & hit_wire_buffer,
						const std::vector<float> & hit_time_vect,
						const std::vector<float> & hit_channel_vect)
{
//Sanity check
if( fVerbose ){
 std::cout << "*************** Buffer Filling Info *****************" << std::endl;
}

  //Loop over the wire plane axes
  for( int iWCAx = 0; iWCAx < fNumber_wire_chambers*2 ; ++iWCAx ){
    //Loop over the tdc labels (possibly repeated) within this trigger and find those with the WCAxis = iWCAx
    for( size_t iTDC = 0; iTDC < tdc_number_vect.size(); ++iTDC ){

      int hit_wire_chamber_axis = int((tdc_number_vect.at(iTDC)-1)/2); //-1 is for tdc index
           // std::cout << "TEST: tdc: " << tdc_number_vect.at(iTDC)-1 << ", WCAx: " << hit_wire_chamber_axis << ", iWCAx: " << iWCAx << std::endl;
      if( hit_wire_chamber_axis == iWCAx ){
	float wire = 0;
	convertToWireNumber( hit_channel_vect.at(iTDC), tdc_number_vect.at(iTDC)-1, wire );
	hit_wire_buffer.at(iWCAx).push_back(wire);
	hit_time_buffer.at(iWCAx).push_back(float(hit_time_vect.at(iTDC)));

	//Sanity Check
	//	if( fVerbose ){
	//  std::cout << "(iTDC,Channel): (" << tdc_number_vect.at(iTDC)-1 << "," << hit_channel_vect.at(trigger_number).at(iTDC) << "), (WCAx,Wire): (" << iWCAx << "," << wire << ")" << ", time: " << hit_time_vect.at(trigger_number).at(iTDC) << std::endl;
	//	}
      }
    }
  }
}
//=====================================================================
//Take in tdc and channel and convert to wire number and axis
//UPDATED
//HITS!!!
void WCHitFinderAlg::convertToWireNumber(int channel,
			 int TDC_index,
			 float & wire)
{
  wire = ((TDC_index+1) % 2) *64 - channel;                  //Want wire number to be zero in middle
}

//=====================================================================
//Use DBSCAN to create clusters of hits that are spatially and temporally close
//UPDATED
//HITS!!!!
void WCHitFinderAlg::createClusters(    std::vector<std::vector<float> > hit_time_buffer,
					std::vector<std::vector<float> > hit_wire_buffer,
					std::vector<std::vector<float> > & cluster_time_buffer,
					std::vector<std::vector<float> > & cluster_wire_buffer)
{

  //Parameter for clusters being too large
  size_t max_hits = 10;

  //Loop through WC Axes
  for( int iWCAx = 0; iWCAx < fNumber_wire_chambers*2; ++iWCAx ){

    //Create hit and scaled hit vectors for use in DBSCAN's clustering
    WCHitList hits;
    WCHitList scaled_hits;
    createHitAndScaledHitVectors(iWCAx, hit_time_buffer, hit_wire_buffer, scaled_hits ); //(fill the above vectors with scaled/hits)
   
    std::vector<WCHitList> cluster_list;
    if( scaled_hits.hits.size() != 0 )
      run_DBSCAN(iWCAx, scaled_hits,cluster_list);
    
    //Loop through clusters and see if they are too big (pancake-like cross-talk)
    for( size_t iClust = 0; iClust < cluster_list.size(); ++iClust ){
      if( cluster_list.at(iClust).hits.size() > max_hits ){ continue; }
      float wire = 9999;
      float time = 9998;
      findLowestTimeHitInCluster( cluster_list.at(iClust), wire, time );
      cluster_time_buffer.at(iWCAx).push_back(time/fTime_bin_scaling);
      cluster_wire_buffer.at(iWCAx).push_back(wire/fWire_scaling);
    }
  }
}
//=====================================================================
//Create hit vectors for convenient use in DBSCAN's clustering
//HITS!!!!!
void WCHitFinderAlg::createHitAndScaledHitVectors( int WCAx_number,
						      const std::vector<std::vector<float> > hit_time_buffer,
						      const std::vector<std::vector<float> > hit_wire_buffer,
						      WCHitList & scaled_hit_list)
{
  //Sanity Check
  if( fVerbose ){
    if( hit_time_buffer.size() != hit_wire_buffer.size() ){ std::cout << "Error: vector size mismatch." << std::endl; }
    if( hit_time_buffer.at(WCAx_number).size() != hit_wire_buffer.at(WCAx_number).size() ){ std::cout << "Error: sub-vector size mismatch." << std::endl; }
  }

  //For each element in the hit time buffer (for each hit)
  for( size_t iHit = 0; iHit < hit_time_buffer.at(WCAx_number).size(); ++iHit ){
    //Create a scaled hit and fill it with time, wire number, and whether it has been visited
    //Also with a hit index and cluster index (-1 means noise cluster)
    WCHit scaled_hit;
    scaled_hit.time = hit_time_buffer.at(WCAx_number).at(iHit)*fTime_bin_scaling;
    scaled_hit.wire = hit_wire_buffer.at(WCAx_number).at(iHit)*fWire_scaling;
    scaled_hit.isVisited = false;
    scaled_hit.hit_index = iHit;
    scaled_hit.cluster_index = -1;
    scaled_hit_list.hits.push_back(scaled_hit);
  }
}
//=====================================================================
//function that runs dbscan algorithm on the hit wire/time set
//Input is ^
//Output is vector of same length as # of hits with labels
//HITS!!!!!
void WCHitFinderAlg::run_DBSCAN(    int WCAx_number,
				    WCHitList scaled_hits, 
				    std::vector<WCHitList> & cluster_list )
{
  //Parameters for algorithim
  float epsilon = fDBSCANEpsilon;
  size_t min_hits = fDBSCANMinHits;

  //Create matrix of neighborhoods of hits, given epsilon
  std::vector<WCHitList> neighborhood_matrix = createNeighborhoodMatrix( scaled_hits, epsilon );
  
  //Cluster counter
  int cluster_counter = 0;
  for( size_t iSH = 0; iSH < scaled_hits.hits.size(); ++iSH ){
    //If hit has been visited, ignore it
    if( scaled_hits.hits.at(iSH).isVisited == true ){ continue; }
    scaled_hits.hits.at(iSH).isVisited = true;
    WCHitList neighbor_hits = regionQuery( neighborhood_matrix, scaled_hits.hits.at(iSH), scaled_hits, epsilon );
    
    //If there aren't enough nearest neighbors, count this hit as noise for now
    if( neighbor_hits.hits.size() < min_hits ){ scaled_hits.hits.at(iSH).cluster_index = -1; }
    
    //Otherwise, create a new cluster from it
    else{ 
      cluster_counter++;
      scaled_hits.hits.at(iSH).cluster_index = cluster_counter-1;
      expandCluster( neighborhood_matrix, scaled_hits.hits.at(iSH), scaled_hits, neighbor_hits, epsilon, min_hits, cluster_counter-1 );
      
    }
  }
  
 
   //Extract cluster info from scaled_hits and put into clusters
   for( int iClust = 0; iClust < cluster_counter; ++iClust ){
     WCHitList cluster;
     cluster_list.push_back(cluster);
   }
   for( size_t iSH = 0; iSH < scaled_hits.hits.size(); ++iSH ){
     if( scaled_hits.hits.at(iSH).cluster_index == -1 ){ 
       continue;
     }   
     cluster_list.at(scaled_hits.hits.at(iSH).cluster_index).hits.push_back(scaled_hits.hits.at(iSH));
   }
}
//=====================================================================
//Creating a matrix of neighbors for use in efficient DBSCAN operation
//HITS!!!!
std::vector<WCHitList> WCHitFinderAlg::createNeighborhoodMatrix( WCHitList scaled_hits,
						 float epsilon )
{
  //Final vector
  std::vector<WCHitList> neighborhoods_vector;

  //Loop through all hits
  for( size_t iHit = 0; iHit < scaled_hits.hits.size() ; ++iHit ){
    //Create a hit list representing all hits (including itself) that are within epsilon of that hit
    WCHitList neighborhood_hits;
    for( size_t iSubHit = 0; iSubHit < scaled_hits.hits.size(); ++iSubHit ){
      float distance = pow(pow(scaled_hits.hits.at(iHit).wire-scaled_hits.hits.at(iSubHit).wire,2) + pow(scaled_hits.hits.at(iHit).time-scaled_hits.hits.at(iSubHit).time,2),0.5);
      if( distance < epsilon ){ neighborhood_hits.hits.push_back(scaled_hits.hits.at(iSubHit)); }
    }
    neighborhoods_vector.push_back(neighborhood_hits);
  }
  return neighborhoods_vector;
}
//=====================================================================
//DBSCAN function - finds hits within epsilon of the central hit
//HITS!!!!!!
WCHitList WCHitFinderAlg::regionQuery( std::vector<WCHitList> neighborhood_matrix,
					  WCHit the_hit,
					  WCHitList & scaled_hits,
					  float epsilon )
{
  //Create a final hit list to return
  WCHitList neighbor_hit_list;

  for( size_t iNH = 0; iNH < neighborhood_matrix.at(the_hit.hit_index).hits.size(); ++iNH ){
    //    if( scaled_hits.hits.at(neighborhood_matrix.at(the_hit.hit_index).hits.at(iNH).hit_index).isVisited == true ){ continue; }
    neighbor_hit_list.hits.push_back(scaled_hits.hits.at(neighborhood_matrix.at(the_hit.hit_index).hits.at(iNH).hit_index));
  }
  return neighbor_hit_list;
}	
//=====================================================================
//DBSCAN function - once a cluster is seeded, reach out from it and find
//other hits that are also in this cluster
//HITS!!!!!
void WCHitFinderAlg::expandCluster( std::vector<WCHitList> neighborhood_matrix,
		    WCHit the_hit,
		    WCHitList & scaled_hits,
		    WCHitList neighbor_hits,
		    float epsilon,
		    size_t min_hits,
		    int cluster_index)
{

   //Loop through neighbors
  for( size_t iNB = 0; iNB < neighbor_hits.hits.size() ; ++iNB ){
    if( neighbor_hits.hits.at(iNB).isVisited == false ){
      //Need to do each setting for the neighbor hit and the scaled hits cluster (scaled hits retains all info)
      neighbor_hits.hits.at(iNB).isVisited = true;
      scaled_hits.hits.at(neighbor_hits.hits.at(iNB).hit_index).isVisited = true;
      WCHitList next_neighbors = regionQuery( neighborhood_matrix, neighbor_hits.hits.at(iNB), scaled_hits, epsilon );
      //If there are enough next-neighbors for this neighbor, append the neighbor hits list
      if( next_neighbors.hits.size() >= min_hits ){
	for( size_t iNN = 0; iNN < next_neighbors.hits.size() ; ++iNN ){
	  //If the hit already exists in the neighbors vector, then continue
	  bool isAlreadyFound = false;
	  for( size_t iHit2 = 0; iHit2 < neighbor_hits.hits.size() ; ++iHit2 ){
	    if( next_neighbors.hits.at(iNN).hit_index == neighbor_hits.hits.at(iHit2).hit_index ){
	      isAlreadyFound = true;
	    }
	  }
	  if( isAlreadyFound == false ){
	    // if( scaled_hits.hits.at(next_neighbors.hits.at(iNN).hit_index).isVisited == false ){
	    neighbor_hits.hits.push_back(next_neighbors.hits.at(iNN));
	  }
	}
      }
    }
    
        //std::cout << "Scaled hits size: " << scaled_hits.hits.size() << ", searched hit index: " << neighbor_hits.hits.at(iNB).hit_index << std::endl;
    
    //If this neighbor hit is not yet part of a cluster, add it
    if( neighbor_hits.hits.at(iNB).cluster_index == -1 ){ 
      neighbor_hits.hits.at(iNB).cluster_index = cluster_index;
      scaled_hits.hits.at(neighbor_hits.hits.at(iNB).hit_index).cluster_index = cluster_index;
    }
    
  }
}  
//=====================================================================
//Finding the one hit with the first time to represent the true
//passing point of the particle
//HITS!!!!
void WCHitFinderAlg::findLowestTimeHitInCluster( WCHitList cluster,
				 float & wire,
				 float & time )
{
  float lowest_time = 9997;
  float lowest_time_index = 9996;
  for( size_t iHit = 0; iHit < cluster.hits.size() ; ++iHit ){
    if( cluster.hits.at(iHit).time < lowest_time ){
      lowest_time = cluster.hits.at(iHit).time;
      lowest_time_index = iHit;
    }
  }
  wire = cluster.hits.at(lowest_time_index).wire;
  time = cluster.hits.at(lowest_time_index).time;
}
//=====================================================================
//NOTE: BE CAREFUL ABOUT REFERENCING STRUCTS - MIGHT NOT REFERENCE STUFF INSIDE? CHECK
//HITS!!!
void WCHitFinderAlg::findGoodHits( std::vector<std::vector<float> > cluster_time_buffer,
				      std::vector<std::vector<float> > cluster_wire_buffer,
				      std::vector<std::vector<WCHitList> > & good_hits)
{
  //Loop through wire chamber axes (remember, 0-7)
  for( int iWCAx = 0; iWCAx < fNumber_wire_chambers*2 ; ++iWCAx ){
    
    //Sanity check
    if(fVerbose){
      if( cluster_time_buffer.at(iWCAx).size() != cluster_wire_buffer.at(iWCAx).size() ){
	std::cout << "Cluster wire/time buffer size mismatch! Error!" << std::endl;
	return;
      }
    }    
    size_t number_clusters = cluster_time_buffer.at(iWCAx).size();   
    

      //std::cout << "number clusters WCAx: " << iWCAx << " is " << number_clusters << std::endl;
 
    //Loop through clusters
    for( size_t iClust = 0; iClust < number_clusters; ++iClust ){
      float time = float(cluster_time_buffer.at(iWCAx).at(iClust));
      float wire = float(cluster_wire_buffer.at(iWCAx).at(iClust));
      //Convert back to wire chamber and axis for good hits. Lack of foresight on my part (as everywhere here)
      int wire_chamber = int(iWCAx/2);
      int axis = iWCAx % 2;
      finalizeGoodHits(wire,time,good_hits.at(wire_chamber).at(axis));
    }
  }
  
  if( fVerbose ){
    std::cout << "Number of good hits in each wire plane: "; 
    for( int iWC = 0; iWC < 4 ; ++iWC ){
      for (int iAx = 0; iAx < 2 ; ++iAx ){
	std::cout << good_hits.at(iWC).at(iAx).hits.size() << ", ";
      }
    }
    std::cout << std::endl;
  }
}  
//=====================================================================
//Finalize the hits and place them into the final good hit list
//HITS!!!!!
void WCHitFinderAlg::finalizeGoodHits(float wire,
					 float time,
					 WCHitList & finalGoodHitList)
{
  //Loop through the existing good hit list
  for( size_t iHit = 0; iHit < finalGoodHitList.hits.size(); ++iHit ){
    WCHit hit = finalGoodHitList.hits.at(iHit);
    //If there are good hits close enough to each other, average them and get rid of the old hit
    if( fabs(hit.wire-wire) < fGoodHitAveragingEps && fabs(hit.time-time) < fGoodHitAveragingEps ){
      float average_wire = (hit.wire+wire)/2;
      float average_time = (hit.time+time)/2;
      finalGoodHitList.hits.erase(finalGoodHitList.hits.begin()+iHit);
      std::cout<<"FINALIZING!!! Wire= "<<wire<<" changing to average wire= "<<average_wire<<" and time= "<<time<<" changing to average time "<<average_time<<std::endl; 
      wire = average_wire;
      time = average_time;
    }
  }
  //Now that we have averaged, push back the final good hit list with the (possibly) averaged hit
  WCHit finalGoodHit;
  finalGoodHit.wire = wire;
  finalGoodHit.time = time;
  finalGoodHitList.hits.push_back(finalGoodHit);
}


