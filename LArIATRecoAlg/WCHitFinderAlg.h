////////////////////////////////////////////////////////////////////////
// This alg takes in TDC information for the LArIAT Wire Chambers     //
// and returns a vector of "good_hits" that can be used in            //
// finding good tracks later, using a track building alg of your      //
// choice.  In effect this is just the hit finding parts of the       //
// code in WCTrackBuildingSlicing_module.cc and WCHitFinderAlg.cxx //
// Author: Greg Pulliam, gkpullia@syr.edu                             //
////////////////////////////////////////////////////////////////////////



#ifndef WCHITFINDERALG_H
#define WCHITFINDERALG_H

//C++ includes
#include <vector>
#include <cmath>
#include <iostream>

//Framework includes
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"


// LArSoft includes
#include "larcore/Geometry/Geometry.h"

//LArIAT include
#include "Utilities/DatabaseUtilityT1034.h"

//ROOT includes
#include <TH1F.h>

//Structs for organizational purposes
struct WCHit
{
  float wire;
  float time;
  bool isVisited;
  int hit_index;
  int cluster_index;
};

struct WCHitList
{
  std::vector<WCHit> hits;
};
//----------------------------------
class WCHitFinderAlg{
 public:
  
//Constructor/destructor
WCHitFinderAlg( fhicl::ParameterSet const& pset );
~WCHitFinderAlg();
  
  
void reconfigure( fhicl::ParameterSet const& pset );
  
int getTrackType(std::vector<std::vector<WCHitList> > & good_hits);
  
void createHits(std::vector<int> tdc_number_vect,
		std::vector<float> hit_channel_vect,
		std::vector<float> hit_time_bin_vect,
		std::vector<std::vector<WCHitList> > & good_hits,
		bool verbose);
		
void findGoodHits( std::vector<std::vector<float> > cluster_time_buffer,
		   std::vector<std::vector<float> > cluster_wire_buffer,
		   std::vector<std::vector<WCHitList> > & good_hits);
  void finalizeGoodHits(float wire,
			float time,
			WCHitList & finalGoodHitList);

  void convertToWireNumber(int channel,
			   int TDC_index,
			   float & wire);
  
  void createClusters( std::vector<std::vector<float> > hit_time_buffer,
		       std::vector<std::vector<float> > hit_wire_buffer,
		       std::vector<std::vector<float> > & cluster_time_buffer,
		       std::vector<std::vector<float> > & cluster_wire_buffer);
  
  void findLowestTimeHitInCluster( WCHitList cluster,
				   float & wire,
				   float & time );
  
  void run_DBSCAN( int WCAx_number,
		   WCHitList scaled_hits, 
		   std::vector<WCHitList> & cluster_list );
  
  std::vector<WCHitList> createNeighborhoodMatrix( WCHitList scaled_hits,
						   float epsilon );
  
  void expandCluster( std::vector<WCHitList> neighborhood_matrix,
		      WCHit the_hit,
		      WCHitList & scaled_hits,
		      WCHitList neighbor_hits,
		      float epsilon,
		      size_t min_hits,
		      int cluster_index);
  
  WCHitList regionQuery( std::vector<WCHitList> neighborhood_matrix,
			 WCHit the_hit,
			 WCHitList & scaled_hits,
			 float epsilon );
  
  void createHitAndScaledHitVectors( int WCAx_number,
				     const std::vector<std::vector<float> > hit_time_buffer,
				     const std::vector<std::vector<float> > hit_wire_buffer,
				     WCHitList & scaled_hit_list);
  
  void fillTimeAndWireBuffers( const std::vector<int> & tdc_number_vect,
			       std::vector<std::vector<float> > & hit_time_buffer,
			       std::vector<std::vector<float> > & hit_wire_buffer,
			       const std::vector<float> & hit_time_vect,
			       const std::vector<float> & hit_channel_vect);
  
  void initializeBuffers( std::vector<std::vector<float> > & hit_time_buffer,
			  std::vector<std::vector<float> > & hit_wire_buffer,
			  std::vector<std::vector<float> > & cluster_time_buffer,
			  std::vector<std::vector<float> > & cluster_wire_buffer );

 private:
  //Hardware constants
  int fNumber_tdcs;
  int fNumber_wire_chambers;

  //Parameters that can be varied
  double fTime_bin_scaling;
  double fWire_scaling;
  double fGoodHitAveragingEps;
  float  fDBSCANEpsilon;
  int    fDBSCANMinHits;
  float  fPrintDisambiguation;
  int    fTrack_Type;

  art::ServiceHandle<util::DatabaseUtilityT1034> fDatabaseUtility;


  //Misc
  bool fVerbose;
  
};


#endif
  
