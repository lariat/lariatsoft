////////////////////////////////////////////////////////////////////////
/// \file  ChannelMapLArIATAlg.cxx
/// \brief Interface to algorithm class for the standar, simplest detector channel mapping
///
/// \version $Id:  $
/// \author  brebel@fnal.gov
////////////////////////////////////////////////////////////////////////

// ART Includes
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "messagefacility/MessageLogger/MessageLogger.h" 

// LArSoft Includes
#include "Geometry/AuxDetGeo.h"
#include "Geometry/AuxDetSensitiveGeo.h"

// LArIATSoft
#include "Geo/AuxDetGeometryCore.h"
#include "Geo/AuxDetChannelMapLArIATAlg.h"

// C++ standard library
#include <ostream> // std::endl


namespace geo{

  //----------------------------------------------------------------------------
  AuxDetChannelMapLArIATAlg::AuxDetChannelMapLArIATAlg(fhicl::ParameterSet const& p)
    : fSorter(geo::AuxDetGeoObjectSorterLArIAT(p))
  {
  }

  //----------------------------------------------------------------------------
  void AuxDetChannelMapLArIATAlg::Initialize( AuxDetGeometryData_t& geodata )
  {
    // start over:
    Uninitialize();
    
    std::vector<geo::AuxDetGeo*>  & adgeo = geodata.auxDets;
    
    // sort the AuxDetGeo objects and map them to names of the detectors
    // Each raw::AuxDetDigit knows the name of the detector it came from and its
    // channel - sometimes channel maps to sensitive volume and some times it doesn't
    fSorter.SortAuxDets(adgeo);    
    //\todo: Uncomment the following line with larsoft v04_09_01
    //for(auto a : adgeo) a->SortSubVolumes(fSorter);

    // map the AuxDetGeo names to their position in the sorted vector
    // Each TOF  detector has 2   channels, and 1  sensitive volume
    // Each AG   detector has 2   channels, and 1  sensitive volume
    // The  MuRS detector has 16  channels, and 16 sensitive volumes
    // Each MWPC detector has 128 channels, and 1  sensitive volume
    // The  Halo detector has 2   channels, and    sensitive volume
    fADNameToGeo.clear();
    fADChannelToSensitiveGeo.clear();

    for(size_t a = 0; a < adgeo.size(); ++a){
      std::string volName(adgeo[a]->TotalVolume()->GetName());

      fADGeoToName[a] = volName;
      fNameToADGeo[volName] = a;

      if(volName.find("TOFUS") != std::string::npos){
	for(size_t c = 0; c < 2; ++c) fADGeoToChannelAndSV[a].push_back(std::make_pair(c, 0));
      }
      else if(volName.find("TOFDS") != std::string::npos){
	for(size_t c = 0; c < 2; ++c) fADGeoToChannelAndSV[a].push_back(std::make_pair(c, 0));
      }
      else if(volName.find("MWPC1") != std::string::npos){
	for(size_t c = 0; c < 128; ++c) fADGeoToChannelAndSV[a].push_back(std::make_pair(c, 0));
      }
      else if(volName.find("MWPC2") != std::string::npos){
	for(size_t c = 0; c < 128; ++c) fADGeoToChannelAndSV[a].push_back(std::make_pair(c, 0));
      }
      else if(volName.find("MWPC3") != std::string::npos){
	for(size_t c = 0; c < 128; ++c) fADGeoToChannelAndSV[a].push_back(std::make_pair(c, 0));
      }
      else if(volName.find("MWPC4") != std::string::npos){
	for(size_t c = 0; c < 128; ++c) fADGeoToChannelAndSV[a].push_back(std::make_pair(c, 0));
      }
      else if(volName.find("AeroGelUS") != std::string::npos){
	for(size_t c = 0; c < 2; ++c) fADGeoToChannelAndSV[a].push_back(std::make_pair(c, 0));
      }	
      else if(volName.find("AeroGelDS") != std::string::npos){
	for(size_t c = 0; c < 2; ++c) fADGeoToChanAndSV[a].push_back(std::make_pair(c, 0));
      }	
      else if(volName.find("MuonRangeStack") != std::string::npos){

	// check that there are the expected 16 sensitive volumes
	if(adgeo[a]->NSensitiveVolume() != 16)
	  throw cet::exception("AuxDetChannelMapLArIATAlg") << "Expected 16 sensitive volumes for MuRS AuxDet, "
							    << volName << " instead have only " 
							    << adgeo[a]->NSensitiveVolume();

	// the sorting should have sorted them by plane going from bottom to top in a plane
	// and that is what we want.
	for(size_t sv = 0; sv < 16; ++sv) fADGeoToChannelAndSV[a].push_back(std::make_pair(sv, sv));	
      }
      else if(volName.find("Halo") != std::string::npos){
	for(size_t c = 0; c < 2; ++c) fADGeoToChannelAndSV[a].push_back(std::make_pair(c, 0));
      }

    } // loop over the AuxDetGeo objects

    return;
  }
   
  //----------------------------------------------------------------------------
  void AuxDetChannelMapLArIATAlg::Uninitialize()
  {
  }

  //----------------------------------------------------------------------------
  uint32_t AuxDetChannelMapLArIATAlg::PositionToAuxDetChannel(double                       const  worldLoc[3],
							      std::vector<geo::AuxDetGeo*> const& auxDets,
							      size_t                            & ad,	    
							      size_t                       	& sv) const
  {
    // set the default to be that we don't find the position 
    // in any of the AuxDets
    uint32_t channel = UINT_MAX;

    // figure out which detector we are in
    size_t ad = 0;
    size_t sv = this->NearestSensitiveAuxDet(worldLoc, auxDets, ad);

    // get the origin of the sensitive volume in the world coordinate system
    double svOrigin[3]    = {0.};
    double localOrigin[3] = {0.};

    auxDets[ad]->SensitiveVolume(sv).LocalToWorld(localOrigin, svOrigin);

    // check to see which AuxDet this position corresponds to
    auto gnItr = fADGeoToName.find(ad);
    if( gnItr != fADGeoToName.end() ){
      
      // get the vector of channel and sensitive volume pairs
      auto csvItr = fADGeoToChannelAndSV.find(ad);
      if( csvItr == fADGeoToChannelAndSV.end() )
	throw cet::exception("AuxDetChannelMapLArIATAlg") << "No entry in channel and sensitive volume"
							  << " map for AuxDet index " << ad << " bail";

      // there are two cases here - the MWPC and everything else
      // the MWPC wires are XXX mm apart
      if( gnItr->second.find("MWPC") != std::string::npos){
	
	// determine the number of wire spacings from the center 
	// of the MWPC to the position in both X and Y directions
	double deltaPos[3] = {worldLoc[0] - svOrigin[0],
			      worldLoc[1] - svOrigin[1],
			      worldLoc[2] - svOrigin[2]};

	double numChanX = 64 + deltaPos[0]/0.1;
	double numChanY = 64 + deltaPos[1]/0.1;

	if(numChanX < 0. || numChanY < 0.) 
	  throw cet::exception("AuxDetChannelMapLArIATAlg") << "position (" 
							    << worldLoc[0] << "," << worldLoc[1] << "," << worldLoc[2] 
							    << ") corresponds to a negative channel number: "
							    << numChanX << " " << numChanY;

	channel = (uint32_t)numChanX;
      } // end if this is an MWPC
      else if(gnItr->second.find("MuonRangeStack") != std::string::npos){
	// each SV in the MuRS has a single channel
	channel = sv;
      }
      else{
	// the remaining aux dets have 2 channels and a single SV.  Figure out which
	// channel the position is closer to
      }

    }// end finding of the detector name

    if(channel == UINT_MAX) 
      throw cet::exception("AuxDetChannelMapLArIATAlg") << "position (" 
							<< worldLoc[0] << "," << worldLoc[1] << "," << worldLoc[2]
							<< ") does not correspond to any AuxDet";

    return channel;
  }

  //----------------------------------------------------------------------------
  const TVector3 AuxDetChannelMapLArIATAlg::AuxDetChannelToPosition(uint32_t                     const& channel,
								    std::string                  const& auxDetName,
								    std::vector<geo::AuxDetGeo*> const& auxDets) const
  {
    double x = 0.;
    double y = 0.;
    double z = 0.;

    // figure out which detector we are in
    size_t ad = UINT_MAX;
    if( fNameToADGeo.count(auxDetName) > 0 )
      ad = fNameToADGeo.find(auxDetName)->second;
    else
      throw cet::exception("AuxDetChannelMapLArIATAlg") << "No AuxDetGeo with name " << auxDetName;

    // get the vector of channel and sensitive volume pairs
    auto csvItr = fADGeoToChannelAndSV.find(ad);
    if( csvItr == fADGeoToChannelAndSV.end() )
      throw cet::exception("AuxDetChannelMapLArIATAlg") << "No entry in channel and sensitive volume"
							<< " map for AuxDet index " << ad << " bail";

    // loop over the vector of channel and sensitive volumes to determine the sensitive volume 
    // for this channel
    // then get the origin of the sensitive volume in the world coordinate system
    double svOrigin[3]    = {0.};
    double localOrigin[3] = {0.};
    for(auto csv : csvItr->second){
      
      if( csv.first == channel ){

	// get the center of the sensitive volume for this channel
	auxDets[ad]->SensitiveVolume(sv).LocalToWorld(localOrigin, svOrigin);

	// there are two cases here - the MWPC and everything else
	// the MWPC wires are XXX mm apart
	if( auxDetName.find("MWPC") != std::string::npos){
	
	} // end if this is an MWPC
	else if(auxDetName.find("MuonRangeStack") != std::string::npos){
	  // each SV in the MuRS has a single channel
	}
	else{
	  // the remaining aux dets have 2 channels and a single SV.  Figure out which
	  // channel the position is closer to
	}
	
	break;
      } // end if using the correct channel
    } // end loop over vector of channels and sensitive volumes

    return TVector3(x, y, z);
  }

} // namespace
