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
    fADNameToGeo.clear();
    fADChannelToSensitiveGeo.clear();

    for(size_t a = 0; a < adgeo.size(); ++a){
      std::string volName(adgeo[a]->TotalVolume()->GetName());

      if(volName.find("TOFUS") != std::string::npos){
	fADNameToGeo["TOFUS"] = a;
	for(size_t c = 0; c < 2; ++c) fADChannelToSensitiveGeo[a].push_back(0);
      }
      else if(volName.find("TOFDS") != std::string::npos){
	fADNameToGeo["TOFDS"] = a;
	for(size_t c = 0; c < 2; ++c) fADChannelToSensitiveGeo[a].push_back(0);
      }
      else if(volName.find("MWPC1") != std::string::npos){
	fADNameToGeo["MWPC1"]= a;
	for(size_t c = 0; c < 128; ++c) fADChannelToSensitiveGeo[a].push_back(0);
      }
      else if(volName.find("MWPC2") != std::string::npos){
	fADNameToGeo["MWPC2"]= a;
	for(size_t c = 0; c < 128; ++c) fADChannelToSensitiveGeo[a].push_back(0);
      }
      else if(volName.find("MWPC3") != std::string::npos){
	fADNameToGeo["MWPC3"]= a;
	for(size_t c = 0; c < 128; ++c) fADChannelToSensitiveGeo[a].push_back(0);
      }
      else if(volName.find("MWPC4") != std::string::npos){
	fADNameToGeo["MWPC4"]= a;
	for(size_t c = 0; c < 128; ++c) fADChannelToSensitiveGeo[a].push_back(0);
      }
      else if(volName.find("AeroGelUS") != std::string::npos){
	fADNameToGeo["AeroGelUS"] = a;
	for(size_t c = 0; c < 2; ++c) fADChannelToSensitiveGeo[a].push_back(0);
      }	
      else if(volName.find("AeroGelDS") != std::string::npos){
	fADNameToGeo["AeroGelDS"] = a;
	for(size_t c = 0; c < 2; ++c) fADChannelToSensitiveGeo[a].push_back(0);
      }	
      else if(volName.find("MuonRangeStack") != std::string::npos){
	fADNameToGeo["MuonRangeStack"] = a;

	// check that there are the expected 16 sensitive volumes
	if(adgeo[a]->NSensitiveVolume() != 16)
	  throw cet::exception("AuxDetChannelMapLArIATAlg") << "Expected 16 sensitive volumes for MuRS AuxDet, "
						      << volName << " instead have only " 
						      << adgeo[a]->NSensitiveVolume();

	// the sorting should have sorted them by plane going from bottom to top in a plane
	// and that is what we want.
	for(size_t sv = 0; sv < 16; ++sv) fADChannelToSensitiveGeo[a].push_back(sv);	
      }
      else if(volName.find("Halo") != std::string::npos){
	fADNameToGeo["Halo"] = a;
	for(size_t c = 0; c < 2; ++c) fADChannelToSensitiveGeo[a].push_back(0);
      }

    } // loop over the AuxDetGeo objects

    return;
  }
   
  //----------------------------------------------------------------------------
  void AuxDetChannelMapLArIATAlg::Uninitialize()
  {
  }


} // namespace
