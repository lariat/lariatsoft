////////////////////////////////////////////////////////////////////////
/// \file  AuxDetGeoObjectSorterLArIAT.cxx
/// \brief Interface to algorithm class for sorting standard geo::XXXGeo objects
///
/// \version $Id:  $
/// \author  brebel@fnal.gov
////////////////////////////////////////////////////////////////////////

#include "Geo/AuxDetGeoObjectSorterLArIAT.h"
#include "Geometry/AuxDetGeo.h"
#include "Geometry/AuxDetSensitiveGeo.h"

namespace geo{

  //----------------------------------------------------------------------------
  // Define sort order for AuxDets in standard configuration
  static bool sortAuxDetLArIAT(const AuxDetGeo* ad1, const AuxDetGeo* ad2)
  {

    // sort using the center of the detector - primary ordering by z,
    // then y, and x
    double c1[3] = {0.};
    double c2[3] = {0.};
    ad1->GetCenter(c1);
    ad2->GetCenter(c2);

    for(int i = 2; i > 0; --i){ 
      if(c1[i] != c2[i]){
	return c1[i] < c2[i];
      }
    }
    
    return c1[0] < c2[0];
   
  }

  //----------------------------------------------------------------------------
  // Define sort order for AuxDetSensitives in standard configuration
  static bool sortAuxDetSensitiveLArIAT(const AuxDetSensitiveGeo* ad1, const AuxDetSensitiveGeo* ad2)
  {

    // sort using the center of the detector - primary ordering by z,
    // then y, and x
    double c1[3] = {0.};
    double c2[3] = {0.};
    ad1->GetCenter(c1);
    ad2->GetCenter(c2);

    for(int i = 2; i > 0; --i){ 
      if(c1[i] != c2[i]){
	return c1[i] < c2[i];
      }
    }
    
    return c1[0] < c2[0];
   
  }

  //----------------------------------------------------------------------------
  AuxDetGeoObjectSorterLArIAT::AuxDetGeoObjectSorterLArIAT(fhicl::ParameterSet const&)
  {
  }

  //----------------------------------------------------------------------------
  AuxDetGeoObjectSorterLArIAT::~AuxDetGeoObjectSorterLArIAT()
  {
  }

  //----------------------------------------------------------------------------
  void AuxDetGeoObjectSorterLArIAT::SortAuxDets(std::vector<geo::AuxDetGeo*> & adgeo) const
  {
    std::sort(adgeo.begin(), adgeo.end(), sortAuxDetLArIAT);
    
    return;
  }

  //----------------------------------------------------------------------------
  void AuxDetGeoObjectSorterLArIAT::SortAuxDetSensitive(std::vector<geo::AuxDetSensitiveGeo*> & adsgeo) const
  {
    std::sort(adsgeo.begin(), adsgeo.end(), sortAuxDetSensitiveLArIAT);
    
    return;
  }

}
