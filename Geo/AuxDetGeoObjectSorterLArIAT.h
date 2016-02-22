////////////////////////////////////////////////////////////////////////
/// \file  GeoObjectSorterLArIAT.h
/// \brief Interface to algorithm class for standard sorting of geo::XXXGeo objects
///
/// \version $Id:  $
/// \author  brebel@fnal.gov
////////////////////////////////////////////////////////////////////////
#ifndef GEO_AUXDETGEOOBJECTSORTERLARIAT_H
#define GEO_AUXDETGEOOBJECTSORTERLARIAT_H

#include <vector>

#include "larcore/Geometry/AuxDetGeoObjectSorter.h"

namespace geo{

  class AuxDetGeoObjectSorterLArIAT : public AuxDetGeoObjectSorter {

  public:

    AuxDetGeoObjectSorterLArIAT(fhicl::ParameterSet const& p);
    ~AuxDetGeoObjectSorterLArIAT();

    void SortAuxDets        (std::vector<geo::AuxDetGeo*>          & adgeo)    const;
    void SortAuxDetSensitive(std::vector<geo::AuxDetSensitiveGeo*> & adsgeo)   const;
    
  private:
    
  };

}

#endif // GEO_AUXDETGEOOBJECTSORTERLARIAT_H
