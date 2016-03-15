////////////////////////////////////////////////////////////////////////
/// \file  GeoObjectSorterLArIAT.h
/// \brief Interface to algorithm class for standard sorting of geo::XXXGeo objects
///
/// \version $Id:  $
/// \author  brebel@fnal.gov
////////////////////////////////////////////////////////////////////////
#ifndef GEO_GEOOBJECTSORTERSTANDARD_H
#define GEO_GEOOBJECTSORTERSTANDARD_H

#include <vector>

#include "larcore/Geometry/GeoObjectSorter.h"

namespace geo{

  class GeoObjectSorterLArIAT : public GeoObjectSorter {

  public:

    GeoObjectSorterLArIAT(fhicl::ParameterSet const& p);
    ~GeoObjectSorterLArIAT();

    void SortAuxDets        (std::vector<geo::AuxDetGeo*>          & adgeo)    const;
    void SortAuxDetSensitive(std::vector<geo::AuxDetSensitiveGeo*> & adsgeo)   const;
    void SortCryostats      (std::vector<geo::CryostatGeo*>        & cgeo)     const;
    void SortTPCs     	    (std::vector<geo::TPCGeo*>      	   & tgeo)     const;
    void SortPlanes   	    (std::vector<geo::PlaneGeo*>    	   & pgeo,	      
		      	     geo::DriftDirection_t     	     const & driftDir) const;
    void SortWires    	    (std::vector<geo::WireGeo*>     	   & wgeo)     const;
    
  private:
    
  };

}

#endif // GEO_GEOOBJECTSORTERLARIAT_H
