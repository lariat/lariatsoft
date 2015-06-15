////////////////////////////////////////////////////////////////////////
/// \file  ChannelMapAlg.h
/// \brief Interface to algorithm class for a specific detector channel mapping
///
/// \version $Id:  $
/// \author  brebel@fnal.gov
////////////////////////////////////////////////////////////////////////
#ifndef GEO_AUXDETCHANNELMAPALG_H
#define GEO_AUXDETCHANNELMAPALG_H

// LArSoft  libraries

// Framework libraries
#include "cetlib/exception.h"

// ROOT libraries
#include "TVector3.h"

// C/C++ standard libraries
#include <vector>
#include <map>
#include <set>


namespace geo{
  
  // forward-declaration from geometry
  class AuxDetGeometryData_t;
  class AuxDetGeo;  
  
  class AuxDetChannelMapAlg{

  public:
    
    virtual ~AuxDetChannelMapAlg() = default;

    virtual void                     Initialize(AuxDetGeometryData_t& geodata) = 0;
    virtual void                     Uninitialize() = 0;
    
    // method returns the entry in the sorted AuxDetGeo vector so that the 
    // Geometry in turn can return that object
    virtual size_t  NearestAuxDet          (const double* point, 
					    std::vector<geo::AuxDetGeo*> const& auxDets) const;
    virtual size_t  NearestSensitiveAuxDet (const double* point, 
					    std::vector<geo::AuxDetGeo*> const& auxDets) const;
    virtual size_t  ChannelToAuxDet        (std::vector<geo::AuxDetGeo*> const& auxDets,
					    std::string                  const& detName,
					    uint32_t                     const& channel) const;
    virtual std::pair<size_t, size_t>  ChannelToSensitiveAuxDet(std::vector<geo::AuxDetGeo*> const& auxDets,
								std::string                  const& detName,
								uint32_t                     const& channel) const;

 protected:

   std::map<std::string, size_t>          fADNameToGeo;             ///< map the names of the dets to the AuxDetGeo objects
   std::map<size_t, std::vector<size_t> > fADChannelToSensitiveGeo; ///< map the AuxDetGeo index to a vector of 
                                                                    ///< indices corresponding to the AuxDetSensitiveGeo index
   
 };
}
#endif // GEO_AUXDETCHANNELMAPALG_H

