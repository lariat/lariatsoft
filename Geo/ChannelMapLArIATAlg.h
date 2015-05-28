////////////////////////////////////////////////////////////////////////
/// \file  ChannelMapLArIATAlg.h
/// \brief Interface to algorithm class for a specific detector channel mapping
///
/// \version $Id:  $
/// \author  brebel@fnal.gov
////////////////////////////////////////////////////////////////////////
#ifndef GEO_CHANNELMAPLARIATALG_H
#define GEO_CHANNELMAPLARIATALG_H

#include <vector>
#include <set>
#include <iostream>

#include "SimpleTypesAndConstants/geo_types.h"
#include "Geometry/ChannelMapAlg.h"
#include "Geo/GeoObjectSorterLArIAT.h"
#include "fhiclcpp/ParameterSet.h"

#include "TVector3.h"

namespace geo{

  class ChannelMapLArIATAlg : public ChannelMapAlg{

  public:

    ChannelMapLArIATAlg(fhicl::ParameterSet const& p);
    ~ChannelMapLArIATAlg();
    
    void                     Initialize( std::vector<geo::CryostatGeo*> & cgeo, 
					 std::vector<geo::AuxDetGeo*>   & adgeo );
    void                     Uninitialize();
    std::vector<WireID>      ChannelToWire(uint32_t channel)           const;
    uint32_t                 Nchannels()                               const;

    double                   WireCoordinate(double YPos, double ZPos,
					    unsigned int    PlaneNo,
					    unsigned int    TPCNo,
					    unsigned int    cstat)     const;

    WireID                   NearestWireID(const TVector3& worldPos,
					   unsigned int    PlaneNo,
					   unsigned int    TPCNo,
					   unsigned int    cstat)      const;
    uint32_t                 PlaneWireToChannel(unsigned int plane,
						unsigned int wire,
						unsigned int tpc,
						unsigned int cstat)    const;
   View_t                    View( uint32_t const channel )            const;
   SigType_t                 SignalType( uint32_t const channel )      const;
   std::set<View_t>  const&  Views()                                   const;
   std::set<PlaneID> const&  PlaneIDs()                                const;

   // methods for the auxiliary detectors

   // method returns the entry in the sorted AuxDetGeo vector so that the 
   // Geometry in turn can return that object
   size_t                   NearestAuxDet         (TVector3 const& point)   const;
   size_t                   NearestSensitiveAuxDet(TVector3 const& point)   const;
   size_t                   NSensitiveAuxDet(std::string const& auxDetName) const;

  private:
    
    unsigned int                                         fNcryostat;      ///< number of cryostats in the detector
    uint32_t                                             fNchannels;      ///< number of channels in the detector
    uint32_t                                             fTopChannel;     ///< book keeping highest channel #
    std::vector<unsigned int>                            fNTPC;           ///< number of TPCs in each cryostat
    std::set<View_t>                                     fViews;          ///< vector of the views present in the detector
    std::set<PlaneID>                                    fPlaneIDs;       ///< vector of the PlaneIDs present in the detector
    std::vector<std::vector<std::vector<float>>>         fFirstWireProj;  ///< Distance (0,0,0) to first wire 	 
                                                                          ///< along orth vector per plane per TPC
    std::vector<std::vector<std::vector<float>>>         fOrthVectorsY;   ///< Unit vectors orthogonal to wires in
    std::vector<std::vector<std::vector<float>>>         fOrthVectorsZ;   ///< each plane - stored as 2 components
                                                                          ///< to avoid having to invoke any bulky
                                                                          ///< TObjects / CLHEP vectors etc	 
    std::vector<std::vector<std::vector<float>>>         fWireCounts;     ///< Number of wires in each plane - for
                                                                          ///< range checking after calculation   
    std::vector<std::vector<unsigned int>>  		 fNPlanes;        ///< Number of planes in each TPC - for
                                                                          ///< range checking after calculation   
    std::vector<std::vector<std::vector<unsigned int>>>  fPlaneBaselines; ///< The number of wires in all the 
                                                                          ///< tpcs and planes up to this one 
                                                                          ///< in the heirachy
    std::vector<std::vector<std::vector<unsigned int>>>  fWiresPerPlane;  ///< The number of wires in this plane 
                                                                          ///< in the heirachy
    geo::GeoObjectSorterLArIAT                           fSorter;         ///< class to sort geo objects
  };


}
#endif // GEO_CHANNELMAPLARIATALG_H

