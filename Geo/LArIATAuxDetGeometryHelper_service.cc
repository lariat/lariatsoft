////////////////////////////////////////////////////////////////////////////////
/// \file LArIATAuxDetGeometryHelper_service.cc
///
/// \version $Id
/// \author  rs@fnal.gov
////////////////////////////////////////////////////////////////////////////////

#include "Geometry/AuxDetChannelMapAlg.h"
#include "Geometry/AuxDetGeometryCore.h"

#include "Geo/LArIATAuxDetGeometryHelper.h"
#include "Geo/AuxDetChannelMapLArIATAlg.h"

#include <memory> // std::make_shared()


namespace lariatgeo
{

  //------------------------------------------------------------------------
  LArIATAuxDetGeometryHelper::LArIATAuxDetGeometryHelper(fhicl::ParameterSet const& pset, 
							 art::ActivityRegistry    &)
  : fPset(pset)
  , fChannelMap()
  {}

  //------------------------------------------------------------------------
  void LArIATAuxDetGeometryHelper::doConfigureAuxDetChannelMapAlg
    (fhicl::ParameterSet const& sortingParameters, geo::AuxDetGeometryCore* geom)
  {
    fChannelMap = std::make_shared<geo::AuxDetChannelMapLArIATAlg>(sortingParameters);
    if(fChannelMap) geom->ApplyChannelMap(fChannelMap);

    return;
  }
  
  //------------------------------------------------------------------------
  LArIATAuxDetGeometryHelper::AuxDetChannelMapAlgPtr_t
  LArIATAuxDetGeometryHelper::doGetAuxDetChannelMapAlg() const
  {
    return fChannelMap;
  }
  
}

DEFINE_ART_SERVICE_INTERFACE_IMPL(lariatgeo::LArIATAuxDetGeometryHelper, geo::AuxDetExptGeoHelperInterface)
