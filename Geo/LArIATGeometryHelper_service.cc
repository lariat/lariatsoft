////////////////////////////////////////////////////////////////////////////////
/// \file LArIATGeometryHelper_service.cc
///
/// \version $Id
/// \author  rs@fnal.gov
////////////////////////////////////////////////////////////////////////////////

#include "Geo/LArIATGeometryHelper.h"
#include "larcore/Geometry/ChannelMapAlg.h"
#include "larcore/Geometry/GeometryCore.h"

#include "Geo/ChannelMapLArIATAlg.h"

#include <memory> // std::make_shared()


namespace lariatgeo
{

  //------------------------------------------------------------------------
  LArIATGeometryHelper::LArIATGeometryHelper(fhicl::ParameterSet const& pset, 
					     art::ActivityRegistry    &)
  : fPset(pset)
  , fChannelMap()
  {}

  //------------------------------------------------------------------------
  void LArIATGeometryHelper::doConfigureChannelMapAlg
    (fhicl::ParameterSet const& sortingParameters, geo::GeometryCore* geom)
  {
    fChannelMap = std::make_shared<geo::ChannelMapLArIATAlg>(sortingParameters);
    if(fChannelMap) geom->ApplyChannelMap(fChannelMap);

    return;
  }
  
  //------------------------------------------------------------------------
  LArIATGeometryHelper::ChannelMapAlgPtr_t
  LArIATGeometryHelper::doGetChannelMapAlg() const
  {
    return fChannelMap;
  }
  
}

DEFINE_ART_SERVICE_INTERFACE_IMPL(lariatgeo::LArIATGeometryHelper, geo::ExptGeoHelperInterface)
