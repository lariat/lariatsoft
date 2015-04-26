////////////////////////////////////////////////////////////////////////////////
/// \file LArIATGeometryHelper_service.cc
///
/// \version $Id
/// \author  rs@fnal.gov
////////////////////////////////////////////////////////////////////////////////

#include "Geo/LArIATGeometryHelper.h"
#include "Geometry/ChannelMapAlg.h"

#include "Geo/ChannelMapLArIATAlg.h"

#include "TString.h"


namespace lariatgeo
{

  //------------------------------------------------------------------------
  LArIATGeometryHelper::LArIATGeometryHelper(fhicl::ParameterSet const& pset, 
					     art::ActivityRegistry    & reg )
  : fPset(pset)
  , fReg(reg)
  , fChannelMap()
  {}

  //------------------------------------------------------------------------
  LArIATGeometryHelper::~LArIATGeometryHelper() throw()
  {}  
  
  //------------------------------------------------------------------------
  void LArIATGeometryHelper::doConfigureChannelMapAlg(const TString                  & detectorName,
						      fhicl::ParameterSet      const & sortingParam,
						      std::vector<geo::CryostatGeo*> & c, 
						      std::vector<geo::AuxDetGeo*>   & ad  )
  {
    fChannelMap = nullptr;
    
    fChannelMap = std::shared_ptr<geo::ChannelMapAlg>(new geo::ChannelMapLArIATAlg(sortingParam));
    if(fChannelMap) fChannelMap->Initialize(c,ad);

    return;
  }
  
  //------------------------------------------------------------------------
  std::shared_ptr<const geo::ChannelMapAlg> LArIATGeometryHelper::doGetChannelMapAlg() const
  {
    return fChannelMap;
  }

}

DEFINE_ART_SERVICE_INTERFACE_IMPL(lariatgeo::LArIATGeometryHelper, geo::ExptGeoHelperInterface)
