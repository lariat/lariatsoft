////////////////////////////////////////////////////////////////////////////////
/// \file LArIATGeometryHelper.h
/// \brief Geometry helper service for LArIAT geometries. 
/// 
/// Handles LArIAT-specific information for the generic Geometry service
/// within LArSoft. Derived from the ExptGeoHelperInterface class
///
/// \verion $Id
/// \author brebel@fnal.gov
////////////////////////////////////////////////////////////////////////////////

#ifndef LArIAT_AUXDETExptGeoHelperInterface_h
#define LArIAT_AUXDETExptGeoHelperInterface_h

#include "larcore/Geometry/AuxDetExptGeoHelperInterface.h"

#include <memory> // std::shared_ptr<>

// Declaration
//
namespace lariatgeo
{
  class LArIATAuxDetGeometryHelper : public geo::AuxDetExptGeoHelperInterface
  {
  public:
  
    LArIATAuxDetGeometryHelper(fhicl::ParameterSet const & pset, 
			       art::ActivityRegistry &);

    /*
      Public interface for ExptGeoHelperInterface (for reference purposes)
      
      Configure, initialize and return the channel map:
      
      void ConfigureChannelMapAlg
        (fhicl::ParameterSet const& sortingParameters, geo::GeometryCore* geom);
      
      Returns null pointer if the initialization failed:
      
      ChannelMapAlgPtr_t GetChannelMapAlg() const;
    */
  
  private:
    
    virtual void doConfigureAuxDetChannelMapAlg
      (fhicl::ParameterSet const& sortingParameters, geo::AuxDetGeometryCore* geom)
      override;
    virtual AuxDetChannelMapAlgPtr_t doGetAuxDetChannelMapAlg() const override;
    
    fhicl::ParameterSet                       fPset;       ///< copy of configuration parameter set
    std::shared_ptr<geo::AuxDetChannelMapAlg> fChannelMap; ///< channel map
  
  };

}
DECLARE_ART_SERVICE_INTERFACE_IMPL(lariatgeo::LArIATAuxDetGeometryHelper, geo::AuxDetExptGeoHelperInterface, LEGACY)

#endif // LArIAT_ExptGeoHelperInterface_h
