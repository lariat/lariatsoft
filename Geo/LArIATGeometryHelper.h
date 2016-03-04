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

#ifndef LArIAT_ExptGeoHelperInterface_h
#define LArIAT_ExptGeoHelperInterface_h

#include "larcore/Geometry/ExptGeoHelperInterface.h"

#include <memory> // std::shared_ptr<>

// Declaration
//
namespace lariatgeo
{
  class LArIATGeometryHelper : public geo::ExptGeoHelperInterface
  {
  public:
  
    LArIATGeometryHelper(fhicl::ParameterSet const & pset, 
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
    
    virtual void doConfigureChannelMapAlg
      (fhicl::ParameterSet const& sortingParameters, geo::GeometryCore* geom)
      override;
    virtual ChannelMapAlgPtr_t doGetChannelMapAlg() const override;
    
    fhicl::ParameterSet                 fPset;       ///< copy of configuration parameter set
    std::shared_ptr<geo::ChannelMapAlg> fChannelMap; ///< channel map
  
  };

}
DECLARE_ART_SERVICE_INTERFACE_IMPL(lariatgeo::LArIATGeometryHelper, geo::ExptGeoHelperInterface, LEGACY)

#endif // LArIAT_ExptGeoHelperInterface_h
