////////////////////////////////////////////////////////////////////////
// Class:       CalWireTPCQuality
// Module Type: filter
// File:        CalWireTPCQuality_module.cc
//
// Generated at Wed Feb 17 11:04:55 2016 by Irene Nutini using artmod
// from cetpkgsupport v1_08_06.
////////////////////////////////////////////////////////////////////////

// ##########################
// ### Framework Includes ###
// ##########################
#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "canvas/Persistency/Common/FindOneP.h" 
#include "canvas/Persistency/Common/Ptr.h" 
#include "canvas/Persistency/Common/PtrVector.h" 
#include "art/Framework/Services/Registry/ServiceHandle.h"

// ###########################
// ### LArIATsoft Includes ###
// ###########################
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "larcoreobj/SimpleTypesAndConstants/RawTypes.h" // raw::ChannelID_t
#include "larcore/Geometry/Geometry.h"	
#include "larcorealg/Geometry/CryostatGeo.h"	
#include "larcorealg/Geometry/TPCGeo.h"	
#include "larcorealg/Geometry/PlaneGeo.h"	
#include "larcorealg/Geometry/WireGeo.h"
#include "lardataobj/RecoBase/Wire.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardata/Utilities/AssociationUtil.h"
// ####################
// ### C++ Includes ###
// ####################
#include <iostream>
#include <memory>

class CalWireTPCQuality;

class CalWireTPCQuality : public art::EDFilter {
public:
  explicit CalWireTPCQuality(fhicl::ParameterSet const & p);
  // The destructor generated by the compiler is fine for classes
  // without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  CalWireTPCQuality(CalWireTPCQuality const &) = delete;
  CalWireTPCQuality(CalWireTPCQuality &&) = delete;
  CalWireTPCQuality & operator = (CalWireTPCQuality const &) = delete;
  CalWireTPCQuality & operator = (CalWireTPCQuality &&) = delete;

  // Required functions.
  bool filter(art::Event & e) override;

  // Selected optional functions.
  void beginJob() override;
  void reconfigure(fhicl::ParameterSet const & p) ;

private:

  // Declare member data here.
std::string     fCalDataModuleLabel;
double 		fnCalWireObjects;

};


CalWireTPCQuality::CalWireTPCQuality(fhicl::ParameterSet const & p)
// :
// Initialize member data here.
{
  // Call appropriate produces<>() functions here.
  this->reconfigure(p);
}

bool CalWireTPCQuality::filter(art::Event & e)
{
  // Implementation of required member function here.
  // ##########################################
   // ### Reading in the Wire List object(s) ###
   // ##########################################   
   art::Handle< std::vector<recob::Wire> > wireVeclHandle;
   std::vector< art::Ptr<recob::Wire> >wireVec;
   
   if(e.getByLabel(fCalDataModuleLabel,wireVeclHandle))
   {art::fill_ptr_vector(wireVec, wireVeclHandle);}
   
   
//std::cout << "Calwire objects expected " << fnCalWireObjects << std::endl;

// ### Reject the event if there is no good CalWire info - calwire objects != 480 ### 
//std::cout << "There are: " << wireVec.size() << " calData wire signals for this event" << std::endl;
if(wireVec.size() <  fnCalWireObjects){return false;}  

else {return true;}
   
   
 
   
   
   
   
}

void CalWireTPCQuality::beginJob()
{
  // Implementation of optional member function here.
}

void CalWireTPCQuality::reconfigure(fhicl::ParameterSet const & p)
{
  // Implementation of optional member function here.
  fCalDataModuleLabel = p.get< std::string  >("CalDataModuleLabel");
  fnCalWireObjects    = p.get< double >("nCalWireObjects", 480.);
}

DEFINE_ART_MODULE(CalWireTPCQuality)
