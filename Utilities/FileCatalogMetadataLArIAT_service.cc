////////////////////////////////////////////////////////////////////////
// Name:  FileCatalogMetadataLArIAT_service.cc.  
//
// Purpose:  Implementation for FileCatalogMetadataLArIAT.
//
// Created:  28-Oct-2014,  H. Greenlee
//
////////////////////////////////////////////////////////////////////////

#include "Utilities/FileCatalogMetadataLArIAT.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/System/FileCatalogMetadata.h"

//--------------------------------------------------------------------
// Constructor.

util::FileCatalogMetadataLArIAT::
FileCatalogMetadataLArIAT(fhicl::ParameterSet const& pset, art::ActivityRegistry& reg)
{
  // Get parameters.

  fFCLName = pset.get<std::string>("FCLName");
  fFCLVersion = pset.get<std::string>("FCLVersion");
  fProjectName = pset.get<std::string>("ProjectName");
  fProjectStage = pset.get<std::string>("ProjectStage");
  fProjectVersion = pset.get<std::string>("ProjectVersion");    

  // Register for callbacks.

  reg.sPostBeginJob.watch(this, &FileCatalogMetadataLArIAT::postBeginJob);
}

//--------------------------------------------------------------------
// PostBeginJob callback.
// Insert per-job metadata via FileCatalogMetadata service.
void util::FileCatalogMetadataLArIAT::postBeginJob()
{
  // Get art metadata service.

  art::ServiceHandle<art::FileCatalogMetadata> mds;

  // Add metadata.

  mds->addMetadata("fclName", fFCLName);
  mds->addMetadata("fclVersion", fFCLVersion);
  mds->addMetadata("lariatProjectName", fProjectName);
  mds->addMetadata("lariatProjectStage", fProjectStage);
  mds->addMetadata("lariatProjectVersion", fProjectVersion);
}

DEFINE_ART_SERVICE(util::FileCatalogMetadataLArIAT)
