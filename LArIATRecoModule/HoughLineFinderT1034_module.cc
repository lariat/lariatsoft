////////////////////////////////////////////////////////////////////////
// Class:       HoughLineFinderT1034
// Module Type: producer
// File:        HoughLineFinderT1034_module.cc
//
// Generated at Thu Jun  4 11:03:49 2015 by Flor De Maria Blaszczyk using artmod
// from cetpkgsupport v1_08_05.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "canvas/Persistency/Common/FindOneP.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileDirectory.h" 
#include "art/Framework/Services/Optional/TFileService.h"
#include "CLHEP/Random/JamesRandom.h"

#include "nutools/RandomUtils/NuRandomService.h"

//C/C++ standard library
#include <algorithm>
#include <vector>
#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include <math.h>
#include <memory>

//LArSoft Includes
#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "larreco/RecoAlg/HoughBaseAlg.h"

//LArIAT Includes
#include "RawDataUtilities/TriggerDigitUtility.h"

class HoughLineFinderT1034;

namespace cluster {
  
  class HoughLineFinderT1034 : public art::EDProducer {
  public:
    explicit HoughLineFinderT1034(fhicl::ParameterSet const & p);
    
    HoughLineFinderT1034(HoughLineFinderT1034 const &) = delete;
    HoughLineFinderT1034(HoughLineFinderT1034 &&) = delete;
    HoughLineFinderT1034 & operator = (HoughLineFinderT1034 const &) = delete;
    HoughLineFinderT1034 & operator = (HoughLineFinderT1034 &&) = delete;


  private:
    void produce(art::Event & e) override;

    //Trigger module here, needed to get the triggers
    std::string fTriggerUtility;
    
    //Clustering module label
    std::string fDBScanModuleLabel;
    
    unsigned int fHoughSeed;

    //Hough transform algorithm    
    HoughBaseAlg fHLAlg;  
    
    CLHEP::HepRandomEngine& fEngine;

  };//Class HoughLineFinderT1034

}//End namespace cluster


namespace cluster {
  //****************************************************************************
  // Get parameters
  //****************************************************************************
  HoughLineFinderT1034::HoughLineFinderT1034(fhicl::ParameterSet const & pset)
    : fTriggerUtility{pset.get< std::string >("TriggerUtility")}
      //Get the clustering module label
    , fDBScanModuleLabel{pset.get< std::string >("DBScanModuleLabel")}
    , fHoughSeed{pset.get< unsigned int >("HoughSeed", 0)}
    , fHLAlg(pset.get< fhicl::ParameterSet >("HoughBaseAlg"))
    , fEngine{art::ServiceHandle<rndm::NuRandomService>()->createEngine(*this, pset, "Seed")}
  {
    produces< std::vector<recob::Cluster> >();
    produces< art::Assns<recob::Cluster, recob::Hit> >();
    produces< art::Assns<raw::Trigger, recob::Cluster > >();
  }
  
  //****************************************************************************
  // Event loop
  //****************************************************************************
  void HoughLineFinderT1034::produce(art::Event & evt)
  {
    //Get the trigger data utility (tdu)
    rdu::TriggerDigitUtility tdu(evt, fTriggerUtility);

    //Read the ClusterList object(s).
    //art::Handle< std::vector<recob::Cluster> > clusterListHandle;
    //evt.getByLabel(fDBScanModuleLabel,clusterListHandle);
    
    //Cluster vector
    std::vector< art::Ptr<recob::Cluster> > clusIn;
    /*for(unsigned int iCL = 0; iCL < clusterListHandle->size(); ++iCL) {

	art::Ptr<recob::Cluster> cluster(clusterListHandle, iCL);
        clusIn.push_back(cluster);
	}*/

    //Output collection of clusters
    auto ccol = std::make_unique<std::vector<recob::Cluster>>();
    //Association between clusters and hits
    auto assn = std::make_unique<art::Assns<recob::Cluster, recob::Hit>>();
    //Association between clusters and triggers
    auto trigclusassn = std::make_unique<art::Assns<raw::Trigger, recob::Cluster>>();

    //Vector of hits associated with the Hough Transform
    std::vector< art::PtrVector<recob::Hit> > clusHitsOut;   

    //Read clusters associated with the trigger
    art::FindManyP<recob::Cluster> ClusterDigits(tdu.EventTriggersPtr(), evt, fDBScanModuleLabel);
   
    //Loop on triggers
    for(size_t itrig = 0; itrig < tdu.NTriggers(); ++itrig) {

      size_t numclus = 0;      

      //Get trigger
      art::Ptr<raw::Trigger> trigger = tdu.EventTriggersPtr()[itrig];
      
      //Get clusters associated with this trigger
      clusIn = ClusterDigits.at(itrig); 

      // If a nonzero random number seed has been provided, overwrite the seed already initialized
      if(fHoughSeed != 0)
      {
        fEngine.setSeed(fHoughSeed,0);
      } 

      numclus = fHLAlg.FastTransform(clusIn, *ccol, clusHitsOut, fEngine, evt, fDBScanModuleLabel);

      MF_LOG_DEBUG("HoughLineClusters") << "found " << numclus << "clusters with HoughBaseAlg"; 

    
      for(size_t i = 0; i < ccol->size(); ++i) {
      
	mf::LogVerbatim("Summary") << ccol->at(i); 
      
	//Associate the hits to this cluster
	//util::CreateAssn(*this, evt, *(ccol.get()), clusHitsOut[i], *(assn.get()), i);
	util::CreateAssn(*this, evt, *ccol, clusHitsOut[i], *assn, i);

	//Associate clusters to trigger
	util::CreateAssn(*this, evt, *ccol, trigger, *trigclusassn);

      }//End loop on cluster collection


    }//End trigger loop

    mf::LogVerbatim("Summary") << std::setfill('-') << std::setw(175) << "-" << std::setfill(' ');
    mf::LogVerbatim("Summary") << "HoughLineFinder Summary:";
    
    
    //Write the cluster collection and associations to the event
    evt.put(std::move(ccol));
    evt.put(std::move(assn));
    evt.put(std::move(trigclusassn));
  }//End event loop
  

}//End namespace cluster

DEFINE_ART_MODULE(cluster::HoughLineFinderT1034)
