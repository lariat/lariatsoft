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
#include "art/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Framework/Core/FindOneP.h"
#include "art/Framework/Core/FindManyP.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileDirectory.h" 
#include "art/Framework/Services/Optional/TFileService.h"
#include "CLHEP/Random/JamesRandom.h"

#include "artextensions/SeedService/SeedService.hh"

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
#include "lardata/RawData/RawDigit.h"
#include "lardata/RecoBase/Cluster.h"
#include "lardata/RecoBase/Hit.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "larreco/RecoAlg/HoughBaseAlg.h"

// ROOT Includes


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

    // Required functions.
    void produce(art::Event & e) override;

    // Selected optional functions.
    void beginJob() override;
    void beginRun(art::Run & r) override;
    void beginSubRun(art::SubRun & sr) override;
    void endJob() override;
    void endRun(art::Run & r) override;
    void endSubRun(art::SubRun & sr) override;
    void reconfigure(fhicl::ParameterSet const & p) override;
    void respondToCloseInputFile(art::FileBlock const & fb) override;
    void respondToCloseOutputFiles(art::FileBlock const & fb) override;
    void respondToOpenInputFile(art::FileBlock const & fb) override;
    void respondToOpenOutputFiles(art::FileBlock const & fb) override;

  private:

    //Trigger module here, needed to get the triggers
    std::string fTriggerUtility;
    
    //Clustering module label
    std::string fDBScanModuleLabel;
    
    unsigned int fHoughSeed;

    //Hough transform algorithm    
    HoughBaseAlg fHLAlg;  
    

  };//Class HoughLineFinderT1034

}//End namespace cluster


namespace cluster {
  //****************************************************************************
  // Get parameters
  //****************************************************************************
  HoughLineFinderT1034::HoughLineFinderT1034(fhicl::ParameterSet const & pset)
  : fHLAlg(pset.get< fhicl::ParameterSet >("HoughBaseAlg"))
  {
    this->reconfigure(pset);
    
    //What this module produces
    produces< std::vector<recob::Cluster> >();
    produces< art::Assns<recob::Cluster, recob::Hit> >();
    produces< art::Assns<raw::Trigger, recob::Cluster > >();
    
    //Create random number engine
    art::ServiceHandle<artext::SeedService>()->createEngine(*this, pset, "Seed");
  }
  
  //****************************************************************************
  // Reconfigure
  //****************************************************************************
  void HoughLineFinderT1034::reconfigure(fhicl::ParameterSet const & p)
  {
    //Get the clustering module label
    fDBScanModuleLabel = p.get< std::string >("DBScanModuleLabel");
    
    //Get trigger utility
    fTriggerUtility = p.get< std::string >("TriggerUtility");
    
    //Configure Hough 
    fHoughSeed = p.get< unsigned int >("HoughSeed", 0);
    fHLAlg.reconfigure(p.get< fhicl::ParameterSet >("HoughBaseAlg"));
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
    std::unique_ptr< std::vector<recob::Cluster> > ccol(new std::vector<recob::Cluster>);
    //Association between clusters and hits
    std::unique_ptr< art::Assns<recob::Cluster, recob::Hit> > assn(new art::Assns<recob::Cluster, recob::Hit>);
    //Association between clusters and triggers
    std::unique_ptr< art::Assns<raw::Trigger, recob::Cluster> > trigclusassn(new art::Assns< raw::Trigger, recob::Cluster > );

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
	art::ServiceHandle<art::RandomNumberGenerator> rng;
	CLHEP::HepRandomEngine &engine = rng->getEngine();
	engine.setSeed(fHoughSeed,0);
      } 

      numclus = fHLAlg.FastTransform(clusIn, *ccol, clusHitsOut, evt, fDBScanModuleLabel);

      LOG_DEBUG("HoughLineClusters") << "found " << numclus << "clusters with HoughBaseAlg"; 

    
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
    
    //return;    
    
  }//End event loop
  

  //****************************************************************************
  //Other methods
  //****************************************************************************
  void HoughLineFinderT1034::beginJob()
  {
    // Implementation of optional member function here.
  }
  
  void HoughLineFinderT1034::beginRun(art::Run & r)
  {
    // Implementation of optional member function here.
  }
  
  void HoughLineFinderT1034::beginSubRun(art::SubRun & sr)
  {
    // Implementation of optional member function here.
  }
  
  void HoughLineFinderT1034::endJob()
  {
    // Implementation of optional member function here.
  }
  
  void HoughLineFinderT1034::endRun(art::Run & r)
  {
    // Implementation of optional member function here.
  }
  
  void HoughLineFinderT1034::endSubRun(art::SubRun & sr)
  {
    // Implementation of optional member function here.
  }
  
  void HoughLineFinderT1034::respondToCloseInputFile(art::FileBlock const & fb)
  {
    // Implementation of optional member function here.
  }
  
  void HoughLineFinderT1034::respondToCloseOutputFiles(art::FileBlock const & fb)
  {
    // Implementation of optional member function here.
  }
  
  void HoughLineFinderT1034::respondToOpenInputFile(art::FileBlock const & fb)
  {
    // Implementation of optional member function here.
  }
  
  void HoughLineFinderT1034::respondToOpenOutputFiles(art::FileBlock const & fb)
  {
    // Implementation of optional member function here.
  }

}//End namespace cluster

namespace cluster{  
  DEFINE_ART_MODULE(HoughLineFinderT1034)
}
