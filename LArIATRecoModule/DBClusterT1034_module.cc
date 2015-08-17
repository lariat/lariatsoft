////////////////////////////////////////////////////////////////////////
// Class:       DBClusterT1034
// Module Type: producer
// File:        DBClusterT1034_module.cc
//
// Generated at Thu Jun  4 10:51:54 2015 by Jonathan Asaadi using artmod
// from cetpkgsupport v1_08_05.
////////////////////////////////////////////////////////////////////////

// ### Framework Includes ###
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Persistency/Common/Ptr.h"
#include "art/Persistency/Common/PtrVector.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"


// ########################
// ### LArSoft includes ###
// ########################
#include "Geometry/Geometry.h"
#include "Geometry/CryostatGeo.h"
#include "Geometry/TPCGeo.h"
#include "Geometry/PlaneGeo.h"
#include "RecoBase/Cluster.h"
#include "RecoBase/Hit.h"
#include "Utilities/AssociationUtil.h"
#include "Filters/ChannelFilter.h"
#include "RecoAlg/DBScanAlg.h"
#include "ClusterFinder/ClusterCreator.h"
#include "RecoAlg/ClusterRecoUtil/StandardClusterParamsAlg.h"
#include "RecoAlg/ClusterParamsImportWrapper.h"


#include <fstream>
#include <cstdlib>
#include "TGeoManager.h"
#include "TH1.h"
#include <vector>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <memory>

// ################################
// ### LArIAT Specific Includes ###
// ################################
#include "RawDataUtilities/TriggerDigitUtility.h"

class TH1F;


class DBClusterT1034;

namespace cluster{

class DBClusterT1034 : public art::EDProducer {
public:
  explicit DBClusterT1034(fhicl::ParameterSet const & p);
  // The destructor generated by the compiler is fine for classes
  // without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  DBClusterT1034(DBClusterT1034 const &) = delete;
  DBClusterT1034(DBClusterT1034 &&) = delete;
  DBClusterT1034 & operator = (DBClusterT1034 const &) = delete;
  DBClusterT1034 & operator = (DBClusterT1034 &&) = delete;

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

  // Declare member data here.
  
  // ################################
  // ### Module that created hits ###
  // ################################
  std::string fhitsModuleLabel;
  
  std::string fTriggerUtility;
  
  //TH1F *fhitwidth;
  //TH1F *fhitwidth_ind_test;  
  //TH1F *fhitwidth_coll_test;  

  
  
  // ####################################################
  // ### Object that implements the DB scan algorithm ###
  // ####################################################
  DBScanAlg fDBScan; 

};

}//<---End cluster namespace



namespace cluster{
// ---------------------------------------------------------------
//			    Parameter Set
// ---------------------------------------------------------------
DBClusterT1034::DBClusterT1034(fhicl::ParameterSet const & pset)
: fDBScan(pset.get< fhicl::ParameterSet >("DBScanAlg"))
// Initialize member data here.
{
   this->reconfigure(pset);
   // ###############################################################
   // ### Calling new things that will be produced by this module ###
   // ###############################################################
   produces< std::vector<recob::Cluster> >();  
   produces< art::Assns<recob::Cluster, recob::Hit>     >();
   produces< art::Assns<raw::Trigger,   recob::Cluster> >();
   
   // TODO: Need an association between clusters and triggers

}

// -------------------------------------------------------------
//		 	    Reconfigure 
// -------------------------------------------------------------
void DBClusterT1034::reconfigure(fhicl::ParameterSet const & p)
{
   fhitsModuleLabel = p.get< std::string >("HitsModuleLabel");
   fTriggerUtility = p.get< std::string >("TriggerUtility");
   fDBScan.reconfigure(p.get< fhicl::ParameterSet >("DBScanAlg"));

}

//---------------------------------------------------------
//			   BeginJob
//---------------------------------------------------------
void DBClusterT1034::beginJob()
{
  // Implementation of optional member function here.
}


// ---------------------------------------------------------
// 			     Event
// ---------------------------------------------------------
void DBClusterT1034::produce(art::Event & evt)
{
   // ###########################################
   // ### Grab the trigger data utility (tdu) ###
   // ###########################################
   rdu::TriggerDigitUtility tdu(evt, fTriggerUtility);
   
   // #########################################################
   // ### Make a collection of clusters to put on the event ###
   // #########################################################  
   std::unique_ptr<std::vector<recob::Cluster> > ccol(new std::vector<recob::Cluster>);
   
   // ###################################################
   // ### Make associations between clusters and hits ###
   // ###################################################
   std::unique_ptr< art::Assns<recob::Cluster, recob::Hit> > assn(new art::Assns<recob::Cluster, recob::Hit>);
   
   // #######################################################
   // ### Make associations between clusters and triggers ###
   // #######################################################
   std::unique_ptr< art::Assns<raw::Trigger, recob::Cluster> > TrigCluAssn(new art::Assns<raw::Trigger, recob::Cluster>);

   // prepare the algorithm to compute the cluster characteristics;
   // we use the "standard" one here; configuration would happen here,
   // but we are using the default configuration for that algorithm
   ClusterParamsImportWrapper<StandardClusterParamsAlg> ClusterParamAlgo;
   
   // ################################
   // ### Calling Geometry Service ###
   // ################################
   art::ServiceHandle<geo::Geometry> geom;
   
    art::Handle< std::vector<recob::Hit> > hitcol;
    evt.getByLabel(fhitsModuleLabel,hitcol);

   
   // ########################################################
   // ### Reading in the hits assocatied with the triggers ###
   // ########################################################
   art::FindManyP<recob::Hit> HitDigits(tdu.EventTriggersPtr(), evt, fhitsModuleLabel);
   
   // ##############################
   // ### Vector for hit objects ###
   // ##############################
   std::vector< art::Ptr<recob::Hit> > allhits;
   
   // ##############################
   // ### Filtering Bad Channels ###
   // ##############################
   filter::ChannelFilter chanFilt;

   
   // ##################################
   // ### Loop over all the triggers ###
   // ##################################
   for(size_t trig = 0; trig < tdu.NTriggers(); trig++)
      { 
      std::cout<<"trigger number = "<<trig<<std::endl;
      
      // === Getting the pointer for this trigger ====
      art::Ptr<raw::Trigger> trigger = tdu.EventTriggersPtr()[trig];
      
      allhits.clear();
      allhits = HitDigits.at(trig);
      
      // Making a map of the geo::PlaneID to vectors of art::Ptr<recob::Hit>
      std::map<geo::PlaneID, std::vector< art::Ptr<recob::Hit> > > planeIDToHits;
      for(size_t i = 0; i < HitDigits.at(trig).size(); ++i){
	planeIDToHits[allhits[i]->WireID().planeID()].push_back(allhits[i]);
	//std::cout<<"trig = "<<trig<<"/"<< tdu.NTriggers() << " " <<allhits[i]->WireID().Plane << " " << allhits[i]->WireID().Wire << " " <<allhits[i]->PeakTime() << std::endl;
	}
      
      
      for(auto & itr : planeIDToHits)
         { std::cout << "Plane ID:	" << itr.first.Plane << std::endl;
	//for(auto &hit: itr.second)
	//{std::cout << itr.first.Plane << " " << hit->WireID().Wire << " "  << hit->PeakTime() << std::endl;}
         //geo::SigType_t sigType = geom->SignalType(itr.first);
	 allhits.resize(itr.second.size());
	 allhits.swap(itr.second);
	 
	 fDBScan.InitScan(allhits, chanFilt.SetOfBadChannels());
	 
	 // #######################################
	 // ### Looping over fDBScan.fbs.size() ###
	 // #######################################
	 std::cout<<"Loop over fps.size"<<std::endl;
	 for(unsigned int j = 0; j < fDBScan.fps.size(); ++j){
	   LOG_VERBATIM("DBClusterT1034") << "j = " << j 
					  << "\nallhits.size() = " << allhits.size()
					  << " , fDBScan.fps.size() = " << fDBScan.fps.size();
	   if(allhits.size() != fDBScan.fps.size()) break;
	   //fhitwidth->Fill(fDBScan.fps[j][1]);
	   
	   //if(sigType == geo::kInduction)  fhitwidth_ind_test->Fill(fDBScan.fps[j][1]);
	   //if(sigType == geo::kCollection) fhitwidth_coll_test->Fill(fDBScan.fps[j][1]);
	   
	   
	 }//<---End j loop
	 
	 // ######################
	 // ### Run Clustering ###
	 // ######################
	 fDBScan.run_cluster();
	
	 // #############################################
	 // ### Looping over fDBScan.fclusters.size() ###
	 // #############################################
	 for(size_t i = 0; i < fDBScan.fclusters.size(); ++i)
	    {
	     
	    art::PtrVector<recob::Hit> clusterHits;
	    double totalQ = 0.;
	    
	    // ######################################################
	    // ### Loop over fDBScan.fpointId_to_clusterId.size() ###
	    // ######################################################
	    for(size_t j = 0; j < fDBScan.fpointId_to_clusterId.size(); ++j)
	       {
	       if(fDBScan.fpointId_to_clusterId[j]==i)
	          {
		  clusterHits.push_back(allhits[j]);
		  totalQ += clusterHits.back()->Integral();
		  
		  }// End if statement
	       }//<--- end j loop
	     
	     // ######################################################  
	     // ### If there are enough hits then fill the cluster ###
	     // ######################################################
	   if (clusterHits.size()>0){
	     const geo::WireID& wireID = clusterHits.front()->WireID();
	     unsigned int sw = wireID.Wire;
	     unsigned int ew = clusterHits.back()->WireID().Wire;
	     
	     // feed the algorithm with all the cluster hits
	     
	     ClusterParamAlgo.ImportHits(clusterHits);
	     
	     // create the recob::Cluster directly in the vector  
	     
	     ClusterCreator cluster(
				    ClusterParamAlgo,                     // algo
				    float(sw),                            // start_wire
				    0.,                                   // sigma_start_wire
				    clusterHits.front()->PeakTime(),      // start_tick
				    clusterHits.front()->SigmaPeakTime(), // sigma_start_tick
				    float(ew),                            // end_wire
				    0.,                                   // sigma_end_wire,
				    clusterHits.back()->PeakTime(),       // end_tick
				    clusterHits.back()->SigmaPeakTime(),  // sigma_end_tick
				    ccol->size(),                         // ID
				    clusterHits.front()->View(),          // view
				    wireID.planeID(),                     // plane
				    recob::Cluster::Sentry                // sentry
				    );
	     
	     ccol->emplace_back(cluster.move());
	     // associate the hits to this cluster
	     util::CreateAssn(*this, evt, *ccol, clusterHits, *assn);
	     util::CreateAssn(*this, evt, *ccol, trigger, *TrigCluAssn);
	     clusterHits.clear(); 
	     
	   }//<---End filling if there are enough hits
	 }//<--- end i loop
	 
	 allhits.clear();
	 
      }//<---End itr auto loop
      
      
      if(ccol->size() == 0){mf::LogWarning("DBCluster") << "No clusters made for this trigger.";}
  

      }//<---End trig loop
   
   mf::LogVerbatim("Summary") << std::setfill('-') << std::setw(175) << "-" << std::setfill(' ');
   mf::LogVerbatim("Summary") << "DBcluster Summary:";
   for(unsigned int i = 0; i<ccol->size(); ++i) mf::LogVerbatim("Summary") << ccol->at(i) ;
   evt.put(std::move(ccol));
   evt.put(std::move(assn));
   evt.put(std::move(TrigCluAssn));
   return;
   
      
   
}// <---End Evt loop



void DBClusterT1034::beginRun(art::Run & r)
{
 
  // Implementation of optional member function here.
}


void DBClusterT1034::beginSubRun(art::SubRun & sr)
{
  // Implementation of optional member function here.
}

void DBClusterT1034::endJob()
{
  // Implementation of optional member function here.
}

void DBClusterT1034::endRun(art::Run & r)
{
  // Implementation of optional member function here.
}

void DBClusterT1034::endSubRun(art::SubRun & sr)
{
  // Implementation of optional member function here.
}



void DBClusterT1034::respondToCloseInputFile(art::FileBlock const & fb)
{
  // Implementation of optional member function here.
}

void DBClusterT1034::respondToCloseOutputFiles(art::FileBlock const & fb)
{
  // Implementation of optional member function here.
}

void DBClusterT1034::respondToOpenInputFile(art::FileBlock const & fb)
{
  // Implementation of optional member function here.
}

void DBClusterT1034::respondToOpenOutputFiles(art::FileBlock const & fb)
{
  // Implementation of optional member function here.
}




DEFINE_ART_MODULE(DBClusterT1034)

}//<---End namespace
