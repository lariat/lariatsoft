////////////////////////////////////////////////////////////////////////
// Class:       ClusterCrawlerLariat
// Module Type: producer
// File:        ClusterCrawlerLariat_module.cc
//
// Generated at Fri Jun  7 09:44:09 2013 by Bruce Baller using artmod 
// from cetpkgsupport v1_02_00.
////////////////////////////////////////////////////////////////////////

// ####################
// ### C++ Includes ###
// ####################
#include <string>
#include <vector>
#include <utility> // std::pair<>
#include <memory> // std::unique_ptr<>
#include <iomanip>
#include <algorithm> // std::max()
#include <functional> // std::mem_fn()
#include <limits> // std::numeric_limits<>

// ##########################
// ### Framework includes ###
// ##########################
#include "fhiclcpp/ParameterSet.h" 
#include "messagefacility/MessageLogger/MessageLogger.h" 
#include "art/Framework/Core/ModuleMacros.h" 
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h" 
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Principal/Handle.h" 
#include "art/Persistency/Common/Ptr.h" 
#include "art/Persistency/Common/PtrVector.h" 
#include "art/Framework/Services/Registry/ServiceHandle.h" 
#include "art/Framework/Services/Optional/TFileService.h" 
#include "art/Framework/Core/FindOneP.h"
#include "art/Utilities/Exception.h"
#include "art/Utilities/InputTag.h"

// ########################
// ### LArSoft Includes ###
// ########################
#include "SimpleTypesAndConstants/geo_types.h"
#include "SimpleTypesAndConstants/RawTypes.h" // raw::ChannelID_t
#include "Geometry/Geometry.h"
#include "Geometry/CryostatGeo.h"
#include "Geometry/TPCGeo.h"
#include "Geometry/PlaneGeo.h"
#include "Filters/ChannelFilter.h"
#include "RawDataUtilities/TriggerDigitUtility.h" 
#include "RawData/RawDigit.h"
#include "RawData/raw.h"
#include "RecoBase/Wire.h"
#include "RecoBase/Hit.h"
#include "RecoBase/Cluster.h"
#include "RecoBase/Vertex.h"
#include "RecoAlg/CCHitFinderAlg.h"
#include "RecoAlg/ClusterCrawlerAlg.h"
#include "Utilities/AssociationUtil.h"
#include "TComplex.h"
#include "TFile.h"
#include "TH2D.h"
#include "TF1.h"

// #include "RecoBase/EndPoint2D.h"
// #include "RecoBaseArt/HitCreator.h" // recob::HitCollectionAssociator
// #include "RecoAlg/ClusterRecoUtil/StandardClusterParamsAlg.h"
// #include "RecoAlg/ClusterParamsImportWrapper.h"


namespace cluster {
  class ClusterCrawlerLariat;
}

class cluster::ClusterCrawlerLariat : public art::EDProducer {

  public:
    explicit ClusterCrawlerLariat(fhicl::ParameterSet const & pset);
    virtual ~ClusterCrawlerLariat();

    void reconfigure(fhicl::ParameterSet const & pset) override;
    void produce(art::Event & evt) override;
    void beginJob();

  private:
    hit::CCHitFinderAlg    fCCHFAlg;            ///< define CCHitFinderAlg object
    cluster::ClusterCrawlerAlg fCCAlg;              ///< define ClusterCrawlerAlg object
    std::string            fCalDataModuleLabel; ///< label of module producing input wires
    std::string            fDigitModuleLabel;   ///< module that made digits                
                                                       
    void Clustering(std::vector< art::Ptr<recob::Wire> > & wireVec, 
		  std::vector<recob::Cluster> & clus, 
		  std::vector< std::vector<recob::Hit> > & clusToHits,
		  std::map< size_t, std::vector<recob::Vertex> > & clusToVertex);

};


namespace cluster {

  ClusterCrawlerLariat::ClusterCrawlerLariat(fhicl::ParameterSet const& pset)
    : fCCHFAlg(pset.get< fhicl::ParameterSet >("CCHitFinderAlg"))
    , fCCAlg  (pset.get< fhicl::ParameterSet >("ClusterCrawlerAlg"))
  {
    mf::LogWarning("ClusterCrawlerLariat") <<
      "\nClusterCrawlerLariat module has been deprecated and will be removed."
      "\nIt is now replaced by HitFinder and LineCluster modules.";
      
    
    this->reconfigure(pset);
    
    produces< std::vector<recob::Cluster> >();  
    produces< std::vector<recob::Vertex> >();
    produces< art::Assns<recob::Wire,    recob::Hit>     >();
    produces< art::Assns<recob::Cluster, recob::Hit>     >();
    produces< art::Assns<recob::Cluster, recob::Vertex>  >();
    produces< art::Assns<raw::Trigger,   recob::Hit>     >();
    produces< art::Assns<raw::Trigger,   recob::Vertex>  >();
    produces< art::Assns<raw::Trigger,   recob::Cluster> >();
  }

  ClusterCrawlerLariat::~ClusterCrawlerLariat()
  {
  }

  void ClusterCrawlerLariat::reconfigure(fhicl::ParameterSet const & pset)
  {
    fCalDataModuleLabel = pset.get< std::string >("CalDataModuleLabel"                 );
    fDigitModuleLabel   = pset.get< std::string >("DigitModuleLabel", "FragmentToDigit");

    fCCAlg  .reconfigure(pset.get< fhicl::ParameterSet >("ClusterCrawlerAlg"));
    fCCHFAlg.reconfigure(pset.get< fhicl::ParameterSet >("CCHitFinderAlg")   );
  }
  
  void ClusterCrawlerLariat::beginJob(){
  }
  
  void ClusterCrawlerLariat::produce(art::Event & evt)
  {
    std::unique_ptr<art::Assns<recob::Wire,    recob::Hit>                     > wh_assn( new art::Assns<recob::Wire,    recob::Hit>                      );
    std::unique_ptr<art::Assns<recob::Cluster, recob::Hit>                     > hc_assn( new art::Assns<recob::Cluster, recob::Hit>                      );
//    std::unique_ptr<art::Assns<recob::Cluster, recob::Vertex, unsigned short>  > cv_assn( new art::Assns<recob::Cluster, recob::Vertex, unsigned short>   );
    std::unique_ptr<art::Assns<recob::Cluster, recob::Vertex>                  > cv_assn( new art::Assns<recob::Cluster, recob::Vertex>                   );
    std::unique_ptr<art::Assns<raw::Trigger,   recob::Hit>                     > th_assn( new art::Assns<raw::Trigger,   recob::Hit>                      );
    std::unique_ptr<art::Assns<raw::Trigger,   recob::Vertex>                  > tv_assn( new art::Assns<raw::Trigger,   recob::Vertex>                   );
    std::unique_ptr<art::Assns<raw::Trigger,   recob::Cluster>                 > tc_assn( new art::Assns<raw::Trigger,   recob::Cluster>                  );

    std::unique_ptr<std::vector<recob::Hit>      > hits     (new std::vector<recob::Hit>);
    std::unique_ptr<std::vector<recob::Cluster>  > clusters (new std::vector<recob::Cluster>);
    std::unique_ptr<std::vector<recob::Vertex>   > vertices (new std::vector<recob::Vertex>);
  
    art::ServiceHandle<geo::Geometry> geo;
    rdu::TriggerDigitUtility tdu(evt, fDigitModuleLabel);

    art::FindManyP<recob::Wire> fmpw(tdu.EventTriggersPtr(), evt, fCalDataModuleLabel);

    std::map< raw::ChannelID_t, art::Ptr<recob::Wire>    > chIDToWire;
    std::vector< art::Ptr<recob::Wire>           > wireVec;
    std::vector< std::vector<recob::Hit>         > clusToHits;
    std::vector< recob::Cluster                  > clus;
    std::map< size_t, std::vector<recob::Vertex> > clusToVertex;
    raw::ChannelID_t chid;

    hits->clear();

    for (size_t t=0; t<tdu.NTriggers(); ++t)                              
    {
       art::Ptr<raw::Trigger> trig = tdu.EventTriggersPtr()[t];            

       // fetch the wires needed by CCHitFinder
       chIDToWire.clear();
       wireVec.clear();
       wireVec = fmpw.at(t); 

       for(size_t w=0; w<wireVec.size(); ++w) 
       {
          chid=wireVec[w]->Channel(); 
          chIDToWire[chid] = wireVec[w];
       }

       ClusterCrawlerLariat::Clustering(wireVec, clus, clusToHits, clusToVertex);

//       size_t firstHitClus = hits->size();
       for(size_t ic = 0; ic < clus.size(); ++ic) 
       {  
          clusters->push_back(clus[ic]);     
          size_t startHitIdx = hits->size();
          for(size_t h = 0; h < clusToHits[ic].size(); ++h)
          {
	     hits->push_back(clusToHits[ic][h]);
             chid=hits->back().Channel();

             // make the wire - hit association
             if(!util::CreateAssn(*this, evt, *hits, chIDToWire[chid], *wh_assn, hits->size()-1))
             {
                throw art::Exception(art::errors::InsertFailure) <<"Failed to associate hit "<< h << " with wire ";
             } // exception

          }   
          size_t endHitIdx = hits->size();


          // make the cluster - vertices association
          for(size_t v = 0; v< clusToVertex[ic].size(); ++v)
          {
   	     vertices->push_back(clusToVertex[ic][v]);
             if(!util::CreateAssn(*this, evt, *clusters, *vertices, *cv_assn, vertices->size()-1, vertices->size(), ic))
             {
                throw art::Exception(art::errors::InsertFailure) <<"Failed to associate vertex "<< vertices->size()-1 << " with cluster "<<ic;
             } // exception
          }

          // make the cluster - hit association
          if(!util::CreateAssn(*this, evt, *clusters, *hits, *hc_assn, startHitIdx, endHitIdx, ic))
          {
             throw art::Exception(art::errors::InsertFailure) <<"Failed to associate hit "<<" with cluster "<<ic;
          } // exception

       } // Loop over clusters

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*
       // make the wire - hit association

       for(size_t h=firstHitClus; h<hits->size(); ++h) 
       {
          chid=hits->at(h).Channel(); 
          if(!util::CreateAssn(*this, evt, *hits, chIDToWire[chid], *wh_assn, h))
          {
             throw art::Exception(art::errors::InsertFailure) <<"Failed to associate hit "<< h << " with wire ";
          } // exception

       }
*/
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

       // make the trigger - cluster association
       for(size_t c = 0; c < clusters->size(); ++c)
       {
          if(!util::CreateAssn(*this, evt, *clusters, trig, *tc_assn, c))
          {
            throw art::Exception(art::errors::InsertFailure) <<"Failed to associate cluster "<< c << " with trigger "<<trig.key();
          } // exception
       }

       // make the trigger - hit association
       for(size_t h = 0; h < hits->size(); ++h)
       {
          if(!util::CreateAssn(*this, evt, *hits, trig, *th_assn, h))
          {
            throw art::Exception(art::errors::InsertFailure) <<"Failed to associate hit "<< h << " with trigger "<<trig.key();
          } // exception
       }

       // make the trigger - vertex association
       for(size_t v = 0; v < vertices->size(); ++v)
       {
          if(!util::CreateAssn(*this, evt, *vertices, trig, *tv_assn, v))
          {
            throw art::Exception(art::errors::InsertFailure) <<"Failed to associate vertex "<< v << " with trigger "<<trig.key();
          } // exception
       }

    } // Loop over triggers

    // move the hit collection and the associations into the event:
    evt.put(std::move(hits));
    evt.put(std::move(clusters));  
    evt.put(std::move(vertices));
    evt.put(std::move(wh_assn));
    evt.put(std::move(hc_assn));
    evt.put(std::move(cv_assn));   
    evt.put(std::move(th_assn));  
    evt.put(std::move(tc_assn)); 
    evt.put(std::move(tv_assn)); 
 } // produce

/////////////////////////////////////////////////////////////////////////
/////***************************************************************/////
/////***************************************************************/////
/////***************************************************************/////
/////////////////////////////////////////////////////////////////////////

 void ClusterCrawlerLariat::Clustering(std::vector< art::Ptr<recob::Wire> > & wireVec, 
				       std::vector<recob::Cluster> & clus, 
				       std::vector< std::vector<recob::Hit> > & clusToHits,
				       std::map< size_t, std::vector<recob::Vertex> > & clusToVertex)
 {  

//example of a std::map
// std::map<size_t, std::vector<reco::Vertex> > clusterToVertices;
// place something in the map, vertex is a recob::Vertex object
// clstr.IDx is the index of the current cluster
// map[clstr.IDx].push_back(vertex);

    std::vector<recob::Wire> wwires;
    wwires.clear();
    for(auto w : wireVec) wwires.push_back(*w);

    clusToHits.clear();
    clus.clear();
    clusToVertex.clear();

    // fetch the wires needed by CCHitFinder
    // find hits in all planes
    fCCHFAlg.RunCCHitFinder(wwires);
    
    // extract the result of the algorithm (it's moved)
    std::vector<recob::Hit> FirstHits = fCCHFAlg.YieldHits();

    // look for clusters in all planes
    fCCAlg.RunCrawler(FirstHits);
    
    // access to the algorithm results
    ClusterCrawlerAlg::HitInCluster_t const& HitInCluster = fCCAlg.GetHitInCluster();

    std::vector<recob::Hit> allHits;
    allHits.clear();
    for(auto hit : fCCAlg.YieldHits() ) allHits.push_back(hit);
    
    std::vector<ClusterCrawlerAlg::ClusterStore> const& tcl = fCCAlg.GetClusters();

    // make 3D vertices
    std::vector<ClusterCrawlerAlg::Vtx3Store> const& Vertx = fCCAlg.GetVertices();

    // make the clusters and associations
    float sumChg=0.0; 
    float sumADC=0.0;
    unsigned int vtxIndex = 0;
    unsigned int nclhits = 0;
    unsigned short plane = 0;
    double xyz[3] = {0, 0, 0};
    clusToHits.resize(tcl.size());
    for(unsigned int icl = 0; icl < tcl.size(); ++icl) 
    {
//       FinalHits.clear();
       ClusterCrawlerAlg::ClusterStore const& clstr = tcl[icl];
       if(clstr.ID < 0) continue;
       sumChg = 0;
       sumADC = 0;
       geo::PlaneID planeID = ClusterCrawlerAlg::DecodeCTP(clstr.CTP);
       plane = planeID.Plane;
       nclhits = clstr.tclhits.size();
       std::vector<unsigned int> clsHitIndices;
       // correct the hit indices to refer to the valid hits that were just added
       for(unsigned int itt = 0; itt < nclhits; ++itt) 
       {  
          clusToHits[icl].push_back(allHits.at(clstr.tclhits[itt]));
          sumChg += clusToHits[icl].back().Integral();
          sumADC += clusToHits[icl].back().SummedADC();
       } // itt

       // get the wire, plane from a hit
       geo::View_t view = clusToHits[icl].front().View();
       clus.emplace_back(
                           (float)clstr.BeginWir,  // Start wire
                           0,                      // sigma start wire
                           clstr.BeginTim,         // start tick
                           0,                      // sigma start tick
                           clstr.BeginChg,         // start charge
                           clstr.BeginAng,         // start angle
                           0,                      // start opening angle (0 for line-like clusters)
                           (float)clstr.EndWir,    // end wire
                           0,                      // sigma end wire
                           clstr.EndTim,           // end tick
                           0,                      // sigma end tick
                           clstr.EndChg,           // end charge
                           clstr.EndAng,           // end angle
                           0,                      // end opening angle (0 for line-like clusters)
                           sumChg,                 // integral
                           0,                      // sigma integral
                           sumADC,                 // summed ADC
                           0,                      // sigma summed ADC
                           nclhits,                // n hits
                           0,                      // wires over hits
                           0,                      // width (0 for line-like clusters)
                           clstr.ID,                  // ID
                           view,                   // view
                           planeID,                // plane
                           recob::Cluster::Sentry  // sentry
                         );

       // make the cluster - endpoint associations
       if(clstr.BeginVtx >= 0) 
       {
          // See if this endpoint is associated with a 3D vertex
          vtxIndex = 0;
          for(ClusterCrawlerAlg::Vtx3Store const& vtx3: Vertx) 
          {  
             // ignore incomplete vertices
             if(vtx3.Ptr2D[0] < 0) continue;
             if(vtx3.Ptr2D[1] < 0) continue;
             if(vtx3.Ptr2D[2] < 0) continue;
             if(vtx3.Ptr2D[plane] == clstr.BeginVtx) 
             {
                xyz[0] = vtx3.X;
                xyz[1] = vtx3.Y;
                xyz[2] = vtx3.Z;
                clusToVertex[icl].emplace_back(xyz, vtxIndex);
                ++vtxIndex;
             } // vertex match
          } // 3D vertices
       } // clstr.BeginVtx >= 0
       if(clstr.EndVtx >= 0) 
       {  
          // See if this endpoint is associated with a 3D vertex
           vtxIndex = 0;
          for(ClusterCrawlerAlg::Vtx3Store const& vtx3: Vertx) 
          {
             // ignore incomplete vertices
             if(vtx3.Ptr2D[0] < 0) continue;
             if(vtx3.Ptr2D[1] < 0) continue;
             if(vtx3.Ptr2D[2] < 0) continue;
             if(vtx3.Ptr2D[plane] == clstr.EndVtx) 
             {
                xyz[0] = vtx3.X;
                xyz[1] = vtx3.Y;
                xyz[2] = vtx3.Z;
                clusToVertex[icl].emplace_back(xyz, vtxIndex);
                ++vtxIndex;
             } // vertex match
          } // 3D vertices
       } // clstr.EndVtx >= 0
    } // icl
        
    // clean up
    fCCAlg.ClearResults();
 } //clustering function

} // namespace


namespace cluster{

  DEFINE_ART_MODULE(ClusterCrawlerLariat)
  
} 
