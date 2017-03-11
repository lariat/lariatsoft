//////////////////////////////////////////////////////////////////////////////////////////
// Class:       ClusterCrawlerT1034
// Module Type: producer
// File:        ClusterCrawlerT1034_module.cc
//
// Generated at Fri Jun  7 09:44:09 2013 by Bruce Baller using artmod 
// from cetpkgsupport v1_02_00.
//
//  ** Modified by Roberto Acciarri to account multiple trigger in one event for LArIAT.
//   
//  acciarri@fnal.gov
//  July 2015
//////////////////////////////////////////////////////////////////////////////////////////

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
#include "canvas/Persistency/Common/Ptr.h" 
#include "canvas/Persistency/Common/PtrVector.h" 
#include "art/Framework/Services/Registry/ServiceHandle.h" 
#include "art/Framework/Services/Optional/TFileService.h" 
#include "canvas/Persistency/Common/FindOneP.h"
#include "canvas/Utilities/Exception.h"
#include "canvas/Utilities/InputTag.h"

// ########################
// ### LArSoft Includes ###
// ########################
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "larcoreobj/SimpleTypesAndConstants/RawTypes.h" // raw::ChannelID_t
#include "larcore/Geometry/Geometry.h"
#include "larcore/Geometry/CryostatGeo.h"
#include "larcore/Geometry/TPCGeo.h"
#include "larcore/Geometry/PlaneGeo.h"
#include "larevt/Filters/ChannelFilter.h"
#include "RawDataUtilities/TriggerDigitUtility.h" 
#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RawData/raw.h"
#include "lardataobj/RecoBase/Wire.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "larreco/RecoAlg/CCHitFinderAlg.h"
#include "larreco/RecoAlg/ClusterCrawlerAlg.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "TComplex.h"
#include "TFile.h"
#include "TH2D.h"
#include "TF1.h"

// #include "RecoBase/EndPoint2D.h"
// #include "RecoBaseArt/HitCreator.h" // recob::HitCollectionAssociator
// #include "RecoAlg/ClusterRecoUtil/StandardClusterParamsAlg.h"
// #include "RecoAlg/ClusterParamsImportWrapper.h"


namespace cluster {
  class ClusterCrawlerT1034;
}

class cluster::ClusterCrawlerT1034 : public art::EDProducer {

  public:
    explicit ClusterCrawlerT1034(fhicl::ParameterSet const & pset);
    virtual ~ClusterCrawlerT1034();

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
		  std::map< size_t, std::vector<recob::Vertex> > & clusToVertex,
                  size_t & lastID);

};


namespace cluster {

  ClusterCrawlerT1034::ClusterCrawlerT1034(fhicl::ParameterSet const& pset)
    : fCCHFAlg(pset.get< fhicl::ParameterSet >("CCHitFinderAlg"))
    , fCCAlg  (pset.get< fhicl::ParameterSet >("ClusterCrawlerAlg"))
  {      
    
    this->reconfigure(pset);

    produces< std::vector<recob::Hit> >();     
    produces< std::vector<recob::Cluster> >();  
    produces< std::vector<recob::Vertex> >();
    produces< art::Assns<recob::Wire,    recob::Hit>     >();
    produces< art::Assns<recob::Cluster, recob::Hit>     >();
    produces< art::Assns<recob::Cluster, recob::Vertex>  >();
    produces< art::Assns<raw::Trigger,   recob::Hit>     >();
    produces< art::Assns<raw::Trigger,   recob::Vertex>  >();
    produces< art::Assns<raw::Trigger,   recob::Cluster> >();
  }

  ClusterCrawlerT1034::~ClusterCrawlerT1034()
  {
  }

  void ClusterCrawlerT1034::reconfigure(fhicl::ParameterSet const & pset)
  {
    fCalDataModuleLabel = pset.get< std::string >("CalDataModuleLabel"                 );
    fDigitModuleLabel   = pset.get< std::string >("DigitModuleLabel", "FragmentToDigit");

    fCCAlg  .reconfigure(pset.get< fhicl::ParameterSet >("ClusterCrawlerAlg"));
    fCCHFAlg.reconfigure(pset.get< fhicl::ParameterSet >("CCHitFinderAlg")   );
  }
  
  void ClusterCrawlerT1034::beginJob(){
  }
  
  void ClusterCrawlerT1034::produce(art::Event & evt)
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

    //clearing vectors before to start. Just to stay on the safe side
    hits->clear();
    clusters->clear();
    vertices->clear();

    size_t startHit = 0;
    size_t startCluster = 0;
    size_t startVertex = 0;
    size_t lastID = 0;

    LOG_VERBATIM("Summary") << "Number of Triggers: " << tdu.NTriggers();
    for (size_t t=0; t<tdu.NTriggers(); ++t)                              
    {
       art::Ptr<raw::Trigger> trig = tdu.EventTriggersPtr()[t];  

       // Skip trigger if empty
       art::PtrVector<raw::RawDigit> rdvec = tdu.TriggerRawDigitsPtr(t);

       LOG_VERBATIM("ClusterCrawlerT1034") << "#####################################################"
                                            << "\n #####################################################"                 
                                            << "\n Trigger Number: " << t << "   Raw Digit vector size: "<< rdvec.size();
       if(!rdvec.size()){mf::LogInfo("ClusterCrawlerT1034") << " Raw Digit vector is empty. Skipping the trigger"; continue;}
       LOG_VERBATIM("ClusterCrawlerT1034") << " ";

       // get the starting index of the hits, clusters and vertices for this trigger
       startHit = hits->size();
       startCluster = clusters->size();
       startVertex = vertices->size();

       // fetch the wires needed by CCHitFinder
       chIDToWire.clear();
       wireVec.clear();
       wireVec = fmpw.at(t); 

       for(size_t w=0; w<wireVec.size(); ++w) 
       {
          chid=wireVec[w]->Channel(); 
          chIDToWire[chid] = wireVec[w];
       }

       ClusterCrawlerT1034::Clustering(wireVec, clus, clusToHits, clusToVertex, lastID);

       for(size_t ic = 0; ic < clus.size(); ++ic) 
       {  
          clusters->push_back(clus[ic]);   

          size_t startHitIdx = hits->size();
          for(size_t h = 0; h < clusToHits[ic].size(); ++h)
          {
	     hits->push_back(clusToHits[ic][h]);
             chid=hits->back().Channel();
              
             // make the wire - hit association
             if(!util::CreateAssn(*this, evt, *hits, chIDToWire[chid], *wh_assn))
             {
                throw art::Exception(art::errors::ProductRegistrationFailure) <<"Failed to associate hit "<< h << " with wire ";
             } // exception

          }   
          size_t endHitIdx = hits->size();

          // make the cluster - vertices association
          for(size_t v = 0; v< clusToVertex[ic].size(); ++v)
          {
   	     vertices->push_back(clusToVertex[ic][v]);
             if(!util::CreateAssn(*this, evt, *clusters, *vertices, *cv_assn, vertices->size()-1, vertices->size()))
             {
                throw art::Exception(art::errors::ProductRegistrationFailure) <<"Failed to associate vertex "<< vertices->size()-1 << " with cluster "<<ic;
             } // exception
          }
 
          // make the cluster - hit association
          if(!util::CreateAssn(*this, evt, *clusters, *hits, *hc_assn, startHitIdx, endHitIdx))
          {
             throw art::Exception(art::errors::ProductRegistrationFailure) <<"Failed to associate hit "<<" with cluster "<<ic;
          } // exception

       } // Loop over clusters

       lastID=clusters->back().ID()+1;

       // make the trigger - cluster association
       for(size_t c = startCluster; c < clusters->size(); ++c)
       {
          if(!util::CreateAssn(*this, evt, *clusters, trig, *tc_assn, c))
          {
            throw art::Exception(art::errors::ProductRegistrationFailure) <<"Failed to associate cluster "<< c << " with trigger "<<trig.key();
          } // exception
       }

       // make the trigger - hit association
       for(size_t h = startHit; h < hits->size(); ++h)
       {
          if(!util::CreateAssn(*this, evt, *hits, trig, *th_assn, h))
          {
            throw art::Exception(art::errors::ProductRegistrationFailure) <<"Failed to associate hit "<< h << " with trigger "<<trig.key();
          } // exception
       }

       // make the trigger - vertex association
       for(size_t v = startVertex; v < vertices->size(); ++v)
       {
          if(!util::CreateAssn(*this, evt, *vertices, trig, *tv_assn, v))
          {
            throw art::Exception(art::errors::ProductRegistrationFailure) <<"Failed to associate vertex "<< v << " with trigger "<<trig.key();
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

 void ClusterCrawlerT1034::Clustering(std::vector< art::Ptr<recob::Wire> > & wireVec, 
				       std::vector<recob::Cluster> & clus, 
				       std::vector< std::vector<recob::Hit> > & clusToHits,
				       std::map< size_t, std::vector<recob::Vertex> > & clusToVertex,
                                       size_t & lastID)
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
    auto const& HitInCluster = fCCAlg.GetinClus();

    std::vector<recob::Hit> allHits;
    allHits.clear();
    for(auto hit : fCCAlg.YieldHits() ) allHits.push_back(hit);

    std::vector<ClusterCrawlerAlg::ClusterStore> const& tcl = fCCAlg.GetClusters();

    // Consistency check
    for(unsigned int icl = 0; icl < tcl.size(); ++icl) 
    {
       ClusterCrawlerAlg::ClusterStore const& clstr = tcl[icl];
       if(clstr.ID < 0) continue;
       geo::PlaneID planeID = ClusterCrawlerAlg::DecodeCTP(clstr.CTP);
       unsigned short plane = planeID.Plane;
       for(unsigned short ii = 0; ii < clstr.tclhits.size(); ++ii) 
       {
          unsigned int iht = clstr.tclhits[ii];
          recob::Hit const& theHit = allHits.at(iht);
          if(theHit.WireID().Plane != plane) 
          {
             std::cout<<"CC: cluster-hit plane mis-match "<<theHit.WireID().Plane<<" "<< plane <<" in cluster "<< clstr.ID <<" WT "<< clstr.BeginWir <<" : "<<(int)clstr.BeginTim << "\n";
             return;
          }
          if(HitInCluster[iht] != clstr.ID) 
          {
             std::cout << "CC: InClus mis-match " << HitInCluster[iht] << " ID " << clstr.ID << " in cluster " << icl << "\n";
             return;
          }
       } // ii
    } // icl

    // make 3D vertices
    std::vector<ClusterCrawlerAlg::Vtx3Store> const& Vertx = fCCAlg.GetVertices();

    // make the clusters and associations
    float sumChg          = 0.0; 
    float sumADC          = 0.0;
    unsigned int idx      =   0;  //index for clusToHits
    unsigned int vtxIndex =   0;
    unsigned int nclhits  =   0;
    unsigned short plane  =   0;
    double xyz[3]   = {0, 0, 0};
    clusToHits.resize(tcl.size());
    LOG_VERBATIM("ClusterCrawlerT1034") << " Number of Clusters in the trigger: " <<tcl.size();
    LOG_VERBATIM("ClusterCrawlerT1034") << " ";
 
    for(unsigned int icl = 0; icl < tcl.size(); ++icl) 
    {
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
          clusToHits[idx].push_back(allHits.at(clstr.tclhits[itt]));
          sumChg += clusToHits[idx].back().Integral();
          sumADC += clusToHits[idx].back().SummedADC();
       } // itt

       // get the wire, plane from a hit
       geo::View_t view = clusToHits[idx].front().View();
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
                           clstr.ID+lastID,        // ID
                           view,                   // view
                           planeID,                // plane
                           recob::Cluster::Sentry  // sentry
                         );

          LOG_VERBATIM("ClusterCrawlerT1034") << " /////////////////////////////////////// "
                                               << "\n /////////////////////////////////////// "
                                               << "\n Cluster ID " << clus.back().ID() << "   Number of hits " << nclhits << "   View " << view
                                               << "\n Cluster Info: Start Wire " << clstr.BeginWir << "   End Wire " << clstr.EndWir
                                               << "\n Cluster Info: Begin Time " << clstr.BeginTim << "   End Time " << clstr.EndTim
                                               << "\n Cluster Info: Charge " << sumChg << "   ADC Sum " << sumADC
                                               << "\n Cluster Info: Begin Charge " << clstr.BeginChg << "   End Charge " << clstr.EndChg
                                               << "\n Cluster Info: Begin Angle " << clstr.BeginAng << "   End Angle " << clstr.EndAng
                                               << "\n Begin Vertex Index " << clstr.BeginVtx <<  "   End Vertex Index " << clstr.EndVtx
                                               << "\n  ";

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
                clusToVertex[idx].emplace_back(xyz, vtxIndex);
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
                clusToVertex[idx].emplace_back(xyz, vtxIndex);
                ++vtxIndex;
             } // vertex match
          } // 3D vertices
       } // clstr.EndVtx >= 0
       idx++;
    } // icl
    // assign the right size to clusToHits
    clusToHits.resize(idx);
    // clean up
    fCCAlg.ClearResults();
 } //clustering function

} // namespace


namespace cluster{

  DEFINE_ART_MODULE(ClusterCrawlerT1034)
  
} 
