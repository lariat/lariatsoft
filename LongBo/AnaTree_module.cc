////////////////////////////////////////////////////////////////////////
// Class:       AnaTree
// Module Type: analyzer
// File:        AnaTree_module.cc
//
// Generated at Sun Mar 24 09:05:02 2013 by Tingjun Yang using artmod
// from art v1_02_06.
////////////////////////////////////////////////////////////////////////
// Framework includes
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h" 
#include "art/Framework/Principal/Event.h" 
#include "fhiclcpp/ParameterSet.h" 
#include "art/Framework/Principal/Handle.h" 
#include "art/Persistency/Common/Ptr.h" 
#include "art/Persistency/Common/PtrVector.h" 
#include "art/Framework/Services/Registry/ServiceHandle.h" 
#include "art/Framework/Services/Optional/TFileService.h" 
#include "art/Framework/Services/Optional/TFileDirectory.h" 
#include "messagefacility/MessageLogger/MessageLogger.h" 

// LArSoft includes
#include "Geometry/Geometry.h"
#include "Geometry/PlaneGeo.h"
#include "Geometry/WireGeo.h"
#include "RecoBase/Hit.h"
#include "RecoBase/Cluster.h"
#include "RecoBase/Track.h"
#include "RecoBase/SpacePoint.h"
#include "Utilities/LArProperties.h"
#include "Utilities/DetectorProperties.h"
#include "Utilities/AssociationUtil.h"
#include "RawData/ExternalTrigger.h"

// ROOT includes
#include "TTree.h"
#include "TTimeStamp.h"

const int kMaxTrack      = 1000;  //maximum number of tracks
const int kMaxHits       = 10000; //maximum number of hits
const int kMaxTrackHits  = 1000;  //maximum number of space points

namespace bo {
  class AnaTree;
}

class bo::AnaTree : public art::EDAnalyzer {
public:
  explicit AnaTree(fhicl::ParameterSet const & p);
  virtual ~AnaTree();

  void analyze(art::Event const & e) override;

  void beginJob() override;
  //void reconfigure(fhicl::ParameterSet const & p) override;

private:

  void ResetVars();
  
  TTree* fTree;
  //run information
  int run;
  int subrun;
  int event;
  double evttime;
  double efield[3];
  int t0;
  int trigtime[16];
  int ntracks_reco;         //number of reconstructed tracks
  double trkvtxx[kMaxTrack];
  double trkvtxy[kMaxTrack];
  double trkvtxz[kMaxTrack];
  double trkendx[kMaxTrack];
  double trkendy[kMaxTrack];
  double trkendz[kMaxTrack];
  double trkstartdcosx[kMaxTrack];
  double trkstartdcosy[kMaxTrack];
  double trkstartdcosz[kMaxTrack];
  double trkenddcosx[kMaxTrack];
  double trkenddcosy[kMaxTrack];
  double trkenddcosz[kMaxTrack];
  int    ntrkhits[kMaxTrack];
  double trkx[kMaxTrack][kMaxTrackHits];
  double trky[kMaxTrack][kMaxTrackHits];
  double trkz[kMaxTrack][kMaxTrackHits];
  double trkpitch[kMaxTrack][3];
  int    nhits;
  int    hit_plane[kMaxHits];
  int    hit_wire[kMaxHits];
  int    hit_channel[kMaxHits];
  double hit_peakT[kMaxHits];
  double hit_charge[kMaxHits];
  double hit_ph[kMaxHits];
  double hit_tstart[kMaxHits];
  double hit_tend[kMaxHits];
  double hit_mult[kMaxHits];
  double hit_goodness[kMaxHits];
  int    hit_trkid[kMaxHits];

  std::string fTrigModuleLabel;
  std::string fHitsModuleLabel;
  std::string fTrackModuleLabel;


};


bo::AnaTree::AnaTree(fhicl::ParameterSet const & pset)
  : EDAnalyzer(pset)
  , fTrigModuleLabel       (pset.get< std::string >("TrigModuleLabel"))
  , fHitsModuleLabel       (pset.get< std::string >("HitsModuleLabel"))
  , fTrackModuleLabel       (pset.get< std::string >("TrackModuleLabel"))
{
}

bo::AnaTree::~AnaTree()
{
  // Clean up dynamic memory and other resources here.
}

void bo::AnaTree::analyze(art::Event const & evt)
{
  // Implementation of required member function here.
  ResetVars();

  art::ServiceHandle<geo::Geometry> geom;
  art::ServiceHandle<util::LArProperties> larprop;
  art::ServiceHandle<util::DetectorProperties> detprop;
  
  run = evt.run();
  subrun = evt.subRun();
  event = evt.id().event();

  art::Timestamp ts = evt.time();
  TTimeStamp tts(ts.timeHigh(), ts.timeLow());
  evttime = tts.AsDouble();


  efield[0] = larprop->Efield(0);
  efield[1] = larprop->Efield(1);
  efield[2] = larprop->Efield(2);

  t0 = detprop->TriggerOffset();

  art::Handle< std::vector<raw::ExternalTrigger> > trigListHandle;
  std::vector<art::Ptr<raw::ExternalTrigger> > triglist;
  if (evt.getByLabel(fTrigModuleLabel,trigListHandle))
    art::fill_ptr_vector(triglist, trigListHandle);

  for (size_t i = 0; i<triglist.size(); ++i){
    trigtime[i] = triglist[i]->GetTrigTime();
  }

  art::Handle< std::vector<recob::Track> > trackListHandle;
  std::vector<art::Ptr<recob::Track> > tracklist;
  if (evt.getByLabel(fTrackModuleLabel,trackListHandle))
    art::fill_ptr_vector(tracklist, trackListHandle);    

  art::Handle< std::vector<recob::Hit> > hitListHandle;
  std::vector<art::Ptr<recob::Hit> > hitlist;
  if (evt.getByLabel(fHitsModuleLabel,hitListHandle))
    art::fill_ptr_vector(hitlist, hitListHandle);


  art::FindManyP<recob::SpacePoint> fmsp(trackListHandle, evt, fTrackModuleLabel);
  art::FindMany<recob::Track>       fmtk(hitListHandle, evt, fTrackModuleLabel);

  //track information
  ntracks_reco=tracklist.size();
  double larStart[3];
  double larEnd[3];
  std::vector<double> trackStart;
  std::vector<double> trackEnd;
  for(size_t i=0; i<tracklist.size();++i){
    trackStart.clear();
    trackEnd.clear();
    memset(larStart, 0, 3);
    memset(larEnd, 0, 3);
    tracklist[i]->Extent(trackStart,trackEnd); 
    tracklist[i]->Direction(larStart,larEnd);
    trkvtxx[i]        = trackStart[0];
    trkvtxy[i]        = trackStart[1];
    trkvtxz[i]        = trackStart[2];
    trkendx[i]        = trackEnd[0];
    trkendy[i]        = trackEnd[1];
    trkendz[i]        = trackEnd[2];
    trkstartdcosx[i]  = larStart[0];
    trkstartdcosy[i]  = larStart[1];
    trkstartdcosz[i]  = larStart[2];
    trkenddcosx[i]    = larEnd[0];
    trkenddcosy[i]    = larEnd[1];
    trkenddcosz[i]    = larEnd[2];
    ntrkhits[i] = fmsp.at(i).size();
    std::vector<art::Ptr<recob::SpacePoint> > spts = fmsp.at(i);
    for (size_t j = 0; j<spts.size(); ++j){
      trkx[i][j] = spts[j]->XYZ()[0];
      trky[i][j] = spts[j]->XYZ()[1];
      trkz[i][j] = spts[j]->XYZ()[2];
    }
    for (int j = 0; j<3; ++j){
      try{
	if (j==0)
	  trkpitch[i][j] = tracklist[i]->PitchInView(geo::kU);
	else if (j==1)
	  trkpitch[i][j] = tracklist[i]->PitchInView(geo::kV);
	else if (j==2)
	  trkpitch[i][j] = tracklist[i]->PitchInView(geo::kZ);
      }
      catch( cet::exception &e){
	mf::LogWarning("AnaTree")<<"caught exeption "<<e<<"\n setting pitch to 0";
	trkpitch[i][j] = 0;
      }
    }
  }

  nhits = hitlist.size();
  for (size_t i = 0; i<hitlist.size(); ++i){
    unsigned int channel = hitlist[i]->Channel();
    geo::WireID wireid = hitlist[i]->WireID();
    hit_plane[i]   = wireid.Plane;
    hit_wire[i]    = wireid.Wire;
    hit_channel[i] = channel;
    hit_peakT[i]   = hitlist[i]->PeakTime();
    hit_charge[i]  = hitlist[i]->Charge();
    hit_ph[i]      = hitlist[i]->Charge(true);
    hit_tstart[i]  = hitlist[i]->StartTime();
    hit_tend[i]    = hitlist[i]->EndTime();
    hit_mult[i]    = hitlist[i]->Multiplicity();
    hit_goodness[i]= hitlist[i]->GoodnessOfFit();
    if (fmtk.at(i).size()!=0){
      hit_trkid[i] = fmtk.at(i)[0]->ID();
    }
  }
  fTree->Fill();

}

void bo::AnaTree::beginJob()
{
  // Implementation of optional member function here.
  art::ServiceHandle<art::TFileService> tfs;
  fTree = tfs->make<TTree>("anatree","analysis tree");
  fTree->Branch("run",&run,"run/I");
  fTree->Branch("subrun",&subrun,"subrun/I");
  fTree->Branch("event",&event,"event/I");
  fTree->Branch("evttime",&evttime,"evttime/D");
  fTree->Branch("efield",efield,"efield[3]/D");
  fTree->Branch("t0",&t0,"t0/I");
  fTree->Branch("trigtime",trigtime,"trigtime[16]/I");
  fTree->Branch("ntracks_reco",&ntracks_reco,"ntracks_reco/I");
  fTree->Branch("trkvtxx",trkvtxx,"trkvtxx[ntracks_reco]/D");
  fTree->Branch("trkvtxy",trkvtxy,"trkvtxy[ntracks_reco]/D");
  fTree->Branch("trkvtxz",trkvtxz,"trkvtxz[ntracks_reco]/D");
  fTree->Branch("trkendx",trkendx,"trkendx[ntracks_reco]/D");
  fTree->Branch("trkendy",trkendy,"trkendy[ntracks_reco]/D");
  fTree->Branch("trkendz",trkendz,"trkendz[ntracks_reco]/D");
  fTree->Branch("trkstartdcosx",trkstartdcosx,"trkstartdcosx[ntracks_reco]/D");
  fTree->Branch("trkstartdcosy",trkstartdcosy,"trkstartdcosy[ntracks_reco]/D");
  fTree->Branch("trkstartdcosz",trkstartdcosz,"trkstartdcosz[ntracks_reco]/D");
  fTree->Branch("trkenddcosx",trkenddcosx,"trkenddcosx[ntracks_reco]/D");
  fTree->Branch("trkenddcosy",trkenddcosy,"trkenddcosy[ntracks_reco]/D");
  fTree->Branch("trkenddcosz",trkenddcosz,"trkenddcosz[ntracks_reco]/D");
  fTree->Branch("ntrkhits",ntrkhits,"ntrkhits[ntracks_reco]/I");
  fTree->Branch("trkx",trkx,"trkx[ntracks_reco][1000]/D");
  fTree->Branch("trky",trky,"trky[ntracks_reco][1000]/D");
  fTree->Branch("trkz",trkz,"trkz[ntracks_reco][1000]/D");
  fTree->Branch("trkpitch",trkpitch,"trkpitch[ntracks_reco][3]/D");
  fTree->Branch("nhits",&nhits,"nhits/I");
  fTree->Branch("hit_plane",hit_plane,"hit_plane[nhits]/I");
  fTree->Branch("hit_wire",hit_wire,"hit_wire[nhits]/I");
  fTree->Branch("hit_channel",hit_channel,"hit_channel[nhits]/I");
  fTree->Branch("hit_peakT",hit_peakT,"hit_peakT[nhits]/D");
  fTree->Branch("hit_charge",hit_charge,"hit_charge[nhits]/D");
  fTree->Branch("hit_ph",hit_ph,"hit_ph[nhits]/D");
  fTree->Branch("hit_tstart",hit_tstart,"hit_tstart[nhits]/D");
  fTree->Branch("hit_tend",hit_tend,"hit_tend[nhits]/D");
  fTree->Branch("hit_mult",hit_mult,"hit_mult[nhits]/D");
  fTree->Branch("hit_goodness",hit_goodness,"hit_goodness[nhits]/D");
  fTree->Branch("hit_trkid",hit_trkid,"hit_trkid[nhits]/I");
}

//void bo::AnaTree::reconfigure(fhicl::ParameterSet const & p)
//{
//  // Implementation of optional member function here.
//}
void bo::AnaTree::ResetVars(){

  run = -99999;
  subrun = -99999;
  event = -99999;
  evttime = -99999;
  for (int i = 0; i<3; ++i){
    efield[i] = -99999;
  }
  t0 = -99999;
  ntracks_reco = -99999;
  for (int i = 0; i < kMaxTrack; ++i){
    trkvtxx[i] = -99999;
    trkvtxy[i] = -99999;
    trkvtxz[i] = -99999;
    trkendx[i] = -99999;
    trkendy[i] = -99999;
    trkendz[i] = -99999;
    trkstartdcosx[i] = -99999;
    trkstartdcosy[i] = -99999;
    trkstartdcosz[i] = -99999;
    trkenddcosx[i] = -99999;
    trkenddcosy[i] = -99999;
    trkenddcosz[i] = -99999;
    ntrkhits[i] = -99999;
    for (int j = 0; j<kMaxTrackHits; ++j){
      trkx[i][j] = -99999;
      trky[i][j] = -99999;
      trkz[i][j] = -99999;
    }
    for (int j = 0; j<3; ++j){
      trkpitch[i][j] = -99999;
    }
  }
  nhits = -99999;
  for (int i = 0; i<kMaxHits; ++i){
    hit_plane[i] = -99999;
    hit_wire[i] = -99999;
    hit_channel[i] = -99999;
    hit_peakT[i] = -99999;
    hit_charge[i] = -99999;
    hit_ph[i] = -99999;
    hit_trkid[i] = -99999;
    hit_tstart[i] = -99999;
    hit_tend[i] = -99999;
    hit_mult[i] = -99999;
    hit_goodness[i] = -99999;
  }
  
}

DEFINE_ART_MODULE(bo::AnaTree)
