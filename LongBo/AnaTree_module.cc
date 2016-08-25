////////////////////////////////////////////////////////////////////////
// Class:       AnaTree
// Module Type: analyzer
// File:        AnaTree_module.cc
//
// Generated at Sun Mar 24 09:05:02 2013 by Tingjun Yang using artmod
// from art v1_02_06.
////////////////////////////////////////////////////////////////////////
// Framework includes
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/View.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/PtrVector.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "canvas/Persistency/Common/FindMany.h"
#include "canvas/Utilities/InputTag.h"
#include "messagefacility/MessageLogger/MessageLogger.h" 
#include "fhiclcpp/ParameterSet.h" 
#include "cetlib/maybe_ref.h"

// LArSoft includes
#include "larcore/Geometry/Geometry.h"
#include "larcore/Geometry/PlaneGeo.h"
#include "larcore/Geometry/WireGeo.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardata/RecoBaseArt/TrackUtils.h" // lar::util::TrackPitchInView()
#include "lardata/DetectorInfoServices/LArPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/Utilities/AssociationUtil.h"
//#include "RawData/ExternalTrigger.h"
#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RawData/raw.h"
#include "larsim/MCCheater/BackTracker.h"
#include "larsimobj/Simulation/SimChannel.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "lardataobj/AnalysisBase/ParticleID.h"
#include "larreco/RecoAlg/TrackMomentumCalculator.h"

// ROOT includes
#include "TTree.h"
#include "TTimeStamp.h"

const int kMaxTrack      = 1000;  //maximum number of tracks
const int kMaxCluster    = 1000;  //maximum number of clusters
const int kMaxHits       = 20000; //maximum number of hits
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
  //int trigtime[16];
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
  double trklen[kMaxTrack];
  double trkmomrange[kMaxTrack];
  double trkmommschi2[kMaxTrack];
  double trkmommsllhd[kMaxTrack];
  int    ntrkhits[kMaxTrack];
  double trkx[kMaxTrack][kMaxTrackHits];
  double trky[kMaxTrack][kMaxTrackHits];
  double trkz[kMaxTrack][kMaxTrackHits];
  double trkpitch[kMaxTrack][2];
  int    trkhits[kMaxTrack][2];
  double trkpida[kMaxTrack][2];
  double trkke[kMaxTrack][2];
  double trkdedx[kMaxTrack][2][1000];
  double trkrr[kMaxTrack][2][1000];
  double trkpitchhit[kMaxTrack][2][1000];
  int    nclus;
  double clustartwire[kMaxCluster];
  double clustarttick[kMaxCluster];
  double cluendwire[kMaxCluster];
  double cluendtick[kMaxCluster];
  int    cluplane[kMaxCluster];
  int    nhits;
  int    hit_plane[kMaxHits];
  int    hit_wire[kMaxHits];
  int    hit_channel[kMaxHits];
  double hit_peakT[kMaxHits];
  double hit_charge[kMaxHits];
  double hit_ph[kMaxHits];
  double hit_tstart[kMaxHits];
  double hit_tend[kMaxHits];
  int    hit_trkid[kMaxHits];
  int    hit_pk[kMaxHits];
  int    hit_t[kMaxHits];
  int    hit_ch[kMaxHits];
  int    hit_fwhh[kMaxHits];
  double hit_rms[kMaxHits];
  double hit_nelec[kMaxHits];
  double hit_energy[kMaxHits];
  int    hit_trkkey[kMaxHits];
  int    hit_clukey[kMaxHits];

  //  std::string fTrigModuleLabel;
  std::string fHitsModuleLabel;
  std::string fClusterModuleLabel;
  std::string fTrackModuleLabel;
  std::string fCalorimetryModuleLabel;
  std::string fParticleIDModuleLabel;
};


bo::AnaTree::AnaTree(fhicl::ParameterSet const & pset)
  : EDAnalyzer(pset)
    //  , fTrigModuleLabel       (pset.get< std::string >("TrigModuleLabel"))
  , fHitsModuleLabel       (pset.get< std::string >("HitsModuleLabel"))
  , fClusterModuleLabel     (pset.get< std::string >("ClusterModuleLabel"))
  , fTrackModuleLabel       (pset.get< std::string >("TrackModuleLabel"))
  , fCalorimetryModuleLabel (pset.get< std::string >("CalorimetryModuleLabel"))
  , fParticleIDModuleLabel  (pset.get< std::string >("ParticleIDModuleLabel"))
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
  //auto const* larprop = lar::providerFrom<detinfo::LArPropertiesService>();
  auto const* detprop = lar::providerFrom<detinfo::DetectorPropertiesService>();
  art::ServiceHandle<cheat::BackTracker> bt;

//  for (size_t i = 0; i<geom->Nplanes(); ++i){
//    for (size_t j = 0; j<geom->Plane(i).Nwires(); ++j){
//      double xyzStart[3];
//      double xyzEnd[3];
//      geom->WireEndPoints(0,0,i,j,xyzStart,xyzEnd);
//      std::cout<<i<<" "<<j<<" "<<xyzStart[0]<<" "<<xyzStart[1]<<" "<<xyzStart[2]<<" "<<xyzEnd[0]<<" "<<xyzEnd[1]<<" "<<xyzEnd[2]<<std::endl;
//    }
//  }
  
  run = evt.run();
  subrun = evt.subRun();
  event = evt.id().event();

  art::Timestamp ts = evt.time();
  TTimeStamp tts(ts.timeHigh(), ts.timeLow());
  evttime = tts.AsDouble();


  // Note: LArProperties::Efield() has moved to DetectorProperties/DetectorPropertiesService
  efield[0] = detprop->Efield(0);
  efield[1] = detprop->Efield(1);
  efield[2] = detprop->Efield(2);

  t0 = detprop->TriggerOffset();

  /*
  art::Handle< std::vector<raw::ExternalTrigger> > trigListHandle;
  std::vector<art::Ptr<raw::ExternalTrigger> > triglist;
  if (evt.getByLabel(fTrigModuleLabel,trigListHandle))
    art::fill_ptr_vector(triglist, trigListHandle);

  for (size_t i = 0; i<triglist.size(); ++i){
    trigtime[i] = triglist[i]->GetTrigTime();
  }
  */

  art::Handle< std::vector<recob::Cluster> > clusterListHandle;
  std::vector<art::Ptr<recob::Cluster> > clusterlist;
  if (evt.getByLabel(fClusterModuleLabel,clusterListHandle))
    art::fill_ptr_vector(clusterlist, clusterListHandle);    


  art::Handle< std::vector<recob::Track> > trackListHandle;
  std::vector<art::Ptr<recob::Track> > tracklist;
  if (evt.getByLabel(fTrackModuleLabel,trackListHandle))
    art::fill_ptr_vector(tracklist, trackListHandle);    

  art::Handle< std::vector<recob::Hit> > hitListHandle;
  std::vector<art::Ptr<recob::Hit> > hitlist;
  if (evt.getByLabel(fHitsModuleLabel,hitListHandle))
    art::fill_ptr_vector(hitlist, hitListHandle);

  art::FindOne<raw::RawDigit>       ford(hitListHandle,   evt, fHitsModuleLabel);
  art::FindManyP<recob::SpacePoint> fmsp(trackListHandle, evt, fTrackModuleLabel);
  art::FindManyP<recob::Track>       fmtk(hitListHandle,   evt, fTrackModuleLabel);
  art::FindManyP<recob::Cluster>     fmc(hitListHandle,   evt, fClusterModuleLabel);
  art::FindManyP<anab::Calorimetry> fmcal(trackListHandle, evt, fCalorimetryModuleLabel);
  art::FindManyP<anab::ParticleID>  fmpid(trackListHandle, evt, fParticleIDModuleLabel);
  std::vector<const sim::SimChannel*> fSimChannels;
  try{
    evt.getView("largeant", fSimChannels);
  }catch (art::Exception const&e){
  }

  //cluster information
  nclus = clusterlist.size();
  for (size_t i = 0; i<clusterlist.size(); ++i){
    clustartwire[i] = clusterlist[i]->StartWire();
    clustarttick[i] = clusterlist[i]->StartTick();
    cluendwire[i] = clusterlist[i]->EndWire();
    cluendtick[i] = clusterlist[i]->EndTick();
    cluplane[i] = clusterlist[i]->Plane().Plane;
  }
  //track information
  trkf::TrackMomentumCalculator trkm;
  trkm.SetMinLength(10); //change the minimal track length requirement to 10 cm
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
    trklen[i]         = tracklist[i]->Length();
    trkmomrange[i]    = trkm.GetTrackMomentum(trklen[i],13);
    trkmommschi2[i]   = trkm.GetMomentumMultiScatterChi2(tracklist[i]);
    trkmommsllhd[i]   = trkm.GetMomentumMultiScatterLLHD(tracklist[i]);
    ntrkhits[i] = fmsp.at(i).size();
    std::vector<art::Ptr<recob::SpacePoint> > spts = fmsp.at(i);
    for (size_t j = 0; j<spts.size(); ++j){
      trkx[i][j] = spts[j]->XYZ()[0];
      trky[i][j] = spts[j]->XYZ()[1];
      trkz[i][j] = spts[j]->XYZ()[2];
    }
    for (int j = 0; j<2; ++j){
      try{
	if (j==0)
	  trkpitch[i][j] = lar::util::TrackPitchInView(*(tracklist[i]), geo::kU);
	else if (j==1)
	  trkpitch[i][j] = lar::util::TrackPitchInView(*(tracklist[i]), geo::kV);
      }
      catch( cet::exception &e){
	mf::LogWarning("AnaTree")<<"caught exeption "<<e<<"\n setting pitch to 0";
	trkpitch[i][j] = 0;
      }
    }
    if (fmcal.isValid()){
      std::vector<art::Ptr<anab::Calorimetry> > calos = fmcal.at(i);
      for (size_t j = 0; j<calos.size(); ++j){
	if (!calos[j]->PlaneID().isValid) continue;
	int pl = calos[j]->PlaneID().Plane;
	if (pl<0||pl>1) continue;
	trkhits[i][pl] = calos[j]->dEdx().size();
	trkke[i][pl] = calos[j]->KineticEnergy();
	for (size_t k = 0; k<calos[j]->dEdx().size(); ++k){
	  if (k>=1000) continue;
	  trkdedx[i][pl][k] = calos[j]->dEdx()[k];
	  trkrr[i][pl][k] = calos[j]->ResidualRange()[k];
	  trkpitchhit[i][pl][k] = calos[j]->TrkPitchVec()[k];
	}
      }
    }
    if (fmpid.isValid()){
      std::vector<art::Ptr<anab::ParticleID> > pids = fmpid.at(i);
      for (size_t j = 0; j<pids.size(); ++j){
	if (!pids[j]->PlaneID().isValid) continue;
	int pl = pids[j]->PlaneID().Plane;
	if (pl<0||pl>1) continue;
	trkpida[i][pl] = pids[j]->PIDA();
      }
    }
  }

  nhits = hitlist.size();
  for (size_t i = 0; i<hitlist.size(); ++i){
    cet::maybe_ref<raw::RawDigit const> rdref(ford.at(i));
    unsigned int channel = hitlist[i]->Channel();
    geo::WireID wireid = hitlist[i]->WireID();
    hit_plane[i]   = wireid.Plane;
    hit_wire[i]    = wireid.Wire;
    hit_channel[i] = channel;
    hit_peakT[i]   = hitlist[i]->PeakTime();
    hit_charge[i]  = hitlist[i]->Integral();
    hit_ph[i]      = hitlist[i]->PeakAmplitude();
    hit_tstart[i]  = hitlist[i]->StartTick();
    hit_tend[i]    = hitlist[i]->EndTick();
    if (fmtk.isValid()){
      if (fmtk.at(i).size()!=0){
	hit_trkid[i] = fmtk.at(i)[0]->ID();
	hit_trkkey[i] = fmtk.at(i)[0].key();
      }
    }
    if (fmc.isValid()){
      if (fmc.at(i).size()!=0){
	hit_clukey[i] = fmc.at(i)[0].key();
      }
    }
    if (hit_plane[i]==1){//collection plane
      if( rdref.isValid() ){
	raw::RawDigit const& rd(rdref.ref());
	int dataSize = rd.Samples();
	short ped = rd.GetPedestal();
	std::vector<short> rawadc(dataSize);
	raw::Uncompress(rd.ADCs(),rawadc,rd.Compression());
	int t0 = hit_peakT[i] - 3*(hit_peakT[i]-hit_tstart[i]);
	if (t0<0) t0 = 0;
	int t1 = hit_peakT[i] + 3*(hit_peakT[i]-hit_tstart[i]);
	if (t1>=dataSize) t1 = dataSize-1;
	hit_pk[i] = -1;
	hit_t[i] = -1;
	for (int j = t0; j<=t1; ++j){
	  if (rawadc[j]-ped>hit_pk[i]){
	    hit_pk[i] = rawadc[j]-ped;
	    hit_t[i] = j;
	  }
	}
	hit_ch[i] = 0;
	hit_fwhh[i] = 0;
	double mean_t = 0;
	double mean_t2 = 0;
	for (int j = t0; j<=t1; ++j){
	  if (rawadc[j]-ped>=0.5*hit_pk[i]){
	    ++hit_fwhh[i];
	  }
	  if (rawadc[j]-ped>=0.1*hit_pk[i]){
	    hit_ch[i] += rawadc[j]-ped;
	    mean_t += j*(rawadc[j]-ped);
	    mean_t2 += j*j*(rawadc[j]-ped);
	  }
	}
	mean_t/=hit_ch[i];
	mean_t2/=hit_ch[i];
	hit_rms[i] = sqrt(mean_t2-mean_t*mean_t);
	if (!evt.isRealData()){
	  //	std::vector<sim::IDE> ides;	
	  //	bt->HitToSimIDEs(hitlist[i], ides);
	  hit_nelec[i] = 0;
	  hit_energy[i] = 0;
	  //	for (size_t j = 0; j<ides.size(); ++j){
	  //	  hit_nelec[i] += ides[j].numElectrons;
	  //	  hit_energy[i] += ides[j].energy;
	  //	}
	  const sim::SimChannel* chan = 0;
	  for(size_t sc = 0; sc < fSimChannels.size(); ++sc){
	    if(fSimChannels[sc]->Channel() == hitlist[i]->Channel()) chan = fSimChannels[sc];
	  }
	  if (chan){
	    const auto & tdcidemap = chan->TDCIDEMap();
	    for(auto mapitr = tdcidemap.begin(); mapitr != tdcidemap.end(); mapitr++){
	      // loop over the vector of IDE objects.
	      
	      const std::vector<sim::IDE> idevec = (*mapitr).second;
	      
	      for(size_t iv = 0; iv < idevec.size(); ++iv){ 

		hit_nelec[i] += idevec[iv].numElectrons;
		hit_energy[i] += idevec[iv].energy;
	      }
	    }
	  }
	}
      //std::cout<<hit_wire[i]<<" "<<hit_peakT[i]<<" "<<hit_ph[i]<<" "<<hit_tend[i]-hit_tstart[i]<<" "<<hit_t[i]<<" "<<hit_pk[i]<<" "<<hit_ch[i]<<" "<<hit_fwhh[i]<<" "<<hit_rms[i]<<" "<<mean_t<<" "<<mean_t2<<std::endl;

      } // if cet::maybe_ref is valid
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
  //fTree->Branch("trigtime",trigtime,"trigtime[16]/I");
  fTree->Branch("nclus",&nclus,"nclus/I");
  fTree->Branch("clustartwire",clustartwire,"clustartwire[nclus]/D");
  fTree->Branch("clustarttick",clustarttick,"clustarttick[nclus]/D");
  fTree->Branch("cluendwire",cluendwire,"cluendwire[nclus]/D");
  fTree->Branch("cluendtick",cluendtick,"cluendtick[nclus]/D");
  fTree->Branch("cluplane",cluplane,"cluplane[nclus]/I");
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
  fTree->Branch("trklen",trklen,"trklen[ntracks_reco]/D");
  fTree->Branch("trkmomrange",trkmomrange,"trkmomrange[ntracks_reco]/D");
  fTree->Branch("trkmommschi2",trkmommschi2,"trkmommschi2[ntracks_reco]/D");
  fTree->Branch("trkmommsllhd",trkmommsllhd,"trkmommsllhd[ntracks_reco]/D");
  fTree->Branch("ntrkhits",ntrkhits,"ntrkhits[ntracks_reco]/I");
  fTree->Branch("trkx",trkx,"trkx[ntracks_reco][1000]/D");
  fTree->Branch("trky",trky,"trky[ntracks_reco][1000]/D");
  fTree->Branch("trkz",trkz,"trkz[ntracks_reco][1000]/D");
  fTree->Branch("trkpitch",trkpitch,"trkpitch[ntracks_reco][2]/D");
  fTree->Branch("trkhits",trkhits,"trkhits[ntracks_reco][2]/I");
  fTree->Branch("trkdedx",trkdedx,"trkdedx[ntracks_reco][2][1000]/D");
  fTree->Branch("trkrr",trkrr,"trkrr[ntracks_reco][2][1000]/D");
  fTree->Branch("trkpitchhit",trkpitchhit,"trkpitchhit[ntracks_reco][2][1000]/D");
  fTree->Branch("trkke",trkke,"trkke[ntracks_reco][2]/D");
  fTree->Branch("trkpida",trkpida,"trkpida[ntracks_reco][2]/D");
  fTree->Branch("nhits",&nhits,"nhits/I");
  fTree->Branch("hit_plane",hit_plane,"hit_plane[nhits]/I");
  fTree->Branch("hit_wire",hit_wire,"hit_wire[nhits]/I");
  fTree->Branch("hit_channel",hit_channel,"hit_channel[nhits]/I");
  fTree->Branch("hit_peakT",hit_peakT,"hit_peakT[nhits]/D");
  fTree->Branch("hit_charge",hit_charge,"hit_charge[nhits]/D");
  fTree->Branch("hit_ph",hit_ph,"hit_ph[nhits]/D");
  fTree->Branch("hit_tstart",hit_tstart,"hit_tstart[nhits]/D");
  fTree->Branch("hit_tend",hit_tend,"hit_tend[nhits]/D");
  fTree->Branch("hit_trkid",hit_trkid,"hit_trkid[nhits]/I");
  fTree->Branch("hit_trkkey",hit_trkkey,"hit_trkkey[nhits]/I");
  fTree->Branch("hit_clukey",hit_clukey,"hit_clukey[nhits]/I");
  fTree->Branch("hit_pk",hit_pk,"hit_pk[nhits]/I");
  fTree->Branch("hit_t",hit_t,"hit_t[nhits]/I");
  fTree->Branch("hit_ch",hit_ch,"hit_ch[nhits]/I");
  fTree->Branch("hit_fwhh",hit_fwhh,"hit_fwhh[nhits]/I");
  fTree->Branch("hit_rms",hit_rms,"hit_rms[nhits]/D");
  fTree->Branch("hit_nelec",hit_nelec,"hit_nelec[nhits]/D");
  fTree->Branch("hit_energy",hit_energy,"hit_energy[nhits]/D");
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
//  for (int i = 0; i < 16; ++i){
//     trigtime[i]=-99999;
//  }
  nclus = -99999;
  for (int i = 0; i < kMaxCluster; ++i){
    clustartwire[i] = -99999;
    clustarttick[i] = -99999;
    cluendwire[i] = -99999;
    cluendtick[i] = -99999;
    cluplane[i] = -99999;
  }
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
    trklen[i] = -99999;
    trkmomrange[i] = -99999;
    trkmommschi2[i] = -99999;
    trkmommsllhd[i] = -99999;
    ntrkhits[i] = -99999;
    for (int j = 0; j<kMaxTrackHits; ++j){
      trkx[i][j] = -99999;
      trky[i][j] = -99999;
      trkz[i][j] = -99999;
    }
    for (int j = 0; j<2; ++j){
      trkpitch[i][j] = -99999;
      trkhits[i][j] = -99999; 
      trkke[i][j] = -99999;
      trkpida[i][j] = -99999;
      for (int k = 0; k<1000; ++k){
	trkdedx[i][j][k] = -99999;
	trkrr[i][j][k] = -99999;
	trkpitchhit[i][j][k] = -99999;
      }
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
    hit_trkkey[i] = -99999;
    hit_clukey[i] = -99999;
    hit_tstart[i] = -99999;
    hit_tend[i] = -99999;
    hit_pk[i] = -99999;
    hit_t[i] = -99999;
    hit_ch[i] = -99999;
    hit_fwhh[i] = -99999;
    hit_rms[i] = -99999;
    hit_nelec[i] = -99999;
    hit_energy[i] = -99999;
  }
  
}

DEFINE_ART_MODULE(bo::AnaTree)
