////////////////////////////////////////////////////////////////////////
// Class:       AnaTreeT1034
// Module Type: analyzer
// File:        AnaTreeT1034_module.cc
//
// Generated at Tue Jul 14 11:21:46 2015 by Roberto Acciarri using artmod
// from cetpkgsupport v1_08_06.
////////////////////////////////////////////////////////////////////////

// Framework includes
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h" 
#include "art/Framework/Principal/Event.h" 
#include "fhiclcpp/ParameterSet.h" 
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Principal/Handle.h" 
#include "art/Persistency/Common/Ptr.h" 
#include "art/Persistency/Common/PtrVector.h" 
#include "art/Framework/Services/Registry/ServiceHandle.h" 
#include "art/Framework/Services/Optional/TFileService.h" 
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "art/Framework/Core/FindOneP.h" 
#include "messagefacility/MessageLogger/MessageLogger.h" 
#include "cetlib/maybe_ref.h"

// LArSoft includes
#include "SimpleTypesAndConstants/geo_types.h"
#include "SimpleTypesAndConstants/RawTypes.h" // raw::ChannelID_t
#include "Geometry/Geometry.h"
#include "Geometry/CryostatGeo.h"
#include "Geometry/TPCGeo.h"
#include "Geometry/PlaneGeo.h"
#include "Geometry/WireGeo.h"
#include "RecoBase/Wire.h"
#include "RecoBase/Hit.h"
#include "RecoBase/Cluster.h"
#include "RecoBase/Track.h"
#include "RecoBase/Vertex.h"
#include "RecoBase/SpacePoint.h"
#include "Utilities/LArProperties.h"
#include "Utilities/DetectorProperties.h"
#include "Utilities/AssociationUtil.h"
#include "RawData/ExternalTrigger.h"
#include "RawData/RawDigit.h"
#include "RawData/raw.h"
#include "RawDataUtilities/TriggerDigitUtility.h" 
#include "MCCheater/BackTracker.h"
#include "Simulation/SimChannel.h"
#include "Filters/ChannelFilter.h"

// ROOT includes
#include "TComplex.h"
#include "TFile.h"
#include "TH2D.h"
#include "TF1.h"
#include "TTree.h"
#include "TTimeStamp.h"

const int kMaxTrack      = 1000;  //maximum number of tracks
const int kMaxHits       = 10000; //maximum number of hits
const int kMaxTrackHits  = 1000;  //maximum number of space points

namespace lariat 
{
   class AnaTreeT1034;
}

class lariat::AnaTreeT1034 : public art::EDAnalyzer 
{
public:
   explicit AnaTreeT1034(fhicl::ParameterSet const & p);
   virtual ~AnaTreeT1034();

   // Required functions.
   void analyze(art::Event const & e) override;

   // Selected optional functions.
   void beginJob();
   void reconfigure(fhicl::ParameterSet const & p);

private:

   void ResetVars();
  
   TTree* fTree;
   //run information
   int run;
   int subrun;
   int event;
   int trig;
   double evttime;
   double efield[3];
   int t0;
//   int trigtime[16];
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
   int    hit_trkid[kMaxHits];
   int    hit_pk[kMaxHits];
   int    hit_t[kMaxHits];
   int    hit_ch[kMaxHits];
   int    hit_fwhh[kMaxHits];
   double hit_rms[kMaxHits];
   double hit_nelec[kMaxHits];
   double hit_energy[kMaxHits];

   std::string fTriggerUtility;     // label for Fragment to digit module
   std::string fHitsModuleLabel;
   std::string fTrackModuleLabel;

};


lariat::AnaTreeT1034::AnaTreeT1034(fhicl::ParameterSet const & pset) : EDAnalyzer(pset)
{
   this->reconfigure(pset);
}

lariat::AnaTreeT1034::~AnaTreeT1034()
{
  // Clean up dynamic memory and other resources here.
}

void lariat::AnaTreeT1034::reconfigure(fhicl::ParameterSet const & pset)
{
   fTriggerUtility       = pset.get< std::string >("TriggerUtility", "FragmentToDigit");
   fHitsModuleLabel      = pset.get< std::string >("HitsModuleLabel");
   fTrackModuleLabel     = pset.get< std::string >("TrackModuleLabel");

   return;
}

void lariat::AnaTreeT1034::analyze(art::Event const & evt)
{
   // Implementation of required member function here.
   ResetVars();

   rdu::TriggerDigitUtility tdu(evt, fTriggerUtility);

   // get services
   art::ServiceHandle<geo::Geometry> geom;
   art::ServiceHandle<util::LArProperties> larprop;
   art::ServiceHandle<util::DetectorProperties> detprop;
   art::ServiceHandle<cheat::BackTracker> bt;
  
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

//  art::Handle< std::vector<raw::ExternalTrigger> > trigListHandle;
//  std::vector<art::Ptr<raw::ExternalTrigger> > triglist;
//  if (evt.getByLabel(fTrigModuleLabel,trigListHandle))
//    art::fill_ptr_vector(triglist, trigListHandle);

//  for (size_t i = 0; i<triglist.size(); ++i){
//    trigtime[i] = triglist[i]->GetTrigTime();
//  }

   std::vector <art::Ptr <recob::Hit>   > hitlist;
   std::vector <art::Ptr <recob::Track> > tracklist;

   art::FindManyP<recob::Hit>   fmh (tdu.EventTriggersPtr(), evt, fHitsModuleLabel);
   art::FindManyP<recob::Track> fmt (tdu.EventTriggersPtr(), evt, fTrackModuleLabel);
  
   for(size_t t = 0; t < tdu.NTriggers(); ++t)        // Loop over triggers                          
   {
      // Skip trigger if empty
      art::PtrVector<raw::RawDigit> rdvec = tdu.TriggerRawDigitsPtr(t);

      LOG_VERBATIM("AnaTreeT1034") << " ";       
      LOG_VERBATIM("AnaTreeT1034") << " ";
      LOG_VERBATIM("AnaTreeT1034") << "Trigger Number: " << t << "   Raw Digit vector size: "<< rdvec.size();
      if(!rdvec.size()){mf::LogInfo("AnaTreeT1034") << " Raw Digit vector is empty. Skipping the trigger"; continue;}                              
      LOG_VERBATIM("AnaTreeT1034") << " ";
      
      trig = t; 

      // get input Hits and Tracks objects per trigger.
      hitlist.clear();
      tracklist.clear();
      hitlist   = fmh.at(t);
      tracklist = fmt.at(t);
/*
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
  art::FindMany<recob::Track>       fmtk(hitListHandle,   evt, fTrackModuleLabel);
*/

//      art::FindOne<raw::RawDigit>       ford (hitlist,   evt, fHitsModuleLabel);
      art::FindManyP<recob::SpacePoint> fmsp (tracklist, evt, fTrackModuleLabel);
      art::FindMany<recob::Track>       fmtk (hitlist,   evt, fTrackModuleLabel);

      std::vector<const sim::SimChannel*> fSimChannels;
      try
      {
         evt.getView("largeant", fSimChannels);
      }
      catch (art::Exception const&e)
      {
      }

      //track information
      ntracks_reco=tracklist.size();
      double larStart[3];
      double larEnd[3];
      std::vector<double> trackStart;
      std::vector<double> trackEnd;
      for(size_t i=0; i<tracklist.size();++i)  // loop over tracks
      {
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
         for (size_t j = 0; j<spts.size(); ++j)
         {
            trkx[i][j] = spts[j]->XYZ()[0];
            trky[i][j] = spts[j]->XYZ()[1];
            trkz[i][j] = spts[j]->XYZ()[2];
         }
         for (int j = 0; j<2; ++j)
         {
            try
            {
	       if (j==0)      trkpitch[i][j] = tracklist[i]->PitchInView(geo::kU);
	       else if (j==1) trkpitch[i][j] = tracklist[i]->PitchInView(geo::kV);
            }
            catch( cet::exception &e)
            {
	       mf::LogWarning("AnaTreeT1034")<<"caught exeption "<<e<<"\n setting pitch to 0";
	       trkpitch[i][j] = 0;
            }
         }
      } // loop over tracks

      //hits information
      nhits = hitlist.size();
      for (size_t i = 0; i<hitlist.size(); ++i)  // loop over hits
      {
 //        cet::maybe_ref<raw::RawDigit const> rdref(ford.at(i));
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
         if (fmtk.at(i).size()!=0)
         {
            hit_trkid[i] = fmtk.at(i)[0]->ID();
         }
         if (hit_plane[i]==1) //collection plane
         {
/*            if( rdref.isValid() )
            {
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
	       for (int j = t0; j<=t1; ++j)
               {
	          if (rawadc[j]-ped>hit_pk[i])
                  {
	             hit_pk[i] = rawadc[j]-ped;
	             hit_t[i] = j;
	          }
	       }
	       hit_ch[i] = 0;
	       hit_fwhh[i] = 0;
	       double mean_t = 0;
	       double mean_t2 = 0;
	       for (int j = t0; j<=t1; ++j)
               {
	          if (rawadc[j]-ped>=0.5*hit_pk[i])
                  {
	             ++hit_fwhh[i];
	          }
	          if (rawadc[j]-ped>=0.1*hit_pk[i])
                  {
	             hit_ch[i] += rawadc[j]-ped;
	             mean_t += j*(rawadc[j]-ped);
	             mean_t2 += j*j*(rawadc[j]-ped);
	          }
	       }
	       mean_t/= hit_ch[i];
	       mean_t2/= hit_ch[i];
	       hit_rms[i] = sqrt(mean_t2-mean_t*mean_t);
	       if (!evt.isRealData())
               {
	          hit_nelec[i] = 0;
	          hit_energy[i] = 0;
 	          const sim::SimChannel* chan = 0;
	          for(size_t sc = 0; sc < fSimChannels.size(); ++sc)
                  {
	             if(fSimChannels[sc]->Channel() == hitlist[i]->Channel()) chan = fSimChannels[sc];
	          }
	          if (chan)
                  {
	             const std::map<unsigned short, std::vector<sim::IDE> >& tdcidemap = chan->TDCIDEMap();
	             for(auto mapitr = tdcidemap.begin(); mapitr != tdcidemap.end(); mapitr++) // loop over the map
                     {
	      	        const std::vector<sim::IDE> idevec = (*mapitr).second;	      
	                for(size_t iv = 0; iv < idevec.size(); ++iv)   // loop over the vector of IDE objects.
                        { 
		           hit_nelec[i] += idevec[iv].numElectrons;
		           hit_energy[i] += idevec[iv].energy;
	                }  // loop over the vector of IDE objects.
	             }  // loop over the map
                  }  // if(chan)
	       } //if (!evt.isRealData())
               //std::cout<<hit_wire[i]<<" "<<hit_peakT[i]<<" "<<hit_ph[i]<<" "<<hit_tend[i]-hit_tstart[i]<<" "<<hit_t[i]<<" "<<hit_pk[i]<<" "<<hit_ch[i]<<" "<<hit_fwhh[i]<<" "<<hit_rms[i]<<" "<<mean_t<<" "<<mean_t2<<std::endl;

            } // if rdref.isValid()*/
         }  //collection plane
      }    // loop over hits

      fTree->Fill();

   } //loop over triggers
}

void lariat::AnaTreeT1034::beginJob()
{
   // Implementation of optional member function here.
   art::ServiceHandle<art::TFileService> tfs;
   fTree = tfs->make<TTree>("anatree","analysis tree");
   fTree->Branch("run",&run,"run/I");
   fTree->Branch("subrun",&subrun,"subrun/I");
   fTree->Branch("event",&event,"event/I");
   fTree->Branch("trig",&trig,"trig/I");
   fTree->Branch("evttime",&evttime,"evttime/D");
   fTree->Branch("efield",efield,"efield[3]/D");
   fTree->Branch("t0",&t0,"t0/I");
//   fTree->Branch("trigtime",trigtime,"trigtime[16]/I");
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
   fTree->Branch("hit_trkid",hit_trkid,"hit_trkid[nhits]/I");
   fTree->Branch("hit_pk",hit_pk,"hit_pk[nhits]/I");
   fTree->Branch("hit_t",hit_t,"hit_t[nhits]/I");
   fTree->Branch("hit_ch",hit_ch,"hit_ch[nhits]/I");
   fTree->Branch("hit_fwhh",hit_fwhh,"hit_fwhh[nhits]/I");
   fTree->Branch("hit_rms",hit_rms,"hit_rms[nhits]/D");
   fTree->Branch("hit_nelec",hit_nelec,"hit_nelec[nhits]/D");
   fTree->Branch("hit_energy",hit_energy,"hit_energy[nhits]/D");
}

void lariat::AnaTreeT1034::ResetVars()
{

   run = -99999;
   subrun = -99999;
   event = -99999;
   trig = -99999;
   evttime = -99999;
   for (int i = 0; i<3; ++i)
   {
      efield[i] = -99999;
   }
   t0 = -99999;
   ntracks_reco = -99999;
   for (int i = 0; i < kMaxTrack; ++i)
   {
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
      for (int j = 0; j<kMaxTrackHits; ++j)
      {
         trkx[i][j] = -99999;
         trky[i][j] = -99999;
         trkz[i][j] = -99999;
      }
      for (int j = 0; j<3; ++j)
      {
         trkpitch[i][j] = -99999;
      }
   }  
   nhits = -99999;
   for (int i = 0; i<kMaxHits; ++i) 
   {
      hit_plane[i] = -99999;
      hit_wire[i] = -99999;
      hit_channel[i] = -99999;
      hit_peakT[i] = -99999;
      hit_charge[i] = -99999;
      hit_ph[i] = -99999;
      hit_trkid[i] = -99999;
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

DEFINE_ART_MODULE(lariat::AnaTreeT1034)
