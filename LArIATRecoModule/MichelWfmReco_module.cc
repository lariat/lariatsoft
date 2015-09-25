////////////////////////////////////////////////////////////////////////
// Class:       MichelWfmReco
// Module Type: producer
// File:        MichelWfmReco_module.cc
//
// This module is used to perform some ID and reconstruction of PMT
// waveforms from stopping/decaying muons (primarily in the Michel 
// trigger sample).
//
// Eventually it may be used to create a new data product related to 
// Michel events, but for now, it does not add to the data file.
//
// Output histograms include:
//  - # hits found in each waveform
//  - amplitude of Michel candidate pulses
//  - integrated charge of Michel candidates (100ns window)
//  - time difference between 1st and 2nd pulse
//    when exactly two pulses are found
//  - information for reconstructed tracks associated with
//    Michel candidate events
//
// Authors: William Foreman, wforeman@uchicago.edu
//
// Generated at Wed Jul 15 13:09:43 2015 by William Foreman using artmod
// from cetpkgsupport v1_08_06.
//
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
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"

//C++ Includes
#include <iostream>
#include <fstream>
#include <vector>
#include <memory>
#include <utility>

//ROOT Includes
#include <TH1F.h>
#include <TH1I.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TTree.h>

// LArSoft Includes
#include "Utilities/AssociationUtil.h"
#include "RawData/TriggerData.h"
#include "RecoBase/Track.h"
#include "AnalysisBase/Calorimetry.h"

//LAriatSoft Includes
#include "RawDataUtilities/TriggerDigitUtility.h"
#include "LArIATRecoAlg/OpHitBuilderAlg.h"
#include "LArIATRecoAlg/TriggerFilterAlg.h"
#include "Utilities/DatabaseUtilityT1034.h"


class MichelWfmReco;

class MichelWfmReco : public art::EDProducer {
public:
  explicit MichelWfmReco(fhicl::ParameterSet const & p);
  // The destructor generated by the compiler is fine for classes
  // without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  MichelWfmReco(MichelWfmReco const &) = delete;
  MichelWfmReco(MichelWfmReco &&) = delete;
  MichelWfmReco & operator = (MichelWfmReco const &) = delete;
  MichelWfmReco & operator = (MichelWfmReco &&) = delete;

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
  
  // Custom functions
  bool IsPointInFiducialVolume(TVector3);

private:

  // Name of the module producing the triggers
  std::string fTriggerUtility;

  // Tunable parameters from fcl
  bool            bUseTriggerFilter;
  bool            bVerbose;
  std::string     fDAQModule;
  std::string     fTrackModule;
  std::string     fTrackCalModule;
  std::string     fInstanceName;
  short           fBaselineWindowLength;
  double          fFiducialMargin_X;
  double          fFiducialMargin_Y;
  double          fFiducialMargin_Z; 

  // Alg objects
  OpHitBuilderAlg   fOpHitBuilderAlg; 
  TriggerFilterAlg  fTrigFiltAlg;

  // Variables/vectors
  bool                  GotETL;
  bool                  flag;
  std::vector<short>    ETL_waveform;
  short                 PostPercentMark;
  int                   NSamples;
  std::vector<TVector3> TrackVertex;
  std::vector<TVector3> TrackEnd;
  int                   MuTrackIndex;
  TVector3              MuCandidateTrackVertex;
  TVector3              MuCandidateTrackEnd;
  TVector3              region_centerpoint;
  double                region_radius;
  std::vector<double>   TrackEnergy;
  

  // Histograms
  TH1I* h_FilterStages;
  TH1I* h_NumOpHits;
  TH2I* h_NumOpHits_vs_NumTracks;
  TH1I* h_NumOpHits_beam;
  TH1I* h_NumOpHits_offbeam;
  TH1F* h_DeltaTime_chargeCut;
  TH1F* h_DeltaTime_stoppingMu;

  TH1F* h_Amplitude;
  TH1F* h_Charge100ns;
  TH1F* h_Charge100ns_stoppingMu;
  TH1F* h_Charge100ns_region;
  TH1F* h_TrackVertex_x;
  TH1F* h_TrackVertex_y;
  TH1F* h_TrackVertex_z;
  //TH3F* h_TrackVertex;
  TH1F* h_TrackEnd_x;
  TH1F* h_TrackEnd_y;
  TH1F* h_TrackEnd_z;
  //TH3F* h_TrackEnd;
  TH2F* h_TrackNode_zx;
  TH2F* h_TrackNode_zy;
  TH2I* h_InFiducialVolume;
  TH1F* h_SecondTrackOffset;

  // Variables for Tree
  Int_t       iEvent = -1;
  Double_t    WaveformBaseline;
  Double_t    WaveformBaselineRMS;
  Double_t    Timestamp;
  Double_t    DeltaTime;
  Int_t       NumHits;
  Double_t    Amplitude;
  Double_t    Charge_100ns; 
  Int_t       NumTracks;
  Int_t       IsSingleStoppingTrack;
  Double_t    StoppingTrackZenithAngle;
  Double_t    MuTrackVertex_x;
  Double_t    MuTrackVertex_y;
  Double_t    MuTrackVertex_z;
  Double_t    MuTrackEnd_x;
  Double_t    MuTrackEnd_y;
  Double_t    MuTrackEnd_z;
  Double_t    MuTrackLength;
  Double_t    MuTrackEnergy;
  
  // TTree info
  TTree* MichelDataTree;
  TBranch* b_iEvent;
  TBranch* b_Timestamp;
  TBranch* b_WaveformBaseline;
  TBranch* b_WaveformBaselineRMS;
  TBranch* b_NumHits;
  TBranch* b_DeltaTime;
  TBranch* b_Amplitude;
  TBranch* b_Charge_100ns;
  TBranch* b_NumTracks;
  TBranch* b_IsSingleStoppingTrack;
  TBranch* b_StoppingTrackZenithAngle;
  TBranch* b_MuTrackVertex_x;
  TBranch* b_MuTrackVertex_y;
  TBranch* b_MuTrackVertex_z;
  TBranch* b_MuTrackEnd_x;
  TBranch* b_MuTrackEnd_y;
  TBranch* b_MuTrackEnd_z;
  TBranch* b_MuTrackLength;
  TBranch* b_MuTrackEnergy;

};



MichelWfmReco::MichelWfmReco(fhicl::ParameterSet const & p)
: fOpHitBuilderAlg(p), fTrigFiltAlg(p)
{

  // Configures the ROOT histograms
  this->reconfigure(p);
  
  // Produces the LArSoft object to be outputted
  // (none for now!) 
 
  // TO DO: make this a fcl parameter 
  region_centerpoint.SetXYZ(23.5,0.,45.);
  region_radius = 5.;
  
}


void MichelWfmReco::produce(art::Event & e)
{
  iEvent++;
  h_FilterStages->Fill(0);
  std::cout<<"---------- Starting event "<<iEvent<<" -----------------\n";
  

  // Initialize variables and vectors
  Timestamp         = -9.; 
  NumHits           = -9;
  NumTracks         = -9;
  IsSingleStoppingTrack   = 0;
  StoppingTrackZenithAngle = -9;
  MuTrackVertex_x   = -99.;
  MuTrackVertex_y   = -99.;
  MuTrackVertex_z   = -99.;
  MuTrackEnd_x      = -99.;
  MuTrackEnd_y      = -99.;
  MuTrackEnd_z      = -99.;
  MuTrackLength     = -9.;
  MuTrackEnergy     = -9.;
  DeltaTime         = -99.;
  Amplitude         = -99.;
  Charge_100ns      = -99.;
  WaveformBaseline  = -99.;
  WaveformBaselineRMS = -9.;
  TrackVertex.clear();
  TrackEnd.clear();
  TrackEnergy.clear();
  
  
  // Filter for the MICHEL trigger pattern (note that for some runs, particularly for 
  // those with the optimized Michel trigger setup, unfortunately the trigger inputs 
  // were not being saved and thus the filter won't work).
  /*
  bool isMichel = 1; 
  if (bUseTriggerFilter){
    std::string myFilter = "+MICHEL";
    isMichel = fTrigFiltAlg.doesTriggerPassFilter( thisTrigger, myFilter ); 
  } else {
    isMichel = 1;
  }
  if ( !isMichel ) return;
  */


  // Define the vector of associations to be saved
  // (none yet!)


  // Get the OpDetPulses; skip event if empty
  art::Handle< std::vector< raw::OpDetPulse >> WaveformHandle;
  e.getByLabel(fDAQModule,fInstanceName,WaveformHandle);
  if( (int)WaveformHandle->size() == 0 ) {
    std::cout<<"No optical detector data found -- skipping the event\n";
    return;
  } else {
    
    // If not empty, store ETL waveform
    GotETL = false;
    for( int ipulse = 0; ipulse < (int)WaveformHandle->size(); ipulse++){
      art::Ptr< raw::OpDetPulse > ThePulsePtr(WaveformHandle,ipulse);
      raw::OpDetPulse ThePulse = *ThePulsePtr;
      if( ThePulse.OpChannel() == 1) {
        ETL_waveform = ThePulse.Waveform();
        NSamples = ETL_waveform.size();
        Timestamp = (double(ThePulse.PMTFrame())*8.)/1.0e09;
        PostPercentMark = short(ThePulse.FirstSample); 

        std::cout<<"ETL pulse recorded: Nsamples = "<<NSamples<<std::endl;
        std::cout<<"   OpDetPulse PMTFrame = "<<ThePulse.PMTFrame()<<
          "  ("Timestamp<<" sec\n";
        std::cout<<"   OpDetPulse FirstSample = "<<ThePulse.FirstSample()<<std::endl;
        GotETL = true;


      }
    }
  }


  // If we somehow didn't get the ETL, skip the event
  if( !GotETL ) return;
  h_FilterStages->Fill(1);
  
  
  // Get the tracks and their associated energy
  art::Handle< std::vector< recob::Track >> TrackHandle;
  e.getByLabel(fTrackModule,TrackHandle);
  art::FindManyP<anab::Calorimetry> fmcal(TrackHandle, e, fTrackCalModule);
  NumTracks = (int)TrackHandle->size();
  std::cout<<"Number of tracks: "<<NumTracks<<"\n";
  
  // Loop through the track list and store their properties
  for( int track_index = 0; track_index < NumTracks; track_index++){

    // Get the recob::Track object and record its endpoint/vertex
    art::Ptr< recob::Track > TheTrackPtr(TrackHandle,track_index);
    recob::Track TheTrack = *TheTrackPtr;
    TrackEnd.push_back(TheTrack.End());
    TrackVertex.push_back(TheTrack.Vertex());
    std::cout<<"  Track vertex: ("
    << TheTrack.Vertex().X() <<"," 
    << TheTrack.Vertex().Y() << "," 
    << TheTrack.Vertex().Z() << ")\n"; 
    std::cout<<"  Track endpoint: ("
    << TheTrack.End().X() <<"," 
    << TheTrack.End().Y() << "," 
    << TheTrack.End().Z() << ")\n"; 
    std::cout<<"  Track length: "<<(TheTrack.Vertex() - TheTrack.End()).Mag()<<std::endl;
    
    std::cout<<"Vertex in fiducial vol?  "<<IsPointInFiducialVolume(TheTrack.Vertex())<<std::endl;
    
    h_InFiducialVolume->Fill(
      IsPointInFiducialVolume(TheTrack.Vertex()),
      IsPointInFiducialVolume(TheTrack.End()));

    h_TrackNode_zx->Fill((float)TheTrack.Vertex().Z(),(float)TheTrack.Vertex().X());
    h_TrackNode_zx->Fill((float)TheTrack.End().Z(),(float)TheTrack.End().X());
    h_TrackNode_zy->Fill((float)TheTrack.Vertex().Z(),(float)TheTrack.Vertex().Y());
    h_TrackNode_zy->Fill((float)TheTrack.End().Z(),(float)TheTrack.End().Y());
    

    // Get measured track energy (if possible)
    if( fmcal.isValid() ){
      std::vector<art::Ptr<anab::Calorimetry> > calos = fmcal.at(track_index);
      // calos[0] is from induction plane
      // calos[1] is from collection plane
      TrackEnergy.push_back(calos[1]->KineticEnergy());
      std::cout<<"  Energy = "<<calos[1]->KineticEnergy()<<"\n"; 
    } else {
      std::cout<<"  Energy = undefined\n";
    }
  
  }  


  // --------------------------------------
  // Higher-level track filtering:
  
  // Require no more than two tracks in the event:
  if( NumTracks <= 2 ) {
    
    flag = false;

    // Cycle through tracks
    for( int i=0; i<NumTracks; i++){
        
      // Since we don't want to blindly trust that the endpoint and 
      // vertex have been assigned correctly, let's say whichever 
      // has a higher Y-value is the vertex (pretty reasonable 
      // assumption for cosmic muons).
      TVector3 vertex;
      TVector3 end;
      if( TrackVertex[i].Y() > TrackEnd[i].Y()  ) { 
        vertex  = TrackVertex[i]; 
        end     = TrackEnd[i];
      } else {
        vertex  = TrackEnd[i];
        end     = TrackVertex[i];
      }
    
      // Fill histograms
      h_TrackEnd_x->Fill((float)end.X());
      h_TrackEnd_y->Fill((float)end.Y());
      h_TrackEnd_z->Fill((float)end.Z());
      //h_TrackEnd  ->Fill(..,..,..); 
      h_TrackVertex_x->Fill((float)vertex.X());
      h_TrackVertex_y->Fill((float)vertex.Y());
      h_TrackVertex_z->Fill((float)vertex.Z());
      //h_TrackVertex  ->Fill(..,..,..);

      // Check that vertex is outside fiducial volume and endpoint is 
      // inside fiducial volume.
      if( (!flag) && (IsPointInFiducialVolume(end)) && (!IsPointInFiducialVolume(vertex))) {
        
        // If there was already a stopping track, this disqualifies 
        // the event.
        if( IsSingleStoppingTrack ){
          IsSingleStoppingTrack = 0;
          flag = true;
          MuTrackVertex_x   = -99.;
          MuTrackVertex_y   = -99.;
          MuTrackVertex_z   = -99.;
          MuTrackEnd_x      = -99.;
          MuTrackEnd_y      = -99.;
          MuTrackEnd_z      = -99.;
          MuTrackLength     = -9.;
          MuTrackEnergy     = -9.;
          StoppingTrackZenithAngle = -9.;
        } else { 
          MuTrackIndex = i;
          MuCandidateTrackVertex  = vertex;
          MuCandidateTrackEnd     = end;
          IsSingleStoppingTrack = 1;
          MuTrackVertex_x = vertex.X();
          MuTrackVertex_y = vertex.Y();
          MuTrackVertex_z = vertex.Z();
          MuTrackEnd_x    = end.X();
          MuTrackEnd_y    = end.Y();
          MuTrackEnd_z    = end.Z();
          MuTrackLength   = (vertex-end).Mag();
          MuTrackEnergy   = TrackEnergy[i];

          TVector3 vert(0.,1.,0.);
          TVector3 tmp = (vertex-end);
          StoppingTrackZenithAngle = tmp.Angle(vert);
        }
      
      }
    } 
  }

  // If there was a stopping "muon" and a second track, check
  // to see if either of its endpoints is close to muon endpoint.
  if( IsSingleStoppingTrack && NumTracks == 2 ) {
    for( int i=0; i<NumTracks; i++){
      if( i != MuTrackIndex ){
        double proximity = std::min( 
          (TrackVertex[i] - MuCandidateTrackEnd).Mag(),
          (TrackEnd[i]   - MuCandidateTrackEnd).Mag());
        h_SecondTrackOffset->Fill(proximity);
      }
    }  
  }



  // ---------------------------------------
  // Analysis:

  // Perform hit-finding/filtering
  std::vector<short> hit_times = fOpHitBuilderAlg.GetHits(ETL_waveform);
  NumHits = hit_times.size();
  std::cout << "We found "<<NumHits<<" hits\n";
  for (int i=0; i<NumHits; i++) std::cout<<"   "<<i<<"    t = "<<hit_times[i]<<"\n";

  // Save the baseline/RMS
  std::vector<double> tmp = fOpHitBuilderAlg.GetBaselineAndRMS(ETL_waveform,0,fBaselineWindowLength);
  WaveformBaseline = tmp[0];
  WaveformBaselineRMS = tmp[1];
  std::cout<<"Baseline: "<<WaveformBaseline<<"  RMS: "<<WaveformBaselineRMS<<"\n";

  // Fill NumOpHits histo (beam vs. offbeam)
  h_NumOpHits ->Fill(NumHits);
  h_NumOpHits_vs_NumTracks->Fill(NumHits,NumTracks);
  if(Timestamp < 5.2)   h_NumOpHits_beam    ->Fill(NumHits);
  if(Timestamp >= 5.3)  h_NumOpHits_offbeam ->Fill(NumHits);
  if(Timestamp >= 5.3)  h_FilterStages      ->Fill(2);
        
  // Event quality control:
  // Require 2 hits (one before PostPercent, one within 1% of PostPercent) 
  if( (NumHits == 2) && (hit_times[0] < PostPercentMark) && ( abs(hit_times[1]-PostPercentMark) <= 0.01*NSamples) && (Timestamp >= 5.3)){
    std::cout << "    passes Michel cut \n";
    h_FilterStages->Fill(3);

    // Measure time difference
    DeltaTime = hit_times[1] - hit_times[0];
        
    // Integral/amplitude of Michel candidate pulse
    std::vector<double> hit_info = fOpHitBuilderAlg.IntegrateHit(ETL_waveform, hit_times[1]);
    Amplitude     = hit_info[0];
    Charge_100ns  = hit_info[1];
    
    h_Amplitude->Fill(Amplitude);
    h_Charge100ns->Fill(Charge_100ns); 
    if( Charge_100ns > 1600. ) h_DeltaTime_chargeCut->Fill(DeltaTime);
  
    if( IsSingleStoppingTrack ) {
      h_FilterStages->Fill(4);
      h_Charge100ns_stoppingMu->Fill(Charge_100ns);
      h_DeltaTime_stoppingMu->Fill(DeltaTime);

      
      if( (MuCandidateTrackEnd-region_centerpoint).Mag() <= region_radius) h_Charge100ns_region->Fill(Charge_100ns);
    }
  }
        
  MichelDataTree->Fill();
    
}


void MichelWfmReco::beginJob()
{
  // Opens up the file service to read information from the ROOT file input
  art::ServiceHandle<art::TFileService> tfs;

  MichelDataTree        = tfs->make<TTree>("MichelDataTree","MichelDataTree");
  b_iEvent              = MichelDataTree->Branch("iEvent",&iEvent,"iEvent/I");
  b_Timestamp           = MichelDataTree->Branch("Timestamp (s)",&Timestamp,"Timestamp/D");
  b_WaveformBaseline    = MichelDataTree->Branch("WaveformBaseline (ADC)",&WaveformBaseline,"WaveformBaseline/D");
  b_WaveformBaselineRMS = MichelDataTree->Branch("WaveformBaselineRMS (ADC)",&WaveformBaselineRMS,"WaveformBaselineRMS/D");
  b_NumHits             = MichelDataTree->Branch("NumHits",&NumHits,"NumHits/I");
  b_DeltaTime           = MichelDataTree->Branch("DeltaTime (ns)",&DeltaTime,"DeltaTime/D");
  b_Amplitude           = MichelDataTree->Branch("Amplitude (mv)",&Amplitude,"Amplitude/D");
  b_Charge_100ns        = MichelDataTree->Branch("Charge_100ns (ADC)",&Charge_100ns,"Charge_100ns/D");
  b_NumTracks           = MichelDataTree->Branch("NumTracks",&NumTracks,"NumTracks/I");
  b_IsSingleStoppingTrack     = MichelDataTree->Branch("IsSingleStoppingTrack",&IsSingleStoppingTrack,"IsStoppingTrack/I");
  b_StoppingTrackZenithAngle = MichelDataTree->Branch("StoppingTrackZenithAngle (rad)",&StoppingTrackZenithAngle,"StoppingTrackZenithAngle/D");
  b_MuTrackVertex_x     = MichelDataTree->Branch("MuTrackVertex_x",&MuTrackVertex_x,"MuTrackVertex_x/D");
  b_MuTrackVertex_y     = MichelDataTree->Branch("MuTrackVertex_y",&MuTrackVertex_y,"MuTrackVertex_y/D");
  b_MuTrackVertex_z     = MichelDataTree->Branch("MuTrackVertex_z",&MuTrackVertex_z,"MuTrackVertex_z/D");
  b_MuTrackEnd_x        = MichelDataTree->Branch("MuTrackEnd_x",&MuTrackEnd_x,"MuTrackEnd_x/D");
  b_MuTrackEnd_y        = MichelDataTree->Branch("MuTrackEnd_y",&MuTrackEnd_y,"MuTrackEnd_y/D");
  b_MuTrackEnd_z        = MichelDataTree->Branch("MuTrackEnd_z",&MuTrackEnd_z,"MuTrackEnd_z/D");
  b_MuTrackLength        = MichelDataTree->Branch("MuTrackLength (cm)",&MuTrackLength,"MuTrackLength/D");
  b_MuTrackEnergy        = MichelDataTree->Branch("MuTrackEnergy (MeV)",&MuTrackEnergy,"MuTrackEnergy/D");
  
  h_FilterStages        = tfs->make<TH1I>("FilterStages", "Sample filtering;filter stage;Events remaining",5,0,5);
  h_FilterStages        ->SetMinimum(0);
  h_NumOpHits           = tfs->make<TH1I>("OpHitsPerEvent", "Optical hits per event", 10, 0, 10);
  h_NumOpHits           ->GetXaxis()->SetTitle("Num hits");
  h_NumOpHits           ->GetYaxis()->SetTitle("Counts");
  h_NumOpHits_beam      = tfs->make<TH1I>("OpHitsPerEvent_beam", "Optical hits per event (beam)",  10, 0, 10);
  h_NumOpHits_beam      ->GetXaxis()->SetTitle("Num hits");
  h_NumOpHits_beam      ->GetYaxis()->SetTitle("Counts");
  h_NumOpHits_offbeam   = tfs->make<TH1I>("OpHitsPerEvent_offbeam", "Optical hits per event (off-beam)", 10, 0, 10);
  h_NumOpHits_offbeam   ->GetXaxis()->SetTitle("Num hits");
  h_NumOpHits_offbeam   ->GetYaxis()->SetTitle("Counts");
  h_NumOpHits_vs_NumTracks = tfs->make<TH2I>("NumOpHits_vs_NumTracks", "Optical hits vs. number of tracks in event", 10, 0,10, 30,0,30);
  h_NumOpHits_vs_NumTracks ->GetXaxis()->SetTitle("Num optical hits");
  h_NumOpHits_vs_NumTracks ->GetYaxis()->SetTitle("Num tracks");
  h_NumOpHits_vs_NumTracks ->SetOption("colz");

  h_DeltaTime_chargeCut = tfs->make<TH1F>("DeltaTime_chargeCutoff", "#Delta t (integral > 1600 ADC)", 700,0., 7000.);
  h_DeltaTime_chargeCut ->GetXaxis()->SetTitle("ns");
  h_DeltaTime_chargeCut ->GetYaxis()->SetTitle("Counts");
  
  h_DeltaTime_stoppingMu = tfs->make<TH1F>("DeltaTime_stoppingMu", "#Delta t (stopping track)", 700,0., 7000.);
  h_DeltaTime_stoppingMu ->GetXaxis()->SetTitle("ns");
  h_DeltaTime_stoppingMu ->GetYaxis()->SetTitle("Counts");

  h_Amplitude           = tfs->make<TH1F>("Amplitude", "Amplitude", 500,   0., 100.);
  h_Amplitude           ->GetXaxis()->SetTitle("Amplitude of Michel PMT pulse [mv]");
  h_Amplitude           ->GetYaxis()->SetTitle("Counts");
  
  h_Charge100ns          = tfs->make<TH1F>("Charge100ns", "Prompt light integral (100ns)",  240,  0., 12000.);
  h_Charge100ns          ->GetXaxis()->SetTitle("Integrated prompt light, 100ns [ADC]");
  h_Charge100ns          ->GetYaxis()->SetTitle("Counts");

  h_Charge100ns_stoppingMu = tfs->make<TH1F>("Charge100ns_stoppingMu","Prompt light integral (100ns) in events with a stopping track",240,0.,12000.);
  h_Charge100ns_stoppingMu ->GetXaxis()->SetTitle("Integrated prompt light, 100ns [ADC]");
  h_Charge100ns_stoppingMu ->GetYaxis()->SetTitle("Counts");

  h_Charge100ns_region    = tfs->make<TH1F>("Charge100ns_region","Prompt light integral (100ns), limited region",240,0.,12000.);
  h_Charge100ns_region    ->GetXaxis()->SetTitle("Integrated prompt light, 100ns [ADC]");
  h_Charge100ns_region    ->GetYaxis()->SetTitle("Counts");

  h_TrackNode_zx          = tfs->make<TH2F>("TrackNode_zx","TrackNode_zx",500,0.,90.,500,0.,47.);
  h_TrackNode_zx          ->GetXaxis()->SetTitle("z [cm]");
  h_TrackNode_zx          ->GetYaxis()->SetTitle("x [cm]");
  h_TrackNode_zx          ->SetMarkerStyle(7);

  h_TrackNode_zy          = tfs->make<TH2F>("TrackNode_zy","TrackNode_zy",500,0.,90.,500,-20.,20.);
  h_TrackNode_zy          ->GetXaxis()->SetTitle("z [cm]");
  h_TrackNode_zy          ->GetYaxis()->SetTitle("y [cm]");
  h_TrackNode_zy          ->SetMarkerStyle(7);

  h_TrackEnd_x            = tfs->make<TH1F>("TrackEnd_x","TrackEnd_x",100,-10.,60.);
  h_TrackEnd_y            = tfs->make<TH1F>("TrackEnd_y","TrackEnd_y",100,-25.,25.);
  h_TrackEnd_z            = tfs->make<TH1F>("TrackEnd_z","TrackEnd_z",100,-10.,100.);
  h_TrackVertex_x         = tfs->make<TH1F>("TrackVertex_x","TrackVertex_x",100,-10.,60.);
  h_TrackVertex_y         = tfs->make<TH1F>("TrackVertex_y","TrackVertex_y",100,-25.,25.);
  h_TrackVertex_z         = tfs->make<TH1F>("TrackVertex_z","TrackVertex_z",100,-10.,100.);
  //h_TrackEnd              = tfs->make<TH3F>("TrackEnd","TrackEnd",  1000,0.,47.,  1000, -20.,20.,  1000, 0.,90.);
  //h_TrackVertex           = tfs->make<TH3F>("TrackVertex","TrackVertex",  1000,0.,47.,  1000, -20.,20.,  1000, 0.,90.);
  
  h_InFiducialVolume      = tfs->make<TH2I>("InFiducialVolume",";Vertex in fiducial volume;End in fiducial volume",2,0.,2.,2,0.,2.);
  h_InFiducialVolume      ->SetOption("colz");
  
  h_SecondTrackOffset     = tfs->make<TH1F>("SecondTrackOffset","Primary / secondary track endpoint offset;cm",100,0.,30.);
}

void MichelWfmReco::beginRun(art::Run & r)
{
  fTrigFiltAlg.loadXMLDatabaseTable( r.run() ); 
}

void MichelWfmReco::beginSubRun(art::SubRun & sr)
{
}

void MichelWfmReco::endJob()
{
}

void MichelWfmReco::endRun(art::Run & r)
{
}

void MichelWfmReco::endSubRun(art::SubRun & sr)
{
}

void MichelWfmReco::reconfigure(fhicl::ParameterSet const & p)
{
  // Pass name of TriggerUtility
  fTriggerUtility         = p.get< std::string >("TriggerUtility","FragmentToDigit");
  bUseTriggerFilter       = p.get< bool >("UseTriggerFilter","false");
  bVerbose                = p.get< bool >("Verbosity","false");
  fDAQModule              = p.get< std::string >("DAQModule","daq");
  fTrackModule            = p.get< std::string >("TrackModule","pmtrack");
  fTrackCalModule         = p.get< std::string >("TrackCalModule","calo");
  fInstanceName           = p.get< std::string >("InstanceName","");
  fBaselineWindowLength   = p.get< short >("BaselineWindowLength",1000);
  fFiducialMargin_X       = p.get< double>("FiducialMargin_X",5.);
  fFiducialMargin_Y       = p.get< double>("FiducialMargin_Y",4.);
  fFiducialMargin_Z       = p.get< double>("FiducialMargin_Z",5.);

}

void MichelWfmReco::respondToCloseInputFile(art::FileBlock const & fb)
{
}

void MichelWfmReco::respondToCloseOutputFiles(art::FileBlock const & fb)
{
}

void MichelWfmReco::respondToOpenInputFile(art::FileBlock const & fb)
{
}

void MichelWfmReco::respondToOpenOutputFiles(art::FileBlock const & fb)
{
}

// Function for determining if a point is inside or outside
// predefined fiducial volume
bool MichelWfmReco::IsPointInFiducialVolume(TVector3 p)
{
  double Lx = 47.;
  double Ly = 40.;
  double Lz = 90.;
  if( (fabs(p.Y()       ) > Ly/2. - fFiducialMargin_Y) ||
      (fabs(p.X()-Lx/2. ) > Lx/2. - fFiducialMargin_X) ||
      (fabs(p.Z()-Lz/2. ) > Lz/2. - fFiducialMargin_Z) )
  {
    return false;
  } else {
    return true;
  }
}

DEFINE_ART_MODULE(MichelWfmReco)
