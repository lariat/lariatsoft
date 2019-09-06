////////////////////////////////////////////////////////////////////////
// Class:       SignalToNoise
// Module Type: filter
// File:        SignalToNoise_module.cc
//
// Simplified module to measure the raw wire signal-to-noise, based on
// analysis work from Silvia Zhang.
//
// This module must be run over reconstructed data.
// The following data products are expected:
//  * raw::Wires
//  * recob::Tracks
//
// W. Foreman
// Aug 2019
////////////////////////////////////////////////////////////////////////

// Framework includes
#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "canvas/Utilities/InputTag.h"
#include "canvas/Persistency/Common/FindOneP.h" 
#include "canvas/Persistency/Common/Ptr.h" 
#include "canvas/Persistency/Common/PtrVector.h" 
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib/maybe_ref.h"
#include "cetlib/exception.h"
#include "cetlib/search_path.h"

// LArSoft includes
#include "lardata/Utilities/AssociationUtil.h"
#include "lardata/DetectorInfo/DetectorPropertiesStandard.h"
#include "lardata/DetectorInfoServices/LArPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardataobj/RawData/OpDetPulse.h"
#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RawData/raw.h"
#include "lardataobj/RecoBase/Track.h"
#include "larevt/Filters/ChannelFilter.h"
#include "LArIATDataProducts/WCTrack.h"

// ROOT includes
#include <TH1.h>
#include <TH2D.h>
#include <TTree.h>
#include <TProfile.h>

// C++ includes
#include <iostream>
#include <memory>
#include <algorithm>

const int NWires         = 480;   // Total number of Wires in TPC, need to mod for 5mm (384 wires)
const int NTTicks        = 3072;  // Number of time ticks in an event, 1 ttick = 128 ns

class SignalToNoise;

class SignalToNoise : public art::EDFilter {
public:
  explicit SignalToNoise(fhicl::ParameterSet const & p);
  
  SignalToNoise(SignalToNoise const &) = delete;
  SignalToNoise(SignalToNoise &&) = delete;
  SignalToNoise & operator = (SignalToNoise const &) = delete;
  SignalToNoise & operator = (SignalToNoise &&) = delete;
  
  bool filter(art::Event & e) override;
  void beginJob() override;
  void beginSubRun(art::SubRun const & sr);
  void reconfigure(fhicl::ParameterSet const & p) override;

private:
  
  TTree*            fTree;
  int               run;                        //run number
  int               subrun;                     //subrun number
  int               event;                      //event
  float             timestamp;                  //timestamp of subrun in unix time (supposedly...)
  raw::ChannelID_t  channel;                    //wire channel, 0-239 induction, 240-480 collection
  int               vPeakval[NWires];           //Vector of the biggest peak on each wire          
  int               wireraw[NWires][NTTicks];   //How we readout raw data                  
  double            pedmean[NWires];            //we calculate the pedestal mean...         
  double            pedrms[NWires];             //and rms        
  
  // Time range in which to look for beam tracks
  // TODO: replace this with WC track matching 
  int               timeTickMin     = 300;      //ignore all waveform info prior to this (PMT-induced noise)
  int               signalRegionT1  = 500;     
  int               signalRegionT2  = 1800;

  // Rudimentary hit-finding on waveforms
  int		    signalheight = 18;          //we define a hit as above pedrms[i]*signalheight for currentsignalwidth time ticks.
  int		    Csignalwidth = 20;          //min signal width on collection
  int		    Isignalwidth = 10;          //min signal width on induction
  int               minWidth[2];

  // module names
  std::string	    fDAQModuleLabel;
  std::string	    fDAQModuleInstanceName;
  std::string       fWCTrackLabel;
  std::string	    fTrackModuleLabel;

  // histograms
  TProfile*	    TProfPed;
  TH2D*             ChitSignals;
  TH2D*             IhitSignals;
  TH1D*             hElectronLifetime;
  TH1D*             hRawPed[2];
  TH1D*             hADC[2];
  TH2D*             hTimeTick[2];
  TH2D*             hTimeTickSignal[2];
  TH1D*             hSignal[2];
  TH1D*             hRMS[2];
  TH1D*             hSNR[2];
  //TH1F*           hWireAdc[2];
  //TH1F*           hWireRms[2];
  
  detinfo::DetectorProperties const* fDetProp;
};

void SignalToNoise::beginJob()
{
  // Implementation of optional member function here.
  art::ServiceHandle<art::TFileService> tfs;

  fTree = tfs->make<TTree>("fTree","analysis tree");
  fTree->Branch("run"      ,&run      ,"run/I");
  fTree->Branch("subrun"   ,&subrun   ,"subrun/I");
  fTree->Branch("event"   ,&event   ,"event/I");
  fTree->Branch("timestamp", &timestamp,"timestamp/I");
  fTree->Branch("pedmean",pedmean,"pedmean[480]/D");
  fTree->Branch("pedrms",pedrms,"pedrms[480]/D");
  fTree->Branch("vPeakval",vPeakval,"vPeakval[480]/I");

  TProfPed = tfs->make<TProfile>("TProfilePed","Pedestal TProfile Mean/RMS", 480, 0,480, -1000, 1000, "s");
  TProfPed->SetAxisRange(-5.,5.,"Y");
  TProfPed->SetStats(0);
  TProfPed->GetXaxis()->SetTitle("Channel Number");
  TProfPed->GetYaxis()->SetTitle("ADC Counts");

  ChitSignals       = tfs->make<TH2D>("CwireSignals","Collection Wire Signals", 200, -100, 100, 1000, -300, 300);  
  ChitSignals->GetXaxis()->SetTitle("Time [TTicks]");
  ChitSignals->GetXaxis()->SetTitle("ADC Counts");

  IhitSignals       = tfs->make<TH2D>("IwireSignals","Induction Wire Signals", 200, -100, 100, 1000, -300, 300);
  IhitSignals->GetXaxis()->SetTitle("Time [TTicks]");
  IhitSignals->GetXaxis()->SetTitle("ADC Counts");
  
  hElectronLifetime       = tfs->make<TH1D>("ElectronLifetime","Electron lifetime from database;Electron lifetime [#mus]",2500,0.,2500.);
 
  hRawPed[0]    = tfs->make<TH1D>("RawPed_0","Induction plane;Raw wire pedestal [ADC]",80,-20,20);
  hRawPed[1]    = tfs->make<TH1D>("RawPed_1","Collection plane;Raw wire pedestal [ADC]",80,-20,20);
  hADC[0]       = tfs->make<TH1D>("ADC_0","Induction plane;Raw wire ADC amplitude [ADC]",200,-100,100);
  hADC[1]       = tfs->make<TH1D>("ADC_1","Collection plane;Raw wire ADC amplitude [ADC]",200,-100,100);
  hTimeTick[0]  = tfs->make<TH2D>("TimeTick_0","Induction plane;Wire number;Time-tick of pulse", 240,0,240,  3072,0,3072);
  hTimeTick[1]  = tfs->make<TH2D>("TimeTick_1","Collection plane;Wire number;Time-tick of pulse",240,240,480,3072,0,3072);
  hTimeTickSignal[0]  = tfs->make<TH2D>("TimeTickSignal_0","Induction plane;Wire number;Time-tick of signal pulse", 240,0,240,  3072,0,3072);
  hTimeTickSignal[1]  = tfs->make<TH2D>("TimeTickSignal_1","Collection plane;Wire number;Time-tick of signal pulse",240,240,480,3072,0,3072);
  hSignal[0]    = tfs->make<TH1D>("Signal_0","Induction plane;Raw wire signal amplitude [ADC]",100,0,500);
  hSignal[1]    = tfs->make<TH1D>("Signal_1","Collection plane;Raw wire signal amplitude [ADC]",100,0,500);
  hRMS[0]       = tfs->make<TH1D>("RMS_0","Induction plane;Raw wire RMS [ADC]",100,0,25);
  hRMS[1]       = tfs->make<TH1D>("RMS_1","Collection plane;Raw wire RMS [ADC]",100,0,25);
  //hWireAdc[1]    = tfs->make<TH1F>("WireAdc_1","Collection plane;Wire baseline [ADC]",200,-10,10);
  //hWireRms[0]    = tfs->make<TH1F>("WireRms_0","Induction plane;Wire RMS [ADC]",200,0.,20);
  //hWireRms[1]    = tfs->make<TH1F>("WireRms_1","Collection plane;Wire RMS [ADC]",200,0.,20);
  hSNR[0]       = tfs->make<TH1D>("SNR_0","Induction plane;Signal-to-noise ratio",100,0,100);
  hSNR[1]       = tfs->make<TH1D>("SNR_1","Collection plane;Signal-to-noise ratio",100,0,100);
 
  hTimeTick[0]  ->SetOption("colz");
  hTimeTick[1]  ->SetOption("colz");
  hTimeTickSignal[0]  ->SetOption("colz");
  hTimeTickSignal[1]  ->SetOption("colz");
  
}


SignalToNoise::SignalToNoise(fhicl::ParameterSet const & p){
  this->reconfigure(p);
  
  // Initialize detprop pointer
  fDetProp = lar::providerFrom<detinfo::DetectorPropertiesService>();
  
  minWidth[0]=Isignalwidth;
  minWidth[1]=Csignalwidth;
}

void SignalToNoise::beginSubRun(art::SubRun const & sr)
{
  // Implementation of optional member function here.
  memset(pedmean, -999, sizeof(pedmean));
  memset(pedrms, -999, sizeof(pedrms));
}

bool SignalToNoise::filter(art::Event & e)
{

  //LOG_VERBATIM("SignalToNoise")
  //<<"---- SignalToNoise -----\n";
  
  // ------------------------------------------------------------------------
  // reset tree variables
  for(int ichannel = 0; ichannel < NWires; ichannel++){ //empty rawadc storage
    for(int ibin = 0; ibin <NTTicks; ibin++) wireraw[ichannel][ibin] = -999;
  }
  memset(vPeakval, -543, sizeof(vPeakval));
  
  run       = e.run();
  subrun    = e.subRun();
  event     = e.id().event();
  timestamp = (float)e.getSubRun().beginTime().value();
 
  // ----------------------------------------------------------------------
  // Get raw digits
  art::Handle< std::vector<raw::RawDigit> > DigitHandle;;
  std::vector<art::Ptr<raw::RawDigit> > digit;
  if(e.getByLabel(fDAQModuleLabel,DigitHandle)) 
    art::fill_ptr_vector(digit, DigitHandle); 
  if( digit.size() == 0 ) return false;

  // ---------------------------------------------------------------------
  // Get TPC tracks
  art::Handle< std::vector<recob::Track> > trackListHandle;
  std::vector<art::Ptr<recob::Track> > tpctrack; 
  if (e.getByLabel(fTrackModuleLabel,trackListHandle)) 
    art::fill_ptr_vector(tpctrack, trackListHandle);
 
  // ---------------------------------------------------------------------
  // Get WC tracks
  art::Handle< std::vector<ldp::WCTrack> > wctrackHandle;
  std::vector<art::Ptr<ldp::WCTrack> > wctrack;
  if(e.getByLabel(fWCTrackLabel, wctrackHandle))
    art::fill_ptr_vector(wctrack, wctrackHandle);
  
  // 1) Filter out events containing exactly 1 TPC track and 1 WC track
  // 2) Loop over all wires and for each one,
  //       i)   Get RMS by masking the raw digit surrounding any "hits"
  //       ii)  Look for max ADC and find the hit "width" to check that it's not noise

  // First, throw out busy events
  if( tpctrack.size() > 3 ) return false;

  // Define signal event as containing 1 WC track and 1 through-going TPC track
  bool isSignalEvent = false;
  if( wctrack.size() == 1 && tpctrack.size() == 1 ) {
    float minLength = 80.;
    if( run > 13421 ) minLength = 50.; // Run-IIIB had smaller active volume
    if( tpctrack.at(0)->Length() >= minLength ) isSignalEvent = true;
  }
  
  art::Handle< std::vector<raw::RawDigit> > digitVecHandle;
  e.getByLabel(fDAQModuleLabel, fDAQModuleInstanceName, digitVecHandle);
  channel = raw::InvalidChannelID; // channel number
  filter::ChannelFilter chanFilt;
  int dataSize = 0;
    
  for(size_t rdIter = 0; rdIter < digitVecHandle->size(); ++rdIter){
      
    // get the reference to the current raw::RawDigit
    art::Ptr<raw::RawDigit> digitVec(digitVecHandle, rdIter);
    channel = digitVec->Channel();
    dataSize = (int)digitVec->Samples();
    std::vector<short> rawadc(dataSize);  // vector holding uncompressed adc values
     
    int plane = 1;
    if( channel < 240 ) plane = 0; 
    // TODO: what is the numbering scheme for Run 3A where there are only 192 wires on each plane?

    // skip bad channels
    if(!chanFilt.BadChannel(channel)) {
      raw::Uncompress(digitVec->ADCs(), rawadc, digitVec->Compression());
 
      // when calculating pedestal and RMS, only use part of the
      // waveform that precedes any large signals to avoid biases
      // due to undershoot and irrelevant post-pulse fluctuations
      //std::vector<bool> mask(dataSize,true);
      int ped_max_i = -9;
      for(int i=0; i<dataSize; i++){
        if( abs(rawadc.at(i)) > 20 ) {
          hTimeTick[plane]->Fill( int(channel), i );
          if( ped_max_i < 0 ) ped_max_i = std::max(0,i-30);
        }
      }

      // require at least 80 samples in the baseline
      if( ped_max_i >= 80 ) {
      
        // calculate pedestal mean
        float sum = 0;
        int n = 0;
        for(int i=0; i<ped_max_i; i++){
          sum += rawadc.at(i);
          n++;
        }
        float ped = -999;
        if(n) ped = sum/float(n);
        hRawPed[plane]->Fill(ped);
   
        // subtract off pedestal and recast waveform in floats, and
        // find the max pulse. also calculate RMS while we're at it.
        float maxAdc    = -999;
        int   maxAdc_i  = -9;
        float sumsq     = 0;
        float rms       = -9;
        std::vector<float> adc(dataSize);
        for(int i=0; i<dataSize; i++) {
          adc.at(i) = float(rawadc.at(i)) - ped;
          hADC[plane]->Fill( adc.at(i) );
          if( i >= timeTickMin && adc.at(i) > maxAdc 
            && i >= signalRegionT1 && i <= signalRegionT2 ) {
            maxAdc = adc.at(i);
            maxAdc_i = i;
          }
          if(i < ped_max_i ) sumsq += adc.at(i)*adc.at(i);
        }
        rms=sqrt(sumsq/(n-1));
       
        // verify that the peak has a realistic width
        bool goodHit = false;
        int pulseWidth = -9;
        int start_i = -9;
        int end_i = -9;
        float thresh = 9999;
        if( rms > 0 ) thresh = std::max(15.,10.*rms);
        if(maxAdc_i > 0 ) {
          for(short i=std::max(maxAdc_i-150,0); i<std::min(maxAdc_i+150,dataSize); i++){
            if( start_i < 0   && adc.at(i) >= thresh ) start_i = i;
            if( start_i >= 0  && adc.at(i) < thresh ) { 
              end_i = i; 
              pulseWidth = end_i - start_i;
              break; 
            }
          }
        }
        if( pulseWidth >= minWidth[plane] ) goodHit = true; 

        // ********************
        pedmean[channel]  = ped;
        pedrms[channel]   = rms;
        hRMS[plane]       ->Fill(rms);
        if( goodHit && isSignalEvent ) {
          hTimeTickSignal[plane]->Fill(int(channel), maxAdc_i );
          hSignal[plane]->Fill(maxAdc);
          hSNR[plane]->Fill( maxAdc / rms );
          hElectronLifetime->Fill(fDetProp->ElectronLifetime());
        }
        // ********************
   
        // find max ADC and check that the width looks OK
      
      }//endif good wire
    
    }//chanfilt

  }//loop over wires

  fTree->Fill();
 
  return true;
  
}


void SignalToNoise::reconfigure(fhicl::ParameterSet const & p)
{
  fDAQModuleLabel       = p.get< std::string >  ("DAQModule","daq");
  fDAQModuleInstanceName= p.get< std::string >  ("DAQInstanceName","");
  fTrackModuleLabel     = p.get< std::string >  ("TrackModuleLabel","pmtracktc");
  fWCTrackLabel         = p.get< std::string >  ("WCTrackLabel","wctrack");
}


DEFINE_ART_MODULE(SignalToNoise)
