////////////////////////////////////////////////////////////////////////
// Class:       WireNoiseFilter
// Module Type: filter
// File:        WireNoiseFilter_module.cc
//
// Calculates the average wire RMS noise on the collection and induction
// planes and filters out events that fall outside of an allowable range.
//
// Sept 2018 -- W. Foreman
////////////////////////////////////////////////////////////////////////

// Framework includes
#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// LArSoft includes
#include "lardataobj/RawData/OpDetPulse.h"
#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RawData/raw.h"
#include "larevt/Filters/ChannelFilter.h"

// ROOT includes
#include <TH1F.h>
#include <TH2F.h>

// C++ includes
#include <iostream>
#include <memory>


class WireNoiseFilter;

class WireNoiseFilter : public art::EDFilter {
public:
  explicit WireNoiseFilter(fhicl::ParameterSet const & p);
  // The destructor generated by the compiler is fine for classes
  // without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  WireNoiseFilter(WireNoiseFilter const &) = delete;
  WireNoiseFilter(WireNoiseFilter &&) = delete;
  WireNoiseFilter & operator = (WireNoiseFilter const &) = delete;
  WireNoiseFilter & operator = (WireNoiseFilter &&) = delete;

  // Required functions.
  bool filter(art::Event & e) override;

  // Selected optional functions.
  void reconfigure(fhicl::ParameterSet const & p) override;

private:
  
  std::string fDAQModuleLabel;
  std::string fDAQModuleInstanceName;
  float       fMinWireRms[2];
  float       fMaxWireRms[2];
  float       fWireAdcThresh[2];
  int         fBaselineSamples;
  
  TH1F*       hMaxSignalPulse[2];
  TH1F*       hWireAdc[2];
  TH1F*       hWireRms[2];
  TH1F*       hAveWireRms[2];
  TH1F*       hAveWireRms_pass[2];
  TH1F*       hEvtPass;
};


WireNoiseFilter::WireNoiseFilter(fhicl::ParameterSet const & p)
// :
// Initialize member data here.
{
  this->reconfigure(p);
  art::ServiceHandle<art::TFileService> tfs;
  hMaxSignalPulse[0]= tfs->make<TH1F>("MaxSignalPulse_0","Induction plane;Wire pulse amplitude [ADC]",200,0,1000);
  hMaxSignalPulse[1]= tfs->make<TH1F>("MaxSignalPulse_1","Collection plane;Wire pulse amplitude [ADC]",200,0,1000);
  hWireAdc[0]    = tfs->make<TH1F>("WireAdc_0","Induction plane;Wire baseline [ADC]",200,-10,10);
  hWireAdc[1]    = tfs->make<TH1F>("WireAdc_1","Collection plane;Wire baseline [ADC]",200,-10,10);
  hWireRms[0]    = tfs->make<TH1F>("WireRms_0","Induction plane;Wire RMS [ADC]",200,0.,20);
  hWireRms[1]    = tfs->make<TH1F>("WireRms_1","Collection plane;Wire RMS [ADC]",200,0.,20);
  hAveWireRms[0]    = tfs->make<TH1F>("AveWireRms_0","Induction plane;Average Wire RMS [ADC]",200,0.,20);
  hAveWireRms[1]    = tfs->make<TH1F>("AveWireRms_1","Collection plane;Average Wire RMS [ADC]",200,0.,20);
  hAveWireRms_pass[0] = tfs->make<TH1F>("AveWireRms_pass_0","Induction plane;Wire RMS [ADC]",200,0.,20);
  hAveWireRms_pass[1] = tfs->make<TH1F>("AveWireRms_pass_1","Collection plane;Wire RMS [ADC]",200,0.,20);
  hEvtPass            = tfs->make<TH1F>("EvtPass","0 = fail, 1 = pass",2,0,2);
}

bool WireNoiseFilter::filter(art::Event & e)
{

  LOG_VERBATIM("WireNoiseFilter")
  <<"---- WireNoiseFilter -----\n";

  // Average wire RMS for each plane (initialize to dummy values)
  float aveRms[2]={-9.,-9.};
  
  // ----------------------------------------------------------------------
  // Check that raw digits exist
  art::Handle< std::vector<raw::RawDigit> > DigitHandle;;
  std::vector<art::Ptr<raw::RawDigit> > digit;
  if(e.getByLabel("daq",DigitHandle))
    {art::fill_ptr_vector(digit, DigitHandle);} 
  
  // ----------------------------------------------------------------------
  // Check wire RMS
  if( digit.size() > 0 ) {
    
    art::Handle< std::vector<raw::RawDigit> > digitVecHandle;
    e.getByLabel(fDAQModuleLabel, fDAQModuleInstanceName, digitVecHandle);
    
    raw::ChannelID_t channel = raw::InvalidChannelID; // channel number
    filter::ChannelFilter chanFilt;
    size_t dataSize = 0;
    
    float sum_rms[2] = {0,0};
    int   nchan[2]   = {0,0};

    for(size_t rdIter = 0; rdIter < digitVecHandle->size(); ++rdIter){
      
      // max pulse amplitude
      float maxpulse = -999.;

      // get the reference to the current raw::RawDigit
      art::Ptr<raw::RawDigit> digitVec(digitVecHandle, rdIter);
      channel = digitVec->Channel();
      dataSize = digitVec->Samples();
      std::vector<short> rawadc(dataSize);  // vector holding uncompressed adc values
      
      int plane = 1;
      if( channel < 240 ) plane = 0; 
      
      // skip bad channels
      if(!chanFilt.BadChannel(channel)) {
        // uncompress the data
        raw::Uncompress(digitVec->ADCs(), rawadc, digitVec->Compression());
        
        
        // calculate RMS assuming pedestal already subtracted off.
        // avoid outliers > [thresh] ADC
        int k = (int)rawadc.size();
        if( fBaselineSamples > 0 && k > fBaselineSamples ) k = fBaselineSamples;
        int nn = 0;
        float sum_sq = 0.;
        float sum_adc = 0.;
        for(int bin = 0; bin < k; ++bin) {
          if( fabs(rawadc[bin]) > maxpulse ) maxpulse = fabs(rawadc[bin]);
          if( fabs(rawadc[bin]) > fWireAdcThresh[plane] ) continue;
          sum_sq  += pow(rawadc[bin],2);
          sum_adc += rawadc[bin];
          nn++;  
        }
        if( nn ) {
          nchan[plane] ++;
          sum_rms[plane] += sqrt( sum_sq / nn );
          hWireAdc[plane] -> Fill( sum_adc / nn );
          hWireRms[plane] -> Fill( sum_sq / nn );
        }
        if( maxpulse > 0 ) hMaxSignalPulse[plane]->Fill(maxpulse);
      }
    } // Done looping over all wires
   
    for(int plane=0; plane<2; plane++){
      if( nchan[plane] > 0 ) aveRms[plane] = sum_rms[plane] / nchan[plane];
      LOG_VERBATIM("WireNoiseFilter")
      <<"   Plane "<<plane<<": "<<aveRms[plane]<<" ADC\n (acceptance range: "<<fMinWireRms[plane]<<"-"<<fMaxWireRms[plane]<<")";
    }// done loop over planes

  }// endif digits.size() > 0

  // -------------------------------------------------------------
  // Pass or fail the event. If any threshold is negative, ignore.
  bool pass = true;
  for(int plane=0; plane<2; plane++){
    if( aveRms[plane]>0 ) {
      hAveWireRms[plane]->Fill( aveRms[plane] );
      if( (fMinWireRms[plane]>0 && aveRms[plane] < fMinWireRms[plane] )
        ||(fMaxWireRms[plane]>0 && aveRms[plane] > fMaxWireRms[plane] ) ){
          pass = false;
      } 
    }
  }
  
  if( pass ) {
    for(int plane=0; plane<2; plane++){
      hAveWireRms_pass[plane]->Fill(aveRms[plane]);
    }
  }
  
  LOG_VERBATIM("WireNoiseFilter")
  <<"--------------------------\n";
  
  hEvtPass->Fill(int(pass));
  return pass;
  
}


void WireNoiseFilter::reconfigure(fhicl::ParameterSet const & p)
{
  fDAQModuleLabel       = p.get< std::string >  ("DAQModule","daq");
  fDAQModuleInstanceName= p.get< std::string >  ("DAQInstanceName","");
  fMinWireRms[0]        = p.get< float >        ("MinWireRmsInd",-9.);
  fMaxWireRms[0]        = p.get< float >        ("MaxWireRmsInd",-9.);
  fMinWireRms[1]        = p.get< float >        ("MinWireRmsCol",-9.);
  fMaxWireRms[1]        = p.get< float >        ("MaxWireRmsCol",-9.);
  fWireAdcThresh[0]     = p.get< float >        ("WireAdcThreshInd",30);
  fWireAdcThresh[1]     = p.get< float >        ("WireAdcThreshCol",30);
  fBaselineSamples      = p.get< int >          ("BaselineSamples",100);
}


DEFINE_ART_MODULE(WireNoiseFilter)
