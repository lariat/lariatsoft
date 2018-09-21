////////////////////////////////////////////////////////////////////////
// Class:       TimestampFilter
// Module Type: filter
// File:        TimestampFilter_module.cc
//
// This filter checks the event timestamp (within spill) and requires
// that it be within some user-defined limit.  Useful way to select beam
// events (1.2 - 5.5 sec).
//
// Option to require event contain non-empty container of TPC wire digits.
//
// Generated at Wed Oct 12 15:23:51 2016 by William Foreman using artmod
// from cetpkgsupport v1_10_02.
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

class TimestampFilter;

class TimestampFilter : public art::EDFilter {
public:
  explicit TimestampFilter(fhicl::ParameterSet const & p);
  // The destructor generated by the compiler is fine for classes
  // without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  TimestampFilter(TimestampFilter const &) = delete;
  TimestampFilter(TimestampFilter &&) = delete;
  TimestampFilter & operator = (TimestampFilter const &) = delete;
  TimestampFilter & operator = (TimestampFilter &&) = delete;

  // Required functions.
  bool filter(art::Event & e) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;
  void reconfigure(fhicl::ParameterSet const & p) override;

private:

  float fT1;
  float fT2;
  bool  fRequireRawDigits;
  std::string fDAQModuleLabel;
  std::string fDAQModuleInstanceName;

  TH1F* hTimestamps;
  TH1F* hTimestamps_pass;
  TH1F* hEvtSelection;

};


TimestampFilter::TimestampFilter(fhicl::ParameterSet const & p)
// :
// Initialize member data here.
{
  this->reconfigure(p);
  art::ServiceHandle<art::TFileService> tfs;
  hTimestamps       = tfs->make<TH1F>("Timestamps",";Time in spill [sec]",300,0,60.);
  hTimestamps_pass  = tfs->make<TH1F>("Timestamps_pass",";Time in spill [sec]",300,0.,60.);
  hEvtSelection     = tfs->make<TH1F>("EvtSelection","",3,0,3);
  hEvtSelection->SetOption("HIST TEXT");
  hEvtSelection->GetXaxis()->SetBinLabel(1,"Total evts");       
  hEvtSelection->GetXaxis()->SetBinLabel(2,"Raw dgts present"); 
  hEvtSelection->GetXaxis()->SetBinLabel(3,"Timestamp cut");  
}

bool TimestampFilter::filter(art::Event & e)
{
  // Set flags
  bool timestampFlag  = true;
  bool rawDigitFlag   = true;

  // ------------------------------------------------------------------------
  // First do timestamp filtering (this only applies to real data)
  if( e.isRealData() ) {
    
    //std::cout<<"TimestampFilter: run "<<e.run()<<", subrun "<<e.subRun()<<", event "<<e.id().event()<<"\n";
    // Get the timestamp (within the spill cycle) from the opdetpulse
    // objects because I don't know where else this info is saved!
    art::Handle< std::vector< raw::OpDetPulse >> opdetHandle;
    e.getByLabel(fDAQModuleLabel, fDAQModuleInstanceName, opdetHandle);
    
    float timeStamp = -1.;
    
    if( (size_t)opdetHandle->size() > 0 ){
      // All we want is the timestamp so just grab the first opdetpulse
      // we can find and call it a day
      art::Ptr< raw::OpDetPulse > ThePulsePtr(opdetHandle,0); 
      raw::OpDetPulse pulse = *ThePulsePtr;
      timeStamp = ((float)pulse.PMTFrame()*8.)/1.0e09;
      //std::cout<<"Timestamp = "<<timeStamp<<" sec\n";
      hTimestamps->Fill(timeStamp);
    }
  
    if( timeStamp < fT1 || timeStamp > fT2 ) {
      timestampFlag = false;
    } else {
      hTimestamps_pass->Fill(timeStamp);
    }
  
  }
    
  // ----------------------------------------------------------------------
  // Check that raw digits exist
  art::Handle< std::vector<raw::RawDigit> > DigitHandle;;
  std::vector<art::Ptr<raw::RawDigit> > digit;
  if(e.getByLabel("daq",DigitHandle))
    {art::fill_ptr_vector(digit, DigitHandle);} 
  //std::cout<<"Number of rawDigits: "<<digit.size()<<"\n";
  if( fRequireRawDigits && digit.size() == 0 ) rawDigitFlag = false;

  
  /*
  // ----------------------------------------------------------------------
  // Check wire RMS
  if( timestampFlag && rawDigitFlag ) {
    
    art::Handle< std::vector<raw::RawDigit> > digitVecHandle;
    e.getByLabel("daq", "", digitVecHandle);
    
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
      
      // only look at collection plane 
      //if( channel < 240 ) continue;
      
      int plane = 1;
      if( channel < 240 ) plane = 0; 
     
      
      // skip bad channels
      if(!chanFilt.BadChannel(channel)) {
        // uncompress the data
        raw::Uncompress(digitVec->ADCs(), rawadc, digitVec->Compression());
        
        // calculate RMS assuming pedestal already subtracted off.
        // avoid outliers > [thresh] ADC
        size_t k = rawadc.size();
        if( k > 300 ) k = 300;
        int nn = 0;
        float sum_sq = 0.;
        for(size_t bin = 0; bin < k; ++bin) {
          if( rawadc[bin] > maxpulse ) maxpulse = rawadc[bin];
          if( fabs(rawadc[bin])>fWireADCThresh ) continue;
          sum_sq += pow(rawadc[bin],2);
          nn++;  
        }
        if( nn ) {
          sum_rms[plane] += sqrt(sum_sq / nn );
          nchan[plane] ++;
        }
        if( maxpulse > 0 ) hMaxSignalPulse[plane]->Fill(maxpulse);
      }
    } // Done looping over all wires
   
    for(int plane=0; plane<2; plane++){
      float aveRms = -9.;
      if( nchan[plane] > 0 ) aveRms = sum_rms[plane] / nchan[plane];
      hAveWireRms[plane]->Fill( aveRms );
      if( plane == 1 ) {
        if( fMaxWireRms > 0. && ( aveRms < 0. || aveRms > fMaxWireRms || aveRms < fMinWireRms ) ) {
          wireRmsFlag = false;
        } else {
          hAveWireRms_pass->Fill( aveRms );
        }
      }
    }
    

  }
  */


  bool passFlag = false;
  hEvtSelection->Fill(0);
  if( rawDigitFlag ) {
    hEvtSelection->Fill(1);
    if( timestampFlag ) {
      hEvtSelection->Fill(2);
      passFlag = true;
      /*
      if( wireRmsFlag ) {
        hEvtSelection->Fill(3);
        passFlag = true;
      }
      */
    }
  }
  
  return passFlag;

}

void TimestampFilter::beginJob()
{
}


void TimestampFilter::endJob()
{
}


void TimestampFilter::reconfigure(fhicl::ParameterSet const & p)
{
  fT1				= p.get< float  >("T1",0);
  fT2				= p.get< float  >("T2",60);
  fDAQModuleLabel               = p.get< std::string >  ("DAQModule","daq");
  fDAQModuleInstanceName        = p.get< std::string >  ("DAQInstanceName","");
  fRequireRawDigits             = p.get< bool >         ("RequireRawDigits",true);
}


DEFINE_ART_MODULE(TimestampFilter)
