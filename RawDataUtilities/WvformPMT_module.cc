////////////////////////////////////////////////////////////////////////
// Class:       WvformPMT
// Module Type: analyzer
// File:        WvformPMT_module.cc
//
// written by Pawel Kryczynski based on the DigitReader_module Generated at Tue Feb 10 15:40:09 2015 by Will Flanagan using artmod
// from cetpkgsupport v1_08_02 and OpDigiAna_module (by Christie Chiu and Ben Jones, MIT, 2012) and CalibrationTPC_Algs by D. Caratelli and WArP analysis code
////////////////////////////////////////////////////////////////////////
#define _USE_MATH_DEFINES
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/PtrVector.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "TVirtualFFT.h"

#include "TF1.h"

#include "lardataobj/RawData/OpDetPulse.h"
#include "lardataobj/RecoBase/OpHit.h"

// ROOT includes
#include "TH1.h"
#include "THStack.h"
#include "TF1.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TArrayD.h"
// C++ Includes
#include <map>
#include <vector>
#include <iostream>
#include <cstring>
#include <sstream>
#include "math.h"
#include <climits>
#include <memory>
class WvformPMT;

class WvformPMT : public art::EDAnalyzer {
public:
  explicit WvformPMT(fhicl::ParameterSet const & p);
  // The destructor generated by the compiler is fine for classes
  // without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  WvformPMT(WvformPMT const &) = delete;
  WvformPMT(WvformPMT &&) = delete;
  WvformPMT & operator = (WvformPMT const &) = delete;
  WvformPMT & operator = (WvformPMT &&) = delete;

  // Required functions.
  void analyze(art::Event const & e) override;

  // Selected optional functions.
  void beginJob() override;
 // void reconfigure(fhicl::ParameterSet const & p) override;
  void endJob() override;
static Double_t twoexp(Double_t *x,Double_t *par);

private:


  // Declare member data here.
    std::string fInputModule;              // Input tag for OpDet collection
    std::string fInstanceName;             // Input tag for OpDet collection
    float fSampleFreq;                     // in MHz
    float fTimeBegin;                      // in us
    float fTimeEnd;                        // in us

   double fFitConst;
   double fExpM1;
   double fExpM2;
   double fExpM3;
   double fExpS1;
   double fExpS1L;
   double fExpS1H;
   double fExpS2;
   double fExpS2L;
   double fExpS2H;
   double fExpS3;
	double fP1;
	double fSigma;
	double fLowLimit;
	double fHighLimit;

   double fFitStart;
   double fFitEnd;
   int fSamples;
    TF1 *fitexp2;
    TF1 *fitconv;

    int  fNumberOfPMTs;
    int  fBaselineCounts;

    char fHistName[50];
  TH1D * avhist [5];
  TH1D * avhist_raw [5];
  TH1D * avhist_fft [5];
  TH1D * fft_scaled [5];
int PMTWvformCount[5];
std::vector<std::vector <double>> bins;
std::vector <double> baseline;
std::vector <double> res;
std::vector <double> ims;
int fEvents;
int pmt1;
std::vector<std::vector<unsigned int>> fOpDetChID;
std::vector<unsigned int> PMTsInRun;
art::ServiceHandle<art::TFileService> tfs;

};

//------------------------------------------------------------------------------

Double_t WvformPMT::twoexp(Double_t *x,Double_t *par)
{
Double_t arg = 0;
arg = x[0];
Double_t twoe =par[0]+ par[1]*TMath::Exp(-arg/par[2]) +par[3]*TMath::Exp(-arg/par[4])+par[5]*TMath::Exp(-arg/par[6])+par[7]*TMath::Exp(-arg/par[8]);
return twoe;
}


//-----------------------------------------------------------------------------
WvformPMT::WvformPMT(fhicl::ParameterSet const & p)
  : EDAnalyzer(p)
{


    fTimeBegin  =  p.get<double>("time_begin",0.0);
    fTimeEnd    =  p.get<double>("time_end",7168.0);
    fSampleFreq =  p.get<double>("sample_freq",1.0);//in samples per ns=GHz


    // Indicate that the Input Module comes from .fcl
    fInputModule = p.get<std::string>("InputModule");
    fInstanceName = p.get<std::string>("InstanceName");
 

   fOpDetChID = p.get< std::vector<std::vector<unsigned int>> >("pmt_channel_ids");
   fNumberOfPMTs = p.get<int>("NumberOfPMTs",5);
   fBaselineCounts = p.get<int>("baseline_counts",90);
   fFitConst = p.get<double>("fit_const",0.0);
   fExpM1 = p.get<double>("fit_expmult1",1.0);
   fExpM2 = p.get<double>("fit_expmult2",1.0);
   fExpM3 = p.get<double>("fit_expmult3",1.0);
   fExpM3 = p.get<double>("fit_expmult4",0.02);
   fExpS1 = p.get<double>("fit_exp_slope1",15.0);
   fExpS1L = p.get<double>("fit_exp_slope1low",1.0);
   fExpS1H = p.get<double>("fit_exp_slope1high",20.0);
   fExpS2 = p.get<double>("fit_exp_slope2",50.0);
   fExpS2L = p.get<double>("fit_exp_slope2low",40.0);
   fExpS2H = p.get<double>("fit_exp_slope2high",60.0);
   fExpS3 = p.get<double>("fit_exp_slope3",3000.0);
   fFitStart = p.get<double>("fit_start",1250.0);
   fFitEnd = p.get<double>("fit_end",5000.0);
   fSamples = p.get<int>("SamplesNumber",7168);

}



//------------------------------------------------------------------------------
void WvformPMT::beginJob()
{

	fEvents=0;
	pmt1=0;
	PMTsInRun.clear();
	for(auto chid: fOpDetChID){
		if(int(chid.size())>0){
			pmt1+=chid.size();
			for(auto ch: chid){
				PMTsInRun.push_back(ch);
			}
		}
	}
	for(int i=0;i <fNumberOfPMTs;++i) PMTWvformCount[i]=0;//zeroing pmt wvform count	
	bins.resize(fNumberOfPMTs);
	baseline.resize(fNumberOfPMTs);
        for(int c=0;c<fNumberOfPMTs;++c) {
            sprintf(fHistName, "OpDet_%i_AvPulse", c);

	    avhist[c]= tfs->make<TH1D>(fHistName, ";t (ns);", 
					     int((fTimeEnd - fTimeBegin) * fSampleFreq), 
					     fTimeBegin, 
					     fTimeEnd);

            sprintf(fHistName, "OpDet_%i_AvPulse_raw", c);

	    avhist_raw[c]= tfs->make<TH1D>(fHistName, ";t (ns);", 
					     int((fTimeEnd - fTimeBegin) * fSampleFreq), 
					     fTimeBegin, 
					     fTimeEnd);
            sprintf(fHistName, "OpDet_%i_AvPulse_fft", c);

	    avhist_fft[c]= tfs->make<TH1D>(fHistName, ";freq;", 
					     (fSamples/2)+1, 
					     0, 
					     (fSamples/2)+1);

            sprintf(fHistName, "OpDet_%i_AvPulse_fft_scaled", c);

	    fft_scaled[c]= tfs->make<TH1D>(fHistName, ";freq,Hz;", 
					     (fSampleFreq*5000/2.), 0, 
					      (fSampleFreq*1000000000/2.));



        bins[c].clear();
        bins[c].resize(int((fTimeEnd - fTimeBegin) * fSampleFreq),0.0);
	  
        }

        

}

//------------------------------------------------------------------------------
void WvformPMT::analyze(art::Event const & evt)
{
	++fEvents;

    art::Handle< std::vector< raw::OpDetPulse > > WaveformHandle;

	evt.getByLabel(fInputModule, fInstanceName, WaveformHandle);
    
	std::cout<<" wvform size "<<WaveformHandle->size()<<std::endl;
    std::vector<std::string> histnames;

    for(unsigned int i = 0; i < WaveformHandle->size(); ++i){ 

			art::Ptr< raw::OpDetPulse > ThePulsePtr(WaveformHandle, i);
			raw::OpDetPulse ThePulse = *ThePulsePtr;
			for(int j =0;j<int(ThePulse.Waveform().size());++j){

				if(int(j)<int(bins[int(ThePulse.OpChannel())].size())) bins[int(ThePulse.OpChannel())].at(j)=bins[int(ThePulse.OpChannel())].at(j)+(double) ThePulse.Waveform()[j];
  		 }
	PMTWvformCount[int(ThePulse.OpChannel())]=PMTWvformCount[int(ThePulse.OpChannel())]+1;
	std::cout<<"pmt "<<int(ThePulse.OpChannel())<<" "<<PMTWvformCount[int(ThePulse.OpChannel())]<<std::endl;


    }

}

//------------------------------------------------------------------------------
void WvformPMT::endJob()
{
  
  for(int c=0;c<fNumberOfPMTs;++c) {
    fitexp2 = new TF1("exp2",twoexp,0,15000,7);
    for(int j =0;j<fBaselineCounts;j++) baseline.at(c)=baseline.at(c)+double(bins[c].at(j));
    baseline.at(c)=baseline.at(c)/fBaselineCounts;
    for(int j =0;j<int(avhist[c]->GetNbinsX());++j) {
      if (PMTWvformCount[c]!=0 && int(j)<int(bins[c].size())) avhist[c]->SetBinContent(j,(-double(bins[c].at(j))+baseline.at(c))/PMTWvformCount[c]);
      if (PMTWvformCount[c]!=0 && int(j)<int(bins[c].size())) avhist_raw[c]->SetBinContent(j,(double(bins[c].at(j))/PMTWvformCount[c]));
    }
      //scaling e.g http://en.wikipedia.org/wiki/Short-time_Fourier_transform, fftÂ from root example
  		double scaling_fft=double(fSampleFreq*1000000000)/double(fSamples);
    TH1* temp_hist=NULL;
    temp_hist=avhist[c]->FFT(temp_hist,"MAG");
    TVirtualFFT *fft = TVirtualFFT::GetCurrentTransform();
    double re1=0.;
    double im1=0.;
    for (int it=0; it<(fSamples/2+1); ++it){
      fft->GetPointComplex(it, re1, im1);
      avhist_fft[c]->SetBinContent(it,temp_hist->GetBinContent(it));
      fft_scaled[c]->SetBinContent(it,scaling_fft*temp_hist->GetBinContent(it));
    }
    
    fitexp2->SetRange(Double_t(fFitStart),Double_t(fFitEnd));
    fitexp2->SetParameter(0,Double_t(fFitConst));
    fitexp2->SetParameter(1,Double_t(fExpM1));
    fitexp2->SetParameter(2,Double_t(fExpS1));
    fitexp2->SetParLimits(2,Double_t(fExpS1L), Double_t(fExpS1H));
    fitexp2->SetParameter(3,Double_t(fExpM2));
    fitexp2->SetParameter(4,Double_t(fExpS2));
    fitexp2->SetParameter(5,Double_t(fExpM3));
    fitexp2->SetParameter(6,Double_t(fExpS3));
    
    
    avhist[c]->Fit("exp2", "REM+");
    delete temp_hist;
    delete fft;
    delete fitexp2;
  }

   
}

DEFINE_ART_MODULE(WvformPMT)

