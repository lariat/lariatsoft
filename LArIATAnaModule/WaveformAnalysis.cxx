#ifndef WAVEFORMANALYSIS_CXX
#define WAVEFORMANALYSIS_CXX

#include "WaveformAnalysis.h"
#include "TH1D.h"
#include "TSpectrum.h"


WaveformAnalysis::WaveformAnalysis(){}
WaveformAnalysis::~WaveformAnalysis(){}

bool WaveformAnalysis::FindPeaks(int npeaks, std::vector <short> fadc, double threshold, std::string polarity){

peaks.clear();
peaks_time.clear();

TSpectrum s(2*npeaks);

std::vector<double> spectrum;

if( polarity== "" || polarity == "negative")
{
  //cout<<"Such negativity"<<endl;
  spectrum.reserve(fadc.size());
  for(int i = 0; i < (int)fadc.size(); ++i) {spectrum.push_back(1./fadc[i]);}
}
else
{
  spectrum.assign(fadc.begin(),fadc.end());
}

wv = new TH1D("wv","",(int)spectrum.size(),0,(double)spectrum.size()-1);
for(int i = 0; i < (int)spectrum.size(); ++i){wv->SetBinContent(i,spectrum.at(i));}
TH1 *hb = s.Background(wv,20,"");

wv_clean = new TH1D("wv_clean","",(int)spectrum.size(),0,(double)spectrum.size()-1);

wv_clean->Add(hb,wv,-1);
/*for(int i = 0; i < hb->GetNbinsX(); i++){std::cout<<1./hb->GetBinContent(i)<<" "<<1./wv_clean->GetBinContent(i)<<std::endl;}*/
double sum_b;
double bsize;

bsize = (double)hb->GetNbinsX();
sum_b = 0;

for (int i = 1; i < bsize; i++) {sum_b += 1./hb->GetBinContent(i);}


baseline = sum_b/bsize;



sum_b = 0;
for (int i = 1; i < bsize; i++) sum_b += (1./hb->GetBinContent(i) - baseline)*(1./hb->GetBinContent(i) - baseline);
rms_baseline = sqrt(sum_b/bsize);
rms_baseline = rms_baseline;
nfound = s.Search(wv_clean,6,"",0.1);
xpeaks = s.GetPositionX();

//std::cout<<"Baseline: "<<baseline<<" "<<rms_baseline<<std::endl;

int xp;
int bin;
double peak;
peaks.reserve((size_t)nfound);
peaks_time.reserve((size_t)nfound);
bool found;

for(int i = 0; i < nfound; i++)
{
   
   xp = (int)xpeaks[i];
   bin = wv->GetXaxis()->FindBin(xp);
   peak =  wv->GetBinContent(bin);
   
   if( polarity== "" || polarity == "negative"){
  // std::cout<<abs(baseline - 1./peak)<<" "<<threshold<<std::endl; 
   if(abs(baseline - 1./peak) > threshold)
   {
     peak = 1./peak;  
     peaks.push_back(abs(peak - baseline));
     peaks_time.push_back(1e-9*xp);

   }
   }
   else{
   if(abs(peak-baseline) > threshold)
   {
     peaks.push_back(abs(peak-baseline));
     peaks_time.push_back(1e-9*xp);
   }
   }
   
}
//  std::cout<<"Waveform Analysis found "<<nfound<<" which passes the threshold!"<<std::endl;
  wv->Delete();
  wv_clean->Delete();
  if(nfound == 0 || peaks.size() == 0){found = 0;}
  else{found = 1;}
  return found;
}

int WaveformAnalysis::GetNPeaks(){return nfound;}

std::vector<short> WaveformAnalysis::GetPeaks(){return peaks;}

std::vector <double> WaveformAnalysis::GetPeaksTime(){return peaks_time;}

double WaveformAnalysis::GetPulseArea(std::vector<short> fadc, double peak_time)
{
  
  peak_position = int(peak_time/1e-9);
  int j;
  j = peak_position;
//  std::cout<<peak_position<<std::endl;
  short sum_rt;
  sum_rt = 0;
  short sum_dt;
  sum_dt = 0;
//  std::cout<<"Let's go back to baseline (RISE)"<<fadc[j+1]<<" "<<sum_rt<<" "<<baseline - 10*rms_baseline<<std::endl;
  while(fadc[j+1] < baseline - 10*rms_baseline)
  {
  //   std::cout<<"Let's go back to baseline (RISE)"<<fadc[j+1]<<" "<<sum_rt<<" "<<baseline - 10*rms_baseline<<std::endl;
     sum_rt = sum_rt+fadc[j];
     j++;
  }
  peak_rise_time = short(j - peak_position);
  j = peak_position;
  while(fadc[j-1] < baseline - 10*rms_baseline)
  {
  //   std::cout<<"Let's go back to baseline (DECAY)"<<fadc[j+1]<<" "<<sum_rt<<" "<< baseline - 10*rms_baseline <<std::endl; 
     sum_dt = sum_dt+fadc[j];
     j--;
  }
    peak_area = sum_dt + sum_rt;
    peak_decay_time = short(peak_position-j);
 return peak_area;

}

short WaveformAnalysis::GetRiseTime(){return peak_rise_time;}

short WaveformAnalysis::GetDecayTime(){return peak_decay_time;}

double WaveformAnalysis::GetBaseline()
{return baseline;}   

double WaveformAnalysis::GetBaselineRMS(){return rms_baseline;}
#endif
