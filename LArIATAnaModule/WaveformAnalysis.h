#ifndef LARIATANAMODULE_WAVEFORMANALYSIS_H
#define LARIATANAMODULE_WAVEFORMANALYSIS_H
#include <vector>
#include <algorithm>
#include <cmath>
#include <iterator>
#include <string>
//ROOT Libraries
#include "TF1.h"
#include "TGraph.h"
#include "TH1.h"
#include "TSpectrum.h"

class WaveformAnalysis{

public:
 WaveformAnalysis();
 ~WaveformAnalysis();
 bool FindPeaks(int, std::vector <short>, double, std::string);
 int GetNPeaks();
 std::vector <short> GetPeaks();
 std::vector <double> GetPeaksTime();
 double GetBaseline();
 double GetBaselineRMS();
 double GetPulseArea(std::vector<short>, double);
 short GetRiseTime();
 short GetDecayTime();
 

private:

 TH1 *wv;
 TH1 *wv_clean;
 Int_t npeaks;
 Int_t peak_position;
 Int_t nfound;
 std::string polarity;
 std::vector <short> peaks;
 std::vector <double> peaks_time;
 double peak_rise_time;
 double peak_decay_time;
 double baseline;
 double rms_baseline;
 Double_t* xpeaks;
 double peak_area;
 
};
 #endif
