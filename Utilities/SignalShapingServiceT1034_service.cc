
////////////////////////////////////////////////////////////////////////
/// \file   SignalShapingServiceT1034_service.cc
/// \author H. Greenlee   (adapted to LArIAT by A. Szelc)
////////////////////////////////////////////////////////////////////////

#include "Utilities/SignalShapingServiceT1034.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib/exception.h"
#include "Geometry/Geometry.h"
#include "Geometry/TPCGeo.h"
#include "Geometry/PlaneGeo.h"
#include "Utilities/DetectorProperties.h"
#include "Utilities/LArProperties.h"
#include "Utilities/LArFFT.h"
#include "TFile.h"

//----------------------------------------------------------------------
// Constructor.
util::SignalShapingServiceT1034::SignalShapingServiceT1034(const fhicl::ParameterSet& pset,
								    art::ActivityRegistry& /* reg */) 
  : fInit(false)
{
  
  //std::cout<<std::endl;
  //std::cout<<"Constructor"<<std::endl;
  //std::cout<<std::endl;
  reconfigure(pset);
}


//----------------------------------------------------------------------
// Destructor.
util::SignalShapingServiceT1034::~SignalShapingServiceT1034()
{
//std::cout<<std::endl;
//std::cout<<"Destructor"<<std::endl;
//std::cout<<std::endl;
}


//----------------------------------------------------------------------
// Reconfigure method.
void util::SignalShapingServiceT1034::reconfigure(const fhicl::ParameterSet& pset)
{
  // Reset initialization flag.
  
  //std::cout<<std::endl;
  //std::cout<<"Reconfigure"<<std::endl;
  //std::cout<<std::endl;
  
  fInit = false;

  // Reset kernels.

  fColSignalShaping.Reset();
  fIndVSignalShaping.Reset();

  // Fetch fcl parameters.

  fADCTicksPerPCAtLowestASICGainSetting = pset.get<double>("ADCTicksPerPCAtLowestASICGainSetting");
  fASICGainInMVPerFC = pset.get<double>("ASICGainInMVPerFC");
  
  // #############################################
  // ### Number of bins for the field response ###
  // #############################################
  fNFieldBins = pset.get<int>("FieldBins");
  fCol3DCorrection = pset.get<double>("Col3DCorrection");
  fInd3DCorrection = pset.get<double>("Ind3DCorrection");
  fColFieldRespAmp = pset.get<double>("ColFieldRespAmp");
  fIndVFieldRespAmp = pset.get<double>("IndVFieldRespAmp");
  fShapeTimeConst = pset.get<std::vector<double> >("ShapeTimeConst");

  fUseFunctionFieldShape= pset.get<bool>("UseFunctionFieldShape");
  fUseSimpleFieldShape = pset.get<bool>("UseSimpleFieldShape");
  if(fUseSimpleFieldShape) {
    fNFieldBins = 300;
  }
  fGetFilterFromHisto= pset.get<bool>("GetFilterFromHisto");
  
  // Construct parameterized collection filter function.
  if(!fGetFilterFromHisto) {

    mf::LogInfo("SignalShapingServiceT1034") << "Getting Filter from .fcl file" ;
    std::string colFilt = pset.get<std::string>("ColFilter");
    std::vector<double> colFiltParams = pset.get<std::vector<double> >("ColFilterParams");
    fColFilterFunc = new TF1("colFilter", colFilt.c_str());
    for(unsigned int i=0; i<colFiltParams.size(); ++i)
      fColFilterFunc->SetParameter(i, colFiltParams[i]);
    
    // Construct parameterized induction filter function.

    std::string indVFilt = pset.get<std::string>("IndVFilter");
    std::vector<double> indVFiltParams = pset.get<std::vector<double> >("IndVFilterParams");
    fIndVFilterFunc = new TF1("indVFilter", indVFilt.c_str());
    for(unsigned int i=0; i<indVFiltParams.size(); ++i)
      fIndVFilterFunc->SetParameter(i, indVFiltParams[i]);
  } else {
  
    std::string histoname = pset.get<std::string>("FilterHistoName");
    mf::LogInfo("SignalShapingServiceT1034") << " using filter from .root file " ;
    //art::ServiceHandle<geo::Geometry> geom;
    int fNPlanes=2;
   
    // constructor decides if initialized value is a path or an environment variable
    std::string fname;
    cet::search_path sp("FW_SEARCH_PATH");
    sp.find_file(pset.get<std::string>("FilterFunctionFname"), fname);
    
    TFile * in=new TFile(fname.c_str(),"READ");
    for(int i=0;i<fNPlanes;i++){
      TH1D * temp=(TH1D *)in->Get(Form(histoname.c_str(),i));
      fFilterHist[i]=new TH1D(Form(histoname.c_str(),i),Form(histoname.c_str(),i),temp->GetNbinsX(),0,temp->GetNbinsX());
      temp->Copy(*fFilterHist[i]); 
    }
   
    in->Close();
   
  }
 
  /////////////////////////////////////
  if(fUseFunctionFieldShape) {

    std::string colField = pset.get<std::string>("ColFieldShape");
    std::vector<double> colFieldParams = pset.get<std::vector<double> >("ColFieldParams");
    fColFieldFunc = new TF1("colField", colField.c_str());
    for(unsigned int i=0; i<colFieldParams.size(); ++i)
      fColFieldFunc->SetParameter(i, colFieldParams[i]);

    // Construct parameterized induction filter function.
    
    
    std::string indVField = pset.get<std::string>("IndVFieldShape");
    std::vector<double> indVFieldParams = pset.get<std::vector<double> >("IndVFieldParams");
    fIndVFieldFunc = new TF1("indVField", indVField.c_str());
    for(unsigned int i=0; i<indVFieldParams.size(); ++i)
      fIndVFieldFunc->SetParameter(i, indVFieldParams[i]);
    // Warning, last parameter needs to be multiplied by the FFTSize, in current version of the code,

  }

}


//----------------------------------------------------------------------
// Accessor for single-plane signal shaper.
const util::SignalShaping&
util::SignalShapingServiceT1034::SignalShaping(unsigned int channel) const
{

  //std::cout<<std::endl;
  //std::cout<<"SignalShaping"<<std::endl;
  //std::cout<<std::endl;
  if(!fInit)
    init();

  // Figure out plane type.
  art::ServiceHandle<geo::Geometry> geom;

  // Return appropriate shaper.

  geo::SigType_t sigtype = geom->SignalType(channel);
  if (sigtype == geo::kInduction)
     {
      return fIndVSignalShaping;
      
      
     }
  else if (sigtype == geo::kCollection)
      {
      return fColSignalShaping;
      }
  else
    {
    std::cout<<"I've got a problem"<<std::endl;
    std::cout<<std::endl;
    throw cet::exception("SignalShapingServiceT1034")<< "can't determine"
                                                          << " SignalType\n";
							  
    }  

  

  return fColSignalShaping;
}


//----------------------------------------------------------------------
// Initialization method.
// Here we do initialization that can't be done in the constructor.
// All public methods should ensure that this method is called as necessary.
void util::SignalShapingServiceT1034::init()
{

  //std::cout<<std::endl;
  //std::cout<<"Initalization"<<std::endl;
  //std::cout<<std::endl;
  
  if(!fInit) {
    fInit = true;
    
    // Do microboone-specific configuration of SignalShaping by providing
    // microboone response and filter functions.

    // Calculate field and electronics response functions.

    SetFieldResponse();
    SetElectResponse();

    // Configure convolution kernels.

    fColSignalShaping.AddResponseFunction(fColFieldResponse);
    fColSignalShaping.AddResponseFunction(fElectResponse);
    fColSignalShaping.SetPeakResponseTime(0.);

    fIndVSignalShaping.AddResponseFunction(fIndVFieldResponse);
    fIndVSignalShaping.AddResponseFunction(fElectResponse);
    fIndVSignalShaping.SetPeakResponseTime(0.);

    // Calculate filter functions.

    SetFilters();

    // Configure deconvolution kernels.
    
    
    fColSignalShaping.AddFilterFunction(fColFilter);
    
    fColSignalShaping.CalculateDeconvKernel();
    fIndVSignalShaping.AddFilterFunction(fIndVFilter);
    fIndVSignalShaping.CalculateDeconvKernel();
  }
}


//----------------------------------------------------------------------
// Calculate LArIAT field response.
//----------------------------------------------------------------------
void util::SignalShapingServiceT1034::SetFieldResponse()
{
  // Get services.

  //std::cout<<std::endl;
  //std::cout<<"SetFieldResponse"<<std::endl;
  //std::cout<<std::endl;
  
  // #############################
  // ### Load Geometry Service ###
  // ############################# 
  art::ServiceHandle<geo::Geometry> geo;
  
  // ################################
  // ### Load Detector Properties ###
  // ################################
  art::ServiceHandle<util::DetectorProperties> detprop;
  
  // ####################################
  // ### Load Liquid Argon Properties ###
  // ####################################
  art::ServiceHandle<util::LArProperties> larp;

  // #################################
  // ### Determine the plane pitch ###
  // #################################
  double xyz1[3] = {0.};
  double xyz2[3] = {0.};
  double xyzl[3] = {0.};
  // should always have at least 2 planes
  geo->Plane(0).LocalToWorld(xyzl, xyz1);
  geo->Plane(1).LocalToWorld(xyzl, xyz2);

  // ### this assumes all planes are equidistant from each other, ###
  // ###             probably not a bad assumption                ###
  double pitch = xyz2[0] - xyz1[0]; ///in cm
  
  // ############################################################
  // ### Resizing the Collection and Induction Field response ###
  // ############################################################
  fColFieldResponse.resize(fNFieldBins, 0.);
  std::cout<<"fColFieldResponse.size() = "<<fColFieldResponse.size()<<std::endl;
  fIndVFieldResponse.resize(fNFieldBins, 0.);

  // set the response for the collection plane first
  // the first entry is 0

  double driftvelocity=larp->DriftVelocity()/1000.;  
  double integral = 0.;
  
  // #################################################################
  // ### If we are going to use the function field shape then this ###
  // ###     variable must have been set true in the fcl file      ###
  // #################################################################
  if(fUseFunctionFieldShape) 
     {
     // ### Get the FFT service ###
     art::ServiceHandle<util::LArFFT> fft;
     // ### Set the size of the FFT, this should be some power 2^n such 
     // ### that the size is greater than the readout window (for LAriAT 3072)
     int signalSize = fft->FFTSize();
    
    std::cout<<"signalSize = "<<signalSize<<std::endl;  
    
    // ### Resize the Field Response vector to the new signal size ###
    fColFieldResponse.resize(signalSize, 0.);
    fIndVFieldResponse.resize(signalSize, 0.);
         
    // #################################################################	 
    // ### Evaluate the Collection/Induction Field Response Function ###
    // ################################################################# 
    for(int i = 0; i < signalSize; i++) 
       {
       //ramp[i]=;
     
       fColFieldResponse[i]=fColFieldFunc->Eval(i);
       //std::cout<<"fColFieldResponse[i] = "<<fColFieldResponse[i]<<std::endl;
      
       integral += fColFieldResponse[i];

       fIndVFieldResponse[i]=fIndVFieldFunc->Eval(i);
       //std::cout<<"fIndVFieldResponse[i] = "<<fIndVFieldResponse[i]<<std::endl;
       }//<---End i loop
     
    // ##########################################################################################
    // ### What does this do? I can't tell why you would do this so I comment it out (JAsaadi)###
    // ##########################################################################################
    for(int i = 0; i < signalSize; ++i)
       {
       fColFieldResponse[i] *= fColFieldRespAmp/integral;
       //std::cout<<"fColFieldResponse[i] = "<<fColFieldResponse[i]<<std::endl;
       
       }
      
    //this might be not necessary if the function definition is not defined in the middle of the signal range  
   // fft->ShiftData(fIndVFieldResponse,signalSize/2.0);

     }//<---End Using Function Field Shape 
    
    // #################################################################################
    // ### Instead of using a function you can always just define a vector of points ###
    // #################################################################################
    else if (fUseSimpleFieldShape) 
       {
       
       // ### Send a message that you are using the hard coded field shapes ###
       mf::LogInfo("SignalShapingServiceT1034") << " using try-2 hard-coded field shapes " ;
       
       // ########################
       // ### Collection Plane ###
       // ########################
       const int nbincPlane = 16;
       double cPlaneResponse[nbincPlane] = 
       {0,               0,               0,   0.02620087336,   0.02620087336, 
        0.04366812227,    0.1310043668,    0.1659388646,    0.1397379913,    0.3711790393, 
        0.06550218341,    0.0480349345,  -0.01310043668, -0.004366812227,               0, 
        0
       };
       
       // ----------------------------------------------
       // --- Again the weird normalization (JAsaadi)---
       // ----------------------------------------------
       for(int i = 1; i < nbincPlane; ++i)
         {
         fColFieldResponse[i] = cPlaneResponse[i];
         integral += fColFieldResponse[i];
         }

       for(int i = 0; i < nbincPlane; ++i)
         {
         //fColFieldResponse[i] *= fColFieldRespAmp/integral;
         fColFieldResponse[i] /= integral;
         }

       // #######################
       // ### Induction Plane ###
       // #######################
       const int nbinvPlane = 20;
       double vPlaneResponse[nbinvPlane] = 
       {0,               0,   0.01090909091,   0.01090909091,   0.01090909091, 
        0.02181818182,   0.03272727273,    0.7636363636,     2.018181818,            2.04, 
        1.090909091,     -1.03861518,    -1.757656458,    -1.757656458,    -1.118508655, 
        -0.2396804261,  -0.07989347537, -0.007989347537,               0,               0
       };
       // ----------------------------------------------
       // --- Again the weird normalization (JAsaadi)---
       // ----------------------------------------------
       for (int i = 0; i < nbinvPlane; ++i) 
         {
         //fIndVFieldResponse[i] = vPlaneResponse[i]*fIndVFieldRespAmp/(nbiniOld);
         fIndVFieldResponse[i] = vPlaneResponse[i]/integral;
         }

       }
    // ##################################################################
    // ### If you want to use the simple field shaping parameters and ###
    // ###     you set all the other things to false you get this     ###
    // ##################################################################
    else {

    //////////////////////////////////////////////////
    mf::LogInfo("SignalShapingServiceT1034") << " using the old field shape " ;
    //int nbinc = TMath::Nint(fCol3DCorrection*(std::abs(pitch))/(driftvelocity*detprop->SamplingRate())); ///number of bins //KP
    int nbinc = TMath::Nint(fCol3DCorrection*(std::abs(pitch))/(driftvelocity*128)); ///number of bins //KP
    
    double integral = 0.;
    for(int i = 1; i < nbinc; ++i){
      fColFieldResponse[i] = fColFieldResponse[i-1] + 1.0;
      integral += fColFieldResponse[i];
    }

    for(int i = 0; i < nbinc; ++i){
      fColFieldResponse[i] *= fColFieldRespAmp/integral;
    }

    // now the induction plane
    
    //int nbini = TMath::Nint(fInd3DCorrection*(std::abs(pitch))/(driftvelocity*detprop->SamplingRate()));//KP
    int nbini = TMath::Nint(fInd3DCorrection*(std::abs(pitch))/(driftvelocity*128));//KP
   

    for(int i = 0; i < nbini; ++i){
      fIndVFieldResponse[i] = fIndVFieldRespAmp/(1.*nbini);
      fIndVFieldResponse[nbini+i] = -fIndVFieldRespAmp/(1.*nbini);
    }

  }
  
  return;
}


//----------------------------------------------------------------------
// Calculate microboone field response.
void util::SignalShapingServiceT1034::SetElectResponse()
{
  //std::cout<<std::endl;
  //std::cout<<"SetElectResponse"<<std::endl;
  //std::cout<<std::endl;
  
  // Get services.
  
  art::ServiceHandle<geo::Geometry> geo;
  art::ServiceHandle<util::DetectorProperties> detprop;
  art::ServiceHandle<util::LArFFT> fft;

  LOG_DEBUG("SignalShapingT1034") << "Setting T1034 electronics response function...";

  int nticks = fft->FFTSize();
  fElectResponse.resize(nticks, 0.);
  std::vector<double> time(nticks,0.);

  //Gain and shaping time variables from fcl file:    
  double Ao = fShapeTimeConst[0];  //gain
  double To = fShapeTimeConst[1];  //peaking time
    
  // this is actually sampling time, in ns
  // mf::LogInfo("SignalShapingT1034") << "Check sampling intervals: " 
  //                                  << fSampleRate << " ns" 
  //                                  << "Check number of samples: " << fNTicks;

  // The following sets the microboone electronics response function in 
  // time-space. Function comes from BNL SPICE simulation of T1034 
  // electronics. SPICE gives the electronics transfer function in 
  // frequency-space. The inverse laplace transform of that function 
  // (in time-space) was calculated in Mathematica and is what is being 
  // used below. Parameters Ao and To are cumulative gain/timing parameters 
  // from the full (ASIC->Intermediate amp->Receiver->ADC) electronics chain. 
  // They have been adjusted to make the SPICE simulation to match the 
  // actual electronics response. Default params are Ao=1.4, To=0.5us. 
  double integral=0.;
  
  for(int i = 0; i < nticks; ++i){

    //convert time to microseconds, to match fElectResponse[i] definition
    //time[i] = (1.*i)*detprop->SamplingRate()*1e-3; 
    time[i] = (1.*i)*128*1e-3;
    fElectResponse[i] = 
      4.31054*exp(-2.94809*time[i]/To)*Ao - 2.6202*exp(-2.82833*time[i]/To)*cos(1.19361*time[i]/To)*Ao
      -2.6202*exp(-2.82833*time[i]/To)*cos(1.19361*time[i]/To)*cos(2.38722*time[i]/To)*Ao
      +0.464924*exp(-2.40318*time[i]/To)*cos(2.5928*time[i]/To)*Ao
      +0.464924*exp(-2.40318*time[i]/To)*cos(2.5928*time[i]/To)*cos(5.18561*time[i]/To)*Ao
      +0.762456*exp(-2.82833*time[i]/To)*sin(1.19361*time[i]/To)*Ao
      -0.762456*exp(-2.82833*time[i]/To)*cos(2.38722*time[i]/To)*sin(1.19361*time[i]/To)*Ao
      +0.762456*exp(-2.82833*time[i]/To)*cos(1.19361*time[i]/To)*sin(2.38722*time[i]/To)*Ao
      -2.6202*exp(-2.82833*time[i]/To)*sin(1.19361*time[i]/To)*sin(2.38722*time[i]/To)*Ao 
      -0.327684*exp(-2.40318*time[i]/To)*sin(2.5928*time[i]/To)*Ao + 
      +0.327684*exp(-2.40318*time[i]/To)*cos(5.18561*time[i]/To)*sin(2.5928*time[i]/To)*Ao
      -0.327684*exp(-2.40318*time[i]/To)*cos(2.5928*time[i]/To)*sin(5.18561*time[i]/To)*Ao
      +0.464924*exp(-2.40318*time[i]/To)*sin(2.5928*time[i]/To)*sin(5.18561*time[i]/To)*Ao;

      integral+=fElectResponse[i];
  }// end loop over time buckets
    
  LOG_DEBUG("SignalShapingT1034") << " Done.";

  // normalize fElectResponse[i], before the convolution
  // Put in overall normalization in a pedantic way:
  // first put in the pulse area per eleectron at the lowest gain setting,
  // then normalize by the actual ASIC gain setting used.
  // This code is executed only during initialization of service,
  // so don't worry about code inefficiencies here.
   for(int i = 0; i < nticks; ++i){
     fElectResponse[i] /= integral;
     fElectResponse[i] *= fADCTicksPerPCAtLowestASICGainSetting*1.6e-7;
     fElectResponse[i] *= fASICGainInMVPerFC/4.7;
   }
  
  
  return;

}


//----------------------------------------------------------------------
// Calculate microboone filter functions.
void util::SignalShapingServiceT1034::SetFilters()
{
  
  //std::cout<<std::endl;
  //std::cout<<"SetFilters"<<std::endl;
  //std::cout<<std::endl;
  // Get services.
  
  art::ServiceHandle<util::DetectorProperties> detprop;
  art::ServiceHandle<util::LArFFT> fft;

  //double ts = detprop->SamplingRate();
  double ts = 128;
  
  //std::cout<<"ts = "<<ts<<std::endl;
  int n = fft->FFTSize() / 2;

  // Calculate collection filter.

  fColFilter.resize(n+1);
  fIndVFilter.resize(n+1);
  
  if(!fGetFilterFromHisto)
  {
  fColFilterFunc->SetRange(0, double(n));

  for(int i=0; i<=n; ++i) {
    //double freq = 500. * i / (ts * n);      // Cycles / microsecond.
    double freq = 1250. * i / (ts * n);      // Cycles / microsecond.
    double f = fColFilterFunc->Eval(freq);
    fColFilter[i] = TComplex(f, 0.);
  }
  

  // Calculate induction filters.
 
 
   fIndVFilterFunc->SetRange(0, double(n));

  for(int i=0; i<=n; ++i) {
    //double freq = 350. * i / (ts * n);      // Cycles / microsecond.
    //double freq = 500. * i / (ts * n);      // Cycles / microsecond.
    double freq = 1000. * i / (ts * n);      // Cycles / microsecond.
    double f = fIndVFilterFunc->Eval(freq);
    fIndVFilter[i] = TComplex(f, 0.);
    }
  
  }
  else
  {
    
    for(int i=0; i<=n; ++i) {
      double f = fFilterHist[1]->GetBinContent(i);  // hardcoded plane numbers. Bad. but detector specific. Leaving for now.
      fColFilter[i] = TComplex(f, 0.);
      double g = fFilterHist[0]->GetBinContent(i);
      fIndVFilter[i] = TComplex(g, 0.);
  
      
    }
  }
  
  fIndVSignalShaping.AddFilterFunction(fIndVFilter);
  fColSignalShaping.AddFilterFunction(fColFilter);
  
}



namespace util {

  DEFINE_ART_SERVICE(SignalShapingServiceT1034)

}
