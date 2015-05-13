////////////////////////////////////////////////////////////////////////
//
// CalWireROI class - variant of CalWire that deconvolves in 
// Regions Of Interest
//
// Originally written by: baller@fnal.gov
//
// Adapted for LArIAT by: jasaadi@fnal.gov
////////////////////////////////////////////////////////////////////////

// ####################
// ### C++ Includes ###
// ####################
#include <string>
#include <vector>
#include <utility> // std::pair<>
#include <memory> // std::unique_ptr<>
#include <iomanip>

// ##########################
// ### Framework includes ###
// ##########################
#include "fhiclcpp/ParameterSet.h" 
#include "messagefacility/MessageLogger/MessageLogger.h" 
#include "art/Framework/Core/ModuleMacros.h" 
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h" 
#include "art/Framework/Principal/Handle.h" 
#include "art/Persistency/Common/Ptr.h" 
#include "art/Persistency/Common/PtrVector.h" 
#include "art/Framework/Services/Registry/ServiceHandle.h" 
#include "art/Framework/Services/Optional/TFileService.h" 
#include "art/Utilities/Exception.h"

// ########################
// ### LArSoft Includes ###
// ########################
#include "Utilities/SignalShapingServiceT1034.h"
#include "Geometry/Geometry.h"
#include "Filters/ChannelFilter.h"
#include "RawData/RawDigit.h"
#include "RawData/raw.h"
#include "RecoBase/Wire.h"
#include "Utilities/LArFFT.h"
#include "TComplex.h"
#include "TFile.h"
#include "TH2D.h"
#include "TF1.h"
#include "RecoBaseArt/WireCreator.h"
#include "Utilities/AssociationUtil.h"


// namespace for calwireROI
namespace caldata {
class CalWireROIT1034 : public art::EDProducer {

  public:
    
    // create calibrated signals on wires. this class runs 
    // an fft to remove the electronics shaping.     
    explicit CalWireROIT1034(fhicl::ParameterSet const& pset); 
    virtual ~CalWireROIT1034();
    
    void produce(art::Event& evt); 
    void beginJob(); 
    void endJob();                 
    void reconfigure(fhicl::ParameterSet const& p);
    
    void reconfFFT(int temp_fftsize);
  private:
    
    // === Module name the made the Digits ====
    std::string  fDigitModuleLabel;
    
    // ===  Spill name  ====
    std::string  fSpillName;  ///< nominal spill is an empty string
                              ///< it is set by the DigitModuleLabel
                              ///< ex.:  "daq:preSpill" for prespill data
    // === Threshold ADC counts for ROI === 
    std::vector<unsigned short> fThreshold; 
    
    // === Minimum width for ROI ===
    unsigned short fMinWid; 
    
    // === Minimum separation between ROI's ===       
    unsigned short fMinSep;
    
    // === Size of FFT for ROI deconvolution ===          
    int fFFTSize; 
    
    // === Pre timing ROI padding size ===
    std::vector<unsigned short> fPreROIPad; ///< ROI padding
    
    // === Post timing ROI padding size ===
    std::vector<unsigned short> fPostROIPad; ///< ROI padding
    
    // === Do baseline subtraction after deconvolution ====
    bool fDoBaselineSub; 
    
    // === Correct the induction plane response ====
    bool fuPlaneRamp;     ///< set true for correct U plane wire response
    int  fSaveWireWF;     ///< Save recob::wire object waveforms
    size_t fEventCount;  ///< count of event processed
    
    // ==== jasaadi: Come back and figure out what we need ======
    //void doDecon(std::vector<float>& holder, 
      //raw::ChannelID_t channel, unsigned int thePlane,
      //std::vector<std::pair<unsigned int, unsigned int>> rois,
      //std::vector<std::pair<unsigned int, unsigned int>> holderInfo,
      //recob::Wire::RegionsOfInterest_t& ROIVec,
      //art::ServiceHandle<util::SignalShapingServiceMicroBooNE>& sss);
    float SubtractBaseline(std::vector<float>& holder, float basePre,
			   float basePost,unsigned int roiStart,unsigned int roiLen,
			   unsigned int dataSize);

    //bool                               fDoBaselineSub_WaveformPropertiesAlg;
    //util::WaveformPropertiesAlg<float> fROIPropertiesAlg;
    //float SubtractBaseline(const std::vector<float>& holder);
    
  protected: 
    
  }; // class CalWireROIT1034

  DEFINE_ART_MODULE(CalWireROIT1034)


//--------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------
CalWireROIT1034::CalWireROIT1034(fhicl::ParameterSet const& pset)
  {
    fSpillName="";
    this->reconfigure(pset);
    if(fSpillName.size()<1) {produces< std::vector<recob::Wire> >();
     produces<art::Assns<raw::RawDigit, recob::Wire>>();
    }
    else { produces< std::vector<recob::Wire> >(fSpillName);
          produces<art::Assns<raw::RawDigit, recob::Wire>>(fSpillName);
    }
  }//<---End CalWireROIT1034 pset
  
//--------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------  
CalWireROIT1034::~CalWireROIT1034()
  {
  }






//--------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------- 
void CalWireROIT1034::reconfigure(fhicl::ParameterSet const& p)
  {
  
    std::vector<unsigned short> vin;
    std::vector<unsigned short> zin;
    
    fDigitModuleLabel = p.get< std::string >          ("DigitModuleLabel", "daq");
    fThreshold        = p.get< std::vector<unsigned short> >   ("Threshold");
    fMinWid           = p.get< unsigned short >       ("MinWid");
    fMinSep           = p.get< unsigned short >       ("MinSep");
    //uin               = p.get< std::vector<unsigned short> >   ("uPlaneROIPad");
    vin               = p.get< std::vector<unsigned short> >   ("vPlaneROIPad");
    zin               = p.get< std::vector<unsigned short> >   ("zPlaneROIPad");
    fDoBaselineSub    = p.get< bool >                 ("DoBaselineSub");
    fuPlaneRamp       = p.get< bool >                 ("uPlaneRamp");
    fFFTSize          = p.get< int  >                 ("FFTSize");
    //fSaveWireWF       = p.get< int >                  ("SaveWireWF");

    //fDoBaselineSub_WaveformPropertiesAlg = p.get< bool >("DoBaselineSub_WaveformPropertiesAlg");
    
    //fPedestalRetrievalAlg.Reconfigure(p.get<fhicl::ParameterSet>("DetPedestalRetrievalAlg"));
    
    if( vin.size() != 2 || zin.size() != 2) {
      throw art::Exception(art::errors::Configuration)
        << "u/v/z plane ROI pad size != 2";
    }

    fPreROIPad.resize(3);
    fPostROIPad.resize(3);
    
    // put the ROI pad sizes into more convenient vectors
    fPreROIPad[0]  = vin[0];
    fPostROIPad[0] = vin[1];
    fPreROIPad[1]  = zin[0];
    fPostROIPad[1] = zin[1];
    //fPreROIPad[2]  = zin[0];
    //fPostROIPad[2] = zin[1];
    
    fSpillName.clear();
    
    size_t pos = fDigitModuleLabel.find(":");
    if( pos!=std::string::npos ) {
      fSpillName = fDigitModuleLabel.substr( pos+1 );
      fDigitModuleLabel = fDigitModuleLabel.substr( 0, pos );
    }
    
    
  }//<---End CalWireROIT1034 reconfigure

  
//--------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------- 
void CalWireROIT1034::beginJob()
  {  

  }

//--------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------- 
void CalWireROIT1034::endJob()
  {  

  }


//--------------------------------------------------------------------------------------------
//					CalWire ROI Event Loop
//-------------------------------------------------------------------------------------------- 
void CalWireROIT1034::produce(art::Event& evt)
  {
  // ####################################
  // ### Loading the geometry service ###
  // ####################################
  art::ServiceHandle<geo::Geometry> geom;
  
  // ###############################
  // ### Loading the FFT service ###
  // ###############################
  art::ServiceHandle<util::LArFFT> fFFT;
  //int transformSize = fFFT->FFTSize();
  
  // ##########################################
  // ### Loading the Signal Shaping Service ###
  // ##########################################
  art::ServiceHandle<util::SignalShapingServiceT1034> sss;
    
  // #########################################
  // ### Make a collection of recob::Wires ###
  // #########################################
  std::unique_ptr<std::vector<recob::Wire> > wirecol(new std::vector<recob::Wire>);
  
  // ##########################################################
  // ### Making an association of wire digits to raw digits ###
  // ##########################################################
  std::unique_ptr<art::Assns<raw::RawDigit,recob::Wire> > WireDigitAssn
    (new art::Assns<raw::RawDigit,recob::Wire>);

  // ####################################
  // ### Read in the RawDigit objects ###
  // ####################################
  art::Handle< std::vector<raw::RawDigit> > digitVecHandle;
  
  // ###########################
  // ### RawDigit getByLabel ###
  // ###########################
  if(fSpillName.size()>0) evt.getByLabel(fDigitModuleLabel, fSpillName, digitVecHandle);
  else evt.getByLabel(fDigitModuleLabel, digitVecHandle);
  
  if (!digitVecHandle->size())  return;
  mf::LogInfo("CalWireROIT1034") << "CalWireT1034:: digitVecHandle size is " << digitVecHandle->size();
  
  // #######################################
  // ### Getting the Dead Channel Filter ###
  // #######################################
  filter::ChannelFilter *chanFilt = new filter::ChannelFilter();
  
  // ##################################################
  // ### Getting the normalization of deconvolution ###
  // ##################################################
  // jasaadi: fix this
  double DeconNorm = 1; //sss->GetDeconNorm();
  
  // ###########################
  // ### Loop over the wires ###
  // ###########################
  wirecol->reserve(digitVecHandle->size());
  
  for(size_t rdIter = 0; rdIter < digitVecHandle->size(); ++rdIter)
     {
     
     // ### Initializing the Channel Number ###
     uint32_t     channel(0); // channel number
     
     // ### Time bin loop variable ###
     unsigned int bin(0);
      
     // vector that will be moved into the Wire object
     recob::Wire::RegionsOfInterest_t ROIVec;
      
     // the starting position and length of each ROI in the packed holder vector
     std::vector<std::pair<unsigned int, unsigned int>> holderInfo;
     // vector of ROI begin and end bins
     std::vector<std::pair<unsigned int, unsigned int>> rois;
      
     // get the reference to the current raw::RawDigit
     art::Ptr<raw::RawDigit> digitVec(digitVecHandle, rdIter);
     channel = digitVec->Channel();
     unsigned int dataSize = digitVec->Samples();
     // vector holding uncompressed adc values
     std::vector<short> rawadc(dataSize);
      
     std::vector<geo::WireID> wids = geom->ChannelToWire(channel);
     unsigned int thePlane = wids[0].Plane;

     // skip bad channels
     if(!chanFilt->BadChannel(channel)) {
        
        // uncompress the data
        raw::Uncompress(digitVec->ADCs(), rawadc, digitVec->Compression());
        // loop over all adc values and subtract the pedestal
	// When we have a pedestal database, can provide the digit timestamp as the third argument of GetPedestalMean
	
	// jasaadi: Need to get the pedestal somehow
        float pdstl = 0;//fPedestalRetrievalAlg.PedMean(channel);
	
	//subtract time-offset added in SimWireMicroBooNE_module
	// Xin, remove the time offset
	//int time_offset = 0.;//sss->FieldResponseTOffset(channel);
        unsigned int roiStart = 0;
	
	// jasaadi: need to get the noise
	double raw_noise = 1;//sss->GetRawNoise(channel);

        // search for ROIs
        for(bin = 1; bin < dataSize; ++bin) {
          float SigVal = fabs(rawadc[bin] - pdstl);
          if(roiStart == 0) {
            // not in a ROI
            // Handle the onset of a ROI differently for the 1st induction plane
            // if it has been modeled correctly
            // if(fuPlaneRamp && thePlane == 0) {
            //   if(rawadc[bin] - rawadc[bin - 1] < sPed) roiStart = bin;
            // } else {
            //   if(SigVal > fThreshold[thePlane]) roiStart = bin;
            // }
	    unsigned int sbin[7];
	    if (bin>=3) {
	      sbin[0] = bin -3;
	      sbin[1] = bin -2;
	      sbin[2] = bin -1;
	    }else if (bin>=2){
	      sbin[0] = 0;
	      sbin[1] = bin-2;
	      sbin[2] = bin-1;
	    }else if (bin>=1){
	      sbin[0] =0;
	      sbin[1] =0;
	      sbin[2] = bin-1;
	    }else if (bin==0) {
	      sbin[0] =0;
	      sbin[1] =0;
	      sbin[2] =0;
	    }
	    sbin[3] = bin ; 
	    sbin[4] = bin + 1; if (sbin[4]>dataSize-1) sbin[4] =dataSize-1;
	    sbin[5] = bin + 2; if (sbin[5]>dataSize-1) sbin[5] =dataSize-1;
	    sbin[6] = bin + 3; if (sbin[6]>dataSize-1) sbin[6] =dataSize-1;
	    float sum = 0;
	    for (int qx = 0; qx!=7;qx++){
	      sum += rawadc[sbin[qx]]-pdstl;
	    }
	    sum = fabs(sum);
	    //std::cout << bin << " " << sum << " " << raw_noise/sqrt(7.)*3. << std::endl;
	    if (sum > raw_noise*sqrt(7.)*6.) roiStart = bin;

	  } else {
            // leaving a ROI?
            if(SigVal < fThreshold[thePlane]) {
              // is the ROI wide enough?
              //unsigned int roiLen = bin - roiStart;
              // if(roiLen > transformSize) {
              //   mf::LogError("CalWireROI")<<"ROI too long "
              //     <<roiLen<<" on plane:wire "<<thePlane<<":"<<theWire;
              //   break;
              // }
              //if(roiLen > fMinWid && roiLen < transformSize) 
	      //std::cout << roiStart << " " << roiLen << std::endl;
	      // if(roiLen > fMinWid) 
	      rois.push_back(std::make_pair(roiStart, bin));
              roiStart = 0;
            }
          } // roiStart test
        } // bin
	// add the last ROI if existed
	if (roiStart!=0){
	  //unsigned int roiLen = dataSize -1 - roiStart;
	   // if(roiLen > fMinWid) 
	   rois.push_back(std::make_pair(roiStart, dataSize-1));
	   roiStart = 0;
	}

	

        // skip deconvolution if there are no ROIs
        if(rois.size() == 0) continue;

        holderInfo.clear();
	//        for(unsigned int bin = 0; bin < holder.size(); ++bin) holder[bin] = 0;
        
	 // pad the ROIs
        for(unsigned int ii = 0; ii < rois.size(); ++ii) {
          // low ROI end
          int low = rois[ii].first - fPreROIPad[thePlane];
          if(low < 0) low = 0;
          rois[ii].first = low;
          // high ROI end
          unsigned int high = rois[ii].second + fPostROIPad[thePlane];
          if(high >= dataSize) high = dataSize-1;
          rois[ii].second = high;
	  
        }

	// if (channel==3218){
	//	std::cout << "Xin " << " " << channel << " " << rois.size() << std::endl;
	//   for(unsigned int ii = 0; ii < rois.size(); ++ii) {
	//     std::cout << rois[ii].first << " " << rois[ii].second << std::endl;
	//   }
	// }

        // merge the ROIs?
        if(rois.size() > 1) {
          // temporary vector for merged ROIs
          std::vector<std::pair<unsigned int, unsigned int>> trois;
          
	  for (unsigned int ii = 0; ii<rois.size();ii++){
	    unsigned int roiStart = rois[ii].first;
	    unsigned int roiEnd = rois[ii].second;

	    // if (channel==806)
	    //   std::cout << "a" << " " << roiStart << " " << roiEnd << std::endl;

	    int flag1 = 1;
	    unsigned int jj=ii+1;
	    while(flag1){	
	      if (jj<rois.size()){
		if(rois[jj].first <= roiEnd  ) {
		  roiEnd = rois[jj].second;
		  ii = jj;
		  jj = ii+1;
		}else{
		  flag1 = 0;
		}
	      }else{
		flag1 = 0;
	      }
	    }
	    
	   // if (channel==806)
	   //   std::cout << "b" << " " << roiStart << " " << roiEnd << std::endl;
	 
	    trois.push_back(std::make_pair(roiStart,roiEnd));	    
	  }
	  
	  rois = trois;
	}
	  
	

	for (unsigned int ir = 0; ir < rois.size(); ++ir) {
	  unsigned int roiLen = rois[ir].second - rois[ir].first + 1;
	  unsigned int roiStart = rois[ir].first;
	  //treat FFT Size
	  // if (channel==806)
	  //   std::cout << roiStart << " " << roiLen << std::endl;
	  

	  int flag =1;
	  float tempPre=0,tempPost=0;
	  std::vector<float> holder;
	  while(flag){
	    
	    unsigned int transformSize = fFFTSize; //current transformsize
	    //if ROI length is longer, take ROI length
	    if (roiLen > transformSize) transformSize = roiLen;
	    
	    // Get signal shaping service.
	    
	    // jasaadi: fix this
	    //sss->SetDecon(transformSize);
	    
	    transformSize = fFFT->FFTSize();
	    // temporary vector of signals
	    holder.resize(transformSize,0);
	    
	    unsigned int hBin = 0;
	    for(unsigned int bin = roiStart; bin < roiStart + holder.size(); ++bin) {
	      if (bin < dataSize){
		holder[hBin] = rawadc[bin]-pdstl;
	      }else{
		holder[hBin] = rawadc[bin-dataSize]-pdstl;
	      }
	      if (bin>=dataSize-1) flag = 0;
	      ++hBin;
	    } // bin

	    //std::cout << channel << " " << roiStart << " " << flag << " " << dataSize << " " << holder.size() << " " << roiLen << std::endl;

	    sss->Deconvolute(channel,holder);
	    for(bin = 0; bin < holder.size(); ++bin) holder[bin]=holder[bin]/DeconNorm;

	    // if (channel==3218){
	    //   for(unsigned int bin = 0; bin <holder.size(); ++bin) {
	    // 	std::cout << bin << " " <<  holder[bin] << std::endl;
	    //   }
	    // }
	    //1. Check Baseline match?
	    // If not, include next ROI(if none, go to the end of signal)
	    // If yes, proceed
	    tempPre=0,tempPost=0;
	    for(unsigned int bin = 0; bin < 20; ++bin) {
	      tempPre  += holder[bin];
	      tempPost += holder[roiLen -1 - bin];
	    }
	    tempPre = tempPre/20.;
	    tempPost = tempPost/20.;
	    
	    // ### jasaadi: fix this
	    double deconNoise = 2;//sss->GetDeconNoise(channel)/sqrt(10.)*4;
	    
	    if (fabs(tempPost-tempPre)<deconNoise){
	      flag = 0;
	    }else{
	      if (tempPre > tempPost && roiStart <= 2){
		//if (tempPre > tempPost){
		flag = 0;
	      }else{
		ir++;
		if (ir<rois.size()){
		  roiLen += 100;
		  if (roiLen >= rois[ir].first - roiStart + 1)
		    roiLen = rois[ir].second - roiStart + 1;
		}else{
		  roiLen += 100;
		  if (roiLen>dataSize-roiStart)
		    roiLen = dataSize - roiStart;
		}
	      }
	    }
	  }
	  


	  // transfer the ROI and start bins into the vector that will be
	  // put into the event
	  std::vector<float> sigTemp;
	  unsigned int bBegin = 0;
	  //unsigned int theROI =ir;
	  unsigned int bEnd = bBegin + roiLen;
	  float basePre = 0., basePost = 0.;

	  float base=0;
	  if(fDoBaselineSub && fPreROIPad[thePlane] > 0 ) {
	    basePre =tempPre;
	    basePost=tempPost;
	    
	    

	    base = SubtractBaseline(holder, basePre,basePost,roiStart,roiLen,dataSize);
	    // if (channel==200)
	    //   std::cout << basePre << " " << basePost << " " << roiStart << " " << roiLen << " " << dataSize << " " << base << std::endl;
	  } // fDoBaselineSub ...
	  
	  
	// jasaadi: FIX ME  
	//  else if(fDoBaselineSub_WaveformPropertiesAlg)
      //{
       //   holder.resize(roiLen);
	  // jasaadi: need to fix this
          //base = fROIPropertiesAlg.GetWaveformPedestal(holder);
     // }


	  for(unsigned int jj = bBegin; jj < bEnd; ++jj) {
	    sigTemp.push_back(holder[jj]-base);
	  } // jj
	        
	  // add the range into ROIVec 
	  ROIVec.add_range(roiStart, std::move(sigTemp));
	  
	}
	




      } // end if not a bad channel 

      // create the new wire directly in wirecol
      wirecol->push_back(recob::WireCreator(std::move(ROIVec),*digitVec).move());
      // add an association between the last object in wirecol
      // (that we just inserted) and digitVec
      if (!util::CreateAssn(*this, evt, *wirecol, digitVec, *WireDigitAssn, fSpillName)) {
        throw art::Exception(art::errors::InsertFailure)
          << "Can't associate wire #" << (wirecol->size() - 1)
          << " with raw digit #" << digitVec.key();
      } // if failed to add association
    //  DumpWire(wirecol->back()); // for debugging
    }
  
  }//<---End evt loop 

// ------------------------------------------------------------------------------------
//				Subtract Baseline
//-------------------------------------------------------------------------------------
float CalWireROIT1034::SubtractBaseline(std::vector<float>& holder, float basePre,
				     float basePost,unsigned int roiStart,
				     unsigned int roiLen,unsigned int dataSize)
  {
    float base=0;

    //can not trust the early part
    if (roiStart < 20 && roiStart + roiLen < dataSize - 20){
      base = basePost;
      // can not trust the later part
    }else if (roiStart >= 20 && roiStart + roiLen >= dataSize - 20){
      base = basePre;
      // can trust both
    }else if (roiStart >= 20 && roiStart + roiLen < dataSize - 20){
      if (fabs(basePre-basePost)<3){
	base = (basePre+basePost)/2.;
      }else{
	if (basePre < basePost){
	  base = basePre;
	}else{
	  base = basePost;
	}
      }
      // can not use both
    }else{
      float min = 0,max=0;
      for (unsigned int bin = 0; bin < roiLen; bin++){
	if (holder[bin] > max) max = holder[bin];
	if (holder[bin] < min) min = holder[bin];
      }
      int nbin = max - min;
      if (nbin!=0){
	TH1F *h1 = new TH1F("h1","h1",nbin,min,max);
	for (unsigned int bin = 0; bin < roiLen; bin++){
	  h1->Fill(holder[bin]);
	}
	float ped = h1->GetMaximum();
	float ave=0,ncount = 0;
	for (unsigned int bin = 0; bin < roiLen; bin++){
	  if (fabs(holder[bin]-ped)<2){
	    ave +=holder[bin];
	    ncount ++;
	  }
	}
	if (ncount==0) ncount=1;
	ave = ave/ncount;
	h1->Delete();
	base = ave;
      }
    }
    
   
    return base;
  }

}//<---End caldata namespace
