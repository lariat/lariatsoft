////////////////////////////////////////////////////////////////////////
//
// CalWireROTT1034 class
//
// brebel@fnal.gov
//
// 11-3-09 Pulled all FFT code out and put into Utilitiess/LArFFT
//
////////////////////////////////////////////////////////////////////////

// C/C++ standard libraries
#include <string>
#include <vector>
#include <utility> // std::move()
#include <memory> // std::unique_ptr<>
// BB temp
#include <fstream>

// ROOT libraries
#include "TComplex.h"

// framework libraries
#include "cetlib/exception.h"
#include "cetlib/search_path.h"
#include "fhiclcpp/ParameterSet.h" 
#include "messagefacility/MessageLogger/MessageLogger.h" 
#include "art/Framework/Core/ModuleMacros.h" 
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h" 
#include "art/Framework/Principal/Handle.h" 
#include "canvas/Persistency/Common/Ptr.h" 
#include "art/Framework/Services/Registry/ServiceHandle.h" 

// LArSoft libraries
#include "larcoreobj/SimpleTypesAndConstants/RawTypes.h" // raw::ChannelID_t
#include "larcore/Geometry/Geometry.h"
#include "larevt/Filters/ChannelFilter.h"
#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RawData/raw.h"
#include "lardataobj/RecoBase/Wire.h"
#include "lardata/ArtDataHelper/WireCreator.h"
#include "lardata/Utilities/LArFFT.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "Utilities/SignalShapingServiceT1034.h"


///creation of calibrated signals on wires
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
    
  private:
    
    int          fDataSize;           ///< size of raw data on one wire
    int          fPostsample;         ///< number of postsample bins
    bool         fDoBaselineSub;      ///< option to perform baseline subtraction using postsample bins
    bool         fAdvancedBaselineSub;///< more sophisticated baseline subtraction (but slower)
    bool         fDoROI;              ///< make ROIs
    std::string  fDigitModuleLabel;   ///< module that made digits
    std::string  fSpillName;          ///< nominal spill is an empty string
                                      ///< it is set by the DigitModuleLabel
                                      ///< ex.:  "daq:preSpill" for prespill data
    unsigned short fPreROIPad;        ///< ROI padding
    unsigned short fPostROIPad;       ///< ROI padding

    int          fSampPrecision;      ///< Limit on number of decimal places for each sample in deconvoluted
                                      ///  signals (recob::Wires) to save disk space post-compression. For example,
                                      ///    = -1 --> 1.35718291... (floating point precision)
                                      ///    = 0  --> 1.
                                      ///    = 1  --> 1.4 (default)
                                      ///    = 2  --> 1.36

    bool                        fDodQdxCalib;          ///< Do we apply wire-by-wire calibration?
    std::string                 fdQdxCalibFileName;    ///< Text file for constants to do wire-by-wire calibration
    std::map<unsigned int, float> fdQdxCalib;          ///< Map to do wire-by-wire calibration, key is channel number, content is correction factor

    void SubtractBaseline(std::vector<float>& holder);    ///< basic baseline subraction function (using postsample bins)
    void SubtractBaselineAdv(std::vector<float>& holder); ///< advanced baseline subtraction function

    
  protected: 
    
  }; // class CalWireROIT1034

  DEFINE_ART_MODULE(CalWireROIT1034)
  
  //-------------------------------------------------
  CalWireROIT1034::CalWireROIT1034(fhicl::ParameterSet const& pset)
  {
    this->reconfigure(pset);

    produces< std::vector<recob::Wire> >(fSpillName);
    produces<art::Assns<raw::RawDigit, recob::Wire>>(fSpillName);
  }
  
  //-------------------------------------------------
  CalWireROIT1034::~CalWireROIT1034()
  {
  }

  //////////////////////////////////////////////////////
  void CalWireROIT1034::reconfigure(fhicl::ParameterSet const& p)
  {
    std::vector<unsigned short> uin;    std::vector<unsigned short> vin;
    std::vector<unsigned short> zin;


    fDigitModuleLabel       = p.get< std::string >("DigitModuleLabel",    "daq");
    fPostsample             = p.get< int >        ("PostsampleBins",      300);
    fDoBaselineSub          = p.get< bool >       ("DoBaselineSub",       true);
    fAdvancedBaselineSub    = p.get< bool >       ("AdvancedBaselineSub", false);
    fDoROI                  = p.get< bool >       ("DoROI",               false);
    fSampPrecision          = p.get< int >        ("SampPrecision",       1);
    uin                     = p.get< std::vector<unsigned short> >   ("PlaneROIPad");
     
    
    // put the ROI pad sizes into more convenient vectors
    fPreROIPad  = uin[0];
    fPostROIPad = uin[1];
    

    fSpillName.clear();
    
    size_t pos = fDigitModuleLabel.find(":");
    if( pos!=std::string::npos ) {
      fSpillName = fDigitModuleLabel.substr( pos+1 );
      fDigitModuleLabel = fDigitModuleLabel.substr( 0, pos );
    }

    //wire-by-wire calibration
    fDodQdxCalib        = p.get< bool >                          ("DodQdxCalib", false);
    if (fDodQdxCalib){
      fdQdxCalibFileName = p.get< std::string >                  ("dQdxCalibFileName");
      std::string fullname;
      cet::search_path sp("FW_SEARCH_PATH");
      sp.find_file(fdQdxCalibFileName, fullname);

      if (fullname.empty()) {
	std::cout << "Input file " << fdQdxCalibFileName << " not found" << std::endl;
	throw cet::exception("File not found");
      }
      else
	std::cout << "Applying wire-by-wire calibration using file " << fdQdxCalibFileName << std::endl;

      std::ifstream inFile(fullname, std::ios::in);
      std::string line;
      
      while (std::getline(inFile,line)) {
	unsigned int channel;
	float        constant;
	std::stringstream linestream(line);
	linestream >> channel >> constant;
	fdQdxCalib[channel] = constant;
	if (channel%100==0) std::cout<<"Channel "<<channel<<" correction factor "<<fdQdxCalib[channel]<<std::endl;
      }
    }

  }

  //-------------------------------------------------
  void CalWireROIT1034::beginJob()
  {  
  }

  //////////////////////////////////////////////////////
  void CalWireROIT1034::endJob()
  {  
  }
  
  //////////////////////////////////////////////////////
  void CalWireROIT1034::produce(art::Event& evt)
  {      
    // get the geometry
    art::ServiceHandle<geo::Geometry> geom;

    // get the FFT service to have access to the FFT size
    art::ServiceHandle<util::LArFFT> fFFT;
    int transformSize = fFFT->FFTSize();

    // Get signal shaping service.
    art::ServiceHandle<util::SignalShapingServiceT1034> sss;
    double DeconNorm = sss->GetDeconNorm();
 
    // make a collection of Wires
    std::unique_ptr<std::vector<recob::Wire> > wirecol(new std::vector<recob::Wire>);
    // ... and an association set
    std::unique_ptr<art::Assns<raw::RawDigit,recob::Wire> > WireDigitAssn
      (new art::Assns<raw::RawDigit,recob::Wire>);
    
    // Read in the digit List object(s). 
    art::Handle< std::vector<raw::RawDigit> > digitVecHandle;
    evt.getByLabel(fDigitModuleLabel, fSpillName, digitVecHandle);

    mf::LogInfo("CalWireROIT1034") << "CalWireROIT1034:: digitVecHandle size is " << digitVecHandle->size();

    // Use the handle to get a particular (0th) element of collection.
    art::Ptr<raw::RawDigit> digitVec0(digitVecHandle, 0);
        
    //  if (digitVec0->Compression() != raw::kZeroSuppression) {
    //   throw art::Exception(art::errors::UnimplementedFeature)
    //	<< "CalGausHFLBNE only supports zero-suppressed raw digit input!";
    //} // if

    unsigned int dataSize = 0; //size of raw data vectors

    raw::ChannelID_t channel = raw::InvalidChannelID; // channel number
    unsigned int bin(0);     // time bin loop variable
    
    filter::ChannelFilter chanFilt;
    
    std::vector<float> holder;                // holds signal data
    std::vector<short> rawadc(transformSize);  // vector holding uncompressed adc values
    std::vector<TComplex> freqHolder(transformSize+1); // temporary frequency data
    

    // loop over all wires
    wirecol->reserve(digitVecHandle->size());
    for(size_t rdIter = 0; rdIter < digitVecHandle->size(); ++rdIter){ // ++ move
      holder.clear();
      
      // get the reference to the current raw::RawDigit
      art::Ptr<raw::RawDigit> digitVec(digitVecHandle, rdIter);
      channel = digitVec->Channel();
      
      // skip bad channels
      if(!chanFilt.BadChannel(channel)) {
        holder.resize(transformSize);
        
	dataSize = digitVec->Samples();
        // uncompress the data
        raw::Uncompress(digitVec->ADCs(), rawadc, digitVec->Compression());
        
        // loop over all adc values. Pedestal is already subtracted
        for(bin = 0; bin < dataSize; ++bin) holder[bin]=rawadc[bin];
        // pad with zeros
        for (bin = dataSize; bin < holder.size(); ++bin) holder[bin] = 0;
        
        // Do deconvolution.
        sss->Deconvolute(channel, holder);
        for(bin = 0; bin < holder.size(); ++bin) holder[bin]=holder[bin]/DeconNorm;
      } // end if not a bad channel
      
      //This restores the DC component to signal removed by the deconvolution.
      if( fDoBaselineSub ){
        SubtractBaseline(holder);
        if( fAdvancedBaselineSub ) SubtractBaselineAdv(holder);
      } 


      // retrieve the wire-by-wire correction factor for this channel
      float constant = 1.0;
      if( fDodQdxCalib ) {
        if( fdQdxCalib.find(channel) != fdQdxCalib.end() )
	  constant = fdQdxCalib[channel];
      }

      // loop over the samples after deconvolution
      for(size_t iholder = 0; iholder < holder.size(); ++iholder) {
        
        // apply wire-by-wire calibration correction
	holder[iholder] *= constant;
        
        // limit precision by truncating digits to save disk space
        if( fSampPrecision >= 0 )
          holder[iholder] = roundf( holder[iholder] * pow(10,fSampPrecision) ) / pow(10,fSampPrecision);
      
      }


      if (fDoROI){
        // work out the ROI
        recob::Wire::RegionsOfInterest_t ROIVec;
        std::vector<std::pair<unsigned int, unsigned int>> holderInfo;
        std::vector<std::pair<unsigned int, unsigned int>> rois;
        
        double max = 0;
        double deconNoise = sss->GetDeconNoise(channel);
        // find out all ROI
        unsigned int roiStart = 0;
        for(bin = 0; bin < dataSize; ++bin) {
          double SigVal = holder[bin];
          if (SigVal > max) max = SigVal;
          if(roiStart == 0) {
            if (SigVal > 3*deconNoise) roiStart = bin; // 3 sigma above noise
          }else{
            if (SigVal < deconNoise){
              rois.push_back(std::make_pair(roiStart, bin));
              roiStart = 0;
            }
          }
        }
        if (roiStart!=0){
          rois.push_back(std::make_pair(roiStart, dataSize-1));
          roiStart = 0;
        }
	
        if(rois.size() == 0) continue;
        holderInfo.clear();
        for(unsigned int ii = 0; ii < rois.size(); ++ii) {
          // low ROI end
          int low = rois[ii].first - fPreROIPad;
          if(low < 0) low = 0;
          rois[ii].first = low;
          // high ROI end
          unsigned int high = rois[ii].second + fPostROIPad;
          if(high >= dataSize) high = dataSize-1;
          rois[ii].second = high;
          
        }
        // merge them
        if(rois.size() >= 1) {
          // temporary vector for merged ROIs
          
          for (unsigned int ii = 0; ii<rois.size();ii++){
            unsigned int roiStart = rois[ii].first;
            unsigned int roiEnd = rois[ii].second;
            
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
            std::vector<float> sigTemp;
            for(unsigned int kk = roiStart; kk < roiEnd; ++kk) {
              sigTemp.push_back(holder[kk]);
            } // jj
            ROIVec.add_range(roiStart, std::move(sigTemp));
         }
        }
	
        // save them
        wirecol->push_back(recob::WireCreator(std::move(ROIVec),*digitVec).move());
        
        // add an association between the last object in wirecol
        // (that we just inserted) and digitVec
        if (!util::CreateAssn(*this, evt, *wirecol, digitVec, *WireDigitAssn, fSpillName)) {
          throw art::Exception(art::errors::ProductRegistrationFailure)
          << "Can't associate wire #" << (wirecol->size() - 1)
          << " with raw digit #" << digitVec.key();
        } // if failed to add association
      }
      else{
        wirecol->push_back(recob::WireCreator(holder,*digitVec).move());
        // add an association between the last object in wirecol
        // (that we just inserted) and digitVec
        if (!util::CreateAssn(*this, evt, *wirecol, digitVec, *WireDigitAssn, fSpillName)) {
          throw art::Exception(art::errors::ProductRegistrationFailure)
          << "Can't associate wire #" << (wirecol->size() - 1)
          << " with raw digit #" << digitVec.key();
        } // if failed to add association
      }
    }
    
    if(wirecol->size() == 0)
      mf::LogWarning("CalWireROIT1034") << "No wires made for this event.";
    
    evt.put(std::move(wirecol), fSpillName);
    evt.put(std::move(WireDigitAssn), fSpillName);
    
    return;
  }



  //========================================================================
  // Basic baseline subtraction using postsample bins 
  void CalWireROIT1034::SubtractBaseline(std::vector<float>& holder)
  {
    if( fPostsample > 0 ) {
      double average=0.0;
      double sum=0.0;
      int n=0;
      for(size_t i = 0; i < (size_t)fPostsample; i++){
        size_t bin = holder.size()-i; 
        if( fabs(holder[bin]) < 20 ){ // avoid outliers
          sum+=holder[bin];
          n++;
        }
      }
      if(n) average=sum/n;
      for(size_t i=0; i < holder.size(); ++i) holder[i]-=average;
    }
  }
  
  //========================================================================
  // Advanced baseline subtraction using the whole waveform. A histogram
  // of all ADC values is created and a Gaussian fit is used to find the mean.
  void CalWireROIT1034::SubtractBaselineAdv(std::vector<float>& holder)
  {
    float min = 0,max=0;
    for (unsigned int bin = 0; bin < holder.size(); bin++){
      if (holder[bin] > max) max = holder[bin];
      if (holder[bin] < min) min = holder[bin];
    }
    int nbin = max - min;
    if (nbin!=0){
      TH1F h1("h1","h1",nbin,min,max);
      for (unsigned int bin = 0; bin < holder.size(); bin++) 
        h1.Fill(holder[bin]);
      float ped = h1.GetMaximum(); // mode
      float rms = h1.GetRMS();
      float ave=0,ncount = 0;
      for (unsigned int bin = 0; bin < holder.size(); bin++){
        if (fabs(holder[bin]-ped)<rms*3.){
          ave +=holder[bin];
          ncount ++;
        }
      }
      if (ncount==0) ncount=1;
      ave = ave/ncount;
      for (unsigned int bin = 0; bin < holder.size(); bin++){
        holder[bin] -= ave;
      }
    }
  }


} // end namespace caldata


