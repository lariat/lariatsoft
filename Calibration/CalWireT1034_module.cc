////////////////////////////////////////////////////////////////////////
//
// CalWireT1034 class
//
// brebel@fnal.gov
//
// 11-3-09 Pulled all FFT code out and put into Utilitiess/LArFFT
//  copied over to 1034 - andrzej.szelc@yale.edu
////////////////////////////////////////////////////////////////////////

#include <string>
#include <vector>
#include <stdint.h>
extern "C" {
#include <sys/types.h>
#include <sys/stat.h>
}
#include "art/Framework/Core/ModuleMacros.h" 
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h" 
#include "art/Framework/Principal/Handle.h" 
#include "art/Persistency/Common/Ptr.h" 
#include "art/Persistency/Common/PtrVector.h" 
#include "art/Framework/Services/Registry/ServiceHandle.h" 
#include "art/Framework/Services/Optional/TFileService.h" 
#include "art/Framework/Services/Optional/TFileDirectory.h" 
#include "fhiclcpp/ParameterSet.h" 
#include "messagefacility/MessageLogger/MessageLogger.h" 
#include "cetlib/exception.h"
#include "cetlib/search_path.h"
#include "larcore/Geometry/Geometry.h"
#include "larevt/Filters/ChannelFilter.h"
#include "lardata/RawData/RawDigit.h"
#include "lardata/RawData/raw.h"
#include "lardata/RecoBase/Wire.h"
#include "lardata/Utilities/LArFFT.h"
#include "TComplex.h"
#include "TFile.h"
#include "TH2D.h"
#include "TF1.h"
#include "lardata/RecoBaseArt/WireCreator.h"
#include "lardata/Utilities/AssociationUtil.h"
// ################################
// ### LArIAT Specific Includes ###
// ################################
#include "RawDataUtilities/TriggerDigitUtility.h"
#include "Utilities/SignalShapingServiceT1034.h"


///creation of calibrated signals on wires
namespace caldata {
  
  class CalWireT1034 : public art::EDProducer {
  
  public:

    // create calibrated signals on wires. this class runs 
    // an fft to remove the electronics shaping.     
    explicit CalWireT1034(fhicl::ParameterSet const& pset); 
    virtual ~CalWireT1034();

    void produce(art::Event& evt); 
    void beginJob(); 
    void endJob();                 
    void reconfigure(fhicl::ParameterSet const& p);

  private:

    int          fDataSize;          ///< size of raw data on one wire
    int          fPostsample;        ///< number of postsample bins
    int          fBaseSampleBins;        ///< number of postsample bins
    float        fBaseVarCut;        ///< baseline variance cut
    std::string  fDigitModuleLabel;  ///< module that made digits
    ///< constants
    // === Size of FFT for ROI deconvolution ===          
    int fFFTSize; 
    std::string  fSpillName;  ///< nominal spill is an empty string
                              ///< it is set by the DigitModuleLabel
                              ///< ex.:  "daq:preSpill" for prespill data
    
    std::string fTriggerUtility;

    void SubtractBaseline(std::vector<float>& holder, int fBaseSampleBins);

  protected: 

  }; // class CalWireT1034
  DEFINE_ART_MODULE(CalWireT1034)

  

  //-------------------------------------------------
  CalWireT1034::CalWireT1034(fhicl::ParameterSet const& pset)
  {
    fSpillName="";
    this->reconfigure(pset);
    produces< std::vector<recob::Wire> >();
    produces<art::Assns<raw::RawDigit, recob::Wire>>();
    produces<art::Assns<raw::Trigger, recob::Wire>>();
    //produces<art::Assns<raw::Trigger, recob::Wire>>();
    
   
  }

  

  //-------------------------------------------------

  CalWireT1034::~CalWireT1034()
  {
  }

  //////////////////////////////////////////////////////

  void CalWireT1034::reconfigure(fhicl::ParameterSet const& p)

  {

    fDigitModuleLabel = p.get< std::string >("DigitModuleLabel", "daq");
    fPostsample       = p.get< int >        ("PostsampleBins");
    fBaseSampleBins   = p.get< int >        ("BaseSampleBins");
    fBaseVarCut       = p.get< int >        ("BaseVarCut");
    fTriggerUtility   = p.get< std::string >("TriggerUtility");
    fFFTSize          = p.get< int  >                 ("FFTSize");
    fSpillName="";

    size_t pos = fDigitModuleLabel.find(":");
    if( pos!=std::string::npos ) {
      fSpillName = fDigitModuleLabel.substr( pos+1 );
      fDigitModuleLabel = fDigitModuleLabel.substr( 0, pos );
    }

  
  }

  //-------------------------------------------------

  void CalWireT1034::beginJob()

  {  

  }

  //////////////////////////////////////////////////////

  void CalWireT1034::endJob()

  {  

  }

  

  //////////////////////////////////////////////////////
  void CalWireT1034::produce(art::Event& evt)
  {
  
      // ################################
      // ### Get the geometry service ###
      // ################################
    art::ServiceHandle<geo::Geometry> geom;
    
    // ###########################
    // ### Get the FFT service ###
    // ###########################
    art::ServiceHandle<util::LArFFT> fFFT;
    
    // === Recall the FFT size ===
    int transformSize = fFFT->FFTSize();
    //transformSize = fFFTSize;
    
    // ############################################
    // ### Get the Signal Shaping Service (sss) ###
    // ############################################
    art::ServiceHandle<util::SignalShapingServiceT1034> sss;
    
    // ####################################
    // ### Making a collection of Wires ###
    // ####################################
    std::unique_ptr<std::vector<recob::Wire> > wirecol(new std::vector<recob::Wire>);
    
    // ################################################
    // ### Making an association of Wires to Digits ###
    // ################################################
    std::unique_ptr<art::Assns<raw::RawDigit,recob::Wire> > WireDigitAssn(new art::Assns<raw::RawDigit,recob::Wire>);
    
    // ##################################################
    // ### Making an association of Wires to Triggers ###
    // ##################################################
    std::unique_ptr<art::Assns<raw::Trigger, recob::Wire> > TrigWireAssn(new art::Assns<raw::Trigger,recob::Wire>);
    
    // ################################################
    // ### Defining the size of the raw data vector ###
    // ################################################
    unsigned int dataSize = 0;
    
    // #########################################################
    // ### Defining the channel number and time bin variable ###
    // #########################################################
    uint32_t     channel(0); // channel number
    unsigned int bin(0);     // time bin loop variable
    
    // ### (Put comments here) ###
    filter::ChannelFilter *chanFilt = new filter::ChannelFilter();
    std::vector<float> holder;                // holds signal data
    std::vector<short> rawadc(transformSize);  // vector holding uncompressed adc values
    std::vector<TComplex> freqHolder(transformSize+1); // temporary frequency data
    
    
    // ###########################################
    // ### Grab the trigger data utility (tdu) ###
    // ###########################################
    rdu::TriggerDigitUtility tdu(evt, fTriggerUtility);
    
    
    // ##############################
    // ### Loop over the triggers ###
    // ##############################
    for(size_t trig = 0; trig < tdu.NTriggers(); trig++)
    {
      //std::cout<<"========================="<<std::endl;
      //std::cout<<"    Trigger Number = "<<trig<<std::endl;
      //std::cout<<"========================="<<std::endl;
    
      // === Getting the pointer for this trigger ===
      art::Ptr<raw::Trigger> trigger = tdu.EventTriggersPtr()[trig];
      art::PtrVector<raw::RawDigit> rdvec = tdu.TriggerRawDigitsPtr(trig);
      
      //std::cout<<"rdvec.size() = "<<rdvec.size()<<std::endl;
      
      // === Put in protection in case ther is no TPC raw digits ===
      if(!rdvec.size()){LOG_INFO("CalWireT1034") << "Problem CalWireT1034:: RawDigiVec size is " << rdvec.size(); continue;}
      
      // === Using the zeroth element to get the number of samples ===
      dataSize = rdvec[0]->Samples();
      
      //std::cout<<"Number of samples = "<<dataSize<<std::endl;
    
      // =======================================================================================
      // === Put in protection in case the transform size is less than the number of samples ===
      // =======================================================================================
      if ( (unsigned int)transformSize < dataSize /*&& trig < 1*/)
      {
        LOG_WARNING("CalWireT1034")<<"FFT size (" << transformSize << ") "
        << "is smaller than the data size (" << dataSize << ") "
        << "\nResizing the FFT now...";
        
        // === Resizing the FFT ====
        fFFT->ReinitializeFFT(dataSize,fFFT->FFTOptions(),fFFT->FFTFitBins());
        
        // === Resetting the transform size ===
        transformSize = fFFT->FFTSize();
        
        LOG_WARNING("CalWireT1034")<<"FFT size is now (" << transformSize << ") "
        << "and should be larger than the data size (" << dataSize << ")";
        
      }//<-end protection
      
      
      LOG_INFO("CalWireT1034") << "Data size is " << dataSize << " and transform size is " << transformSize;
      
      // === Protection in case sample size is not some multiple of the number of sample bins ===
      if(fBaseSampleBins > 0 && dataSize % fBaseSampleBins != 0) {LOG_ERROR("CalWireT1034")<<"Set BaseSampleBins modulo dataSize= "<<dataSize;}
      
      // === Reserve the right amount of data in the vector ===
      wirecol->reserve(rdvec.size());
      
      // ###########################
      // ### Loop over the wires ###
      // ###########################
      for(size_t rdItr = 0; rdItr < rdvec.size(); ++rdItr)
      {
        //std::cout<<"Looping over wire # "<<rdItr<<std::endl;
        
        holder.clear();
        
        // === Grab the current channel ===
        channel = rdvec[rdItr]->Channel();
        
        //std::cout<<"Channel = "<<channel<<std::endl;
      
        // === Only deconvolute non-"Bad" channels ====
        if(!chanFilt->BadChannel(channel))
        {
          // === Resize the holder vector and pad with zeros ====
          holder.resize(transformSize, 0.);
          
          // === Uncompress the data ===
          raw::Uncompress(rdvec[rdItr]->ADCs(), rawadc, rdvec[rdItr]->Compression());
          
          // === Getting the pedestal value ====
          // 06/23/15 brebel: as of this date we presubtract the pedestal
          // when making RawDigits to mesh with changes in LArSoft v04_13_00
          float pdstl = 0.;//rdvec[rdItr]->GetPedestal();
                           //std::cout<<"pedestal = "<<pdstl<<std::endl;
          
          // === Loop over all adc values and subtract the pedestal ===
          for(bin = 0; bin < dataSize; ++bin)
          {
            holder[bin]=(rawadc[bin]-pdstl);
          }
          
          // =================================
          // === Perform the Deconvolution ===
          // =================================
          
          /*std::cout<<std::endl;
           std::cout<<"Getting Ready to Deconvolute a wire"<<std::endl;
           std::cout<<std::endl;*/
          
          sss->Deconvolute(channel, holder);
          
          /*std::cout<<std::endl;
           std::cout<<std::endl;
           std::cout<<"Done Deconvolution of a wire"<<std::endl;
           std::cout<<std::endl;
           std::cout<<std::endl;*/
          
        } // end if not a bad channel
        
        // === Resize the holder vector ====
        holder.resize(dataSize,1e-5);
        
        //This restores the DC component to signal removed by the deconvolution.
        if(fPostsample)
        {
          float average=0.0;
          for(bin=0; bin < (unsigned short)fPostsample; ++bin){average += holder[holder.size()-1+bin];}
          
          average = average / (float)fPostsample;
          
          
          for(bin = 0; bin < holder.size(); ++bin) {holder[bin]-=average;}
        }
        
        // adaptive baseline subtraction
        if(fBaseSampleBins) SubtractBaseline(holder, fBaseSampleBins);
        // Make a single ROI that spans the entire data size
        wirecol->push_back(recob::WireCreator(holder,*rdvec[rdItr]).move());
        
        // add an association between the last object in wirecol
        // (that we just inserted) and digitVec
        // util::CreateAssn(*this, evt, *wirecol, rdvec[rdItr], *WireDigitAssn);
        if (!util::CreateAssn(*this, evt, *wirecol, rdvec[rdItr], *WireDigitAssn))
        {
          throw art::Exception(art::errors::InsertFailure) << "Can't associate wire #" << (wirecol->size() - 1)
          << " with raw digit #" << rdvec[rdItr].key();
        } // if failed to add association
        //util::CreateAssn(*this, evt, *wirecol, trigger, *TrigWireAssn);
        // add an association between the last object in wirecol
        // (that we just inserted) and digitVec
        if (!util::CreateAssn(*this, evt, *wirecol, trigger, *TrigWireAssn))
        {
          throw art::Exception(art::errors::InsertFailure) << "Can't associate wire #" << (wirecol->size() - 1)
          << " with Trigger #" << trigger.key();
        } // if failed to add association      
        
        
      }//<---End looping over wires
      
      // ###########################################################
      // ### Send a warning if no wires were made for this event ###
      // ###########################################################
      if(wirecol->size() == 0){LOG_WARNING("CalWireT1034") << "No wires made for this event.";}
      
      
    }//<---End trig loop
    
    evt.put(std::move(wirecol));
    evt.put(std::move(WireDigitAssn));
    evt.put(std::move(TrigWireAssn));
    delete chanFilt;
    return;
    
  }//<---End Event loop





  void CalWireT1034::SubtractBaseline(std::vector<float>& holder,
				      int fBaseSampleBins)
  {
    // subtract baseline using linear interpolation between regions defined
    // by the datasize and fBaseSampleBins
    // number of points to characterize the baseline
    unsigned short nBasePts = holder.size() / fBaseSampleBins;
    // the baseline offset vector
    std::vector<float> base;
    for(unsigned short ii = 0; ii < nBasePts; ++ii) base.push_back(0.);
    // find the average value in each region, using values that are
    // similar
    float fbins = fBaseSampleBins;
    unsigned short nfilld = 0;
    for(unsigned short ii = 0; ii < nBasePts; ++ii) {

      unsigned short loBin = ii * fBaseSampleBins;

      unsigned short hiBin = loBin + fBaseSampleBins;

      float ave = 0.;

      float sum = 0.;

      for(unsigned short bin = loBin; bin < hiBin; ++bin) {

        ave += holder[bin];

        sum += holder[bin] * holder[bin];

      } // jj

      ave = ave / fbins;

      float var = (sum - fbins * ave * ave) / (fbins - 1.);

      // Set the baseline for this region if the variance is small

      if(var < fBaseVarCut) {

        base[ii] = ave;

        ++nfilld;

      }

    } // ii

    // fill in any missing points if there aren't too many missing

    if(nfilld < nBasePts && nfilld > nBasePts / 2) {

      bool baseOK = true;

      // check the first region

      if(base[0] == 0) {

        unsigned short ii1 = 0;

        for(unsigned short ii = 1; ii < nBasePts; ++ii) {

          if(base[ii] != 0) {

            ii1 = ii;

            break;

          }

        } // ii

        unsigned short ii2 = 0;

        for(unsigned short ii = ii1 + 1; ii < nBasePts; ++ii) {

          if(base[ii] != 0) {

            ii2 = ii;

            break;

          }

        } // ii

        // failure

        if(ii2 > 0) {

          float slp = (base[ii2] - base[ii1]) / (float)(ii2 - ii1);

          base[0] = base[ii1] - slp * ii1;

        } else {

          baseOK = false;

        }

      } // base[0] == 0

      // check the last region

      if(baseOK && base[nBasePts] == 0) {

        unsigned short ii2 = 0;

        for(unsigned short ii = nBasePts - 1; ii > 0; --ii) {

          if(base[ii] != 0) {

            ii2 = ii;

            break;

          }

        } // ii

        baseOK = false; // assume the worst, hope for better

        unsigned short ii1 = 0;

        if (ii2 >= 1) {

          for(unsigned short ii = ii2 - 1; ii > 0; --ii) {

            if(base[ii] != 0) {

              ii1 = ii;

              baseOK = true;

              break;

            } // if base[ii]

          } // ii

        } // if ii2

        if (baseOK) {

          float slp = (base[ii2] - base[ii1]) / (float)(ii2 - ii1);

          base[nBasePts] = base[ii2] + slp * (nBasePts - ii2);

        }

      } // baseOK && base[nBasePts] == 0

      // now fill in any intermediate points

      for(unsigned short ii = 1; ii < nBasePts - 1; ++ii) {

        if(base[ii] == 0) {

          // find the next non-zero region

          for(unsigned short jj = ii + 1; jj < nBasePts; ++jj) {

            if(base[jj] != 0) {

              float slp = (base[jj] - base[ii - 1]) / (jj - ii + 1);

              base[ii] = base[ii - 1] + slp;

              break;

            }

          } // jj

        } // base[ii] == 0

      } // ii

    } // nfilld < nBasePts

    // interpolate and subtract

    float slp = (base[1] - base[0]) / (float)fBaseSampleBins;

    // bin offset to the origin (the center of the region)

    unsigned short bof = fBaseSampleBins / 2;

    unsigned short lastRegion = 0;

    for(unsigned short bin = 0; bin < holder.size(); ++bin) {

      // in a new region?

      unsigned short region = bin / fBaseSampleBins;

      if(region > lastRegion) {

        // update the slope and offset

        slp = (base[region] - base[lastRegion]) / (float)fBaseSampleBins;

        bof += fBaseSampleBins;

        lastRegion = region;

      }

      holder[bin] -= base[region] + (bin - bof) * slp;

    }

  } // SubtractBaseline

  

} // end namespace caldata

