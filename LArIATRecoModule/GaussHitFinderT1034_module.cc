#ifndef GAUSHITFINDER_H
#define GAUSHITFINDER_H

////////////////////////////////////////////////////////////////////////
//
// GaussHitFinder class
//
// jaasaadi@syr.edu
//
//  This algorithm is designed to find hits on wires after deconvolution.
// -----------------------------------
// This algorithm is based on the FFTHitFinder written by Brian Page, 
// Michigan State University, for the ArgoNeuT experiment.
// 
//
// The algorithm walks along the wire and looks for pulses above threshold
// The algorithm then attempts to fit n-gaussians to these pulses where n
// is set by the number of peaks found in the pulse
// If the Chi2/NDF returned is "bad" it attempts to fit n+1 gaussians to
// the pulse. If this is a better fit it then uses the parameters of the 
// Gaussian fit to characterize the "hit" object 
//
// To use this simply include the following in your producers:
// gaushit:     @local::microboone_gaushitfinder
// gaushit:	@local::argoneut_gaushitfinder
////////////////////////////////////////////////////////////////////////


// C/C++ standard library
#include <algorithm> // std::accumulate()
#include <vector>
#include <string>
#include <utility> // std::move()


// Framework includes
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/EDProducer.h"
#include "canvas/Persistency/Common/FindOneP.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Optional/TFileService.h"



// LArSoft Includes
#include "larcoreobj/SimpleTypesAndConstants/RawTypes.h" // raw::ChannelID_t
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/CryostatGeo.h"
#include "larcorealg/Geometry/TPCGeo.h"
#include "larcorealg/Geometry/PlaneGeo.h"
#include "lardataobj/RecoBase/Wire.h"
#include "lardataobj/RecoBase/Hit.h"
//#include "RecoBaseArt/HitCreator.h"
//#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/Utilities/AssociationUtil.h"

// ROOT Includes
#include "TGraphErrors.h"
#include "TH1D.h"
#include "TDecompSVD.h"
#include "TMath.h"
#include "TF1.h"
#include "TTree.h"
#include "TStopwatch.h"

// ################################
// ### LArIAT Specific Includes ###
// ################################
#include "RawDataUtilities/TriggerDigitUtility.h"

namespace hit{
  class GausHitFinder : public art::EDProducer {
    
  public:
    
    explicit GausHitFinder(fhicl::ParameterSet const& pset); 
    virtual ~GausHitFinder();
         
    void produce(art::Event& evt); 
    void beginJob(); 
    void endJob(); 
    void reconfigure(fhicl::ParameterSet const& p);                

  private:
  
    // ### This function will fit N-Gaussians to at TH1D where N is set ###
    // ###            by the number of peaks found in the pulse         ###
    void FitGaussians(std::vector<float> SignalVector, std::vector<int> PeakTime, int nGauss, int StartTime, int EndTime);

    // ### This function will fit N+1-Gaussians to at TH1D where N+1 is set ###
    // ###     by the number of peaks found plus one more in the pulse      ###
    void ReFitGaussians(std::vector<float> ReSignalVector,std::vector<int> RePeakTime, int mGauss, int ReStartTime, int ReEndTime);
    
    void FillOutHitParameterVector(const std::vector<double>& input,
				   std::vector<double>& output);

    // load the fit into a temp vector
    std::vector<double> tempGaus1Par;
    std::vector<double> tempGaus1ParError;
    double Chi2PerNDF = 0.;
    int NDF = -1;
    
    // load the fit into a temp vector
    std::vector<double> tempGaus2Par;
    std::vector<double> tempGaus2ParError;
    double Chi2PerNDF2 = 0.;
    int NDF2 = -1;
    
    
    double threshold              = 0.;  // minimum signal size for id'ing a hit
    double fitWidth               = 0.;  // hit fit width initial value
    double minWidth		  = 0 ;  // hit minimum width
    std::string     fCalDataModuleLabel;

    std::vector<double> fMinSig;    ///<signal height threshold
    std::vector<double> fInitWidth; ///<Initial width for fit
    std::vector<double> fMinWidth;  ///<Minimum hit width
    
    int             fMaxMultiHit;   ///<maximum hits for multi fit 
    int             fAreaMethod;    ///<Type of area calculation  
    std::vector<double> fAreaNorms; ///<factors for converting area to same units as peak height 
    int		    fTryNplus1Fits; ///<whether we will (0) or won't (1) try n+1 fits
    double	    fChi2NDFRetry;  ///<Value at which a second n+1 Fit will be tried
    double	    fChi2NDF;       ///maximum Chisquared / NDF allowed for a hit to be saved
    
    TH1F* fFirstChi2;
    TH1F* fChi2;
    
    std::string fTriggerUtility;   ///Module name which made the TriggerDigit Block
		
  protected: 
    
  
  }; // class GausHitFinder
  

//-------------------------------------------------
//-------------------------------------------------
GausHitFinder::GausHitFinder(fhicl::ParameterSet const& pset)
{
  this->reconfigure(pset);
  
  // let HitCollectionCreator declare that we are going to produce
  // hits and associations with wires and raw digits
  // (with no particular product label)
  //recob::HitCollectionCreator::declare_products(*this);
  produces< std::vector<recob::Hit> >();
  
  // ### Produce an association between hits and Triggers ###
  produces<art::Assns<raw::Trigger, recob::Hit>>();
  
} // GausHitFinder::GausHitFinder()


//-------------------------------------------------
//-------------------------------------------------
  GausHitFinder::~GausHitFinder()
{
  // Clean up dynamic memory and other resources here.

}
  
//-------------------------------------------------
//-------------------------------------------------
void GausHitFinder::FillOutHitParameterVector(const std::vector<double>& input,
					      std::vector<double>& output)
{
  if(input.size()==0)
    throw std::runtime_error("GausHitFinder::FillOutHitParameterVector ERROR! Input config vector has zero size.");

   art::ServiceHandle<geo::Geometry> geom;
   const unsigned int N_PLANES = geom->Nplanes();

  if(input.size()==1)
    output.resize(N_PLANES,input[0]);
  else if(input.size()==N_PLANES)
    output = input;
  else
    throw std::runtime_error("GausHitFinder::FillOutHitParameterVector ERROR! Input config vector size !=1 and !=N_PLANES.");
      
}
						
  
//-------------------------------------------------
//-------------------------------------------------
void GausHitFinder::reconfigure(fhicl::ParameterSet const& p)
{
  // Implementation of optional member function here.
  fCalDataModuleLabel = p.get< std::string  >("CalDataModuleLabel");
  fTriggerUtility = p.get< std::string >("TriggerUtility");
  FillOutHitParameterVector(p.get< std::vector<double> >("MinSig"),fMinSig);
  FillOutHitParameterVector(p.get< std::vector<double> >("InitWidth"),fInitWidth);
  FillOutHitParameterVector(p.get< std::vector<double> >("MinWidth"),fMinWidth);
  FillOutHitParameterVector(p.get< std::vector<double> >("AreaNorms"),fAreaNorms);

  
  fMaxMultiHit        = p.get< int          >("MaxMultiHit");
  fAreaMethod         = p.get< int          >("AreaMethod");
  fTryNplus1Fits      = p.get< int          >("TryNplus1Fits");
  fChi2NDFRetry       = p.get< double       >("Chi2NDFRetry");
  fChi2NDF            = p.get< double       >("Chi2NDF");
  
}  

//-------------------------------------------------
//-------------------------------------------------
void GausHitFinder::beginJob()
{
// get access to the TFile service
art::ServiceHandle<art::TFileService> tfs;
   
    
// ======================================
// === Hit Information for Histograms ===
fFirstChi2	= tfs->make<TH1F>("fFirstChi2", "#chi^{2}", 10000, 0, 5000);
fChi2	        = tfs->make<TH1F>("fChi2", "#chi^{2}", 10000, 0, 5000);


}

//-------------------------------------------------
//-------------------------------------------------
void GausHitFinder::endJob()
{

}

//  This algorithm uses the fact that deconvolved signals are very smooth 
//  and looks for hits as areas between local minima that have signal above 
//  threshold.
//-------------------------------------------------
void GausHitFinder::produce(art::Event& evt)
{
//==================================================================================================  
   
   TH1::AddDirectory(kFALSE);
   
   // Instantiate and Reset a stop watch
    //TStopwatch StopWatch;
    //StopWatch.Reset();
     
   // #########################################################
   // ### List of useful variables used throughout the code ###
   // #########################################################
   
   //double minWidth               = 0.;  //minimum hit width
   std::vector<int> startTimes;         // stores time of 1st local minimum
   std::vector<int> maxTimes;    	// stores time of local maximum    
   std::vector<int> endTimes;    	// stores time of 2nd local minimum
   
   std::vector <int> peakT;
   
   std::vector<double> StartTime;      // stores the start time of the hit
   std::vector<double> StartTimeError; // stores the error assoc. with the start time of the hit
   std::vector<double> EndTime;	       // stores the end time of the hit
   std::vector<double> EndTimeError;   // stores the error assoc. with the end time of the hit
   std::vector<double> RMS;            // stores the sigma parameter of the gaussian
   std::vector<double> MeanPosition;   // stores the peak time position of the hit
   std::vector<double> MeanPosError;   // stores the error assoc. with thte peak time of the hit
   std::vector<double> SumADC;         // stores the total charge from sum of ADC counts
   std::vector<double> Charge;         // stores the total charge assoc. with the hit (from the fit)
   std::vector<double> ChargeError;    // stores the error on the charge
   std::vector<double> Amp;	       // stores the amplitude of the hit
   std::vector<double> AmpError;       // stores the error assoc. with the amplitude
   std::vector<double> NumOfHits;      // stores the multiplicity of the hit
   std::vector<double> FitGoodness;    // stores the Chi2/NDF of the hit
   std::vector<int>    FitNDF;         // stores the NDF of the hit
   std::vector<double>  hitSig;
   
   // ###################################
   // ### Calling Detector Properties ###
   // ###################################
   //auto const* detprop = lar::providerFrom<detinfo::DetectorPropertiesService>();
   
   // ################################
   // ### Calling Geometry service ###
   // ################################
   art::ServiceHandle<geo::Geometry> geom;

   // ###############################################
   // ### Making a ptr vector to put on the event ###
   // ###############################################
   std::unique_ptr<std::vector<recob::Hit> > hcol(new std::vector<recob::Hit>);
   
   // ##########################################
   // ### Reading in the Wire List object(s) ###
   // ##########################################   
   std::vector< art::Ptr<recob::Wire> >wireVecHandle;
   
   // Signal Type (Collection or Induction)
   geo::SigType_t sigType;
     
   // #################################################
   // ### Making an association of Triggers to Hits ###
   // #################################################
   std::unique_ptr<art::Assns<raw::Trigger, recob::Hit> > TrigHitAssn(new art::Assns<raw::Trigger,recob::Hit>);
   
   // Channel Number
   raw::ChannelID_t channel = raw::InvalidChannelID;
   
   // ###########################################
   // ### Grab the trigger data utility (tdu) ###
   // ###########################################
   rdu::TriggerDigitUtility tdu(evt, fTriggerUtility);
   art::PtrVector<raw::Trigger> const& EventTriggersPtr = tdu.EventTriggersPtr();

   
   // #################################################################
   // ### Reading in the RawDigit associated with these wires, too  ###
   // #################################################################
   //art::FindManyP<raw::RawDigit> RawDigits
     //(wireVecHandle, evt, fCalDataModuleLabel);
   
   //std::cout<<"RawDigits(tdu.EventTriggersPtr(), evt, fTriggerUtility);"<<std::endl;  
   art::FindManyP<raw::RawDigit> RawDigits(tdu.EventTriggersPtr(), evt, fTriggerUtility);
   art::FindManyP<recob::Wire>   CalWireDigits(tdu.EventTriggersPtr(), evt, fCalDataModuleLabel);
   
   // ##############################
   // ### Loop over the triggers ###
   // ##############################
   size_t startHit = 0;
   for(size_t trig = 0; trig < tdu.NTriggers(); trig++)
      {

      // get the starting index of the hits for this trigger
      startHit = hcol->size();
      
      // === Getting the pointer for this trigger ===
      //art::Ptr<raw::Trigger> trigger = tdu.EventTriggersPtr()[trig];
     art::Ptr<raw::Trigger> theTrigger = (EventTriggersPtr[trig]);

     art::PtrVector<raw::RawDigit> rdvec = tdu.TriggerRawDigitsPtr(trig);

      if(!rdvec.size()){mf::LogInfo("GausHitFinder") << "Problem GaussHitFinder:: RawDigiVec size is " << rdvec.size(); continue;}
      wireVecHandle.clear();
      wireVecHandle = CalWireDigits.at(trig);
      std::cout<<"wireVecHandle.size()"<<wireVecHandle.size()<<std::endl;
      //##############################
      //### Looping over the wires ###
      //############################## 
      for(size_t wireIter = 0; wireIter < wireVecHandle.size(); ++wireIter)
         {
         // -----------------------------------------------------------
         // -- Clearing variables at the start of looping over wires --
         // -----------------------------------------------------------
         startTimes.clear();
         maxTimes.clear();
         endTimes.clear();
      
         StartTime.clear();
         StartTimeError.clear();
         EndTime.clear();
         EndTimeError.clear();
         RMS.clear();
         MeanPosition.clear();
         MeanPosError.clear();
         Amp.clear();
         AmpError.clear();
         NumOfHits.clear();
         hitSig.clear();
         SumADC.clear();
         Charge.clear();
         ChargeError.clear();
         tempGaus1Par.clear();
      
      
      
         // ####################################
         // ### Getting this particular wire ###
         // ####################################
         //art::Ptr<recob::Wire> wire(wireVecHandle[wireIter], wireIter);
	 art::Ptr<recob::Wire> wire = wireVecHandle[wireIter];
	 
         std::vector< art::Ptr<raw::RawDigit> > rawdigits = RawDigits.at(trig);
      
         // --- Setting Channel Number and Signal type ---
         channel = wireVecHandle[wireIter]->Channel();
	 sigType = geom->SignalType(channel);
         // ----------------------------------------------------------
         // -- Setting the appropriate signal widths and thresholds --
         // --    for the right plane.      --
         // ----------------------------------------------------------

         threshold     = fMinSig.at(wireVecHandle[wireIter]->View());
         fitWidth      = fInitWidth.at(wireVecHandle[wireIter]->View());
         minWidth      = fMinWidth.at(wireVecHandle[wireIter]->View());	
	
         // #################################################
         // ### Getting a vector of signals for this wire ###
         // #################################################
         std::vector<float> signal(wireVecHandle[wireIter]->Signal());
      
         // ##########################################################
         // ### Making an iterator for the time ticks of this wire ###
         // ##########################################################
         std::vector<float>::iterator timeIter;  	    // iterator for time bins
      
         // current time bin
         int time             = 0;
      
         // current start time
         int minTimeHolder    = 0; 
       
         // Flag for whether a peak > threshold has been found
         bool maxFound        = false;                    
         //std::cout<<std::endl;
         //std::cout<<" Wire # "<<wireIter<<std::endl;
         //std::cout<<"===================="<<std::endl; 
      
         // ### Starting the StopWatch ###
         //StopWatch.Start(); 
          
         // ##################################
         // ### Looping over Signal Vector ###
         // ##################################
         for(timeIter = signal.begin(); timeIter+2<signal.end(); timeIter++)
      	    {
	 
	    // ##########################################################
      	    // ###                LOOK FOR A MINIMUM                  ###
      	    // ### Testing if the point timeIter+1 is at a minimum by ###
      	    // ###  checking if timeIter and timeIter+2 are greater   ###
      	    // ###   and if it is then we add this to the endTimes    ###
      	    // ##########################################################
	    if(*(timeIter+1) < *timeIter && *(timeIter+1) < *(timeIter+2))
	       {
	       //--- Note: We only keep the a minimum if we've already ---
	       //---          found a point above threshold            ---
	       if(maxFound)
	          {
	          endTimes.push_back(time+1);
	          maxFound = false;
	          //keep these in case new hit starts right away
	          minTimeHolder = time+2; 
	          }
	       else 
	          {minTimeHolder = time+1;} 
	  
      	       }//<---End Checking if this is a minimum
	 
	 
            // ########################################################	
            // ### Testing if the point timeIter+1 is a maximum and ###
            // ###  if it and is above threshold then we add it to  ###
            // ###                  the startTime                   ###
            // ########################################################
            //if not a minimum, test if we are at a local maximum 
            //if so, and the max value is above threshold, add it and proceed.
	    else if(*(timeIter+1) > *timeIter && *(timeIter+1) > *(timeIter+2) && *(timeIter+1) > threshold)
	       { 
	       maxFound = true;
	       maxTimes.push_back(time+1);
	       startTimes.push_back(minTimeHolder);          
               }
	
	
            time++;
            }//<---End timeIter loop
       
         // ###########################################################    
         // ### If there was no minimum found before the end, but a ### 
         // ###  maximum was found then add an end point to the end ###
         // ###########################################################
         while( maxTimes.size() > endTimes.size() )
            { endTimes.push_back(signal.size()-1); }//<---End maxTimes.size > endTimes.size
	
         // ####################################################
         // ### If no startTime hit was found skip this wire ###
         // ####################################################
         if( startTimes.size() == 0 ){continue;}

      
      
         // #####################################
         // ### Variables used for each pulse ###
         // #####################################
         int startT = 0;
         int endT = 0;
         double totSig = 0;
         peakT.clear();

         // #######################################################
         // ### Lets loop over the pulses we found on this wire ###
         // #######################################################
         for(size_t num = 0; num < maxTimes.size(); num++)
            {
	    peakT.clear();
	 
	    // ### Setting the start, peak, and end time of the pulse ###
	    startT = startTimes[num];
	    endT   = endTimes[num];
	    peakT.push_back(maxTimes[num]);
	 
	    // ###       See if we want to merge pulses together          ###
	    // ### First check if we have more than one pulse on the wire ###
 	    if(endTimes.size() > 1 )
	       {

	       // ###      If the start time of the next pulse is one wire away and    ###
	       // ### the height of that start point is 0.5* the thereshold then merge ###
	       if(num < (maxTimes.size() - 1) && startTimes[num+1] - endTimes[num]  < 2 && signal[endTimes[num]] > threshold/2)
	          {
	          endT = endTimes[num+1];
	          peakT.push_back(maxTimes[num+1]);
	       
	          num++;
	          }//<---Checking adjacent pulses

	       }//<---End checking if there is more than one pulse on the wire

	    // ### Putting in a protection in case things went wrong ###
	    if(startT - endT == 0){startT = 0; endT = 9600;}
	    if(startT > endT){endT = 9999;}
	 
	    int nFill = 0;
	    // #######################################################
	    // ### Clearing the parameter vector for the new pulse ###
	    // #######################################################
	    tempGaus1Par.clear();
	 
	    // === Setting the number of Gaussians to try ===
	    int nGausForFit = peakT.size();
	    int nGausReFit = nGausForFit + 1;
	 
	    nFill = nGausForFit;
	    // ##################################################
	    // ### Calling the function for fitting Gaussians ###
	    // ##################################################
	    FitGaussians(signal, peakT, nGausForFit, startT, endT);
	    fFirstChi2->Fill(Chi2PerNDF);
	 
	    // #######################################################
	    // ### Clearing the parameter vector for the new pulse ###
	    // #######################################################
	    tempGaus2Par.clear();
            bool secondTry = false;
	    // #####################################################
	    // ### Trying extra gaussians for an initial bad fit ###
	    // #####################################################
	    if( (Chi2PerNDF > (2*fChi2NDFRetry) && fTryNplus1Fits == 0 && nGausForFit == 1)|| 
	        (Chi2PerNDF > (fChi2NDFRetry) && fTryNplus1Fits == 0 && nGausForFit > 1))
	       {
	    
	       nFill = nGausReFit;
	       // #########################################################
	       // ### Calling the function for re-fitting n+1 Gaussians ###
	       // #########################################################
	       ReFitGaussians(signal, peakT, nGausReFit, startT, endT);
	    
	       secondTry = true;
	       }
	    
	    
	    // ### Creating temperary parameters for storing hits ###
	    double tempChi2NDF = 0;
	    int tempNDF = -1;
	    std::vector<double> tempPar;
            std::vector<double> tempParError;
	 
	    // #########################################################
	    // ### Getting the appropriate parameter into the vector ###
	    // #########################################################
	    if(secondTry && Chi2PerNDF2 < Chi2PerNDF)
	       {
	       nFill = nGausReFit;
	       tempChi2NDF = Chi2PerNDF2;
	       tempNDF = NDF2;
	       // ### Filling the vector of information
	       for(size_t ab = 0; ab < tempGaus2Par.size(); ab++)
	          {
	          tempPar.push_back( tempGaus2Par[ab] );
                  tempParError.push_back( tempGaus2ParError[ab] );  
		  
		  
	          }//<---End ab loop
	       }//<---End using n+1 Gaussian Fit
	       
	    else
	       {
	       nFill = nGausForFit;
	       tempChi2NDF = Chi2PerNDF;
	       tempNDF = NDF;
	       // ### Filling the vector of information
	       for(size_t ab = 0; ab < tempGaus1Par.size(); ab++)
	          {
	          tempPar.push_back( tempGaus1Par[ab] );
                  tempParError.push_back( tempGaus1ParError[ab] );  
		  
		  
	          }//<---End ab loop

	       }//<---End using n Gaussian Fit
	    fChi2->Fill(tempChi2NDF);
	 
	    // #################################################
	    // ### Clearing the variables for this new pulse ###
	    // #################################################  
	    totSig = 0;
	    StartTime.clear();
	    StartTimeError.clear();
	    EndTime.clear();
	    EndTimeError.clear();
	    RMS.clear();
	    MeanPosition.clear();
	    MeanPosError.clear();
	    Amp.clear();
	    AmpError.clear();
	    NumOfHits.clear();
	    hitSig.clear();
	    Charge.clear();
	    SumADC.clear();
	    ChargeError.clear();
	    FitGoodness.clear();
	    FitNDF.clear();
	 
	    int numHits = 0;
	    // #################################################
	    // ### Recording all the relevant hit parameters ###
	    // #################################################
	 
	    for(int cc = 0; cc < nFill; cc++)
	       {
	       // #############################
	       // ### Skip poorly made hits ###
	       // #############################
	       // Peak Below threshold
	       if(tempPar[3*cc] < threshold){continue;}
	       // Poor Chi2
	       if(tempChi2NDF > fChi2NDF){continue;} 
	    
	       // ### Start Time ###
	       StartTime.push_back( tempPar[(3*cc)+1] - tempPar[(3*cc)+2] ); //<---(Mean - Width)
	       StartTimeError.push_back(TMath::Sqrt( (tempParError[(3*cc)+1]*tempParError[(3*cc)+1]) + 
		   			             (tempParError[(3*cc)+2]*tempParError[(3*cc)+2])) );
	    
	       // ### End Time ###	
	       EndTime.push_back( tempPar[(3*cc)+1] + tempPar[(3*cc)+2] ); //<---(Mean + Width)
	       EndTimeError.push_back( TMath::Sqrt( (tempParError[(3*cc)+1]*tempParError[(3*cc)+1]) + 
					            (tempParError[(3*cc)+2]*tempParError[(3*cc)+2])) );
						 
	       // ### Mean Time ###
	       MeanPosition.push_back( tempPar[(3*cc)+1] );
	       MeanPosError.push_back( tempParError[(3*cc)+1] );
	    
	       // ### Width ###
	       RMS.push_back( tempPar[(3*cc)+2] );
	    
	       // ### Amplitude ###
	       Amp.push_back( tempPar[3*cc] );
	       AmpError.push_back( tempParError[3*cc] );
	    
	       // ### Number of hits in this pulse ###
	       int mul = nFill;
	       NumOfHits.push_back( mul );
	    
	       // ### Chi^2 / NDF ###
	       FitGoodness.push_back( tempChi2NDF );
	       FitNDF.push_back( tempNDF );
	    
	       // ### Charge ###
	       hitSig.resize(endT - startT);
	    
	       // ##################################
	       // ### Integral Method for charge ###
	       // ##################################
	       for(int sigPos = 0; sigPos<(endT - startT); sigPos++)//<---Loop over the size (endT - startT)
	          { 
	          hitSig[sigPos] = tempPar[3*cc]*TMath::Gaus(sigPos+startT,tempPar[(3*cc)+1], tempPar[(3*cc)+2]);
	          totSig+=hitSig[(int)sigPos];
	          }//<---End Signal postion loop
	     
	       // ###################################################### 
	       // ### Getting the total charge using the area method ###
	       // ######################################################
	       if(fAreaMethod) 
	          {
		    totSig = std::sqrt(2*TMath::Pi())*tempPar[3*cc]*tempPar[(3*cc)+2]/fAreaNorms[(size_t)(wireVecHandle[wireIter]->View())];
	          }//<---End Area Method
	    
	       Charge.push_back( totSig );
	       ChargeError.push_back( TMath::Sqrt(TMath::Pi())*(tempParError[(3*cc)+0]*tempParError[(3*cc)+2]+
	  			      tempParError[(3*cc)+2]*tempParError[(3*cc)+0]) );   //estimate from area of Gaussian
	    
	    
	       // ### Sum of ADC counts ### FIXME I don't think this copes well with overlaps
	       SumADC.push_back(
	         std::accumulate(signal.begin() + startT, signal.begin() + endT, 0.)
	         );
	    
	       numHits++;
	       }//<---End cc loop
	          
            tempPar.clear();
	    tempParError.clear();
	 
	    // get the WireID for this hit
	    std::vector<geo::WireID> wids = geom->ChannelToWire(channel);
	    // for now, just take the first option returned from ChannelToWire
	    geo::WireID wid = wids[0];
	    geo::View_t view = wire->View();
	    
	    raw::TDCtick_t ST = startT;
	    
	    raw::TDCtick_t ET = endT;
	    // ### Recording each hit in the pulse ###
          	
	    for(size_t dd = 0; dd < MeanPosition.size(); dd++)
	       { 
		 recob::Hit hit(
		   channel, 		//raw::ChannelID_t        channel
		   ST,	    		//raw::TDCtick_t          start_tick
		   ET,			//raw::TDCtick_t          end_tick
		   MeanPosition[dd],	//float                   peak_time
		   MeanPosError[dd],	//float                   sigma_peak_time
		   RMS[dd],
		   Amp[dd],
		   AmpError[dd],
		   SumADC[dd],
		   Charge[dd],
		   ChargeError[dd],
		   NumOfHits[dd],
		   dd,
		   FitGoodness[dd],
		   FitNDF[dd] ,
		   view,
		   sigType,
		   wid
		   );
		 
		 //hcol.emplace_back(hit.move(), wire, rawdigits[wireIter]);
		 hcol->push_back(hit);
	       }
	    
	    }//<---End num loop
	 

         }//<---End looping over all the wires

      // make sure we don't have the wires from the previous trigger
      wireVecHandle.clear();


      // ######################################################
      // ### Creating association between hits and triggers ###
      // ######################################################
      for(size_t h = startHit; h < hcol->size(); ++h)
         {
	 if(!util::CreateAssn(*this, evt, *hcol, theTrigger, *TrigHitAssn, h))
	 {throw art::Exception(art::errors::ProductRegistrationFailure) <<"Failed to associate hit "<< h << " with trigger "<<theTrigger.key();} // exception

         }//<---End h loop
	 
      if(hcol->size() == 0){mf::LogWarning("GaussHitFinder") << "No hits made for this trigger.";}
      // move the hit collection and the associations into the event
      //hcol.put_into(evt);
      
      }//<---End loop over trigger
   

   evt.put(std::move(hcol));
   evt.put(std::move(TrigHitAssn));
   return;
//==================================================================================================  
// End of the event  
   
  
} // End of produce() 

// --------------------------------------------------------------------------------------------
// Fit Gaussians
// --------------------------------------------------------------------------------------------
void hit::GausHitFinder::FitGaussians(std::vector<float> SignalVector, std::vector<int> PeakTime, int nGauss, int StartTime, int EndTime)
   {
   std::string eqn = "gaus(0)";  // string for equation for gaus fit
   std::stringstream numConv;
   double amplitude(0); //fit parameters
   
   tempGaus1Par.clear();
   tempGaus1ParError.clear();
   
   
   int size = EndTime - StartTime;
   // #############################################
   // ### If size < 0 then set the size to zero ###
   // #############################################
   if(EndTime - StartTime < 0){size = 0;}
   //std::cout<<"size = "<<size<<std::endl;
   // --- TH1D HitSignal ---
   TH1F hitSignal("hitSignal","",std::max(size,1),StartTime,EndTime);
   hitSignal.Sumw2();
   
   // #############################
   // ### Filling the histogram ###
   // #############################
   for(int aa = StartTime; aa < EndTime; aa++)
      {
      //std::cout<<"Time Tick = "<<aa<<", ADC = "<<signal[aa]<<std::endl;
      hitSignal.Fill(aa,SignalVector[aa]);

      
      if(EndTime > 10000){break;} // FIXME why?
      }//<---End aa loop

   // ############################################
   // ### Building TFormula for basic Gaussian ###
   // ############################################
   eqn = "gaus(0)";
   for(int i = 3; i < (nGauss)*3; i+=3)
      {
      eqn.append("+gaus(");
      numConv.str("");
      numConv << i;
      eqn.append(numConv.str());
      eqn.append(")");
      }   
   //std::cout<<"Fit 1 Check 4"<<std::endl;
   // ---------------------------------	 
   // --- TF1 function for GausHit  ---
   // ---------------------------------
   TF1 Gaus("Gaus",eqn.c_str(),0,std::max(size,1));
   
   // ### Setting the parameters for the Gaussian Fit ###
   for(int bb = 0; bb < nGauss; bb++)
      {
      amplitude = SignalVector[PeakTime[bb]];
      //std::cout<<"fitWidth = "<<fitWidth<<", amplitude = "<<SignalVector[PeakTime[bb]]<<std::endl;
      Gaus.SetParameter(3*bb,amplitude);
      Gaus.SetParameter(1+(3*bb), PeakTime[bb]);
      Gaus.SetParameter(2+(3*bb), fitWidth);
      Gaus.SetParLimits(3*bb, 0.0, 3.0*amplitude);
      Gaus.SetParLimits(1+(3*bb), StartTime , EndTime);
      Gaus.SetParLimits(2+(3*bb), 0.0, 15.0*fitWidth);
      }//<---End bb loop
 
   // ####################################################
   // ### PERFORMING THE TOTAL GAUSSIAN FIT OF THE HIT ###
   // ####################################################
   
   //hitSignal.Fit(&Gaus,"QNRWIB","", StartTime, EndTime);
   //hitSignal.Fit(&Gaus,"QNRWB","", StartTime, EndTime);
   //hitGraph.Fit(&Gaus,"QNB","",StartTime, EndTime);
   
   //std::cout<<"StartTime = "<<StartTime<<" , EndTime = "<<EndTime<<std::endl;
   try
      { hitSignal.Fit(&Gaus,"QNRWB","", StartTime, EndTime);}
   catch(...)
      {mf::LogWarning("GausHitFinder") << "Fitter failed finding a hit";}
   
   // ##################################################
   // ### Getting the fitted parameters from the fit ###
   // ##################################################
   Chi2PerNDF = (Gaus.GetChisquare() / Gaus.GetNDF());
   NDF = Gaus.GetNDF();
   
   for(unsigned short ipar = 0; ipar < (3 * nGauss); ++ipar)
      {
      tempGaus1Par.push_back( Gaus.GetParameter(ipar) );
      tempGaus1ParError.push_back( Gaus.GetParError(ipar) );
	 
      }//<---End ipar loop 
   
   Gaus.Delete();
   hitSignal.Delete();
   }//<----End FitGaussians



// --------------------------------------------------------------------------------------------
// Fit Gaussians
// --------------------------------------------------------------------------------------------
void hit::GausHitFinder::ReFitGaussians(std::vector<float> ReSignalVector, std::vector<int> RePeakTime, int mGauss, int ReStartTime, int ReEndTime)
   {
   tempGaus2Par.clear();
   tempGaus2ParError.clear();
   
   std::string eqn2        = "gaus(0)";  // string for equation for gaus fit
   std::stringstream numConv2;
   
   double amplitude2(0);        //fit parameters
   
   int size2 = ReEndTime - ReStartTime;
   
   // #############################################
   // ### If size < 0 then set the size to zero ###
   // #############################################
   if(size2 < 0){size2 = 0;}
   
   // --- TH1D HitSignal ---
   TH1F hitSignal2("hitSignal2","",std::max(size2,1),ReStartTime,ReEndTime);
   hitSignal2.Sumw2();
   
   // #############################
   // ### Filling the histogram ###
   // #############################
   for(int aa = ReStartTime; aa < ReEndTime; aa++)
      {
      hitSignal2.Fill(aa,ReSignalVector[aa]);
      
      if(ReEndTime > 10000){break;}
      }//<---End aa loop
   
   // ############################################
   // ### Building TFormula for basic Gaussian ###
   // ############################################
   eqn2 = "gaus(0)";
   for(int i = 3; i < (mGauss)*3; i+=3)
      {
      eqn2.append("+gaus(");
      numConv2.str("");
      numConv2 << i;
      eqn2.append(numConv2.str());
      eqn2.append(")");
      }
   
   // --- TF1 function for GausHit  ---
   TF1 Gaus2("Gaus2",eqn2.c_str(),0,std::max(size2,1));
	 
   // ### Setting the parameters for the Gaussian Fit ###
   for(int bb = 0; bb < mGauss; bb++)
      {
      float TrialPeak = 0;
      if(bb<mGauss -1){TrialPeak = RePeakTime[bb];}
      else{TrialPeak = RePeakTime[0]+(bb*5);}
      
      if(TrialPeak > ReEndTime){TrialPeak = ReEndTime -1;}
      amplitude2 = 0.5* ReSignalVector[TrialPeak];
      //std::cout<<"fitWidth = "<<fitWidth<<", amplitude2 = "<<ReSignalVector[TrialPeak]<<std::endl;
      Gaus2.SetParameter(3*bb,amplitude2);
      Gaus2.SetParameter(1+(3*bb), TrialPeak);
      Gaus2.SetParameter(2+(3*bb), fitWidth);
      Gaus2.SetParLimits(3*bb, 0.0, 10.0*amplitude2);
      Gaus2.SetParLimits(1+(3*bb), ReStartTime , ReEndTime);
      Gaus2.SetParLimits(2+(3*bb), 0.0, 10.0*fitWidth);
      }//<---End bb loop
      
      
   // ####################################################
   // ### PERFORMING THE TOTAL GAUSSIAN FIT OF THE HIT ###
   // ####################################################
   //hitSignal2.Fit(&Gaus2,"QNRWIB","", ReStartTime, ReEndTime);
   //hitSignal2.Fit(&Gaus2,"QNRWB","", ReStartTime, ReEndTime);
   //std::cout<<"ReStartTime = "<<ReStartTime<<" , ReEndTime = "<<ReEndTime<<std::endl;
   try
      { hitSignal2.Fit(&Gaus2,"QNRW","", ReStartTime, ReEndTime);}
      
   catch(...)
      {mf::LogWarning("GausHitFinder") << "Fitter failed finding a hit";}
   
   // ##################################################
   // ### Getting the fitted parameters from the fit ###
   // ##################################################

   Chi2PerNDF2 = (Gaus2.GetChisquare() / Gaus2.GetNDF());
   NDF2 = Gaus2.GetNDF();
      
   for(unsigned short ipar = 0; ipar < (3 * mGauss); ++ipar)
      {
      tempGaus2Par.push_back( Gaus2.GetParameter(ipar) );
      tempGaus2ParError.push_back( Gaus2.GetParError(ipar) );
	 
      }//<---End ipar loop 
   
   Gaus2.Delete();
   hitSignal2.Delete();
   }//<---End ReFitGaussians



  DEFINE_ART_MODULE(GausHitFinder)

} // end of hit namespace
#endif // GAUSHITFINDER_H
