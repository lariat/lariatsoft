//////////////////////////////////////////////////////////////////////////
// Class:       XSAnalysis
// Module Type: analyzer
// File:        XSAnalysis_module.cc
//
// Generated at Wed Aug  5 12:52:45 2015 by Ryan Linehan using artmod
// from cetpkgsupport v1_08_06.
// Readapted by Irene Nutini for reading calorimetry from the TPC tracks and with an updated XS calculation routine 
// Cleaned up by Elena Gramellini:
// To Do:
// [ x ] Take external cut configuration
// [ x ] Avoid if repetition
// [   ] Code clean up
// [ x ] From Producer (WHY?!?!) to Analyzer
// [   ] Decent comments
// [   ] Modify to include MC switch
//       [   ] Exclude WC information
// [ x ] Read Associations
// [   ] Write an Exception class like SelectionTool/ERTool/Base/ERException.cxx 
////////////////////////////////////////////////////////////////////////

// ########################
// ### LArSoft includes ###
// ########################
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Hit.h"

// #######################
// ### LArIAT includes ###
// #######################
#include "Utilities/DatabaseUtilityT1034.h"
#include "LArIATRecoAlg/TriggerFilterAlg.h"
#include "LArIATDataProducts/WCTrack.h"

// ##########################
// ### Framework includes ###
// ##########################
#include "canvas/Utilities/InputTag.h"
#include "art/Framework/Core/EDAnalyzer.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/FindOneP.h" 
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art_root_io/TFileService.h"
#include "art_root_io/TFileDirectory.h"
#include "art/Framework/Services/Registry/ServiceHandle.h" 
#include "canvas/Persistency/Common/Ptr.h" 
#include "canvas/Persistency/Common/PtrVector.h" 
#include "canvas/Utilities/InputTag.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib/maybe_ref.h"
#include "fhiclcpp/ParameterSet.h" 

// ########################
// ### LArSoft includes ###
// ########################
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "larcoreobj/SimpleTypesAndConstants/RawTypes.h" 
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/CryostatGeo.h"
#include "larcorealg/Geometry/TPCGeo.h"
#include "larcorealg/Geometry/PlaneGeo.h"
#include "larcorealg/Geometry/WireGeo.h"
#include "lardataobj/RecoBase/Wire.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardata/DetectorInfoServices/LArPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "lardataobj/RawData/ExternalTrigger.h"
#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RawData/raw.h"
#include "larsim/MCCheater/BackTracker.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "larevt/Filters/ChannelFilter.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "lardataobj/AnalysisBase/ParticleID.h"
#include "LArIATDataProducts/WCTrack.h"
#include "RawDataUtilities/TriggerDigitUtility.h"

// #####################
// ### ROOT includes ###
// #####################
#include <TH1F.h>
#include <TF1.h>
#include <TH2F.h>
#include <TGraph.h>

// ####################
// ### C++ includes ###
// ####################
#include <map>
#include <memory>


class XSAnalysis;

class XSAnalysis : public art::EDAnalyzer {
public:
  explicit XSAnalysis(fhicl::ParameterSet const & p);
  // The destructor generated by the compiler is fine for classes
  // without bare pointers or other resource use.
  
  // Plugins should not be copied or assigned.
  XSAnalysis(XSAnalysis const &) = delete;
  XSAnalysis(XSAnalysis      &&) = delete;
  XSAnalysis & operator = (XSAnalysis const &) = delete;
  XSAnalysis & operator = (XSAnalysis      &&) = delete;
  
  // Required functions.
  void analyze(art::Event const & e) override;
  
  // Selected optional functions.
  void beginJob   (                ) override;
  void endJob     (                ) override;
  void reconfigure(fhicl::ParameterSet const & p) ;

  void calcXS (float wcTrackMomentum,
	       std::vector<art::Ptr<anab::Calorimetry> > caloTPCTrack);
  
			
  
/*
  PSEUDO CODE
  findTheTPCTrack gives back the id of TPC track  
  calcXS calculates the XS
*/


private:
  // Declare member data here.

  //////////////////
  //Producers' names
  //////////////////
  std::string      fSlicerSourceLabel;
  std::string      fWCTrackBuilderLabel;
  std::string      fCalorimetryModuleLabel;
  std::string      fTPCTrackHandleLabel;
  std::string      fNonStoppingTPCTrackBuildLabel;

  //double fupstreamZPosition; // unused
  // The parameters we'll read from the .fcl file.
  //  bool  fAmIData;                /// amIData True if real Data, False if MC
  bool  fUseNonStoppingOnly;

  /////////////
  //Parameters
  /////////////
  int   run;			//<---Run Number
  int   subrun;			//<---SubRun Number
  int   event;			//<---Event Number
  
  int   fNEvents;               /// Doxygen comment
  int   fEventCounter;          /// Doxygen comment
  float fEpsUpstream;           /// Max z of most upstream point of a track to be tagged as a primary 
  float fRad2Deg;               /// Doxygen comment
  bool  fVerbose;		/// Doxygen comment
  bool  fXSVerbose;		/// Doxygen comment
  float fRecombinationConstant; /// Doxygen comment
  float fMass;                  /// Mass of the particle we want to calculate the XS for
  float fInitialLostEn;         /// Energy Lost by the particle at in the Ar before the TPC
  float fSlabThick;             /// Doxygen comment
  int   ntracks_reco;
  int   nwctrks;

  
  
  //################################################
  //### Histogram declaration for whole analysis ###
  //################################################
  //////////////////////////////////////////////////////////////
  // ~~*~~ Cross section histos - the important ones!!! ~~*~~ //
  //////////////////////////////////////////////////////////////
  TH1F*  fInteractions;                     /// Doxygen comment
  TH1F*  fIncident;                         /// Doxygen comment
  TH1F*  fCrossSection;                     /// Doxygen comment

  ///////////////////
  //Checks Histograms
  ///////////////////
  TH1F*  fInitialEnergy;                    /// Doxygen comment
  TH2F*	 fdEdxVSrSel;			    /// Doxygen comment
  TH1F*	 fEdep;			            /// Doxygen comment

  TH1F*  fTPCTrackEventNumbers;             /// Doxygen comment  
  TH1F*  fWCTrackEventNumbers;	            /// Doxygen comment

  // === Storing the tracks Calorimetry Information
  //int    trkhits[2]; // unused
  //double trkke[2]; // unused
  //double trkdedx[2][1000]; // unused and redefined elsewhere
  //double trkrr[2][1000]; // unused and redefined elsewhere
  //double trkpitchhit[2][1000]; // unused and redefined elsewhere

  //TH1F  *fTotNTrack; // unused
  //TH1F  *fStoppingTrack; // unused
  //TH1F  *fDifferenceTrack; // unused
  //TH1F  *fNTrackPassedCuts; // unused
  
  //TH1F *fdEdx; // unused
  //TH2F *fdEdxRange; // unused
  //TH1F *fFitPar; // unused

  //double fLowLimitStop; // unused
  //double fUpLimitStop; // unused



};


XSAnalysis::XSAnalysis(fhicl::ParameterSet const & p)
  :
  EDAnalyzer(p)
// Initialize member data here.
{
  // Call appropriate produces<>() functions here.
  this->reconfigure( p );
  fEventCounter = 1;
  fRad2Deg = 180.0/3.141592654;
}

void XSAnalysis::analyze(art::Event const & e)
{
  /**
     This functions expects events 
     with 1 wcTrack
   */

  //#####################################
  //### Let's start by retrieving all ###
  //### the reconstructed objects     ###
  //### we will need for the analysis ###
  //#####################################
  
 
  //Retrieving the identifiers
  run = e.run();          /// Run     Number
  subrun = e.subRun();    /// Sub-Run Number
  event = e.id().event(); /// Event   Number
  float wcTrackMomentum;
  
  
  ///Retrieving the WCtracks from the sliced event
  
  art::Handle< std::vector<ldp::WCTrack> > wcTrackHandle;        
  std::vector<art::Ptr<ldp::WCTrack> >     wctrack;                /// Define wctrack as a pointer to ldp::WCTrack
  if(!e.getByLabel(fWCTrackBuilderLabel, wcTrackHandle)) return;   /// If there are no wire chamber tracks for the right label, return 
  
  art::fill_ptr_vector(wctrack, wcTrackHandle);			   /// Filling the wctrack from the wcTrackHandle 	       
  nwctrks = wctrack.size();	
  if ( nwctrks !=1 ) return;                                       /// If there are more than 1 wire chamber tracks, return 
  
  // Take the wcTrack momentum
  wcTrackMomentum = wctrack[0]->Momentum();
  
  
  
  ///Retrieving the TPC tracks from the sliced event
  if(fUseNonStoppingOnly) {
    auto tpcTrackHandle = e.getValidHandle<art::PtrVector<recob::Track> >(fNonStoppingTPCTrackBuildLabel);
    if( !tpcTrackHandle.isValid() ) return;
    
    auto tracklist =  *tpcTrackHandle;   
    ntracks_reco = tracklist.max_size();                                    /// Store the number of tracks per event   
    if ( !ntracks_reco ){ return; }       /// If there are no TPC tracks in general, return
    ///Fill hist info with total number of TPC and WC tracks per event
    fTPCTrackEventNumbers->SetBinContent(fEventCounter,ntracks_reco);
    fWCTrackEventNumbers->SetBinContent(fEventCounter,nwctrks);
    
    
    // Get the rigth calorimtry: the calorimetry is associated with the TPC tracks
    // not with the Stopping Tracks. We'll use the key of the non stopping track to
    // understand which calo obj we need and a mock tpc handle to retreive the right calo handler.
    art::Handle< std::vector<recob::Track> > mocktpcTrackHandle;
    if(!e.getByLabel(fTPCTrackHandleLabel, mocktpcTrackHandle)) return;  
    ///Association between Calorimetry objects and Tracks
    art::FindManyP<anab::Calorimetry> fmcal(mocktpcTrackHandle, e, fCalorimetryModuleLabel);
    
    for ( auto const& thisTrack : tracklist )
      {
	unsigned int indexTPCTrack = thisTrack.key();
	// Get the calorimetry associated with the right tpcTrack
	std::vector<art::Ptr<anab::Calorimetry> > caloTPCTrack = fmcal.at(indexTPCTrack);//here the track index is needed
	//// Calculate the XS
	calcXS(wcTrackMomentum, caloTPCTrack); 
      }
    
                  

  }else{
    art::Handle< std::vector<recob::Track> > tpcTrackHandle;        
    std::vector<art::Ptr<recob::Track> >     tracklist;                /// Define tracklist as a pointer to recob::Track
    if(!e.getByLabel(fTPCTrackHandleLabel, tpcTrackHandle)) return;   /// If there are no tpc chamber tracks for the right label, return 
    
    art::fill_ptr_vector(tracklist, tpcTrackHandle);			   /// Filling the wctrack from the tpcTrackHandle 	      
    
    ntracks_reco = tracklist.size();                                    /// Store the number of tracks per event   
    if ( !ntracks_reco ){ return; }       /// If there are no TPC tracks in general, return
    ///Fill hist info with total number of TPC and WC tracks per event
    fTPCTrackEventNumbers->SetBinContent(fEventCounter,ntracks_reco);
    
    fWCTrackEventNumbers->SetBinContent(fEventCounter,nwctrks);
    
    art::FindManyP<anab::Calorimetry> fmcal(tpcTrackHandle, e, fCalorimetryModuleLabel);
    
    for ( auto const& thisTrack : tracklist )
      {
	unsigned int indexTPCTrack = thisTrack.key();
	// Get the calorimetry associated with the right tpcTrack
	std::vector<art::Ptr<anab::Calorimetry> > caloTPCTrack = fmcal.at(indexTPCTrack);//here the track index is needed
	//// Calculate the XS
	calcXS(wcTrackMomentum, caloTPCTrack); 
      }
  }
  
  
  fEventCounter++;
    
  // If you feel like talking, tell me all the info you've got
  if(fVerbose)
    {
      std::cout<<std::endl;
      std::cout<<"========================================="<<std::endl;
      std::cout<<"Run = "<<run<<", SubRun = "<<subrun<<", Evt = "<<event<<std::endl;
      std::cout<<"Number of  WC tracks " << nwctrks      <<std::endl;
      std::cout<<"Number of TPC tracks " << ntracks_reco <<std::endl;
      std::cout<<"========================================="<<std::endl;
      std::cout<<std::endl;
    } 
  
}

//===========================================================================================
					
//===========================================================================================
void XSAnalysis::beginJob()
{
  // Implementation of optional member function here.
  art::ServiceHandle<art::TFileService> tfs;

  /////////////// DEFINITIONS ////////////////
  fInteractions = tfs->make<TH1F>("Interactions","Number of pion interactions",20,0.,1000.);
  fIncident     = tfs->make<TH1F>("Incident"    ,"Number of pions incident"   ,20,0.,1000.);
  fCrossSection = tfs->make<TH1F>("CrossSection","Cross Section"              ,20,0.,1000.);  

  fInitialEnergy  = tfs->make<TH1F>("InitialEnergy"            ,"Initial Kinetic Energy of Particle in TPC"                     ,180,0,1800);
  fdEdxVSrSel     = tfs->make<TH2F>("dEdxVSresRangeNonStopping", "dEdxVSresRange for the tpc track selected  not stopping"      , 100,0.,100., 100,0.,20.);    
  fEdep           = tfs->make<TH1F>("EnergyAtTrkEnd"           , "energy at end track (with neg dEdx correction and selections)", 100, 0., 1000.);
  
  fTPCTrackEventNumbers             = tfs->make<TH1F>("TPCTrackEventNumbers","Number of TPC tracks as a function of event",fNEvents,0,fNEvents);
  fWCTrackEventNumbers              = tfs->make<TH1F>("WCTrackEventNumbers" ,"Number of  WC tracks as a function of event",fNEvents,0,fNEvents);

  
  ///////////////// LABELING //////////////////
  fInteractions->GetXaxis()->SetTitle("Kinetic Energy (Mev)");
  fInteractions->GetYaxis()->SetTitle("Number of Interactions");
  fIncident->GetXaxis()->SetTitle("Kinetic Energy (MeV)");
  fIncident->GetYaxis()->SetTitle("Number of Incident Particles");
  fCrossSection->GetXaxis()->SetTitle("Kinetic Energy, MeV");
  fCrossSection->GetYaxis()->SetTitle("Cross Section, barns");

  fInitialEnergy->GetXaxis()->SetTitle("Kinetic Energy (MeV)");
  fInitialEnergy->GetYaxis()->SetTitle("Counts");

  fTPCTrackEventNumbers->GetXaxis()->SetTitle("Event number");
  fTPCTrackEventNumbers->GetYaxis()->SetTitle("Number of TPC Tracks");
  fWCTrackEventNumbers->GetXaxis()->SetTitle("Event number");
  fWCTrackEventNumbers->GetYaxis()->SetTitle("Number of WC Tracks");

  
}

void XSAnalysis::endJob()
{
  /// Create the cross section from the incident and interaction plots
  float rho        = 1400;     ///kg/m^3
  float molar_mass = 39.9;     ///g/mol
  float g_per_kg   = 1000;     /// Conversion factor
  float avogadro   = 6.02e+23; ///number/mol
  float number_density = rho*g_per_kg/molar_mass*avogadro;
  float slab_width = fSlabThick/100.;//in m

  std::cout << "Now start the xs calculation" << std::endl;
  
  for( int iBin = 1; iBin <= fInteractions->GetNbinsX(); ++iBin ){
    float crossSection=0.;
    if( fIncident->GetBinContent(iBin) == 0 ) continue; /// We don't have enought stat to say anything for that bin
    else{
      crossSection= fInteractions->GetBinContent(iBin)/fIncident->GetBinContent(iBin)/number_density/slab_width;
      crossSection /= 1e-28; //To put this into barns
      fCrossSection->SetBinContent(iBin,crossSection);

      //Now, let's take care of the erros
      float numError = pow(fInteractions->GetBinContent(iBin),0.5);
      float denomError = pow(fIncident->GetBinContent(iBin),0.5);
      float totalError = 0.;
      if(fInteractions->GetBinContent(iBin) == 0) {continue;} //totalError=0.;}
      else{ 
	totalError = fInteractions->GetBinContent(iBin)/fIncident->GetBinContent(iBin)*pow((pow(numError/fInteractions->GetBinContent(iBin),2)
		     + pow(denomError/fIncident->GetBinContent(iBin),2)),0.5)/number_density/slab_width/1e-28;
       }
      fCrossSection->SetBinError(iBin,totalError); 
    }
    std::cout << iBin 
	      << " " << fIncident    ->GetBinContent(iBin) 
	      << " " << fInteractions->GetBinContent(iBin) 
	      << " " << fCrossSection->GetBinContent(iBin)
	      << " " << fCrossSection->GetBinError(iBin)  << std::endl;
  }  
}

void XSAnalysis::reconfigure(fhicl::ParameterSet const & p)
{
  // Implementation of optional member function here.
  fSlicerSourceLabel             = p.get<std::string>("SourceLabel"           , "SlicerInput"    );
  fWCTrackBuilderLabel           = p.get<std::string>("WCTrackBuilderLabel"   , "wctrack"        );
  fNonStoppingTPCTrackBuildLabel = p.get<std::string>("TPCTrackBuilderLabel"  , "NonStoppingTrk" );
  fCalorimetryModuleLabel        = p.get<std::string>("CalorimetryModuleLabel", "calo"           );
  fTPCTrackHandleLabel           = p.get<std::string>("TPCTrackHandle"        , "pmtrack"           );
  fNEvents                 = p.get< int   >("NumEvents");
  fEpsUpstream             = p.get< float >("UpstreamEpsilon");
  fRecombinationConstant   = p.get< float >("RecombinationFactor");
  fMass                    = p.get< float >("Mass"             ,493.677);
  fInitialLostEn           = p.get< float >("InitialLostEn"    ,  8.6  );
  fSlabThick  		   = p.get< float >("ThinSlabThickness");//in cm
  //fAmIData                 = p.get< bool  >("AmIData"            , true );
  fUseNonStoppingOnly      = p.get< bool  >("UseNonStoppingOnly" , false );
  fVerbose                 = p.get< bool  >("Verbose"            , false);
  fXSVerbose               = p.get< bool  >("XSVerbose"          , true );
}

void XSAnalysis::calcXS (float momentum,
			 std::vector<art::Ptr<anab::Calorimetry> > caloTPCTrack)
{
  /** This fuction takes one wcTrack to define the initial energy and 
      the calorimetry associated with one tpcTrack and calculated the XS
  */
  
  ///Get momentum from WCTrack and convert it to energy UNDER THE ASSUMPTION THAT THIS IS A KAON
  ///Also subtract 8.6 MeV for first approximation of TPC entry (MeV This needs checking!!!)
  float entryTPCEnergyLoss = fInitialLostEn; /// MeV This needs changing!!!
  float kineticEnergy = pow(momentum*momentum+fMass*fMass,0.5)-fMass;
  kineticEnergy -= entryTPCEnergyLoss; 
  // You have calculated the right initial energy. Good Boy.
  
  ///Fill initial energy information
  fInitialEnergy->Fill(kineticEnergy);
  if (fXSVerbose) std::cout << "Initial Kinetic energy " << kineticEnergy << std::endl;


  //Vectors with calo info of the matched tpc track
  double trkdedx[1000]={0.};
  double trkrr[1000]={0.};
  double trkpitchhit[1000]={0.};

    
  art::Ptr<anab::Calorimetry> theInterestingCaloObj;
  //Keep the calorimetry only on the collection plane
  for(const auto &tmpCalo : caloTPCTrack){
    if (!tmpCalo->PlaneID().Plane) theInterestingCaloObj = tmpCalo;
  }
  
  //std::cout << "Number of calo hits on collection plane for this track " << theInterestingCaloObj->dEdx().size() << std::endl;
  // Note: the NCalo can change if the track is thru going. NEEDS CHECKING!!!!  
  size_t NCalo =  theInterestingCaloObj->dEdx().size();
  int    negPitch=0;
  int    negdEdx=0;
  
  for (size_t k = 0; k<theInterestingCaloObj->dEdx().size(); ++k){
    // ### Recording the dE/dX information for this calo point along the track in this plane ###
    trkdedx[k] = theInterestingCaloObj->dEdx()[k];
    
    if(theInterestingCaloObj->TrkPitchVec()[k]< 0.) negPitch++; //std::cout << "Negative pitch " << std::endl;
    if(theInterestingCaloObj->dEdx()[k]< 0.) negdEdx++;
    // ### Recording the residual range for this calo point along the track in this plane ###
    trkrr[k] = theInterestingCaloObj->ResidualRange()[k];
    
    // ### Recording the pitch of this calo point along the track in this plane ###
    trkpitchhit[k] = theInterestingCaloObj->TrkPitchVec()[k];
    //std::cout << "deltaE per each point " << (theInterestingCaloObj->dEdx()[k])*(theInterestingCaloObj->TrkPitchVec()[k]) << std::endl;
  }
  
  
  std::cout << " Now start filling histos for xs calc" << std::endl;
  
  size_t caloK = 0;
  float totalLostEnergy = 0.;
  while(caloK < NCalo)       
    {
      fIncident->Fill(kineticEnergy);
      fdEdxVSrSel->Fill(trkrr[caloK],trkdedx[caloK]);
      if(caloK == NCalo-1) 
	{
	  fInteractions->Fill(kineticEnergy);
	} 
      
      float dEdxTemp=0.;
      dEdxTemp=trkdedx[caloK];
      float lostEnergy = dEdxTemp*trkpitchhit[caloK]*fRecombinationConstant;
      totalLostEnergy += lostEnergy;
      kineticEnergy -= lostEnergy; 
      caloK++;
    }
  fEdep->Fill(kineticEnergy);
  
}

		


DEFINE_ART_MODULE(XSAnalysis)
