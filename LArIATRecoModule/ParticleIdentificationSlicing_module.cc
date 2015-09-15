////////////////////////////////////////////////////////////////////////
// Class:       ParticleIdentificationSlicing
// Module Type: producer
// File:        ParticleIdentificationSlicing_module.cc
//
// Generated at Tue Jul 21 10:59:53 2015 by Ryan Linehan using artmod
// from cetpkgsupport v1_08_06.
////////////////////////////////////////////////////////////////////////

//Framework includes
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Framework/Core/FindManyP.h"
#include "art/Framework/Services/Optional/TFileService.h"

//Lariatsoft includes
#include "LArIATDataProducts/WCTrack.h"
#include "LArIATDataProducts/TOF.h"
#include "LArIATDataProducts/MuonRangeStackHits.h"
#include "LArIATDataProducts/AGCounter.h"
#include "RawDataUtilities/TriggerDigitUtility.h"

//C++ includes
#include <vector>
#include <memory>
#include <iostream>

//ROOT includes
#include <TH2F.h>
#include <TH1F.h>
#include <TF1.h>

//LArIAT Includes
#include "Utilities/DatabaseUtilityT1034.h"

class ParticleIdentificationSlicing;

class ParticleIdentificationSlicing : public art::EDProducer {
public:
  explicit ParticleIdentificationSlicing(fhicl::ParameterSet const & p);
  // The destructor generated by the compiler is fine for classes
  // without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  ParticleIdentificationSlicing(ParticleIdentificationSlicing const &) = delete;
  ParticleIdentificationSlicing(ParticleIdentificationSlicing &&) = delete;
  ParticleIdentificationSlicing & operator = (ParticleIdentificationSlicing const &) = delete;
  ParticleIdentificationSlicing & operator = (ParticleIdentificationSlicing &&) = delete;

  // Required functions.
  void produce(art::Event & e) override;

  // Selected optional functions.
  void beginJob() override;
  void beginRun(art::Run & r) override;
  void beginSubRun(art::SubRun & sr) override;
  void endJob() override;
  void endRun(art::Run & r) override;
  void endSubRun(art::SubRun & sr) override;
  void reconfigure(fhicl::ParameterSet const & p) override;
  void respondToCloseInputFile(art::FileBlock const & fb) override;
  void respondToCloseOutputFiles(art::FileBlock const & fb) override;
  void respondToOpenInputFile(art::FileBlock const & fb) override;
  void respondToOpenOutputFiles(art::FileBlock const & fb) override;

  void doThePionMuonSeparation( float thePenetrationDepth,
				float reco_momentum,
				std::vector<float> & pion_muon_likelihood_ratios );

  bool isThereAGoodMuRSTrack( art::Handle< std::vector<ldp::MuonRangeStackHits> > & MuRSColHandle,
			      int & thePenetrationDepth );
  void setPriors( float momentum );
  void pullPriorsFromTable( float momentum );
  void queryDataBaseForMagnetAndEnergy();

  /*
  void getActivePriors( std::string runSetting );
  void getActivePriorsDefault();
  */

private:

  //Run parameters and constants
  float fDistanceTraveled;           //Path length of particle, in meters
  float fSpeedOfLight;                
  float fMagnetSetting;
  float fEnergySetting;
  float fPolaritySetting;
  int   fRun;
  int   fSubRun;
  float fPathLength;
  int   fTriggerWindowFirstTick;
  int   fTriggerWindowLastTick;

  //Database Utility for getting run/subrun beam setting info
  art::ServiceHandle<util::DatabaseUtilityT1034> fDatabaseUtility;

  //Priors storage
  float fMuonPrior;
  float fPionPrior;
  float fKaonPrior;
  float fProtonPrior;
  
  //Gathering information from histograms
  std::vector<float> fPenetrationDepthInfo;
  

  // Declare member data here.
  std::string       fWCTrackModuleLabel;
  std::string       fTOFModuleLabel;
  std::string       fMuRSModuleLabel;
  bool              fVerbose;
  bool              fPlotHistograms;


  //Histogramming and distribution generation
  //parameters and hists
  TH2F*             fPzVsTOF;
  TH1F*             fNTOF;
  TH1F*             fPz;
  TH1F*             fY_Kink;
  TH1F*             fX_Dist;
  TH1F*             fY_Dist;
  TH1F*             fZ_Dist;
  TH1F*             fX_Face_Dist;
  TH1F*             fY_Face_Dist;
  TH1F*             fTheta_Dist;
  TH1F*             fPhi_Dist;
  TH1F*             fTOF;
  TH1F*             fParticleMass;
  TH1F*             fTOFType;
  TH1F*             fMagnetSettingHist;
  TH1F*             fEnergySettingHist;
  TH1F*             fBDiffHist;
  TH1F*             fEDiffHist;

  std::vector<size_t>  fKaonRun;
  std::vector<size_t>  fKaonSubRun;
  std::vector<size_t>  fKaonEvent;

  //Gaussian parameters for mass fits
  float             fPiMuMassMean;
  float             fPiMuMassSigma;
  float             fKaonMassMean;
  float             fKaonMassSigma;
  float             fProtonMassMean;
  float             fProtonMassSigma;

  //Priors for likelihood analysis
  std::map<std::string,float> fPionPriorMap;
  std::map<std::string,float> fMuonPriorMap;
  std::map<std::string,float> fKaonPriorMap;
  std::map<std::string,float> fProtonPriorMap;
  float fPionActivePrior;
  float fMuonActivePrior;
  float fKaonActivePrior;
  float fProtonActivePrior;

  //Misc Parameters
  float fMaxMomentumForPID;  //MeV/c
  float fPiMuLRThreshold;
  
};


ParticleIdentificationSlicing::ParticleIdentificationSlicing(fhicl::ParameterSet const & p)
// :
// Initialize member data here.
{
  // Call appropriate produces<>() functions here.
  this->reconfigure(p);
  
}

void ParticleIdentificationSlicing::produce(art::Event & e)
{
  //  std::cout << "Size of penetrationdepthinfo: " << fPenetrationDepthInfo.size() << std::endl;

  //Get the collection of WCTracks produced by the WCTrackBuilder module
  art::Handle< std::vector<ldp::WCTrack> > WCTrackColHandle;
  e.getByLabel(fWCTrackModuleLabel,WCTrackColHandle);
  
  //Get the collection of TOF objects produced by the TOF module
  art::Handle< std::vector<ldp::TOF> > TOFColHandle;
  e.getByLabel(fTOFModuleLabel,TOFColHandle);

  //Get the collection of MuonRangeStackHits objects produced by the MuonRangeStackHitsBuilder module
  art::Handle< std::vector<ldp::MuonRangeStackHits> > MuRSColHandle;
  e.getByLabel(fMuRSModuleLabel,MuRSColHandle);

  //Assume that there is only one good WCTrack. Identify the WCTrack's momentum.
  if( WCTrackColHandle->size() != 1 ) return;
  float reco_momentum = WCTrackColHandle->at(0).Momentum();
  setPriors( reco_momentum );
  
  //This particle ID work is done under the assumption that there is only one good WCTrack.
  //We set our prior knowledge of particle type (the information about the numbers of
  //pions, protons, etc. in the beam) as a function of the particle momentum. 
  


  /*
  std::cout << "Produce has commenced." << std::endl;

  if( fVerbose ){
    std::cout << "NumTracks in the event: " << WCTrackColHandle->size() << std::endl;
    std::cout << "NumTOF Objects in the event: " << TOFColHandle->size() << std::endl;
    if( TOFColHandle->size() > 0 )
      std::cout << "NumTOFs for first event TOF objct: " << TOFColHandle->at(0).NTOF() << std::endl;
  }

  ///////////////// DISTRIBUTION GENERATION ///////////////////

  //Quick-and-dirty method for matching TOF and WCTracks: good WCTracks only have one track per trigger.
  //I think (?) that only one TOF comes out of the TimeOfFlight module, but just in case that's not true, we'll require
  //that only one TOF can correspond to one trigger
  
  if( TOFColHandle->size() > 1 ) std::cout << "More than one TOF object per event? This is weird." << std::endl;
  if( TOFColHandle->size() > 0 ){
    if( fPlotHistograms )
      fTOFType->Fill(TOFColHandle->at(0).NTOF());
  }

  //Temporary solution - need to address cases where there is more than 1 TOF per event
  if( TOFColHandle->size() > 0 ){
    if((WCTrackColHandle->size() == 1 && TOFColHandle->at(0).NTOF() == 1)){
      
       //Identify Kaons
      if((WCTrackColHandle->at(0).Momentum() < 700 && WCTrackColHandle->at(0).Momentum() > 500 ) &&
	 (TOFColHandle->at(0).SingleTOF(0) < 45 && TOFColHandle->at(0).SingleTOF(0) > 37) ){
	fKaonRun.push_back(size_t(e.run()));
	fKaonSubRun.push_back(size_t(e.subRun()));
	fKaonEvent.push_back(size_t(e.event()));
      }
			  	

      if( fPlotHistograms ){
	fPzVsTOF->Fill(WCTrackColHandle->at(0).Momentum(),TOFColHandle->at(0).SingleTOF(0));
	fPz->Fill(WCTrackColHandle->at(0).Momentum());
	fTOF->Fill(TOFColHandle->at(0).SingleTOF(0));
	fY_Kink->Fill(WCTrackColHandle->at(0).YKink());
	fX_Dist->Fill(WCTrackColHandle->at(0).DeltaDist(0));
	fY_Dist->Fill(WCTrackColHandle->at(0).DeltaDist(1));
	fZ_Dist->Fill(WCTrackColHandle->at(0).DeltaDist(2));
	fX_Face_Dist->Fill(WCTrackColHandle->at(0).XYFace(0));
	fY_Face_Dist->Fill(WCTrackColHandle->at(0).XYFace(1));
	fTheta_Dist->Fill(WCTrackColHandle->at(0).Theta());
	fPhi_Dist->Fill(WCTrackColHandle->at(0).Phi());
	

	if(WCTrackColHandle->at(0).Momentum() < 1000 ){
	  float mass = pow(pow((c*(TOFColHandle->at(0).SingleTOF(0)-11)*1e-9/distance_traveled)*WCTrackColHandle->at(0).Momentum(),2)-pow(WCTrackColHandle->at(0).Momentum(),2),0.5);
	  fParticleMass->Fill(mass);
	  std::cout << "mass: " << mass << std::endl;
	}
      }
    }
  }
  */
  ///////////////// DISTRIBUTION GENERATION ENDING ///////////////////  
  ///////////////// PARTICLE ID Pz vs. TOF/////////////////////////////
  //We take the calculated mass of the particle and find the value of 3
  //p.d.f.s at that mass. Each p.d.f. is a gaussian with the parameters
  //defined in the class definition above, and represents a particle
  //flavor (pi/mu, kaon, proton). Each of these probabilities is equal
  //to the likelihood of the associated particle being the correct one
  //given the mass. We normalize these three probabilities and store them
  //for future use.
  //NOTE: If I can get my mitts on MC beam information, we can also weight
  //these with priors (relative beam composition of pions, protons, kaons,
  //and improve our likelihood estimates)

  //Indices after filling:
  //0: proton likelihood ratio
  //1: kaon likelihood ratio
  //2: pimu likelihood ratio
  //Note that the size might still be zero after the
  //following function if cuts aren't passed. Careful!
  std::vector<float> proton_kaon_pimu_likelihood_ratios;
  std::vector<float> pion_muon_likelihood_ratios;
  doThePiMu_Proton_KaonSeparation( WCTrackColHandle,
				   TOFColHandle,
				   proton_kaon_pimu_likelihood_ratios );
  ///////////////// PARTICLE ID ENDING Pz vs. TOF /////////////////////

  //Now we check to see if the pimu likelihood ratio is above some 
  //parametrized threshold. If it is, run pi/mu separation on it.
  if( proton_kaon_pimu_likelihood_ratios.at(2) > fPiMuLRThreshold ){
    
    ///////////////// PARTICLE ID Pion Vs. Muon /////////////////////////
    // We start with the assumption that there is only one good WCTrack,
    // the one that triggered the event. The trigger waveform times for
    // the wire chambers will therefore have hits at around time tick
    // 136. We then cut on track time, saying that the MuRS track must be
    // found within a reasonable distance of tick 136. This very roughly
    // matches the WCTrack to the MuRS track.
    //
    
    MuRSTrack theGoodMuRSTrack;
    int thePenetrationDepth = 9989;
    bool goodMuRS = isThereAGoodMuRSTrack( MuRSColHandle, thePenetrationDepth );
    if( goodMuRS ){
      doThePionMuonSeparation( thePenetrationDepth, 
			       reco_momentum,
			       pion_muon_likelihood_ratios );
    }
    
  }

  //Now fill in the AuxDetParticle with the likelihoods
  
  



}

//============================================================================================
//Finding the likelihood ratio of pions and muons
void ParticleIdentificationSlicing::doThePionMuonSeparation( float thePenetrationDepth,
							     float reco_momentum,
							     std::vector<float> & pion_muon_likelihood_ratios )
{
  /*
  float P_depth_given_pion = getProbabilityOfDepthGivenPionAtPunchThrough( reco_momentum, thePenetrationDepth );
  float P_pion = getProbabilityOfPionAtPunchThrough( reco_momentum );
  float P_depth_given_muon = getProbabilityOfDepthGivenMuonAtPunchThrough( reco_momentum, thePenetrationDepth );
  float P_muon = getProbabilityOfMuonAtPunchThrough( reco_momentum );
  
  float P_pion_given_depth = (P_depth_given_pion*P_pion)/(P_depth_given_pion*P_pion+P_depth_given_muon*P_muon);
  float P_muon_given_depth = (P_depth_given_muon*P_muon)/(P_depth_given_pion*P_pion+P_depth_given_muon*P_muon);
  
  pion_muon_likelihood_ratios.push_back(P_pion_given_depth,P_muon_given_depth);*/
}


//============================================================================================
//Finding whether there is a MuRSTrack at the appropriate time
bool ParticleIdentificationSlicing::isThereAGoodMuRSTrack( art::Handle< std::vector<ldp::MuonRangeStackHits> > & MuRSColHandle,
							   int & thePenetrationDepth )
{
  //Loop through the MuRS objects
  int counter = 0;
  for( size_t iMuRS; iMuRS < MuRSColHandle->size(); ++iMuRS ){
    ldp::MuonRangeStackHits theMuRS = MuRSColHandle->at(iMuRS);
    for( size_t iTrack = 0; iTrack < theMuRS.NTracks(); ++iTrack ){
      float theTrackTime = theMuRS.GetArrivalTime(iTrack);
      if( theTrackTime >= fTriggerWindowFirstTick && theTrackTime <= fTriggerWindowLastTick ){
	thePenetrationDepth = theMuRS.GetPenetrationDepth(iTrack);
	counter++;
      }
    }
    
  }

  if( counter == 1 ) return true;
  else if( counter == 0 ) return false;
  else{
    std::cout << "I'm not sure what to do here. There were 2+ good murs tracks. Returning true." << std::endl;
    return true;
  }
}

//============================================================================================
//Setting prior information on the likelihood
//of different particle ID hypotheses
void ParticleIdentificationSlicing::setPriors(float momentum)
{
  //We have a magnet and energy setting. Use these with momentum to identify:
  //1. Which priors histogram file's associated vector to use (based on B field and Energy )
  //2. Which momentum range to draw our priors from
  //Now dig through the prior table (which is at the moment unfinished) using the momentum,
  //magnet settings, and energy settings
  pullPriorsFromTable( momentum );
}


//============================================================================================
//Extracts the prior information about beamline populations from a table built from MC simulation
void ParticleIdentificationSlicing::pullPriorsFromTable( float momentum )
{
  //Make sure to use the fMag and fEnergy settings in this!!!



  
  //By the end, I should have a set of scalar priors for this specific event (electron, muon, pion, proton, kaon, etc.)
  

}

//============================================================================================
//Gets the settings and rounds them to the nearest discrete setting that we use
//(example: 99.7 Amps would round to 100 A)
void ParticleIdentificationSlicing::queryDataBaseForMagnetAndEnergy()
{
  fMagnetSetting = std::stod(fDatabaseUtility->GetIFBeamValue("mid_f_mc7an",fRun,fSubRun));
  fEnergySetting = std::stod(fDatabaseUtility->GetIFBeamValue("mid_f_mcenrg",fRun,fSubRun));
  if( fMagnetSetting > 0 ) fPolaritySetting = 1;
  else if( fMagnetSetting < 0 ) fPolaritySetting = -1;
  else fPolaritySetting = 0;

  //Known settings
  int magSettingSize = 5;
  int energySettingSize = 4;
  float magSettings[magSettingSize] = {40,50,60,80,100};
  float energySettings[energySettingSize] = {8,16,32,64};
  
  //Loop through and find the closest settings
  float theTrueMagSetting = 9996;
  float theLowestDifference = 9995;
  float theBDiff = 0;
  for( int iMag = 0; iMag < magSettingSize; ++iMag ){
    if( fabs(magSettings[iMag] - fabs(fMagnetSetting)) < theLowestDifference ){
      theLowestDifference = fabs(magSettings[iMag]-fabs(fMagnetSetting));
      theTrueMagSetting = magSettings[iMag];
      theBDiff = fabs(theTrueMagSetting - fabs(fMagnetSetting) );
    }
  }
  float theTrueEnergySetting = 9994;
  theLowestDifference = 9994;
  float theEDiff = 0;
  for( int iEn = 0; iEn < energySettingSize; ++iEn ){
    if( fabs(energySettings[iEn] - fabs(fEnergySetting)) < theLowestDifference ){
      theLowestDifference = fabs(energySettings[iEn]-fabs(fEnergySetting));
      theTrueEnergySetting = energySettings[iEn];
      theEDiff = fabs(theTrueEnergySetting - fabs(fEnergySetting));
    }
  }
  
  //Sanity checks
  fMagnetSettingHist->Fill(fMagnetSetting);
  fEnergySettingHist->Fill(fEnergySetting);
  fBDiffHist->Fill(theBDiff);
  fEDiffHist->Fill(theEDiff);


  //Assign these true values 
  fMagnetSetting = theTrueMagSetting;
  fEnergySetting = theTrueEnergySetting;

}

/*
//============================================================================================
//Setting the active prior variables from the prior maps set in setPriors() function
void ParticleIdentificationSlicing::getActivePriors( std::string runSetting )
{
  fPionActivePrior = fPionPriorMap.at( runSetting );
  fMuonActivePrior = fMuonPriorMap.at( runSetting );
  fKaonActivePrior = fKaonPriorMap.at( runSetting );
  fProtonActivePrior = fProtonPriorMap.at( runSetting );
}


//============================================================================================
//Setting the active priors to 1 for defualt  
void ParticleIdentificationSlicing::getActivePriorsDefault()
{
  fPionActivePrior = 1;
  fMuonActivePrior = 1;
  fKaonActivePrior = 1;
  fProtonActivePrior = 1;
}
*/


//============================================================================================  
void ParticleIdentifiactionSlicing::doThePiMu_Proton_KaonSeparation( art::Handle< std::vector<ldp::WCTrack> > WCTrackColHandle,
								     art::Handle< std::vector<ldp::TOF> > TOFColHandle,
								     std::vector<float> & proton_kaon_pimu_likelihood_ratios );
{
  if( TOFColHandle->size() > 0 ){
    if((WCTrackColHandle->size() == 1 && TOFColHandle->at(0).NTOF() == 1)){

      //Finding the mass
      if(WCTrackColHandle->at(0).Momentum() < fMaxMomentumForPID ){
	float mass = pow(pow((c*(TOFColHandle->at(0).SingleTOF(0)-fMissingNanoseconds)*1e-9/fDistanceTraveled)*WCTrackColHandle->at(0).Momentum(),2)-pow(WCTrackColHandle->at(0).Momentum(),2),0.5);
      }
      else return;   //If mass is too high, it's hard to disambiguate protons, pi/mu, and kaons

      //Finding values of pdf for mass given proton, kaon, pi/mu distributions
      float proton_prob = fProtonActivePrior*(1/pow(2*3.1415926,0.5)/fProtonMassSigma)*exp(-0.5*pow((mass-fProtonMassMean)/fProtonMassSigma,2));
      float kaon_prob = fKaonActivePrior*(1/pow(2*3.1415926,0.5)/fKaonMassSigma)*exp(-0.5*pow((mass-fKaonMassMean)/fKaonMassSigma,2));
      float pimu_prob = (fPionActivePrior+fMuonActivePrior)*(1/pow(2*3.1415926,0.5)/fPiMuMassSigma)*exp(-0.5*pow((mass-fPiMuMassMean)/fPiMuMassSigma,2));
      
      //These ^ are likelihoods, so find the likelihood ratio of each to the total
      float proton_likelihood = proton_prob/(proton_prob+kaon_prob+pimu_prob);
      float kaon_likelihood = kaon_prob/(proton_prob+kaon_prob+pimu_prob);
      float pimu_likelihood = pimu_prob/(proton_prob+kaon_prob+pimu_prob);
      
      proton_kaon_pimu_likelihood_ratios.push_back(proton_likelihood);
      proton_kaon_pimu_likelihood_ratios.push_back(kaon_likelihood);
      proton_kaon_pimu_likelihood_ratios.push_back(pimu_likelihood);
    }
  }
}

//============================================================================================  
void ParticleIdentificationSlicing::beginJob()
{
  // Implementation of optional member function here.
  if( fPlotHistograms ){
    art::ServiceHandle<art::TFileService> tfs;
    fPzVsTOF = tfs->make<TH2F>("PzVsTOF","Pz vs. Time of Flight",160,0,1600,60,20,80);
    fNTOF = tfs->make<TH1F>("NTOF","Number of TOF values per TOF object",10,0,10);
    fPz = tfs->make<TH1F>("Reco_Pz","Reconstructed momentum",180,0,1800);
    fTOF = tfs->make<TH1F>("Reco_TOF","Reconstructed Time of Flight",70,20,90);
    fY_Kink = tfs->make<TH1F>("Y_Kink","Angle between US/DS tracks in Y direction (degrees)",100,-5*3.1415926/180,5*3.141592654/180);
    fX_Dist = tfs->make<TH1F>("X_Dist","X distance between US/DS tracks at midplane (mm)",120,-60,60);
    fY_Dist = tfs->make<TH1F>("Y_Dist","Y distance between US/DS tracks at midplane (mm)",120,-60,60);
    fZ_Dist = tfs->make<TH1F>("Z_Dist","Z distance between US/DS tracks at midplane (mm)",120,-60,60);
    fX_Face_Dist = tfs->make<TH1F>("X_Face","X Location of Track's TPC Entry (mm)",800,-200,600);
    fY_Face_Dist = tfs->make<TH1F>("Y_Face","Y Location of Track's TPC Entry (mm)",800,-400,400);
    fTheta_Dist = tfs->make<TH1F>("Theta","Track Theta (w.r.t. TPC Z axis), (radians),",100,0,0.2);
    fPhi_Dist = tfs->make<TH1F>("Phi","Track Phi (w.r.t. TPC X axis), (radians)",100,0,6.28318);
    fTOFType = tfs->make<TH1F>("TOFType","Number of valid times of flight per event",10,0,10);
    fParticleMass = tfs->make<TH1F>("Mass","Particle Mass",50,0,2000);
    fMagnetSettingHist = tfs->make<TH1F>("MagnetSetting","Magnet Setting",1,0,-1);
    fEnergySettingHist = tfs->make<TH1F>("EnergySetting","Energy Setting",1,0,-1);
    fBDiffHist = tfs->make<TH1F>("BDiff","Difference between B field setting and closest discrete one",1,0,-1);
    fEDiffHist = tfs->make<TH1F>("EDiff","Difference between Energy setting and closest discrete one",1,0,-1);

    
    fPz->GetXaxis()->SetTitle("Reconstructed momentum (MeV/c)");
    fPz->GetYaxis()->SetTitle("Tracks per 10 MeV/c");
    fY_Kink->GetXaxis()->SetTitle("Reconstructed y_kink (radians)");
    fY_Kink->GetYaxis()->SetTitle("Tracks per 0.000872 radians");
    fX_Dist->GetXaxis()->SetTitle("X distance between US and DS track ends");
    fX_Dist->GetYaxis()->SetTitle("Tracks per 1 mm");
    fY_Dist->GetXaxis()->SetTitle("Y distance between US and DS track ends");
    fY_Dist->GetYaxis()->SetTitle("Tracks per 1 mm");
    fZ_Dist->GetXaxis()->SetTitle("Z distance between US and DS track ends");
    fZ_Dist->GetYaxis()->SetTitle("Tracks per 1 mm");
    fX_Face_Dist->GetXaxis()->SetTitle("X (mm)");
    fX_Face_Dist->GetYaxis()->SetTitle("Tracks per 1 mm");
    fY_Face_Dist->GetXaxis()->SetTitle("Y (mm)");
    fY_Face_Dist->GetYaxis()->SetTitle("Tracks per 1 mm");
    fTheta_Dist->GetXaxis()->SetTitle("Theta (radians)");
    fTheta_Dist->GetYaxis()->SetTitle("Tracks per .002 radians");
    fPhi_Dist->GetXaxis()->SetTitle("Phi (radians)");
    fPhi_Dist->GetYaxis()->SetTitle("Tracks per 0.0628 radians");
    fParticleMass->GetXaxis()->SetTitle("Particle Mass (MeV/c^2)");
    fParticleMass->GetYaxis()->SetTitle("Counts");
  }

 
  //Setting the distance traveled by a particle via the geometry service
  fDistanceTraveled = 6.7; //Temporarily hacked together

}

void ParticleIdentificationSlicing::beginRun(art::Run & r)
{
  // getActivePriorsDefault();
}

void ParticleIdentificationSlicing::beginSubRun(art::SubRun & sr)
{
  // Implementation of optional member function here.
  fRun = sr.run();
  fSubRun = sr.subRun();

  //Draw information from database on the beam settings
  queryDataBaseForMagnetAndEnergy();
}

void ParticleIdentificationSlicing::endJob()
{
  // Implementation of optional member function here.
  for( size_t iKaon = 0; iKaon < fKaonRun.size() ; ++iKaon ){
    std::cout << "Kaon Run: " << fKaonRun.at(iKaon) << ", SubRun: " << fKaonSubRun.at(iKaon) << ", Event: " << fKaonEvent.at(iKaon) << std::endl;
  }
  
  //Fitting the mass histogram in the appropriate regions to get the
  //parameters used for likelihood estimation of particle flavor
  TF1 * f1 = new TF1("f1","[2]/pow(2*3.14159265,0.5)*exp(-pow((x-[0])/[1],2)/2)",0,400);
  TF1 * f2 = new TF1("f2","[2]/pow(2*3.14159265,0.5)*exp(-pow((x-[0])/[1],2)/2)",400,700);
  TF1 * f3 = new TF1("f3","[2]/pow(2*3.14159265,0.5)*exp(-pow((x-[0])/[1],2)/2)",700,1200);
  f1->SetParameter(0,140);
  f1->SetParameter(1,50);
  f1->SetParameter(2,50);
  f2->SetParameter(0,490);
  f2->SetParameter(1,20);
  f2->SetParameter(2,5);
  f3->SetParameter(0,940);
  f3->SetParameter(1,100);
  f3->SetParameter(2,20);
  std::cout << "Fitting to Pi/Mu peak: " << std::endl;
  fParticleMass->Fit("f1","R");
  std::cout << "Fitting to Kaon peak: " << std::endl;
  fParticleMass->Fit("f2","R");
  std::cout << "Fitting to Proton peak: " << std::endl;
  fParticleMass->Fit("f3","R");
  


}

void ParticleIdentificationSlicing::endRun(art::Run & r)
{
  // Implementation of optional member function here.
}

void ParticleIdentificationSlicing::endSubRun(art::SubRun & sr)
{
  // Implementation of optional member function here.
}

void ParticleIdentificationSlicing::reconfigure(fhicl::ParameterSet const & p)
{
  // Implementation of optional member function here.
  fPathLength                   =p.get< float >("TrajectoryPathLength",6.7);
  fSpeedOfLight                 = 3e+8;
  fTriggerWindowFirstTick       =p.get< int >("TriggerWindowFirstTick",135);
  fTriggerWindowLastTick        =p.get< int >("TriggerWindowLastTick",138);

  fPenetrationDepthInfo         =p.get<std::vector<float> >("PDepthInfo");

  fWCTrackModuleLabel           =p.get< std::string >("WCTrackModuleLabel");
  fTOFModuleLabel               =p.get< std::string >("TOFModuleLabel");
  fMuRSModuleLabel              =p.get< std::string >("MuRSModuleLabel");
  fVerbose                      =p.get< bool >("Verbose");
  fPlotHistograms               =p.get< bool >("PlotHistograms");

  fPiMuMassMean                 =p.get< float >("PiMuMassMean",1.892e+2);
  fPiMuMassSigma                =p.get< float >("PiMuMassSigma",4.646e+1);
  fKaonMassMean                 =p.get< float >("KaonMassMean",5.56e+2);
  fKaonMassSigma                =p.get< float >("KaonMassSigma",8.39e+1);
  fProtonMassMean               =p.get< float >("ProtonMassMean",9.16e+2);
  fProtonMassSigma              =p.get< float >("ProtonMassSigma",9.85e+1);

  fMaxMomentumForPID            =p.get< float >("MaxMomentumForPID",1000);
  fPiMuLRThreshold              =p.get< float >("PiMuLikelihoodRatioThreshold",0.5);
  

  
}

void ParticleIdentificationSlicing::respondToCloseInputFile(art::FileBlock const & fb)
{
  // Implementation of optional member function here.
}

void ParticleIdentificationSlicing::respondToCloseOutputFiles(art::FileBlock const & fb)
{
  // Implementation of optional member function here.
}

void ParticleIdentificationSlicing::respondToOpenInputFile(art::FileBlock const & fb)
{
  // Implementation of optional member function here.
}

void ParticleIdentificationSlicing::respondToOpenOutputFiles(art::FileBlock const & fb)
{
  // Implementation of optional member function here.
}

DEFINE_ART_MODULE(ParticleIdentificationSlicing)
