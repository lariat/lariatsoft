////////////////////////////////////////////////////////////////////////
// Class:       ParticleFilter
// Module Type: filter
// File:        ParticleFilter_module.cc
//
// Generated at Thur Dec 10 2015 by Irene Nutini using artmod
// from cetpkgsupport v1_08_06.
//Following KaonFilter - PiMuFilter structure
//Using BeamlinePID object to identify and select different particle species
//From the fcl file of the filter you can choose the PDG of the particle you want to select and filter
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "LArIATDataProducts/WCTrack.h"
#include "LArIATDataProducts/TOF.h"
#include "LArIATDataProducts/AuxDetParticleID.h"
#include <TH2F.h>
#include "art/Framework/Services/Optional/TFileService.h"

#include <iostream>
#include <memory>

class ParticleFilter;

class ParticleFilter : public art::EDFilter {
public:
  explicit ParticleFilter(fhicl::ParameterSet const & p);
  // The destructor generated by the compiler is fine for classes
  // without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  ParticleFilter(ParticleFilter const &) = delete;
  ParticleFilter(ParticleFilter &&) = delete;
  ParticleFilter & operator = (ParticleFilter const &) = delete;
  ParticleFilter & operator = (ParticleFilter &&) = delete;

  // Required functions.
  bool filter(art::Event & e) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;
  void reconfigure(fhicl::ParameterSet const & p) override;
  void respondToCloseInputFile(art::FileBlock const & fb) override;
  void respondToCloseOutputFiles(art::FileBlock const & fb) override;
  void respondToOpenInputFile(art::FileBlock const  &fb) override;
  void respondToOpenOutputFiles(art::FileBlock const & fb) override;

private:

  // Declare member data here.
  TH2F* fParticlePzVsTOF;
  TH2F* fPzVsTOF;
  
  std::string fParticleIDModuleLabel;
  std::string fWCTrackModuleLabel;
  std::string fTOFModuleLabel;
  
  double fParticlePDG;
  double fParticleProbCutOff;

};


// ########################################
// ### Calling the reconfigure function ###
// ########################################
ParticleFilter::ParticleFilter(fhicl::ParameterSet const & p)
// :
// Initialize member data here.
{
  // Call appropriate produces<>() functions here.
  this->reconfigure(p);

}

// ############################
// ### Reconfigure Function ###
// ############################
void ParticleFilter::reconfigure(fhicl::ParameterSet const & p)
{
  // Implementation of optional member function here.
  fParticleProbCutOff    = p.get<   float   >("ParticleProbabilityThreshold",0.5);
  fParticleIDModuleLabel = p.get<std::string>("ParticleIDModuleLabel");
  fWCTrackModuleLabel    = p.get<std::string>("WCTrackModuleLabel");
  fTOFModuleLabel        = p.get<std::string>("TOFModuleLabel");
  fParticlePDG           = p.get<   float   >("ParticlePDG",211); //default value is for piplus

}


// ##################
// ### Event Loop ###
// ##################
bool ParticleFilter::filter(art::Event & e)
{
  std::cout<<"Puppaaaaaaaaaaaaaaaaaa"<<"\n";
  // Implementation of required member function here.
  //Retrieving the Particle IDs from the event record
  art::Handle< std::vector<ldp::AuxDetParticleID> > particleIDCol;
  e.getByLabel(fParticleIDModuleLabel,particleIDCol);
  
  
  //Get the collection of WCTracks produced by the WCTrackBuilder module
  art::Handle< std::vector<ldp::WCTrack> > WCTrackColHandle;
  e.getByLabel(fWCTrackModuleLabel,WCTrackColHandle);
  
  
  //Get the collection of TOF objects produced by the TOF module
  art::Handle< std::vector<ldp::TOF> > TOFColHandle;
  e.getByLabel(fTOFModuleLabel,TOFColHandle);

  std::cout << "@@@@@@@@@Selected Possible Kaon Run/Subrun/Event: " << e.run() << "/" << e.subRun() << "/" << e.event() << std::endl;
  

  // #########################################################
  // ## If there is no ParticleID object then return false ###
  // #########################################################
  if(!particleIDCol->size()) return false;
  
  //Finding best-guess Particles
  double pdg_temp = 0;
  pdg_temp=fParticlePDG;

  std::cout<<"Particle PDG "<<pdg_temp<<"\n";
  
  // ################################################################
  // ### Filling the TOF vs WC Track Momentum Histo prior to cuts ###
  // ################################################################ 
  fPzVsTOF->Fill(WCTrackColHandle->at(0).Momentum(),TOFColHandle->at(0).SingleTOF(0));

  
  // ###################################
  // ### Identifying Kaon Candidates ###
  // ###################################
  if(pdg_temp == 321 || pdg_temp == -321){
    if( particleIDCol->at(0).PDGCode() == 321 || particleIDCol->at(0).PDGCode() == -321 ){
      if( particleIDCol->at(0).KaonProbability() > fParticleProbCutOff ){
	std::cout << "Selected Possible Kaon Run/Subrun/Event: " << e.run() << "/" << e.subRun() << "/" << e.event() << std::endl;
	fParticlePzVsTOF->Fill(WCTrackColHandle->at(0).Momentum(),TOFColHandle->at(0).SingleTOF(0));
	return true;
      }
    }
    return false;
  } 
  
  // ###################################
  // ### Identifying Pion Candidates ###
  // ###################################
  else if(pdg_temp == 211 || pdg_temp == -211){
    if( particleIDCol->at(0).PDGCode() == 211 || particleIDCol->at(0).PDGCode() == -211 ){
      if( particleIDCol->at(0).PionProbability() > fParticleProbCutOff ){
	std::cout << "Selected Possible Pion Run/Subrun/Event: " << e.run() << "/" << e.subRun() << "/" << e.event() << std::endl;
	fParticlePzVsTOF->Fill(WCTrackColHandle->at(0).Momentum(),TOFColHandle->at(0).SingleTOF(0));
	return true;
      }
    }
    return false;
  }
  
  // ###################################
  // ### Identifying Muon Candidates ###
  // ###################################
  else if(pdg_temp == 13 || pdg_temp == -13){
    if( particleIDCol->at(0).PDGCode() == 13 || particleIDCol->at(0).PDGCode() == -13 ){
      if( particleIDCol->at(0).MuonProbability() > fParticleProbCutOff ){
	std::cout << "Selected Possible Muon Run/Subrun/Event: " << e.run() << "/" << e.subRun() << "/" << e.event() << std::endl;
	fParticlePzVsTOF->Fill(WCTrackColHandle->at(0).Momentum(),TOFColHandle->at(0).SingleTOF(0));
	return true;
      }
    }
    return false;
  }
  
  // ########################################
  // ### Identifying Muon/Pion Candidates ### (Note: This is a degenerate case)
  // ########################################
  else if(pdg_temp == 21113){

    
    auto WC = *WCTrackColHandle;
    auto TOF = *TOFColHandle;
  
    if (TOF.size() < 1) return false;
    
    // ### LOOP OVER THE WCTRACKS AND TOF OBJECTS ###
    
    
    if(WC[0].Momentum() > 100 && WC[0].Momentum() < 1500 && 
       TOF[0].SingleTOF(0) > 10 && TOF[0].SingleTOF(0) < 25)
      { return true;}
	
    	
	
    
    //PDG -> Pi: 211, Mu: 13, so PiMu is 21113
    //Use PiMu pdg flag due to lack of good MuRS Tracks (as it is now)  
    std::cout<< "###################################################### "<<std::endl ;
    std::cout<< "Looking At Every Event..... Run/Subrun/Event: " << e.run() << "/" << e.subRun() << "/" << e.event() << std::endl;
    for (unsigned int i = 0; i<  particleIDCol->size();i++)  
      {
	std::cout<<i<<" Particle 2113? "<<particleIDCol->at(i).PDGCode() <<"\n";
	std::cout<<i<<" Particle p  ? "<<particleIDCol->at(i).ProtonProbability() <<" cutOff "<<fParticleProbCutOff <<"\n";
	std::cout<<i<<" Particle mu ? "<<particleIDCol->at(i).MuonProbability() <<" cutOff "<<fParticleProbCutOff <<"\n";
	std::cout<<i<<" Particle pi ? "<<particleIDCol->at(i).PionProbability() <<" cutOff "<<fParticleProbCutOff <<"\n";
	std::cout<<i<<" Particle pimu? "<<particleIDCol->at(i).PiMuProbability() <<" cutOff "<<fParticleProbCutOff <<"\n";
	std::cout<<i<<" Particle K? "<<particleIDCol->at(i).KaonProbability() <<" cutOff "<<fParticleProbCutOff <<"\n";

      }
    std::cout<< "###################################################### "<<(*particleIDCol).size()<<std::endl ;
    std::cout<< "###################################################### "<<particleIDCol->size()<<std::endl ;

    if( particleIDCol->at(0).PDGCode() == 21113 ){ 
      if( particleIDCol->at(0).PiMuProbability() > fParticleProbCutOff ){
	if(!WCTrackColHandle->size()) return false;
	if(!TOFColHandle->size()) return false;
	std::cout << "Selected Possible Pi/Mu Run/Subrun/Event: " << e.run() << "/" << e.subRun() << "/" << e.event() << std::endl;
	fParticlePzVsTOF->Fill(WCTrackColHandle->at(0).Momentum(),TOFColHandle->at(0).SingleTOF(0));
	std::cout<< "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% "<<std::endl ;
	return true;
      }
    }
    std::cout << "@@@@@@@@@8\n";
  

  return false;
  }
  
  // #####################################
  // ### Identifying Proton Candidates ###
  // #####################################
  else if(pdg_temp == 2212 || pdg_temp == -2212){
    if( particleIDCol->at(0).PDGCode() == 2212 || particleIDCol->at(0).PDGCode() == -2212 ){
      if( particleIDCol->at(0).ProtonProbability() > fParticleProbCutOff ){
	std::cout << "Selected Possible Proton Run/Subrun/Event: " << e.run() << "/" << e.subRun() << "/" << e.event() << std::endl;
	fParticlePzVsTOF->Fill(WCTrackColHandle->at(0).Momentum(),TOFColHandle->at(0).SingleTOF(0));
	std::cout << "@@@@@@@@@9\n";
	return true;
      }
    }
    std::cout << "@@@@@@@@@10\n";
    return false;
  }
  
  // ###########################################
  // ###    If you don't satisfy any of      ###
  // ### these, then the event does not pass ###
  // ###########################################
  else return false; 
  
  
}



// ##########################
// ### Begin Job Function ###
// ##########################
void ParticleFilter::beginJob()
{
  // #######################################
  // ### Defining histograms to be saved ###
  // #######################################
  art::ServiceHandle<art::TFileService> tfs;
  fPzVsTOF = tfs->make<TH2F>("PzVsTOF","Pz Vs. TOF (All) ",160,0,1600,70,10,80);
  fParticlePzVsTOF = tfs->make<TH2F>("ParticlePzVsTOF","Particle Pz Vs. TOF",160,0,1600,70,10,80);  //that's for the selected particles

}

void ParticleFilter::endJob()
{
  // Implementation of optional member function here.
}



void ParticleFilter::respondToCloseInputFile(art::FileBlock const & fb)
{
  // Implementation of optional member function here.
}

void ParticleFilter::respondToCloseOutputFiles(art::FileBlock const & fb)
{
  // Implementation of optional member function here.
}

void ParticleFilter::respondToOpenInputFile(art::FileBlock const  &fb)
{
  // Implementation of optional member function here.
}

void ParticleFilter::respondToOpenOutputFiles(art::FileBlock const & fb)
{
  // Implementation of optional member function here.
}

DEFINE_ART_MODULE(ParticleFilter)
