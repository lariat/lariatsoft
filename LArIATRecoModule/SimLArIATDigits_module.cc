////////////////////////////////////////////////////////////////////////
// Class:       SimLArIATDigits
// Module Type: producer
// File:        SimLArIATDigits_module.cc
//
// Generated at Tue May 31 13:53:13 2016 by Greg Pulliam using artmod
// from cetpkgsupport v1_08_06.
////////////////////////////////////////////////////////////////////////
#include <vector>
#include <string>
#include <algorithm>
#include <sstream>
#include <fstream>
#include <bitset>
#include <iostream>

extern "C" {
#include <sys/types.h>
#include <sys/stat.h>
}

#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "lardata/Utilities/LArFFT.h"
#include "lardata/RawData/RawDigit.h"
#include "lardata/RawData/raw.h"
#include "lardata/RawData/TriggerData.h"
#include "lardata/DetectorInfoServices/LArPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/DetectorInfoServices/DetectorClocksServiceStandard.h" // special (see below)
#include "Utilities/SignalShapingServiceT1034.h"
#include "larcore/Geometry/Geometry.h"
#include "larsim/Simulation/sim.h"
#include "larsim/Simulation/SimChannel.h"
#include "larsim/Simulation/AuxDetSimChannel.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"

#include "TMath.h"
#include "TComplex.h"
#include "TString.h"
#include "TH2.h"
#include "TH1D.h"
#include "TFile.h"
#include "TRandom.h"
#include "TTree.h"
#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/RandGaussQ.h"

const int kMaxDet=50;
const int kMaxIDE=1000;

class SimLArIATDigits;

class SimLArIATDigits : public art::EDProducer {
public:
  explicit SimLArIATDigits(fhicl::ParameterSet const & p);
  // The destructor generated by the compiler is fine for classes
  // without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  SimLArIATDigits(SimLArIATDigits const &) = delete;
  SimLArIATDigits(SimLArIATDigits &&) = delete;
  SimLArIATDigits & operator = (SimLArIATDigits const &) = delete;
  SimLArIATDigits & operator = (SimLArIATDigits &&) = delete;

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

private:
 // Declare member data here.
  //std::vector<raw::AuxDetDigit> const& MakeWCDigits(art::Event & evt);  
  std::string fG4ModuleLabel; 
  void ResetVars(); //Function to reset all values for next run
  int numSimChannels;  //Number of Aux Dets 
  int numIDEs[kMaxDet];  //Number of IDEs per Aux Det
 //The ID number for a particular Aux Det
  //Storing positon and momentum of Geant particle//
  
  
  double entryx[kMaxDet][kMaxIDE];
  double entryy[kMaxDet][kMaxIDE];
  double entryz[kMaxDet][kMaxIDE];
  double exitx[kMaxDet][kMaxIDE];
  double exity[kMaxDet][kMaxIDE];
  double exitz[kMaxDet][kMaxIDE];
  double exitmomx[kMaxDet][kMaxIDE];
  double exitmomy[kMaxDet][kMaxIDE];
  double exitmomz[kMaxDet][kMaxIDE];
  double AuxDetID[kMaxDet];
  int iterarray[kMaxDet];
  double TOFangle[kMaxIDE];
  double TrackID[kMaxDet][kMaxIDE]; 
  double Energy[kMaxDet][kMaxIDE]; 
  //double enterx;
  //double entery;
  //double enterz;
  int ID; //The AuxDetID() for an AuxDetSimChannel
  TTree* fTree; 
  TH2D* XZHit;
};


SimLArIATDigits::SimLArIATDigits(fhicl::ParameterSet const & p)
// :
// Initialize member data here.
{
  // Call appropriate produces<>() functions here.
  this->reconfigure(p);
}

void SimLArIATDigits::produce(art::Event & e)
{
  ResetVars();
  
  art::Handle< std::vector<sim::AuxDetSimChannel> > AuxDetHandle;
  //std::vector<sim::AuxDetSimChannel*> AuxDetCollection;
  e.getByLabel(fG4ModuleLabel, AuxDetHandle);
  numSimChannels=AuxDetHandle->size();
  std::cout<<"numSimChannel: "<<numSimChannels<<std::endl;
  int iter=0;
  // for(size_t i=0; i<numSimChannels; ++i){
   for(std::vector<sim::AuxDetSimChannel>::const_iterator auxiter = AuxDetHandle->begin(); auxiter!=AuxDetHandle->end(); ++auxiter){
     const sim::AuxDetSimChannel & aux = *auxiter;
     AuxDetID[iter]=aux.AuxDetID();
     ID=aux.AuxDetID();
     iterarray[iter]=iter;
     std::vector<sim::AuxDetIDE> SimIDE=aux.AuxDetIDEs();
     numIDEs[iter]=SimIDE.size();  
     std::cout<<"For Sim Channel: "<<iter<<", there are "<<SimIDE.size()<<" IDEs. AuxDetID: "<<aux.AuxDetID()<<std::endl;
     for(size_t nIDE=0; nIDE<SimIDE.size(); ++nIDE){
       sim::AuxDetIDE TheIDE=SimIDE[nIDE];
       entryx[iter][nIDE]=TheIDE.entryX;
       entryy[iter][nIDE]=TheIDE.entryY;
       entryz[iter][nIDE]=TheIDE.entryZ;
       exitx[iter][nIDE]=TheIDE.exitX;
       exity[iter][nIDE]=TheIDE.exitY;
       exitz[iter][nIDE]=TheIDE.exitZ;
       TrackID[iter][nIDE]=TheIDE.trackID; 
       exitmomx[iter][nIDE]=TheIDE.exitMomentumX;
       exitmomy[iter][nIDE]=TheIDE.exitMomentumY;
       exitmomz[iter][nIDE]=TheIDE.exitMomentumZ;
       Energy[iter][nIDE]=TheIDE.energyDeposited;
       if(iter==0){
         TOFangle[nIDE]=180/(3.141593)*tan(TheIDE.exitMomentumX/TheIDE.exitMomentumZ);
       } 
       std::cout<<"ID: "<<ID<<" nIDE: "<<nIDE<<" Total IDEs: "<<SimIDE.size()<<std::endl;
       //enterx=TheIDE.entryX;
       //entery=TheIDE.entryY;
       //enterz=TheIDE.entryZ;
       XZHit->Fill(TheIDE.entryZ,TheIDE.entryX); 
     }
     ++iter;
   }
fTree->Fill();  
}

void SimLArIATDigits::beginJob()
{
  // Implementation of optional member function here.
  art::ServiceHandle<art::TFileService> tfs;
  fTree=tfs->make<TTree>("GeantTree","Geant Tree");
  fTree->Branch("numSimChannels",&numSimChannels,"numSimChannels/I");
  fTree->Branch("numIDEs",numIDEs,"numIDEs[numSimChannels]/I");
  fTree->Branch("entryx",entryx,"entryx[numSimChannels][1000]/D");
  fTree->Branch("entryy",entryy,"entryy[numSimChannels][1000]/D");
  fTree->Branch("entryz",entryz,"entryz[numSimChannels][1000]/D");
  fTree->Branch("exitx",exitx,"exitx[numSimChannels][1000]/D");
  fTree->Branch("exity",exity,"exity[numSimChannels][1000]/D");
  fTree->Branch("exitz",exitz,"exitz[numSimChannels][1000]/D");
  fTree->Branch("exitmomx",exitmomx,"exitmomx[numSimChannels][1000]/D");
  fTree->Branch("exitmomy",exitmomy,"exitmomy[numSimChannels][1000]/D");
  fTree->Branch("exitmomz",exitmomz,"exitmomz[numSimChannels][1000]/D");
  fTree->Branch("TrackID",TrackID,"TrackID[numSimChannels][1000]/D");
  fTree->Branch("Energy",Energy,"Energy[numSimChannels][1000]/D");
  //fTree->Branch("enterx",enterx,"enterx/D");
  //fTree->Branch("entery",entery,"entery/D");
  //fTree->Branch("enterz",enterz,"enterz/D");
  fTree->Branch("AuxDetID",AuxDetID,"AuxDetID[numSimChannels]/D");
  fTree->Branch("iterarray",iterarray,"iter[numSimChannels]/I");
  fTree->Branch("TOFangle",TOFangle,"TOFangle[numIDEs]/D");
  XZHit = tfs->make<TH2D>("XZHit", "XZHit", 1500,-1000,500,1500,-150,150);  
  
  
}
void SimLArIATDigits::ResetVars()
{
  numSimChannels=-9999;
  for(int i=0; i<kMaxDet; ++i){
    AuxDetID[i]=-9999;
    iterarray[i]=-9999;
    numIDEs[i]=-9999;
    for(int j=0; j<kMaxIDE; ++j){
      entryx[i][j]=-9999;
      entryy[i][j]=-9999;
      entryz[i][j]=-9999;
      exitx[i][j]=-9999;
      exity[i][j]=-9999;
      exitz[i][j]=-9999;
      exitmomx[i][j]=-9999;
      exitmomy[i][j]=-9999;
      exitmomz[i][j]=-9999;
      TrackID[i][j]=-9999;
      TOFangle[j]=-9999;
      Energy[i][j]=-9999;
    }
  }
}
void SimLArIATDigits::beginRun(art::Run & r)
{
  // Implementation of optional member function here.
}

void SimLArIATDigits::beginSubRun(art::SubRun & sr)
{
  // Implementation of optional member function here.
}

void SimLArIATDigits::endJob()
{
  // Implementation of optional member function here.
}

void SimLArIATDigits::endRun(art::Run & r)
{
  // Implementation of optional member function here.
}

void SimLArIATDigits::endSubRun(art::SubRun & sr)
{
  // Implementation of optional member function here.
}

void SimLArIATDigits::reconfigure(fhicl::ParameterSet const & p)
{
  // Implementation of optional member function here.
  fG4ModuleLabel = p.get<std::string>("G4ModuleLabel");
}

void SimLArIATDigits::respondToCloseInputFile(art::FileBlock const & fb)
{
  // Implementation of optional member function here.
}

void SimLArIATDigits::respondToCloseOutputFiles(art::FileBlock const & fb)
{
  // Implementation of optional member function here.
}

void SimLArIATDigits::respondToOpenInputFile(art::FileBlock const & fb)
{
  // Implementation of optional member function here.
}

void SimLArIATDigits::respondToOpenOutputFiles(art::FileBlock const & fb)
{
  // Implementation of optional member function here.
}

DEFINE_ART_MODULE(SimLArIATDigits)
