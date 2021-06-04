////////////////////////////////////////////////////////////////////////
// Class:       TOFSimDigits
// Module Type: producer
// File:        TOFSimDigits_module.cc
//
// Generated at Wed Jul 20 14:54:59 2016 by Lucas Mendes Santos using artmod
// from cetpkgsupport v1_10_01.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art_root_io/TFileService.h"
#include "art_root_io/TFileDirectory.h"

//LArSoft libraries

#include "lardata/Utilities/LArFFT.h"
#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RawData/AuxDetDigit.h"
#include "lardataobj/RawData/raw.h"
#include "lardataobj/RawData/TriggerData.h"
#include "RawDataUtilities/FragmentToDigitAlg.h"
#include "lardata/DetectorInfoServices/LArPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/DetectorInfoServices/DetectorClocksServiceStandard.h" // special (see below)
#include "Utilities/SignalShapingServiceT1034.h"
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/AuxDetGeo.h"
#include "larcore/Geometry/AuxDetGeometry.h"
#include "lardataobj/Simulation/sim.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "lardataobj/Simulation/AuxDetSimChannel.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"

#include "TF1.h"
#include "TTree.h"

#include <memory>

class TOFSimDigits;

class TOFSimDigits : public art::EDProducer {
public:
  explicit TOFSimDigits(fhicl::ParameterSet const & p);
  // The destructor generated by the compiler is fine for classes
  // without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  TOFSimDigits(TOFSimDigits const &) = delete;
  TOFSimDigits(TOFSimDigits &&) = delete;
  TOFSimDigits & operator = (TOFSimDigits const &) = delete;
  TOFSimDigits & operator = (TOFSimDigits &&) = delete;
  void  reconfigure(fhicl::ParameterSet const & p) ;
  void beginJob() override;
  // Required functions.
  void produce(art::Event & e) override;
  std::string fG4ModuleLabel;

private:
 
  TTree *fTree;
  int toft0, toftf;

  double usentryT, usexitT;
  double dsentryT, dsexitT;

  // Declare member data here.

};

void TOFSimDigits::reconfigure(fhicl::ParameterSet const & p)
{
  // Implementation of optional member function here.
  
  fG4ModuleLabel = p.get<std::string>("G4ModuleLabel");
  
}

void TOFSimDigits::beginJob()
{
  art::ServiceHandle<art::TFileService> tfs;
  fTree = tfs->make<TTree>("tofana","Sim TOF");
  fTree->Branch("usentryT", &usentryT,"usentryT/D");
  fTree->Branch("usexitT", &usexitT,"usexitT/D");
  fTree->Branch("dsentryT", &usentryT,"usentryT/D");
  fTree->Branch("dseexitT", &usexitT,"usexitT/D");
  fTree->Branch("simt0", &toft0, "simt0/I");
  fTree->Branch("simtf", &toftf, "simtf/I");
}


TOFSimDigits::TOFSimDigits(fhicl::ParameterSet const & p)
: EDProducer(p)
// Initialize member data here.
{
  // Call appropriate produces<>() functions here.
  this->reconfigure(p);
  produces< std::vector<raw::AuxDetDigit>   >();
}

void TOFSimDigits::produce(art::Event & e)
{
  // Implementation of required member function here.

  art::Handle< std::vector<sim::AuxDetSimChannel> > AuxDetHandle;
  //std::vector<sim::AuxDetSimChannel*> AuxDetCollection;
  e.getByLabel(fG4ModuleLabel, AuxDetHandle);
  std::unique_ptr< std::vector<raw::AuxDetDigit>> tofdigits(new std::vector<raw::AuxDetDigit>);
  int numSimChannels=AuxDetHandle->size();
  std::cout<<"numSimChannel: "<<numSimChannels<<std::endl;
  //int iter=0;
  int ID;

  //double energy;
  //std::string detName;
  //short channel;

  //Generating a example waveform
  std::vector <short> fadc;
  TF1 *wvsim = new TF1("wvsim", "[0] -[1]*TMath::Landau(x,[2],[3])",0,28672);
  
  std::vector <short> fadc2;
  TF1 *wvsim2 = new TF1("wvsim2", "[0] -[1]*TMath::Landau(x,[2],[3])",0,28672);

  std::vector <short> fadc3;
  std::vector <short> fadc4;

  double t0;
  double tf;



   if(AuxDetHandle->size() > 0){
   std::cout<<AuxDetHandle->size()<<std::endl;   
   for(std::vector<sim::AuxDetSimChannel>::const_iterator auxiter = AuxDetHandle->begin(); auxiter!=AuxDetHandle->end(); ++auxiter){
     const sim::AuxDetSimChannel & aux = *auxiter;
    
     ID=aux.AuxDetID();
     //art::ServiceHandle<geo::Geometry> adGeoServ;
     
     std::vector<sim::AuxDetIDE> SimIDE=aux.AuxDetIDEs();
     if(ID == 0 || ID == 6){
     switch(ID)
     {
       case 0:
	//detName = "TOFUS";
	//channel = 5;
	
	if(SimIDE.size()>0){ 
  	usentryT = (double)SimIDE.at(0).entryT;
       	usexitT = (double)SimIDE.at(0).exitT;   
 	t0 = (usentryT + usexitT)/2;

        toft0 = (int)t0;
	std::cout<<usentryT<<" "<<usexitT<<" TOFUS time: "<<t0<<std::endl;
        wvsim->SetParameters(929.25, 1323.58, t0 + 8800, 2.10069);
        wvsim2->SetParameters(924.43, 996.51, t0 + 8800, 2.11423);
       
        for(int i = 0; i < 28672; i++){fadc.push_back((short)wvsim->Eval(i));  fadc2.push_back((short)wvsim2->Eval(i));}    		
        tofdigits->push_back(raw::AuxDetDigit(5, fadc,  "TOFUS",0));
        tofdigits->push_back(raw::AuxDetDigit(5, fadc2, "TOFUS",0));

        }
	else toft0 = 0;
        
        
        break;
   
      case 6:
      	if(SimIDE.size() > 0){
        dsentryT = (double)SimIDE.at(0).entryT;
       	dsexitT = (double)SimIDE.at(0).exitT;
	tf = (dsexitT + dsentryT)/2;
	std::cout<<dsentryT<<" "<<dsexitT<<" TOFDS time: "<<tf<<std::endl;
        wvsim->SetParameters(919.446, 2502.53, tf + 8800, 2.54303);
        wvsim2->SetParameters(920.706, 2264.09, tf + 8800, 2.22591);
        for(int i = 0; i < 28672; i++){fadc3.push_back((short)wvsim->Eval(i));
        fadc4.push_back((short)wvsim2->Eval(i));}
        tofdigits->push_back(raw::AuxDetDigit(6, fadc3, "TOFDS",0));
        tofdigits->push_back(raw::AuxDetDigit(6, fadc4, "TOFDS",0));}
        else{ tf = 0.;}
        std::cout<<tf<<" "<<t0<<std::endl;

        break;

    
    }//switch
   }//if ID
   
   }//for
   }//if AuxDetDigit->size()
  if(tofdigits->size()>3){

     if(toftf - toft0 > 0){
     
     //
     short fHitThreshold;
     int j = 0;
	    while(j < (int)fadc3.size())
	    {
		fHitThreshold = fadc3[j+1] - fadc3[j];
		if(fHitThreshold < -40){ toftf = j; break;}	
		else j++;
	    }
     j = 0;
	   while(j < (int)fadc.size())
	    {
		fHitThreshold = fadc[j+1] - fadc[j];
		if(fHitThreshold < -40){ toft0 = j; break;}	
		else j++;
	    }

      fTree->Fill();
  }  



  }
  e.put(std::move(tofdigits));
}

DEFINE_ART_MODULE(TOFSimDigits)
