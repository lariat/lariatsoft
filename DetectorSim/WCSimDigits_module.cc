////////////////////////////////////////////////////////////////////////
// Class:       WCSimDigits
// Module Type: producer
// File:        WCSimDigits_module.cc
//
// Generated at Wed Sep 28 11:16:11 2016 by Greg Pulliam using artmod
// from cetpkgsupport v1_10_02.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "lardataobj/RawData/AuxDetDigit.h"
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/AuxDetGeo.h"
#include "larcore/Geometry/AuxDetGeometry.h"
#include "lardataobj/Simulation/sim.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "lardataobj/Simulation/AuxDetSimChannel.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include <memory>
#include "TMath.h"
#include "TComplex.h"
#include "TString.h"
#include "TH2.h"
#include "TH1F.h"
#include "TFile.h"
#include "TRandom.h"
#include "TTree.h"

class WCSimDigits;

class WCSimDigits : public art::EDProducer {
public:
  explicit WCSimDigits(fhicl::ParameterSet const & p);
  // The destructor generated by the compiler is fine for classes
  // without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  WCSimDigits(WCSimDigits const &) = delete;
  WCSimDigits(WCSimDigits &&) = delete;
  WCSimDigits & operator = (WCSimDigits const &) = delete;
  WCSimDigits & operator = (WCSimDigits &&) = delete;

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
  std::vector<raw::AuxDetDigit> MakeWCDigits(sim::AuxDetIDE & TheIDE, int  AuxDetID);
  void ResetVars();
  int neutrals[9]={22,14,2112,310,130,311,12,16,111}; //Neutrally charged particles that shouldn't make digits
  bool IsNeutral;
  TTree* fTree;
  int WCx1[100], WCx2[100], WCx3[100], WCx4[100], WCy1[100], WCy2[100], WCy3[100], WCy4[100]; //Array with wire hit for each Digit. Hopefully no WC is hit by >100 particles.
  float wc1pz[100],wc2pz[100],wc3pz[100],wc4pz[100],wc1theta[100],wc2theta[100],wc3theta[100],wc4theta[100]; //Arrays for P_z and Theta_xz for each digit.
  int WC1PDG[100], WC2PDG[100], WC3PDG[100], WC4PDG[100]; // PDG of particle passing through WC
  float wc1t[100],wc2t[100],wc3t[100],wc4t[100]; //time of digit made in WC
  int WC1Count, WC2Count,WC3Count,WC4Count, WC1PDGCount, WC2PDGCount, WC3PDGCount, WC4PDGCount; //A counter for how often a WC is hit.
  bool pickytrack, highyield, WC1Hit, WC2Hit, WC3Hit, WC4Hit; //Some bools if WCs have any hits, and if the digits made allow for a picky track or a high yield track.
private:

 std::string fG4ModuleLabel; 
 std::string fPickyTracks;

};


WCSimDigits::WCSimDigits(fhicl::ParameterSet const & p)
// :
// Initialize member data here.
{
  this->reconfigure(p);
  produces<std::vector<raw::AuxDetDigit> >();
}

void WCSimDigits::produce(art::Event & e)
{
  std::unique_ptr<std::vector<raw::AuxDetDigit> > MWPCDigits(new std::vector<raw::AuxDetDigit> );
  
   ResetVars(); //Reset tree variables  
   
    //Get the ADSCs
  art::Handle< std::vector<sim::AuxDetSimChannel> > AuxDetHandle;
  std::vector<sim::AuxDetSimChannel> AuxDetColl;
  
  //Get MCParticle info
  art::Handle<std::vector<simb::MCParticle> > g4_part;
  e.getByLabel(fG4ModuleLabel, g4_part);
  int numG4=(int)g4_part->size();
   
  //if(!e.getByLabel(fG4ModuleLabel,AuxDetHandle)){e.put(std::move(DigitIDEAssn)); return;} //If no sim information found, put an empty assn and return.
  e.getByLabel(fG4ModuleLabel, AuxDetHandle);
  
  //Loop over the ADSC
  for(std::vector<sim::AuxDetSimChannel>::const_iterator auxiter = AuxDetHandle->begin(); auxiter!=AuxDetHandle->end(); ++auxiter){
    const sim::AuxDetSimChannel & aux = *auxiter;
    
     
    //Only want WCs, which have AuxDet of 1-4
    if(aux.AuxDetID()>0 && aux.AuxDetID()<5){
      std::vector<sim::AuxDetIDE> SimIDE=aux.AuxDetIDEs();
      
     
    
      //Loop over IDEs, which are individual particles passsing through the AuxDet
      for(size_t nIDE=0; nIDE<SimIDE.size(); ++nIDE){
        sim::AuxDetIDE TheIDE=SimIDE[nIDE];
	// Assuming trigger happens around t=0, as it would be if the event occured exactly as in G4BL. WC timing only reads out for ~300 bins before to ~724 bins after the trigger. 
	// Only attempt to make digits if the time would be in that window.
	if(TheIDE.exitT<-300*1.177 || TheIDE.exitT>724*1.177) continue; 
	int IDETrackID=TheIDE.trackID;
	IsNeutral=false;
	
	
	//Loop over the MCParticles for this event, and try to match an IDE TrackID to an MCParticle
	for(int i=0; i<numG4; ++i){
	  if((double)g4_part->at(i).TrackId()==IDETrackID){
	    int G4PDG=g4_part->at(i).PdgCode();
	    
	    
	    //Now loop through neutrals and skip this IDE if it was made by a neutral.
	    for(int nNeu=0; nNeu<8; ++nNeu){
	      if(fabs(G4PDG)==neutrals[nNeu]){IsNeutral=true; break;}
	    }//Neutral loop
	    
	    //If we dont have a neutral, make a WCDigit
	    if(!IsNeutral){
	       std::cout<<"Making Digit with PDG: "<<G4PDG<<std::endl;
	       int AuxDetId=aux.AuxDetID();
	       std::cout<<"Using "<<AuxDetId<<std::endl;
	       if(AuxDetId==1){
	         WC1Hit=true;
		 WC1PDG[WC1PDGCount]=G4PDG;
		 wc1pz[WC1PDGCount]=TheIDE.exitMomentumZ;
		 wc1theta[WC1PDGCount]=atan(TheIDE.exitMomentumX/TheIDE.exitMomentumZ);
		 wc1t[WC1PDGCount]=TheIDE.exitT;
		 ++WC1PDGCount;
	       }
	       if(AuxDetId==2){
	         WC2Hit=true;
		 WC2PDG[WC2PDGCount]=G4PDG;
		 wc2pz[WC2PDGCount]=TheIDE.exitMomentumZ;
		 wc2theta[WC2PDGCount]=atan(TheIDE.exitMomentumX/TheIDE.exitMomentumZ);
		 wc2t[WC2PDGCount]=TheIDE.exitT;		 
		 ++WC2PDGCount;		 
	       }
	       if(AuxDetId==3){
	         WC3Hit=true;
		 WC3PDG[WC3PDGCount]=G4PDG;
		 wc3pz[WC3PDGCount]=TheIDE.exitMomentumZ;
		 wc3theta[WC3PDGCount]=atan(TheIDE.exitMomentumX/TheIDE.exitMomentumZ);
		 wc3t[WC3PDGCount]=TheIDE.exitT;		 		 
		 ++WC3PDGCount;		 
	       }
	       if(AuxDetId==4){
	         WC4Hit=true;
		 WC4PDG[WC4PDGCount]=G4PDG;
		 wc4pz[WC4PDGCount]=TheIDE.exitMomentumZ;
		 wc4theta[WC4PDGCount]=atan(TheIDE.exitMomentumX/TheIDE.exitMomentumZ);	
		 wc4t[WC4PDGCount]=TheIDE.exitT;	 
		 ++WC4PDGCount;		 
	       }
	       std::vector<raw::AuxDetDigit> SimWCDigits=MakeWCDigits(TheIDE, AuxDetId);
	       for(size_t Digititer=0; Digititer<SimWCDigits.size(); ++Digititer){
	         (*MWPCDigits).push_back(SimWCDigits[Digititer]);
	       }//Digititer
	    }//IsNeutral
	  }//If matched particle	  
	}// numG4 loop
      } //nIDE loop
    }  //If ADSC is a WC
  }//ADSC loop
  if(WC1Hit && WC4Hit && (WC2Hit || WC3Hit)){highyield=true;}
  if(WC1Hit && WC2Hit && WC3Hit && WC4Hit && (*MWPCDigits).size()==8){pickytrack=true;}
  std::cout<<"Putting "<<(*MWPCDigits).size()<<" digits on the event"<<std::endl;
  //Put the Digits on the event.
  e.put(std::move(MWPCDigits));
  fTree->Fill();
}//produce


std::vector<raw::AuxDetDigit> WCSimDigits::MakeWCDigits(sim::AuxDetIDE & TheIDE, int AuxDetID)
{


  
  std::vector<raw::AuxDetDigit> MWPCDigits;
  //Use AuxDetGeo to get the center of the WC, given the AuxDetID
  art::ServiceHandle<geo::Geometry> adGeoServ;
  geo::AuxDetGeo const& adg=adGeoServ->AuxDet(AuxDetID);
  double centerarray[3]={0.0,0.0,0.0};
  adg.GetCenter(centerarray,0);
  
  //Each Digit needs a string for the WC the digit is for. Make a vector of the strings, and we can get them when we need them.
  std::vector<std::string> WCnames;
  WCnames.resize(4);
  WCnames[0]="MWPC1";
  WCnames[1]="MWPC2";
  WCnames[2]="MWPC3";
  WCnames[3]="MWPC4";
  
  //Start with the XY position of the hit in TPC coordinate frame
  double xy_tpc[2]={.5*(TheIDE.exitX+TheIDE.entryX),.5*(TheIDE.exitY+TheIDE.entryY)}; 
  
  //Move the origin to the center of the WC
  double xy_wcCent[2]={xy_tpc[0]-centerarray[0],xy_tpc[1]-centerarray[1]}; 
  
  //The WCs are rotated at 3 and 13 degrees. Therefore, the z position of the hit is tan(theta)*x. Make tan(theta) array 
  double rotarray[4]={tan(13*3.14159265359/180),tan(13*3.14159265359/180),tan(3*3.14159265359/180),tan(3*3.14159265359/180)};
  
  //Using the distance formula d^2=x^2+z^2, we can find how far along the WC the hit is.
  double tempx_dist=pow(pow(xy_wcCent[0],2)+pow(rotarray[AuxDetID-1]*xy_wcCent[0],2),.5); 
  
  //There is an ambiguity whether tempx_dist is +ive or -ive. Using the hit in the WC frame we can determine that and add/subtract half a WC length
  
  if(xy_wcCent[0]<0){tempx_dist=tempx_dist+6.4;} // The hit would be on the right side of the WC, so add on half the size of the WC to make up for the fact that the origin is halfway along the x axis
  if(xy_wcCent[0]>0){tempx_dist=6.4-tempx_dist;} //Left side of the WC, so subtract off half length.  Makes the extreme edge of the WC x=0
  
  
  double tempy_dist=6.4-xy_wcCent[1]; //No distance formula necessary, just moving the origin in y.
  
  //The WC are 12.8cm by 12.8cm. If the hit is not inside the wc, don't make digits.
  
  if(tempx_dist>0 && tempx_dist<12.8 && tempy_dist>0 && tempy_dist<12.8){
  //Now we have the hit in the WC frame, with the center at the top left corner (From beam perspective).  Also, convert from mm to cm.  Wires are 1mm apart, so working in mm makes sense.
    double xy_wc[2]={tempx_dist*10, tempy_dist*10};
  
  //Initialize the x wire and y wire for each particle, and the timebin of the hit.
    int xwire, ywire;
    double timebin;
    
  //The WC have 1mm spaced wires, and I assume the first wire is .5mm from the edge of the frame. The wires that gets hit will be the closest wires. Hits in the edge .5mm gaps get associated to the first or last wire
    if(xy_wc[0]<=.5) {xwire=0;}
    else if(xy_wc[0]>=127.5) {xwire=127;}
    else {xwire=round(xy_wc[0]-.5);}
  
    if(xy_wc[1]<=.5){ywire=0;}
    else if(xy_wc[1]>=127.5){ywire=127;}
    else {ywire=round(xy_wc[1]-.5);}
  
    timebin=TheIDE.exitT/1.177; //WCs readout in 1.117 ns bins
    if(timebin<0){std::cout<<"Time check "<<TheIDE.exitT<<std::endl;}
  //We have all the information necessary for a digit. Now we create the 4 variables that go into a digit. 
  //To know the digit is made from simulation, I'm going to store the TrackID as the timestamp.  
  
    unsigned long long TimestampfromTrackID=TheIDE.trackID;
    std::vector<short> timebinvector;
    timebinvector.push_back(timebin);
  
  //Each IDE makes 2 digits, one each for the x and y wire

  
    MWPCDigits.push_back(raw::AuxDetDigit(static_cast<unsigned short> (xwire),
	 				     timebinvector,
					     WCnames[AuxDetID-1],
					     static_cast <unsigned long long> (TimestampfromTrackID)));
  
    MWPCDigits.push_back(raw::AuxDetDigit(static_cast<unsigned short> (ywire+128),
	 				     timebinvector,
					     WCnames[AuxDetID-1],
					     static_cast <unsigned long long> (TimestampfromTrackID)));
  if(AuxDetID==1)
  {
    WCx1[WC1Count]=xwire;
    WCy1[WC1Count]=ywire;
    ++WC1Count;
  }
    if(AuxDetID==2)
  {
    WCx2[WC2Count]=xwire;
    WCy2[WC2Count]=ywire;
    ++WC2Count;
  }
    if(AuxDetID==3)
  {
    WCx3[WC3Count]=xwire;
    WCy3[WC3Count]=ywire;
    ++WC3Count;
  }
    if(AuxDetID==4)
  {
    WCx4[WC4Count]=xwire;
    WCy4[WC4Count]=ywire;
    ++WC4Count;
  }

  }
  //Fill the Tree with the wire hits for this IDE
  					     
  return MWPCDigits;
					     
}
void WCSimDigits::beginJob()
{
  art::ServiceHandle<art::TFileService> tfs;
  fTree=tfs->make<TTree>("WCDigitTree", "WCDigitTree");
  fTree->Branch("WCx1",&WCx1,"WCx1[100]/I");
  fTree->Branch("WCx2",&WCx2,"WCx2[100]/I");
  fTree->Branch("WCx3",&WCx3,"WCx3[100]/I");
  fTree->Branch("WCx4",&WCx4,"WCx4[100]/I");
  fTree->Branch("WCy1",&WCy1,"WCy1[100]/I");
  fTree->Branch("WCy2",&WCy2,"WCy2[100]/I");
  fTree->Branch("WCy3",&WCy3,"WCy3[100]/I");
  fTree->Branch("WCy4",&WCy4,"WCy4[100]/I");
  fTree->Branch("pickytrack", &pickytrack, "pickytrack/O"); 
  fTree->Branch("highyield", &highyield, "highyield/O");
  fTree->Branch("WC1Hit", &WC1Hit,"WC1Hit/O");
  fTree->Branch("WC2Hit", &WC2Hit,"WC2Hit/O");
  fTree->Branch("WC3Hit", &WC3Hit,"WC3Hit/O");
  fTree->Branch("WC4Hit", &WC4Hit,"WC4Hit/O");
  fTree->Branch("WC1PDG", &WC1PDG,"WC1PDG[100]/I");
  fTree->Branch("WC2PDG", &WC2PDG,"WC2PDG[100]/I");
  fTree->Branch("WC3PDG", &WC3PDG,"WC3PDG[100]/I");
  fTree->Branch("WC4PDG", &WC4PDG,"WC4PDG[100]/I");
  fTree->Branch("wc1pz", &wc1pz, "wc1pz[100]/F");
  fTree->Branch("wc1theta", &wc1theta, "wc1theta[100]/F");
  fTree->Branch("wc2pz", &wc2pz, "wc2pz[100]/F");
  fTree->Branch("wc2theta", &wc2theta, "wc2theta[100]/F");
  fTree->Branch("wc3pz", &wc3pz, "wc3pz[100]/F");
  fTree->Branch("wc3theta", &wc3theta, "wc3theta[100]/F");
  fTree->Branch("wc4pz", &wc4pz, "wc4pz[100]/F");
  fTree->Branch("wc4theta", &wc4theta, "wc4theta[100]/F"); 
  fTree->Branch("wc1t", &wc1t, "wc1t[100]/F");
  fTree->Branch("wc2t", &wc2t, "wc2t[100]/F");     
  fTree->Branch("wc3t", &wc3t, "wc3t[100]/F");     
  fTree->Branch("wc4t", &wc4t, "wc4t[100]/F");          
}
void WCSimDigits::ResetVars()
{
  for(int i=0; i<100; ++i)
  {
    WCx1[i]=-99999;
    WCx2[i]=-99999;
    WCx3[i]=-99999;
    WCx4[i]=-99999;
    WCy1[i]=-99999;
    WCy2[i]=-99999;
    WCy3[i]=-99999;
    WCy4[i]=-99999;
    WC1PDG[i]=-99999;
    WC2PDG[i]=-99999;
    WC3PDG[i]=-99999;
    WC4PDG[i]=-99999;
    wc1pz[i]=-99999;
    wc1theta[i]=-99999;
    wc2pz[i]=-99999;
    wc2theta[i]=-99999;
    wc3pz[i]=-99999;
    wc3theta[i]=-99999;
    wc4pz[i]=-99999;
    wc4theta[i]=-99999; 
    wc1t[i]=-99999;
    wc2t[i]=-99999; 
    wc3t[i]=-99999; 
    wc4t[i]=-99999;            
  }
  highyield=false;
  pickytrack=false;
  WC1Hit=false;
  WC2Hit=false;
  WC3Hit=false;
  WC4Hit=false;
  WC1Count=0;
  WC2Count=0;
  WC3Count=0;
  WC4Count=0;
  WC1PDGCount=0;
  WC2PDGCount=0;
  WC3PDGCount=0;
  WC4PDGCount=0;
 
  
}
void WCSimDigits::beginRun(art::Run & r)
{
  
}

void WCSimDigits::beginSubRun(art::SubRun & sr)
{
  // Implementation of optional member function here.
}

void WCSimDigits::endJob()
{
  // Implementation of optional member function here.
}

void WCSimDigits::endRun(art::Run & r)
{
  // Implementation of optional member function here.
}

void WCSimDigits::endSubRun(art::SubRun & sr)
{
  // Implementation of optional member function here.
}

void WCSimDigits::reconfigure(fhicl::ParameterSet const & p)
{
  fG4ModuleLabel = p.get<std::string>("G4ModuleLabel");
  fPickyTracks= p.get<std::string>("PickyTracks");
  std::cout<<"Picky Tracks is "<<fPickyTracks<<std::endl;
}

void WCSimDigits::respondToCloseInputFile(art::FileBlock const & fb)
{
  // Implementation of optional member function here.
}

void WCSimDigits::respondToCloseOutputFiles(art::FileBlock const & fb)
{
  // Implementation of optional member function here.
}

void WCSimDigits::respondToOpenInputFile(art::FileBlock const & fb)
{
  // Implementation of optional member function here.
}

void WCSimDigits::respondToOpenOutputFiles(art::FileBlock const & fb)
{
  // Implementation of optional member function here.
}

DEFINE_ART_MODULE(WCSimDigits)
