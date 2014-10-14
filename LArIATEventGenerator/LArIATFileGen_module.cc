////////////////////////////////////////////////////////////////////////
/// \file  LArIATFileGen.cxx
/// \brief Generator for muons from a file.
///
/// Module designed to produce a set list of particles for a MC event
///
/// \version $Id: LArIATFileGen.cxx,v 1.4 2010/03/29 09:54:01 echurch Exp $
/// \author  echurch@fnal.gov
////////////////////////////////////////////////////////////////////////
// C++ includes.
#include <string>
#include <cmath>
#include <memory>

// Framework includes
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// nutools includes
#include "SimulationBase/MCTruth.h"
#include "SimulationBase/MCParticle.h"

// lar includes
#include "Geometry/Geometry.h"
#include "SummaryData/RunData.h"
#include "Utilities/DetectorProperties.h"

#include "TVector3.h"
#include "TDatabasePDG.h"

#include "art/Framework/Core/ModuleMacros.h"


#include <vector>
#include <string>
#include "art/Framework/Core/EDProducer.h"
#include "TFile.h"
#include "TTree.h"
#include "TTreeIndex.h"
#include "TEntryList.h"

#include <stdio.h>
#include <iostream>
#include <fstream>

namespace simb { class MCTruth; }

namespace evgen {

    
  class EventFrame{
    public:
    EventFrame(std::vector< int > pdg, std::vector< TLorentzVector > pos, std::vector< TLorentzVector > mom)
    {
     fPdg=pdg; 
     fPos=pos;
     fMom=mom;
    };
    
    ~EventFrame()
    {
      fPdg.clear();
      fPos.clear();
      fMom.clear();
    };
    
    std::vector< int > fPdg; 
    std::vector< TLorentzVector > fPos;
    std::vector< TLorentzVector > fMom;
    
    
    
  };
  
  
  /// module to produce single or multiple specified particles in the detector
  class LArIATFileGen : public art::EDProducer {

  public:
    explicit LArIATFileGen(fhicl::ParameterSet const& pset);
    virtual ~LArIATFileGen();

    // This is called for each event.
    void produce(art::Event& evt);
    void beginJob();
    void beginRun(art::Run& run);
    void endJob();

  private:

    void ReadEvents(simb::MCTruth &mct);        

    
    
    
    int                 fSeed;           // random number seed    
    std::vector<int>    fPDG;           
    std::vector<double> fXYZ_Off;           
    
    
    
    
    
    
    
    std::string fFileName;
    std::string fMuonsFileType;
    std::string fTreeName;
    std::vector<std::string> fBranchNames; 
    int fEventSpillOffset;  // Where in file to start.
    int fEventsPerSpill;
    int fEventsPerFile; 
    int fEventFileOffset; 
    std::string fTriggerCut;
    bool fUseTrigger; 
    
    
    
    ifstream *fMuonFile;
    TFile *fMuonFileR;
    TTree *TNtuple;
    TTree *TNtupleSmall;
    unsigned int countFile;

    Float_t xtmp, ytmp, ztmp;
    Float_t pxtmp, pytmp, pztmp;
    Float_t charge;
    Float_t         E;
    Float_t         costheta;
    Float_t         phi;
    Float_t         xdet;
    Float_t         ydet;
    Float_t         zdet;

    
    
    int EventID,TrackID;
    TBranch * b_eventid, * b_trackid;
    TEntryList *elist; //used for triggered events
    
    
    //Det1
    double xDet1,yDet1,zDet1,tDet1,PxDet1,PyDet1,PzDet1,PDGidDet1,ParentIDDet1;bool TrackPresentDet1;
    
    TBranch *b_xDet1,*b_yDet1,*b_zDet1,*b_tDet1,*b_PxDet1,*b_PyDet1,*b_PzDet1,*b_PDGidDet1,*b_ParentIDDet1,*b_TrackPresentDet1;
    //Det2
    double xDet2,yDet2,zDet2,tDet2,PxDet2,PyDet2,PzDet2,PDGidDet2,ParentIDDet2;bool TrackPresentDet2;
    TBranch *b_xDet2,*b_yDet2,*b_zDet2,*b_tDet2,*b_PxDet2,*b_PyDet2,*b_PzDet2,*b_PDGidDet2,*b_ParentIDDet2,*b_TrackPresentDet2;
    //Det3
    double xDet3,yDet3,zDet3,tDet3,PxDet3,PyDet3,PzDet3,PDGidDet3,ParentIDDet3;bool TrackPresentDet3;
    TBranch *b_xDet3,*b_yDet3,*b_zDet3,*b_tDet3,*b_PxDet3,*b_PyDet3,*b_PzDet3,*b_PDGidDet3,*b_ParentIDDet3,*b_TrackPresentDet3;
    //Det4
    double xDet4,yDet4,zDet4,tDet4,PxDet4,PyDet4,PzDet4,PDGidDet4,ParentIDDet4;bool TrackPresentDet4;
    TBranch *b_xDet4,*b_yDet4,*b_zDet4,*b_tDet4,*b_PxDet4,*b_PyDet4,*b_PzDet4,*b_PDGidDet4,*b_ParentIDDet4,*b_TrackPresentDet4;
    //TOFus
    double xTOFus,yTOFus,zTOFus,tTOFus,PxTOFus,PyTOFus,PzTOFus,PDGidTOFus,ParentIDTOFus;bool TrackPresentTOFus;
    TBranch *b_xTOFus,*b_yTOFus,*b_zTOFus,*b_tTOFus,*b_PxTOFus,*b_PyTOFus,*b_PzTOFus,*b_PDGidTOFus,*b_ParentIDTOFus,*b_TrackPresentTOFus;
    //TODdsVert
    double xTODdsVert,yTODdsVert,zTODdsVert,tTODdsVert,PxTODdsVert,PyTODdsVert,PzTODdsVert,PDGidTODdsVert,ParentIDTODdsVert;bool TrackPresentTODdsVert;
    TBranch *b_xTODdsVert,*b_yTODdsVert,*b_zTODdsVert,*b_tTODdsVert,*b_PxTODdsVert,*b_PyTODdsVert,*b_PzTODdsVert,*b_PDGidTODdsVert,*b_ParentIDTODdsVert,*b_TrackPresentTODdsVert;
    //TODdsHorz
    double xTODdsHorz,yTODdsHorz,zTODdsHorz,tTODdsHorz,PxTODdsHorz,PyTODdsHorz,PzTODdsHorz,PDGidTODdsHorz,ParentIDTODdsHorz;bool TrackPresentTODdsHorz;
    TBranch *b_xTODdsHorz,*b_yTODdsHorz,*b_zTODdsHorz,*b_tTODdsHorz,*b_PxTODdsHorz,*b_PyTODdsHorz,*b_PzTODdsHorz,*b_PDGidTODdsHorz,*b_ParentIDTODdsHorz,*b_TrackPresentTODdsHorz;
    //Halo
    double xHalo,yHalo,zHalo,tHalo,PxHalo,PyHalo,PzHalo,PDGidHalo,ParentIDHalo;bool TrackPresentHalo;
    TBranch *b_xHalo,*b_yHalo,*b_zHalo,*b_tHalo,*b_PxHalo,*b_PyHalo,*b_PzHalo,*b_PDGidHalo,*b_ParentIDHalo,*b_TrackPresentHalo;
    //HaloHole
    double xHaloHole,yHaloHole,zHaloHole,tHaloHole,PxHaloHole,PyHaloHole,PzHaloHole,PDGidHaloHole,ParentIDHaloHole;bool TrackPresentHaloHole;
    TBranch *b_xHaloHole,*b_yHaloHole,*b_zHaloHole,*b_tHaloHole,*b_PxHaloHole,*b_PyHaloHole,*b_PzHaloHole,*b_PDGidHaloHole,*b_ParentIDHaloHole,*b_TrackPresentHaloHole;
    //TiWindow
    double xTiWindow,yTiWindow,zTiWindow,tTiWindow,PxTiWindow,PyTiWindow,PzTiWindow,PDGidTiWindow,ParentIDTiWindow;bool TrackPresentTiWindow;
    TBranch *b_xTiWindow,*b_yTiWindow,*b_zTiWindow,*b_tTiWindow,*b_PxTiWindow,*b_PyTiWindow,*b_PzTiWindow,*b_PDGidTiWindow,*b_ParentIDTiWindow,*b_TrackPresentTiWindow;
    
    std::vector<  EventFrame  > fEventFrames;
    int fEventCounter;

  };
  
  
  

  
  
  
} // namespace



namespace evgen{

  //____________________________________________________________________________
  LArIATFileGen::LArIATFileGen(fhicl::ParameterSet const& pset)
     : fSeed(314159)					                          
     , fFileName         (pset.get< std::string     	     	>("FileName")         	)      
     , fMuonsFileType    (pset.get< std::string             	>("MuonsFileType")    	)
     , fTreeName         (pset.get< std::string   	     	>("TreeName")         	)      
     , fBranchNames      (pset.get< std::vector<std::string> 	>("BranchNames")      	)
     , fEventSpillOffset (pset.get<int				>("EventSpillOffset")	)
     , fEventsPerSpill (pset.get<int				>("EventsPerSpill")	)
     , fEventsPerFile (pset.get<int				>("EventsPerFile")	)
     , fEventFileOffset (pset.get<int				>("EventFileOffset")	) 
     , fTriggerCut         (pset.get< std::string     	     	>("TriggerCut")       	)  
     , fUseTrigger(pset.get<     bool                		>("UseTrigger")		)
  {

    produces< std::vector<simb::MCTruth> >();
    produces< sumdata::RunData, art::InRun >();

  }

  //____________________________________________________________________________
  void LArIATFileGen::beginJob()
  {
    countFile=fEventSpillOffset;
    mf::LogInfo("LArIATFileGen : starting at event ")  << countFile*fEventsPerSpill <<std::endl;
    
  
    if (fMuonsFileType.compare("root")==0) 
      {
	std::cout << "LArIATFileGen: You have chosen to read muons from Root File " << fFileName << std::endl; 
      }
    else 
      {
	std::cout << "LArIATFileGen: You must specify one of source/text/root file to read for muons."<< std::endl; 
        return ; 
      }

      
    fEventCounter=0;
  //Create subtree of the     
      
    fMuonFileR = new TFile(fFileName.c_str(),"READ");
    TNtuple = (TTree*)(fMuonFileR->Get(fTreeName.c_str()));

    //Set Branches for all detectors to extract the data per spill.
    
    TNtuple->SetBranchAddress("EventID", &EventID, &b_eventid);
    TNtuple->SetBranchAddress("TrackID", &TrackID, &b_trackid);
    
    //Det1
    TNtuple->SetBranchAddress("xDet1", &xDet1, &b_xDet1);
    TNtuple->SetBranchAddress("yDet1", &yDet1, &b_yDet1);
    TNtuple->SetBranchAddress("zDet1", &zDet1, &b_zDet1);
    TNtuple->SetBranchAddress("tDet1", &tDet1, &b_tDet1);
    TNtuple->SetBranchAddress("PxDet1", &PxDet1, &b_PxDet1);
    TNtuple->SetBranchAddress("PyDet1", &PyDet1, &b_PyDet1);
    TNtuple->SetBranchAddress("PzDet1", &PzDet1, &b_PzDet1);
    TNtuple->SetBranchAddress("PDGidDet1", &PDGidDet1, &b_PDGidDet1);
    TNtuple->SetBranchAddress("ParentIDDet1", &ParentIDDet1, &b_ParentIDDet1);
    TNtuple->SetBranchAddress("TrackPresentDet1", &TrackPresentDet1, &b_TrackPresentDet1);
    
    //Det2
    TNtuple->SetBranchAddress("xDet2", &xDet2, &b_xDet2);
    TNtuple->SetBranchAddress("yDet2", &yDet2, &b_yDet2);
    TNtuple->SetBranchAddress("zDet2", &zDet2, &b_zDet2);
    TNtuple->SetBranchAddress("tDet2", &tDet2, &b_tDet2);
    TNtuple->SetBranchAddress("PxDet2", &PxDet2, &b_PxDet2);
    TNtuple->SetBranchAddress("PyDet2", &PyDet2, &b_PyDet2);
    TNtuple->SetBranchAddress("PzDet2", &PzDet2, &b_PzDet2);
    TNtuple->SetBranchAddress("PDGidDet2", &PDGidDet2, &b_PDGidDet2);
    TNtuple->SetBranchAddress("ParentIDDet2", &ParentIDDet2, &b_ParentIDDet2);
    TNtuple->SetBranchAddress("TrackPresentDet2", &TrackPresentDet2, &b_TrackPresentDet2);
    
    //Det3
    TNtuple->SetBranchAddress("xDet3", &xDet3, &b_xDet3);
    TNtuple->SetBranchAddress("yDet3", &yDet3, &b_yDet3);
    TNtuple->SetBranchAddress("zDet3", &zDet3, &b_zDet3);
    TNtuple->SetBranchAddress("tDet3", &tDet3, &b_tDet3);
    TNtuple->SetBranchAddress("PxDet3", &PxDet3, &b_PxDet3);
    TNtuple->SetBranchAddress("PyDet3", &PyDet3, &b_PyDet3);
    TNtuple->SetBranchAddress("PzDet3", &PzDet3, &b_PzDet3);
    TNtuple->SetBranchAddress("PDGidDet3", &PDGidDet3, &b_PDGidDet3);
    TNtuple->SetBranchAddress("ParentIDDet3", &ParentIDDet3, &b_ParentIDDet3);
    TNtuple->SetBranchAddress("TrackPresentDet3", &TrackPresentDet3, &b_TrackPresentDet3);
    
    //Det4
    TNtuple->SetBranchAddress("xDet4", &xDet4, &b_xDet4);
    TNtuple->SetBranchAddress("yDet4", &yDet4, &b_yDet4);
    TNtuple->SetBranchAddress("zDet4", &zDet4, &b_zDet4);
    TNtuple->SetBranchAddress("tDet4", &tDet4, &b_tDet4);
    TNtuple->SetBranchAddress("PxDet4", &PxDet4, &b_PxDet4);
    TNtuple->SetBranchAddress("PyDet4", &PyDet4, &b_PyDet4);
    TNtuple->SetBranchAddress("PzDet4", &PzDet4, &b_PzDet4);
    TNtuple->SetBranchAddress("PDGidDet4", &PDGidDet4, &b_PDGidDet4);
    TNtuple->SetBranchAddress("ParentIDDet4", &ParentIDDet4, &b_ParentIDDet4);
    TNtuple->SetBranchAddress("TrackPresentDet4", &TrackPresentDet4, &b_TrackPresentDet4);
    
    //TOFus
    TNtuple->SetBranchAddress("xTOFus", &xTOFus, &b_xTOFus);
    TNtuple->SetBranchAddress("yTOFus", &yTOFus, &b_yTOFus);
    TNtuple->SetBranchAddress("zTOFus", &zTOFus, &b_zTOFus);
    TNtuple->SetBranchAddress("tTOFus", &tTOFus, &b_tTOFus);
    TNtuple->SetBranchAddress("PxTOFus", &PxTOFus, &b_PxTOFus);
    TNtuple->SetBranchAddress("PyTOFus", &PyTOFus, &b_PyTOFus);
    TNtuple->SetBranchAddress("PzTOFus", &PzTOFus, &b_PzTOFus);
    TNtuple->SetBranchAddress("PDGidTOFus", &PDGidTOFus, &b_PDGidTOFus);
    TNtuple->SetBranchAddress("ParentIDTOFus", &ParentIDTOFus, &b_ParentIDTOFus);
    TNtuple->SetBranchAddress("TrackPresentTOFus", &TrackPresentTOFus, &b_TrackPresentTOFus);
    
    //TODdsVert
    TNtuple->SetBranchAddress("xTODdsVert", &xTODdsVert, &b_xTODdsVert);
    TNtuple->SetBranchAddress("yTODdsVert", &yTODdsVert, &b_yTODdsVert);
    TNtuple->SetBranchAddress("zTODdsVert", &zTODdsVert, &b_zTODdsVert);
    TNtuple->SetBranchAddress("tTODdsVert", &tTODdsVert, &b_tTODdsVert);
    TNtuple->SetBranchAddress("PxTODdsVert", &PxTODdsVert, &b_PxTODdsVert);
    TNtuple->SetBranchAddress("PyTODdsVert", &PyTODdsVert, &b_PyTODdsVert);
    TNtuple->SetBranchAddress("PzTODdsVert", &PzTODdsVert, &b_PzTODdsVert);
    TNtuple->SetBranchAddress("PDGidTODdsVert", &PDGidTODdsVert, &b_PDGidTODdsVert);
    TNtuple->SetBranchAddress("ParentIDTODdsVert", &ParentIDTODdsVert, &b_ParentIDTODdsVert);
    TNtuple->SetBranchAddress("TrackPresentTODdsVert", &TrackPresentTODdsVert, &b_TrackPresentTODdsVert);
    
    //TODdsHorz
    TNtuple->SetBranchAddress("xTODdsHorz", &xTODdsHorz, &b_xTODdsHorz);
    TNtuple->SetBranchAddress("yTODdsHorz", &yTODdsHorz, &b_yTODdsHorz);
    TNtuple->SetBranchAddress("zTODdsHorz", &zTODdsHorz, &b_zTODdsHorz);
    TNtuple->SetBranchAddress("tTODdsHorz", &tTODdsHorz, &b_tTODdsHorz);
    TNtuple->SetBranchAddress("PxTODdsHorz", &PxTODdsHorz, &b_PxTODdsHorz);
    TNtuple->SetBranchAddress("PyTODdsHorz", &PyTODdsHorz, &b_PyTODdsHorz);
    TNtuple->SetBranchAddress("PzTODdsHorz", &PzTODdsHorz, &b_PzTODdsHorz);
    TNtuple->SetBranchAddress("PDGidTODdsHorz", &PDGidTODdsHorz, &b_PDGidTODdsHorz);
    TNtuple->SetBranchAddress("ParentIDTODdsHorz", &ParentIDTODdsHorz, &b_ParentIDTODdsHorz);
    TNtuple->SetBranchAddress("TrackPresentTODdsHorz", &TrackPresentTODdsHorz, &b_TrackPresentTODdsHorz);
    
    //Halo
    TNtuple->SetBranchAddress("xHalo", &xHalo, &b_xHalo);
    TNtuple->SetBranchAddress("yHalo", &yHalo, &b_yHalo);
    TNtuple->SetBranchAddress("zHalo", &zHalo, &b_zHalo);
    TNtuple->SetBranchAddress("tHalo", &tHalo, &b_tHalo);
    TNtuple->SetBranchAddress("PxHalo", &PxHalo, &b_PxHalo);
    TNtuple->SetBranchAddress("PyHalo", &PyHalo, &b_PyHalo);
    TNtuple->SetBranchAddress("PzHalo", &PzHalo, &b_PzHalo);
    TNtuple->SetBranchAddress("PDGidHalo", &PDGidHalo, &b_PDGidHalo);
    TNtuple->SetBranchAddress("ParentIDHalo", &ParentIDHalo, &b_ParentIDHalo);
    TNtuple->SetBranchAddress("TrackPresentHalo", &TrackPresentHalo, &b_TrackPresentHalo);
    
    //HaloHole
    TNtuple->SetBranchAddress("xHaloHole", &xHaloHole, &b_xHaloHole);
    TNtuple->SetBranchAddress("yHaloHole", &yHaloHole, &b_yHaloHole);
    TNtuple->SetBranchAddress("zHaloHole", &zHaloHole, &b_zHaloHole);
    TNtuple->SetBranchAddress("tHaloHole", &tHaloHole, &b_tHaloHole);
    TNtuple->SetBranchAddress("PxHaloHole", &PxHaloHole, &b_PxHaloHole);
    TNtuple->SetBranchAddress("PyHaloHole", &PyHaloHole, &b_PyHaloHole);
    TNtuple->SetBranchAddress("PzHaloHole", &PzHaloHole, &b_PzHaloHole);
    TNtuple->SetBranchAddress("PDGidHaloHole", &PDGidHaloHole, &b_PDGidHaloHole);
    TNtuple->SetBranchAddress("ParentIDHaloHole", &ParentIDHaloHole, &b_ParentIDHaloHole);
    TNtuple->SetBranchAddress("TrackPresentHaloHole", &TrackPresentHaloHole, &b_TrackPresentHaloHole);
    
    //TiWindow
    TNtuple->SetBranchAddress("xTiWindow", &xTiWindow, &b_xTiWindow);
    TNtuple->SetBranchAddress("yTiWindow", &yTiWindow, &b_yTiWindow);
    TNtuple->SetBranchAddress("zTiWindow", &zTiWindow, &b_zTiWindow);
    TNtuple->SetBranchAddress("tTiWindow", &tTiWindow, &b_tTiWindow);
    TNtuple->SetBranchAddress("PxTiWindow", &PxTiWindow, &b_PxTiWindow);
    TNtuple->SetBranchAddress("PyTiWindow", &PyTiWindow, &b_PyTiWindow);
    TNtuple->SetBranchAddress("PzTiWindow", &PzTiWindow, &b_PzTiWindow);
    TNtuple->SetBranchAddress("PDGidTiWindow", &PDGidTiWindow, &b_PDGidTiWindow);
    TNtuple->SetBranchAddress("ParentIDTiWindow", &ParentIDTiWindow, &b_ParentIDTiWindow);
    TNtuple->SetBranchAddress("TrackPresentTiWindow", &TrackPresentTiWindow, &b_TrackPresentTiWindow);
 
    
    TFile * testfile=new TFile("out.root","RECREATE");
    TNtupleSmall = TNtuple->CloneTree(0);    
    Long64_t nentries = TNtuple->GetEntries();
    art::ServiceHandle<geo::Geometry> geom;
    //Create small TTree for sorting.
    // check that neiter of these numbers goes out of bounds.
    //nentries = (nentries > fEventNumberOffset+fEventsPerSpill ) ? fEventNumberOffset+fEventsPerSpill : nentries;
    std::cout << " creating a subtree for sorting - this make take a while. from evt: " <<  fEventSpillOffset*fEventsPerSpill << " to: " << nentries << std::endl;
    
  //   int counter =0;
    long start_event=fEventSpillOffset*fEventsPerSpill;
    for (Long64_t i=0;i<nentries; i++) {
      TNtuple->GetEntry(i);
    //  if(counter++ < 20)
//	std::cout << " EventID " << EventID << " " <<  TrackID << " " << tHalo << std::endl;
      if (EventID > start_event   && EventID < start_event+fEventsPerSpill) TNtupleSmall->Fill();
      if(EventID > start_event+fEventsPerSpill)
	break;
   }
    
    TNtupleSmall->Write();
    
    std::cout << "old tree size, new tree size: " << TNtuple->GetEntries() << " " << TNtupleSmall->GetEntries() << std::endl;
    
    art::ServiceHandle<util::DetectorProperties> detprop;

    
    const double frametime=detprop->SamplingRate()*detprop->NumberTimeSamples()/1000.; 
    int frameindex=1;
    
    std::vector< int > locPDG;
    std::vector< TLorentzVector > locXYZ;
    std::vector< TLorentzVector > locMOM;
    
    
    TDatabasePDG  pdgt;
    double m = 0.;
    std::cout << " building index on time on the subtree - this make also take a while " << std::endl;
    TNtupleSmall->BuildIndex("tHalo*1000000");  //sorting on integers, need microseconds
    TTreeIndex *index = (TTreeIndex*)TNtupleSmall->GetTreeIndex();
   
  // if(fUseTrigger)
  // {
    TNtupleSmall->Draw(">>elist", fTriggerCut.data(), "entrylist");
    elist = (TEntryList*)gDirectory->Get("elist");
    Long64_t ListNEntries=elist->GetN();
    if(fUseTrigger)
      std::cout << ListNEntries << " entries in Event List " << std::endl;
    double last_tHalo_trig=0;
    
    //if running with trigger then loop only on selected entries, i.e. ListNEntries
    // ir running without trigger then loop on all events and save them by frame.
    Long64_t MainIndexEntries = (fUseTrigger) ? ListNEntries : index->GetN() - 1;
    //loop over selected trigger entries.
    for (int i = 0;  i<MainIndexEntries ; i++)
      {
	int tree_entry; 
	
	if(fUseTrigger)
	  tree_entry=elist->GetEntry(i);
	else
	  tree_entry=TNtupleSmall->LoadTree( index->GetIndex()[i] );
	  
	TNtupleSmall->GetEntry(tree_entry);
	
	// std::cout <<  "  Selected Trigger event, i " << i << " " << EventID << " " << TrackID << " present in TrackPresentTOFus" << TrackPresentTOFus << " tHalo " <<  tHalo*1000 << " tree_entry " << tree_entry << " curr trig time and distance: "<< last_tHalo_trig*1000 << " "<< (tHalo- last_tHalo_trig)*1000 << " frtime " << frametime/1000000 << " " << last_tHalo_trig +frametime/1000000<<  std::endl;
	
 	if(tHalo<last_tHalo_trig+frametime/1000000) // cannot retrigger on same frame.
 	  continue;
	
	if(fEventFrames.size() >  (unsigned int)fEventsPerFile)  // break if too many events.
	  break;
	
	if(fUseTrigger)
	  last_tHalo_trig=tHalo;
	else
	  {
	   while(tHalo > frameindex*frametime/1000000)  
	      {frameindex++; }
	   last_tHalo_trig=frameindex*frametime/1000000;
	  }
	//std::cout << "cut: " << Form("tHalo > %lf && tHalo < %lf",last_tHalo_trig-frametime/1000000,last_tHalo_trig+frametime/1000000) << std::endl;
	TNtupleSmall->Draw(">>loclist",Form("tHalo > %lf && tHalo < %lf",last_tHalo_trig-frametime/1000000,last_tHalo_trig+frametime/1000000),"entrylist");
	
	TEntryList * loclist = (TEntryList*)gDirectory->Get("loclist");
	 Long64_t ListlocNEntries=loclist->GetN();
	
	// std::cout << " found" << ListlocNEntries << " other particles in frame " << std::endl;
	 
	locPDG.clear();
	locXYZ.clear();
	locMOM.clear(); 
	 
	for (int xx = 0;  xx<ListlocNEntries ; xx++)
	  {int tree_entry=loclist->GetEntry(xx);
	  TNtupleSmall->GetEntry(tree_entry);
	//  std::cout <<  "  Selected Frame event " << EventID << " " << TrackID << " present in TrackPresentTOFus" << TrackPresentTOFus << " tHalo " << tHalo << " tree_entry " << tree_entry << std::endl;
	
	 double x,y,Px,Py,Pz,PDG ;//z,
      
	if(TrackPresentHaloHole)
	  {
	  x=xHaloHole;y=yHaloHole;//z=zHaloHole;
	  Px=PxHaloHole;Py=PyHaloHole;Pz=PzHaloHole;	
	  PDG=PDGidHaloHole;
	  }
	else
	  {
	  x=xHalo;y=yHalo;//z=zHalo;
	  Px=PxHalo;Py=PyHalo;Pz=PzHalo;	
	  PDG=PDGidHalo;
	  }
      
        
	TParticlePDG* pdgp = pdgt.GetParticle(PDG);
	if (pdgp) m = pdgp->Mass();
	
	double Ener=std::sqrt((Px*Px+Py*Py+Pz*Pz)/1000./1000.+m*m);

	//      titatnium window position wrt. TPCActive (based on gdml file wdisk)
      // in gdml	    (-0,435 ,0.2 ,-84.54646643)
      // to transfer into LArSOFT coordinates, X changes by HalfWidth
      // Z changes by HalfLength. 
      // geom->CryostatHalfHeight()*0.01,geom->CryostatLength()
         
	locPDG.push_back(PDGidHalo);
 	locXYZ.push_back(TLorentzVector(x/10.-0.435+geom->DetHalfWidth(),y/10.+0.2,-84.54646+geom->DetLength()/2.-1,(tHalo-last_tHalo_trig)*1e9)); //convert to ns
	locMOM.push_back(TLorentzVector(Px/1000.,Py/1000.,Pz/1000.,Ener));
	  
	//first back track to find particles that are one frame before:
	//while(tHalo > tHalo_trig-frametime/1000000)
	 // {
	  } //end for loop inside of frame

	if(locPDG.size() && locXYZ.size() && locMOM.size() )  // skipping non needed files
	  {  fEventFrames.push_back(EventFrame(locPDG,locXYZ,locMOM)); 
	  }
	  
      } //end for loop on triggered particles
	
//    }
//    else{ 
//     std::cout << " looping on subtree. " << std::endl;
//     for(Long64_t i = 0; i<index->GetN() - 1; ++i ) {
//       Long64_t local = TNtupleSmall->LoadTree( index->GetIndex()[i] );
//       TNtupleSmall->GetEntry(local);
//    
//       
//      if(fEventFrames.size() >  fEventsPerFile)
// 	break;
//       
//      
//        if(tHalo > frameindex*frametime/1000000)   // want in seconds not us
//        {
// 	if(locPDG.size() && locXYZ.size() && locMOM.size() && frameindex>fEventFileOffset)  // skipping non needed files
// 	  {  fEventFrames.push_back(EventFrame(locPDG,locXYZ,locMOM));   }
// 	while(tHalo > frameindex*frametime/1000000)  
// 	  {frameindex++; }
// 	locPDG.clear();
// 	locXYZ.clear();
// 	locMOM.clear();
//        } // end if tHalo > frameIndex
//      
//      
//      //pick up particle info
//       
//       double x,y,z,Px,Py,Pz,PDG ;
//       
//       if(TrackPresentHaloHole)
// 	{
// 	x=xHaloHole;y=yHaloHole;z=zHaloHole;
// 	Px=PxHaloHole;Py=PyHaloHole;Pz=PzHaloHole;	
// 	PDG=PDGidHaloHole;
// 	}
//       else
// 	{
// 	x=xHalo;y=yHalo;z=zHalo;
// 	Px=PxHalo;Py=PyHalo;Pz=PzHalo;	
// 	PDG=PDGidHalo;
// 	}
//       
//           
//       TParticlePDG* pdgp = pdgt.GetParticle(PDG);
//       if (pdgp) m = pdgp->Mass();
//       double Ener=std::sqrt((Px*Px+Py*Py+Pz*Pz)/1000./1000.+m*m);
//     
//     //      titatnium window position wrt. TPCActive (based on gdml file wdisk)
//     // in gdml	    (-0,435 ,0.2 ,-84.54646643)
//     // to transfer into LArSOFT coordinates, X changes by HalfWidth
//      // Z changes by HalfLength. 
//      // geom->CryostatHalfHeight()*0.01,geom->CryostatLength()
//          
//       locPDG.push_back(PDGidHalo);
//       locXYZ.push_back(TLorentzVector(x/10.-0.435+geom->DetHalfWidth(),y/10.+0.2,-84.54646+geom->DetLength()/2.-1,(tHalo*1000000-(frameindex-1)*frametime)*1e3));
//       locMOM.push_back(TLorentzVector(Px/1000.,Py/1000.,Pz/1000.,Ener));
//     	
//       
//      
//       }
//      } // if !useTrigger
     
     if (fMuonsFileType.compare("root")==0) 
      fMuonFileR->Close();
      testfile->Close();
    
  }

  //____________________________________________________________________________
  void LArIATFileGen::endJob()
  {
  //  if (fMuonsFileType.compare("root")==0) 
  //    fMuonFileR->Close();
  }

  //____________________________________________________________________________
  void LArIATFileGen::beginRun(art::Run& run)
  {

    // grab the geometry object to see what geometry we are using
    art::ServiceHandle<geo::Geometry> geo;

    std::unique_ptr<sumdata::RunData> runcol(new sumdata::RunData(geo->DetectorName()));

    run.put(std::move(runcol));

    return;
  }

  //____________________________________________________________________________
  LArIATFileGen::~LArIATFileGen()
  {
  }

  //____________________________________________________________________________
  void LArIATFileGen::produce(art::Event& evt)
  {

   // check that the file is still good
    if((unsigned int)fEventCounter>=fEventFrames.size())
    { throw cet::exception("LArIATFileGen") << "ran out of events in fEventFrames. Exiting "; }
    
    
    ///unique_ptr allows ownership to be transferred to the art::Event after the put statement
    std::unique_ptr< std::vector<simb::MCTruth> > truthcol(new std::vector<simb::MCTruth>);

    simb::MCTruth truth;
    truth.SetOrigin(simb::kSingleParticle);

      
    for(unsigned int ix=0;ix<fEventFrames[fEventCounter].fPdg.size();ix++)
      {
	
      std::cout << " adding particle:"  << fEventFrames[fEventCounter].fPdg[ix] << " "<< fEventFrames[fEventCounter].fPos[ix].X() << " "<< fEventFrames[fEventCounter].fPos[ix].Y() <<" "<< fEventFrames[fEventCounter].fPos[ix].Z() <<" time "<< fEventFrames[fEventCounter].fPos[ix].T()  << " energy: " << fEventFrames[fEventCounter].fMom[ix].T() << std::endl;
	

      int trackid = -1*(ix+1); // set track id to -i as these are all primary particles and have id <= 0
      std::string primary("primary");
      simb::MCParticle part(trackid, fEventFrames[fEventCounter].fPdg[ix], primary);
      part.AddTrajectoryPoint(fEventFrames[fEventCounter].fPos[ix], fEventFrames[fEventCounter].fMom[ix]);
      
      //       std::cout << "add the particle to the primary" << std::endl;
      
      truth.Add(part);
    

      }
//     std::cout << "put mctruth into the vector" << std::endl;
    truthcol->push_back(truth);

//     std::cout << "add vector to the event " << truthcol->size() << std::endl;
    evt.put(std::move(truthcol));
    fEventCounter++;
    return;
  }

  //____________________________________________________________________________
  void LArIATFileGen::ReadEvents(simb::MCTruth &mct) 
  {

//     std::cout << "size of particle vector is " << fPDG.size() << std::endl;

    ///every event will have one of each particle species in the fPDG array
    for (unsigned int i=0; i<fPDG.size(); ++i) {
      
      // Choose momentum
      //double p = 0.0;
      double m(0.108);

      TVector3 x;
      TVector3 p;
      Double_t q = 0.;
      Int_t pdgLocal;
      
      if (fMuonsFileType.compare("text")==0)
	{

	  std::string line;
	  getline(*fMuonFile,line);
	  if (!fMuonFile->good())
	    {
	      std::cout << "LArIATFileGen: Problem reading muon file line ...."<< countFile << ". Perhaps you've exhausted the events in " << fFileName << std::endl; exit(0);
	    }
	  else
	    {
	      //	  std::cout << "LArIATFileGen: getline() gives "<< line << " for event " << countFile << std::endl; 
	    }
	  countFile++;
	  
	  LOG_DEBUG("LArIATFileGen: countFile is ") << countFile <<std::endl;
	  char * cstr, *ptok;
	  
      // Split this line into tokens
	  cstr = new char [line.size()+1];
	  strcpy (cstr, line.c_str());
	  // cstr now contains a c-string copy of str
	  ptok=strtok (cstr,"*");
	  unsigned int fieldCount = 0;
	  unsigned int posIndex = 0;
	  unsigned int pIndex = 0;
	  while (ptok!=NULL)
	    {
	      // std::cout << ptok << std::endl;
	      ptok=strtok(NULL,"*");
	      if (fieldCount==9 || fieldCount==10 || fieldCount==11)
		{
		  p[pIndex] = atof(ptok); pIndex++;
		  //   std::cout << ptok << std::endl;
		}
	      if (fieldCount==6 || fieldCount==7 || fieldCount==8)
		{
		  x[posIndex] = atof(ptok); 
		  // make the z axis point up for x, as with p
		  if (posIndex==2) {x[posIndex] = -1.0*x[posIndex];}
		  posIndex++;
		}
	      if (fieldCount==12)
		{
		  q = atof(ptok); 
		}
	      fieldCount++;
	    }
	  
	  delete[] cstr;  

	}
      else if (fMuonsFileType.compare("root")==0) // from root file
	{
	    /*
	      // Don't use this yet. Keep the specific branch-by-branch identification.
	      for (unsigned int ii=0;ii<fBranchNames.size();ii++)
	      {
	       TNtuple->SetBranchAddress(fBranchNames[ii], x+ii); 
	      }
	    */
	  //	  TNtuple->ResetBranchAddresses();
	  TNtuple->GetEntry(countFile);
	  //TNtuple->Show(countFile);

	  x.SetXYZ(xdet,ydet,-zdet); // as with txt file, make z point up.
	  // Watch for units change to mm in Modern JdJ files!!
	  // This is for pre Spring-2012 JdJ Ntuples.
	  p.SetXYZ(pxtmp,pytmp,pztmp);
	  q = charge;

	  countFile++;

	} // End read.
      
      static TDatabasePDG  pdgt;
      pdgLocal = -q*fPDG[i];
      
      TParticlePDG* pdgp = pdgt.GetParticle(pdgLocal);
      if (pdgp) m = pdgp->Mass();


      //       std::cout << "set the position "<<std::endl;
      // This gives coordinates at the center of the 300mx300m plate that is 3m above top of 
      // cavern. Got these by histogramming deJong's xdet,ydet,zdet.
      const double cryoGap = 15.0; 
      x[0] -= fXYZ_Off[0];
      x[1] -= fXYZ_Off[1];
      x[2] -= fXYZ_Off[2]; // 3 for plate height above top of cryostat.
      // Now, must rotate to TPC coordinates. Let's orient TPC axis along z axis,
      // Cosmics, mostly going along deJong's +z axis must be going along TPC -y axis.
      x.RotateX(-M_PI/2);
      p.RotateX(-M_PI/2);
      //add vector of the position of the center of the point between Cryostats
      // level with top. (To which I've added 3m - in above code - in height.)
      // This is referenced from origin at center-right of first cryostat.
      art::ServiceHandle<geo::Geometry> geom;
      TVector3 off3(geom->CryostatHalfWidth()*0.01,geom->CryostatHalfHeight()*0.01,geom->CryostatLength()*0.01+cryoGap*0.01/2.0) ;
      x += off3;

      TLorentzVector pos(x[0]*100.0, x[1]*100.0, x[2]*100.0, 0.0);
      TLorentzVector pvec(p[0]*1000.0,p[1]*1000.0,p[2]*1000.0,std::sqrt(p.Mag2()*1000.0*1000.0+m*m));
      std::cout << "x[m] and p [TeV] are " << std::endl;
      x.Print();
      p.Print();
    
      int trackid = -1*(i+1); // set track id to -i as these are all primary particles and have id <= 0
      std::string primary("primary");
      simb::MCParticle part(trackid, pdgLocal, primary);
      part.AddTrajectoryPoint(pos, pvec);
      
      //       std::cout << "add the particle to the primary" << std::endl;
      
      mct.Add(part);

    }//end loop over particles
    
    return;
  }

  DEFINE_ART_MODULE(LArIATFileGen)

}//end namespace evgen
////////////////////////////////////////////////////////////////////////
