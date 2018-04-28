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
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCParticle.h"

// lar includes
#include "larcore/Geometry/Geometry.h"
#include "larcoreobj/SummaryData/RunData.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"

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


// ##helper class to avoid hardcoding branches. For now, for lack of a smarter solution am assuiming that TrackPresent is the last variable.
// ## all things will go to hell if it's not.






namespace evgen {

  
  
  class det_info{
    
public:
    
  det_info(std::vector< std::string > var_names, std::string detname);
  void SetBranches(TTree * TNtuple,std::vector< std::string > var_names);
  

  
  std::vector< TBranch *> b_obj;
  std::vector< TString > branchname;
  double x,y,z,t,Px,Py,Pz,PDGid,ParentID;
  bool TrackPresent;
  int nvars;
  
};

det_info::det_info(std::vector< std::string > var_names, std::string detname)
{
 nvars =var_names.size();
 b_obj.resize(nvars);
 branchname.resize(nvars);
// value.resize(nvars); 
 for(int i=0;i<nvars;i++)
 {
   branchname[i]=var_names[i]+detname;
 }
   
 TrackPresent=false;
}
  
  
void det_info::SetBranches(TTree * TNtuple,std::vector< std::string > var_names)
{
  for(unsigned int br_id=0;br_id<var_names.size();br_id++)
      {
      if(var_names[br_id]=="x")	
        TNtuple->SetBranchAddress(branchname[br_id], &x, &b_obj[br_id]);
      else if(var_names[br_id]=="y")	 
	TNtuple->SetBranchAddress(branchname[br_id], &y, &b_obj[br_id]);
      else if(var_names[br_id]=="z")	 
	TNtuple->SetBranchAddress(branchname[br_id], &z, &b_obj[br_id]);
      else if(var_names[br_id]=="t")	 
	TNtuple->SetBranchAddress(branchname[br_id], &t, &b_obj[br_id]);
      else if(var_names[br_id]=="Px")	
        TNtuple->SetBranchAddress(branchname[br_id], &Px, &b_obj[br_id]);
      else if(var_names[br_id]=="Py")	 
	TNtuple->SetBranchAddress(branchname[br_id], &Py, &b_obj[br_id]);
      else if(var_names[br_id]=="Pz")	 
	TNtuple->SetBranchAddress(branchname[br_id], &Pz, &b_obj[br_id]);
      else if(var_names[br_id]=="PDGid")	 
	TNtuple->SetBranchAddress(branchname[br_id], &PDGid, &b_obj[br_id]);
       else if(var_names[br_id]=="TrackPresent")	 
	TNtuple->SetBranchAddress(branchname[br_id], &TrackPresent, &b_obj[br_id]);
   
      }
      
  
  
}
  
  
  
    
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
  
  

#if defined __clang__
  #pragma clang diagnostic push
  #pragma clang diagnostic ignored "-Wunused-private-field"
#endif
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
    std::vector<std::string> fDetectorNames; 
    double fTiWindowBeamx;
    double fTiWindowBeamy;
    double fTiWindowBeamz;
    double fTiWindowTPCx;
    double fTiWindowTPCy;
    double fTiWindowTPCz;
    int fEventSpillOffset;  // Where in file to start.
    int fEventsPerSpill;
    int fEventsPerFile; 
    int fEventFileOffset; 
    std::string fTriggerCut;
    std::string fInnerDetName;
    std::string fOuterDetName;
    std::string fTimeDetName;
    bool fUseTrigger; 
  
    
    std::ifstream *fMuonFile;
    TFile *fMuonFileR;
    TTree *TNtuple;
    TTree *TNtupleSmall;
    TTree *TNtupleSorted;
    unsigned int countFile;

    //Float_t xtmp, ytmp, ztmp; // unused
    Float_t pxtmp, pytmp, pztmp;
    Float_t charge;
    //Float_t         E; // unused
    //Float_t         costheta; // unused
    //Float_t         phi; // unused
    Float_t         xdet;
    Float_t         ydet;
    Float_t         zdet;

    
    
    int EventID,TrackID;
    TBranch * b_eventid, * b_trackid;
    TEntryList *elist; //used for triggered events
    
    std::vector < det_info > detectors;
//     //Det1
//     double xDet1,yDet1,zDet1,tDet1,PxDet1,PyDet1,PzDet1,PDGidDet1,ParentIDDet1;bool TrackPresentDet1;
//     
//     TBranch *b_xDet1,*b_yDet1,*b_zDet1,*b_tDet1,*b_PxDet1,*b_PyDet1,*b_PzDet1,*b_PDGidDet1,*b_ParentIDDet1,*b_TrackPresentDet1;
//     //Det2
//     double xDet2,yDet2,zDet2,tDet2,PxDet2,PyDet2,PzDet2,PDGidDet2,ParentIDDet2;bool TrackPresentDet2;
//     TBranch *b_xDet2,*b_yDet2,*b_zDet2,*b_tDet2,*b_PxDet2,*b_PyDet2,*b_PzDet2,*b_PDGidDet2,*b_ParentIDDet2,*b_TrackPresentDet2;
//     //Det3
//     double xDet3,yDet3,zDet3,tDet3,PxDet3,PyDet3,PzDet3,PDGidDet3,ParentIDDet3;bool TrackPresentDet3;
//     TBranch *b_xDet3,*b_yDet3,*b_zDet3,*b_tDet3,*b_PxDet3,*b_PyDet3,*b_PzDet3,*b_PDGidDet3,*b_ParentIDDet3,*b_TrackPresentDet3;
//     //Det4
//     double xDet4,yDet4,zDet4,tDet4,PxDet4,PyDet4,PzDet4,PDGidDet4,ParentIDDet4;bool TrackPresentDet4;
//     TBranch *b_xDet4,*b_yDet4,*b_zDet4,*b_tDet4,*b_PxDet4,*b_PyDet4,*b_PzDet4,*b_PDGidDet4,*b_ParentIDDet4,*b_TrackPresentDet4;
//     //TOFus
//     double xTOFus,yTOFus,zTOFus,tTOFus,PxTOFus,PyTOFus,PzTOFus,PDGidTOFus,ParentIDTOFus;bool TrackPresentTOFus;
//     TBranch *b_xTOFus,*b_yTOFus,*b_zTOFus,*b_tTOFus,*b_PxTOFus,*b_PyTOFus,*b_PzTOFus,*b_PDGidTOFus,*b_ParentIDTOFus,*b_TrackPresentTOFus;
//     //TODdsVert
//     double xTODdsVert,yTODdsVert,zTODdsVert,tTODdsVert,PxTODdsVert,PyTODdsVert,PzTODdsVert,PDGidTODdsVert,ParentIDTODdsVert;bool TrackPresentTODdsVert;
//     TBranch *b_xTODdsVert,*b_yTODdsVert,*b_zTODdsVert,*b_tTODdsVert,*b_PxTODdsVert,*b_PyTODdsVert,*b_PzTODdsVert,*b_PDGidTODdsVert,*b_ParentIDTODdsVert,*b_TrackPresentTODdsVert;
//     //TODdsHorz
//     double xTODdsHorz,yTODdsHorz,zTODdsHorz,tTODdsHorz,PxTODdsHorz,PyTODdsHorz,PzTODdsHorz,PDGidTODdsHorz,ParentIDTODdsHorz;bool TrackPresentTODdsHorz;
//     TBranch *b_xTODdsHorz,*b_yTODdsHorz,*b_zTODdsHorz,*b_tTODdsHorz,*b_PxTODdsHorz,*b_PyTODdsHorz,*b_PzTODdsHorz,*b_PDGidTODdsHorz,*b_ParentIDTODdsHorz,*b_TrackPresentTODdsHorz;
//     //Halo
//     double xHalo,yHalo,zHalo,tHalo,PxHalo,PyHalo,PzHalo,PDGidHalo,ParentIDHalo;bool TrackPresentHalo;
//     TBranch *b_xHalo,*b_yHalo,*b_zHalo,*b_tHalo,*b_PxHalo,*b_PyHalo,*b_PzHalo,*b_PDGidHalo,*b_ParentIDHalo,*b_TrackPresentHalo;
//     //HaloHole
//     double xHaloHole,yHaloHole,zHaloHole,tHaloHole,PxHaloHole,PyHaloHole,PzHaloHole,PDGidHaloHole,ParentIDHaloHole;bool TrackPresentHaloHole;
//     TBranch *b_xHaloHole,*b_yHaloHole,*b_zHaloHole,*b_tHaloHole,*b_PxHaloHole,*b_PyHaloHole,*b_PzHaloHole,*b_PDGidHaloHole,*b_ParentIDHaloHole,*b_TrackPresentHaloHole;
//     //TiWindow
//     double xTiWindow,yTiWindow,zTiWindow,tTiWindow,PxTiWindow,PyTiWindow,PzTiWindow,PDGidTiWindow,ParentIDTiWindow;bool TrackPresentTiWindow;
//     TBranch *b_xTiWindow,*b_yTiWindow,*b_zTiWindow,*b_tTiWindow,*b_PxTiWindow,*b_PyTiWindow,*b_PzTiWindow,*b_PDGidTiWindow,*b_ParentIDTiWindow,*b_TrackPresentTiWindow;
    
    std::vector<  EventFrame  > fEventFrames;
    int fEventCounter;

  };
#if defined __clang__
  #pragma clang diagnostic pop
#endif
  
  
  

  
  
  
} // namespace



namespace evgen{

  //____________________________________________________________________________
  LArIATFileGen::LArIATFileGen(fhicl::ParameterSet const& pset)
     : fSeed(314159)					                          
     , fFileName         (pset.get< std::string     	     	>("FileName")         	)      
     , fMuonsFileType    (pset.get< std::string             	>("MuonsFileType")    	)
     , fTreeName         (pset.get< std::string   	     	>("TreeName")         	)      
     , fBranchNames      (pset.get< std::vector<std::string> 	>("BranchNames")      	)
     , fDetectorNames    (pset.get< std::vector<std::string> 	>("DetectorNames")     	)
     , fTiWindowBeamx    (pset.get<double			>("TiWindowBeamx")     	)
     , fTiWindowBeamy    (pset.get<double			>("TiWindowBeamy")     	)
     , fTiWindowBeamz    (pset.get<double			>("TiWindowBeamz")     	)
     , fTiWindowTPCx     (pset.get<double			>("TiWindowTPCx")     	) 
     , fTiWindowTPCy     (pset.get<double			>("TiWindowTPCy")     	) 
     , fTiWindowTPCz     (pset.get<double			>("TiWindowTPCz")     	) 
     , fEventSpillOffset (pset.get<int				>("EventSpillOffset")	)
     , fEventsPerSpill (pset.get<int				>("EventsPerSpill")	)
     , fEventsPerFile (pset.get<int				>("EventsPerFile")	)
     , fEventFileOffset (pset.get<int				>("EventFileOffset")	) 
     , fTriggerCut         (pset.get< std::string     	     	>("TriggerCut")       	)
     , fInnerDetName       (pset.get< std::string     	 	>("InnerDetName")      	)  
     , fOuterDetName       (pset.get< std::string     	     	>("OuterDetName")      	)  
     , fTimeDetName        (pset.get< std::string     	     	>("TimeDetName")      	) 
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
    
     // constructor decides if initialized value is a path or an environment variable
    //std::string fname;
    //cet::search_path sp("FW_SEARCH_PATH");
    //sp.find_file(fFileName, fname);
    //fMuonFileR = new TFile(fname.c_str(),"READ");
    fMuonFileR = new TFile(fFileName.c_str(),"READ");
    
    TNtuple = (TTree*)(fMuonFileR->Get(fTreeName.c_str()));

    //Set Branches for all detectors to extract the data per spill.
    
    TNtuple->SetBranchAddress("EventID", &EventID, &b_eventid);
    TNtuple->SetBranchAddress("TrackID", &TrackID, &b_trackid);
    
    //reserve memory to avoid vector memory location shifting, found out by Dr E. Church, esq.
    detectors.reserve(fDetectorNames.size());
    
    for (unsigned int detId=0;detId<fDetectorNames.size();detId++ )
    {
     std::cout << " setting detector branches for " << detId << " " << fDetectorNames[detId] << std::endl;
      
     detectors.push_back( det_info(fBranchNames,fDetectorNames[detId]));
     detectors[detId].SetBranches(TNtuple,fBranchNames);
     std::cout << " setting detector branches for " << detId << " " << fDetectorNames[detId] << " " << detectors[detId].branchname[0]  << " detinfo address " << &detectors[0] << std::endl;
      
    }
     
    TFile * testfile=new TFile("out.root","RECREATE");
    TNtupleSmall = TNtuple->CloneTree(0);    
    Long64_t nentries = TNtuple->GetEntries();
    art::ServiceHandle<geo::Geometry> geom;

    //Create small TTree for sorting.
    // check that neiter of these numbers goes out of bounds.
    //nentries = (nentries > fEventNumberOffset+fEventsPerSpill ) ? fEventNumberOffset+fEventsPerSpill : nentries;
   
    
    long start_event=fEventSpillOffset*fEventsPerSpill;
     std::cout << " creating a subtree for sorting - this make take a while. from evt: " <<  fEventSpillOffset*fEventsPerSpill << " to: " << start_event+fEventsPerSpill << std::endl;
    for (Long64_t i=0;i<nentries; i++) 
    {
      TNtuple->GetEntry(i);

      if (EventID > start_event   && EventID < start_event+fEventsPerSpill) TNtupleSmall->Fill();
      if(EventID > start_event+fEventsPerSpill) break;
    }
    
    TNtupleSmall->Write();
    
    std::cout << "old tree size, new tree size: " << TNtuple->GetEntries() << " " << TNtupleSmall->GetEntries() << std::endl;
    
    auto const* detprop = lar::providerFrom<detinfo::DetectorPropertiesService>();
    
    const double frametime=detprop->SamplingRate()*detprop->NumberTimeSamples()/1000.; 
    int frameindex=1;
    
    std::vector< int > locPDG;
    std::vector< TLorentzVector > locXYZ;
    std::vector< TLorentzVector > locMOM;
        
    TDatabasePDG  pdgt;
    double m = 0.;
    int time_index=0;  // Index of the detector time sorting is based on
    int time_tag=0;    // index for the position of t variable in fBranchNames

    std::cout << " building index on time on the subtree - this make also take a while " << std::endl;

    TNtupleSorted = TNtupleSmall->CloneTree(0); 
       
    for(unsigned int ix=0;ix<fDetectorNames.size();ix++)
    {
//     std::cout << "Detector " << fDetectorNames[ix] << "   Index: " << ix << std::endl;
     if(fDetectorNames[ix]==fTimeDetName) time_index=ix;
    }
    
    for(unsigned int it=0;it<fBranchNames.size();it++)
    {
//     std::cout << "Variable " << fBranchNames[it] << "   Index: " << it << std::endl;
     if(fBranchNames[it]=="t") time_tag=it;
    }

//    std::cout << "Time Detector Index " << time_index << std::endl;
//    std::cout << "Variable Index " << time_tag << std::endl;
    std::cout << "Selected branch for time ordering: " << (detectors[time_index].branchname[time_tag]).Data() << std::endl;
    
    TNtupleSmall->BuildIndex(Form("%s*10000000",(detectors[time_index].branchname[time_tag]).Data()));  //sorting on integers, need microseconds
    TTreeIndex *index = (TTreeIndex*)TNtupleSmall->GetTreeIndex();

    for(int qq=0; qq<TNtupleSmall->GetEntries(); qq++)
    {
       int tree_entry;
       tree_entry=TNtupleSmall->LoadTree( index->GetIndex()[qq] );
       TNtupleSmall->GetEntry(tree_entry);
       TNtupleSorted->Fill();
//       std::cout << "Step: " << qq << "   Event ID: " << EventID << " Time: " << detectors[time_index].t << std::endl;       
    }
   
   std::cout << "Size of sorted Ntuple: " << TNtupleSorted->GetEntries() << std::endl; 

    TNtupleSorted->Draw(">>elist", fTriggerCut.data(), "entrylist");
    elist = (TEntryList*)gDirectory->Get("elist");
    Long64_t ListNEntries=elist->GetN();
    if(fUseTrigger)
    std::cout <<  "Number of entries in Event List " << ListNEntries << std::endl;  

    double last_tHalo_trig=0;
        
    int INNER_DET_INDEX=0;  // position in the detectors vector of the inner detector
    int OUTER_DET_INDEX=0;  // position in the detectors vector of the outer detector
    double x_shift=0;       // distance between TPC and beam reference frames in the x axis
    double y_shift=0;       // distance between TPC and beam reference frames in the x axis
    double z_shift=0;       // distance between TPC and beam reference frames in the x axis
/*    double x_inner_shift=0; //x position of the center of inner detector in Jason's coordinate frame
    double x_outer_shift=0; //x position of the center of outer detector in Jason's coordinate frame
    double y_inner_shift=0; //y position of the center of inner detector in Jason's coordinate frame
    double y_outer_shift=0; //y position of the center of outer detector in Jason's coordinate frame
    double z_inner_shift=0; //z position of the center of inner detector in Jason's coordinate frame
    double z_outer_shift=0; //z position of the center of outer detector in Jason's coordinate frame*/
    
    for(unsigned int ix=0;ix<fDetectorNames.size();ix++)
    {
     if(fDetectorNames[ix]==fInnerDetName) 
      INNER_DET_INDEX=ix;
     if(fDetectorNames[ix]==fOuterDetName) 
      OUTER_DET_INDEX=ix;
    }
  
    x_shift=fTiWindowTPCx-fTiWindowBeamx;
    y_shift=fTiWindowTPCy-fTiWindowBeamy;
    z_shift=fTiWindowTPCz-fTiWindowBeamz;

    //if running with trigger then loop only on selected entries, i.e. ListNEntries
    // ir running without trigger then loop on all events and save them by frame.
    Long64_t MainIndexEntries = (fUseTrigger) ? ListNEntries : index->GetN() - 1;
     
    //loop over selected trigger entries.
    for (int i = 0;  i<MainIndexEntries ; i++)
      {
	int tree_entry; 
	
	if(fUseTrigger)  tree_entry=elist->GetEntry(i);
	else   tree_entry=TNtupleSorted->LoadTree( index->GetIndex()[i] );
	  
	TNtupleSorted->GetEntry(tree_entry);

        std::cout << std::endl;
        std::cout <<  "  Selected Trigger event # " << i << " Event ID: " << EventID << " Event Time: " << detectors[time_index].t << std::endl;
        	
 	if(detectors[time_index].t<last_tHalo_trig+frametime/1000000) // cannot retrigger on same frame.
 	  continue;
	
	if(fEventFrames.size() >=  (unsigned int)fEventsPerFile)  // break if too many events.
	  break;
	
	if(fUseTrigger)
        {
	  last_tHalo_trig=detectors[time_index].t;
        }
	else
	{
	   while(detectors[time_index].t > frameindex*frametime/1000000)  
	      {frameindex++; }
	   last_tHalo_trig=frameindex*frametime/1000000;
	}

//	std::cout << "Time cut: " << Form("%s > %lf && %s < %lf",(detectors[time_index].branchname[time_tag]).Data(),last_tHalo_trig-frametime/1000000,(detectors[time_index].branchname[time_tag]).Data(),last_tHalo_trig+frametime/1000000) << std::endl;

	TNtupleSorted->Draw(">>loclist",Form("%s > %lf && %s < %lf",(detectors[time_index].branchname[time_tag]).Data(),last_tHalo_trig-frametime/1000000,(detectors[time_index].branchname[time_tag]).Data(),last_tHalo_trig+frametime/1000000),"entrylist");
	
	TEntryList * loclist = (TEntryList*)gDirectory->Get("loclist");
        Long64_t ListlocNEntries=loclist->GetN();
	
	std::cout << "Found " << ListlocNEntries << " particles in frame " << std::endl;
	 
	locPDG.clear();
	locXYZ.clear();
	locMOM.clear(); 
	 
	for (int xx = 0;  xx<ListlocNEntries ; xx++)
	{ 
          int tree_entry=loclist->GetEntry(xx);
	  TNtupleSorted->GetEntry(tree_entry);
//	  std::cout <<  "  Selected Frame event " << EventID << " Time:  " << detectors[time_index].t << " tree_entry " << tree_entry << std::endl;
	
	  double x,y,z,Px,Py,Pz,PDG ; //z
         //double x_shift,y_shift,z_shift;
         //double x_extra_shift,y_extra_shift; // shift between the x(y) position of TiWindow and x(y) position of the inner or outer detector considered
      
	 if(detectors[INNER_DET_INDEX].TrackPresent)
         {
	    x=detectors[INNER_DET_INDEX].x;y=detectors[INNER_DET_INDEX].y;z=detectors[INNER_DET_INDEX].z;//z=zHalo;
	    Px=detectors[INNER_DET_INDEX].Px;Py=detectors[INNER_DET_INDEX].Py;Pz=detectors[INNER_DET_INDEX].Pz;	
	    PDG=detectors[INNER_DET_INDEX].PDGid;
//          x_shift=x_inner_shift;
//          y_shift=y_inner_shift;
//          z_shift=z_inner_shift;
//            std::cout << "Track in inner detector " << std::endl;
          }
	  else
	  {
	    x=detectors[OUTER_DET_INDEX].x;y=detectors[OUTER_DET_INDEX].y;z=detectors[OUTER_DET_INDEX].z;//z=zHaloHole;
	    Px=detectors[OUTER_DET_INDEX].Px;Py=detectors[OUTER_DET_INDEX].Py;Pz=detectors[OUTER_DET_INDEX].Pz;	
	    PDG=detectors[OUTER_DET_INDEX].PDGid;
//          x_shift=x_outer_shift;
//          y_shift=y_outer_shift;
//          z_shift=z_outer_shift;
	    std::cout << "Track in outer detector " <<OUTER_DET_INDEX << " "<< fDetectorNames[OUTER_DET_INDEX] << " " << fOuterDetName << " PDG " <<detectors[OUTER_DET_INDEX].PDGid << std::endl;
	  }
      
/*        x_extra_shift=x_shift+112.43; // x_shift - Ti window x position. To account for
        y_extra_shift=y_shift-0.0; // y_shift - Ti window y position
        
       std::cout << " X Shift: " << x_shift << "   Y Shift: " << y_shift << "   Z Shift: " << z_shift << std::endl;  
        std::cout << " X Extra Shift: " << x_extra_shift << "   Y Extra Shift: " << y_extra_shift << std::endl;*/        
 
	  TParticlePDG* pdgp = pdgt.GetParticle(PDG);
	  if (pdgp) m = pdgp->Mass();
	
	  double Ener=std::sqrt((Px*Px+Py*Py+Pz*Pz)/1000./1000.+m*m);

      // titatnium window position wrt. TPCActive (based on gdml file wdisk)
      // in gdml (0,035 ,0.0 ,-84.54646643)
      // to transfer into LArSOFT coordinates, 
      // X changes by HalfWidth + difference btw TiWindow and TPCactive + difference btw TiWindow and Jason coordinate system's center
      // Z changes by distance btw TPCactive and Tiwindow + distance btw TiWindow and inner or outer detector considered. 
      // geom->CryostatHalfHeight()*0.01,geom->CryostatLength()

	  locPDG.push_back(PDG);
//	locXYZ.push_back(TLorentzVector(x/10.-0.435+geom->DetHalfWidth(),y/10.+0.2,-84.54646+geom->DetLength()/2.-1,(detectors[INNER_DET_INDEX].t-last_tHalo_trig)*1e9)); //convert to ns
//        locXYZ.push_back(TLorentzVector(x/10.-x_shift+x_extra_shift+geom->DetHalfWidth()+0.035,y/10.+y_shift+y_extra_shift,-39.546-(716.792-z_shift),(detectors[INNER_DET_INDEX].t-last_tHalo_trig)*1e9)); //convert to ns
          locXYZ.push_back(TLorentzVector(x/10.+x_shift,y/10.+y_shift,z/10.+z_shift,(detectors[time_index].t-last_tHalo_trig)*1e9)); //convert to ns
	  locMOM.push_back(TLorentzVector(Px/1000.,Py/1000.,Pz/1000.,Ener));
	  
	//first back track to find particles that are one frame before:
	//while(tHalo > tHalo_trig-frametime/1000000)
	 // {
	  } //end for loop inside of frame

	if(locPDG.size() && locXYZ.size() && locMOM.size() )  // skipping non needed files
	{  
          fEventFrames.push_back(EventFrame(locPDG,locXYZ,locMOM)); 
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
     
     if (fMuonsFileType.compare("root")==0) {
      fMuonFileR->Close(); }
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
