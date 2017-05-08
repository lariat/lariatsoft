////////////////////////////////////////////////////////////////////////
// Class:       MassSelectionFilter
// Module Type: filter
// File:        MassSelectionFilter_module.cc
//
// Generated at Tue Sep  6 19:32:40 2016 by Elena Gramellini using artmod
// from cetpkgsupport v1_10_02.
//
//  With this filter I want to select a sample of perfect protons
//  [ x ] Keep events with only 1 TOF
//  [ x ] Keep events with only 1 WC Track
//  [ x ] Keep events with max  4 TPC Track in the first 14 cm
//  [ x ] Keep events with 1+ TPC point in the first 2 cm
//
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include <memory>
#include <TH1F.h>
#include <TH2F.h>

#include "art/Framework/Services/Optional/TFileService.h"
#include "LArIATDataProducts/TOF.h"
#include "LArIATDataProducts/WCTrack.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/SpacePoint.h"

#include "canvas/Persistency/Common/FindOneP.h" 
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "canvas/Persistency/Common/Ptr.h" 
#include "canvas/Persistency/Common/PtrVector.h" 


#include "larcore/Geometry/Geometry.h"
#include "lardata/Utilities/AssociationUtil.h"



class MassSelectionFilter;

class MassSelectionFilter : public art::EDFilter {
public:
  explicit MassSelectionFilter(fhicl::ParameterSet const & p);
  // The destructor generated by the compiler is fine for classes
  // without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  MassSelectionFilter(MassSelectionFilter const &) = delete;
  MassSelectionFilter(MassSelectionFilter &&) = delete;
  MassSelectionFilter & operator = (MassSelectionFilter const &) = delete;
  MassSelectionFilter & operator = (MassSelectionFilter &&) = delete;

  void reconfigure(fhicl::ParameterSet const & p) override;
  void beginJob() override;
  void endJob() override;
  // Required functions.
  bool filter(art::Event & e) override;
  int nPassed = 0 ;
  int nFailed = 0 ;
  
private:

  // Declare member data here.
  //---------- Filter Parameters ----------
  std::string fTOFModuleLabel  ;
  size_t      fnTOFObjects     ;
 
  std::string fWCTrackLabel    ;
  size_t      fNumberWCTrack;

  std::string fTrackModuleLabel;  
  double fUpstreamZPosition;  
  double fnTracksUpstream;    
  double fUpstreamZThreshold; 
  double fMassUpperLimit;  
  double fMassLowerLimit;
        
  art::ServiceHandle<geo::Geometry> fGeo;  

  //---------- Histos ----------
  TH1F* hBeamlineMassAfterCut;
  TH1F* hBeamlineMassBeforeCut;
  TH2F* hTOFVsMomAfterCut;
  TH2F* hTOFVsMomBeforeCut;

};

// ---------------------- Begin Job ---------------------------
void MassSelectionFilter::endJob()
{
  std::cout<<"&&& Number of events in Mass filter = "<<nPassed + nFailed<<"; passed = "<<nPassed<<"; failed = "<<nFailed<<std::endl;
}

// ---------------------- Begin Job ---------------------------
void MassSelectionFilter::beginJob()
{

  // Implementation of optional member function here.
  art::ServiceHandle<art::TFileService> tfs;
  hTOFVsMomAfterCut      = tfs->make<TH2F>("hTOFvsMomAfterCuts"     , "hTOFvsMomAfterCuts; WC Momentum [MeV/c];TOF [ns];", 200, 0, 2000, 100, 0, 100); 
  hTOFVsMomBeforeCut     = tfs->make<TH2F>("hTOFvsMomBeforeCuts"    , "hTOFvsMomBeforeCuts;WC Momentum [MeV/c];TOF [ns];", 200, 0, 2000, 100, 0, 100); 
  hBeamlineMassAfterCut  = tfs->make<TH1F>("BeamlineMassAfterCut" ,"BeamLine Mass After Cuts" ,400,0,2000);  
  hBeamlineMassBeforeCut = tfs->make<TH1F>("BeamlineMassBeforeCut","BeamLine Mass Before Cuts",400,0,2000);  
  
}


MassSelectionFilter::MassSelectionFilter(fhicl::ParameterSet const & p)
// :
// Initialize member data here.
{
  this->reconfigure(p);
  // Call appropriate produces<>() functions here.
}

bool MassSelectionFilter::filter(art::Event & evt)
{
  
  // ####################################################
  // ### Getting the Time of Flight (TOF) Information ###
  // ####################################################
  art::Handle< std::vector<ldp::TOF> > TOFColHandle;
  std::vector<art::Ptr<ldp::TOF> > tof;

  

  if(!evt.getByLabel(fTOFModuleLabel,TOFColHandle)) return false;
  art::fill_ptr_vector(tof, TOFColHandle);
   
  if(tof.size() < 1 || tof.size() > 1)      {nFailed++; return false; }
  //  if (!tof[0]->NTOF()) {nFailed++; return false; }


  // #############################################
  // ### Getting the Momentum (WC) Information ###
  // #############################################
  
  art::Handle< std::vector<ldp::WCTrack> > wctrackHandle;
  std::vector<art::Ptr<ldp::WCTrack> > wctrack;

  if(!evt.getByLabel(fWCTrackLabel, wctrackHandle)){nFailed++; return false; } 
  art::fill_ptr_vector(wctrack, wctrackHandle);

  if (wctrack.size() > 1 || wctrack.size() < 1){nFailed++; return false; }


   
  double tofObject[100];            //<---The TOF calculated (in ns?) for this TOF object                                                                                           
  float fDistanceTraveled = 6.652; 
  double reco_momo = wctrack[0]->Momentum();

  
  size_t tof_counter = 0;
  for(size_t i = 0; i < tof.size(); i++) {
    size_t number_tof = tof[i]->NTOF();
    for (size_t tof_idx = 0; tof_idx < number_tof; ++tof_idx) {
      tofObject[tof_counter] =  tof[i]->SingleTOF(tof_idx);
      ++tof_counter;
    } // loop over TOF
    
  }//<---End tof_count loop

  double reco_tof = tofObject[0];
  
  float mass = reco_momo*pow(reco_tof*0.299792458*0.299792458*reco_tof/(fDistanceTraveled*fDistanceTraveled) - 1 ,0.5);
 

  if(std::isnan(mass)) return false;

  hTOFVsMomBeforeCut->Fill(reco_momo,reco_tof);
  if(!std::isnan(mass)) hBeamlineMassBeforeCut->Fill(mass);     

  
  if (mass < fMassLowerLimit ) {nFailed++; return false; }
  if (mass > fMassUpperLimit ) {nFailed++; return false; }



  if(!std::isnan(mass)) hBeamlineMassAfterCut->Fill(mass);  
  hTOFVsMomAfterCut->Fill(reco_momo,reco_tof);
  nPassed++;
  return true;   
}


void MassSelectionFilter::reconfigure(fhicl::ParameterSet const & p)
{                                                                                                   
  fTOFModuleLabel = p.get< std::string >("TOFModuleLabel");
  fnTOFObjects    = p.get<   size_t    >("nTOFObjects", 1);
 
  fWCTrackLabel   = p.get< std::string >("WCTrackLabel");
  fNumberWCTrack  = p.get<  size_t     >("NumberWCTrack",   1);

  fMassUpperLimit = p.get<   double    >("MassUpperLimit", 1000.0);
  fMassLowerLimit = p.get<   double    >("MassLowerLimit",  900.0);
}

DEFINE_ART_MODULE(MassSelectionFilter)


