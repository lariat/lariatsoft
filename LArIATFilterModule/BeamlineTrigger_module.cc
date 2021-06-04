////////////////////////////////////////////////////////////////////////
// Class:       BeamlineTrigger
// Module Type: filter
// File:        BeamlineTrigger_module.cc
//
// Generated at Sun Mar 13 15:42:13 2016 by Irene Nutini using artmod
// from cetpkgsupport v1_08_06.
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
#include "art_root_io/TFileService.h"


// ###########################
// ### LArIATsoft Includes ###
// ###########################
#include "LArIATDataProducts/TOF.h"
#include "LArIATDataProducts/WCTrack.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "larcoreobj/SimpleTypesAndConstants/RawTypes.h" // raw::ChannelID_t
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/CryostatGeo.h"
#include "larcorealg/Geometry/TPCGeo.h"
#include "larcorealg/Geometry/PlaneGeo.h"
#include "larcorealg/Geometry/WireGeo.h"
#include "lardataobj/RecoBase/Wire.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardata/Utilities/AssociationUtil.h"

// ####################
// ### C++ Includes ###
// ####################
#include <iostream>
#include <memory>
#include <TH1F.h>
#include <TH2F.h>

class BeamlineTrigger;

class BeamlineTrigger : public art::EDFilter {
public:
  explicit BeamlineTrigger(fhicl::ParameterSet const & p);
  // The destructor generated by the compiler is fine for classes
  // without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  BeamlineTrigger(BeamlineTrigger const &) = delete;
  BeamlineTrigger(BeamlineTrigger &&) = delete;
  BeamlineTrigger & operator = (BeamlineTrigger const &) = delete;
  BeamlineTrigger & operator = (BeamlineTrigger &&) = delete;

  // Required functions.
  bool filter(art::Event & e) override;

  // Selected optional functions.
  void beginJob() override;
  void reconfigure(fhicl::ParameterSet const & p) ;

private:

  // Declare member data here.
    std::string fTOFModuleLabel;		// Name of the producer that made the TOF objects
    double fnTOFObjects;
    
    std::string fWCTrackModuleLabel;
    double fminNumberWCTrack;
    double fWCLowerBound; //<---Lowest WC value we can have
    double fWCUpperBound; //<---Upper WC value we can have
    
    std::string     fCalDataModuleLabel;
    double 		fnCalWireObjects;
    
    
    TH1F* fReco_Pz;
    TH1F* fTOFBeforeCut;
    TH2F* fPiMuPzVsTOF;
    TH2F* fPzVsTOF;


};


BeamlineTrigger::BeamlineTrigger(fhicl::ParameterSet const & p)
: EDFilter(p)
// Initialize member data here.
{
  // Call appropriate produces<>() functions here.
    this->reconfigure(p);
}

bool BeamlineTrigger::filter(art::Event & e)
{
  // Implementation of required member function here.
    
    // ##########################################
    // ### Reading in the Wire List object(s) ###
    // ##########################################
    art::Handle< std::vector<recob::Wire> > wireVeclHandle;
    std::vector< art::Ptr<recob::Wire> >wireVec;
    
    if(e.getByLabel(fCalDataModuleLabel,wireVeclHandle))
    {art::fill_ptr_vector(wireVec, wireVeclHandle);}
    
    
    //std::cout << "Calwire objects expected " << fnCalWireObjects << std::endl;
    bool GoodTPCValue = false;
    // ### Reject the event if there is no good CalWire info - calwire objects != 480 ###
    //std::cout << "There are: " << wireVec.size() << " calData wire signals for this event" << std::endl;
    if(wireVec.size() <  fnCalWireObjects){return false;}
    
    else { GoodTPCValue =true;}
  
    // ####################################################
    // ### Getting the Time of Flight (TOF) Information ###
    // ####################################################
    art::Handle< std::vector<ldp::TOF> > TOFColHandle;
    std::vector<art::Ptr<ldp::TOF> > tof;
    
    if(e.getByLabel(fTOFModuleLabel,TOFColHandle))
    {art::fill_ptr_vector(tof, TOFColHandle);}
    
    // ######################################
    // ### Get the collection of WCTracks ###
    // ######################################
    art::Handle< std::vector<ldp::WCTrack> > wctrackHandle;
    std::vector<art::Ptr<ldp::WCTrack> > wctrack;
    
    if(e.getByLabel(fWCTrackModuleLabel, wctrackHandle))
    {art::fill_ptr_vector(wctrack, wctrackHandle);}
    
    
    // ### Set a boolian to false for the event ###
    bool GoodTOFValue = false;
    // ### Reject the event if there is no TOF info ###
    if(tof.size() <  fnTOFObjects){return false;}
	 std::cout << "Ntofs " << TOFColHandle->at(0).NTOF() << std::endl;
    if( TOFColHandle->at(0).NTOF() < 1 ) {
		 std::cout << "Empty tof value created " << std::endl;
		 return false;}
    // ################################
    // ### Looping over TOF objects ###
    // ################################
    else{
        for(size_t i = 0; i < tof.size(); i++)
       {
        for (size_t tof_idx = 0; tof_idx < tof[i]->NTOF(); ++tof_idx)
        {
           std::cout << "Tof value " << tof[i]->SingleTOF(tof_idx) << std::endl;
			   fTOFBeforeCut->Fill( tof[i]->SingleTOF(tof_idx) );
        }//<---End tof_idx loop
	  }//<---End i loop
	  GoodTOFValue = true;
     }
    
    bool GoodWCValue =false;
    // ###     If the number of WCTracks in the event       ###
    // ### is less than the min number set, skip this event ###
    if(wctrack.size() < fminNumberWCTrack){return false;}
    
    for(size_t wct_count = 0; wct_count < wctrack.size(); wct_count++)
    {
        // ### Filling histogram with Wire Chamber Track momentum Pz and check momentum range ###
		  std::cout << "wctrack momentum " << wctrack[wct_count]->Momentum() << std::endl;
        fReco_Pz->Fill(wctrack[wct_count]->Momentum());
        if(wctrack[wct_count]->Momentum() > fWCLowerBound && wctrack[wct_count]->Momentum() < fWCUpperBound)
        {GoodWCValue = true;}
        
    }
    
    fPzVsTOF->Fill(wctrackHandle->at(0).Momentum(),TOFColHandle->at(0).SingleTOF(0));
    
    if(GoodTOFValue && GoodWCValue && GoodTPCValue){
        fPiMuPzVsTOF->Fill(wctrackHandle->at(0).Momentum(),TOFColHandle->at(0).SingleTOF(0));
        return true;}
    else {return false;}
}

    
    
    


void BeamlineTrigger::beginJob()
{
  // Implementation of optional member function here.
   art::ServiceHandle<art::TFileService> tfs;
   fTOFBeforeCut = tfs->make<TH1F>("TOFBeforeCut", "TOF(ns)", 100, 0, 100);
   fReco_Pz = tfs->make<TH1F>("Reco_Pz","Reconstructed momentum in XZ plane", 180, 0, 1800);
   fReco_Pz->GetXaxis()->SetTitle("Reconstructed momentum (MeV/c)");
   fReco_Pz->GetYaxis()->SetTitle("Tracks per 10 MeV/c");
   fPzVsTOF = tfs->make<TH2F>("PzVsTOF","Pz Vs. TOF (All) ",160,0,1600,70,10,80);
   fPiMuPzVsTOF = tfs->make<TH2F>("PiMuPzVsTOF","PiMu Pz Vs. TOF",160,0,1600,70,10,80); 
}

void BeamlineTrigger::reconfigure(fhicl::ParameterSet const & p)
{
  // Implementation of optional member function here.
    fCalDataModuleLabel = p.get< std::string  >("CalDataModuleLabel");
    fnCalWireObjects    = p.get< double >("nCalWireObjects", 480.);
    // Implementation of optional member function here.
    fTOFModuleLabel 		= p.get< std::string >("TOFModuleLabel");
    fnTOFObjects			= p.get< double >("nTOFObjects", 1.0);
    fWCTrackModuleLabel           = p.get<  std::string  >("WCTrackModuleLabel");
    fminNumberWCTrack   		= p.get<     double    >("minNumberWCTrack", 1.0);
    fWCLowerBound                 = p.get< double >("WCLowerBound", 0.0);
    fWCUpperBound                 = p.get< double >("WCUpperBound", 2000.0);
}

DEFINE_ART_MODULE(BeamlineTrigger)
