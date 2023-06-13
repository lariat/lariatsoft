////////////////////////////////////////////////////////////////////////
// Class:       ParticleDump
// Module Type: analyzer
// File:        ParticleDump_module.cc
//
// Generated at Tue Jul 14 11:21:46 2015 by Roberto Acciarri using artmod
// from cetpkgsupport v1_08_06.
////////////////////////////////////////////////////////////////////////

// ##########################
// ### Framework includes ###
// ##########################
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/PtrVector.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"
#include "art_root_io/TFileDirectory.h"
#include "canvas/Persistency/Common/FindOneP.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
//#include "cetlib/maybe_ref.h"

// ########################
// ### LArSoft includes ###
// ########################
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "larcoreobj/SimpleTypesAndConstants/RawTypes.h" // raw::ChannelID_t
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/CryostatGeo.h"
#include "larcorealg/Geometry/TPCGeo.h"
#include "larcorealg/Geometry/PlaneGeo.h"
#include "larcorealg/Geometry/WireGeo.h"
#include "lardataobj/RecoBase/Wire.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/TrackHitMeta.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardata/ArtDataHelper/TrackUtils.h" // lar::util::TrackPitchInView()
#include "lardata/DetectorInfoServices/LArPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardata/Utilities/AssociationUtil.h"

//#include "RawData/ExternalTrigger.h"
#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RawData/raw.h"
#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "larevt/Filters/ChannelFilter.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "lardataobj/AnalysisBase/ParticleID.h"
#include "larreco/RecoAlg/TrackMomentumCalculator.h"
#include "LArIATRecoAlg/BeamlineMassAlg.h"
#include "LArIATDataProducts/WCTrack.h"
#include "LArIATDataProducts/TOF.h"
#include "LArIATDataProducts/AGCounter.h"
#include "RawDataUtilities/TriggerDigitUtility.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/EndPoint2D.h"
#include "lardataobj/MCBase/MCShower.h"
#include "lardataobj/MCBase/MCStep.h"
#include "larreco/Calorimetry/CalorimetryAlg.h"

// #####################
// ### ROOT includes ###
// #####################
#include "TComplex.h"
#include "TFile.h"
#include "TH2D.h"
#include "TF1.h"
#include "TTree.h"
#include "TTimeStamp.h"

namespace lariat
{
  class ParticleDump;
}

class lariat::ParticleDump : public art::EDAnalyzer
{
public:
  explicit ParticleDump(fhicl::ParameterSet const & p);
  virtual ~ParticleDump();

  // Required functions.
  void analyze(art::Event const & e) override;

  // Selected optional functions.
  //void beginJob() override;
  //void reconfigure(fhicl::ParameterSet const & p);

private:
 
  // === Get path length of MCParticle in TPC ===
  float TrajLengthInTpcAV(const simb::MCParticle*);
  bool  IsPointInTpcAV(const simb::MCParticle*, int);
  bool  IsPointInTpcAV(TVector3&);

};


lariat::ParticleDump::ParticleDump(fhicl::ParameterSet const & pset)
  : EDAnalyzer(pset)
{
  //this->reconfigure(pset);
}

lariat::ParticleDump::~ParticleDump()
{
}

//void lariat::ParticleDump::reconfigure(fhicl::ParameterSet const & pset)
//{
//}

void lariat::ParticleDump::analyze(art::Event const & evt)
{
  
  bool isData = (bool)evt.isRealData();
  if( isData ) return;

 printf("ParticleDump: looking at event %i\n",(int)evt.event());

  // #######################################
  // ### Get potentially useful services ###
  // #######################################
  art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;
  const sim::ParticleList& plist = pi_serv->ParticleList();
  
  std::cout<<"Looping through particles: "<<plist.size()<<"\n";

  for(size_t p = 0; p < plist.size(); ++p) {
    auto& particle = plist.Particle(p);
    
    int     last        = particle->NumberTrajectoryPoints() - 1;
    TVector3 loc0       = particle->Position(0).Vect();
    TVector3 locf       = particle->Position(last).Vect();
    //if( !IsPointInTpcAV(loc0) && !IsPointInTpcAV(locf) ) return;
    
    int     trackId     = particle->TrackId();
    int     mother      = particle->Mother();
    int     pdg         = particle->PdgCode();
    int     ndaughters  = particle->NumberDaughters();
    double  dL          = (particle->Position(0).Vect()-particle->Position(last).Vect()).Mag();
    float   KE0         = 1e3*(particle->E(0)-particle->Mass());
    float   T0          = 1e-3*particle->T(0);
    std::string proc    = particle->Process().c_str();
    std::string endProc = particle->EndProcess().c_str();
    float fracPx        = particle->Px()/particle->P();
    float fracPy        = particle->Py()/particle->P();
    float fracPz        = particle->Pz()/particle->P();
    
    printf("  %5i PDG: %7i, dL: %6.1fcm, XYZ: (%7.1f,%7.1f,%7.1f),  dir: (%5.2f, %5.2f, %5.2f), KE0: %7.3f, T0: %7.2f us, moth: %5i, %16.16s -->%16.16s, nD: %i\n",
      trackId,
      pdg,
      dL,
      loc0.X(),loc0.Y(),loc0.Z(),
      fracPx, fracPy, fracPz,
      KE0,
      T0, 
      mother,
      proc.c_str(),endProc.c_str(),
      ndaughters
    );
 

  }
  

}


//void lariat::ParticleDump::beginJob()
//{
//
//}

bool lariat::ParticleDump::IsPointInTpcAV(TVector3& v){
  // Get active volume boundary.
  art::ServiceHandle<geo::Geometry> geom;
  float   xmin = 0. +1e-8;
  float   xmax = 2.0 * geom->DetHalfWidth() - 1e-8;
  float   ymin = -geom->DetHalfHeight() + 1e-8;
  float   ymax = geom->DetHalfHeight()  - 1e-8;
  float   zmin = 0. +1e-8;
  float   zmax = geom->DetLength() - 1e-8;
  if(   (v.X()>=xmin && v.X()<=xmax)
    &&  (v.Y()>=ymin && v.Y()<=ymax)
    &&  (v.Z()>=zmin && v.Z()<=zmax) ) 
    return true;
  else 
    return false;
}

bool lariat::ParticleDump::IsPointInTpcAV(const simb::MCParticle* part, int i){
  TVector3 p(part->Vx(i),part->Vy(i), part->Vz(i));
  return IsPointInTpcAV(p);
}


float lariat::ParticleDump::TrajLengthInTpcAV(const simb::MCParticle* part){
  float L = 0.; 
  int nPts = part->NumberTrajectoryPoints();
  for(int i = 1; i < nPts; ++i) {
    // if both this and previous point are in AV, add to total traj length
    TVector3 pA(part->Vx(i-1),part->Vy(i-1),part->Vz(i-1));
    TVector3 pB(part->Vx(i),part->Vy(i),part->Vz(i));
    if( IsPointInTpcAV(pA) && IsPointInTpcAV(pB) ) {
      L += (pB-pA).Mag();
    }
  }
  return L;
}


DEFINE_ART_MODULE(lariat::ParticleDump)
