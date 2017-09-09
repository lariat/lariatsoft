////////////////////////////////////////////////////////////////////////
// Class:       WCQualityFilter
// Module Type: filter
// File:        WCQualityFilter_module.cc
//
// Created by Daniel Smith, dansmith@bu.edu
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
#include <TH3F.h>
#include <TVector3.h>

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



class WCQualityFilter;

class WCQualityFilter : public art::EDFilter {
public:
  explicit WCQualityFilter(fhicl::ParameterSet const & p);
  // The destructor generated by the compiler is fine for classes
  // without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  WCQualityFilter(WCQualityFilter const &) = delete;
  WCQualityFilter(WCQualityFilter &&) = delete;
  WCQualityFilter & operator = (WCQualityFilter const &) = delete;
  WCQualityFilter & operator = (WCQualityFilter &&) = delete;

  void reconfigure(fhicl::ParameterSet const & p) override;

  bool insideImagPipe(std::vector<double> pos);

  std::vector<double> projToZ(std::vector<double> hit0, std::vector<double> hit1, double zpos);

  void beginJob() override;
  void endJob() override;
  // Required functions.
  bool filter(art::Event & e) override;
  
private:

  // Declare member data here.
  //---------- Filter Parameters ----------
  std::string fTOFModuleLabel;
  std::string fWCTrackLabel;

  bool bCreateMassPlots;
  double fMidPlaneCut;
  double fWC4ProjCut;       


  art::ServiceHandle<geo::Geometry> fGeo;  

  //---------- Histos ----------
  TH1F* hBeamlineMassOriginal;
  TH1F* hBeamlineMassAfterQuality;
  TH1F* hBeamlineMassRejected;

  TH2F* hTOFVsMomOriginal;
  TH2F* hTOFVsMomAfterQuality;
  TH2F* hTOFVsMomRejected;

  TH1F* hMomOriginal;
  TH1F* hMomAfterQuality;
  TH1F* hMomRejected;

  TH2F* hProjVsRealWC;
  TH1F* hRadDist;

  TH2F* hMidPlane;
  TH1F* hRadDistMidPlane;

  TH3F* hMomoVsProjXVsZ;
  TH3F* hMomoVsProjYVsZ;

};

// ---------------------- Begin Job ---------------------------
void WCQualityFilter::endJob()
{

}

// ---------------------- Begin Job ---------------------------
void WCQualityFilter::beginJob()
{

  // Implementation of optional member function here.
  art::ServiceHandle<art::TFileService> tfs;

  if(bCreateMassPlots) {
    hTOFVsMomOriginal = tfs->make<TH2F>("hTOFVsMomOriginal", "hTOFVsMomOriginal; WC Momentum [MeV/c];TOF [ns];", 1000, 0, 2000, 500, 0, 100); 
    hTOFVsMomAfterQuality     = tfs->make<TH2F>("hTOFvsMomAfterQualitys"    , "hTOFvsMomAfterQualitys;WC Momentum [MeV/c];TOF [ns];", 1000, 0, 2000, 500, 0, 100); 
    hTOFVsMomRejected = tfs->make<TH2F>("hTOFvsMomRejected"    , "hTOFvsMomRejected;WC Momentum [MeV/c];TOF [ns];", 1000, 0, 2000, 500, 0, 100); 

    hBeamlineMassOriginal = tfs->make<TH1F>("BeamlineMassOriginal" ,"Original BeamLine Mass" ,400,0,2000);  
    hBeamlineMassAfterQuality = tfs->make<TH1F>("BeamlineMassAfterQuality","BeamLine Mass After Quality Cuts",400,0,2000);  
    hBeamlineMassRejected = tfs->make<TH1F>("BeamlineMassRejected","BeamLine Mass Rejected Events",400,0,2000);  
  }

  hMomOriginal = tfs->make<TH1F>("MomOriginal" ,"Original Beamline Momentum" ,400,0,2000);  
  hMomAfterQuality = tfs->make<TH1F>("MomAfterQuality" ,"Beamline Momentum After Quality Cuts" ,400,0,2000);  
  hMomRejected = tfs->make<TH1F>("MomRejected" ,"Rejected beamline momentum" ,400,0,2000);  

  hProjVsRealWC = tfs->make<TH2F>("hProjVsRealWC","hProjVsRealWC", 100, -40.0, 40.0, 200, -40.0, 40.0);  
  hRadDist = tfs->make<TH1F>("hRadDist","hRadDist", 200, 0.0, 100.0);  

  hMidPlane = tfs->make<TH2F>("hMidPlaneDiff","hMidPlaneDiff", 100, -10.0, 10.0, 200, -10.0, 10.0);  
  hRadDistMidPlane  = tfs->make<TH1F>("hRadDistMid","hRadDistmid", 200, 0.0, 10.0);  

  hMomoVsProjXVsZ = tfs->make<TH3F>("hMomoVsProjXVsZ","hMomoVsProjXVsZ", 400, 0.0, 2000.0, 200, 0.0, 80.0, 80, -800.0, 0.0);
  hMomoVsProjYVsZ = tfs->make<TH3F>("hMomoVsProjYVsZ","hMomoVsProjYVsZ", 400, 0.0, 2000.0, 200, -20.0, 40.0, 80, -800.0, 0.0);

}


WCQualityFilter::WCQualityFilter(fhicl::ParameterSet const & p)
// :
// Initialize member data here.
{
  this->reconfigure(p);
  // Call appropriate produces<>() functions here.
}

bool WCQualityFilter::filter(art::Event & evt)
{
 
  art::Handle< std::vector<ldp::TOF> > TOFColHandle;
  std::vector<art::Ptr<ldp::TOF> > tof;  

  if(bCreateMassPlots) {
    // Getting the Time of Flight (TOF) Information
    if(!evt.getByLabel(fTOFModuleLabel,TOFColHandle)) { return false; }
    art::fill_ptr_vector(tof, TOFColHandle);   
    if(tof.size() != 1) { return false; }
    if(tof[0]->NTOF() != 1) { return false; }
  }

  // Getting the Momentum (WC) Information  
  art::Handle< std::vector<ldp::WCTrack> > wctrackHandle;
  std::vector<art::Ptr<ldp::WCTrack> > wctrack;

  if(!evt.getByLabel(fWCTrackLabel, wctrackHandle)){ return false; } 
  art::fill_ptr_vector(wctrack, wctrackHandle);

  int nGoodWC = 0;

  double reco_momo = -1.0;
  double reco_tof = -1.0;
  float mass = -1.0;

  for(size_t iWC = 0; iWC < wctrack.size(); iWC++) {
    if(wctrack[iWC]->NHits() != 8) { continue; } // require a hit on each wc plane

    //
    // Calculating the mass
    //

    float fDistanceTraveled = 6.652; 
    reco_momo = wctrack[iWC]->Momentum();
    hMomOriginal->Fill(reco_momo);

    if(bCreateMassPlots) {
      double tofObject[1];            // The TOF calculated (in ns) for this TOF object   
      tofObject[0] =  tof[0]->SingleTOF(0);

      reco_tof = tofObject[0];  
      mass = reco_momo*pow(reco_tof*0.299792458*0.299792458*reco_tof/(fDistanceTraveled*fDistanceTraveled) - 1 ,0.5);

      hTOFVsMomOriginal->Fill(reco_momo,reco_tof);
      hBeamlineMassOriginal->Fill(mass);
    }


    // 
    // Here starts my own code
    // 
    // First, I create an downstream beamline track 
    //   project it to WC4, and make sure that it matches the real hit
    // Then, I project the upstream and downstream beamline tracks to the midplane to make sure they match
    // Lastly, I check every 10 cm alone the entire beamline if the particle passes through matter
    //   calls up the function insideImagPipe(hit) which checks, using some geometry, if the hit passes into matter

    std::vector<double> wcHit0 {wctrack[iWC]->HitPosition(0, 0), wctrack[iWC]->HitPosition(0, 1), wctrack[iWC]->HitPosition(0, 2)};
    std::vector<double> wcHit1 {wctrack[iWC]->HitPosition(1, 0), wctrack[iWC]->HitPosition(1, 1), wctrack[iWC]->HitPosition(1, 2)};
    std::vector<double> wcHit2 {wctrack[iWC]->HitPosition(2, 0), wctrack[iWC]->HitPosition(2, 1), wctrack[iWC]->HitPosition(2, 2)};
    std::vector<double> wcHit3 {wctrack[iWC]->HitPosition(3, 0), wctrack[iWC]->HitPosition(3, 1), wctrack[iWC]->HitPosition(3, 2)};
      
    // 
    // Project to WC4 and check that it matches well
    // 

    // First, project downstream to midplane @ -437.97 (not using any angular corrections)
    std::vector<double> midUp = projToZ(wcHit0, wcHit1, -437.97);
    // Then use this point and WC3 to project up to WC4
    std::vector<double> ProjDown = projToZ(midUp, wcHit2, -95.0);
    // Requires some corrections because magnets are not the same
    ProjDown[0] -= tan(1.32 * TMath::Pi() / 180.0) * (-95.0 - -437.97);   // Corrections for the magnet irregularity 

    hProjVsRealWC->Fill(ProjDown[0] - wctrack[iWC]->HitPosition(3, 0), ProjDown[1] - wctrack[iWC]->HitPosition(3, 1));
    hRadDist->Fill(TMath::Sqrt(pow(ProjDown[0] - wctrack[iWC]->HitPosition(3, 0), 2.0) + 
			       pow(ProjDown[1] - wctrack[iWC]->HitPosition(3, 1), 2.0)));
  
    // The actual cut - make sure that hits match within 8 cm, decided from the hRadDist plot
    if(TMath::Sqrt(pow(ProjDown[0] - wctrack[iWC]->HitPosition(3, 0), 2.0) + 
		   pow(ProjDown[1] - wctrack[iWC]->HitPosition(3, 1), 2.0)) > fWC4ProjCut) { 
      if(bCreateMassPlots) {
	hTOFVsMomRejected->Fill(reco_momo,reco_tof);
	hBeamlineMassRejected->Fill(mass);     
      }	

      hMomRejected->Fill(reco_momo);
      continue;

    }


    // 
    // Project to midplane and check that upstream and downstream match well
    // 

    // upstream to midplane
    std::vector<double> resultUp = projToZ(wcHit0, wcHit1, -437.97);
    // downstream to midplane
    std::vector<double> resultDown = projToZ(wcHit2, wcHit3, -437.97);
    // Small (probably meaningless) correction
    resultDown[0] -= tan(1.32 * TMath::Pi() / 180.0) * (-339.57 - -437.97);   // Corrections for the magnet irregularity 

    hMidPlane->Fill(resultUp[0] - resultDown[0], resultUp[1] - resultDown[1]);
    hRadDistMidPlane->Fill(TMath::Sqrt(pow(resultUp[0] - resultDown[0] + 0.75, 2.0) + pow(resultUp[1] - resultDown[1], 2.0)));

    // The actual cut - make sure that tracks match within 3 cm, decided from the hRadDistMidPlane plot
    // NOTE : 0.75 added to the X coordinates to make sure that the 2D guassian is centered at (0, 0). 
    //   the reason it is off origin is because of 1) Magnets aren't identical or 2) Something else?
    if(TMath::Sqrt(pow(resultUp[0] - resultDown[0]  + 0.75, 2.0) + pow(resultUp[1] - resultDown[1], 2.0)) > fMidPlaneCut) { 
      if(bCreateMassPlots) {
	hTOFVsMomRejected->Fill(reco_momo,reco_tof);
	hBeamlineMassRejected->Fill(mass);     
      }
      hMomRejected->Fill(reco_momo);

      continue;

    }


    // 
    // Check every 10 cm if particle passed through matter
    // 

    bool bad = false;

    for(float iZ = -800.0; iZ < 0.0; iZ += 10.0) {
      std::vector<double> hito;

      // For every 10cm, project the beamline track to that point
      if(iZ < -437.97) {
	hito = projToZ(wcHit0, wcHit1, iZ);  // upstream
      } else { 
	// Assuming I can trust this point now that it has passed quality controls
	//std::vector<double> midUp = projToZ(wcHit0, wcHit1, -437.97);
	//hito = projToZ(midUp, wcHit2, iZ);
	//hito[0] -= tan(1.32 * TMath::Pi() / 180.0) * (iZ - -437.97);   // Corrections for the magnet irregularity 

	hito = projToZ(wcHit2, wcHit3, iZ); // downstream      
      }

      // Then, run it through the function that checks if it passed through any matter in the beamline
      //   if it did, reject the event
      if(!insideImagPipe(hito)) { 

	if(bCreateMassPlots) {
	  hTOFVsMomRejected->Fill(reco_momo,reco_tof);
	  hBeamlineMassRejected->Fill(mass);     
	}

	hMomRejected->Fill(reco_momo);

	bad = true;
	break;
      }

      hMomoVsProjXVsZ->Fill(wctrack[iWC]->Momentum(), hito[0], iZ);
      hMomoVsProjYVsZ->Fill(wctrack[iWC]->Momentum(), hito[1], iZ); 
    }

    if(bad) { continue; }
    nGoodWC += 1;
  }


  if(nGoodWC != 1) { return false; } // if we didn't find a single good WC, eject that stuff


  // 
  // This is the end of WC quality control. Thank you for visiting.
  //

  if(bCreateMassPlots) {
    hTOFVsMomAfterQuality->Fill(reco_momo,reco_tof);
    hBeamlineMassAfterQuality->Fill(mass);      
  }

  hMomAfterQuality->Fill(reco_momo);

  // If the event makes it all the way here, congrats you get to graduate now
  return true;
  
}


void WCQualityFilter::reconfigure(fhicl::ParameterSet const & p)
{                                                                                                   
  fTOFModuleLabel = p.get< std::string >("TOFModuleLabel");
  fWCTrackLabel   = p.get< std::string >("WCTrackLabel");

  bCreateMassPlots = p.get<bool>("CreateMassPlots", true);
  fMidPlaneCut = p.get<double>("MidPlaneCut", 3.0);
  fWC4ProjCut = p.get<double>("WC4ProjCut", 8.0);

}

bool WCQualityFilter::insideImagPipe(std::vector<double> pos)
{
 
  // 
  // This function uses the hit passed in (pos) and checks if it within
  //   the matter of 1) Magnet 1 2) Magnet or 3) the Downstream collimator
  // It does this using some complicated vector geometry. A description of the math involved can be found here: 
  // https://math.stackexchange.com/questions/1472049/check-if-a-point-is-inside-a-rectangular-shaped-area-3d
  //

  // Sanity check
  if(pos.size() != 3) { return false; }

  // All WC size {19.1, 19.1, 2.5}
  // WC1: center {107.95, 2.30, -686.13}, 13deg
  // WC2: center {71.47, 2.34, -536.36}, 13deg  
  // WC3: center {40.87, 2.043, -339.57}, 3deg
  // WC4: center {28.10, 0.30, -95.76}, 3 deg

  std::vector<double> c; // Center of box
  std::vector<double> vol; // volume of box
  double rot; // rotation of box

  // All numbers found in various gdml files. Sure there is a better way of doing this
  // If the particle is in the z-location of a given object, 
  //   that object's center point, volume and rotation are loaded up
  // These numbers should represent the Apeture of each of the objects (volume of apeture).
  //   if they don't, then my bad and better numbers are always welcome!

  if(-504.08 < pos[2]  and pos[2] < -444.98) { // Magnet 1
    c = {55.44, 2.29, -474.54};
    vol = {31.75, 14.22, 59.10};
    rot = 10.5 * TMath::Pi() / 180.0;

  } else if(-434.97 < pos[2] and pos[2] < -375.87) { // Magnet 2
   c = {46.51, 2.33, -405.42};
   vol = {31.75, 14.22, 59.10};
   rot = 5.5 * TMath::Pi() / 180.0;

  } else if(-260.0 < pos[2] and pos[2] < -170.0) { // Downstream collimator

    c = {34.49, 1.17, -217.67};

    // NOTE It was really hard to get an idea of where the apeture of the collimator was
    // More percise numbers are welcome if you can find them! 
    vol = {15.0, 15.0, 100.0}; //16.0, 243.81}; 

    rot = 3.0 * TMath::Pi() / 180.0;

  } else {
    return true; // Default true
  }  

  // From the center, volume and rotation loaded up, the fancy vector geomtry comes
  //  into play to make sure that our point z is within this rotated rectangular prism

  // Box corners 
  // {c0 + v0 / 2, TMath::Cos(10.5) (c1 + v1 / 2), c2 - v2 / 2}
  // {c0 + v0 / 2, TMath::Cos(10.5) (c1 - v1 / 2), c2 - v2 / 2}
  // {c0 - v0 / 2, TMath::Cos(10.5) (c1 + v1 / 2), c2 - v2 / 2}
  // {c0 - v0 / 2, TMath::Cos(10.5) (c1 - v1 / 2), c2 - v2 / 2}

  // {c0 + v0 / 2, TMath::Cos(10.5) (c1 + v1 / 2), c2 + v2 / 2}
  // {c0 + v0 / 2, TMath::Cos(10.5) (c1 - v1 / 2), c2 + v2 / 2}
  // {c0 - v0 / 2, TMath::Cos(10.5) (c1 + v1 / 2), c2 + v2 / 2}
  // {c0 - v0 / 2, TMath::Cos(10.5) (c1 - v1 / 2), c2 + v2 / 2}

  TVector3 P1, P2, P4, P5; // Points needed to find if a third point is within the box
  P1.SetXYZ(TMath::Cos(rot)*(c[0] + vol[0] / 2.0), (c[1] + vol[1] / 2.0), c[2] - vol[2] / 2.0);
  P2.SetXYZ(TMath::Cos(rot)*(c[0] + vol[0] / 2.0), (c[1] + vol[1] / 2.0), c[2] + vol[2] / 2.0);
  P4.SetXYZ(TMath::Cos(rot)*(c[0] - vol[0] / 2.0), (c[1] + vol[1] / 2.0), c[2] - vol[2] / 2.0);
  P5.SetXYZ(TMath::Cos(rot)*(c[0] + vol[0] / 2.0), (c[1] - vol[1] / 2.0), c[2] - vol[2] / 2.0);

  TVector3 u = (P1 - P4).Cross(P1 - P5);
  TVector3 v = (P1 - P2).Cross(P1 - P5);
  TVector3 w = (P1 - P2).Cross(P1 - P4);

  TVector3 x(pos[0], pos[1], pos[2]);

  bool verbose = false; // My super cheap verbose variable 

  if(verbose) {
    std::cout << 
      "u.x " << u.Dot(x) << " " << 
      "v.x " << v.Dot(x) << " " << 
      "w.x " << w.Dot(x) << " " << std::endl;
    std::cout << 
      "u.p1 " << u.Dot(P1) << " " << 
      "v.p1 " << v.Dot(P1) << " " << 
      "w.p1 " << w.Dot(P1) << " " << std::endl;
    std::cout << 
      "u.p2 " << u.Dot(P2) << " " << 
      "v.p4 " << v.Dot(P4) << " " << 
      "w.p5 " << w.Dot(P5) << " " << std::endl;
  }

  // If this long statement is true, then the point is within the prism. 
  // if this long statement is false, then the point is Outside of the prism

  return (((u.Dot(P1) < u.Dot(x) and u.Dot(x) < u.Dot(P2)) or (u.Dot(P1) > u.Dot(x) and u.Dot(x) > u.Dot(P2))) and
  	  ((v.Dot(P1) < v.Dot(x) and v.Dot(x) < v.Dot(P4)) or (v.Dot(P1) > v.Dot(x) and v.Dot(x) > v.Dot(P4))) and
  	  ((w.Dot(P1) < w.Dot(x) and w.Dot(x) < w.Dot(P5)) or (w.Dot(P1) > w.Dot(x) and w.Dot(x) > w.Dot(P5))));

}


std::vector<double> WCQualityFilter::projToZ(std::vector<double> hit0, std::vector<double> hit1, double zpos) {

  //
  // This powerful little piece of code is what is used to find the point at z = zpos along 
  //   the line created by hit0 and hit1. Uses the parameterized vector form of a line
  //
  // <x, y, z> = <sx, sy, sz> * t + <startx, starty, startz>
  // sx, sy, sz are all slopes
  //
  // (z - startz) / sz = t
  // x = sx * t + startx
  // y = sx * t + srarty
  //
      
  double sx = hit1[0] - hit0[0];
  double sy = hit1[1] - hit0[1];
  double sz = hit1[2] - hit0[2];

  double t = (zpos - hit0[2]) / sz;

  std::vector<double> resulto {sx * t + hit0[0], sy * t + hit0[1], zpos};

  return resulto;

}


DEFINE_ART_MODULE(WCQualityFilter)


