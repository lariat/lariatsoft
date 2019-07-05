////////////////////////////////////////////////////////////////////////
// Class:       WCQualityProducer
// Module Type: filter
// File:        WCQualityProducer_module.cc
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
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Track.h"

#include "canvas/Persistency/Common/FindOneP.h" 
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "canvas/Persistency/Common/Ptr.h" 
#include "canvas/Persistency/Common/PtrVector.h" 


#include "larcore/Geometry/Geometry.h"
#include "lardata/Utilities/AssociationUtil.h"



class WCQualityProducer;

class WCQualityProducer : public art::EDProducer {
public:
  explicit WCQualityProducer(fhicl::ParameterSet const & p);
  // The destructor generated by the compiler is fine for classes
  // without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  WCQualityProducer(WCQualityProducer const &) = delete;
  WCQualityProducer(WCQualityProducer &&) = delete;
  WCQualityProducer & operator = (WCQualityProducer const &) = delete;
  WCQualityProducer & operator = (WCQualityProducer &&) = delete;

  void reconfigure(fhicl::ParameterSet const & p) override;

  bool CheckUpstreamMagnetAperture(std::vector<double> hit1, std::vector<double> hit2);
  bool CheckDownstreamMagnetAperture(std::vector<double> hit1, std::vector<double> hit2);
  bool CheckDownstreamCollimatorAperture(std::vector<double> hit1, std::vector<double> hit2);
  void projToZ(double hit0[3], double hit1[3], double (&result)[3] , double zpos);
  std::vector<bool> CheckUSMagApertures(double hit0[3], double hit1[3]);
  std::vector<bool> CheckDSMagApertures(double hit0[3], double hit1[3]);
  std::vector<bool> CheckDSColApertures(double hit0[3], double hit1[3]);
  void beginJob() override;
  void endJob() override;
  // Required functions.
  void produce(art::Event & evt) override;
  
private:

  // Declare member data here.
  //---------- Filter Parameters ----------
  std::string fTOFModuleLabel;
  std::string fWCTrackLabel;

  bool UseMidplaneCut;
  bool UseCollimatorCut;
  bool ApplyMassCut;
  double fMidPlaneCut;
  double fMassLowerLimit;
  double fMassUpperLimit;      
  bool IsThisMC;
  bool ApplyWC4Fiducialization;

  double fDataMidDiffXOffset, fDataMidDiffYOffset, fDataXCut, fDataYCut;
  double mid1Hit[3]={-9999,-9999,-9999};
  double mid2Hit[3]={-9999,-9999,-9999};
  double fWC4FiducialCutLow, fWC4FiducialCutHigh;
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

  TH2F* hMidPlaneAfterQuality;
  TH1F* hRadDistMidPlaneAfterQuality;
  TH3F* hMomoVsProjXVsZ;
  TH3F* hMomoVsProjYVsZ;
  
  
  bool ExtrapolateToMP;  //Using WC1, WC2 project to midplane. Use WC3, WC4, project to Midplane. Are those points close? Used with fMidplaneCut.

  //Bools that track hit apertures.
  bool Magnet1ApertureCheck;
  bool Magnet2ApertureCheck;
  bool DSColApertureCheck;
  bool WC4FidCheck;
  bool KeepTheEvent; //Depending on which Checks you want to use, the final boolean that combines these checks to decide if the event is good.
  bool Mag1USPassX, Mag1USPassY, Mag1DSPassY, Mag1Pass, Mag2USPassY, Mag2DSPassX, Mag2DSPassY, Mag2Pass, DSColPassX, DSColPassY, DSColPass, AperturePass;  //All the aperture check booleans.
  // === Storing Time of Flight information === 
  
    
  // I HATE TO DO THIS: I'M HARDCODING THE POSITIONS OF THE MAGNET CENTERS 
  // CAUSE I DON'T KNOW HOW TO FETCH non-AuxDet pieces in the gdml
  double NDB1_Center[3] = { 55.801, 5.048, -472.218};
  double NDB2_Center[3] = { 47.114, 4.726, -403.146};
  double Magnets_Mid[3] = { (NDB1_Center[0]+NDB2_Center[0])/2., (NDB1_Center[1]+NDB2_Center[1])/2., (NDB1_Center[2]+NDB2_Center[2])/2.};
  
  
  //I hate to do this too, but I have to hard code the corners of the magnetic aperatures and the center of the aperature faces.
  double NDB1_FFace_center[3]={61.186,    5.048,    -501.273};  //XYZ of mag1 upstream face
  double NDB1_BFace_center[3]={50.416,    5.048,    -443.163};  //XYZ of mag1 downstream face
  
  double NDB1_FFace_XZ_slope=5.3955; //using the center and the expected rotation around the Y axis (-10.5deg, nominally), find the equation of the line of the face, in the XZ plane
  double NDB1_FFace_XZ_intercept=2765.814;

  double NDB1_BFace_XZ_slope=5.3955;
  double NDB1_BFace_XZ_intercept=2441.509;  
  
 
    //magnet 2
  double NDB2_FFace_center[3]={49.946,    4.726,    -432.56};  //XYZ of mag2 upstream face
  double NDB2_BFace_center[3]={44.282,    4.726,    -373.732};  //XYZ of mag2 downstream face
  
  double NDB2_FFace_XZ_slope=10.3854;  //Same logic as Magnet 1, but with a -5.5deg rotation.
  double NDB2_FFace_XZ_intercept=4542.253;
  
  double NDB2_BFace_XZ_slope=10.3854;
  double NDB2_BFace_XZ_intercept=3925.637;
  
  
    //DS Collimator
  // Because we expect WC3/4 to be at the same rotation as the DSCol, the planes defining the WCs and the faces of the col are all parallel. 
  //we only need to check at one face, so I use the front face. Yay simplified geometry!
  double DSCol_FFace_center[3]={38.404,    3.323,    -293.256};
  
  double DSCol_FFace_XZ_slope=19.0811;  //Same logic, but with a 3deg rotation.
  double DSCol_FFace_XZ_intercept=5634.062;  
  
   //The position of each corner of the aperture. F or B for Front or Back face. top or bottom (relative to gravity), left or right relative to beam (left = TPC cathode side)
  double NDB1_F_top_left[3]=     {71.879,    12.158,    -449.291};
  double NDB1_F_top_right[3]=    {50.493,    12.158,    -503.255};
  double NDB1_F_bottom_left[3]=  {71.879,    -2.062,    -499.291};
  double NDB1_F_bottom_right[3]= {50.493,    -2.062,    -503.255};
  
  double NDB1_B_top_left[3]=     {61.109,    12.158,    -441.181};    
  double NDB1_B_top_right[3]=    {39.723,    12.158,    -445.145};
  double NDB1_B_bottom_left[3]=  {61.109,    -2.062,    -441.181};
  double NDB1_B_bottom_right[3]= {39.723,    -2.062,    -445.145};
  

  double NDB2_F_top_left[3]=     {60.771,    11.836,    -431.518};
  double NDB2_F_top_right[3]=    {39.121,    11.836,    -433.602};
  double NDB2_F_bottom_left[3]=  {60.771,    -2.384,    -431.518};
  double NDB2_F_bottom_right[3]= {39.121,    -2.384,    -433.602};
  
  double NDB2_B_top_left[3]=     {55.107,    11.836,    -372.690};   
  double NDB2_B_top_right[3]=    {33.457,    11.833,    -374.774};
  double NDB2_B_bottom_left[3]=  {55.107,    -2.384,    -372.690};
  double NDB2_B_bottom_right[3]= {33.457,    -2.384,    -374.774};
  //DS Collimator
  // Because we expect WC3/4 to be at the same rotation as the DSCol, the planes defining the WCs and the faces of the col are all parallel. 
  //we only need to check at one face, so I only need the front face. Yay simplified geometry!

  
  double DSCol_F_top_left[3]=     {45.096,    11.273,    -292.905};
  double DSCol_F_top_right[3]=    {31.712,    11.273,    -293.607};
  double DSCol_F_bottom_left[3]=  {45.096,    -4.627,    -292.905};
  double DSCol_F_bottom_right[3]= {31.712,    -4.627,    -293.607}; 
/*   //For each "collimator", the bounds of the face of both aperatures [xlow_frontface, xhigh_frontface, xlow_backface, xhigh_backface], similarly for y. In cm, in TPC coordinates. Taken from survey.
  double xboundMagnet1[4]={45.74, 75.52, 35.09, 64.87};
  double yboundMagnet1[4]={-13.12, 13.59, -13.16, 13.55};
  
  double xboundMagnet2[4]={34.88, 65.12, 27.90, 58.15};
  double yboundMagnet2[4]={-13.10, 13.61, -13.08, 13.63};
 
  double xboundDSCol[4]={30.33, 45.42, 26.65, 41.73};
  double yboundDSCol[4]={-15.70, 14.91, -15.53, 15.08};
  // Z Position of the center of the aperatures of each collimator, found by taking the average of the z bounds of the aperature. [zcent_US, zcent_DS]
  double zcentMagnet1[2] = { (-501.95-494.98)/2, (-449.49-456.46)/2};
  double zcentMagnet2[2] = { (-432.04-427.50)/2, (-381.27-385.81)/2};
  double zcentDSCol[2]   = { (-296.67-297.36)/2, (-205.94-206.63)/2};
  double Keepcount=0; */
};

// ---------------------- Begin Job ---------------------------
void WCQualityProducer::endJob()
{

}

// ---------------------- Begin Job ---------------------------
void WCQualityProducer::beginJob()
{

  // Implementation of optional member function here.
  art::ServiceHandle<art::TFileService> tfs;

  if(ApplyMassCut &&!IsThisMC) {
    hTOFVsMomOriginal = tfs->make<TH2F>("hTOFVsMomOriginal", "hTOFVsMomOriginal; WC Momentum [MeV/c];TOF [ns];", 1000, 0, 2000, 500, 0, 100); 
    hTOFVsMomAfterQuality     = tfs->make<TH2F>("hTOFvsMomAfterQualitys"    , "hTOFvsMomAfterQualitys;WC Momentum [MeV/c];TOF [ns];", 1000, 0, 2000, 500, 0, 100); 
    hTOFVsMomRejected = tfs->make<TH2F>("hTOFvsMomRejected"    , "hTOFvsMomRejected;WC Momentum [MeV/c];TOF [ns];", 1000, 0, 2000, 500, 0, 100); 

    hBeamlineMassOriginal = tfs->make<TH1F>("BeamlineMassOriginal" ,"Original BeamLine Mass" ,800,-2000,2000);  
    hBeamlineMassAfterQuality = tfs->make<TH1F>("BeamlineMassAfterQuality","BeamLine Mass After Quality Cuts",800,-2000,2000);  
    hBeamlineMassRejected = tfs->make<TH1F>("BeamlineMassRejected","BeamLine Mass Rejected Events",800,-2000,2000);  
  }

  hMomOriginal = tfs->make<TH1F>("MomOriginal" ,"Original Beamline Momentum" ,100,0,2000);  
  hMomAfterQuality = tfs->make<TH1F>("MomAfterQuality" ,"Beamline Momentum After Quality Cuts" ,100,0,2000);  
  hMomRejected = tfs->make<TH1F>("MomRejected" ,"Rejected beamline momentum" ,100,0,2000);  

  hProjVsRealWC = tfs->make<TH2F>("hProjVsRealWC","hProjVsRealWC", 100, -40.0, 40.0, 200, -40.0, 40.0);  
  hRadDist = tfs->make<TH1F>("hRadDist","hRadDist", 200, 0.0, 100.0);  

  hMidPlane = tfs->make<TH2F>("hMidPlaneDiff","hMidPlaneDiff", 200, -50.0, 50.0, 200, -50.0, 50.0);  
  hRadDistMidPlane  = tfs->make<TH1F>("hRadDistMid","hRadDistmid", 200, 0.0, 100.0);  
  hMidPlaneAfterQuality = tfs->make<TH2F>("hMidPlaneDiffAfterQuality","hMidPlaneDiff", 200, -50.0, 50.0, 200, -50.0, 50.0);  
  hRadDistMidPlaneAfterQuality  = tfs->make<TH1F>("hRadDistMidAfterQuality","hRadDistmid", 200, 0.0, 100.0); 
  hMomoVsProjXVsZ = tfs->make<TH3F>("hMomoVsProjXVsZ","hMomoVsProjXVsZ", 400, 0.0, 2000.0, 200, 0.0, 80.0, 80, -800.0, 0.0);
  hMomoVsProjYVsZ = tfs->make<TH3F>("hMomoVsProjYVsZ","hMomoVsProjYVsZ", 400, 0.0, 2000.0, 200, -20.0, 40.0, 80, -800.0, 0.0);

}


WCQualityProducer::WCQualityProducer(fhicl::ParameterSet const & p)
// :
// Initialize member data here.
{
  this->reconfigure(p);

  // Call appropriate produces<>() functions here.

  // produce recob::PFParticle objects
  produces< std::vector< recob::PFParticle > >();
}

void WCQualityProducer::produce(art::Event & evt)
{

  
  // This I can fetch from the Geo: much better!
  double USTOF_Center[3];
  double WC1_Center[3];
  double WC2_Center[3];
  double WC3_Center[3];
  double WC4_Center[3];
  double DSTOF_Center[3];
    

  for( size_t iDet = 0; iDet < fGeo->NAuxDets() ; ++iDet ){
    geo::AuxDetGeo const& anAuxDetGeo = fGeo->AuxDet(iDet);
    std::string detName = anAuxDetGeo.Name();
    if( detName == "volAuxDetTOFUS")        anAuxDetGeo.GetCenter(USTOF_Center); 
    if( detName == "volAuxDetSensitiveWC1") anAuxDetGeo.GetCenter(WC1_Center); 
    if( detName == "volAuxDetSensitiveWC2") anAuxDetGeo.GetCenter(WC2_Center); 
    if( detName == "volAuxDetSensitiveWC3") anAuxDetGeo.GetCenter(WC3_Center); 
    if( detName == "volAuxDetSensitiveWC4") anAuxDetGeo.GetCenter(WC4_Center); 
    if( detName == "volAuxDetTOFDS")        anAuxDetGeo.GetCenter(DSTOF_Center); 
  }

 double  centerTOFsGeoDist = TMath::Sqrt((USTOF_Center[0] - DSTOF_Center[0])*(USTOF_Center[0] - DSTOF_Center[0]) + 
				  (USTOF_Center[1] - DSTOF_Center[1])*(USTOF_Center[1] - DSTOF_Center[1]) + 
				  (USTOF_Center[2] - DSTOF_Center[2])*(USTOF_Center[2] - DSTOF_Center[2]) ) ;
				  
				    
  ExtrapolateToMP=true;  
  

  Magnet1ApertureCheck=true;
  Magnet2ApertureCheck=true;
  DSColApertureCheck=true;
  KeepTheEvent=true;
  WC4FidCheck=true;
  art::Handle< std::vector<ldp::TOF> > TOFColHandle;
  std::vector<art::Ptr<ldp::TOF> > tof;  

  if(ApplyMassCut && !IsThisMC) {
    // Getting the Time of Flight (TOF) Information
    if(!evt.getByLabel(fTOFModuleLabel,TOFColHandle)) { KeepTheEvent=false; }
    art::fill_ptr_vector(tof, TOFColHandle);   
    if(tof.size() != 1) { KeepTheEvent=false; }
    else if(tof[0]->NTOF() != 1) { KeepTheEvent=false; }
  }

  // Getting the Momentum (WC) Information  
  art::Handle< std::vector<ldp::WCTrack> > wctrackHandle;
  std::vector<art::Ptr<ldp::WCTrack> > wctrack;

  if(!evt.getByLabel(fWCTrackLabel, wctrackHandle)){ KeepTheEvent=false; } 
  art::fill_ptr_vector(wctrack, wctrackHandle);
  
  //int nGoodWC = 0;
  double reco_momo = -1.0;
  double reco_tof = -1.0;
  float mass = -1.0;
  if(KeepTheEvent && wctrack.size()==0){KeepTheEvent=false;}
  if(KeepTheEvent){  
  for(size_t iWC = 0; iWC < wctrack.size(); iWC++) {
    if(KeepTheEvent && wctrack[iWC]->NHits() != 8) { KeepTheEvent=false; } // require a hit on each wc plane

    //
    // Calculating the mass
    //
   
    
    reco_momo = wctrack[iWC]->Momentum();
    hMomOriginal->Fill(reco_momo);
    if(KeepTheEvent && ApplyMassCut && !IsThisMC){    //Cant check TOF until I know it exists. So breaking into two if statements.
      if(tof[0]->NTOF()==1){
        double tofObject[1];            // The TOF calculated (in ns) for this TOF object   
        tofObject[0] =  tof[0]->SingleTOF(0);

        reco_tof = tofObject[0];  
	
	double fDistanceTraveled = centerTOFsGeoDist/100; //Convert from cm to m.
	double radical = reco_tof*0.299792458*0.299792458*reco_tof/(fDistanceTraveled*fDistanceTraveled) - 1;
	if (radical<0)
	{
        mass = -reco_momo*pow(-radical,0.5);
	}
	
	if (radical>0)
	{
        mass = reco_momo*pow(radical,0.5);
	}
	if(KeepTheEvent && ApplyMassCut && !IsThisMC) {
          hTOFVsMomOriginal->Fill(reco_momo,reco_tof);
          hBeamlineMassOriginal->Fill(mass);
        }
        if(KeepTheEvent && ApplyMassCut && !IsThisMC && (mass<=fMassLowerLimit || mass>=fMassUpperLimit)){KeepTheEvent=false;}  
      }
    }
    // 
    // Here starts my own code
    // 
    // First, I create an downstream beamline track 


    double WC1Hits[3]= {wctrack[iWC]->HitPosition(0, 0), wctrack[iWC]->HitPosition(0, 1), wctrack[iWC]->HitPosition(0, 2)};
    double WC2Hits[3]= {wctrack[iWC]->HitPosition(1, 0), wctrack[iWC]->HitPosition(1, 1), wctrack[iWC]->HitPosition(1, 2)};
    double WC3Hits[3]= {wctrack[iWC]->HitPosition(2, 0), wctrack[iWC]->HitPosition(2, 1), wctrack[iWC]->HitPosition(2, 2)};
    double WC4Hits[3]= {wctrack[iWC]->HitPosition(3, 0), wctrack[iWC]->HitPosition(3, 1), wctrack[iWC]->HitPosition(3, 2)};
    std::cout<<WC4_Center[0]-fWC4FiducialCutLow<<" "<<WC4_Center[0]+fWC4FiducialCutHigh<<std::endl;
    if (WC4Hits[0]<(WC4_Center[0]-fWC4FiducialCutLow) || WC4Hits[0]>(WC4_Center[0]+fWC4FiducialCutHigh)){WC4FidCheck=false;}
     


    projToZ(WC1Hits, WC2Hits, mid1Hit      , Magnets_Mid [2]);
    projToZ(WC3Hits, WC4Hits, mid2Hit      , Magnets_Mid [2]);


    hMidPlane->Fill(10*(mid1Hit[0] - mid2Hit[0]), 10*(mid1Hit[1] - mid2Hit[1])); //Fill this in mm, because I expect an offset of 3mm and 0mm respectively
    hRadDistMidPlane->Fill(TMath::Sqrt(pow(10*(mid1Hit[0] - mid2Hit[0]), 2.0) + pow(10*(mid1Hit[1] -    mid2Hit[1]), 2.0)));

    if((fabs((mid1Hit[0] - mid2Hit[0])-fDataMidDiffXOffset))>fDataXCut || (fabs((mid1Hit[1] - mid2Hit[1])-fDataMidDiffYOffset))>fDataYCut)
      ExtrapolateToMP=false;

    
  
  
  std::vector<bool> Mag1Bools=CheckUSMagApertures(WC1Hits, WC2Hits);  //{USX, USY, DSY, Total}
  std::vector<bool> Mag2Bools=CheckDSMagApertures(WC3Hits, WC4Hits);  //{USY, DSX, DSY, Total}
  std::vector<bool> DSColBools=CheckDSColApertures(WC3Hits, WC4Hits); //{X,Y,total}
   
   
  
   
  Mag1USPassX=Mag1Bools[0];
  Mag1USPassY=Mag1Bools[1]; 
  Mag1DSPassY=Mag1Bools[2]; 
  Mag1Pass=Mag1Bools[3]; 
  
  Mag2USPassY=Mag2Bools[0]; 
  Mag2DSPassX=Mag2Bools[1]; 
  Mag2DSPassY=Mag2Bools[2]; 
  Mag2Pass=Mag2Bools[3]; 
  
  DSColPassX=DSColBools[0]; 
  DSColPassY=DSColBools[1]; 
  DSColPass=DSColBools[2]; 
  }
  if (Mag1Pass && Mag2Pass && DSColPass){AperturePass=true;}
  else {AperturePass=false;}

}
  if(UseMidplaneCut && !ExtrapolateToMP){KeepTheEvent=false;}
  if(UseCollimatorCut && !AperturePass ){KeepTheEvent=false;}
  if(ApplyWC4Fiducialization && !WC4FidCheck){ std::cout<<"Rejecting: "<<wctrack[0]->HitPosition(3, 0)<<std::endl;KeepTheEvent=false;}
  if(KeepTheEvent){
    if(ApplyMassCut && !IsThisMC) {
      hTOFVsMomAfterQuality->Fill(reco_momo,reco_tof);
      hBeamlineMassAfterQuality->Fill(mass);    
    }
          hMidPlaneAfterQuality->Fill(10*(mid1Hit[0] - mid2Hit[0]), 10*(mid1Hit[1] - mid2Hit[1]));
      hRadDistMidPlaneAfterQuality->Fill(TMath::Sqrt(pow(10*(mid1Hit[0] - mid2Hit[0]), 2.0) + pow(10*(mid1Hit[1] -    mid2Hit[1]), 2.0)));  
    hMomAfterQuality->Fill(reco_momo);
  }  
  if(!KeepTheEvent)
  {
    if (ApplyMassCut && !IsThisMC)
    {
     	hTOFVsMomRejected->Fill(reco_momo,reco_tof);
	hBeamlineMassRejected->Fill(mass);    
    }
    hMomRejected->Fill(reco_momo);
  }

  // apply WC quality flag if this is MC
  if (IsThisMC) KeepTheEvent = true;

  // If the event makes it all the way here, congrats you get to graduate now

  //return KeepTheEvent;

  //-------------------------------------------------------------------
  // point to a collection of PFParticle objects
  //-------------------------------------------------------------------
  std::vector< recob::PFParticle > pfparticle_vector;

  std::vector< size_t > pfparticle_daughter_indices;

  if (KeepTheEvent)
  {
    pfparticle_vector.emplace_back(211, 0, 0, pfparticle_daughter_indices);
  }
  //else
  //{
  //  pfparticle_vector.emplace_back(-1, 0, 0, pfparticle_daughter_indices);
  //}

  // convert PFParticle vector to unique_ptr
  std::unique_ptr< std::vector< recob::PFParticle > >
      pfparticle_collection(
          new std::vector< recob::PFParticle >(std::move(pfparticle_vector)));

  // put collections into event
  evt.put(std::move(pfparticle_collection));

  return;

}

//==================================================================================================
void WCQualityProducer::reconfigure(fhicl::ParameterSet const & p)
{
  fTOFModuleLabel = p.get< std::string >("TOFModuleLabel");
  fWCTrackLabel   = p.get< std::string >("WCTrackLabel");

  fMidPlaneCut = p.get<double>("MidPlaneCut", 3.0);
  fDataXCut = p.get<double>("DataMPXCut", 1.5);
  fDataYCut = p.get<double>("DataMPYCut", 1.5);
  fDataMidDiffXOffset = p.get<double>("DataMPXOffset", 0.3);
  fDataMidDiffYOffset = p.get<double>("DataMPYOffset", 0);
  UseMidplaneCut = p.get<bool>("ApplyMidplaneCut", true);
  UseCollimatorCut = p.get<bool>("ApplyCollimatorCut", true);
  ApplyMassCut=p.get<bool>("ApplyMassCut",false);
  fMassLowerLimit=p.get<double>("LowerMassLimit",0);
  fMassUpperLimit=p.get<double>("UpperMassLimit",2000);
  IsThisMC = p.get< bool >("IsThisMC", false);
  fWC4FiducialCutLow=p.get<double>("WC4LowCut",5.715);
  fWC4FiducialCutHigh=p.get<double>("WC4HighCut",4.285);
  ApplyWC4Fiducialization=p.get<bool>("ApplyWC4Fiducialization",true);


}

std::vector<bool> WCQualityProducer::CheckUSMagApertures(double hit0[3], double hit1[3])
{

bool Mag1USPassX, Mag1USPassY, Mag1DSPassY, Mag1Pass;
//Front face checks
//Get slope-intercept form of track
double XZ_slope=(hit1[0]-hit0[0])/(hit1[2]-hit0[2]);
double XZ_intercept= hit0[0]-XZ_slope*hit0[2];

double z_intersect_us=(XZ_intercept-NDB1_FFace_XZ_intercept)/(NDB1_FFace_XZ_slope-XZ_slope);

//We now have the z coordinate where the track intersects the infinite plane of the US aperture. Get the 3D intersection point.
double USFFHit[3];
projToZ(hit0, hit1, USFFHit, z_intersect_us);

//And check whether this intersection point is in the finite bounds of the magnet

//US X Bound check
if(USFFHit[0]<NDB1_F_top_left[0] && USFFHit[0]>NDB1_F_top_right[0]) {Mag1USPassX=true;}
else{Mag1USPassX=false;}

//US Y Bound Check
if(USFFHit[1]<NDB1_F_top_left[1] && USFFHit[1]>NDB1_F_bottom_left[1]) {Mag1USPassY=true;}
else{Mag1USPassY=false;}  


//Same process for the DS aperture  
double z_intersect_ds=(XZ_intercept-NDB1_BFace_XZ_intercept)/(NDB1_BFace_XZ_slope-XZ_slope);
double USBFHit[3];
projToZ(hit0, hit1, USBFHit, z_intersect_ds); 

//We only check Y here, as there is a bend in X.
if(USBFHit[1]<NDB1_B_top_left[1] && USBFHit[1]>NDB1_B_bottom_left[1]) {Mag1DSPassY=true;}
else{Mag1DSPassY=false;} 

//lastly, get the overall boolean that for this magnet

if(Mag1USPassX && Mag1USPassY && Mag1DSPassY){Mag1Pass=true;}
else{Mag1Pass=false;}

std::vector<bool> result;
result.push_back(Mag1USPassX);
result.push_back(Mag1USPassY);
result.push_back(Mag1DSPassY);
result.push_back(Mag1Pass);
return result;
}

std::vector<bool> WCQualityProducer::CheckDSMagApertures(double hit0[3], double hit1[3])
{
bool Mag2DSPassX, Mag2DSPassY, Mag2USPassY, Mag2Pass;
//Back face checks
//Get slope-intercept form of track
double XZ_slope=(hit1[0]-hit0[0])/(hit1[2]-hit0[2]);
double XZ_intercept= hit0[0]-XZ_slope*hit0[2];

double z_intersect_ds=(XZ_intercept-NDB2_BFace_XZ_intercept)/(NDB2_BFace_XZ_slope-XZ_slope);

//We now have the z coordinate where the track intersects the infinite plane of the DS aperture. Get the 3D intersection point.
double DSBFHit[3];
projToZ(hit0, hit1, DSBFHit, z_intersect_ds);

//And check whether this intersection point is in the finite bounds of the magnet

//DS X Bound check
if(DSBFHit[0]<NDB2_B_top_left[0] && DSBFHit[0]>NDB2_B_top_right[0]) {Mag2DSPassX=true;}
else{Mag2DSPassX=false;}

//US Y Bound Check
if(DSBFHit[1]<NDB2_B_top_left[1] && DSBFHit[1]>NDB2_B_bottom_left[1]) {Mag2DSPassY=true;}
else{Mag2DSPassY=false;}  


//Same process for the US aperture  
double z_intersect_us=(XZ_intercept-NDB2_FFace_XZ_intercept)/(NDB2_FFace_XZ_slope-XZ_slope);
double DSFFHit[3];
projToZ(hit0, hit1, DSFFHit, z_intersect_us); 

//We only check Y here, as there is a bend in X.
if(DSFFHit[1]<NDB2_F_top_left[1] && DSFFHit[1]>NDB2_F_bottom_left[1]) {Mag2USPassY=true;}
else{Mag2USPassY=false;} 

//lastly, get the overall boolean that for this magnet

if(Mag2DSPassX && Mag2DSPassY && Mag2USPassY){Mag2Pass=true;}
else{Mag2Pass=false;}

std::vector<bool> result;
result.push_back(Mag2DSPassX);
result.push_back(Mag2DSPassY);
result.push_back(Mag2USPassY);
result.push_back(Mag2Pass);
return result;
}


std::vector<bool> WCQualityProducer::CheckDSColApertures(double hit0[3], double hit1[3])
{
bool DSColPassX, DSColPassY, DSColPass;
//Front face checks
//Get slope-intercept form of track
double XZ_slope=(hit1[0]-hit0[0])/(hit1[2]-hit0[2]);
double XZ_intercept= hit0[0]-XZ_slope*hit0[2];

double z_intersect_us=(XZ_intercept-DSCol_FFace_XZ_intercept)/(DSCol_FFace_XZ_slope-XZ_slope);

//We now have the z coordinate where the track intersects the infinite plane of the US aperture. Get the 3D intersection point.
double DSColFFHit[3];
projToZ(hit0, hit1, DSColFFHit, z_intersect_us);

//And check whether this intersection point is in the finite bounds of the DSCol

//US X Bound check
if(DSColFFHit[0]<DSCol_F_top_left[0] && DSColFFHit[0]>DSCol_F_top_right[0]) {DSColPassX=true;}
else{DSColPassX=false;}

//US Y Bound Check
if(DSColFFHit[1]<DSCol_F_top_left[1] && DSColFFHit[1]>DSCol_F_bottom_left[1]) {DSColPassY=true;}
else{DSColPassY=false;}  

if(DSColPassX && DSColPassY){DSColPass=true;}
else{DSColPass=false;}
std::vector<bool> result;
result.push_back(DSColPassX);
result.push_back(DSColPassY);
result.push_back(DSColPass);

return result;
}



void  WCQualityProducer::projToZ(double hit0[3], double hit1[3], double (&result)[3] , double zpos){
  //
  // This code is used to find the point at z = zpos along 
  // the line created by hit0 and hit1. Uses the parameterized vector form of a line
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

  result[0] = sx * t + hit0[0]; 
  result[1] = sy * t + hit0[1]; 
  result[2] = zpos;

}


DEFINE_ART_MODULE(WCQualityProducer)


