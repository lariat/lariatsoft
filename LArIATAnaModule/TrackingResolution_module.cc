////////////////////////////////////////////////////////////////////////
// Class:       TrackingResolution
// Module Type: analyzer
// File:        TrackingResolution_module.cc
//
// Generated at Tue June 29 11:21:46 2017 by Elena Gramellini
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
#include "art/Framework/Services/Optional/TFileService.h" 
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "canvas/Persistency/Common/FindOneP.h" 
#include "canvas/Persistency/Common/FindManyP.h"
#include "messagefacility/MessageLogger/MessageLogger.h" 
//#include "cetlib/maybe_ref.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "larcoreobj/SimpleTypesAndConstants/RawTypes.h" // raw::ChannelID_t
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/CryostatGeo.h"
#include "larcorealg/Geometry/TPCGeo.h"
#include "larcorealg/Geometry/PlaneGeo.h"
#include "larcorealg/Geometry/WireGeo.h"
#include "larsim/MCCheater/BackTrackerService.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "larevt/Filters/ChannelFilter.h"
#include "RawDataUtilities/TriggerDigitUtility.h"
// ########################
// ### LArSoft includes ###
// ########################
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "larcoreobj/SimpleTypesAndConstants/RawTypes.h" // raw::ChannelID_t
#include "larcore/Geometry/Geometry.h"
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
#include "lardata/Utilities/AssociationUtil.h"
#include "LArIATDataProducts/WCTrack.h"


#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RawData/raw.h"

// #####################
// ### ROOT includes ###
// #####################
#include "TComplex.h"
#include "TFile.h"
#include "TH2D.h"
#include "TMath.h"
#include "TF1.h"
#include "TGraph2D.h"
#include "TTree.h"
#include "TTimeStamp.h"
#include "TCanvas.h"
#include "TF2.h"
#include "TH1.h"
#include "Math/Functor.h"
#include "TPolyLine3D.h"
#include "Math/Vector3D.h"
#include "Fit/Fitter.h"

#include <vector>
#include <array>
#include "TLinearFitter.h"



#include <map>


struct SumDistance2 {
  // the TGraph is a data member of the object
  TGraph2D *fGraph;
  bool first = true;
  // function Object to be minimized
  SumDistance2(TGraph2D *g) : fGraph(g) {}
  
  // calculate distance line-point
  double distance2(double x,double y,double z, const double *p) {
    // distance line point is D= | (xp-x0) cross  ux |
    // where ux is direction of line and x0 is a point in the line (like t = 0)
    ROOT::Math::XYZVector xp(x,y,z);
    ROOT::Math::XYZVector x0(p[0], p[2], 0. );
    ROOT::Math::XYZVector x1(p[0] + p[1], p[2] + p[3], 1. );
    ROOT::Math::XYZVector u = (x1-x0).Unit();
    double d2 = ((xp-x0).Cross(u)).Mag2();
    return d2;
  }

  // implementation of the function to be minimized
  double operator() (const double *par) {
    assert(fGraph != 0);
    double * x = fGraph->GetX();
    double * y = fGraph->GetY();
    double * z = fGraph->GetZ();
    int npoints = fGraph->GetN();
    double sum = 0;
    /*
    double myMeasure   = 0.;
    double distanceAvg = 0.; 

    for (int i  = 0; i < npoints; ++i) {
      double d = distance2(x[i],y[i],z[i],par);                                                                                                                                                                   
      distanceAvg += TMath::Sqrt(d);
    }
    
    distanceAvg /= npoints ; 
    
    for (int i  = 0; i < npoints; ++i) {
      double d2  = distance2(x[i],y[i],z[i],par);                                                                                                                                        
      double std =  (TMath::Sqrt(d2) - distanceAvg)*(TMath::Sqrt(d2) - distanceAvg);     
      double thisDev = d2/std;
      myMeasure += thisDev;
    }

    myMeasure /=npoints;
    return myMeasure;
    */
    
    // The original idea
    double distanceAvg = 0;
    double distanceStd = 0;
    
    std::cout<<"Hello world"<<std::endl;
    
    for (int i  = 0; i < npoints; ++i) {
      double d = distance2(x[i],y[i],z[i],par);
      distanceAvg += TMath::Sqrt(d);
      sum += d;
    }
    
    distanceAvg /= npoints ;
    
    for (int i  = 0; i < npoints; ++i)      
      {
	double d = TMath::Sqrt(distance2(x[i],y[i],z[i],par));
	distanceStd += ( d - distanceAvg ) * ( d - distanceAvg ) ;
      }
    
    distanceStd /= npoints - 1 ;
    distanceStd = TMath::Sqrt(distanceStd);
    
    double chi2 = sum/((npoints-1)*distanceStd); 
    //if (first)    std::cout << "Total Initial distance square = " << sum  << std::endl;
    //std::cout << "Chi Square -----------------------------------> " << chi2 << std::endl;
    first = false;   
    return chi2;
    
  }
};



namespace lariat 
{
  class TrackingResolution;
}

class lariat::TrackingResolution : public art::EDAnalyzer 
{
public:
  explicit TrackingResolution(fhicl::ParameterSet const & p);
  virtual ~TrackingResolution();

  // Required functions.
  void analyze(art::Event const & e) override;

  // Selected optional functions.
  void   beginJob();
  void   reconfigure(fhicl::ParameterSet const & p);
  
  //void   line(double t, const double *p, double &x, double &y, double &z);
  double calculateChi2(std::vector<const double*> points, std::vector<double>& parFit );
  double Alpha(std::vector<double> parLine1, std::vector<double> parLine2 );
  //double fitMe(std::vector<double*> points, std::vector<double>& parFit );
  
private:
  std::string fTrackModuleLabel; 
  std::string fWCTrackLabel; 
  std::string fWC2TPCModuleLabel;
 
  int    f_nSpt_Gap;
  int    f_nSpt_Min_Half;
  double f_chi2_Cut;  

  TH1D*  h_Chi2_Full;
  TH1D*  h_Chi2_1Half;
  TH1D*  h_Chi2_2Half;

  TH1D*  h_NSpt_Full;
  TH1D*  h_NSpt_1Half;
  TH1D*  h_NSpt_2Half;

  TH1D*  h_DeltaAlpha;
  TH1D*  h_Alpha_1Half;
  TH1D*  h_Alpha_2Half;

  TTree* fTree; 
  int    nSpt_Full;
  int    nSpt_1Half;
  int    nSpt_2Half;
  float  chi2_Full;
  float  chi2_1Half;
  float  chi2_2Half;
  float  delta_Alpha;
  float  alpha_1Half;
  float  alpha_2Half;


};


lariat::TrackingResolution::TrackingResolution(fhicl::ParameterSet const & pset) : EDAnalyzer(pset)
{
  this->reconfigure(pset);
}

lariat::TrackingResolution::~TrackingResolution()
{
}

/*
// define the parametric line equation
void lariat::TrackingResolution::line(double t, const double *p, double &x, double &y, double &z) {
   // a parametric line is define from 6 parameters but 4 are independent
   // x0,y0,z0,z1,y1,z1 which are the coordinates of two points on the line
   // can choose z0 = 0 if line not parallel to x-y plane and z1 = 1;
   x = p[0] + p[1]*t;
   y = p[2] + p[3]*t;
   z = t;
}
*/

double lariat::TrackingResolution::calculateChi2(std::vector<const double*> points, std::vector<double>& parFit )
{
  TGraph2D * gr = new TGraph2D();
  // Fill the 2D graph
  size_t sillyIt = 0;
  for (auto s : points)  {gr->SetPoint(sillyIt,s[0],s[1],s[2]); sillyIt++;}
  ROOT::Fit::Fitter  fitter; 
  // make the functor objet 
  SumDistance2 sdist(gr);
  ROOT::Math::Functor fcn(sdist,4); 
  // set the function and the initial parameter values  
  double pStart[4] = {1,1,1,1};  
  fitter.SetFCN(fcn,pStart);
  // set step sizes different than default ones (0.3 times parameter values) 
  for (int i = 0; i < 4; ++i) fitter.Config().ParSettings(i).SetStepSize(0.01);  
  bool ok = fitter.FitFCN();

  if (!ok) { 
    //Error("line3Dfit","Line3D Fit failed");  
    return -100000; 
  }

  const ROOT::Fit::FitResult & result = fitter.Result();
  std::cout << "Total final chi2  " << result.MinFcnValue() << std::endl;  

  result.Print(std::cout);  
  // get fit parameters
  const double * parFit_bogus = result.GetParams();
  parFit.clear();
  for (size_t i = 0; i < 4; i++) parFit.push_back( parFit_bogus[i] );
  
  return result.MinFcnValue();

}

double lariat::TrackingResolution::Alpha(std::vector<double> parLine1, std::vector<double> parLine2 )
{
  
  if (parLine1.size() != 4 ) throw cet::exception("Geometry")    /// <== ... throw exception
			       << "cannot find right number of parameters line1"<<std::endl;
  if (parLine2.size() != 4 ) throw cet::exception("Geometry")    /// <== ... throw exception
			       << "cannot find right number of parameters line2"<<std::endl;

  double numerator   = parLine1[1]*parLine2[1] + parLine1[3]*parLine2[3] + 1;
  double denominator = (parLine1[1]*parLine1[1] + parLine1[3]*parLine1[3] + 1) * (parLine2[1]*parLine2[1] + parLine2[3]*parLine2[3] + 1);
  denominator = TMath::Sqrt(denominator);
  double positiveCosAlpha = numerator/denominator;
  double alpha = TMath::ACos(positiveCosAlpha);

  if (alpha > 0 )
    if (alpha <    TMath::Pi()/2. ) return alpha;
    else  return TMath::Pi() - alpha;
  else
    if (alpha > -1*TMath::Pi()/2. && alpha < 0) return -1*alpha;
    else return TMath::Pi() + alpha;
}

void lariat::TrackingResolution::reconfigure(fhicl::ParameterSet const & pset)
{
  
  f_nSpt_Gap        = pset.get<   int  >("nSpt_Gap");
  f_nSpt_Min_Half   = pset.get<   int  >("nSpt_Min_Half");
  f_chi2_Cut        = pset.get< double >("chi2_Cut");

  fTrackModuleLabel = pset.get< std::string >("TrackModuleLabel"      , "pmtrack");
  fWCTrackLabel     = pset.get< std::string >("WCTrackLabel"          , "wctrack"  );
  fWC2TPCModuleLabel= pset.get< std::string >("WC2TPCModuleLabel"     , "wctracktpctrackmatch");

  return;
  
}

void lariat::TrackingResolution::analyze(art::Event const & evt)
{
  art::Handle< std::vector<ldp::WCTrack> > wctrackHandle;   // Container of wc tracks  
  std::vector<art::Ptr<ldp::WCTrack> >     wctrack;         // Vector of wc tracks  
  art::Handle< std::vector<recob::Track> > trackListHandle; // Container of reco tracks
  std::vector<art::Ptr<recob::Track> >     tracklist;       // Vector of wc tracks  

  // Find which recontructed track is associated with the WC
  int matchedRecoTrkKey = -99999;
  if(!evt.getByLabel(fWCTrackLabel, wctrackHandle)) return;
  if(!evt.getByLabel(fTrackModuleLabel,trackListHandle)) return;

  // Fill the container of reco tracks
  art::fill_ptr_vector(wctrack, wctrackHandle);
  art::fill_ptr_vector(tracklist, trackListHandle);
  
  art::FindOneP<recob::Track> fWC2TPC(wctrackHandle, evt, fWC2TPCModuleLabel); 
  if (fWC2TPC.isValid())   
    {
      for (unsigned int indexAssn = 0; indexAssn < fWC2TPC.size(); ++indexAssn ) 
	{
	  cet::maybe_ref<recob::Track const> trackWC2TPC(*fWC2TPC.at(indexAssn)); 
	  if (!trackWC2TPC) continue;   
	  recob::Track const& aTrack(trackWC2TPC.ref()); 
	  matchedRecoTrkKey = aTrack.ID(); // This checks out OK   
	}
    }// if there's an associated track
  // Now we know what ID identifies the reco track we're interested in.
  if (matchedRecoTrkKey == -99999) return;
  

  //Find space points
  art::FindManyP<recob::SpacePoint> fmsp(trackListHandle, evt, fTrackModuleLabel); 

  //Loop on the reco tracks
  for (auto recoTrk : tracklist ) 
    {
      //we keep only the matched one
      if (recoTrk->ID() != matchedRecoTrkKey) continue;
      //This is the reco track I want!!!!
      
      //Work flow
      // 1) Get SpacePoints
      // 2) Calculate Full Chi2 and nSpt_Full
      // 2.5) Skip if nSpt < 2*f_nSpt_Min_Half+f_nSpt_Gap
      // 3) Order by Z position
      // 4) Split in 2
      // 4.1) Remove last  f_nSpt_Gap/2 from end       of first half
      // 4.2) Remove first f_nSpt_Gap/2 from beginning of second half
      // 5) Fit halves --> chi2_halves
      // 6) Skip if chi2_1half or chi2_2half is shit 
      // 7) Calculate alpha and delta alpha

      // Let's start
      // 1) Get SpacePoints
      std::vector<art::Ptr<recob::SpacePoint> > const& spts = fmsp.at(recoTrk.key()); 
      // 3) Order by Z position
      std::map<double, const double*> orderedSpts;
      std::vector<const double*>     spts_Full;
      std::vector<const double*>     spts_1Half;
      std::vector<const double*>     spts_2Half;
      
      std::vector<double>     par_Full;
      std::vector<double>     par_1Half;
      std::vector<double>     par_2Half;

      
      for(size_t iSpts=0; iSpts<spts.size() ;++iSpts)
        {
	  auto currentSpacePoint =  spts.at(iSpts)->XYZ();
	  orderedSpts[currentSpacePoint[2]]= currentSpacePoint ;
	  spts_Full.push_back(currentSpacePoint);
	}
      
      // 2) Calculate Full Chi2 and nSpt_Full
      // std::cout<<"Calculating the Chi2 for full track -------------> "<< evt.run()<<" "<<evt.subRun()<<" "<<evt.event()<<"\n";
      nSpt_Full = spts.size();
      if (nSpt_Full < 3 ) continue;
      chi2_Full = calculateChi2(spts_Full, par_Full);
      //std::cout<<"-------------------------------------------------> Full "<< chi2_Full <<"\n";     
      if (chi2_Full < -9999. ) continue;
      // 2.5) Skip if nSpt < 2*f_nSpt_Min_Half+f_nSpt_Gap
      if ( nSpt_Full < 2*f_nSpt_Min_Half+f_nSpt_Gap ) 
	{
	  orderedSpts.clear();   spts_Full.clear();
	  spts_1Half.clear();    spts_2Half.clear();
	  continue;
	}

      // 4) Split in 2
      // 4.1) Remove last  f_nSpt_Gap/2 from end       of first half
      // 4.2) Remove first f_nSpt_Gap/2 from beginning of second half
      std::map<double,const double*>::iterator it_fistHalf   = orderedSpts.begin();
      std::map<double,const double*>::iterator it_secondHalf = orderedSpts.begin();
      size_t stop_halves = (size_t) ( (float)(spts.size() - f_nSpt_Gap)/2.);
      for (size_t i = 0; i <  stop_halves              ; i++) it_fistHalf++;
      for (size_t i = 0; i <  stop_halves+f_nSpt_Gap   ; i++) it_secondHalf++;
      
      for (std::map<double,const double*>::iterator it=orderedSpts.begin() ; it!=it_fistHalf      ; ++it) spts_1Half.push_back(it->second);
      for (std::map<double,const double*>::iterator it=it_secondHalf       ; it!=orderedSpts.end(); ++it) spts_2Half.push_back(it->second);
      
      // 5) Fit halves --> chi2_halves
      chi2_1Half = calculateChi2(spts_1Half,par_1Half); 
      std::cout<<"-------------------------------------------------> 1Half "<< chi2_1Half <<"\n";     
      if (chi2_1Half < -9999. ) continue;
      chi2_2Half = calculateChi2(spts_2Half,par_2Half); 
      std::cout<<"-------------------------------------------------> 2Half "<< chi2_2Half <<"\n";     
      if (chi2_2Half < -9999. ) continue;
      nSpt_1Half = spts_1Half.size();
      nSpt_2Half = spts_2Half.size();
      // 6) Skip if chi2_1half or chi2_2half is shit 
      //if ( chi2_1Half > 4 || chi2_1Half < 0.001 ) continue;
      //if ( chi2_2Half > 4 || chi2_2Half < 0.001 ) continue;
      // 7) Calculate alpha and delta alpha
      alpha_1Half = Alpha(par_Full, par_1Half);
      alpha_2Half = Alpha(par_Full, par_2Half);

      delta_Alpha = alpha_1Half - alpha_2Half;



      fTree->Fill();
      h_Chi2_Full ->Fill(chi2_Full);
      h_Chi2_1Half->Fill(chi2_1Half);
      h_Chi2_2Half->Fill(chi2_2Half);
      
      h_NSpt_Full ->Fill(nSpt_Full);
      h_NSpt_1Half->Fill(nSpt_1Half);
      h_NSpt_2Half->Fill(nSpt_2Half);
      
      h_DeltaAlpha ->Fill(delta_Alpha*180./TMath::Pi());
      h_Alpha_1Half->Fill(alpha_1Half*180./TMath::Pi());
      h_Alpha_2Half->Fill(alpha_2Half*180./TMath::Pi());
       
    }// Loop on tracks
  

}


void lariat::TrackingResolution::beginJob()
{
  art::ServiceHandle<art::TFileService> tfs;
  h_Chi2_Full    = tfs->make<TH1D>("h_Chi2_Full"   ,"h_Chi2_Full"    ,400,-200,200);
  h_Chi2_1Half   = tfs->make<TH1D>("h_Chi2_1Half"  ,"h_Chi2_1Half"   ,400,-200,200);
  h_Chi2_2Half   = tfs->make<TH1D>("h_Chi2_2Half"  ,"h_Chi2_2Half"   ,400,-200,200);
  h_NSpt_Full    = tfs->make<TH1D>("h_NSpt_Full"   ,"h_NSpt_Full"    ,500,   0,500);
  h_NSpt_1Half   = tfs->make<TH1D>("h_NSpt_1Half"  ,"h_NSpt_1Half"   ,500,   0,500);
  h_NSpt_2Half   = tfs->make<TH1D>("h_NSpt_2Half"  ,"h_NSpt_2Half"   ,500,   0,500);               			                                   
  h_DeltaAlpha   = tfs->make<TH1D>("h_DeltaAlpha"  ,"h_DeltaAlpha"   ,720, -90,90);
  h_Alpha_1Half  = tfs->make<TH1D>("h_Alpha_1Half" ,"h_Alpha_1Half"  ,720,   0,180);
  h_Alpha_2Half  = tfs->make<TH1D>("h_Alpha_2Half" ,"h_Alpha_2Half"  ,720,   0,180);

  fTree = tfs->make<TTree>("trackResTree","analysis tree");
  fTree->Branch("nSpt_Full"   ,&nSpt_Full   ,"nSpt_Full/I"  );
  fTree->Branch("nSpt_1Half"  ,&nSpt_1Half  ,"nSpt_1Half/I"  );
  fTree->Branch("nSpt_2Half"  ,&nSpt_2Half  ,"nSpt_2Half/I"  );

  fTree->Branch("chi2_Full"    ,&chi2_Full    ,"chi2_Full/F"    );
  fTree->Branch("chi2_1Half"   ,&chi2_1Half   ,"chi2_1Half/F"   );
  fTree->Branch("chi2_2Half"   ,&chi2_2Half   ,"chi2_2Half/F"   );
  fTree->Branch("delta_Alpha"  ,&delta_Alpha  ,"delta_Alpha/F"  );
  fTree->Branch("alpha_1Half"  ,&alpha_1Half  ,"alpha_1Half/F"  );
  fTree->Branch("alpha_2Half"  ,&alpha_2Half  ,"alpha_2Half/F"  );


}



DEFINE_ART_MODULE(lariat::TrackingResolution)
