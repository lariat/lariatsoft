////////////////////////////////////////////////////////////////////////
// Class:       WCTrackAnalyzer
// Module Type: analyzer
// File:        WCTrackAnalyzer_module.cc
//
// Generated at Tue Jun 23 18:58:00 2015 by Matthew Smylie using artmod
// from cetpkgsupport v1_08_06.
////////////////////////////////////////////////////////////////////////

// LArSoft includes
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/AuxDetGeo.h"
//add AuxDetGeoSensitive here when relevant
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "lardataobj/RawData/AuxDetDigit.h"

//lariatsoft includes
#include "LArIATDataProducts/WCTrack.h"
#include "RawDataUtilities/TriggerDigitUtility.h"

// Framework includes
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "fhiclcpp/ParameterSet.h"
#include "cetlib/exception.h"

// ROOT includes. Note: To look up the properties of the ROOT classes,
// use the ROOT web site; e.g.,
// <http://root.cern.ch/root/html532/ClassIndex.html>
#include "TH1.h"
#include "TH2.h"
#include "TTree.h"
#include "TString.h"

// C++ Includes
#include <map>
#include <vector>
#include <algorithm>
#include <iostream>
#include <string>
#include <cmath>
#include <string>
#include <sstream>

namespace wct {

  class WCTrackAnalyzer;

  class WCTrackAnalyzer : public art::EDAnalyzer {
  public:
    explicit WCTrackAnalyzer(fhicl::ParameterSet const & p);
    // The destructor generated by the compiler is fine for classes
    // without bare pointers or other resource use.

    // Plugins should not be copied or assigned.
    WCTrackAnalyzer(WCTrackAnalyzer const &) = delete;
    WCTrackAnalyzer(WCTrackAnalyzer &&) = delete;
    WCTrackAnalyzer & operator = (WCTrackAnalyzer const &) = delete;
    WCTrackAnalyzer & operator = (WCTrackAnalyzer &&) = delete;

    // Required functions.
    void analyze(art::Event const & e) override;

    // Selected optional functions.
    void beginJob() override;
    void reconfigure(fhicl::ParameterSet const & p) ;


  private:

    // Declare member data here.

    std::string fWCTrackLabel; // The name of the producer that made tracks through the MWPCs
    std::string fWCDigLabel; //The name of the producer that created AuxDetDigits with MWPC TDC data
    std::string fTriggerUtility;

    // Pointers to the 1D histograms we'll create. 
    TH1D* fPHist; //reconstructed momentum
    TH1D* fXFaceHist;  //X coord of entry into TPC
    TH1D* fYFaceHist;  //Y coord of entry into TPC
    TH1D* fX_Dist;     //X distance between upstream and downstream tracks
    TH1D* fY_Dist;     //ditto.
    TH1D* fZ_Dist;
    TH1D* fXWireHist;  //coord in terms of wire number
    TH1D* fYWireHist;
    TH1D* fXAxisHist;  //coord in terms of units.
    TH1D* fYAxisHist;
    TH1D* fNTracks;    //number of tracks (per spill)
    TH1D* fWCTDCs;     //histogram of TDC counts from WC AuxDetDigits
    TH1D* fY_Kink;     //angle in Y between upstream and downstream tracks.
    TH1D* fTheta_Dist; //angle of track w.r.t. z axis 
    TH1D* fPhi_Dist;   //angle of track w.r.t. x axis
    TH1D* fNTracksPerTrig; //number of tracks per trigger

    //Pointers to the 2D histograms we'll create.
    TH2D* fXYFaceHist;
    std::vector<TH2D*> fXYWireHist;
    std::vector<TH2D*> fXYAxisHist;

    //histogram limits to avoid recompiling to reconfigure histograms
    size_t fNWC; //number of wire chambers
    double fWCScale; //conversion from WC wires to cm;
    double fPLo;
    double fPHi;
    double fPBins;
    double fXLo;
    double fXHi;
    double fXBins;
    double fYLo;
    double fYHi;
    double fYBins;
    double fNBins;
    double fNLo;  
    double fNHi;
    double fTimeBins;
    double fTimeLo;
    double fTimeHi;
    double fZBins;
    double fZHi;
    double fZLo;
    double fKinkBins;
    double fKinkHi;
    double fKinkLo;
    double fThetaBins;
    double fThetaHi;
    double fThetaLo;
    double fPhiBins;
    double fPhiHi;
    double fPhiLo;
    double fFaceBins;
    double fFaceHi;
    double fFaceLo;

    art::ServiceHandle<geo::Geometry> fGeo;


  };


  WCTrackAnalyzer::WCTrackAnalyzer(fhicl::ParameterSet const & p)
    :
    EDAnalyzer(p)  // ,
    // More initializers here.
  {
    this->reconfigure(p);
  }

  void WCTrackAnalyzer::analyze(art::Event const & event)
  {
    // Implementation of required member function here.

    //Read in the track object made by the producer.
    art::Handle< std::vector<ldp::WCTrack> > trackHandle;
    event.getByLabel(fWCTrackLabel, trackHandle);

    size_t nTracks = trackHandle->size();
    fNTracks->Fill(nTracks); //fill histogram of number of tracks per event
    if(nTracks)
    {
      mf::LogWarning("TrackList") << nTracks << " tracks in event " << event.id().event() << "\n";
    }

    //Get number of tracks associated with each trigger.
    art::Handle< std::vector<raw::Trigger> > trigHandle;
    event.getByLabel(fTriggerUtility, trigHandle);
    art::FindManyP< ldp::WCTrack > fmt(trigHandle,event,fWCTrackLabel);
    for (size_t t=0; t < trigHandle->size(); ++t) {
      fNTracksPerTrig->Fill(fmt.at(t).size());
    }

    for(const auto& track : (*trackHandle)) //trackHandle works somewhat like a pointer to a vector of tracks, so dereference the handle to loop over
    //the vector, then use each "track" as a ldp::WCTrack
    {
      fPHist->Fill(track.Momentum());
      fXFaceHist->Fill(track.XYFace(0)); //indices are 0 for x and 1 for y according to header for WCTrack
      fYFaceHist->Fill(track.XYFace(1));
      fXYFaceHist->Fill(track.XYFace(0), track.XYFace(1));
      fY_Kink->Fill(track.YKink());
      fTheta_Dist->Fill(track.Theta());
      fPhi_Dist->Fill(track.Phi());
      fX_Dist->Fill(track.DeltaDist(0));
      fY_Dist->Fill(track.DeltaDist(1));
      fZ_Dist->Fill(track.DeltaDist(2));
      for(size_t chIt = 0; 2*chIt+1 < track.NHits(); ++chIt)
      {
        fXWireHist->Fill(track.HitWire(2*chIt));
        fYWireHist->Fill(track.HitWire(2*chIt+1));
        fXAxisHist->Fill(track.WC(2*chIt));
        fYAxisHist->Fill(track.WC(2*chIt+1));
        
        //Fill 2D histograms
        fXYWireHist.at(chIt)->Fill(track.HitWire(2*chIt), track.HitWire(2*chIt+1));
        fXYAxisHist.at(chIt)->Fill(track.WC(2*chIt), track.WC(2*chIt+1));
      } //end loop over hit number
    } //end loop over tracks

    //Begin WC AuxDetDigit analysis
    art::Handle<std::vector<raw::AuxDetDigit> > auxDetHandle;
    event.getByLabel(fWCDigLabel, auxDetHandle);

    for(const auto& dig : (*auxDetHandle))
    {
      if(TString(dig.AuxDetName()).Contains("MWPC"))
      {
        for(size_t digIt = 0; digIt < dig.NADC(); ++digIt)
        {
          if(dig.ADC(digIt)) fWCTDCs->Fill(digIt); //actually TDC digits
        }
      }
    }

    return;
  }

  void WCTrackAnalyzer::beginJob()
  {
    // Implementation of optional member function here.
    

    // Access ART's TFileService, which will handle creating and writing
    // histograms and n-tuples for us. 
    art::ServiceHandle<art::TFileService> tfs;

    // Define the histograms. Putting semi-colons around the title
    // causes it to be displayed as the x-axis label if the histogram
    // is drawn.
    fPHist           = tfs->make<TH1D>("momentum", "Momentum of WC Tracks", fPBins, fPLo, fPHi);

    fXFaceHist       = tfs->make<TH1D>("X_Face","X Location of Track's TPC Entry (mm)", fFaceBins, fFaceLo, fFaceHi);

    fYFaceHist       = tfs->make<TH1D>("Y_Face","Y Location of Track's TPC Entry (mm)", fFaceBins, fFaceLo, fFaceHi);    

    fXWireHist       = tfs->make<TH1D>("xwire", "X Wire Values",fXBins, fXLo, fXHi);

    fYWireHist       = tfs->make<TH1D>("ywire", "Y Wire Values",fYBins, fYLo, fYHi);

    fXAxisHist       = tfs->make<TH1D>("xaxis", "X Axis Values",fWCScale*fXBins, fWCScale*fXLo, fWCScale*fXHi);
   
    fYAxisHist       = tfs->make<TH1D>("yaxis", "Y Axis Values", fWCScale*fYBins, fWCScale*fYLo, fWCScale*fYHi);

    fY_Kink          = tfs->make<TH1D>("Y_Kink","Angle between US/DS tracks in Y direction (radians)",fKinkBins,fKinkLo*3.1415926/180,fKinkHi*3.141592654/180);

    fX_Dist          = tfs->make<TH1D>("X_Dist","X distance between US/DS tracks at midplane (mm)",fXBins,fXLo,fXHi);

    fY_Dist          = tfs->make<TH1D>("Y_Dist","Y distance between US/DS tracks at midplane (mm)",fYBins,fYLo,fYHi); //120,-60,60);

    fZ_Dist          = tfs->make<TH1D>("Z_Dist","Z distance between US/DS tracks at midplane (mm)",fZBins,fZLo,fZHi);

    fTheta_Dist      = tfs->make<TH1D>("Theta","Track Theta (w.r.t. TPC Z axis), (radians),",fThetaBins,fThetaLo,fThetaHi);

    fPhi_Dist        = tfs->make<TH1D>("Phi","Track Phi (w.r.t. TPC X axis), (radians)",fPhiBins,fPhiLo,fPhiHi);

    fNTracks         = tfs->make<TH1D>("ntracks", "Number of WC Tracks in Spill", fNBins, fNLo, fNHi);

    fWCTDCs          = tfs->make<TH1D>("wctdcs", "TDC Counts from Wire Chambers", fTimeBins, fTimeLo, fTimeHi);

    fNTracksPerTrig  = tfs->make<TH1D>("nTracksPerTrig","Number of WC Tracks in Trigger",fNBins,fNLo,fNHi);

    //Enter 2D Histograms
    fXYFaceHist      = tfs->make<TH2D>("xyface", "XY Positions on TPC Face", fFaceBins, fFaceLo, fFaceHi, fFaceBins, fFaceLo, fFaceHi);



    for(size_t wcIt = 0; wcIt < fNWC; ++wcIt)
    {
      std::stringstream wcSs;
      wcSs << " for WC " << wcIt; 
      std::stringstream wcName;
      wcName << "wc" << wcIt;

      TH2D* wireHist = tfs->make<TH2D>(("xywire"+wcName.str()).c_str(), ("XY Wire Values" + wcSs.str()).c_str()
                                            ,fXBins, fXLo, fXHi
                                            , fYBins, fYLo, fYHi);

      TH2D* axisHist = tfs->make<TH2D>(("xyaxis"+wcName.str()).c_str(), ("XY Values" + wcSs.str()).c_str()
                                            ,fWCScale*fXBins, fWCScale*fXLo, fWCScale*fXHi, fWCScale*fYBins, fWCScale*fYLo, fWCScale*fYHi);

      wireHist->GetXaxis()->SetTitle("X Wire");
      wireHist->GetYaxis()->SetTitle("Y Wire");

      axisHist->GetXaxis()->SetTitle("X (mm)");
      axisHist->GetYaxis()->SetTitle("Y (mm)");

      fXYWireHist.push_back(wireHist);
      fXYAxisHist.push_back(axisHist);
    }

 
    //Label axes of histograms.


    fPHist->GetXaxis()->SetTitle("Reconstructed momentum (MeV/c)");

    fPHist->GetYaxis()->SetTitle("Tracks per 10 MeV/c");

    fY_Kink->GetXaxis()->SetTitle("Reconstructed y_kink (radians)");

    fY_Kink->GetYaxis()->SetTitle("Tracks per 0.000872 radians");

    fX_Dist->GetXaxis()->SetTitle("X distance between US and DS track ends");

    fX_Dist->GetYaxis()->SetTitle("Tracks per 1 mm");

    fY_Dist->GetXaxis()->SetTitle("Y distance between US and DS track ends");

    fY_Dist->GetYaxis()->SetTitle("Tracks per 1 mm");

    fZ_Dist->GetXaxis()->SetTitle("Z distance between US and DS track ends");

    fZ_Dist->GetYaxis()->SetTitle("Tracks per 1 mm");

    fXFaceHist->GetXaxis()->SetTitle("X (mm)");

    fXFaceHist->GetYaxis()->SetTitle("Tracks per 1 mm");

    fYFaceHist->GetXaxis()->SetTitle("Y (mm)");

    fYFaceHist->GetYaxis()->SetTitle("Tracks per 1 mm");

    fXYFaceHist->GetXaxis()->SetTitle("X (mm)");

    fXYFaceHist->GetYaxis()->SetTitle("Y (mm)");

    fTheta_Dist->GetXaxis()->SetTitle("Theta (radians)");

    fTheta_Dist->GetYaxis()->SetTitle("Tracks per .002 radians");

    fPhi_Dist->GetXaxis()->SetTitle("Phi (radians)");

    fPhi_Dist->GetYaxis()->SetTitle("Tracks per 0.0628 radians");

 

  }

 
  void WCTrackAnalyzer::reconfigure(fhicl::ParameterSet const & p)
  {
    // Implementation of optional member function here.
    
    // Read parameters from the .fcl file. The names in the arguments
    // to p.get<TYPE> must match names in the .fcl file.
    fWCTrackLabel = p.get< std::string >("WCTrackLabel");
    fWCDigLabel   = p.get< std::string >("WCDigLabel");
    fTriggerUtility = "FragmentToDigit";
    fWCScale = p.get< double      >("WCScale");
    fNWC     = p.get< size_t      >("NWC");
    fPBins   = p.get< double      >("MomBins");
    fPLo     = p.get< double      >("MomLo");
    fPHi     = p.get< double      >("MomHi");
    fXBins   = p.get< double      >("XBins");
    fXLo     = p.get< double      >("XLo");
    fXHi     = p.get< double      >("XHi");
    fYBins   = p.get< double      >("YBins"); 
    fYLo     = p.get< double      >("YLo");
    fYHi     = p.get< double      >("YHi");
    fZBins   = p.get< double      >("ZBins"); 
    fZLo     = p.get< double      >("ZLo");
    fZHi     = p.get< double      >("ZHi"); 
    fNBins   = p.get< double      >("NBins");
    fNLo     = p.get< double      >("NLo");
    fNHi     = p.get< double      >("NHi");
    fTimeBins= p.get< double      >("TimeBins");
    fTimeLo  = p.get< double      >("TimeLo");
    fTimeHi  = p.get< double      >("TimeHi");
    fKinkHi  = p.get< double      >("KinkHi");
    fKinkLo  = p.get< double      >("KinkLo");
    fKinkBins  = p.get< double    >("KinkBins");
    fThetaBins = p.get< double    >("ThetaBins");
    fThetaHi  = p.get< double     >("ThetaHi");
    fThetaLo  = p.get< double     >("ThetaLo");
    fPhiLo  = p.get< double       >("PhiLo");
    fPhiHi  = p.get< double       >("PhiHi");
    fPhiBins  = p.get< double     >("PhiBins");
    fFaceBins = p.get< double     >("FaceBins");
    fFaceHi = p.get< double       >("FaceHi");
    fFaceLo = p.get< double       >("FaceLo");
    
    
    return;
  }

  

  DEFINE_ART_MODULE(WCTrackAnalyzer)

} //end namespace
