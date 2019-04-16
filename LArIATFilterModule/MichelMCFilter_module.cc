////////////////////////////////////////////////////////////////////////
// Class:       MichelMCFilter
// Module Type: filter
// File:        MichelMCFilter_module.cc
//
// This module is intended to be used in series during the cosmic muon
// Monte Carlo generation stage.  It requires that the generated muon
// enters or stops inside the TPC.
//
// Generated at Tue Sep  6 15:40:14 2016 by William Foreman using artmod
// from cetpkgsupport v1_10_02.
////////////////////////////////////////////////////////////////////////

// Framework includes
#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

//ROOT Includes
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TH1F.h"
#include "TH2F.h"

// LArSoft Includes
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "larsim/Simulation/LArG4Parameters.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "lardataobj/RawData/TriggerData.h"

#include <memory>

class MichelMCFilter;

class MichelMCFilter : public art::EDFilter {
public:
  explicit MichelMCFilter(fhicl::ParameterSet const & p);

  // Plugins should not be copied or assigned.
  MichelMCFilter(MichelMCFilter const &) = delete;
  MichelMCFilter(MichelMCFilter &&) = delete;
  MichelMCFilter & operator = (MichelMCFilter const &) = delete;
  MichelMCFilter & operator = (MichelMCFilter &&) = delete;

  // Required functions.
  bool filter(art::Event & e) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;
  void reconfigure(fhicl::ParameterSet const & p) override;
  bool IsPointInFiducialVolume(TVector3,double,double,double);
  bool IsPointInFiducialVolume(TLorentzVector,double,double,double);

private:
   
  geo::GeometryCore     *fGeometry;

  // Filter requirements
  bool		      fRequireStoppingMuon;
  bool		      fRequireEnteringMuon;
  bool		      fRequireDecayElectron;
  float               fMinDecayTime;
  float               fMaxDecayTime;

  // Producer label from gen stage
  std::string         fSimProducerLabel;
  
  // TPC dimensions
  double	      TPC_Range_X[2];
  double	      TPC_Range_Y[2];
  double              TPC_Range_Z[2];

  // Diagnostic histograms
  TH1F* hEventCount;
  TH1F* hMuMomentum;
  TH1F* hMuMomentum_pass;
  TH1F* hMuDecayTime;
  TH1F* hMuDecayTime_pass;
  TH1F* hStoppingMuEndptX;
  TH1F* hStoppingMuEndptY;
  TH1F* hStoppingMuEndptZ;
  TH2F* hStoppingMuEndptZX;
  TH2F* hStoppingMuEndptZY;

  // Counters
  int numMuPlus;
  int numMuMinus;
  int numMuPlusDecay;
  int numMuMinusDecay;
  int numMuPlusStop;
  int numMuMinusStop;

};


MichelMCFilter::MichelMCFilter(fhicl::ParameterSet const & p)
{
  // Read in fhicl parameters
  this->reconfigure(p);
  
  // Get a pointer to the geometry service provider
  fGeometry = &*(art::ServiceHandle<geo::Geometry>());
  TPC_Range_X[0]  =  0.00;
  TPC_Range_X[1]  =  2.*fGeometry->DetHalfWidth();
  TPC_Range_Y[0]  =  -1.*fGeometry->DetHalfHeight();
  TPC_Range_Y[1]  =  fGeometry->DetHalfHeight();
  TPC_Range_Z[0]  =  0.00;
  TPC_Range_Z[1]  =  fGeometry->DetLength();
  
  std::cout
  <<"MichelMCFilter: using TPC dimensions \n"
  <<"  x = ("<<TPC_Range_X[0]<<", "<<TPC_Range_X[1]<<")\n"
  <<"  y = ("<<TPC_Range_Y[0]<<", "<<TPC_Range_Y[1]<<")\n"
  <<"  x = ("<<TPC_Range_Z[0]<<", "<<TPC_Range_Z[1]<<")\n";

  numMuPlus=0;
  numMuMinus=0;
  numMuPlusDecay=0;
  numMuMinusDecay=0;
  numMuPlusStop=0;
  numMuMinusStop=0;

}

bool MichelMCFilter::filter(art::Event & e)
{

  //std::cout<<"Beginning event filter.\n";
  hEventCount->Fill(0);

  bool muonStopsInTPC = false;
  bool muonEntersTPC = false;
  bool isDecayElectron = false;
  bool isPrimaryElectron = false;
  bool outFlag = false;
  int muID = -99;
  int muCharge = 0;
  float muStopTime = -999.;
  float muMomentum = -99.;
  float decayTime = -999.;

  // Get MCParticles from event
  art::Handle< std::vector<simb::MCParticle> > particleHandle;
  e.getByLabel(fSimProducerLabel, particleHandle);
 
  // Loop through particles and look for primary muon
  for( auto const& particle : (*particleHandle) )
  {
    if( particle.Process() == "primary" && abs(particle.PdgCode()) == 13 ) {
      int Npts = particle.NumberTrajectoryPoints();
      int last = Npts - 1;
      muMomentum  = particle.P(0)*1000.;
      muID        = particle.TrackId();
      muCharge    = -1.*(particle.PdgCode()/abs(particle.PdgCode()));
      muStopTime  = particle.T(last);
      //std::cout<<"Found primary muon, charge = "<<muCharge<<"\n";
      
      hMuMomentum->Fill( muMomentum );
      
      // Add to counter
      if( muCharge == -1 ) numMuMinus++;
      if( muCharge == 1 ) numMuPlus++;
 
    
      // Does it stop in TPC?
      if( IsPointInFiducialVolume( particle.Position(last), 0., 0., 0.) ) {
        
        // Add to counter
        if( muCharge == -1 )  numMuMinusStop++;
        if( muCharge == 1 )   numMuPlusStop++;
        
        muonStopsInTPC = true;

        hStoppingMuEndptX->Fill( particle.Position(last).X() );
        hStoppingMuEndptY->Fill( particle.Position(last).Y() );
        hStoppingMuEndptZ->Fill( particle.Position(last).Z() );
        hStoppingMuEndptZX->Fill( particle.Position(last).Z(), particle.Position(last).X() );
        hStoppingMuEndptZY->Fill( particle.Position(last).Z(), particle.Position(last).Y() );

          
      }
      
      
      for(int ii=0; ii<Npts; ii++)
	if( IsPointInFiducialVolume( particle.Position(ii), 0., 0., 0.) ) muonEntersTPC = true;
      
    }
   
    // -------------------------------------------------------------------------- 
    // Look for decay electron
    if( abs(particle.PdgCode()) == 11 && particle.Mother() == muID && muonStopsInTPC ) {
      
      // since the Geant4 "Process" param is useless in case of mu-, 
      // make some clever cuts to ensure we really have a Michel:
      float dT = particle.T(0) - muStopTime; 
      if( -1.*particle.PdgCode()/abs(particle.PdgCode()) == muCharge && dT > 1.0 ) {
        decayTime  = dT;
        if( muCharge == 1 ) numMuPlusDecay++;
        if( muCharge == -1 ) numMuMinusDecay++;
        isDecayElectron = true;
        std::cout<<"--> Found decay electron\n";
        hMuDecayTime->Fill( decayTime );
      }
    }
    
    // -------------------------------------------------------------------------- 
    // Look for a primary electron shower
    if( particle.Process() == "primary" && abs(particle.PdgCode()) == 11 ) {
      
        isPrimaryElectron = true;
        std::cout<<"--> Found primary electron shower\n";
    }

  }

  //std::cout<<"Does it enter TPC?   "<<muonEntersTPC<<"\n";
  //std::cout<<"Does it stop in TPC? "<<muonStopsInTPC<<"\n";

  if(  ( (!fRequireStoppingMuon) || (fRequireStoppingMuon && muonStopsInTPC) )
    && ( (!fRequireEnteringMuon) || (fRequireEnteringMuon && muonEntersTPC)  ) 
    && ( (!fRequireDecayElectron) || (fRequireDecayElectron && isDecayElectron)  ) 
    && ( (fMinDecayTime <  0. || decayTime > fMinDecayTime ) )
    && ( (fMaxDecayTime <  0. || decayTime < fMaxDecayTime ) )
    )
    outFlag = true;

  // If this is just an isolated electron that was simulated, 
  // pass it no matter what!
  if( isPrimaryElectron ) outFlag = true;

  if(outFlag) {
    hEventCount->Fill(1);
    hMuMomentum_pass->Fill( muMomentum );
    hMuDecayTime_pass->Fill( decayTime );
  }
  std::cout<<"Passing event? "<<outFlag<<"\n";
  return outFlag;

}

void MichelMCFilter::beginJob()
{
  // Open the TFileService
  art::ServiceHandle<art::TFileService> tfs;
  hEventCount = tfs->make<TH1F>("EventCount","MichelMCFilter event count",2,0,2);
  hEventCount ->GetXaxis()->SetBinLabel(1,"Generated");  
  hEventCount ->GetXaxis()->SetBinLabel(2,"Passed");  
  hMuMomentum         = tfs->make<TH1F>("MuMomentum",     "Initial muon momentum;Momentum [MeV/c]",100.,0.,1000.);
  hMuMomentum_pass    = tfs->make<TH1F>("MuMomentum_pass","Initial muon momentum (passing filter);Momentum [MeV/c]",100.,0.,1000.);
  hMuDecayTime        = tfs->make<TH1F>("MuDecayTime",      "Muon decay time [ns]",375,0,7500);
  hMuDecayTime_pass   = tfs->make<TH1F>("MuDecayTime_pass", "Muon decay time (passing filter) [ns]",375,0,7500);
  hStoppingMuEndptX  = tfs->make<TH1F>("StoppingMuEndptX","Stopping muon endpoint;X [cm]",55,-5.,50.);
  hStoppingMuEndptY  = tfs->make<TH1F>("StoppingMuEndptY","Stopping muon endpoint;Y [cm]",50,-25.,25.);
  hStoppingMuEndptZ  = tfs->make<TH1F>("StoppingMuEndptZ","Stopping muon endpoint;Z [cm]",100,-5.,95.);
  hStoppingMuEndptZX = tfs->make<TH2F>("StoppingMuEndptZX","Stopping muon endpoint;Z [cm];X [cm]",100,-5.,95.,55,-5.,50.);
    hStoppingMuEndptZX->SetOption("colz");
  hStoppingMuEndptZY = tfs->make<TH2F>("StoppingMuEndptZY","Stopping muon endpoint;Z [cm];Y [cm]",100,-5.,95.,50,-25.,25.);
    hStoppingMuEndptZY->SetOption("colz");
}

void MichelMCFilter::endJob()
{
  std::cout
  <<"============================================\n"
  <<"MichelMCFilter\n"
  <<"  total mu+          = "<<numMuPlus<<"\n" 
  <<"  -> stops in TPC    = "<<numMuPlusStop<<"\n"
  <<"  -> michel decays   = "<<numMuPlusDecay<<"\n"
  <<"  total mu-          = "<<numMuMinus<<"\n"
  <<"  -> stops in TPC    = "<<numMuMinusStop<<"\n"
  <<"  -> michel decays   = "<<numMuMinusDecay<<"\n"
  <<"============================================\n";
}


void MichelMCFilter::reconfigure(fhicl::ParameterSet const & p)
{
  fSimProducerLabel       = p.get< std::string >  ("SimProducer","largeant");
  fRequireStoppingMuon	  = p.get< bool >	  ("RequireStoppingMuon",true);
  fRequireEnteringMuon	  = p.get< bool >	  ("RequireEnteringMuon",true);
  fRequireDecayElectron	  = p.get< bool >	  ("RequireDecayElectron",true);
  fMinDecayTime           = p.get< float >        ("MinDecayTime",-9.);
  fMaxDecayTime           = p.get< float >        ("MaxDecayTime",-9.);
}

// Function for determining if a point is inside or outside
// predefined fiducial volume (originally from MichelMCAna_module.cc)
bool MichelMCFilter::IsPointInFiducialVolume(TVector3 p, double fX, double fY, double fZ)
{
  if(    (p.X() >= TPC_Range_X[0] + fX) && (p.X() <= TPC_Range_X[1] - fX)
      && (p.Y() >= TPC_Range_Y[0] + fY) && (p.Y() <= TPC_Range_Y[1] - fY)
      && (p.Z() >= TPC_Range_Z[0] + fZ) && (p.Z() <= TPC_Range_Z[1] - fZ) )
  {
    return 1;
  } else {
    return 0;
  }
}

bool MichelMCFilter::IsPointInFiducialVolume(TLorentzVector pl, double fX, double fY, double fZ)
{
  TVector3 p(pl.X(),pl.Y(),pl.Z());
  return IsPointInFiducialVolume(p,fX,fY,fZ);
}



DEFINE_ART_MODULE(MichelMCFilter)
