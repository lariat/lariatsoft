#ifndef BEAMLINEMASSALG_H
#define BEAMLINEMASSALG_H

// C++ includes
#include <iostream>
#include <memory>
#include <math.h>

// Art includes
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// LArSoft includes
#include "larcore/Geometry/Geometry.h"
#include "lardataobj/RecoBase/Track.h"
#include "nusimdata/SimulationBase/MCParticle.h"

// LArIATSoft includes
#include "LArIATDataProducts/WCTrack.h"
#include "LArIATDataProducts/TOF.h"

//--------------------------------------------
class BeamlineMassAlg{
 public:
  
  //Constructor/destructor
  BeamlineMassAlg( fhicl::ParameterSet const& pset );
  BeamlineMassAlg();
  ~BeamlineMassAlg();

  void    reconfigure( fhicl::ParameterSet const& pset );

  // setters
  void    SetLength( float L ){ fLength = L; }
  void    SetLightTravelTime( float Tc){ fTc = Tc; }
  void    SetTOFModuleLabel(std::string label){ fTOFModuleLabel = label;}
  void    SetWCTrackModuleLabel(std::string label){ fWCTrackModuleLabel = label;}

  // getters
  float   GetMass( float tof, float p);
  float   GetMass( art::Event const &evt);
  float   GetTOF( art::Event const &evt);
  float   GetWCTrackMomentum( art::Event const &evt);
  
 private:
   
  std::string fWCTrackModuleLabel;
  std::string fTOFModuleLabel;
  float fLength;
  float fTc;
  
};

#endif
