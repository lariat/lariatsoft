/////////////////////////////////////////////////////////////////////
// BeamlineMassAlg
// 
// Contains the function which calculates particle mass using beamline 
// information like TOF and momentum. The event must have already
// undergone TOF and WC track reconstruction and therefore must have:
//  - ldp::TOF
//  - ldp::WCTrack
//
//  Note that there are two different (but ultimately equivalent) 
//  methods of doing the calculation: 
//    (a) using some assumed total path length of particles as they 
//        travel from the upstream TOF to the downstream TOF ("L");
//    (b) using the time a particle traveling at the speed of light 
//        would take to travel between the TOF paddles ("Tc").
//
//  Historically, we in LArIAT have used the "Length" as the reference
//  input here, with L~6.6m. More recent calibration work on the other
//  hand has made use of method (b), where Tc is used as input instead. 
//  Since Tc = L/c, the two methods produce the same results as long as 
//  Tc and L are set to be compatible.
//
//  ----------------------------------------------------------------
//
//  To include in your module, add the header file at the top:
//    #include "LArIATRecoAlg/BeamlineMassAlg.h"
//
//  then declare a private mass calculation alg in your constructor:
//    BeamlineMassAlg   fMassAlg;
//
//  then in your code, to get mass, simply do:
//    fMassAlg->GetMass( evt );
//    or
//    fMassAlg->GetMass(TOF, momentum); 
//
//  where "evt" is the art::Event. You can set L or Tc by doing this:
//    fMassAlg->SetLength( <new L value>);
//    fMassAlg->SetLightTravelTime( <new Tc value> );
//
//  By default, L is used in the calculation, *unless* the user sets
//  the light travel time as shown above, in which case the time is
//  used.
//
//  W. Foreman, May 2020
//  wforeman @ iit.edu
//
/////////////////////////////////////////////////////////////////////
#include "LArIATRecoAlg/BeamlineMassAlg.h"

//--------------------------------------------------------------
//Constructor
BeamlineMassAlg::BeamlineMassAlg( fhicl::ParameterSet const& pset )
{
  this->reconfigure(pset);
}

BeamlineMassAlg::BeamlineMassAlg( )
{
  fTOFModuleLabel       = "tof";
  fWCTrackModuleLabel   = "wctrack";
  fLength               = 6.685;
  fTc                   = -1;
}


//--------------------------------------------------------------  
//Destructor
BeamlineMassAlg::~BeamlineMassAlg()
{
}

//--------------------------------------------------------------
void BeamlineMassAlg::reconfigure( fhicl::ParameterSet const& pset ){
  fTOFModuleLabel       = pset.get< std::string > ("TOFModuleLabel","tof");
  fWCTrackModuleLabel   = pset.get< std::string > ("WCTrackModuleLabel","wctrack");
  fLength               = pset.get< float >       ("Length", 6.685 );      // m
  fTc                   = pset.get< float >       ("LightTravelTime", -1); // ns (22.3ns)
}

//---------------------------------------------------------------
// Calculate particle mass [MeV/c^2] from beamline momentum and TOF
float BeamlineMassAlg::GetMass( float TOF, float P ) {
  float m = -999.;
  float c = 0.299792458;    // m/ns
  float tc = fLength / c;
  if( fTc > 0 ) tc = fTc;
  float radical = (TOF*TOF)/(tc*tc) - 1.;
  if( radical >= 0. ) m = P * sqrt( radical );
  else                m = -P * sqrt( -radical ); 
  return m;
}

float BeamlineMassAlg::GetMass( art::Event const & evt ){
  return GetMass( GetTOF(evt), GetWCTrackMomentum(evt) );
}

//---------------------------------------------------------------
// grabs the first saved TOF (assumes prior TOF filtering was done)
float BeamlineMassAlg::GetTOF(art::Event const & e){
  float t = -999.;
  art::Handle< std::vector<ldp::TOF> > tofHandle;
  std::vector<art::Ptr<ldp::TOF> > tof;
  if(e.getByLabel(fTOFModuleLabel,tofHandle)) {
    art::fill_ptr_vector(tof, tofHandle);
    if( tof.size() > 0 && tof[0]->NTOF() > 0 ) t = tof[0]->SingleTOF(0);
  }
  return t; 
}

//---------------------------------------------------------------
// grabs the first saved wctrack (assumes prior filtering was done)
float BeamlineMassAlg::GetWCTrackMomentum(art::Event const & e){
  float p = -999.;
  art::Handle< std::vector<ldp::WCTrack> > wctrackHandle; 
  std::vector<art::Ptr<ldp::WCTrack> > wctrack; 
  if(e.getByLabel(fWCTrackModuleLabel, wctrackHandle)) {
    art::fill_ptr_vector(wctrack, wctrackHandle);
    if( wctrack.size() > 0 ) p = wctrack[0]->Momentum();
  } 
  return p;
}
