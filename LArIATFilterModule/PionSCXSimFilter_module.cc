//////////////////////////////////////////////////////////////
// Name:      PionSCXSimFilter
// Date:      24 September 2015
// Author:    Everybody is an author!
//////////////////////////////////////////////////////////////

#ifndef PionSCXSimFilter_Module
#define PionSCXSimFilter_Module

// Framework includes
#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// LArSoft includes
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "nusimdata/SimulationBase/MCParticle.h"

// ROOT includes
#include "TLorentzVector.h"

// C++ includes
#include <cmath>
#include <map>
#include <memory>
#include <vector>

namespace PionSCXSimFilter {

  //-----------------------------------------------------------------------
  //-----------------------------------------------------------------------
  // class definition

  class PionSCXSimFilter : public art::EDFilter {

   public:

    // standard constructor and destructor for an art module
    explicit PionSCXSimFilter(fhicl::ParameterSet const& parameterSet);
    virtual ~PionSCXSimFilter();

    // this method reads in any parameters from the .fcl files
    virtual void reconfigure(fhicl::ParameterSet const& parameterSet) ;

    // the filter routine, called once per event
    bool filter(art::Event & event) override;

   private:

    // returns true if point (x, y, z) is inside TPC
    bool isInsideTPC_(double const& x,
                      double const& y,
                      double const& z);

    // returns distance between two 3D points (x0, y0, z0) and (x1, y1, z1)
    double distance_(double const& x0,
                     double const& y0,
                     double const& z0,
                     double const& x1,
                     double const& y1,
                     double const& z1);

    // returns 3D spatial distance between two four-vectors
    double distance3D_(double const (& aXYZT)[4],
                       double const (& bXYZT)[4]);

    // returns 4D distance between two four-vectors
    double distance4D_(double const (& aXYZT)[4],
                       double const (& bXYZT)[4]);

    // get four-vectors of start/end positions and momenta from particle
    void getFourVectors_(simb::MCParticle const& particle,
                         double               (& startXYZT)[4],
                         double               (& endXYZT)[4],
                         double               (& startPE)[4],
                         double               (& endPE)[4]);

    // get length of a particle's trajectory
    double getTrajectoryLength_(simb::MCParticle const& particle);

    // name of the producer that tracked simulated particles through
    // the detector
    std::string fSimulationProducerLabel;

    // pointer to Geometry provider
    geo::GeometryCore const* fGeometry;

  };


  //-----------------------------------------------------------------------
  // constructor
  PionSCXSimFilter::PionSCXSimFilter(fhicl::ParameterSet const& parameterSet)
  : EDFilter(parameterSet)
  {
    // get a pointer to the geometry service provider
    fGeometry = &*(art::ServiceHandle<geo::Geometry>());

    // reconfigure parameters
    this->reconfigure(parameterSet);
  }


  //-----------------------------------------------------------------------
  // destructor
  PionSCXSimFilter::~PionSCXSimFilter() {}


  //-----------------------------------------------------------------------
  void PionSCXSimFilter::reconfigure(fhicl::ParameterSet const& parameterSet)
  {
    // read parameters from the .fcl file
    fSimulationProducerLabel = parameterSet.get< std::string >("SimulationLabel", "largeant");
  }


  //-----------------------------------------------------------------------
  bool PionSCXSimFilter::filter(art::Event & event)
  {

    // initialize variables for primary particle
    int primaryTrackID = 0;
    int primaryPDGCode = 0;
    bool primaryChargedPion = false;
    bool primaryTrackStopsInTPC = false;

    // declare primary start/end positions and momenta arrays
    double primaryStartXYZT[4];
    double primaryEndXYZT[4];
    double primaryStartPE[4];
    double primaryEndPE[4];

    // initialize variables for neutral pion(s)
    int numberNeutralPionsInEvent = 0;  // counter for neutral pions in event
    int numberSCXNeutralPions = 0;      // counter for neutral pions from SCX interaction point
    int numberSCXNeutrons = 0;          // counter for neutrons from SCX interaction point
    int numberSCXProtons = 0;           // counter for protons from SCX interaction point
    int numberSCXPhotonConversions = 0; // counter for photons from neutral pion decay that convert inside TPC

    // declare vector of neutral pion track IDs
    std::vector<int> neutralPionTrackIDs;

    // get all the simulated particles for the event
    art::ValidHandle< std::vector<simb::MCParticle> >
      particleHandle = event.getValidHandle< std::vector<simb::MCParticle> > (fSimulationProducerLabel);

    // map of pointers to MCParticle objects
    std::map< int, const simb::MCParticle* > particleMap;

    // loop through vector of MCParticle objects
    for (auto const& particle : (*particleHandle))
    {

      // get track ID
      int trackID = particle.TrackId();

      // add the address of the MCParticle to the map, with the track ID
      // as the key
      particleMap[trackID] = &particle;

      // if particle is primary and charged pion
      if (particle.Process() == "primary" &&
          std::abs(particle.PdgCode()) == 211)
      {

        // change flag
        primaryChargedPion = true;

        // get primary track ID
        primaryTrackID = trackID;

        // get primary PDG code
        primaryPDGCode = particle.PdgCode();

        // get primary start and end 4-positions and 4-momenta
        this->getFourVectors_(particle,
                              primaryStartXYZT,
                              primaryEndXYZT,
                              primaryStartPE,
                              primaryEndPE);

        // check to see if primary charged pion track stops inside TPC
        primaryTrackStopsInTPC = this->isInsideTPC_(primaryEndXYZT[0],
                                                    primaryEndXYZT[1],
                                                    primaryEndXYZT[2]);

        // return false if primary track does not stop inside TPC
        if (!primaryTrackStopsInTPC)
          return false;

      } // if primary and charged pion

      // if neutral pion, increment neutral pion counter
      if (particle.PdgCode() == 111)
        ++numberNeutralPionsInEvent;

    } // loop over all particles in the event

    // return false if primary is not a charged pion or if there are no
    // neutral pions in the event
    if (!primaryChargedPion || numberNeutralPionsInEvent < 1)
      return false;

    // loop through vector of MCParticle objects
    for (auto const& particle : (*particleHandle))
    {

      // get track ID
      int trackID = particle.TrackId();

      // get PDG code
      int pdgCode = particle.PdgCode();

      // skip if particle is not a child of the primary
      if (particle.Mother() != primaryTrackID)
        continue;

      // declare start/end positions and momenta arrays
      double startXYZT[4];
      double endXYZT[4];
      double startPE[4];
      double endPE[4];

      // get start end 4-positions and 4-momenta
      this->getFourVectors_(particle, startXYZT, endXYZT, startPE, endPE);

      // get 4D distance between the end position of the primary charged
      // pion track and the start position of the neutral pion
      double distance = this->distance4D_(primaryEndXYZT, startXYZT);

      // if 4D distance is 0
      if (distance == 0)
      {

        // get trajectory length of particle
        double trajectoryLength = this->getTrajectoryLength_(particle);

        mf::LogVerbatim("PionSCXSimFilter")
          << "\n/////////////////////////////////////////////////////"
          << "\nTrack ID:          " << trackID
          << "\nPDG code:          " << pdgCode
          << "\nProcess:           " << particle.Process()
          << "\nTrajectory length: " << trajectoryLength << " cm"
          << "\n/////////////////////////////////////////////////////\n";

        if (pdgCode == 111)
        {
          // keep track of neutral pions that spawn at SCX
          // interaction point
          neutralPionTrackIDs.push_back(trackID);
          // increment neutral pion counter for neutral pions
          // spawning at SCX interaction point
          ++numberSCXNeutralPions;
        } // if neutral pion
        else if (pdgCode == 2112)
        {
          // increment neutron counter for neutrons spawning
          // at SCX interaction point
          ++numberSCXNeutrons;
        } // if neutron
        else if (pdgCode == 2212)
        {
          // increment proton counter for protons spawning
          // at SCX interaction point
          ++numberSCXProtons;
        } // if proton

      } // if 4D distance is 0

    } // loop over all particles in the event

    // return false if there are no neutral pions from SCX
    if (numberSCXNeutralPions < 1)
      return false;

    // loop over neutral pions in the event
    for (auto const& trackID : neutralPionTrackIDs)
    {

      // get particle from map of particles using track ID as the key
      auto const& particle = particleMap[trackID];

      // get number of daughter particles
      int numberDaughters = particle->NumberDaughters();

      for (int i = 0; i < numberDaughters; ++i)
      {

        // get daughter track ID
        int daughterTrackID = particle->Daughter(i);

        // get daughter particle from map of particles using daughter 
        // track ID as the key
        auto const& daughterParticle = particleMap[daughterTrackID];

        if (daughterParticle->PdgCode() == 22 &&
            daughterParticle->Process() == "Decay")
        {

          // declare start/end positions and momenta arrays
          double startXYZT[4];
          double endXYZT[4];
          double startPE[4];
          double endPE[4];

          // get start and end 4-positions and 4-momenta
          this->getFourVectors_(*daughterParticle,
                                startXYZT,
                                endXYZT,
                                startPE,
                                endPE);

          // check to see if photon converts inside TPC
          bool photonConvertsInTPC = this->isInsideTPC_(endXYZT[0],
                                                        endXYZT[1],
                                                        endXYZT[2]);

          // if photon converts inside TPC, increment counter
          // for photons from neutral pion decay that convert
          // inside TPC
          if (photonConvertsInTPC)
            ++numberSCXPhotonConversions;

        } // if particle is a photon from neutral pion decay

      } // loop over daughter particles of neutral pion

    } // loop over neutral pions in the event

    mf::LogVerbatim("PionSCXSimFilter")
      << "\n/////////////////////////////////////////////////////"
      << "\nPrimary PDG code:                        " << primaryPDGCode
      << "\nNumber of neutral pions in event:        " << numberNeutralPionsInEvent
      << "\nNumber of SCX neutral pions:             " << numberSCXNeutralPions
      << "\nNumber of SCX neutrons:                  " << numberSCXNeutrons
      << "\nNumber of SCX protons:                   " << numberSCXProtons
      << "\nNumber of SCX photon conversions in TPC: " << numberSCXPhotonConversions
      << "\n/////////////////////////////////////////////////////\n";

    return true;

  } // PionSCXSimFilter::filter()


  //-----------------------------------------------------------------------
  bool PionSCXSimFilter::isInsideTPC_(double const& x,
                                      double const& y,
                                      double const& z)
  {

    const double length =      fGeometry->DetLength();
    const double width  = 2. * fGeometry->DetHalfWidth();
    const double height = 2. * fGeometry->DetHalfHeight();

    if (x >= 0.         &&
        x <= width      &&
        y >= -height/2. &&
        y <= height/2.  &&
        z >= 0.         &&
        z <= length)
    {
      return true;
    }
    else
    {
      return false;
    }

    return false;

  } // PionSCXSimFilter::isInsideTPC_()


  //-----------------------------------------------------------------------
  double PionSCXSimFilter::distance_(double const& x0,
                                     double const& y0,
                                     double const& z0,
                                     double const& x1,
                                     double const& y1,
                                     double const& z1)
  {

    return std::sqrt((x1 - x0)*(x1 - x0) +
                     (y1 - y0)*(y1 - y0) +
                     (z1 - z0)*(z1 - z0));

  } // PionSCXSimFilter::distance_()


  //-----------------------------------------------------------------------
  double PionSCXSimFilter::distance3D_(double const (& aXYZT)[4],
                                       double const (& bXYZT)[4])
  {

    double const& x0 = aXYZT[0];
    double const& y0 = aXYZT[1];
    double const& z0 = aXYZT[2];

    double const& x1 = bXYZT[0];
    double const& y1 = bXYZT[1];
    double const& z1 = bXYZT[2];

    return std::sqrt((x1 - x0)*(x1 - x0) +
                     (y1 - y0)*(y1 - y0) +
                     (z1 - z0)*(z1 - z0));

  } // PionSCXSimFilter::distance3D_()


  //-----------------------------------------------------------------------
  double PionSCXSimFilter::distance4D_(double const (& aXYZT)[4],
                                       double const (& bXYZT)[4])
  {

    double const& x0 = aXYZT[0];
    double const& y0 = aXYZT[1];
    double const& z0 = aXYZT[2];
    double const& t0 = aXYZT[3];

    double const& x1 = bXYZT[0];
    double const& y1 = bXYZT[1];
    double const& z1 = bXYZT[2];
    double const& t1 = bXYZT[3];

    return std::sqrt((x1 - x0)*(x1 - x0) +
                     (y1 - y0)*(y1 - y0) +
                     (z1 - z0)*(z1 - z0) +
                     (t1 - t0)*(t1 - t0));

  } // PionSCXSimFilter::distance4D_()


  //-----------------------------------------------------------------------
  void PionSCXSimFilter::getFourVectors_(simb::MCParticle const& particle,
                                         double               (& startXYZT)[4],
                                         double               (& endXYZT)[4],
                                         double               (& startPE)[4],
                                         double               (& endPE)[4])
  {

    // get the number of trajectory points of particle
    size_t numberTrajectoryPoints = particle.NumberTrajectoryPoints();

    // get the index of the last trajectory point
    int last = numberTrajectoryPoints - 1;

    // get start/end positions and momenta of particle
    const TLorentzVector & positionStart = particle.Position(0);
    const TLorentzVector & positionEnd   = particle.Position(last);
    const TLorentzVector & momentumStart = particle.Momentum(0);
    const TLorentzVector & momentumEnd   = particle.Momentum(last);

    // fill start/end positions and momenta arrays
    positionStart.GetXYZT(startXYZT);
    positionEnd.GetXYZT(endXYZT);
    momentumStart.GetXYZT(startPE);
    momentumEnd.GetXYZT(endPE);

    return;

  } // PionSCXSimFilter::getFourVectors_()


  //-----------------------------------------------------------------------
  double PionSCXSimFilter::getTrajectoryLength_(simb::MCParticle const& particle)
  {

    // initialize trajectory length
    double length = 0;

    // get the number of trajectory points of particle
    size_t numberTrajectoryPoints = particle.NumberTrajectoryPoints();

    // get the index of the last trajectory point
    size_t last = numberTrajectoryPoints - 1;

    for (size_t i = 0; i < last; ++i)
    {

      // initialize position arrays
      double aXYZT[4];
      double bXYZT[4];

      // get this trajectory point and the next one
      const TLorentzVector & positionA = particle.Position(i);
      const TLorentzVector & positionB = particle.Position(i+1);

      // fill position arrays
      positionA.GetXYZT(aXYZT);
      positionB.GetXYZT(bXYZT);

      // add segment length to total length
      length += this->distance3D_(aXYZT, bXYZT);

    } // loop over trajectory points

    return length;

  } // PionSCXSimFilter::getTrajectoryLength_()


  DEFINE_ART_MODULE(PionSCXSimFilter)

} // namespace PionSCXSimFilter

#endif // PionSCXSimFilter_module
