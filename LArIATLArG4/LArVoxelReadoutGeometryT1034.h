////////////////////////////////////////////////////////////////////////
/// \file  LArVoxelReadoutGeometryT1034.h
/// \brief Define the "parallel" geometry that's seen by the LAr Voxels.
///
/// \author  seligman@nevis.columbia.edu
////////////////////////////////////////////////////////////////////////
///
/// This class defines the parallel geometry that will be divided into
/// the three-dimensional voxels for the detector read-out.
/// 
/// Why define a parallel geometry?  Here are some reasons:
///
/// - The regular LAr TPC is one large volume of liquid argon.  When
///   Geant4 does its physics modeling, it can be unconstrained in
///   step size by the voxels.  Only for readout would the steps be
///   sub-divided.
///
/// - There may be more than one kind of readout, depending on a
///   detector's instrumentation (e.g., OpDets in addition to the wire
///   planes).  It's possible that the voxelization appropriate for
///   the wire planes may not be an appropriate readout for the other
///   readouts.  Geant4 allows the construction of multiple parallel
///   readouts, so this mechanism is relatively easy to extend for
///   each type of readout.

#ifndef LArIATLArG4_LArVoxelReadoutGeometryT1034_h
#define LArIATLArG4_LArVoxelReadoutGeometryT1034_h

//#include "larsim/LArG4/LArVoxelReadout.h"
#include "larcore/Geometry/Geometry.h"
#include "Geant4/G4VUserParallelWorld.hh"
#include "Geant4/G4String.hh"
#include "Geant4/G4Transform3D.hh"
#include "Geant4/G4VPhysicalVolume.hh"
#include "Geant4/G4UserLimits.hh"

#include "LArIATLArG4/LArVoxelReadoutT1034.h"

// Forward declarations
namespace CLHEP { class HepRandomEngine; }

namespace larg4 {

  class LArVoxelReadoutGeometryT1034 : public G4VUserParallelWorld
  {
  public:
    
    /// Collection of all it takes to set up this object.
    struct Setup_t {
      
      /// Set up data for `LArVoxelReadoutT1034`.
      larg4::LArVoxelReadoutT1034::Setup_t readoutSetup;
      
    }; // struct Setup_t
    
    
    /// Constructor: sets up all its LArVoxelReadoutT1034 instances.
    LArVoxelReadoutGeometryT1034(const G4String name, Setup_t const& setupData);

    /// The key method in this class; creates a parallel world view of
    /// those volumes relevant to the LAr voxel readout.  Required of
    /// any class that inherits from G4VUserParallelWorld
    virtual void Construct();

  private:

    G4VPhysicalVolume* FindNestedVolume(G4VPhysicalVolume* mother,
					G4Transform3D&     motherTransform,
					G4Transform3D&     daughterTransform,
					std::string&       daughterName,
					unsigned int       expectedNum);

    art::ServiceHandle<geo::Geometry> fGeo;       ///< Handle to the geometry    
    std::unique_ptr<G4UserLimits>     fStepLimit; ///< G4 doesn't handle memory management, 
                                                  ///< so we have to
    
    /// Data for `LArVoxelReadoutT1034` setup.
    larg4::LArVoxelReadoutT1034::Setup_t fReadoutSetupData;
    
  };

} // namespace larg4

#endif // LArIATLArG4_LArVoxelReadoutGeometryT1034_h
