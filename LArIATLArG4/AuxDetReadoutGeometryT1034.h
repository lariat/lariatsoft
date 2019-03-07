////////////////////////////////////////////////////////////////////////
/// \file AuxDetReadoutGeometryT1034.h
/// \brief Define the "parallel" geometry that's seen by the AuxDet.
/// \author miceli@fnal.gov, talion@fnal.gov
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

#ifndef LArIATLArG4_AuxDetReadoutGeometryT1034_h
#define LArIATLArG4_AuxDetReadoutGeometryT1034_h

#include "larcore/Geometry/Geometry.h"
#include "Geant4/G4VUserParallelWorld.hh"
#include "Geant4/G4String.hh"
#include "Geant4/G4Transform3D.hh"
#include "Geant4/G4VPhysicalVolume.hh"

namespace larg4 {

  class AuxDetReadoutGeometryT1034 : public G4VUserParallelWorld
  {
  public:
    /// Constructor and destructor.
    AuxDetReadoutGeometryT1034( const G4String name = "AuxDetReadoutGeometry" );
    virtual ~AuxDetReadoutGeometryT1034();

    /// The key method in this class; creates a parallel world view of
    /// those volumes relevant to the auxiliary detector readout.
    /// Required of  any class that inherits from G4VUserParallelWorld
    virtual void Construct();

  private:

    void FindAndMakeAuxDet(std::vector<const G4VPhysicalVolume*>& path,
			   unsigned int depth,
			   G4Transform3D DepthToWorld);    

    void FindAndMakeAuxDetSensitive(std::vector<const G4VPhysicalVolume*>& path,
				    unsigned int depth,
				    G4Transform3D DepthToWorld);    

    art::ServiceHandle<geo::Geometry> fGeo;             ///< Handle to the geometry    
    uint32_t                          fNumSensitiveVol; ///< number of sensitive volumes

  };

} // namespace larg4

#endif // LArIATLArG4_AuxDetReadoutGeometryT1034_h
