////////////////////////////////////////////////////////////////////////
/// \file  AuxDetReadoutGeometryT1034.cxx
/// \brief Define the "parallel" geometry that's seen by the LAr Voxels.
/// \author miceli@fnal.gov, talion@fnal.gov
////////////////////////////////////////////////////////////////////////

#include "nutools/G4Base/DetectorConstruction.h"

// Framework includes
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// LArSoft includes
//#include "larsim/LArG4/AuxDetReadoutGeometry.h"
//#include "larsim/LArG4/AuxDetReadout.h"

// LArIATSoft includes
#include "LArIATLArG4/AuxDetReadoutGeometryT1034.h"
#include "LArIATLArG4/AuxDetReadoutT1034.h"

// G4 includes
#include "Geant4/G4PVPlacement.hh"
#include "Geant4/G4PVReplica.hh"
#include "Geant4/G4LogicalVolume.hh"
#include "Geant4/G4VisAttributes.hh"
#include "Geant4/G4VSolid.hh"
#include "Geant4/G4Box.hh"
#include "Geant4/G4Tubs.hh"
#include "Geant4/G4ThreeVector.hh"
#include "Geant4/G4RotationMatrix.hh"
#include "Geant4/G4VSensitiveDetector.hh"
#include "Geant4/G4SDManager.hh"
#include "Geant4/G4Material.hh"
#include "Geant4/G4Point3D.hh"
#include "Geant4/globals.hh"

#include <vector>
#include <cmath>

namespace larg4 {

  // Constructor and destructor.
  AuxDetReadoutGeometryT1034::AuxDetReadoutGeometryT1034(const G4String name)
    : G4VUserParallelWorld(name)
    , fNumSensitiveVol(0)
  {}
  AuxDetReadoutGeometryT1034::~AuxDetReadoutGeometryT1034() 
  {}

  ////////////////////////////////////////////////////////////////////
  void AuxDetReadoutGeometryT1034::Construct()
  {

    // We want to find all of the AuxDets that the Geometry service would
    // find and make each one a G4 sensitive detector.

    //  Call initial case of a function that will rucursively run through
    //  the volume tree down to max depth. Start at 0 depth with World,
    //  where the initial translation and rotation should be 0 as well
    unsigned int MaxDepth = 8; // should be plenty
    std::vector<const G4VPhysicalVolume*> path(MaxDepth);
    path[0] = g4b::DetectorConstruction::GetWorld();
    G4Transform3D InitTransform( path[0]->GetObjectRotationValue(),    
				 path[0]->GetObjectTranslation()    );

    // first try to make sensitive volumes, if those are not in 
    // the gdml file (ie file was made before the introduction of 
    // sensitive volumes) fall back to looking for the aux dets
    this->FindAndMakeAuxDetSensitive(path, 0, InitTransform);
    if(fNumSensitiveVol < 1)
      this->FindAndMakeAuxDet(path, 0, InitTransform);
    

    return;
  }// end Construct
  
  //---------------------------------------------------------------
  void AuxDetReadoutGeometryT1034::FindAndMakeAuxDetSensitive(std::vector<const G4VPhysicalVolume*>& path,
							 unsigned int depth,
							 G4Transform3D DepthToWorld)
  {

    G4LogicalVolume* LogicalVolumeAtDepth = path[depth]->GetLogicalVolume();

    std::string volName(path[depth]->GetName());
    if( volName.find("volAuxDet") != std::string::npos && 
	volName.find("Sensitive") != std::string::npos){

      // find world coordinate of the AuxDet origin in cm
      G4Point3D local(0., 0., 0.);
      G4Point3D world = DepthToWorld * local; // G4 works in mm
      double worldPos[3] = { world.x() / CLHEP::cm, world.y() / CLHEP::cm, world.z() / CLHEP::cm };

      size_t adNum = 0;
      size_t svNum = 0;
      fGeo->FindAuxDetSensitiveAtPosition(worldPos, adNum, svNum);
      //  N.B. This name is expected by code in LArG4:
      std::string SDName = "AuxDetSD_AuxDet" + std::to_string(adNum) + "_" + std::to_string(svNum); 
      AuxDetReadoutT1034* adReadout = new larg4::AuxDetReadoutT1034(SDName, adNum, svNum);

      LOG_DEBUG("AuxDetReadoutGeometryT1034") << "found" << path[depth]->GetName() 
					 << ", number " << adNum << ":" << svNum;

      // Tell Geant4's sensitive-detector manager about the AuxDetReadoutT1034 class
      (G4SDManager::GetSDMpointer())->AddNewDetector(adReadout);
      LogicalVolumeAtDepth->SetSensitiveDetector(adReadout);      
      return;
    }

    // Explore the next layer down -- unless it is a very deep geometry,
    // recursion should end before exception is thrown.
    unsigned int deeper = depth+1;
    if(deeper >= path.size()){
      throw cet::exception("AuxDetReadoutGeometryT1034") << "exceeded maximum TGeoNode depth\n";
    }

    // Note that there will be nd different branches off of path[depth] 
    G4int nd = LogicalVolumeAtDepth->GetNoDaughters();
    for(int d = 0; d < nd; ++d){

      // get the physvol daughter in the logicalvol
      path[deeper] = LogicalVolumeAtDepth->GetDaughter(d);

      // keep track of the transform to world coordinates for PositionToAuxDet
      G4Transform3D DeeperToMother( path[deeper]->GetObjectRotationValue(),    
				    path[deeper]->GetObjectTranslation()    );
      G4Transform3D DeeperToWorld = DepthToWorld * DeeperToMother;
      this->FindAndMakeAuxDetSensitive(path, deeper, DeeperToWorld);
    }

  }

  //---------------------------------------------------------------
  void AuxDetReadoutGeometryT1034::FindAndMakeAuxDet(std::vector<const G4VPhysicalVolume*>& path,
						unsigned int depth,
						G4Transform3D DepthToWorld)
  {

    G4LogicalVolume* LogicalVolumeAtDepth = path[depth]->GetLogicalVolume();

    std::string volName(path[depth]->GetName());
    if( volName.find("volAuxDet") != std::string::npos ){

      // find world coordinate of the AuxDet origin in cm
      G4Point3D local(0., 0., 0.);
      G4Point3D world = DepthToWorld * local; // G4 works in mm
      double worldPos[3] = { world.x() / CLHEP::cm, world.y() / CLHEP::cm, world.z() / CLHEP::cm };

      unsigned int adNum;
      fGeo->PositionToAuxDet(worldPos, adNum);
      //  N.B. This name is expected by code in LArG4:
      std::string SDName = "AuxDetSD_AuxDet" + std::to_string(adNum) + "_0"; 
      AuxDetReadoutT1034* adReadout = new larg4::AuxDetReadoutT1034(SDName, adNum, 0);

      LOG_DEBUG("AuxDetReadoutGeometryT1034") << "found" << path[depth]->GetName() 
					 << ", number " << adNum << ":0";

      // Tell Geant4's sensitive-detector manager about the AuxDetReadoutT1034 class
      (G4SDManager::GetSDMpointer())->AddNewDetector(adReadout);
      LogicalVolumeAtDepth->SetSensitiveDetector(adReadout);      
      return;
    }

    // Explore the next layer down -- unless it is a very deep geometry,
    // recursion should end before exception is thrown.
    unsigned int deeper = depth+1;
    if(deeper >= path.size()){
      throw cet::exception("AuxDetReadoutGeometryT1034") << "exceeded maximum TGeoNode depth\n";
    }

    // Note that there will be nd different branches off of path[depth] 
    G4int nd = LogicalVolumeAtDepth->GetNoDaughters();
    for(int d = 0; d < nd; ++d){

      // get the physvol daughter in the logicalvol
      path[deeper] = LogicalVolumeAtDepth->GetDaughter(d);

      // keep track of the transform to world coordinates for PositionToAuxDet
      G4Transform3D DeeperToMother( path[deeper]->GetObjectRotationValue(),    
				    path[deeper]->GetObjectTranslation()    );
      G4Transform3D DeeperToWorld = DepthToWorld * DeeperToMother;
      this->FindAndMakeAuxDet(path, deeper, DeeperToWorld);
    }

  }


} // namespace larg4
