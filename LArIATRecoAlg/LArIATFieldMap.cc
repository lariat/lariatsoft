#include "LArIATFieldMap.hh"
#include <vector>
#include "TMath.h"
#include "cetlib/exception.h"
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

#include "Geant4/G4SystemOfUnits.hh"

// Defined coordinates of the field map
#define HEIGHT 196
#define WIDTH 3
#define DEPTH 7

LArIATFieldMap::LArIATFieldMap(std::string const & fieldMap) {

  // Takes the field map file, and loads it into an array
  // Array is going to have indices that is easy to go from index to xyz cm space. 

  std::ifstream infile;

  long i = 0;
  long j = 0;

  fieldMapper.resize(numX);
  for(size_t i = 0; i < numX; i++) {
    fieldMapper[i].resize(numY);
    for(size_t j = 0; j < numY; j++) {
      fieldMapper[i][j].resize(numZ);
      for(size_t k = 0; k < numZ; k++) {
	fieldMapper[i][j][k].resize(numB);
      }
    }
  }
    

  double temp_values[7];

  char cNum[10];
  infile.open (fieldMap, std::ifstream::in);
  if (infile.is_open()) {
    if(bVerbose) { std::cout << "Loading the B-field map. " << std::endl; }
    while (infile.good()) {
      if(j % 1000000 == 0) { if(bVerbose) { std::cout << j << " lines scanned " << std::endl; } }
      j += 1;

      infile.getline(cNum, 256, ' ');

      temp_values[i] = std::atof(cNum);

      i += 1;

      if(i == 7) { 
	i = 0; 

	//std::cout << temp_values[0] << " " << temp_values[1] << " " << temp_values[2] << " " 
	//<< temp_values[3] << " " << temp_values[4] << " " << temp_values[5] << " " 
	//<< temp_values[6] << " " << std::endl; 

	size_t arrayValue[3];
	double posMag1[3] = {temp_values[1], temp_values[2], temp_values[3]};
	SpaceToArrayCords(posMag1, arrayValue);

	//std::cout << " array indices " << arrayValue[0] << " " << arrayValue[1] << " " << arrayValue[2] << std::endl;

	if(arrayValue[0] < numX and arrayValue[1] < numY and arrayValue[2] < numZ) {
	  fieldMapper[arrayValue[0]][arrayValue[1]][arrayValue[2]][0] = temp_values[4];
	  fieldMapper[arrayValue[0]][arrayValue[1]][arrayValue[2]][1] = temp_values[5];
	  fieldMapper[arrayValue[0]][arrayValue[1]][arrayValue[2]][2] = temp_values[6];
	} 

      } 
         

    }
    infile.close();
  } else { 
    throw cet::exception("LArIATFieldMap")
      << "Could not load the geometry file.\n";
  }

}

LArIATFieldMap::~LArIATFieldMap() {
  
}

void LArIATFieldMap::TransformToMagCords(double og[3], double trans_mag1[3], double trans_mag2[3]) const {

  // Going to create a translation matrix
  // and a rotation matrix  
  double og_temp[3] = {og[0], og[1], og[2]};

  // Inverse translation
  og_temp[0] -= NDB1_center[0];
  og_temp[1] -= NDB1_center[1];
  og_temp[2] -= NDB1_center[2];

  // an inverse rotation matrix 
  trans_mag1[0] =  og_temp[0] * cos(NDB1_rot[1] * TMath::Pi() / 180.0) + og_temp[2] * sin(NDB1_rot[1] * TMath::Pi() / 180.0);
  trans_mag1[1] =  og_temp[1];
  trans_mag1[2] = -og_temp[0] * sin(NDB1_rot[1] * TMath::Pi() / 180.0) + og_temp[2] * cos(NDB1_rot[1] * TMath::Pi() / 180.0);

  // Move into correct units
  trans_mag1[0] /= 1000.0;
  trans_mag1[1] /= 1000.0;
  trans_mag1[2] /= 1000.0;

  double og_temp2[3] = {og[0], og[1], og[2]};

  // Translate into the other magnet
  og_temp2[0] -= NDB2_center[0];
  og_temp2[1] -= NDB2_center[1];
  og_temp2[2] -= NDB2_center[2];

  // Rotation
  trans_mag2[0] =  og_temp2[0] * cos(NDB2_rot[1] * TMath::Pi() / 180.0) + og_temp2[2] * sin(NDB2_rot[1] * TMath::Pi() / 180.0);
  trans_mag2[1] =  og_temp2[1];
  trans_mag2[2] = -og_temp2[0] * sin(NDB2_rot[1] * TMath::Pi() / 180.0) + og_temp2[2] * cos(NDB2_rot[1] * TMath::Pi() / 180.0);

  // Units
  trans_mag2[0] /= 1000.0;
  trans_mag2[1] /= 1000.0;
  trans_mag2[2] /= 1000.0;

}

void LArIATFieldMap::GetFieldValue (const G4double local[4], G4double *B)	const
{

  // Return B field from a given location

  // Code first transforms the local variable into magnet coordinates
  //  then checks what the field map says is the magnetic field there
 
  double cords_to_pass[3] = {local[0], local[1], local[2]};

  double posMag1[3];
  double posMag2[3];
  TransformToMagCords(cords_to_pass, posMag1, posMag2);

  size_t arrayValue[3];

  //std::cout << local[0] << " " << local[1] << " " << local[2] << " " << std::endl;
  //std::cout << posMag1[0] << " " << posMag1[1] << " " << posMag1[2] << " " << std::endl;
  //std::cout << posMag2[0] << " " << posMag2[1] << " " << posMag2[2] << " " << std::endl;

  // check which magnet the particle is closer to
  if(TMath::Sqrt(pow(posMag1[2], 2)) < TMath::Sqrt(pow(posMag2[2], 2))) {     
    double passer1[3] = {posMag1[0], abs(posMag1[1]), posMag1[2]};
    SpaceToArrayCords(passer1, arrayValue);
  } else {
    double passer2[3] = {posMag2[0], abs(posMag2[1]), posMag2[2]};
    SpaceToArrayCords(passer2, arrayValue);
    posMag1[0] = posMag2[0]; // cheaty so I don't need if statements later
    posMag1[1] = posMag2[1]; // cheaty so I don't need if statements later
    posMag1[2] = posMag2[2]; // cheaty so I don't need if statements later
  }  	 
  
  
  if(arrayValue[0]+1 < numX and arrayValue[1]+1 < numY and arrayValue[2]+1 < numZ) {

    B[0] = -(fBField / 0.31056) * fieldMapper[arrayValue[0]][arrayValue[1]][arrayValue[2]][0] * tesla;
    B[1] = -(fBField / 0.31056) * fieldMapper[arrayValue[0]][arrayValue[1]][arrayValue[2]][1] * tesla;
    B[2] = -(fBField / 0.31056) * fieldMapper[arrayValue[0]][arrayValue[1]][arrayValue[2]][2] * tesla;

  } else {
    B[0] = 0.0;
    B[1] = 0.0;
    B[2] = 0.0;
  }


  //std::cout << local[0] << " " << local[1] << " "<< local[2] << " " << B[0] << " " << B[1] << " " << B[2] << std::endl;

}

void LArIATFieldMap::SpaceToArrayCords(double og[3], size_t transformed[3]) const
{

  // The Field map is an array with coordinates:

  // [-0.195, 0.195] at 0.002 increments, or 196 entries
  // [ 0.0, 0.06] at 0.002 increments, or 34 entries
  // [-0.8, 0.8] at 0.002 increments, or 800 entries

  transformed[0] = int((og[0] + 0.195) / 0.002);
  transformed[1] = int((og[1] + 0.0) / 0.002);
  transformed[2] = int((og[2] + 0.8) / 0.002);

}
