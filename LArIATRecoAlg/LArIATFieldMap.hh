#ifndef LARIATRECOALG_SMITHFELD_HH
#define LARIATRECOALG_SMITHFELD_HH

#include "Geant4/globals.hh"
#include <vector>
#include "Geant4/G4Types.hh"
#include "Geant4/G4MagneticField.hh"

class LArIATFieldMap : public G4MagneticField 
{

public:
  explicit LArIATFieldMap(std::string const & fieldMap);

  virtual ~LArIATFieldMap();

  virtual void GetFieldValue(const G4double [4], G4double *B) const;

  void SetBField(double cur_B) { fBField = abs(cur_B); }

  bool bVerbose = true;

private:

  void TransformToMagCords(double og[3], double trans_mag1[3], double trans_mag2[3]) const;
  void SpaceToArrayCords(double og[3], size_t transformed[3]) const;

  //double NDBN_vol[3] = {317.5, 142.2, 591.0}; // unused

  double NDB1_center[3] = {553.06, 21.6165, -4757.22}; //553.31, 22.86, -4745.33};
  double NDB1_rot[3] = {0.0, 10.5, 0.0};

  double NDB2_center[3] = {465.03, 25.17, -4066.5}; // 465.05, 23.34, -4054.22};
  double NDB2_rot[3] = {0.0, 5.5, 0.0};

  std::vector< std::vector<std::vector<std::vector<double> > > >  fieldMapper;
  const size_t numX = 196;
  const size_t numY = 34;
  const size_t numZ = 801;
  const size_t numB = 3;

  double fBField = 0.0;

};

#endif 
