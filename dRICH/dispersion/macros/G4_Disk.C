#ifndef MACRO_G4_Disk_C
#define MACRO_G4_Disk_C

#include <GlobalVariables.C>

#include <g4detectors/PHG4CylinderSubsystem.h>
#include <g4main/PHG4Reco.h>

#include <fun4all/Fun4AllServer.h>

#include <Geant4/G4SystemOfUnits.hh>
#include <Geant4/G4Material.hh>

#include <algorithm>
#include <cmath>
#include <vector>

using namespace std;

namespace mySim
{
  double dRICH_width = 80;
  pair<double, double> dRICH_eta_coverage = {1, 4};
  double dRICH_z = 240;
  int nLayers = 100;
}

using namespace mySim;

double etaToR(double eta, double z)
{
  double radius = 2*z*atan(exp(-1*eta));
  return radius;
}

void DiskInit()
{
  BlackHoleGeometry::max_radius = max(BlackHoleGeometry::max_radius, etaToR(dRICH_eta_coverage.first, dRICH_z));
  BlackHoleGeometry::min_z = min(BlackHoleGeometry::min_z, -1.);
  BlackHoleGeometry::max_z = max(BlackHoleGeometry::max_z, dRICH_z + 0.5*dRICH_width); 
}

double Disk(PHG4Reco *g4Reco, double max_radius)
{
  PHG4CylinderSubsystem *cyl = nullptr;

  double diskWidth = dRICH_width/nLayers;
  pair<double, double> dRICH_r_coverage(etaToR(dRICH_eta_coverage.second, dRICH_z),
                                        etaToR(dRICH_eta_coverage.first, dRICH_z));

  for (int i = 0; i < nLayers; ++i)
  {
    double z_position = dRICH_z - 0.5*dRICH_width + (i+0.5)*diskWidth;

    if (Enable::VERBOSITY > 2)
    {
      std::cout << "Layer " << i << std::endl;
      std::cout << "Disk width: = " << diskWidth << std::endl;
      std::cout << "IR = " << dRICH_r_coverage.first << ", OR = " << dRICH_r_coverage.second << std::endl;
      std::cout << "Z position = " << z_position << std::endl;
    }

    cyl = new PHG4CylinderSubsystem("disk", i);
    cyl->Verbosity(0);
    cyl->set_double_param("place_z", z_position);
    cyl->set_double_param("radius", dRICH_r_coverage.first);
    cyl->set_double_param("length", 0.5*diskWidth);
    cyl->set_string_param("material", "G4_AIR");
    cyl->set_color(0.4, 0.4, 0.4, 0.4);
    cyl->set_double_param("thickness", (dRICH_r_coverage.second - dRICH_r_coverage.first));
    cyl->SuperDetector("MySuperDetector");
    cyl->OverlapCheck(false);
    cyl->SetActive();
    g4Reco->registerSubsystem(cyl);
  }

  return dRICH_r_coverage.first;
}

#endif
