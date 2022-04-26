#ifndef MACRO_G4SETUP_C
#define MACRO_G4SETUP_C

#include <GlobalVariables.C>

#include <G4_BlackHole.C>
#include <G4_Magnet.C>
#include "G4_Disk.C"
#include <G4_World.C>
#include <G4_Pipe.C>

#include <g4main/PHG4Reco.h>
#include <g4main/PHG4TruthSubsystem.h>

#include <phfield/PHFieldConfig.h>

#include <g4decayer/EDecayType.hh>

#include <fun4all/Fun4AllServer.h>

R__LOAD_LIBRARY(libg4detectors.so)

using namespace mySim;

void G4Init()
{
  DiskInit();
  MagnetFieldInit();
  BlackHoleInit();
}

int G4Setup()
{
  //---------------
  // Fun4All server
  //---------------
  Fun4AllServer *se = Fun4AllServer::instance();

  PHG4Reco *g4Reco = new PHG4Reco();
  g4Reco->set_rapidity_coverage(dRICH_eta_coverage.first);

  WorldInit(g4Reco);

  double fieldstrength;
  istringstream stringline(G4MAGNET::magfield);
  stringline >> fieldstrength;
  if (stringline.fail())
  {  // conversion to double fails -> we have a string
    if (G4MAGNET::magfield.find("sphenix3dbigmapxyz") != string::npos ||
        G4MAGNET::magfield == "CDB")
    {
      g4Reco->set_field_map(G4MAGNET::magfield, PHFieldConfig::Field3DCartesian);
    }
    else
    {
      g4Reco->set_field_map(G4MAGNET::magfield, PHFieldConfig::kField2D);
    }
  }
  else
  {
    g4Reco->set_field(fieldstrength);  // use const soleniodal field
  }
  g4Reco->set_field_rescale(G4MAGNET::magfield_rescale);

  double radius = 0.;

  Disk(g4Reco, radius);
  BlackHole(g4Reco, radius);

  PHG4TruthSubsystem *truth = new PHG4TruthSubsystem();
  g4Reco->registerSubsystem(truth);
  // finally adjust the world size in case the default is too small
  WorldSize(g4Reco, radius);
  
  se->registerSubsystem(g4Reco);
  return 0;
}

#endif
