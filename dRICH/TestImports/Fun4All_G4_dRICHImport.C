#ifndef FUN4ALL_G4_MOMENTUM_C
#define FUN4ALL_G4_MOMENTUM_C

#include "DisplayOn.C"
#include <G4_Input.C>

#include <g4detectors/PHG4CylinderSubsystem.h>

#include <g4trackfastsim/PHG4TrackFastSim.h>
#include <g4trackfastsim/PHG4TrackFastSimEval.h>

#include <g4main/PHG4ParticleGenerator.h>
#include <g4main/PHG4Reco.h>
#include <g4main/PHG4TruthSubsystem.h>

#include <fun4all/Fun4AllDstOutputManager.h>
#include <fun4all/Fun4AllDummyInputManager.h>
#include <fun4all/Fun4AllInputManager.h>
#include <fun4all/Fun4AllOutputManager.h>
#include <fun4all/Fun4AllServer.h>
#include <fun4all/SubsysReco.h>
#include <g4detectors/PHG4GDMLSubsystem.h>

#include <phool/recoConsts.h>

#include <cmath>

R__LOAD_LIBRARY(libfun4all.so)
R__LOAD_LIBRARY(libg4testbench.so)
R__LOAD_LIBRARY(libg4detectors.so)
R__LOAD_LIBRARY(libg4trackfastsim.so)

int Fun4All_G4_dRICHImport(const int nEvents = -1 // negative value run a Geant4 GUI for event display
    )
{
  ///////////////////////////////////////////
  // Make the Server
  //////////////////////////////////////////
  Fun4AllServer *se = Fun4AllServer::instance();
  se->Verbosity(0);

  recoConsts *rc = recoConsts::instance();
  // if you want to use a fixed seed for reproducible results
  //  rc->set_IntFlag("RANDOMSEED", 12345);

//  // PHG4ParticleGenerator generates particle
//  // distributions in eta/phi/mom range
//  PHG4ParticleGenerator *gen = new PHG4ParticleGenerator("PGENERATOR");
//  gen->set_name("pi-");
//  gen->set_vtx(0, 0, 0);
//  gen->set_eta_range(2, 2.5);            // around midrapidity
//  gen->set_mom_range(10, 10);                  // fixed 4 GeV/c
//  gen->set_phi_range(0., 90. / 180. * M_PI);  // 0-90 deg
//  se->registerSubsystem(gen);


  Input::PYTHIA6 = true;
  //-----------------
  // Initialize the selected Input/Event generation
  //-----------------
  InputInit();

  INPUTGENERATOR::Pythia6->set_config_file(string(getenv("CALIBRATIONROOT")) + "/Generators/phpythia6_ep.cfg");

  // register all input generators with Fun4All
  InputRegister();


  PHG4Reco *g4Reco = new PHG4Reco();
  g4Reco->set_field(1.5);  // 1.5 T constant solenoidal field

  PHG4GDMLSubsystem* gdml = new PHG4GDMLSubsystem("dRICH");
//  gdml->set_string_param("GDMPath", "drich.gdml");
  gdml->set_string_param("GDMPath", "drich_only_check.gdml");
  gdml->set_string_param("TopVolName", "DRICH");
  gdml->set_double_param("place_z", 262.);
  gdml->OverlapCheck(true);
  g4Reco->registerSubsystem(gdml);

  PHG4TruthSubsystem *truth = new PHG4TruthSubsystem();
  g4Reco->registerSubsystem(truth);

  se->registerSubsystem(g4Reco);

  if (nEvents >= 0)
  {
    se->run(nEvents);
    // finish job - close and save output files
    se->End();
    std::cout << "All done" << std::endl;

    // cleanup - delete the server and exit
    delete se;
    gSystem->Exit(0);
  }
  else
  {
    // negative value run a Geant4 GUI for event display
    QTGui();

  }
  return 0;
}

PHG4ParticleGenerator *get_gen(const char *name = "PGENERATOR")
{
  Fun4AllServer *se = Fun4AllServer::instance();
  PHG4ParticleGenerator *pgun = (PHG4ParticleGenerator *) se->getSubsysReco(name);
  return pgun;
}

#endif
