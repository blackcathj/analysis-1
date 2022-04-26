#ifndef MACRO_FUN4ALLPHOTONDISPERSION_C
#define MACRO_FUN4ALLPHOTONDISPERSION_C

#include <GlobalVariables.C>

#include <DisplayOn.C>
#include <G4Setup.C>
#include <G4_Global.C>
#include <G4_Input.C>

#include <fun4all/Fun4AllServer.h>

#include <phool/PHRandomSeed.h>
#include <phool/recoConsts.h>

#include <getvectors/getVectors.h>

R__LOAD_LIBRARY(libfun4all.so)
R__LOAD_LIBRARY(libgetvectors.so)

using namespace mySim;

int Fun4All_photonDispersion(
    const int nEvents = 1,
    const string particle = "pi+",
    const string productionNumber = "00000",
    const string revisionNumber = "000")
{
  //The next set of lines figures out folder revisions, file numbers etc
  string outDir = "./";
  if (outDir.substr(outDir.size() - 1, 1) != "/") outDir += "/";
  outDir += particle + "/" + revisionNumber + "/";
  string outputFileName = "outputData_" + particle + "_" + revisionNumber + "_" + productionNumber + ".root";

  string outputRecoDir = outDir + "/inReconstruction/";
  string makeDirectory = "mkdir -p " + outputRecoDir;
  system(makeDirectory.c_str());
  string outputRecoFile = outputRecoDir + outputFileName;

  Fun4AllServer *se = Fun4AllServer::instance();
  se->Verbosity(1);

  //Opt to print all random seed used for debugging reproducibility. Comment out to reduce stdout prints.
  PHRandomSeed::Verbosity(1);

  //===============
  // Input options
  //===============
  // verbosity setting (applies to all input managers)
  Enable::VERBOSITY = 2;
  // First enable the input generators
  Input::SIMPLE = true;
  // Input::SIMPLE_NUMBER = 2; // if you need 2 of them
  // Input::SIMPLE_VERBOSITY = 1;

  InputInit();

  //--------------
  // Set generator specific options
  //--------------
  // can only be set after InputInit() is called

  // Simple Input generator:
  // if you run more than one of these Input::SIMPLE_NUMBER > 1
  // add the settings for other with [1], next with [2]...
  INPUTGENERATOR::SimpleEventGenerator[0]->add_particles(particle, 100);
  INPUTGENERATOR::SimpleEventGenerator[0]->set_vertex_distribution_function(PHG4SimpleEventGenerator::Gaus,
                                                                            PHG4SimpleEventGenerator::Gaus,
                                                                            PHG4SimpleEventGenerator::Gaus);
  INPUTGENERATOR::SimpleEventGenerator[0]->set_vertex_distribution_mean(0., 0., 0.);
  INPUTGENERATOR::SimpleEventGenerator[0]->set_vertex_distribution_width(0., 0., 0.);
  INPUTGENERATOR::SimpleEventGenerator[0]->set_eta_range(dRICH_eta_coverage.first, dRICH_eta_coverage.second);
  INPUTGENERATOR::SimpleEventGenerator[0]->set_phi_range(-M_PI, M_PI);
  INPUTGENERATOR::SimpleEventGenerator[0]->set_p_range(30., 30.);
  // register all input generators with Fun4All
  InputRegister();

  // turn the display on (default off)
  Enable::DISPLAY = false;

  G4MAGNET::magfield_rescale = 1.;  // make consistent with expected Babar field strength of 1.4T

  // Initialize the selected subsystems
  Enable::PIPE = true;

  G4Init();
  G4Setup();
  Global_Reco();

  //-----------------
  // Global Vertexing
  //-----------------

  if (Enable::GLOBAL_RECO && Enable::GLOBAL_FASTSIM)
  {
    cout << "You can only enable Enable::GLOBAL_RECO or Enable::GLOBAL_FASTSIM, not both" << endl;
    gSystem->Exit(1);
  }
  if (Enable::GLOBAL_RECO)
  {
    Global_Reco();
  }
  else if (Enable::GLOBAL_FASTSIM)
  {
    Global_FastSim();
  }

  InputManagers();

  getVectors *dRICHCalc = new getVectors();
  dRICHCalc->setOutputName(outputRecoFile);
  se->registerSubsystem(dRICHCalc);
  //-----------------
  // Event processing
  //-----------------
  if (Enable::DISPLAY)
  {
    DisplayOn();

    gROOT->ProcessLine("Fun4AllServer *se = Fun4AllServer::instance();");
    gROOT->ProcessLine("PHG4Reco *g4 = (PHG4Reco *) se->getSubsysReco(\"PHG4RECO\");");

    cout << "-------------------------------------------------" << endl;
    cout << "You are in event display mode. Run one event with" << endl;
    cout << "se->run(1)" << endl;
    cout << "Run Geant4 command with following examples" << endl;
    gROOT->ProcessLine("displaycmd()");

    return 0;
  }

  // if we use a negative number of events we go back to the command line here
  if (nEvents < 0)
  {
    return 0;
  }
  // if we run the particle generator and use 0 it'll run forever
  if (nEvents == 0)
  {
    cout << "using 0 for number of events is a bad idea when using particle generators" << endl;
    cout << "it will run forever, so I just return without running anything" << endl;
    return 0;
  }

  se->run(nEvents);

  //-----
  // Exit
  //-----

  se->End();

  ifstream file(outputRecoFile.c_str());
  if (file.good())
  {
    string moveOutput = "mv " + outputRecoFile + " " + outDir;
    system(moveOutput.c_str());
  }

  std::cout << "All done" << std::endl;
  delete se;

  gSystem->Exit(0);
  return 0;
}
#endif
