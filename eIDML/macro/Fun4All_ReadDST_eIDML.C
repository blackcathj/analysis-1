#ifndef MACRO_FUN4ALL_READDST_C
#define MACRO_FUN4ALL_READDST_C

#include <G4_Input.C>
#include <GlobalVariables.C>

#include <eidml/eIDMLInterface.h>

#include <TROOT.h>
#include <fun4all/Fun4AllDstOutputManager.h>
#include <fun4all/Fun4AllOutputManager.h>
#include <fun4all/Fun4AllServer.h>

#include <phool/recoConsts.h>

R__LOAD_LIBRARY(libfun4all.so)
// Use libg4eicdst for campaign 2 DSTs
R__LOAD_LIBRARY(libg4eicdst.so)
// Use libg4dst for campaign 1 DSTs
//R__LOAD_LIBRARY(libg4dst.so)
R__LOAD_LIBRARY(libeidml.so)

int Fun4All_ReadDST_eIDML(const int nEvents = 2,
                          //                          const string &inputFile = "singleElectron.lst"  //
                          const string &inputFile = "singleMuonPlus.prop.7.lst"  //
)

{
  //---------------
  // Fun4All server
  //---------------
  Fun4AllServer *se = Fun4AllServer::instance();
  se->Verbosity(1);

  // just if we set some flags somewhere in this macro
  recoConsts *rc = recoConsts::instance();

  Input::READHITS = true;
//  INPUTREADHITS::filename[0] = inputFile;
  INPUTREADHITS::listfile[0] = inputFile;

//  {
//    eIDMLInterface *anaTutorial = new eIDMLInterface(
//        {"BECAL", "HCALIN", "HCALOUT"},
//        inputFile + "_BECAL_" + to_string(nEvents) + ".root");
////    anaTutorial->Verbosity(1);
//    anaTutorial->setEtaRange(-2, 2);
//    se->registerSubsystem(anaTutorial);
//  }

//  {
//    eIDMLInterface *anaTutorial = new eIDMLInterface({"EEMC"}, inputFile + "_EEMC_" + to_string(nEvents) + ".root");
////    anaTutorial->Verbosity(1);
//    anaTutorial->setEtaRange(-4, 1);
//    se->registerSubsystem(anaTutorial);
//  }
//
  {
    eIDMLInterface *anaTutorial = new eIDMLInterface({"FEMC"
      ,"LFHCAL_0"
      ,"LFHCAL_1"
      ,"LFHCAL_2"
      ,"LFHCAL_3"
      ,"LFHCAL_4"
      ,"LFHCAL_5"
      ,"LFHCAL_6"
    }, inputFile + "_FEMC_" + to_string(nEvents) + ".root");
    anaTutorial->Verbosity(2);
    anaTutorial->setEtaRange(1, 4);
    se->registerSubsystem(anaTutorial);
  }

  InputManagers();

  se->run(nEvents);

  se->End();

  delete se;
  std::cout << "All done processing" << std::endl;
  gSystem->Exit(0);
  return 0;
}
#endif
