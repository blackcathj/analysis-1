#ifndef MACRO_FUN4ALL_READDST_C
#define MACRO_FUN4ALL_READDST_C


#include <GlobalVariables.C>
#include <G4_Input.C>

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



int Fun4All_ReadDST_eIDML(const int nEvents = 2000000,
    const string& inputFile = "singleElectron.lst"//
//        const string& inputFile = "singlePion.lst"//
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
  INPUTREADHITS::listfile[0] = inputFile;

  eIDMLInterface *anaTutorial = new eIDMLInterface("BECAL", inputFile + "_BECAL.root");
//  anaTutorial->setMinJetPt(3.);
  anaTutorial->Verbosity(1);
//  anaTutorial->analyzeTracks(true);
//  anaTutorial->analyzeClusters(true);
//  anaTutorial->analyzeJets(true);
//  anaTutorial->analyzeTruth(false);
  se->registerSubsystem(anaTutorial);

  InputManagers();

  se->run(nEvents);

  se->End();
  
  delete se;
  std::cout << "All done processing" << std::endl;
  gSystem->Exit(0);
  return 0;

}
#endif
