// $Id: $

/*!
 * \file Fun4All_FileSplit.C
 * \brief 
 * \author Jin Huang <jhuang@bnl.gov>
 * \version $Revision:   $
 * \date $Date: $
 */

#include <fun4all/Fun4AllDummyInputManager.h>
#include <fun4all/Fun4AllInputManager.h>
#include <fun4all/Fun4AllServer.h>
#include <g4main/ReadEICFiles.h>
#include <phhepmc/Fun4AllHepMCOutputManager.h>

R__LOAD_LIBRARY(libfun4all.so)
R__LOAD_LIBRARY(libg4testbench.so)

void Fun4All_FileSplit(const int nEvents = 2000,
                       const int nSkip = 2000,
                       const std::string &InputFile = "ep_noradcor.18x100_run001.root",
                       const std::string &OutputFile = "ep_noradcor.18x100_run001.HepMC.dat")
{
  //---------------
  // Fun4All server
  //---------------
  Fun4AllServer *se = Fun4AllServer::instance();
  se->Verbosity(0);

  ReadEICFiles *EICFileReader = new ReadEICFiles();
  EICFileReader->OpenInputFile(InputFile);
  EICFileReader->SetFirstEntry(nSkip);
  EICFileReader->set_embedding_id(0);
  se->registerSubsystem(EICFileReader);

  Fun4AllInputManager *in = new Fun4AllDummyInputManager("JADE");
  se->registerInputManager(in);

  Fun4AllHepMCOutputManager *out = new Fun4AllHepMCOutputManager("HepMCOut", OutputFile);
  se->registerOutputManager(out);

  se->run(nEvents);
}
