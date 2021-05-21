#include "MyNewAnalysis.h"

#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/PHCompositeNode.h>
#include <phool/getClass.h>

MyNewAnalysis::MyNewAnalysis(const std::string &name)
  : SubsysReco(name)
{}


MyNewAnalysis::~MyNewAnalysis()
{

}


int MyNewAnalysis::Init(PHCompositeNode *topNode)
{

  return Fun4AllReturnCodes::EVENT_OK;
}

int MyNewAnalysis::process_event(PHCompositeNode *topNode)
{
 
  
  return Fun4AllReturnCodes::EVENT_OK;
}

int MyNewAnalysis::End(PHCompositeNode *topNode)
{

  return Fun4AllReturnCodes::EVENT_OK;
}
