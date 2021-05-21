#ifndef MYNEWANALYSIS_H__
#define MYNEWANALYSIS_H__

#include <fun4all/SubsysReco.h>

/// Definition of this analysis module class
class MyNewAnalysis : public SubsysReco
{
 public:
  /// Constructor
  MyNewAnalysis(const std::string &name = "MyNewAnalysis");

  // Destructor
  virtual ~MyNewAnalysis();

  /// SubsysReco initialize processing method
  int Init(PHCompositeNode *);

  /// SubsysReco event processing method
  int process_event(PHCompositeNode *);

  /// SubsysReco end processing method
  int End(PHCompositeNode *);

 private:
  
};

#endif
