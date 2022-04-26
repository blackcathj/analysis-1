// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef GETVECTORS_H
#define GETVECTORS_H

#include <fun4all/SubsysReco.h>

#include <g4main/PHG4Hit.h>
#include <g4main/PHG4HitContainer.h>
#include <g4main/PHG4Particle.h>
#include <g4main/PHG4TruthInfoContainer.h>

#include <TBranch.h>
#include <TFile.h>
#include <TTree.h>

#include <string>

class PHCompositeNode;
class GlobalVertexMap;
class PHG4HitContainer;
class PHG4Hit;
class PHG4TruthInfoContainer;
class PHG4Particle;

namespace CLHEP {
class Hep3Vector;
}

class getVectors : public SubsysReco {
public:
  getVectors(const std::string &name = "getVectors");

  virtual ~getVectors();

  /** Called during initialization.
      Typically this is where you can book histograms, and e.g.
      register them to Fun4AllServer (so they can be output to file
      using Fun4AllServer::dumpHistos() method).
   */
  int Init(PHCompositeNode *topNode) override;

  /** Called for each event.
      This is where you do the real work.
   */
  int process_event(PHCompositeNode *topNode) override;

  /// Called at the end of all processing.
  int End(PHCompositeNode *topNode) override;

  void setOutputName(const std::string &name) { m_outfile_name = name; }

private:
  PHG4HitContainer *m_hit_container;
  PHG4TruthInfoContainer *m_truth_info;
  PHG4Particle *m_g4particle;

  void initializeBranches();

  std::string m_outfile_name;
  TFile *m_outfile;
  TTree *m_tree;

  unsigned int m_event_number = 0;
  int m_pdg_id = 0;
  int m_track_id = 0;
  int m_barcode = 0;
  float m_particle_px = -99;
  float m_particle_py = -99;
  float m_particle_pz = -99;
  float m_particle_pt = -99;
  float m_particle_p = -99;
  float m_particle_eta = -99;
  float m_particle_phi = -99;

  float m_delta_phi = -99;
  float m_std_dev = -99;
};

#endif // GETVECTORS_H
