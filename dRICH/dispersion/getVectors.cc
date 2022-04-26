#include "getVectors.h"

#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/PHCompositeNode.h>
#include <phool/getClass.h>

#include <CLHEP/Vector/ThreeVector.h>

#include <math.h>
#include <float.h>
#include <numeric>

#include <g4main/PHG4HitContainer.h>
//____________________________________________________________________________..
getVectors::getVectors(const std::string &name)
    : SubsysReco(name)
    , m_truth_info(nullptr)
    , m_g4particle(nullptr)
    , m_outfile_name("outputData.root")
    , m_outfile(nullptr)
    , m_tree(nullptr)
{}

//____________________________________________________________________________..
getVectors::~getVectors() {}

//____________________________________________________________________________..
int getVectors::Init(PHCompositeNode *topNode)
{
  assert(topNode);
  initializeBranches();
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int getVectors::process_event(PHCompositeNode *topNode)
{
  m_truth_info = findNode::getClass<PHG4TruthInfoContainer>(topNode, "G4TruthInfo");
  if (!m_truth_info) 
  {
    std::cout << "DecayFinder: Missing node G4TruthInfo" << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  CLHEP::Hep3Vector *momentum_hit = NULL;
  CLHEP::Hep3Vector *the_first_hit = NULL;

  PHG4HitContainer* hit_container = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_MySuperDetector");
  if (!hit_container)
  {
    std::cout << __FILE__ << ": Missing node G4TruthInfo" << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  PHG4TruthInfoContainer::ConstRange range = m_truth_info->GetParticleRange();
  for (PHG4TruthInfoContainer::ConstIterator iter = range.first;
       iter != range.second; ++iter) 
  {
    m_g4particle = iter->second;
    if (m_g4particle->get_parent_id() == 0)
    {
      const CLHEP::Hep3Vector momentum_particle( m_g4particle->get_px()
                                               , m_g4particle->get_py()
                                               , m_g4particle->get_pz());
  
      m_pdg_id = m_g4particle->get_pid();
      m_track_id = m_g4particle->get_track_id();
      m_barcode = m_g4particle->get_barcode();
      m_particle_eta = momentum_particle.eta();
      m_particle_px = momentum_particle.x();
      m_particle_py = momentum_particle.y();
      m_particle_pz = momentum_particle.z();
      m_particle_pt = momentum_particle.perp();
      m_particle_p = momentum_particle.mag();
      m_particle_phi = momentum_particle.phi();

      std::vector<double> delta_phi;
      //Has to be done this way as truth eval functions are sPHENIX specific
      PHG4HitContainer::ConstRange hit_range = hit_container->getHits();
      double shortest_distance = DBL_MAX;
      for (PHG4HitContainer::ConstIterator hit_iter = hit_range.first;
           hit_iter != hit_range.second; ++hit_iter)
      {
        if (hit_iter->second->get_trkid() == m_g4particle->get_track_id())
        {
          double hit_distance = sqrt(pow(hit_iter->second->get_x(0), 2) + pow(hit_iter->second->get_y(0), 2) + pow(hit_iter->second->get_z(0), 2));
          if (hit_distance < shortest_distance)
          {   
            the_first_hit = new CLHEP::Hep3Vector( hit_iter->second->get_px(0)
                                                 , hit_iter->second->get_py(0)
                                                 , hit_iter->second->get_pz(0));
            shortest_distance = hit_distance;
          }
        }
      }
      for (PHG4HitContainer::ConstIterator hit_iter = hit_range.first;
           hit_iter != hit_range.second; ++hit_iter)
      {
        if (hit_iter->second->get_trkid() == m_g4particle->get_track_id())
        {
            
           momentum_hit = new CLHEP::Hep3Vector( hit_iter->second->get_px(0)
                                               , hit_iter->second->get_py(0)
                                               , hit_iter->second->get_pz(0));

           const CLHEP::Hep3Vector momentum_first_hit = (*the_first_hit);
           double hit_diff = momentum_hit->angle(momentum_first_hit);
           //double hit_diff = momentum_hit->angle(momentum_particle);
           delta_phi.push_back(hit_diff);
         }
      }
  
      double sum = std::accumulate(delta_phi.begin(), delta_phi.end(), 0.0);
      double mean = sum / delta_phi.size();
      double square_sum = std::inner_product(delta_phi.begin(), delta_phi.end(), delta_phi.begin(), 0.0);
      m_delta_phi = sqrt(square_sum / delta_phi.size());
      m_std_dev = sqrt(square_sum / delta_phi.size() - mean * mean);
  
      m_tree->Fill();
      delta_phi.clear();
    }
  }

  ++m_event_number;

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int getVectors::End(PHCompositeNode *topNode) 
{
  m_outfile->Write();
  m_outfile->Close();
  delete m_outfile;

  assert(topNode);

  return Fun4AllReturnCodes::EVENT_OK;
}

void getVectors::initializeBranches()
{
  m_outfile = new TFile(m_outfile_name.c_str(), "RECREATE");
  delete m_tree;
  m_tree = new TTree("dRICH", "dRICH");
  m_tree->OptimizeBaskets();
  m_tree->SetAutoSave(-5e6); // Save the output file every 5MB

  m_tree->Branch("EventNumber", &m_event_number, "EventNumber/I");
  m_tree->Branch("particle_PDG_ID", &m_pdg_id, "particle_PDG_ID/I");
  m_tree->Branch("particle_track_ID", &m_track_id, "particle_track_ID/I");
  m_tree->Branch("particle_barcode", &m_barcode, "particle_barcode/I");
  m_tree->Branch("particle_eta", &m_particle_eta, "particle_eta/F");
  m_tree->Branch("particle_px", &m_particle_px, "particle_px/F");
  m_tree->Branch("particle_py", &m_particle_py, "particle_py/F");
  m_tree->Branch("particle_pz", &m_particle_pz, "particle_pz/F");
  m_tree->Branch("particle_pt", &m_particle_pt, "particle_pt/F");
  m_tree->Branch("particle_p", &m_particle_p, "particle_p/F");
  m_tree->Branch("particle_phi", &m_particle_phi, "particle_phi/F");
  m_tree->Branch("delta_phi", &m_delta_phi, "delta_phi/F");
  m_tree->Branch("std_dev", &m_std_dev, "std_dev/F");
}
