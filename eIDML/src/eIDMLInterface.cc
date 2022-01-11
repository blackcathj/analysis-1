#include "eIDMLInterface.h"

/// Cluster/Calorimeter includes
#include <calobase/RawCluster.h>
#include <calobase/RawClusterContainer.h>
#include <calobase/RawClusterUtility.h>
#include <calobase/RawTower.h>
#include <calobase/RawTowerContainer.h>
#include <calobase/RawTowerGeom.h>
#include <calobase/RawTowerGeomContainer.h>
#include <calotrigger/CaloTriggerInfo.h>

/// Jet includes
#include <g4jets/Jet.h>
#include <g4jets/JetMap.h>

/// Tracking includes
#include <calobase/RawCluster.h>
#include <calobase/RawClusterContainer.h>
#include <calobase/RawTower.h>
#include <calobase/RawTowerContainer.h>
#include <calobase/RawTowerGeomContainer.h>
#include <g4vertex/GlobalVertex.h>
#include <g4vertex/GlobalVertexMap.h>
#include <trackbase_historic/SvtxTrack.h>
#include <trackbase_historic/SvtxTrackMap.h>
#include <trackbase_historic/SvtxTrack_FastSim.h>
#include <trackbase_historic/SvtxVertex.h>
#include <trackbase_historic/SvtxVertexMap.h>
//PID includes
#include <eicpidbase/EICPIDParticle.h>
#include <eicpidbase/EICPIDParticleContainer.h>
/// Truth evaluation includes
#include <g4eval/JetEvalStack.h>
#include <g4eval/SvtxEvalStack.h>

/// HEPMC truth includes
#include <HepMC/GenEvent.h>
#include <HepMC/GenVertex.h>
#include <phhepmc/PHHepMCGenEvent.h>
#include <phhepmc/PHHepMCGenEventMap.h>

/// Fun4All includes
#include <fun4all/Fun4AllHistoManager.h>
#include <fun4all/Fun4AllReturnCodes.h>
#include <g4main/PHG4Hit.h>
#include <g4main/PHG4Particle.h>
#include <g4main/PHG4TruthInfoContainer.h>
#include <phool/PHCompositeNode.h>
#include <phool/getClass.h>

/// ROOT includes
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TMath.h>
#include <TNtuple.h>
#include <TTree.h>
#include <TVector3.h>

/// C++ includes
#include <cassert>
#include <cmath>
#include <sstream>
#include <string>

using namespace std;

/**
 * This class demonstrates the construction and use of an analysis module 
 * within the EIC Fun4All framework. It is intended to show how to 
 * obtain physics objects from the analysis tree, and then write them out
 * to a ROOT tree and file for further offline analysis.  
 */

/**
 * Constructor of module
 */
eIDMLInterface::eIDMLInterface(const std::vector<std::string> &names, const std::string &filename)
  : SubsysReco("eIDMLInterface_" + names[0])
  , _calo_names(names)
  //  , m_TTree_Tower_dEta(m_sizeTowerPatch * m_sizeTowerPatch, 0)
  //  , m_TTree_Tower_dPhi(m_sizeTowerPatch * m_sizeTowerPatch, 0)
  //  , m_TTree_Tower_iEta_patch(m_sizeTowerPatch * m_sizeTowerPatch, 0)
  //  , m_TTree_Tower_iPhi_patch(m_sizeTowerPatch * m_sizeTowerPatch, 0)
  //  , m_TTree_Tower_E(m_sizeTowerPatch * m_sizeTowerPatch, 0)
  , m_outfilename(filename)
  , m_hm(nullptr)
  , m_minjetpt(5.0)
  , m_mincluspt(0.25)
  , m_analyzeTracks(true)
  , m_analyzeClusters(true)
  , m_analyzeJets(true)
  , m_analyzeTruth(false)
{
  /// Initialize variables and trees so we don't accidentally access
  /// memory that was never allocated
  initializeVariables();
  initializeTrees();
}

/**
 * Destructor of module
 */
eIDMLInterface::~eIDMLInterface()
{
}

/**
 * Initialize the module and prepare looping over events
 */
int eIDMLInterface::Init(PHCompositeNode *topNode)
{
  if (Verbosity() > 5)
  {
    cout << "Beginning Init in eIDMLInterface" << endl;
  }

  m_outfile = new TFile(m_outfilename.c_str(), "RECREATE");

  m_phi_h = new TH1D("phi_h", ";Counts;#phi [rad]", 50, -6, 6);
  m_eta_phi_h = new TH2F("phi_eta_h", ";#eta;#phi [rad]", 10, -1, 1, 50, -6, 6);

  return 0;
}

/**
 * Main workhorse function where each event is looped over and 
 * data from each event is collected from the node tree for analysis
 */
int eIDMLInterface::process_event(PHCompositeNode *topNode)
{
  if (Verbosity() > 5)
  {
    cout << "Beginning process_event in eIDMLInterface" << endl;
  }

  /// G4 truth particle node
  PHG4TruthInfoContainer *truthinfo = findNode::getClass<PHG4TruthInfoContainer>(topNode, "G4TruthInfo");
  SvtxTrackMap *trackmap = findNode::getClass<SvtxTrackMap>(topNode, "TrackMap");

  if (!truthinfo)
  {
    cout << PHWHERE
         << "PHG4TruthInfoContainer node is missing, can't collect G4 truth particles"
         << endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  /// Get the primary particle range
  PHG4TruthInfoContainer::Range range = truthinfo->GetPrimaryParticleRange();
  /// Loop over the G4 truth (stable) particles
  for (PHG4TruthInfoContainer::ConstIterator iter = range.first;
       iter != range.second;
       ++iter)
  {
    initializeVariables();

    /// Get this truth particle
    const PHG4Particle *truth = iter->second;

    /// Get this particles momentum, etc.
    m_truthpx = truth->get_px();
    m_truthpy = truth->get_py();
    m_truthpz = truth->get_pz();
    m_truthp = sqrt(m_truthpx * m_truthpx + m_truthpy * m_truthpy + m_truthpz * m_truthpz);
    m_truthenergy = truth->get_e();

    m_truthpt = sqrt(m_truthpx * m_truthpx + m_truthpy * m_truthpy);

    m_truthphi = atan2(m_truthpy, m_truthpx);

    m_trutheta = atanh(m_truthpz / m_truthenergy);

    // eta cut
    if (m_trutheta < m_etaRange.first or m_trutheta > m_etaRange.second) continue;

    /// Check for nans
    if (m_trutheta != m_trutheta)
      m_trutheta = -99;
    m_truthpid = truth->get_pid();

    SvtxTrack_FastSim *track = nullptr;

    if (Verbosity() > 1) cout << __PRETTY_FUNCTION__ << "TRACKmap size " << trackmap->size() << std::endl;
    for (SvtxTrackMap::ConstIter track_itr = trackmap->begin();
         track_itr != trackmap->end();
         track_itr++)
    {
      //std::cout << "TRACK * " << track_itr->first << std::endl;
      SvtxTrack_FastSim *temp = dynamic_cast<SvtxTrack_FastSim *>(track_itr->second);
      if (!temp)
      {
        if (Verbosity() > 1)
        {
          std::cout << "PHG4TrackFastSimEval::fill_track_tree - ignore track that is not a SvtxTrack_FastSim:";
          track_itr->second->identify();
        }
        continue;
      }
      if (Verbosity() > 1) cout << __PRETTY_FUNCTION__ << " PARTICLE!" << std::endl;

      if ((temp->get_truth_track_id() - truth->get_track_id()) == 0)
      {
        track = temp;

        break;
      }
    }

    if (track)
    {
      // track matched

      /// Get the reconstructed track info
      m_tr_px = track->get_px();
      m_tr_py = track->get_py();
      m_tr_pz = track->get_pz();
      m_tr_p = sqrt(m_tr_px * m_tr_px + m_tr_py * m_tr_py + m_tr_pz * m_tr_pz);

      m_tr_pt = sqrt(m_tr_px * m_tr_px + m_tr_py * m_tr_py);

      m_tr_phi = track->get_phi();
      m_tr_eta = track->get_eta();

      m_charge = track->get_charge();
      m_chisq = track->get_chisq();
      m_ndf = track->get_ndf();
      m_dca = track->get_dca();
      m_tr_x = track->get_x();
      m_tr_y = track->get_y();
      m_tr_z = track->get_z();

      for (const std::string detector : _calo_names)
      {
        //  const std::string detector(_calo_name);

        std::string towernodename = "TOWER_CALIB_" + detector;
        // Grab the towers
        RawTowerContainer *towers = findNode::getClass<RawTowerContainer>(topNode,
                                                                          towernodename);
        if (!towers)
        {
          std::cout << PHWHERE << ": Could not find node " << towernodename
                    << std::endl;
          return Fun4AllReturnCodes::ABORTRUN;
        }
        std::string towergeomnodename = "TOWERGEOM_" + detector;
        RawTowerGeomContainer *towergeom = findNode::getClass<RawTowerGeomContainer>(
            topNode, towergeomnodename);
        if (!towergeom)
        {
          std::cout << PHWHERE << ": Could not find node " << towergeomnodename
                    << std::endl;
          return Fun4AllReturnCodes::ABORTRUN;
        }

        bool has_projection(false);
        // find projections
        for (SvtxTrack::ConstStateIter trkstates = track->begin_states();
             trkstates != track->end_states();
             ++trkstates)
        {
          if (Verbosity() > 1) cout << __PRETTY_FUNCTION__ << " checking " << trkstates->second->get_name() << endl;
          if (trkstates->second->get_name().compare(detector) == 0)
          {
            if (Verbosity() > 1) cout << __PRETTY_FUNCTION__ << " found " << trkstates->second->get_name() << endl;
            has_projection = true;

            // setting the projection (xyz and pxpypz)
            for (int i = 0; i < 3; i++)
            {
              m_CaloDataMap[detector].m_TTree_proj_vec[i] = trkstates->second->get_pos(i);
              m_CaloDataMap[detector].m_TTree_proj_p_vec[i] = trkstates->second->get_mom(i);
            }
            // fourth element is the path length
            m_CaloDataMap[detector].m_TTree_proj_vec[3] = trkstates->first;
          }
        }  //       for (SvtxTrack::ConstStateIter trkstates = track->begin_states();

        // projection match to calorimeter
        if (has_projection)
        {
          TVector3 vec_proj(m_CaloDataMap[detector].m_TTree_proj_vec[0],
                            m_CaloDataMap[detector].m_TTree_proj_vec[1],
                            m_CaloDataMap[detector].m_TTree_proj_vec[2]);

          const double eta_proj = vec_proj.Eta();
          const double phi_proj = vec_proj.Phi();

          double min_tower_r2 = 1000;
          RawTowerDefs::keytype central_tower_key = -1;
          const RawTowerGeom *central_tower(nullptr);
          int maxBinPhi = 0;
          int minBinPhi = 1000;

          RawTowerGeomContainer::ConstRange range = towergeom->get_tower_geometries();
          /// Loop over the G4 truth (stable) particles
          for (RawTowerGeomContainer::ConstIterator titer = range.first;
               titer != range.second;
               ++titer)
          {
            //          const int bineta = RawTowerDefs::decode_index1(titer->first);
            const int binphi = RawTowerDefs::decode_index2(titer->first);
            if (maxBinPhi < binphi) maxBinPhi = binphi;
            if (minBinPhi > binphi) minBinPhi = binphi;

            const RawTowerGeom *tower_geom = titer->second;
            assert(tower_geom);

            TVector3 vec_tower(
                tower_geom->get_center_x(),
                tower_geom->get_center_y(),
                tower_geom->get_center_z());

            const double deta = eta_proj - vec_tower.Eta();
            double dphi = phi_proj - vec_tower.Phi();

            if (dphi > M_PI) dphi -= 2 * M_PI;
            if (dphi < -M_PI) dphi += 2 * M_PI;

            const double r2 = dphi * dphi + deta * deta;

            if (r2 < min_tower_r2)
            {
              min_tower_r2 = r2;
              central_tower_key = titer->first;
              central_tower = titer->second;
            }

          }  //         for (RawTowerGeomContainer::ConstIterator titer = range.first;

          if (central_tower == nullptr) continue;
          if (central_tower_key < 0) continue;

          if (Verbosity() > 1)
          {
            cout << __PRETTY_FUNCTION__ << " found tower " << central_tower_key << ": ";
            cout << " min_tower_r2 =  " << min_tower_r2;
            cout << " decode_index1 =  " << RawTowerDefs::decode_index1(central_tower_key);
            cout << " decode_index2 =  " << RawTowerDefs::decode_index2(central_tower_key);
            cout << " minBinPhi =  " << minBinPhi;
            cout << " maxBinPhi =  " << maxBinPhi;
            central_tower->identify();
          }

          // print tower patch
          const int central_tower_eta = RawTowerDefs::decode_index1(central_tower_key);
          const int central_tower_phi = RawTowerDefs::decode_index2(central_tower_key);
          const int size_half_tower_patch = (m_CaloDataMap[detector].m_sizeTowerPatch - 1) / 2;  // 7x7
          size_t tower_index_patch = 0;

          m_CaloDataMap[detector].m_centralTowerBinEta = central_tower_eta;
          m_CaloDataMap[detector].m_centralTowerBinPhi = central_tower_phi;

          for (int ieta_patch = -size_half_tower_patch;
               ieta_patch <= +size_half_tower_patch;
               ++ieta_patch)
          {
            const int bin_eta = central_tower_eta + ieta_patch;

            for (int iphi_patch = -size_half_tower_patch;
                 iphi_patch <= +size_half_tower_patch;
                 ++iphi_patch)
            {
              assert(tower_index_patch < m_CaloDataMap[detector].m_TTree_Tower_E.size());

              int bin_phi = central_tower_phi + iphi_patch;

              if (bin_phi < minBinPhi) bin_phi = bin_phi - minBinPhi + maxBinPhi + 1;
              if (bin_phi > maxBinPhi) bin_phi = bin_phi - maxBinPhi + minBinPhi - 1;

              m_CaloDataMap[detector].m_TTree_Tower_iEta_patch[tower_index_patch] = ieta_patch;
              m_CaloDataMap[detector].m_TTree_Tower_iPhi_patch[tower_index_patch] = iphi_patch;

              if (bin_eta > 4095 or bin_phi > 4095 or bin_eta < 0 or bin_phi < 0)
              {
                if (Verbosity())
                {
                  cout << __PRETTY_FUNCTION__ << " invalid tower geom " << central_tower_key << ": ";
                  cout << " bin_eta =  " << bin_eta;
                  cout << " bin_phi =  " << bin_phi;
                  cout << " central_tower_eta =  " << central_tower_eta;
                  cout << " central_tower_phi =  " << central_tower_phi;
                  cout << " central_tower_key =  " << central_tower_key;
                  cout << " min_tower_r2 =  " << min_tower_r2;
                  cout << " decode_index1 =  " << RawTowerDefs::decode_index1(central_tower_key);
                  cout << " decode_index2 =  " << RawTowerDefs::decode_index2(central_tower_key);
                  cout << " minBinPhi =  " << minBinPhi;
                  cout << " maxBinPhi =  " << maxBinPhi;
                  central_tower->identify();
                }

                return Fun4AllReturnCodes::ABORTEVENT;
                //              continue;
              }
              RawTowerDefs::keytype tower_key = RawTowerDefs::encode_towerid(
                  towergeom->get_calorimeter_id(), bin_eta, bin_phi);
              const RawTowerGeom *tower_geom = towergeom->get_tower_geometry(tower_key);
              if (tower_geom)
              {
                TVector3 vec_tower(
                    tower_geom->get_center_x(),
                    tower_geom->get_center_y(),
                    tower_geom->get_center_z());

                const double deta = eta_proj - vec_tower.Eta();
                double dphi = phi_proj - vec_tower.Phi();

                m_CaloDataMap[detector].m_TTree_Tower_dEta[tower_index_patch] = deta;
                m_CaloDataMap[detector].m_TTree_Tower_dPhi[tower_index_patch] = dphi;

                if (Verbosity() > 2)
                {
                  cout << __PRETTY_FUNCTION__ << " process tower geom " << tower_key << ": ";
                  cout << " ieta_patch =  " << ieta_patch;
                  cout << " iphi_patch =  " << iphi_patch;
                  cout << " bin_eta =  " << bin_eta;
                  cout << " bin_phi =  " << bin_phi;
                  cout << " deta =  " << deta;
                  cout << " dphi =  " << dphi;
                  tower_geom->identify();
                }

                const RawTower *tower = towers->getTower(tower_key);

                if (tower)
                {
                  const double energy = tower->get_energy();

                  if (abs(ieta_patch) <= 1 and abs(iphi_patch) <= 1) m_CaloDataMap[detector].m_E3x3 += energy;
                  if (abs(ieta_patch) <= 2 and abs(iphi_patch) <= 2) m_CaloDataMap[detector].m_E5x5 += energy;
                  if (abs(ieta_patch) <= 3 and abs(iphi_patch) <= 3) m_CaloDataMap[detector].m_E7x7 += energy;

                  m_CaloDataMap[detector].m_TTree_Tower_E[tower_index_patch] = energy;

                  if (Verbosity() > 2)
                  {
                    cout << __PRETTY_FUNCTION__ << " process tower " << tower_key << ": ";
                    cout << " ieta_patch =  " << ieta_patch;
                    cout << " iphi_patch =  " << iphi_patch;
                    cout << " bin_eta =  " << bin_eta;
                    cout << " bin_phi =  " << bin_phi;
                    cout << " energy =  " << energy;
                    tower->identify();
                  }

                }  //               if (tower)

                ++tower_index_patch;
              }  //             if (tower_geom)

            }  //           for (int iphi_patch = central_tower_phi - size_half_tower_patch;

          }  //         for (int ieta_patch = central_tower_eta - size_half_tower_patch;

          //        const int bineta = RawTowerDefs::decode_index1(central_tower_key);
          //        const int binphi = RawTowerDefs::decode_index2(central_tower_key);

        }  //       if (has_projection)

      }  // iterate detectors

    }  //     if (track)

    /// Fill the g4 truth tree
    assert(m_truthtree);
    m_truthtree->Fill();
  }  /// Loop over the G4 truth (stable) particles

  return Fun4AllReturnCodes::EVENT_OK;
}

/**
 * End the module and finish any data collection. Clean up any remaining
 * loose ends
 */
int eIDMLInterface::End(PHCompositeNode *topNode)
{
  if (Verbosity() > 1)
  {
    cout << "Ending eIDMLInterface analysis package" << endl;
  }

  /// Change to the outfile
  m_outfile->cd();
  m_truthtree->Write();
  /// Write out any other histograms
  m_phi_h->Write();
  m_eta_phi_h->Write();

  /// Write and close the outfile
  m_outfile->Write();
  m_outfile->Close();

  /// Clean up pointers and associated histos/trees in TFile
  delete m_outfile;

  if (Verbosity() > 1)
  {
    cout << "Finished eIDMLInterface analysis package" << endl;
  }

  return 0;
}

/**
 * This method gets all of the HEPMC truth particles from the node tree
 * and stores them in a ROOT TTree. The HEPMC truth particles are what, 
 * for example, directly comes out of PYTHIA and thus gives you all of
 * the associated parton information
 */
void eIDMLInterface::getHEPMCTruth(PHCompositeNode *topNode)
{
  /// Get the node from the node tree
  PHHepMCGenEventMap *hepmceventmap = findNode::getClass<PHHepMCGenEventMap>(topNode, "PHHepMCGenEventMap");

  /// If the node was not properly put on the tree, return
  if (!hepmceventmap)
  {
    cout << PHWHERE
         << "HEPMC event map node is missing, can't collected HEPMC truth particles"
         << endl;
    return;
  }

  /// Could have some print statements for debugging with verbosity
  if (Verbosity() > 1)
  {
    cout << "Getting HEPMC truth particles " << endl;
  }

  /// You can iterate over the number of events in a hepmc event
  /// for pile up events where you have multiple hard scatterings per bunch crossing
  for (PHHepMCGenEventMap::ConstIter eventIter = hepmceventmap->begin();
       eventIter != hepmceventmap->end();
       ++eventIter)
  {
    /// Get the event
    PHHepMCGenEvent *hepmcevent = eventIter->second;

    if (hepmcevent)
    {
      /// Get the event characteristics, inherited from HepMC classes
      HepMC::GenEvent *truthevent = hepmcevent->getEvent();
      if (!truthevent)
      {
        cout << PHWHERE
             << "no evt pointer under phhepmvgeneventmap found "
             << endl;
        return;
      }

      /// Get the parton info
      HepMC::PdfInfo *pdfinfo = truthevent->pdf_info();

      /// Get the parton info as determined from HEPMC
      m_partid1 = pdfinfo->id1();
      m_partid2 = pdfinfo->id2();
      m_x1 = pdfinfo->x1();
      m_x2 = pdfinfo->x2();

      /// Are there multiple partonic intercations in a p+p event
      m_mpi = truthevent->mpi();

      /// Get the PYTHIA signal process id identifying the 2-to-2 hard process
      m_process_id = truthevent->signal_process_id();

      if (Verbosity() > 2)
      {
        cout << " Iterating over an event" << endl;
      }
      /// Loop over all the truth particles and get their information
      for (HepMC::GenEvent::particle_const_iterator iter = truthevent->particles_begin();
           iter != truthevent->particles_end();
           ++iter)
      {
        /// Get each pythia particle characteristics
        m_truthenergy = (*iter)->momentum().e();
        m_truthpid = (*iter)->pdg_id();

        m_trutheta = (*iter)->momentum().pseudoRapidity();
        m_truthphi = (*iter)->momentum().phi();
        m_truthpx = (*iter)->momentum().px();
        m_truthpy = (*iter)->momentum().py();
        m_truthpz = (*iter)->momentum().pz();
        m_truthpt = sqrt(m_truthpx * m_truthpx + m_truthpy * m_truthpy);

        /// Fill the truth tree
        m_hepmctree->Fill();
        m_numparticlesinevent++;
      }
    }
  }
}

/**
 * This function collects the truth PHG4 stable particles that
 * are produced from PYTHIA, or some other event generator. These
 * are the stable particles, e.g. there are not any (for example)
 * partons here.
 */
void eIDMLInterface::getPHG4Truth(PHCompositeNode *topNode)
{
  /// G4 truth particle node
  PHG4TruthInfoContainer *truthinfo = findNode::getClass<PHG4TruthInfoContainer>(topNode, "G4TruthInfo");

  if (!truthinfo)
  {
    cout << PHWHERE
         << "PHG4TruthInfoContainer node is missing, can't collect G4 truth particles"
         << endl;
    return;
  }

  /// Get the primary particle range
  PHG4TruthInfoContainer::Range range = truthinfo->GetPrimaryParticleRange();

  /// Loop over the G4 truth (stable) particles
  for (PHG4TruthInfoContainer::ConstIterator iter = range.first;
       iter != range.second;
       ++iter)
  {
    /// Get this truth particle
    const PHG4Particle *truth = iter->second;

    /// Get this particles momentum, etc.
    m_truthpx = truth->get_px();
    m_truthpy = truth->get_py();
    m_truthpz = truth->get_pz();
    m_truthp = sqrt(m_truthpx * m_truthpx + m_truthpy * m_truthpy + m_truthpz * m_truthpz);
    m_truthenergy = truth->get_e();

    m_truthpt = sqrt(m_truthpx * m_truthpx + m_truthpy * m_truthpy);

    m_truthphi = atan2(m_truthpy, m_truthpx);

    m_trutheta = atanh(m_truthpz / m_truthenergy);
    /// Check for nans
    if (m_trutheta != m_trutheta)
      m_trutheta = -99;
    m_truthpid = truth->get_pid();

    /// Fill the g4 truth tree
    m_truthtree->Fill();
  }
}

/**
 * This method gets the tracks as reconstructed from the tracker. It also
 * compares the reconstructed track to its truth track.
 */
void eIDMLInterface::getTracks(PHCompositeNode *topNode)
{
  /// Tracks node
  SvtxTrackMap *trackmap = findNode::getClass<SvtxTrackMap>(topNode, "TrackMap");
  EICPIDParticleContainer *pidcontainer = findNode::getClass<EICPIDParticleContainer>(topNode, "EICPIDParticleMap");

  if (Verbosity() > 1 and pidcontainer == nullptr)
  {
    cout << "EICPIDParticleContainer named EICPIDParticleMap does not exist. Skip saving PID info" << endl;
  }

  if (!trackmap)
  {
    cout << PHWHERE
         << "TrackMap node is missing, can't collect tracks"
         << endl;
    return;
  }

  /// Get the range for primary tracks
  PHG4TruthInfoContainer *truthinfo = findNode::getClass<PHG4TruthInfoContainer>(topNode, "G4TruthInfo");

  if (Verbosity() > 1)
  {
    cout << "Get the tracks" << endl;
  }
  for (SvtxTrackMap::Iter iter = trackmap->begin();
       iter != trackmap->end();
       ++iter)
  {
    SvtxTrack *track = iter->second;

    /// Get the reconstructed track info
    m_tr_px = track->get_px();
    m_tr_py = track->get_py();
    m_tr_pz = track->get_pz();
    m_tr_p = sqrt(m_tr_px * m_tr_px + m_tr_py * m_tr_py + m_tr_pz * m_tr_pz);

    m_tr_pt = sqrt(m_tr_px * m_tr_px + m_tr_py * m_tr_py);

    // Make some cuts on the track to clean up sample
    if (m_tr_pt < 0.5)
      continue;

    m_tr_phi = track->get_phi();
    m_tr_eta = track->get_eta();

    m_charge = track->get_charge();
    m_chisq = track->get_chisq();
    m_ndf = track->get_ndf();
    m_dca = track->get_dca();
    m_tr_x = track->get_x();
    m_tr_y = track->get_y();
    m_tr_z = track->get_z();

    /// Ensure that the reco track is a fast sim track
    SvtxTrack_FastSim *temp = dynamic_cast<SvtxTrack_FastSim *>(iter->second);
    if (!temp)
    {
      if (Verbosity() > 0)
        std::cout << "Skipping non fast track sim object..." << std::endl;
      continue;
    }

    /// Get truth track info that matches this reconstructed track
    PHG4Particle *truthtrack = truthinfo->GetParticle(temp->get_truth_track_id());
    if (truthtrack)
    {
      m_truth_is_primary = truthinfo->is_primary(truthtrack);

      m_truthtrackpx = truthtrack->get_px();
      m_truthtrackpy = truthtrack->get_py();
      m_truthtrackpz = truthtrack->get_pz();
      m_truthtrackp = sqrt(m_truthtrackpx * m_truthtrackpx + m_truthtrackpy * m_truthtrackpy + m_truthtrackpz * m_truthtrackpz);

      m_truthtracke = truthtrack->get_e();

      m_truthtrackpt = sqrt(m_truthtrackpx * m_truthtrackpx + m_truthtrackpy * m_truthtrackpy);
      m_truthtrackphi = atan2(m_truthtrackpy, m_truthtrackpx);
      m_truthtracketa = atanh(m_truthtrackpz / m_truthtrackp);
      m_truthtrackpid = truthtrack->get_pid();
    }

    // match to PIDparticles
    if (pidcontainer)
    {
      // EICPIDParticle are index the same as the tracks
      const EICPIDParticle *pid_particle =
          pidcontainer->findEICPIDParticle(track->get_id());

      if (pid_particle)
      {
        // top level log likelihood sums.
        // More detailed per-detector information also available at  EICPIDParticle::get_LogLikelyhood(EICPIDDefs::PIDCandidate, EICPIDDefs::PIDDetector)
        m_tr_pion_loglikelihood = pid_particle->get_SumLogLikelyhood(EICPIDDefs::PionCandiate);
        m_tr_kaon_loglikelihood = pid_particle->get_SumLogLikelyhood(EICPIDDefs::KaonCandiate);
        m_tr_proton_loglikelihood = pid_particle->get_SumLogLikelyhood(EICPIDDefs::ProtonCandiate);
      }
    }

    m_tracktree->Fill();
  }
}

/**
 * Method that gets the truth jets and stores them in a tree
 */
void eIDMLInterface::getTruthJets(PHCompositeNode *topNode)
{
  if (Verbosity() > 1)
  {
    cout << "get the truth jets" << endl;
  }

  /// Get the truth jet node
  JetMap *truth_jets = findNode::getClass<JetMap>(topNode, "AntiKt_Truth_r05");

  /// Get reco jets associated to truth jets to study e.g. jet efficiencies
  JetMap *reco_jets = findNode::getClass<JetMap>(topNode, "AntiKt_Tower_r05");
  if (!m_jetEvalStack)
  {
    m_jetEvalStack = new JetEvalStack(topNode, "AntiKt_Tower_r05",
                                      "AntiKt_Truth_r05");
  }
  m_jetEvalStack->next_event(topNode);
  JetTruthEval *trutheval = m_jetEvalStack->get_truth_eval();

  if (!truth_jets)
  {
    cout << PHWHERE
         << "Truth jet node is missing, can't collect truth jets"
         << endl;
    return;
  }

  /// Iterate over the truth jets
  for (JetMap::Iter iter = truth_jets->begin();
       iter != truth_jets->end();
       ++iter)
  {
    Jet *truthJet = iter->second;

    m_truthjetpt = truthJet->get_pt();

    std::set<PHG4Particle *> truthjetcomp =
        trutheval->all_truth_particles(truthJet);
    int ntruthconstituents = 0;
    //loop over the constituents of the truth jet
    for (std::set<PHG4Particle *>::iterator iter2 = truthjetcomp.begin();
         iter2 != truthjetcomp.end();
         ++iter2)
    {
      //get the particle of the truthjet
      PHG4Particle *truthpart = *iter2;
      if (!truthpart)
      {
        cout << "no truth particles in the jet??" << endl;
        break;
      }

      ntruthconstituents++;
    }

    if (ntruthconstituents < 3)
      continue;
    /// Only collect truthjets above the _minjetpt cut
    if (m_truthjetpt < m_minjetpt)
      continue;

    m_truthjeteta = truthJet->get_eta();
    m_truthjetpx = truthJet->get_px();
    m_truthjetpy = truthJet->get_py();
    m_truthjetpz = truthJet->get_pz();
    m_truthjetphi = truthJet->get_phi();
    m_truthjetp = truthJet->get_p();
    m_truthjetenergy = truthJet->get_e();

    m_recojetpt = 0;
    m_recojetid = 0;
    m_recojetpx = 0;
    m_recojetpy = 0;
    m_recojetpz = 0;
    m_recojetphi = 0;
    m_recojetp = 0;
    m_recojetenergy = 0;
    m_dR = -99;
    float closestjet = 9999;
    /// Iterate over the reconstructed jets
    for (JetMap::Iter recoIter = reco_jets->begin();
         recoIter != reco_jets->end();
         ++recoIter)
    {
      const Jet *recoJet = recoIter->second;
      m_recojetpt = recoJet->get_pt();
      if (m_recojetpt < m_minjetpt - m_minjetpt * 0.4)
        continue;

      m_recojeteta = recoJet->get_eta();
      m_recojetphi = recoJet->get_phi();

      if (Verbosity() > 1)
      {
        cout << "matching by distance jet" << endl;
      }

      float dphi = m_recojetphi - m_truthjetphi;
      if (dphi > TMath::Pi())
        dphi -= TMath::Pi() * 2.;
      if (dphi < -1 * TMath::Pi())
        dphi += TMath::Pi() * 2.;

      float deta = m_recojeteta - m_truthjeteta;
      /// Determine the distance in eta phi space between the reconstructed
      /// and truth jets
      m_dR = sqrt(pow(dphi, 2.) + pow(deta, 2.));

      /// If this truth jet is closer than the previous truth jet, it is
      /// closer and thus should be considered the truth jet
      if (m_dR < truth_jets->get_par() && m_dR < closestjet)
      {
        // Get reco jet characteristics
        m_recojetid = recoJet->get_id();
        m_recojetpx = recoJet->get_px();
        m_recojetpy = recoJet->get_py();
        m_recojetpz = recoJet->get_pz();
        m_recojetphi = recoJet->get_phi();
        m_recojetp = recoJet->get_p();
        m_recojetenergy = recoJet->get_e();
      }
    }

    /// Fill the truthjet tree
    m_truthjettree->Fill();
  }
}

/**
 * Get the reconstructed jets and store them in a tree
 */
void eIDMLInterface::getReconstructedJets(PHCompositeNode *topNode)
{
  /// Get the reconstructed tower jets
  JetMap *reco_jets = findNode::getClass<JetMap>(topNode, "AntiKt_Tower_r05");
  /// Get the truth jets
  JetMap *truth_jets = findNode::getClass<JetMap>(topNode, "AntiKt_Truth_r05");

  if (!m_jetEvalStack)
  {
    m_jetEvalStack = new JetEvalStack(topNode, "AntiKt_Tower_r05",
                                      "AntiKt_Truth_r05");
  }
  m_jetEvalStack->next_event(topNode);
  JetRecoEval *recoeval = m_jetEvalStack->get_reco_eval();
  if (!reco_jets)
  {
    cout << PHWHERE
         << "Reconstructed jet node is missing, can't collect reconstructed jets"
         << endl;
    return;
  }

  if (Verbosity() > 1)
  {
    cout << "Get all Reco Jets" << endl;
  }

  /// Iterate over the reconstructed jets
  for (JetMap::Iter recoIter = reco_jets->begin();
       recoIter != reco_jets->end();
       ++recoIter)
  {
    Jet *recoJet = recoIter->second;
    m_recojetpt = recoJet->get_pt();
    if (m_recojetpt < m_minjetpt)
      continue;

    m_recojeteta = recoJet->get_eta();

    // Get reco jet characteristics
    m_recojetid = recoJet->get_id();
    m_recojetpx = recoJet->get_px();
    m_recojetpy = recoJet->get_py();
    m_recojetpz = recoJet->get_pz();
    m_recojetphi = recoJet->get_phi();
    m_recojetp = recoJet->get_p();
    m_recojetenergy = recoJet->get_e();

    if (Verbosity() > 1)
    {
      cout << "matching by distance jet" << endl;
    }

    /// Set the matched truth jet characteristics to 0
    m_truthjetid = 0;
    m_truthjetp = 0;
    m_truthjetphi = 0;
    m_truthjeteta = 0;
    m_truthjetpt = 0;
    m_truthjetenergy = 0;
    m_truthjetpx = 0;
    m_truthjetpy = 0;
    m_truthjetpz = 0;

    Jet *truthjet = recoeval->max_truth_jet_by_energy(recoJet);
    if (truthjet)
    {
      m_truthjetid = truthjet->get_id();
      m_truthjetp = truthjet->get_p();
      m_truthjetpx = truthjet->get_px();
      m_truthjetpy = truthjet->get_py();
      m_truthjetpz = truthjet->get_pz();
      m_truthjeteta = truthjet->get_eta();
      m_truthjetphi = truthjet->get_phi();
      m_truthjetenergy = truthjet->get_e();
      m_truthjetpt = sqrt(m_truthjetpx * m_truthjetpx + m_truthjetpy * m_truthjetpy);
    }

    /// Check to make sure the truth jet node is available
    else if (truth_jets)
    {
      /// Match the reconstructed jet to the closest truth jet in delta R space
      /// Iterate over the truth jets
      float closestjet = 9999;
      for (JetMap::Iter truthIter = truth_jets->begin();
           truthIter != truth_jets->end();
           ++truthIter)
      {
        const Jet *truthJet = truthIter->second;

        float thisjetpt = truthJet->get_pt();
        if (thisjetpt < m_minjetpt)
          continue;

        float thisjeteta = truthJet->get_eta();
        float thisjetphi = truthJet->get_phi();

        float dphi = m_recojetphi - thisjetphi;
        if (dphi > TMath::Pi())
          dphi -= TMath::Pi() * 2.;
        if (dphi < -1. * TMath::Pi())
          dphi += TMath::Pi() * 2.;

        float deta = m_recojeteta - thisjeteta;
        /// Determine the distance in eta phi space between the reconstructed
        /// and truth jets
        m_dR = sqrt(pow(dphi, 2.) + pow(deta, 2.));

        /// If this truth jet is closer than the previous truth jet, it is
        /// closer and thus should be considered the truth jet
        if (m_dR < reco_jets->get_par() && m_dR < closestjet)
        {
          m_truthjetid = -9999;
          m_truthjetp = truthJet->get_p();
          m_truthjetphi = truthJet->get_phi();
          m_truthjeteta = truthJet->get_eta();
          m_truthjetpt = truthJet->get_pt();
          m_truthjetenergy = truthJet->get_e();
          m_truthjetpx = truthJet->get_px();
          m_truthjetpy = truthJet->get_py();
          m_truthjetpz = truthJet->get_pz();
          closestjet = m_dR;
        }
      }
    }
    m_recojettree->Fill();
  }
}

/**
 * This method gets clusters from the EMCal and stores them in a tree. It
 * also demonstrates how to get trigger emulator information. Clusters from
 * other containers can be obtained in a similar way (e.g. clusters from
 * the HCal, etc.)
 */
void eIDMLInterface::getEMCalClusters(PHCompositeNode *topNode)
{
  /// Get the raw cluster container
  /// Note: other cluster containers exist as well. Check out the node tree when
  /// you run a simulation, for example look for the node CLUSTER_EEMC
  RawClusterContainer *clusters = findNode::getClass<RawClusterContainer>(topNode, "CLUSTER_BECAL");

  if (!clusters)
  {
    cout << PHWHERE
         << "EMCal cluster node is missing, can't collect EMCal clusters"
         << endl;
    return;
  }

  /// Get the global vertex to determine the appropriate pseudorapidity of the clusters
  GlobalVertexMap *vertexmap = findNode::getClass<GlobalVertexMap>(topNode, "GlobalVertexMap");
  if (!vertexmap)
  {
    cout << "eIDMLInterface::getEmcalClusters - Fatal Error - GlobalVertexMap node is missing. Please turn on the do_global flag in the main macro in order to reconstruct the global vertex." << endl;
    assert(vertexmap);  // force quit

    return;
  }

  if (vertexmap->empty())
  {
    cout << "eIDMLInterface::getEmcalClusters - Fatal Error - GlobalVertexMap node is empty. Please turn on the do_global flag in the main macro in order to reconstruct the global vertex." << endl;
    return;
  }

  GlobalVertex *vtx = vertexmap->begin()->second;
  if (vtx == nullptr)
    return;

  /// Trigger emulator
  CaloTriggerInfo *trigger = findNode::getClass<CaloTriggerInfo>(topNode, "CaloTriggerInfo");

  /// Can obtain some trigger information if desired
  if (trigger)
  {
    m_E_4x4 = trigger->get_best_EMCal_4x4_E();
  }
  RawClusterContainer::ConstRange begin_end = clusters->getClusters();
  RawClusterContainer::ConstIterator clusIter;

  /// Loop over the EMCal clusters
  for (clusIter = begin_end.first;
       clusIter != begin_end.second;
       ++clusIter)
  {
    /// Get this cluster
    const RawCluster *cluster = clusIter->second;

    /// Get cluster characteristics
    /// This helper class determines the photon characteristics
    /// depending on the vertex position
    /// This is important for e.g. eta determination and E_T determination
    CLHEP::Hep3Vector vertex(vtx->get_x(), vtx->get_y(), vtx->get_z());
    CLHEP::Hep3Vector E_vec_cluster = RawClusterUtility::GetECoreVec(*cluster, vertex);
    m_clusenergy = E_vec_cluster.mag();
    m_cluseta = E_vec_cluster.pseudoRapidity();
    m_clustheta = E_vec_cluster.getTheta();
    m_cluspt = E_vec_cluster.perp();
    m_clusphi = E_vec_cluster.getPhi();

    m_phi_h->Fill(m_clusphi);
    m_eta_phi_h->Fill(m_clusphi, m_cluseta);

    if (m_cluspt < m_mincluspt)
      continue;

    m_cluspx = m_cluspt * cos(m_clusphi);
    m_cluspy = m_cluspt * sin(m_clusphi);
    m_cluspz = sqrt(m_clusenergy * m_clusenergy - m_cluspx * m_cluspx - m_cluspy * m_cluspy);

    //fill the cluster tree with all emcal clusters
    m_clustertree->Fill();
  }
}

/**
 * This function puts all of the tree branch assignments in one place so as to not
 * clutter up the eIDMLInterface::Init function.
 */
void eIDMLInterface::initializeTrees()
{
  m_outfile = new TFile();
  m_phi_h = new TH1F();
  m_eta_phi_h = new TH2F();
  //  m_recojettree = new TTree("jettree", "A tree with reconstructed jets");
  //  m_recojettree->Branch("m_recojetpt", &m_recojetpt, "m_recojetpt/D");
  //  m_recojettree->Branch("m_recojetid", &m_recojetid, "m_recojetid/I");
  //  m_recojettree->Branch("m_recojetpx", &m_recojetpx, "m_recojetpx/D");
  //  m_recojettree->Branch("m_recojetpy", &m_recojetpy, "m_recojetpy/D");
  //  m_recojettree->Branch("m_recojetpz", &m_recojetpz, "m_recojetpz/D");
  //  m_recojettree->Branch("m_recojetphi", &m_recojetphi, "m_recojetphi/D");
  //  m_recojettree->Branch("m_recojeteta", &m_recojeteta, "m_recojeteta/D");
  //  m_recojettree->Branch("m_recojetenergy", &m_recojetenergy, "m_recojetenergy/D");
  //  m_recojettree->Branch("m_truthjetid", &m_truthjetid, "m_truthjetid/I");
  //  m_recojettree->Branch("m_truthjetp", &m_truthjetp, "m_truthjetp/D");
  //  m_recojettree->Branch("m_truthjetphi", &m_truthjetphi, "m_truthjetphi/D");
  //  m_recojettree->Branch("m_truthjeteta", &m_truthjeteta, "m_truthjeteta/D");
  //  m_recojettree->Branch("m_truthjetpt", &m_truthjetpt, "m_truthjetpt/D");
  //  m_recojettree->Branch("m_truthjetenergy", &m_truthjetenergy, "m_truthjetenergy/D");
  //  m_recojettree->Branch("m_truthjetpx", &m_truthjetpx, "m_truthjetpx/D");
  //  m_recojettree->Branch("m_truthjetpy", &m_truthjetpy, "m_truthjyetpy/D");
  //  m_recojettree->Branch("m_truthjetpz", &m_truthjetpz, "m_truthjetpz/D");
  //  m_recojettree->Branch("m_dR", &m_dR, "m_dR/D");
  //
  //  m_truthjettree = new TTree("truthjettree", "A tree with truth jets");
  //  m_truthjettree->Branch("m_truthjetid", &m_truthjetid, "m_truthjetid/I");
  //  m_truthjettree->Branch("m_truthjetp", &m_truthjetp, "m_truthjetp/D");
  //  m_truthjettree->Branch("m_truthjetphi", &m_truthjetphi, "m_truthjetphi/D");
  //  m_truthjettree->Branch("m_truthjeteta", &m_truthjeteta, "m_truthjeteta/D");
  //  m_truthjettree->Branch("m_truthjetpt", &m_truthjetpt, "m_truthjetpt/D");
  //  m_truthjettree->Branch("m_truthjetenergy", &m_truthjetenergy, "m_truthjetenergy/D");
  //  m_truthjettree->Branch("m_truthjetpx", &m_truthjetpx, "m_truthjetpx/D");
  //  m_truthjettree->Branch("m_truthjetpy", &m_truthjetpy, "m_truthjetpy/D");
  //  m_truthjettree->Branch("m_truthjetpz", &m_truthjetpz, "m_truthjetpz/D");
  //  m_truthjettree->Branch("m_dR", &m_dR, "m_dR/D");
  //  m_truthjettree->Branch("m_recojetpt", &m_recojetpt, "m_recojetpt/D");
  //  m_truthjettree->Branch("m_recojetid", &m_recojetid, "m_recojetid/I");
  //  m_truthjettree->Branch("m_recojetpx", &m_recojetpx, "m_recojetpx/D");
  //  m_truthjettree->Branch("m_recojetpy", &m_recojetpy, "m_recojetpy/D");
  //  m_truthjettree->Branch("m_recojetpz", &m_recojetpz, "m_recojetpz/D");
  //  m_truthjettree->Branch("m_recojetphi", &m_recojetphi, "m_recojetphi/D");
  //  m_truthjettree->Branch("m_recojeteta", &m_recojeteta, "m_recojeteta/D");
  //  m_truthjettree->Branch("m_recojetenergy", &m_recojetenergy, "m_recojetenergy/D");
  //
  //  m_tracktree = new TTree("tracktree", "A tree with svtx tracks");
  //
  //  m_hepmctree = new TTree("hepmctree", "A tree with hepmc truth particles");
  //  m_hepmctree->Branch("m_partid1", &m_partid1, "m_partid1/I");
  //  m_hepmctree->Branch("m_partid2", &m_partid2, "m_partid2/I");
  //  m_hepmctree->Branch("m_x1", &m_x1, "m_x1/D");
  //  m_hepmctree->Branch("m_x2", &m_x2, "m_x2/D");
  //  m_hepmctree->Branch("m_mpi", &m_mpi, "m_mpi/I");
  //  m_hepmctree->Branch("m_process_id", &m_process_id, "m_process_id/I");
  //  m_hepmctree->Branch("m_truthenergy", &m_truthenergy, "m_truthenergy/D");
  //  m_hepmctree->Branch("m_trutheta", &m_trutheta, "m_trutheta/D");
  //  m_hepmctree->Branch("m_truthphi", &m_truthphi, "m_truthphi/D");
  //  m_hepmctree->Branch("m_truthpx", &m_truthpx, "m_truthpx/D");
  //  m_hepmctree->Branch("m_truthpy", &m_truthpy, "m_truthpy/D");
  //  m_hepmctree->Branch("m_truthpz", &m_truthpz, "m_truthpz/D");
  //  m_hepmctree->Branch("m_truthpt", &m_truthpt, "m_truthpt/D");
  //  m_hepmctree->Branch("m_numparticlesinevent", &m_numparticlesinevent, "m_numparticlesinevent/I");
  //  m_hepmctree->Branch("m_truthpid", &m_truthpid, "m_truthpid/I");

  m_truthtree = new TTree("T", "A tree one enetry for a truth g4 particles matched with reco objects");
  m_truthtree->Branch("m_truthenergy", &m_truthenergy, "m_truthenergy/D");
  m_truthtree->Branch("m_truthp", &m_truthp, "m_truthp/D");
  m_truthtree->Branch("m_truthpx", &m_truthpx, "m_truthpx/D");
  m_truthtree->Branch("m_truthpy", &m_truthpy, "m_truthpy/D");
  m_truthtree->Branch("m_truthpz", &m_truthpz, "m_truthpz/D");
  m_truthtree->Branch("m_truthpt", &m_truthpt, "m_truthpt/D");
  m_truthtree->Branch("m_truthphi", &m_truthphi, "m_truthphi/D");
  m_truthtree->Branch("m_trutheta", &m_trutheta, "m_trutheta/D");
  m_truthtree->Branch("m_truthpid", &m_truthpid, "m_truthpid/I");

  m_truthtree->Branch("m_tr_px", &m_tr_px, "m_tr_px/D");
  m_truthtree->Branch("m_tr_py", &m_tr_py, "m_tr_py/D");
  m_truthtree->Branch("m_tr_pz", &m_tr_pz, "m_tr_pz/D");
  m_truthtree->Branch("m_tr_p", &m_tr_p, "m_tr_p/D");
  m_truthtree->Branch("m_tr_pt", &m_tr_pt, "m_tr_pt/D");
  m_truthtree->Branch("m_tr_phi", &m_tr_phi, "m_tr_phi/D");
  m_truthtree->Branch("m_tr_eta", &m_tr_eta, "m_tr_eta/D");
  m_truthtree->Branch("m_charge", &m_charge, "m_charge/I");
  m_truthtree->Branch("m_chisq", &m_chisq, "m_chisq/D");
  m_truthtree->Branch("m_ndf", &m_ndf, "m_ndf/I");
  m_truthtree->Branch("m_dca", &m_dca, "m_dca/D");
  m_truthtree->Branch("m_tr_x", &m_tr_x, "m_tr_x/D");
  m_truthtree->Branch("m_tr_y", &m_tr_y, "m_tr_y/D");
  m_truthtree->Branch("m_tr_z", &m_tr_z, "m_tr_z/D");

  const string xyzt[] = {"x", "y", "z", "t"};

  for (std::string _calo_name : _calo_names)
  {
    for (int i = 0; i < 4; i++)
    {
      string bname = _calo_name + "_proj_" + xyzt[i];
      string bdef = bname + "/F";

      // fourth element is the path length
      if (i == 3)
      {
        bdef = _calo_name + "_proj_path_length" + "/F";
      }

      m_truthtree->Branch(bname.c_str(), &m_CaloDataMap[_calo_name].m_TTree_proj_vec[i], bdef.c_str());
    }

    for (int i = 0; i < 3; i++)
    {
      string bname = _calo_name + "_proj_p" + xyzt[i];
      string bdef = bname + "/F";
      m_truthtree->Branch(bname.c_str(), &m_CaloDataMap[_calo_name].m_TTree_proj_p_vec[i], bdef.c_str());
    }

    //  static const int nTowerInPatch = m_sizeTowerPatch * m_sizeTowerPatch;
    m_truthtree->Branch((_calo_name + "_Tower_E3x3").c_str(), &m_CaloDataMap[_calo_name].m_E3x3, (_calo_name + "_Tower_E3x3/F").c_str());
    m_truthtree->Branch((_calo_name + "_Tower_E5x5").c_str(), &m_CaloDataMap[_calo_name].m_E5x5, (_calo_name + "_Tower_E5x5/F").c_str());
    m_truthtree->Branch((_calo_name + "_Tower_E7x7").c_str(), &m_CaloDataMap[_calo_name].m_E7x7, (_calo_name + "_Tower_E7x7/F").c_str());
    m_truthtree->Branch((_calo_name + "_centralTowerBinEta").c_str(), &m_CaloDataMap[_calo_name].m_centralTowerBinEta, (_calo_name + "_centralTowerBinEta/I").c_str());
    m_truthtree->Branch((_calo_name + "_centralTowerBinPhi").c_str(), &m_CaloDataMap[_calo_name].m_centralTowerBinPhi, (_calo_name + "_centralTowerBinPhi/I").c_str());
    m_truthtree->Branch((_calo_name + "_nTowerInPatch").c_str(), &m_CaloDataMap[_calo_name].nTowerInPatch, (_calo_name + "_nTowerInPatch/I").c_str());
    m_truthtree->Branch((_calo_name + "_Tower_dEta").c_str(), m_CaloDataMap[_calo_name].m_TTree_Tower_dEta.data(), (_calo_name + "_Tower_dEta[" + _calo_name + "_nTowerInPatch]/F").c_str());
    m_truthtree->Branch((_calo_name + "_Tower_dPhi").c_str(), m_CaloDataMap[_calo_name].m_TTree_Tower_dPhi.data(), (_calo_name + "_Tower_dPhi[" + _calo_name + "_nTowerInPatch]/F").c_str());
    m_truthtree->Branch((_calo_name + "_Tower_iEta_patch").c_str(), m_CaloDataMap[_calo_name].m_TTree_Tower_iEta_patch.data(), (_calo_name + "_Tower_iEta_patch[" + _calo_name + "_nTowerInPatch]/I").c_str());
    m_truthtree->Branch((_calo_name + "_Tower_iPhi_patch").c_str(), m_CaloDataMap[_calo_name].m_TTree_Tower_iPhi_patch.data(), (_calo_name + "_Tower_iPhi_patch[" + _calo_name + "_nTowerInPatch]/I").c_str());
    m_truthtree->Branch((_calo_name + "_Tower_E").c_str(), m_CaloDataMap[_calo_name].m_TTree_Tower_E.data(), (_calo_name + "_Tower_E[" + _calo_name + "_nTowerInPatch]/F").c_str());

    //  m_clustertree = new TTree("clustertree", "A tree with emcal clusters");
    //  m_clustertree->Branch("m_clusenergy", &m_clusenergy, "m_clusenergy/D");
    //  m_clustertree->Branch("m_cluseta", &m_cluseta, "m_cluseta/D");
    //  m_clustertree->Branch("m_clustheta", &m_clustheta, "m_clustheta/D");
    //  m_clustertree->Branch("m_cluspt", &m_cluspt, "m_cluspt/D");
    //  m_clustertree->Branch("m_clusphi", &m_clusphi, "m_clusphi/D");
    //  m_clustertree->Branch("m_cluspx", &m_cluspx, "m_cluspx/D");
    //  m_clustertree->Branch("m_cluspy", &m_cluspy, "m_cluspy/D");
    //  m_clustertree->Branch("m_cluspz", &m_cluspz, "m_cluspz/D");
    //  m_clustertree->Branch("m_E_4x4", &m_E_4x4, "m_E_4x4/D");
  }  //  for (std::string _calo_name : _calo_names)
}

void eIDMLInterface::CaloData::initializeVariables()
{
  std::fill(m_TTree_proj_vec.begin(), m_TTree_proj_vec.end(), -9999);
  std::fill(m_TTree_proj_p_vec.begin(), m_TTree_proj_p_vec.end(), -9999);

  std::fill(m_TTree_Tower_dEta.begin(), m_TTree_Tower_dEta.end(), 0);
  std::fill(m_TTree_Tower_dPhi.begin(), m_TTree_Tower_dPhi.end(), 0);
  std::fill(m_TTree_Tower_iEta_patch.begin(), m_TTree_Tower_iEta_patch.end(), 0);
  std::fill(m_TTree_Tower_iPhi_patch.begin(), m_TTree_Tower_iPhi_patch.end(), 0);
  std::fill(m_TTree_Tower_E.begin(), m_TTree_Tower_E.end(), 0);

  m_centralTowerBinEta = -9999;
  m_centralTowerBinPhi = -9999;
  m_E3x3 = 0;
  m_E5x5 = 0;
  m_E7x7 = 0;
}
/**
 * This function initializes all of the member variables in this class so that there
 * are no variables that might not be set before e.g. writing them to the output
 * trees. 
 */
void eIDMLInterface::initializeVariables()
{
  m_partid1 = -99;
  m_partid2 = -99;
  m_x1 = -99;
  m_x2 = -99;
  m_mpi = -99;
  m_process_id = -99;
  m_truthenergy = -99;
  m_trutheta = -99;
  m_truthphi = -99;
  m_truthp = -99;
  m_truthpx = -99;
  m_truthpy = -99;
  m_truthpz = -99;
  m_truthpt = -99;
  m_numparticlesinevent = -99;
  m_truthpid = -99;

  m_tr_px = -99;
  m_tr_py = -99;
  m_tr_pz = -99;
  m_tr_p = -99;
  m_tr_pt = -99;
  m_tr_phi = -99;
  m_tr_eta = -99;
  m_charge = -99;
  m_chisq = -99;
  m_ndf = -99;
  m_dca = -99;
  m_tr_x = -99;
  m_tr_y = -99;
  m_tr_z = -99;
  m_tr_pion_loglikelihood = -99;
  m_tr_kaon_loglikelihood = -99;
  m_tr_proton_loglikelihood = -99;

  m_truthtrackpy = -99;
  m_truthtrackpz = -99;
  m_truthtrackp = -99;
  m_truthtracke = -99;
  m_truthtrackpt = -99;
  m_truthtrackphi = -99;
  m_truthtracketa = -99;
  m_truthtrackpid = -99;

  m_recojetpt = -99;
  m_recojetid = -99;
  m_recojetpx = -99;
  m_recojetpy = -99;
  m_recojetpz = -99;
  m_recojetphi = -99;
  m_recojetp = -99;
  m_recojetenergy = -99;
  m_recojeteta = -99;
  m_truthjetid = -99;
  m_truthjetp = -99;
  m_truthjetphi = -99;
  m_truthjeteta = -99;
  m_truthjetpt = -99;
  m_truthjetenergy = -99;
  m_truthjetpx = -99;
  m_truthjetpy = -99;
  m_truthjetpz = -99;
  m_dR = -99;

  for (std::string _calo_name : _calo_names)
  {
    m_CaloDataMap[_calo_name].initializeVariables();
  }
}
