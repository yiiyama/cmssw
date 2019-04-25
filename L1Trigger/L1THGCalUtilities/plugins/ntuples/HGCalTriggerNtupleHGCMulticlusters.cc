#include "DataFormats/L1THGCal/interface/HGCalMulticluster.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "L1Trigger/L1THGCal/interface/HGCalTriggerGeometryBase.h"
#include "L1Trigger/L1THGCalUtilities/interface/HGCalTriggerNtupleBase.h"
#include "L1Trigger/L1THGCal/interface/backend/HGCalTriggerClusterIdentificationBase.h"

#include "SimDataFormats/CaloAnalysis/interface/CaloParticleFwd.h"
#include "SimDataFormats/CaloAnalysis/interface/CaloParticle.h"
#include "SimDataFormats/CaloAnalysis/interface/SimCluster.h"
#include "DataFormats/Math/interface/LorentzVector.h"

class HGCalTriggerNtupleHGCMulticlusters : public HGCalTriggerNtupleBase
{

  public:
    HGCalTriggerNtupleHGCMulticlusters(const edm::ParameterSet& conf);
    ~HGCalTriggerNtupleHGCMulticlusters() override {}
    void initialize(TTree&, const edm::ParameterSet&, edm::ConsumesCollector&&) final;
    void fill(const edm::Event& e, const edm::EventSetup& es) final;

  private:
    void clear() final;

    edm::EDGetToken multiclusters_token_;
    edm::EDGetToken caloParticlesToken_;

    std::unique_ptr<HGCalTriggerClusterIdentificationBase> id_;

    bool matchSimTrack_;

    int cl3d_n_ ;
    std::vector<uint32_t> cl3d_id_;
    std::vector<float> cl3d_pt_;
    std::vector<float> cl3d_energy_;
    std::vector<float> cl3d_eta_;
    std::vector<float> cl3d_phi_;
    std::vector<int> cl3d_clusters_n_;
    std::vector<std::vector<uint32_t>> cl3d_clusters_id_;
    // cluster shower shapes
    std::vector<int> cl3d_showerlength_;
    std::vector<int> cl3d_coreshowerlength_;
    std::vector<int> cl3d_firstlayer_;
    std::vector<int> cl3d_maxlayer_;
    std::vector<float> cl3d_seetot_;
    std::vector<float> cl3d_seemax_;
    std::vector<float> cl3d_spptot_;
    std::vector<float> cl3d_sppmax_;
    std::vector<float> cl3d_szz_;
    std::vector<float> cl3d_srrtot_;
    std::vector<float> cl3d_srrmax_;
    std::vector<float> cl3d_srrmean_;
    std::vector<float> cl3d_emaxe_;
    std::vector<float> cl3d_bdteg_;
    std::vector<int> cl3d_quality_;

    std::vector<std::vector<int>> cl3d_simtracks_pid_;
    std::vector<std::vector<int>> cl3d_simtracks_eventId_;
    std::vector<std::vector<float>> cl3d_simtracks_pt_;
    std::vector<std::vector<float>> cl3d_simtracks_eta_;
    std::vector<std::vector<float>> cl3d_simtracks_phi_;
};

DEFINE_EDM_PLUGIN(HGCalTriggerNtupleFactory,
    HGCalTriggerNtupleHGCMulticlusters,
    "HGCalTriggerNtupleHGCMulticlusters" );

HGCalTriggerNtupleHGCMulticlusters::
HGCalTriggerNtupleHGCMulticlusters(const edm::ParameterSet& conf):HGCalTriggerNtupleBase(conf)
{
}

void
HGCalTriggerNtupleHGCMulticlusters::
initialize(TTree& tree, const edm::ParameterSet& conf, edm::ConsumesCollector&& collector)
{
  multiclusters_token_ = collector.consumes<l1t::HGCalMulticlusterBxCollection>(conf.getParameter<edm::InputTag>("Multiclusters"));
  id_ = std::unique_ptr<HGCalTriggerClusterIdentificationBase>{ HGCalTriggerClusterIdentificationFactory::get()->create("HGCalTriggerClusterIdentificationBDT") };
  id_->initialize(conf.getParameter<edm::ParameterSet>("EGIdentification")); 

  std::string prefix(conf.getUntrackedParameter<std::string>("Prefix", "cl3d"));

  std::string bname;
  auto withPrefix([&prefix, &bname](char const* vname)->char const* {
      bname = prefix + "_" + vname;
      return bname.c_str();
    });

  tree.Branch(withPrefix("n"), &cl3d_n_, (prefix + "_n/I").c_str());
  tree.Branch(withPrefix("id"), &cl3d_id_);
  tree.Branch(withPrefix("pt"), &cl3d_pt_);
  tree.Branch(withPrefix("energy"), &cl3d_energy_);
  tree.Branch(withPrefix("eta"), &cl3d_eta_);
  tree.Branch(withPrefix("phi"), &cl3d_phi_);
  tree.Branch(withPrefix("clusters_n"), &cl3d_clusters_n_);
  tree.Branch(withPrefix("clusters_id"), &cl3d_clusters_id_);
  tree.Branch(withPrefix("showerlength"), &cl3d_showerlength_);
  tree.Branch(withPrefix("coreshowerlength"), &cl3d_coreshowerlength_);
  tree.Branch(withPrefix("firstlayer"), &cl3d_firstlayer_);
  tree.Branch(withPrefix("maxlayer"), &cl3d_maxlayer_);
  tree.Branch(withPrefix("seetot"), &cl3d_seetot_);
  tree.Branch(withPrefix("seemax"), &cl3d_seemax_);
  tree.Branch(withPrefix("spptot"), &cl3d_spptot_);
  tree.Branch(withPrefix("sppmax"), &cl3d_sppmax_);
  tree.Branch(withPrefix("szz"), &cl3d_szz_);
  tree.Branch(withPrefix("srrtot"), &cl3d_srrtot_);
  tree.Branch(withPrefix("srrmax"), &cl3d_srrmax_);
  tree.Branch(withPrefix("srrmean"), &cl3d_srrmean_);
  tree.Branch(withPrefix("emaxe"), &cl3d_emaxe_);
  tree.Branch(withPrefix("bdteg"), &cl3d_bdteg_);
  tree.Branch(withPrefix("quality"), &cl3d_quality_);

  matchSimTrack_ = conf.getUntrackedParameter<bool>("MatchSimTrack", false);
  if (matchSimTrack_) {
    caloParticlesToken_ = collector.consumes<CaloParticleCollection>(conf.getParameter<edm::InputTag>("CaloParticles"));

    std::cout << "matching calo" << std::endl;

    tree.Branch(withPrefix("simtracks_pid"), &cl3d_simtracks_pid_);
    tree.Branch(withPrefix("simtracks_eventId"), &cl3d_simtracks_eventId_);
    tree.Branch(withPrefix("simtracks_pt"), &cl3d_simtracks_pt_);
    tree.Branch(withPrefix("simtracks_eta"), &cl3d_simtracks_eta_);
    tree.Branch(withPrefix("simtracks_phi"), &cl3d_simtracks_phi_);
  }
}

void
HGCalTriggerNtupleHGCMulticlusters::
fill(const edm::Event& e, const edm::EventSetup& es)
{
  // retrieve clusters 3D
  edm::Handle<l1t::HGCalMulticlusterBxCollection> multiclusters_h;
  e.getByToken(multiclusters_token_, multiclusters_h);
  const l1t::HGCalMulticlusterBxCollection& multiclusters = *multiclusters_h;

  std::vector<std::set<uint32_t>> constituentIds(multiclusters.size(0));

  clear();
  for(auto cl3d_itr=multiclusters.begin(0); cl3d_itr!=multiclusters.end(0); cl3d_itr++)
  {
    cl3d_n_++;
    cl3d_id_.emplace_back(cl3d_itr->detId());
    // physical values 
    cl3d_pt_.emplace_back(cl3d_itr->pt());
    cl3d_energy_.emplace_back(cl3d_itr->energy());
    cl3d_eta_.emplace_back(cl3d_itr->eta());
    cl3d_phi_.emplace_back(cl3d_itr->phi());
    cl3d_clusters_n_.emplace_back(cl3d_itr->constituents().size());
    cl3d_showerlength_.emplace_back(cl3d_itr->showerLength());
    cl3d_coreshowerlength_.emplace_back(cl3d_itr->coreShowerLength());
    cl3d_firstlayer_.emplace_back(cl3d_itr->firstLayer());
    cl3d_maxlayer_.emplace_back(cl3d_itr->maxLayer());
    cl3d_seetot_.emplace_back(cl3d_itr->sigmaEtaEtaTot());
    cl3d_seemax_.emplace_back(cl3d_itr->sigmaEtaEtaMax());
    cl3d_spptot_.emplace_back(cl3d_itr->sigmaPhiPhiTot());
    cl3d_sppmax_.emplace_back(cl3d_itr->sigmaPhiPhiMax());
    cl3d_szz_.emplace_back(cl3d_itr->sigmaZZ());
    cl3d_srrtot_.emplace_back(cl3d_itr->sigmaRRTot());
    cl3d_srrmax_.emplace_back(cl3d_itr->sigmaRRMax());
    cl3d_srrmean_.emplace_back(cl3d_itr->sigmaRRMean());
    cl3d_emaxe_.emplace_back(cl3d_itr->eMax()/cl3d_itr->energy());
    cl3d_bdteg_.emplace_back(id_->value(*cl3d_itr));
    cl3d_quality_.emplace_back(cl3d_itr->hwQual());

    // Retrieve indices of trigger cells inside cluster
    cl3d_clusters_id_.emplace_back(cl3d_itr->constituents().size());
    std::transform(cl3d_itr->constituents_begin(), cl3d_itr->constituents_end(),
        cl3d_clusters_id_.back().begin(), [](const std::pair<uint32_t,edm::Ptr<l1t::HGCalCluster>>& id_cl){return id_cl.second->detId();}
        );

    if (matchSimTrack_) {
      for (auto&& id_cl_itr(cl3d_itr->constituents_begin()); id_cl_itr != cl3d_itr->constituents_end(); ++id_cl_itr) {
        auto& cluster(*id_cl_itr->second);
        for (auto&& id_cell_itr(cluster.constituents_begin()); id_cell_itr != cluster.constituents_end(); ++id_cell_itr)
          constituentIds[cl3d_n_ - 1].insert(id_cell_itr->second->detId());
      }
    }
  }

  if (matchSimTrack_) {
    cl3d_simtracks_pid_.resize(cl3d_n_);
    cl3d_simtracks_eventId_.resize(cl3d_n_);
    cl3d_simtracks_pt_.resize(cl3d_n_);
    cl3d_simtracks_eta_.resize(cl3d_n_);
    cl3d_simtracks_phi_.resize(cl3d_n_);

    edm::Handle<CaloParticleCollection> caloParticlesHandle;
    e.getByToken(caloParticlesToken_, caloParticlesHandle);
    CaloParticleCollection const& caloParticles(*caloParticlesHandle);

    // retrieve geometry
    edm::ESHandle<HGCalTriggerGeometryBase> geometryHandle;
    es.get<CaloGeometryRecord>().get(geometryHandle);
    auto const& geometry(*geometryHandle);

    // collect all trigger cells belonging to SimClusters in the CaloParticles
    std::vector<std::vector<std::set<uint32_t>>> clusterConstituents(caloParticles.size());

    auto hasOverlap([](std::set<uint32_t> set1, std::set<uint32_t> set2)->bool {
        auto&& itr1(set1.begin());
        auto&& itr2(set2.begin());
        auto&& end1(set1.end());
        auto&& end2(set2.end());

        while (itr1 != end1 && itr2 != end2) {
          if (*itr1 < *itr2)
            ++itr1;
          else if (*itr2 < *itr1)
            ++itr2;
          else
            return true;
        }

        return false;
      });

    unsigned iP(0);
    for (auto const& caloParticle : caloParticles) {
      clusterConstituents[iP].resize(caloParticle.simClusters().size());

      unsigned iC(0);
      for (auto const& simClusterRef : caloParticle.simClusters()) {
        auto const& simCluster(*simClusterRef);

        for (auto const& hAndF : simCluster.hits_and_fractions()) {
          uint32_t tcId;
          try {
            tcId = geometry.getTriggerCellFromCell(hAndF.first);
          }
          catch (cms::Exception const& ex) {
            edm::LogError("CaloTruthCellsProducer") << ex.what();
            continue;
          }

          clusterConstituents[iP][iC].insert(tcId);
        }

        ++iC;
      }

      ++iP;
    }

    for (int icl3d(0); icl3d != cl3d_n_; ++icl3d) {
      for (iP = 0; iP != caloParticles.size(); ++iP) {
        // Find SimClusters that overlap with the multicluster.
        // If only one overlaps, use the SimTrack of the given SimCluster
        // If all overlap, use the SimTrack of the CaloParticle
        // If multiple but not all overlap, flag it unknown unless it's e+e- (then call it a photon)

        auto const& caloParticle(caloParticles[iP]);
        auto const& simClusters(caloParticle.simClusters());

        std::vector<unsigned> overlapping;
        bool atLeastOne(false);
        bool all(true);

        for (unsigned iC(0); iC != simClusters.size(); ++iC) {
          if (hasOverlap(clusterConstituents[iP][iC], constituentIds[icl3d])) {
            overlapping.push_back(iC);
            atLeastOne = true;
          }
          else
            all = false;
        }

        if (!atLeastOne)
          continue;

        int eventId(caloParticle.g4Tracks().at(0).eventId().event());

        int pid(0);
        math::XYZTLorentzVector momentum(0., 0., 0., 0.);

        if (all) {
          pid = caloParticle.g4Tracks().at(0).type();
          momentum = caloParticle.g4Tracks().at(0).momentum();
        }
        else if (overlapping.size() == 1) {
          pid = simClusters.at(overlapping[0])->g4Tracks().at(0).type();
          momentum = simClusters.at(overlapping[0])->g4Tracks().at(0).momentum();
        }
        else {
          if (overlapping.size() == 2) {
            int pdgId1(simClusters.at(overlapping[0])->g4Tracks().at(0).type());
            int pdgId2(simClusters.at(overlapping[1])->g4Tracks().at(0).type());
            if (pdgId1 * pdgId2 == -121)
              pid = 22;
          }

          for (unsigned iC : overlapping)
            momentum += simClusters.at(iC)->g4Tracks().at(0).momentum();
        }

        cl3d_simtracks_pid_[icl3d].push_back(pid);
        cl3d_simtracks_eventId_[icl3d].push_back(eventId);
        cl3d_simtracks_pt_[icl3d].push_back(momentum.pt());
        cl3d_simtracks_eta_[icl3d].push_back(momentum.eta());
        cl3d_simtracks_phi_[icl3d].push_back(momentum.phi());
      }
    }
  }
}


void
HGCalTriggerNtupleHGCMulticlusters::
clear()
{
  cl3d_n_ = 0;
  cl3d_id_.clear();
  cl3d_pt_.clear();
  cl3d_energy_.clear();
  cl3d_eta_.clear();
  cl3d_phi_.clear();
  cl3d_clusters_n_.clear();
  cl3d_clusters_id_.clear();
  cl3d_showerlength_.clear();
  cl3d_coreshowerlength_.clear();
  cl3d_firstlayer_.clear();
  cl3d_maxlayer_.clear();
  cl3d_seetot_.clear();
  cl3d_seemax_.clear();
  cl3d_spptot_.clear();
  cl3d_sppmax_.clear();
  cl3d_szz_.clear();
  cl3d_srrtot_.clear();
  cl3d_srrmax_.clear();
  cl3d_srrmean_.clear();
  cl3d_emaxe_.clear();
  cl3d_bdteg_.clear();
  cl3d_quality_.clear();

  cl3d_simtracks_pid_.clear();
  cl3d_simtracks_eventId_.clear();
  cl3d_simtracks_pt_.clear();
  cl3d_simtracks_eta_.clear();
  cl3d_simtracks_phi_.clear();
}




