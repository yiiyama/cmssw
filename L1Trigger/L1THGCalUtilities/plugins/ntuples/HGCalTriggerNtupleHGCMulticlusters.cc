#include "DataFormats/L1THGCal/interface/HGCalMulticluster.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "L1Trigger/L1THGCal/interface/HGCalTriggerGeometryBase.h"
#include "L1Trigger/L1THGCalUtilities/interface/HGCalTriggerNtupleBase.h"
#include "L1Trigger/L1THGCal/interface/backend/HGCalTriggerClusterIdentificationBase.h"

#include "SimDataFormats/CaloHit/interface/PCaloHit.h"
#include "SimDataFormats/Vertex/interface/SimVertex.h"
#include "SimDataFormats/Track/interface/SimTrack.h"

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
    edm::EDGetToken simTracksToken_;
    edm::EDGetToken simVerticesToken_;
    edm::EDGetToken simHitsTokenEE_;
    edm::EDGetToken simHitsTokenHEfront_;
    edm::EDGetToken simHitsTokenHEback_;

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

    std::vector<int> cl3d_simtrack_pid_;
    std::vector<float> cl3d_simtrack_pt_;
    std::vector<float> cl3d_simtrack_eta_;
    std::vector<float> cl3d_simtrack_phi_;
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
    simTracksToken_ = collector.consumes<std::vector<SimTrack>>(conf.getParameter<edm::InputTag>("SimTracks"));
    simVerticesToken_ = collector.consumes<std::vector<SimVertex>>(conf.getParameter<edm::InputTag>("SimVertices"));
    simHitsTokenEE_ = collector.consumes<std::vector<PCaloHit>>(conf.getParameter<edm::InputTag>("SimHitsEE"));
    simHitsTokenHEfront_ = collector.consumes<std::vector<PCaloHit>>(conf.getParameter<edm::InputTag>("SimHitsHEfront"));
    simHitsTokenHEback_ = collector.consumes<std::vector<PCaloHit>>(conf.getParameter<edm::InputTag>("SimHitsHEback"));

    tree.Branch(withPrefix("simtrack_pid"), &cl3d_simtrack_pid_);
    tree.Branch(withPrefix("simtrack_pt"), &cl3d_simtrack_pt_);
    tree.Branch(withPrefix("simtrack_eta"), &cl3d_simtrack_eta_);
    tree.Branch(withPrefix("simtrack_phi"), &cl3d_simtrack_phi_);
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
  }

  if (matchSimTrack_) {
    auto haveOverlap([](std::set<uint32_t> set1, std::set<uint32_t> set2)->bool {
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

    cl3d_simtrack_pid_.resize(cl3d_n_);
    cl3d_simtrack_pt_.resize(cl3d_n_);
    cl3d_simtrack_eta_.resize(cl3d_n_);
    cl3d_simtrack_phi_.resize(cl3d_n_);

    // loop over the hits, fill simTrackToHits map
    // loop over 3D clusters, collect matching tracks
    // trace the track to parent

    // retrieve geometry
    edm::ESHandle<HGCalTriggerGeometryBase> geometryHandle;
    es.get<CaloGeometryRecord>().get(geometryHandle);
    auto& geometry(*geometryHandle);

    edm::Handle<std::vector<SimTrack>> simTracksHandle;
    e.getByToken(simTracksToken_, simTracksHandle);
    auto const& simTracks(*simTracksHandle);

    std::unordered_map<int, SimTrack const*> tracks;
    for (auto& track : simTracks)
      tracks.emplace(track.trackId(), &track);

    edm::Handle<std::vector<SimVertex>> simVerticesHandle;
    e.getByToken(simVerticesToken_, simVerticesHandle);
    auto const& simVertices(*simVerticesHandle);

    std::unordered_map<int, SimVertex const*> vertices;
    for (auto& vertex : simVertices)
      vertices.emplace(vertex.vertexId(), &vertex);

    std::unordered_map<int, std::set<uint32_t>> simTrackToCells;

    for (auto* token : {&simHitsTokenEE_, &simHitsTokenHEfront_, &simHitsTokenHEback_}) {
      edm::Handle<std::vector<PCaloHit>> handle;
      e.getByToken(*token, handle);
      auto const& simhits(*handle);

      for (auto const& simhit : simhits) {
        if (simhit.geantTrackId() <= 0)
          continue;

        uint32_t tcId;
        try {
          tcId = geometry.getTriggerCellFromCell(simhit.id());
        }
        catch (cms::Exception const& ex) {
          edm::LogError("CaloTruthCellsProducer") << ex.what();
          continue;
        }

        simTrackToCells[simhit.geantTrackId()].insert(tcId);
      }
    }

    unsigned icl(0);
    for(auto cl3d_itr=multiclusters.begin(0); cl3d_itr!=multiclusters.end(0); ++cl3d_itr, ++icl) {
      std::set<uint32_t> cellIds;

      for (auto&& clItr(cl3d_itr->constituents_begin()); clItr != cl3d_itr->constituents_end(); ++clItr) {
        auto const& cluster(*clItr->second);
        for (auto&& cellItr(cluster.constituents_begin()); cellItr != cluster.constituents_end(); ++cellItr)
          cellIds.insert(cellItr->first);
      }

      std::vector<int> trackIds;

      for (auto&& th : simTrackToCells) {
        if (haveOverlap(cellIds, th.second))
          trackIds.push_back(th.first);
      }

      if (trackIds.size() == 0) {
        cl3d_simtrack_pid_[icl] = 0;
        cl3d_simtrack_pt_[icl] = 0.;
        cl3d_simtrack_eta_[icl] = 0.;
        cl3d_simtrack_phi_[icl] = 0.;
      }
      else if (trackIds.size() == 1) {
        auto& track(*tracks[trackIds[0]]);
        cl3d_simtrack_pid_[icl] = track.type();
        auto& momentum(track.momentum());
        cl3d_simtrack_pt_[icl] = momentum.pt();
        cl3d_simtrack_eta_[icl] = momentum.eta();
        cl3d_simtrack_phi_[icl] = momentum.phi();
      }
      else {
        std::vector<std::vector<int>> histories(trackIds.size());

        for (unsigned iT(0); iT != trackIds.size(); ++iT) {
          auto tItr(tracks.find(trackIds[iT]));
          if (tItr == tracks.end())
            continue;
          
          auto const* track(tItr->second);
          while (!track->noVertex()) {
            auto& vertex(*vertices[track->vertIndex()]);
            if (vertex.noParent())
              break;

            histories[iT].push_back(vertex.parentIndex());
            tItr = tracks.find(vertex.parentIndex());
            if (tItr == tracks.end())
              break;
            
            track = tItr->second;
          }
        }

        unsigned iHist(1);
        int ancestor(-1);
        while (true) {
          if (histories[0].size() < iHist)
            break;

          int candidate(histories[0][histories[0].size() - iHist]);

          unsigned iT(1);
          for (; iT != histories.size(); ++iT) {
            if (histories[iT].size() < iHist)
              break;

            if (histories[iT][histories[iT].size() - iHist] != candidate)
              break;
          }
          if (iT == histories.size()) {
            ancestor = candidate;
            ++iHist;
          }
          else
            break;
        }

        if (ancestor == -1) {
          cl3d_simtrack_pid_[icl] = 0;
          cl3d_simtrack_pt_[icl] = 0.;
          cl3d_simtrack_eta_[icl] = 0.;
          cl3d_simtrack_phi_[icl] = 0.;
        }
        else {
          auto& track(*tracks[ancestor]);
          cl3d_simtrack_pid_[icl] = track.type();
          auto& momentum(track.momentum());
          cl3d_simtrack_pt_[icl] = momentum.pt();
          cl3d_simtrack_eta_[icl] = momentum.eta();
          cl3d_simtrack_phi_[icl] = momentum.phi();
        }
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

  cl3d_simtrack_pid_.clear();
  cl3d_simtrack_pt_.clear();
  cl3d_simtrack_eta_.clear();
  cl3d_simtrack_phi_.clear();
}
