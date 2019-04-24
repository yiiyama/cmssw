#include "DataFormats/L1THGCal/interface/HGCalMulticluster.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "L1Trigger/L1THGCal/interface/HGCalTriggerGeometryBase.h"
#include "L1Trigger/L1THGCalUtilities/interface/HGCalTriggerNtupleBase.h"
#include "L1Trigger/L1THGCal/interface/backend/HGCalTriggerClusterIdentificationBase.h"
#include "RecoLocalCalo/HGCalRecAlgos/interface/RecHitTools.h"

#include "FastSimulation/Event/interface/FSimEvent.h"
#include "DataFormats/GeometrySurface/interface/Plane.h"
#include "TrackingTools/TrajectoryParametrization/interface/CurvilinearTrajectoryError.h"
#include "TrackPropagation/RungeKutta/interface/defaultRKPropagator.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"

// Taken from HGCalAnalysis
class SimpleTrackPropagator {
public:
  SimpleTrackPropagator(MagneticField const&, double hgcalSurface);
  bool propagate(math::XYZTLorentzVectorD const& momentum, math::XYZTLorentzVectorD const& position,
                 float charge, math::XYZVectorD& coords) const;

  double abszTarget() const { return abszTarget_; }

private:
  MagneticField const& field_;
  Plane::PlanePointer targetPlaneForward_;
  Plane::PlanePointer targetPlaneBackward_;
  float abszTarget_;
  CurvilinearTrajectoryError err_;
  defaultRKPropagator::Product prod_;
};

class HGCalTriggerNtupleHGCMulticlusters : public HGCalTriggerNtupleBase
{

  public:
    HGCalTriggerNtupleHGCMulticlusters(const edm::ParameterSet& conf);
    ~HGCalTriggerNtupleHGCMulticlusters() override { delete propagator_; }
    void initialize(TTree&, const edm::ParameterSet&, edm::ConsumesCollector&&) final;
    void beginRun(const edm::Run& r, const edm::EventSetup& es) override;
    void fill(const edm::Event& e, const edm::EventSetup& es) final;

  private:
    void clear() final;

    edm::EDGetToken multiclusters_token_;
    edm::EDGetToken simTracksToken_;
    edm::EDGetToken simVerticesToken_;

    std::unique_ptr<HGCalTriggerClusterIdentificationBase> id_;

    bool matchSimTrack_;
    FSimEvent* fsimEvent_;
    SimpleTrackPropagator* propagator_;

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
    std::vector<std::vector<float>> cl3d_simtracks_pt_;
    std::vector<std::vector<float>> cl3d_simtracks_eta_;
    std::vector<std::vector<float>> cl3d_simtracks_phi_;
};

DEFINE_EDM_PLUGIN(HGCalTriggerNtupleFactory,
    HGCalTriggerNtupleHGCMulticlusters,
    "HGCalTriggerNtupleHGCMulticlusters" );

SimpleTrackPropagator::SimpleTrackPropagator(MagneticField const& _field, double _hgcalSurface) :
  field_(_field),
  targetPlaneForward_(Plane::build(Plane::PositionType(0, 0, std::abs(_hgcalSurface)), Plane::RotationType())),
  targetPlaneBackward_(Plane::build(Plane::PositionType(0, 0, -std::abs(_hgcalSurface)), Plane::RotationType())),
  abszTarget_(std::abs(_hgcalSurface)),
  prod_(&_field, alongMomentum, 5.e-5)
{
  ROOT::Math::SMatrixIdentity id;
  AlgebraicSymMatrix55 C(id);
  C *= 0.001;
  err_ = CurvilinearTrajectoryError(C);
}

bool
SimpleTrackPropagator::propagate(math::XYZTLorentzVectorD const& _momentum,
                                 math::XYZTLorentzVectorD const& _position,
                                 float _charge,
                                 math::XYZVectorD& _output) const
{
  double px(_momentum.px());
  double py(_momentum.py());
  double pz(_momentum.pz());
  double x(_position.x());
  double y(_position.y());
  double z(_position.z());

  typedef TrajectoryStateOnSurface TSOS;
  GlobalPoint startingPosition(x, y, z);
  GlobalVector startingMomentum(px, py, pz);
  Plane::PlanePointer startingPlane(Plane::build(Plane::PositionType(x, y, z), Plane::RotationType()));

  TSOS startingStateP(GlobalTrajectoryParameters(startingPosition, startingMomentum, _charge, &field_),
                      err_,
                      *startingPlane);

  TSOS trackStateP(prod_.propagator.propagate(startingStateP, *(pz > 0 ? targetPlaneForward_ : targetPlaneBackward_)));

  if (trackStateP.isValid()) {
    auto&& gp(trackStateP.globalPosition());
    _output.SetCoordinates(gp.x(), gp.y(), gp.z());
    return true;
  }

  _output = math::XYZVectorD();
  return false;
}


HGCalTriggerNtupleHGCMulticlusters::
HGCalTriggerNtupleHGCMulticlusters(const edm::ParameterSet& conf):HGCalTriggerNtupleBase(conf), propagator_(nullptr)
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

    tree.Branch(withPrefix("simtracks_pid"), &cl3d_simtracks_pid_);
    tree.Branch(withPrefix("simtracks_pt"), &cl3d_simtracks_pt_);
    tree.Branch(withPrefix("simtracks_eta"), &cl3d_simtracks_eta_);
    tree.Branch(withPrefix("simtracks_phi"), &cl3d_simtracks_phi_);

    fsimEvent_ = new FSimEvent(conf.getUntrackedParameterSet("SimTrackFilter"));
  }
}

void
HGCalTriggerNtupleHGCMulticlusters::
beginRun(const edm::Run&, const edm::EventSetup& es)
{
  if (matchSimTrack_) {
    edm::ESHandle<HepPDT::ParticleDataTable> pdt;
    es.getData(pdt);
    fsimEvent_->initializePdt(&(*pdt));

    edm::ESHandle<MagneticField> magfield;
    es.get<IdealMagneticFieldRecord>().get(magfield);

    hgcal::RecHitTools recHitTools;
    recHitTools.getEventSetup(es);

    delete propagator_;
    propagator_ = new SimpleTrackPropagator(*magfield, recHitTools.getPositionLayer(1).z());
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

  // retrieve geometry
  edm::ESHandle<HGCalTriggerGeometryBase> geometry;
  es.get<CaloGeometryRecord>().get(geometry);

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
    // approximate bounds
    double const outerRadius(160.);
    double const innerRadius(25.);

    // vector of propagated positions and simtracks
    std::vector<std::pair<math::XYZVectorD, FSimTrack const*>> positionsAndTracks;

    edm::Handle<std::vector<SimTrack>> simTracksHandle;
    e.getByToken(simTracksToken_, simTracksHandle);
    auto const& simTracks(*simTracksHandle);

    edm::Handle<std::vector<SimVertex>> simVerticesHandle;
    e.getByToken(simVerticesToken_, simVerticesHandle);
    auto const& simVertices(*simVerticesHandle);

    fsimEvent_->fill(simTracks, simVertices);

    for (unsigned iT(0); iT != fsimEvent_->nTracks(); ++iT) {
      auto const& track(fsimEvent_->track(iT));

      // skip tracks that start after HGCal
      if (std::abs(track.vertex().position().z()) > propagator_->abszTarget())
        continue;

      // skip tracks that end before HGCal
      if (!track.noEndVertex() && std::abs(track.endVertex().position().z()) < propagator_->abszTarget())
        continue;

      math::XYZVectorD positionOnSurface;
      if (!propagator_->propagate(track.momentum(), track.vertex().position(), track.charge(), positionOnSurface))
        continue;

      double rho(positionOnSurface.Rho());
      if (rho < outerRadius && rho > innerRadius)
        positionsAndTracks.emplace_back(positionOnSurface, &track);
    }

    // now match the tracks to clusters
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
  cl3d_simtracks_pt_.clear();
  cl3d_simtracks_eta_.clear();
  cl3d_simtracks_phi_.clear();
}




