#include <memory>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"
#include "DataFormats/Common/interface/AssociationMap.h"
#include "DataFormats/Common/interface/OneToMany.h"
#include "DataFormats/Common/interface/Ref.h"
#include "DataFormats/L1THGCal/interface/HGCalTriggerCell.h"
#include "DataFormats/L1THGCal/interface/HGCalMulticluster.h"
#include "DataFormats/ForwardDetId/interface/HGCalDetId.h"
#include "SimDataFormats/CaloAnalysis/interface/CaloParticleFwd.h"
#include "SimDataFormats/CaloAnalysis/interface/CaloParticle.h"
#include "SimDataFormats/CaloAnalysis/interface/SimCluster.h"
#include "L1Trigger/L1THGCal/interface/HGCalTriggerGeometryBase.h"
#include "L1Trigger/L1THGCal/interface/backend/HGCalShowerShape.h"
#include "L1Trigger/L1THGCal/interface/backend/HGCalClusteringDummyImpl.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "Geometry/HcalCommonData/interface/HcalDDDRecConstants.h"
#include "Geometry/HcalCommonData/interface/HcalHitRelabeller.h"
#include "SimDataFormats/CaloTest/interface/HGCalTestNumbering.h"

class CaloTruthCellsProducer : public edm::stream::EDProducer<> {
public:
  explicit CaloTruthCellsProducer(edm::ParameterSet const&);
  ~CaloTruthCellsProducer();

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void produce(edm::Event&, edm::EventSetup const&) override;
  void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

  std::unordered_map<uint32_t, double> makeHitMap(edm::Event const&, edm::EventSetup const&) const;

  typedef edm::AssociationMap<edm::OneToMany<CaloParticleCollection, l1t::HGCalTriggerCellBxCollection>> CaloToCellsMap;

  bool makeCellsCollection_;
  edm::EDGetTokenT<CaloParticleCollection> caloParticlesToken_;
  edm::EDGetTokenT<l1t::HGCalTriggerCellBxCollection> triggerCellsToken_;
  edm::EDGetTokenT<std::vector<PCaloHit>> simHitsTokenEE_;
  edm::EDGetTokenT<std::vector<PCaloHit>> simHitsTokenHEfront_;
  edm::EDGetTokenT<std::vector<PCaloHit>> simHitsTokenHEback_;

  HGCalClusteringDummyImpl dummyClustering_;
  HGCalShowerShape showerShape_;

  // for input with old (<= V8) geometry, sim and reco hits use different ID schemes
  // automatically detected in beginLuminosityBlock through the geometry object
  bool simHitsUseRecoDetIds_{true};
};

CaloTruthCellsProducer::CaloTruthCellsProducer(edm::ParameterSet const& config) :
  makeCellsCollection_(config.getParameter<bool>("makeCellsCollection")),
  caloParticlesToken_(consumes<CaloParticleCollection>(config.getParameter<edm::InputTag>("caloParticles"))),
  triggerCellsToken_(consumes<l1t::HGCalTriggerCellBxCollection>(config.getParameter<edm::InputTag>("triggerCells"))),
  simHitsTokenEE_(consumes<std::vector<PCaloHit>>(config.getParameter<edm::InputTag>("simHitsEE"))),
  simHitsTokenHEfront_(consumes<std::vector<PCaloHit>>(config.getParameter<edm::InputTag>("simHitsHEfront"))),
  simHitsTokenHEback_(consumes<std::vector<PCaloHit>>(config.getParameter<edm::InputTag>("simHitsHEback"))),
  dummyClustering_(config.getParameterSet("dummyClustering"))
{
  produces<CaloToCellsMap>();
  produces<l1t::HGCalClusterBxCollection>();
  produces<l1t::HGCalMulticlusterBxCollection>();
  if (makeCellsCollection_)
    produces<l1t::HGCalTriggerCellBxCollection>();
}

CaloTruthCellsProducer::~CaloTruthCellsProducer()
{
}

void
CaloTruthCellsProducer::produce(edm::Event& event, edm::EventSetup const& setup)
{
  edm::Handle<CaloParticleCollection> caloParticlesHandle;
  event.getByToken(caloParticlesToken_, caloParticlesHandle);
  auto const& caloParticles(*caloParticlesHandle);

  edm::Handle<l1t::HGCalTriggerCellBxCollection> triggerCellsHandle;
  event.getByToken(triggerCellsToken_, triggerCellsHandle);
  auto const& triggerCells(*triggerCellsHandle);

  dummyClustering_.eventSetup(setup);
  showerShape_.eventSetup(setup);

  edm::ESHandle<HGCalTriggerGeometryBase> geometryHandle;
  setup.get<CaloGeometryRecord>().get(geometryHandle);
  auto const& geometry(*geometryHandle);

  std::unordered_map<uint32_t, CaloParticleRef> tcToCalo;

  std::unordered_map<uint32_t, double> hitToEnergy(makeHitMap(event, setup)); // cellId -> sim energy
  std::unordered_map<uint32_t, std::pair<double, double>> tcToEnergies; // tcId -> {total sim energy, fractioned sim energy}

  // used later to order multiclusters
  std::map<int, CaloParticleRef> orderedCaloRefs;

  for (unsigned iP(0); iP != caloParticles.size(); ++iP) {
    auto const& caloParticle(caloParticles.at(iP));
    if (caloParticle.g4Tracks().at(0).eventId().event() != 0) // pileup
      continue;

    CaloParticleRef ref(caloParticlesHandle, iP);

    SimClusterRefVector const& simClusters(caloParticle.simClusters());
    for (auto const& simCluster : simClusters) {
      for (auto const& hAndF : simCluster->hits_and_fractions()) {
        DetId hitId(hAndF.first);
        uint32_t tcId;
        try {
          tcId = geometry.getTriggerCellFromCell(hitId);
        }
        catch (cms::Exception const& ex) {
          edm::LogError("CaloTruthCellsProducer") << ex.what();
          continue;
        }

        tcToCalo.emplace(tcId, ref);

        double cellE(hitToEnergy[hAndF.first]);
        tcToEnergies[tcId].first += cellE;
        tcToEnergies[tcId].second += cellE * hAndF.second;
      }
    }

    // ordered by the gen particle index
    int genIndex(caloParticle.g4Tracks().at(0).genpartIndex() - 1);
    if (genIndex >= 0) // < 0 shouldn't happen
      orderedCaloRefs[genIndex] = ref;
  }

  auto outMap(std::make_unique<CaloToCellsMap>(caloParticlesHandle, triggerCellsHandle));
  std::unique_ptr<l1t::HGCalTriggerCellBxCollection> outCollection;
  if (makeCellsCollection_)
    outCollection = std::move(std::make_unique<l1t::HGCalTriggerCellBxCollection>());

  typedef edm::Ptr<l1t::HGCalTriggerCell> TriggerCellPtr;
  typedef edm::Ptr<l1t::HGCalCluster> ClusterPtr;

  // ClusteringDummyImpl only considers BX 0, so we dump all cells to one vector
  std::vector<TriggerCellPtr> triggerCellPtrs;

  // loop through all bunch crossings
  for (int bx(triggerCells.getFirstBX()); bx <= triggerCells.getLastBX(); ++bx) {
    for (auto&& cItr(triggerCells.begin(bx)); cItr != triggerCells.end(bx); ++cItr) {
      auto const& cell(*cItr);

      auto mapElem(tcToCalo.find(cell.detId()));
      if (mapElem == tcToCalo.end())
        continue;

      outMap->insert(mapElem->second, edm::Ref<l1t::HGCalTriggerCellBxCollection>(triggerCellsHandle, triggerCells.key(cItr)));

      if (makeCellsCollection_) {
        auto& simEnergies(tcToEnergies.at(cell.detId()));
        if (simEnergies.first > 0.) {
          outCollection->push_back(bx, cell);
          (*outCollection)[outCollection->size() - 1].setMipPt(cell.mipPt() * simEnergies.second / simEnergies.first);
        }
      }

      triggerCellPtrs.emplace_back(triggerCellsHandle, triggerCells.key(cItr));
    }
  }

  event.put(std::move(outMap));
  if (makeCellsCollection_)
    event.put(std::move(outCollection));

  auto outClusters(std::make_unique<l1t::HGCalClusterBxCollection>());

  auto sortCellPtrs([](TriggerCellPtr const& lhs, TriggerCellPtr const& rhs)->bool {
      return lhs->mipPt() > rhs->mipPt();
    });

  std::sort(triggerCellPtrs.begin(), triggerCellPtrs.end(), sortCellPtrs);
  dummyClustering_.clusterizeDummy(triggerCellPtrs, *outClusters);

  std::unordered_map<unsigned, std::vector<unsigned>> caloToClusterIndices;
  for (unsigned iC(0); iC != outClusters->size(); ++iC) {
    auto const& cluster((*outClusters)[iC]);
    // cluster detId and cell detId are identical
    auto caloRef(tcToCalo.at(cluster.detId()));
    caloToClusterIndices[caloRef.key()].push_back(iC);
  }

  auto&& clustersHandle(event.put(std::move(outClusters)));

  auto outMulticlusters(std::make_unique<l1t::HGCalMulticlusterBxCollection>());

  for (auto const& ocr : orderedCaloRefs) {
    auto const& ref(ocr.second);
    
    if (ref.isNull()) // shouldn't happen
      continue;

    auto const& caloParticle(*ref);

    l1t::HGCalMulticluster multicluster;
    
    for (unsigned iC : caloToClusterIndices[ref.key()]) {
      ClusterPtr clPtr(clustersHandle, iC);
      multicluster.addConstituent(clPtr, true, 1.);
    }

    // Set the gen particle index as the DetId
    multicluster.setDetId(caloParticle.g4Tracks().at(0).genpartIndex() - 1);

    auto&& centre(multicluster.centre());
    math::PtEtaPhiMLorentzVector multiclusterP4(multicluster.sumPt(), centre.eta(), centre.phi(), 0.);
    multicluster.setP4(multiclusterP4);

    showerShape_.fillShapes(multicluster, geometry);

    // not setting the quality flag
    // multicluster.setHwQual(id_->decision(multicluster));
    // fill H/E
    multicluster.saveHOverE();            
    
    outMulticlusters->push_back(0, multicluster);
  }

  event.put(std::move(outMulticlusters));
}

void
CaloTruthCellsProducer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const& setup)
{
  // check geometry version

  edm::ESHandle<CaloGeometry> geom;
  setup.get<CaloGeometryRecord>().get(geom);
  auto const* eegeom(geom->getSubdetectorGeometry(DetId::HGCalEE, ForwardSubdetector::ForwardEmpty));

  if (eegeom)
    simHitsUseRecoDetIds_ = true; // >= V9 geometry
  else
    simHitsUseRecoDetIds_ = false; // <= V8 geometry
}

std::unordered_map<uint32_t, double>
CaloTruthCellsProducer::makeHitMap(edm::Event const& event, edm::EventSetup const& setup) const
{
  std::unordered_map<uint32_t, double> hitMap; // cellId -> sim energy

  if (simHitsUseRecoDetIds_) {
    edm::ESHandle<HGCalTriggerGeometryBase> geometry;
    setup.get<CaloGeometryRecord>().get(geometry);

    for (auto* token : {&simHitsTokenEE_, &simHitsTokenHEfront_, &simHitsTokenHEback_}) {
      edm::Handle<std::vector<PCaloHit>> handle;
      event.getByToken(*token, handle);
      auto const& simhits(*handle);

      for (auto const& simhit : simhits)
        hitMap.emplace(simhit.id(), simhit.energy());
    }
  }
  else {
    // Conversion code from SimGeneral/CaloAnalysis/plugins/CaloTruthAccumulator.cc

    edm::ESHandle<CaloGeometry> calogeom;
    setup.get<CaloGeometryRecord>().get(calogeom);

    auto const* bhgeom(static_cast<HcalGeometry const*>(calogeom->getSubdetectorGeometry(DetId::Hcal, HcalEndcap)));
    HcalDDDRecConstants const* hcddd(bhgeom->topology().dddConstants());

    auto const* eegeom(static_cast<HGCalGeometry const*>(calogeom->getSubdetectorGeometry(DetId::Forward, HGCEE)));
    HGCalDDDConstants const* hgeddd(&eegeom->topology().dddConstants());

    auto const* fhgeom(static_cast<HGCalGeometry const*>(calogeom->getSubdetectorGeometry(DetId::Forward, HGCHEF)));
    HGCalDDDConstants const* hghddd(&fhgeom->topology().dddConstants());

    for (auto* token : {&simHitsTokenEE_, &simHitsTokenHEfront_}) {
      edm::Handle<std::vector<PCaloHit>> handle;
      event.getByToken(*token, handle);
      auto const& simhits(*handle);

      for (auto const& simhit : simhits) {
        // skip simhits with bad barcodes
        if (simhit.geantTrackId() == 0)
          continue;

        int subdet, layer, cell, sec, subsec, zp;
        HGCalTestNumbering::unpackHexagonIndex(simhit.id(), subdet, zp, layer, sec, subsec, cell);

        HGCalDDDConstants const* ddd{nullptr};
        switch (subdet) {
        case HGCEE:
          ddd = hgeddd;
          break;
        case HGCHEF:
          ddd = hghddd;
          break;
        default:
          throw cms::Exception("LogicError")
            << "Invalid unpacked hexagon index: subdet = " << subdet;
        }

        // last argument is topology->detectorType() which is identically false
        std::pair<int, int> recoLayerCell(ddd->simToReco(cell, layer, sec, false));
        cell = recoLayerCell.first;
        layer = recoLayerCell.second;
        // skip simhits with non-existant layers
        if (layer == -1)
          continue;

        uint32_t hitId(HGCalDetId(ForwardSubdetector(subdet), zp, layer, subsec, sec, cell).rawId());

        hitMap.emplace(hitId, simhit.energy());
      }
    }

    {
      edm::Handle<std::vector<PCaloHit>> handle;
      event.getByToken(simHitsTokenHEback_, handle);
      auto const& simhits(*handle);

      for (auto const& simhit : simhits) {
        // Using HCAL DetId
        HcalDetId hcalId(HcalHitRelabeller::relabel(simhit.id(), hcddd));
        if (hcalId.subdet() != HcalEndcap)
          continue;
      
        hitMap.emplace(hcalId.rawId(), simhit.energy());
      }
    }
  }

  return hitMap;
}

void
CaloTruthCellsProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

DEFINE_FWK_MODULE(CaloTruthCellsProducer);
