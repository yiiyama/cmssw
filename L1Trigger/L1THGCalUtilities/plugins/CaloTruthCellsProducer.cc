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

class CaloTruthCellsProducer : public edm::stream::EDProducer<> {
public:
  explicit CaloTruthCellsProducer(edm::ParameterSet const&);
  ~CaloTruthCellsProducer();

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  virtual void produce(edm::Event&, edm::EventSetup const&) override;

  typedef edm::AssociationMap<edm::OneToMany<CaloParticleCollection, l1t::HGCalTriggerCellBxCollection>> CaloToCellsMap;

  bool makeCellsCollection_;
  edm::EDGetTokenT<CaloParticleCollection> caloParticlesToken_;
  edm::EDGetTokenT<l1t::HGCalTriggerCellBxCollection> triggerCellsToken_;

  HGCalClusteringDummyImpl dummyClustering_;
  HGCalShowerShape showerShape_;
};

CaloTruthCellsProducer::CaloTruthCellsProducer(edm::ParameterSet const& config) :
  makeCellsCollection_(config.getParameter<bool>("makeCellsCollection")),
  caloParticlesToken_(consumes<CaloParticleCollection>(config.getParameter<edm::InputTag>("caloParticles"))),
  triggerCellsToken_(consumes<l1t::HGCalTriggerCellBxCollection>(config.getParameter<edm::InputTag>("triggerCells"))),
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
  auto& caloParticles(*caloParticlesHandle);

  edm::Handle<l1t::HGCalTriggerCellBxCollection> triggerCellsHandle;
  event.getByToken(triggerCellsToken_, triggerCellsHandle);
  auto& triggerCells(*triggerCellsHandle);

  dummyClustering_.eventSetup(setup);
  showerShape_.eventSetup(setup);

  edm::ESHandle<HGCalTriggerGeometryBase> geometryHandle;
  setup.get<CaloGeometryRecord>().get(geometryHandle);
  auto& geometry(*geometryHandle);

  std::unordered_map<uint32_t, CaloParticleRef> tcToCalo;

  // used later to order multiclusters
  std::map<int, CaloParticleRef> orderedCaloRefs;

  for (unsigned iP(0); iP != caloParticles.size(); ++iP) {
    auto& caloParticle(caloParticles.at(iP));
    if (caloParticle.g4Tracks().at(0).eventId().event() != 0) // pileup
      continue;

    CaloParticleRef ref(caloParticlesHandle, iP);

    SimClusterRefVector const& simClusters(caloParticle.simClusters());
    for (auto& simCluster : simClusters) {
      for (auto& hAndF : simCluster->hits_and_fractions()) {
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
      auto& cell(*cItr);

      auto mapElem(tcToCalo.find(cell.detId()));
      if (mapElem == tcToCalo.end())
        continue;

      outMap->insert(mapElem->second, edm::Ref<l1t::HGCalTriggerCellBxCollection>(triggerCellsHandle, triggerCells.key(cItr)));

      if (makeCellsCollection_)
        outCollection->push_back(bx, cell);

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
    auto& cluster((*outClusters)[iC]);
    // cluster detId and cell detId are identical
    auto caloRef(tcToCalo.at(cluster.detId()));
    caloToClusterIndices[caloRef.key()].push_back(iC);
  }

  auto&& clustersHandle(event.put(std::move(outClusters)));

  auto outMulticlusters(std::make_unique<l1t::HGCalMulticlusterBxCollection>());

  for (auto& ocr : orderedCaloRefs) {
    auto& ref(ocr.second);
    
    if (ref.isNull()) // shouldn't happen
      continue;

    auto& caloParticle(*ref);

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
CaloTruthCellsProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

DEFINE_FWK_MODULE(CaloTruthCellsProducer);
