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
#include "DataFormats/ForwardDetId/interface/HGCalDetId.h"
#include "SimDataFormats/CaloAnalysis/interface/CaloParticleFwd.h"
#include "SimDataFormats/CaloAnalysis/interface/CaloParticle.h"
#include "SimDataFormats/CaloAnalysis/interface/SimCluster.h"
#include "L1Trigger/L1THGCal/interface/HGCalTriggerGeometryBase.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
// #include "Geometry/HcalCommonData/interface/HcalDDDRecConstants.h"
// #include "Geometry/HcalCommonData/interface/HcalHitRelabeller.h"
// #include "SimDataFormats/CaloTest/interface/HGCalTestNumbering.h"

class CaloTruthCellsProducer : public edm::stream::EDProducer<> {
public:
  explicit CaloTruthCellsProducer(edm::ParameterSet const&);
  ~CaloTruthCellsProducer();

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  virtual void beginStream(edm::StreamID) override;
  virtual void produce(edm::Event&, edm::EventSetup const&) override;
  virtual void endStream() override;

  virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
  virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
  virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
  virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

  // std::map<uint32_t, double> makeHitMap(edm::Event const&, edm::EventSetup const&) const;

  typedef edm::AssociationMap<edm::OneToMany<CaloParticleCollection, l1t::HGCalTriggerCellBxCollection>> CaloToCellsMap;

  edm::EDGetTokenT<CaloParticleCollection> caloParticlesToken_;
  edm::EDGetTokenT<l1t::HGCalTriggerCellBxCollection> triggerCellsToken_;
  // // three simhits tokens, depedent on old or new geometry
  // edm::EDGetTokenT<std::vector<PCaloHit>> simHitsTokenEE_;
  // edm::EDGetTokenT<std::vector<PCaloHit>> simHitsTokenHEfront_;
  // edm::EDGetTokenT<std::vector<PCaloHit>> simHitsTokenHEback_;

  // // for input with old (<= V8) geometry, sim and reco hits use different ID schemes
  // bool needSimToRecoMapping_{false};
};

CaloTruthCellsProducer::CaloTruthCellsProducer(edm::ParameterSet const& _config) :
  caloParticlesToken_(consumes<CaloParticleCollection>(_config.getParameter<edm::InputTag>("caloParticles"))),
  triggerCellsToken_(consumes<l1t::HGCalTriggerCellBxCollection>(_config.getParameter<edm::InputTag>("triggerCells")))
  // simHitsTokenEE_(consumes<std::vector<PCaloHit>>(_config.getParameter<edm::InputTag>("simHitsEE"))),
  // simHitsTokenHEfront_(consumes<std::vector<PCaloHit>>(_config.getParameter<edm::InputTag>("simHitsHEfront"))),
  // simHitsTokenHEback_(consumes<std::vector<PCaloHit>>(_config.getParameter<edm::InputTag>("simHitsHEback")))
{
  produces<CaloToCellsMap>();
}

CaloTruthCellsProducer::~CaloTruthCellsProducer()
{
}

void
CaloTruthCellsProducer::produce(edm::Event& _event, edm::EventSetup const& _setup)
{
  edm::Handle<CaloParticleCollection> caloParticlesHandle;
  _event.getByToken(caloParticlesToken_, caloParticlesHandle);
  auto& caloParticles(*caloParticlesHandle);

  edm::Handle<l1t::HGCalTriggerCellBxCollection> triggerCellsHandle;
  _event.getByToken(triggerCellsToken_, triggerCellsHandle);
  auto& triggerCells(*triggerCellsHandle);

  // std::map<uint32_t, double> hitMap(makeHitMap(_event, _setup)); // cellId -> sim energy
  // std::map<uint32_t, double> tcMap; // tcId -> total sim energy

  edm::ESHandle<HGCalTriggerGeometryBase> geometry;
  _setup.get<CaloGeometryRecord>().get(geometry);

  std::map<uint32_t, CaloParticleRef> tcToCalo;

  for (unsigned iP(0); iP != caloParticles.size(); ++iP) {
    auto& caloParticle(caloParticles.at(iP));
    if (caloParticle.g4Tracks().at(0).eventId().event() != 0) // pileup
      continue;

    CaloParticleRef ref(caloParticlesHandle, iP);

    SimClusterRefVector const& simClusters(caloParticle.simClusters());
    for (auto& simCluster : simClusters) {
      for (auto& hAndF : simCluster->hits_and_fractions()) {
        DetId hitId(hAndF.first);
        HGCalDetId tcId;
        try {
          tcId = HGCalDetId(geometry->getTriggerCellFromCell(hitId));
        }
        catch (cms::Exception const& ex) {
          edm::LogError("CaloTruthCellsProducer") << ex.what();
          continue;
        }

        // tcMap[tcId] += hitMap[hitId.rawId()] * hAndF.first;
        tcToCalo.emplace(tcId, ref);
      }
    }
  }

  auto outMap(std::make_unique<CaloToCellsMap>(caloParticlesHandle, triggerCellsHandle));

  // loop through all bunch crossings
  for (unsigned iC(0); iC != triggerCells.size(); ++iC) {
    auto& cell(triggerCells[iC]);

    auto mapElem(tcToCalo.find(cell.detId()));
    if (mapElem == tcToCalo.end())
      continue;

    edm::Ref<l1t::HGCalTriggerCellBxCollection> ref(triggerCellsHandle, iC);

    outMap->insert(mapElem->second, ref);
  }

  _event.put(std::move(outMap));
}

void
CaloTruthCellsProducer::beginStream(edm::StreamID)
{
}

void
CaloTruthCellsProducer::endStream() {
}
void
CaloTruthCellsProducer::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
 
void
CaloTruthCellsProducer::endRun(edm::Run const&, edm::EventSetup const&)
{
}
 
void
CaloTruthCellsProducer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const& _setup)
{
  // check geometry version

  // edm::ESHandle<CaloGeometry> geom;
  // _setup.get<CaloGeometryRecord>().get(geom);
  // HGCalGeometry const* eegeom(static_cast<const HGCalGeometry*>(geom->getSubdetectorGeometry(DetId::HGCalEE,ForwardSubdetector::ForwardEmpty)));

  // if (eegeom)
  //   needSimToRecoMapping_ = false;
  // else
  //   needSimToRecoMapping_ = true;
}

void
CaloTruthCellsProducer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// void
// CaloTruthCellsProducer::makeHitMap(edm::Event const& _event, edm::EventSetup const& _setup) const
// {
//   std::map<uint32_t, double> hitMap; // cellId -> sim energy

//   if (needSimToRecoMapping_) {
//     // Conversion code from SimGeneral/CaloAnalysis/plugins/CaloTruthAccumulator.cc

//     edm::ESHandle<CaloGeometry> calogeom;
//     _setup.get<CaloGeometryRecord>().get(calogeom);

//     HcalGeometry const* bhgeom(static_cast<HcalGeometry const*>(calogeom->getSubdetectorGeometry(DetId::Hcal, HcalEndcap)));
//     HcalDDDRecConstants const* hcddd(bhgeom->topology().dddConstants());

//     HGCalGeometry const* eegeom(static_cast<HGCalGeometry const*>(calogeom->getSubdetectorGeometry(DetId::Forward, HGCEE)));
//     HGCalDDDConstants const* hgeddd(&eegeom->topology().dddConstants());

//     HGCalGeometry const* fhgeom(static_cast<HGCalGeometry const*>(calogeom->getSubdetectorGeometry(DetId::Forward, HGCHEF)));
//     HGCalDDDConstants const* hghddd(&fhgeom->topology().dddConstants());

//     for (auto* token : {&simHitsTokenEE_, &simHitsTokenHEfront_}) {
//       edm::Handle<std::vector<PCaloHit>> handle;
//       _event.getByToken(*token, handle);
//       const std::vector<PCaloHit>& simhits(*handle);

//       for (auto& simhit : simhits) {
//         // skip simhits with bad barcodes
//         if (simhit.geantTrackId() == 0)
//           continue;

//         int subdet, layer, cell, sec, subsec, zp;
//         HGCalTestNumbering::unpackHexagonIndex(simhit.id(), subdet, zp, layer, sec, subsec, cell);

//         HGCalDDDConstants const* ddd{nullptr};
//         switch (subdet) {
//         case HGCEE:
//           ddd = hgeddd;
//           break;
//         case HGCHEF:
//           ddd = hghddd;
//           break;
//         default:
//           throw cms::Exception("LogicError")
//             << "Invalid unpacked hexagon index: subdet = " << subdet;
//         }

//         // last argument is topology->detectorType() which is identically false
//         std::pair<int, int> recoLayerCell(ddd->simToReco(cell, layer, sec, false));
//         cell = recoLayerCell.first;
//         layer = recoLayerCell.second;
//         // skip simhits with non-existant layers
//         if (layer == -1)
//           continue;

//         uint32_t hitId(HGCalDetId(ForwardSubdetector(subdet), zp, layer, subsec, sec, cell).rawId());

//         hitMap.emplace(hitId, simhit.energy());
//       }
//     }

//     {
//       edm::Handle<std::vector<PCaloHit>> handle;
//       _event.getByToken(simHitsTokenHEback_, handle);
//       const std::vector<PCaloHit>& simhits(*handle);

//       for (auto& simhit : simhits) {
//         // Using HCAL DetId
//         HcalDetId hcalId(HcalHitRelabeller::relabel(simhit.id(), hcddd));
//         if (hcalId.subdet() != HcalEndcap)
//           continue;
      
//         hitMap.emplace(hcalId.rawId(), simhit.energy());
//       }
//     }
//   }
//   else {
//     edm::ESHandle<HGCalTriggerGeometryBase> geometry;
//     _setup.get<CaloGeometryRecord>().get(geometry);

//     for (auto* token : {&simHitsTokenEE_, &simHitsTokenHEfront_, &simHitsTokenHEback_}) {
//       edm::Handle<std::vector<PCaloHit>> handle;
//       _event.getByToken(*token, handle);
//       auto& simhits(*handle);

//       for (auto& simhit : simhits)
//         hitMap.emplace(simhit.id(), simhit.energy());
//     }
//   }

//   return hitMap;
// }

void
CaloTruthCellsProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

DEFINE_FWK_MODULE(CaloTruthCellsProducer);
