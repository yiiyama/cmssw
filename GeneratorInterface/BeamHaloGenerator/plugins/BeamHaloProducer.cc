#include "GeneratorInterface/BeamHaloGenerator/interface/PYR.h"

#include "FWCore/Framework/interface/Run.h"
#include "FWCore/Framework/interface/LuminosityBlock.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/one/EDProducer.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ServiceRegistry/interface/RandomEngineSentry.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenRunInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"

#include "CLHEP/Random/RandomEngine.h"

#include "HepMC/IO_HEPEVT.h"
#include "HepMC/HEPEVT_Wrapper.h"
#include "HepMC/WeightContainer.h"
#include "HepMC/GenEvent.h"

#include <limits>
#include <iostream>
#include <array>
#include <vector>
#include <string>
#include <memory>

extern "C" {
  // fortran subroutines to be called from C++
  void bhsetparam_(int* iparam, float* fparam, const char* cparam, int length);
  void ki_bhg_init_(long& seed);
  void ki_bhg_fill_(int& iret, float& weight);
  void ki_bhg_stat_(int& iret);

  struct PRATES {
    double R_MU_P;
    double R_PI_P;
    double R_KA_P;
    double R_PROT;
    double R_NEUT;
    double R_MU_M;
    double R_PI_M;
    double R_KA_M;
    double NPRIME;
    double NSTACK;
  } prates_;
}

namespace CLHEP {
  class HepRandomEngine;
}

// Needs to be one module because we use FORTRAN globals - should port all from FORTRAN to C++ at some point
class BeamHaloProducer : public edm::one::EDProducer<edm::EndRunProducer, edm::one::WatchLuminosityBlocks, edm::one::SharedResources> {
 public:
  BeamHaloProducer(edm::ParameterSet const&);
  ~BeamHaloProducer();

  void setRandomEngine(CLHEP::HepRandomEngine*);

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

 private:
  void produce(edm::Event&, edm::EventSetup const&) override;
  void endRunProduce(edm::Run&, edm::EventSetup const&) override;
  void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
  void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override {}

  bool isInitialized_;
  double minR2_;
  double maxR2_;
};

BeamHaloProducer::BeamHaloProducer(edm::ParameterSet const& pset) :
  isInitialized_{false}
{
  // -- from bhgctrl.inc

  // param given in cm, HEPMC uses mm
  auto minR{pset.getParameter<double>("minR")};
  auto maxR{pset.getParameter<double>("maxR")};
  minR2_ = minR * minR * 100.;
  maxR2_ = maxR * maxR * 100.;

  std::array<int, 8> iparam{{
    pset.getParameter<int>("GENMOD"),
    pset.getParameter<int>("LHC_B1"),
    pset.getParameter<int>("LHC_B2"),
    pset.getParameter<int>("IW_MUO"),
    pset.getParameter<int>("IW_HAD"),
    std::numeric_limits<int>::max(),
    pset.getUntrackedParameter<int>("OFFSET", 0),
    pset.getParameter<int>("shift_bx")
  }};

  std::array<float, 4> fparam{{
    float(pset.getParameter<double>("EG_MIN")),
    float(pset.getParameter<double>("EG_MAX")),
    float(pset.getParameter<double>("BXNS")),
    float(pset.getUntrackedParameter<double>("W0", 1.0)),
  }};
  
  auto cparam{pset.getUntrackedParameter<std::string>("G3FNAME", "input.txt")};

  bhsetparam_(iparam.data(), fparam.data(), cparam.c_str(), cparam.length());

  produces<edm::HepMCProduct>("unsmeared");
  produces<GenEventInfoProduct>();
  produces<GenRunInfoProduct, edm::InRun>();
  produces<std::vector<double>, edm::InRun>("prates");

  usesResource("BeamHaloProducer");
}

BeamHaloProducer::~BeamHaloProducer()
{
  int iret{0};
  ki_bhg_stat_(iret);
}

void
BeamHaloProducer::setRandomEngine(CLHEP::HepRandomEngine* v)
{
  // Global defined in PYR.h
  _BeamHalo_randomEngine = v;
}

void
BeamHaloProducer::beginLuminosityBlock(edm::LuminosityBlock const& lumi, edm::EventSetup const&)
{
  if (!isInitialized_) {
    isInitialized_ = true;
    edm::RandomEngineSentry<BeamHaloProducer> randomEngineSentry{this, lumi.index()};

    // -- initialisation
    long seed{0}; // This seed is not actually used
    ki_bhg_init_(seed);
  }
}

void
BeamHaloProducer::produce(edm::Event& event, edm::EventSetup const&)
{
  edm::RandomEngineSentry<BeamHaloProducer> randomEngineSentry{this, event.streamID()};

  HepMC::HEPEVT_Wrapper wrapper;

  double sumWeight{0.};
  double eventWeight{1.};
  while (true) {
    int iret{0};
    float weight{0.};
    ki_bhg_fill_(iret, weight);
    sumWeight += weight;

    if (wrapper.number_entries() < 1)
      iret = -1;
    else if (minR2_ < maxR2_) {
      double R2(wrapper.x(1) * wrapper.x(1) + wrapper.y(1) * wrapper.y(1));
      if (R2 >= minR2_ && R2 <= maxR2_)
        eventWeight = weight / sumWeight;
      else if (iret >= 0)
        continue;
    }

    // Throw an exception if call_ki_bhg_fill(...) fails.  Use the EventCorruption
    // exception since it maps onto SkipEvent which is what we want to do here.

    if (iret < 0)
      throw edm::Exception(edm::errors::EventCorruption)
        << "BeamHaloProducer: function call_ki_bhg_fill returned " << iret;

    break;
  }

  // deleted in the HepMCProduct destructor
  auto genEvent{std::make_unique<HepMC::GenEvent>()};

  for (int idx{1}; idx <= wrapper.number_entries(); ++idx) {
    auto* vtx{new HepMC::GenVertex(HepMC::FourVector(wrapper.x(idx), wrapper.y(idx), wrapper.z(idx), wrapper.t(idx)))};
    HepMC::FourVector momentum{wrapper.px(idx), wrapper.py(idx), wrapper.pz(idx), wrapper.e(idx)};
    auto* part{new HepMC::GenParticle(momentum, wrapper.id(idx), wrapper.status(idx))};
    vtx->add_particle_out(part);
    genEvent->add_vertex(vtx);
  }

  genEvent->set_event_number(event.id().event());

  auto& weights{genEvent->weights()};
  weights.push_back(eventWeight);

  auto genEventInfo{std::make_unique<GenEventInfoProduct>(genEvent.get())};
  event.put(std::move(genEventInfo));

  auto product{std::make_unique<edm::HepMCProduct>()};
  product->addHepMCData(genEvent.release());
  event.put(std::move(product), "unsmeared");
}

void
BeamHaloProducer::endRunProduce(edm::Run &run, edm::EventSetup const&)
{
  // just create an empty product to keep the EventContent definitions happy
  // later on we might put the info into the run info that this is a PGun

  auto genRunInfo{std::make_unique<GenRunInfoProduct>()};
  run.put(std::move(genRunInfo));

  auto prates{std::make_unique<std::vector<double>>()};
  (*prates) = {
    prates_.R_MU_P * 2.,
    prates_.R_PI_P * 2.,
    prates_.R_KA_P * 2.,
    prates_.R_PROT * 2.,
    prates_.R_NEUT * 2.,
    prates_.R_MU_M * 2.,
    prates_.R_PI_M * 2.,
    prates_.R_KA_M * 2.
  };
  run.put(std::move(prates), "prates");
}

void
BeamHaloProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions)
{
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;

  desc.add<int>("GENMOD", 1); // mode (0 = "rate only", 1 = "MC generator", 2 = "read from file", 3 = ?
  desc.add<int>("LHC_B1", 1); // beam 1 on
  desc.add<int>("LHC_B2", 1); // beam 2 on
  desc.add<int>("IW_MUO", 1); // generate muons
  desc.add<int>("IW_HAD", 1); // generate hadrons
  desc.addUntracked<int>("OFFSET", 0); // number of particles to read when GENMOD = 3
  desc.add<int>("shift_bx", 0); // bunch crossing index wrt in-time

  desc.add<double>("EG_MIN", 10.); // minimum particle energy (GeV)
  desc.add<double>("EG_MAX", 14000.); // maximum particle energy (GeV)
  desc.add<double>("BXNS", 25.); // time between two bunch crossings (ns)
  desc.addUntracked<double>("W0", 1.0); // external per second nomrmalization

  desc.addUntracked<std::string>("G3FNAME", "input.txt"); // file name for GENMODE 3

  desc.add<double>("minR", 0.); // minimum transverse distance from origin (cm)
  desc.add<double>("maxR", 1000.); // maximum transverse distance from origin (cm)

  descriptions.add("beamHaloProducer", desc);
}

DEFINE_FWK_MODULE(BeamHaloProducer);
