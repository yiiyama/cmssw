import FWCore.ParameterSet.Config as cms

from L1Trigger.L1THGCalUtilities.caloTruthCellsProducer_cfi import *

from L1Trigger.L1THGCal.hgcalTriggerGeometryESProducer_cfi import *
from L1Trigger.L1THGCal.hgcalConcentrator_cff import *
from L1Trigger.L1THGCal.hgcalBackEndLayer1_cff import *
from L1Trigger.L1THGCal.hgcalBackEndLayer2_cff import *
from L1Trigger.L1THGCal.hgcalTowerMap_cff import *
from L1Trigger.L1THGCal.hgcalTower_cff import *

hgcalTruthConcentratorProducer = hgcalConcentratorProducer.clone(
    InputTriggerCells = cms.InputTag('caloTruthCellsProducer')
)

hgcalTruthBackEndLayer1Producer = hgcalBackEndLayer1Producer.clone(
    InputTriggerCells = cms.InputTag('hgcalTruthConcentratorProducer:HGCalConcentratorProcessorSelection')
)

hgcalTruthBackEndLayer2Producer = hgcalBackEndLayer2Producer.clone(
    InputCluster = cms.InputTag('hgcalTruthBackEndLayer1Producer:HGCalBackendLayer1Processor2DClustering')
)

hgcalTruthTowerMapProducer = hgcalTowerMapProducer.clone(
    InputTriggerCells = cms.InputTag('caloTruthCellsProducer')
)

hgcalTruthTowerProducer = hgcalTowerProducer.clone(
    InputTowerMaps = cms.InputTag('hgcalTruthTowerMapProducer:HGCalTowerMapProcessor')
)

caloTruthCells = cms.Sequence(
    caloTruthCellsProducer *
    hgcalTruthConcentratorProducer *
    hgcalTruthBackEndLayer1Producer *
    hgcalTruthBackEndLayer2Producer *
    hgcalTruthTowerMapProducer *
    hgcalTruthTowerProducer
)

from L1Trigger.L1THGCalUtilities.hgcalTriggerNtuples_cfi import *

ntuple_triggercells.caloParticlesToCells = cms.InputTag('caloTruthCellsProducer')
ntuple_triggercells.FillTruthMap = cms.bool(True)

ntuple_triggercells_truth = ntuple_triggercells.clone(
    TriggerCells = cms.InputTag('caloTruthCellsProducer'),
    MultiClusters = cms.InputTag('hgcalTruthBackEndLayer2Producer:HGCalBackendLayer1Processor2DClustering'),
    Prefix = cms.untracked.string('tctruth'),
    FillTruthMap = cms.bool(False)
)

ntuple_clusters_truth = ntuple_clusters.clone(
    Clusters = cms.InputTag('hgcalTruthBackEndLayer1Producer:HGCalBackendLayer1Processor2DClustering'),
    Prefix = cms.untracked.string('cltruth')
)

ntuple_multiclusters_truth = ntuple_multiclusters.clone(
    Multiclusters = cms.InputTag('hgcalTruthBackEndLayer2Producer:HGCalBackendLayer2Processor3DClustering'),
    Prefix = cms.untracked.string('cl3dtruth')
)

ntuple_towers_truth = ntuple_towers.clone(
    Towers = cms.InputTag('hgcalTruthTowerProducer:HGCalTowerProcessor'),
    Prefix = cms.untracked.string('towertruth')
)

hgcalTriggerNtuplizer.Ntuples.append(ntuple_triggercells_truth)
hgcalTriggerNtuplizer.Ntuples.append(ntuple_multiclusters_truth)
hgcalTriggerNtuplizer.Ntuples.append(ntuple_towers_truth)
