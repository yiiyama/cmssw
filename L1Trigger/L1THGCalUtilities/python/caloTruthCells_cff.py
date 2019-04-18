import FWCore.ParameterSet.Config as cms

from L1Trigger.L1THGCalUtilities.caloTruthCellsProducer_cfi import *
from L1Trigger.L1THGCalUtilities.hgcalTriggerNtuples_cfi import *

caloTruthCells = cms.Sequence(caloTruthCellsProducer)

# Add the "full truth" multiclusters (each cluster made with all cells belonging to one gen particle) to the ntuples

ntuple_multiclusters_fulltruth = ntuple_multiclusters.clone(
    Multiclusters = cms.InputTag('caloTruthCellsProducer:HGCalBackendLayer2Processor3DClustering'),
    Prefix = cms.untracked.string('cl3dfulltruth')
)
hgcalTriggerNtuplizer.Ntuples.append(ntuple_multiclusters_fulltruth)

# If caloTruthCellsProducer.makeCellsCollection is True, can run the clustering algorithm over output cells too

if caloTruthCellsProducer.makeCellsCollection:
    ## cluster and tower sequence

    from L1Trigger.L1THGCal.hgcalTriggerGeometryESProducer_cfi import *
    from L1Trigger.L1THGCal.hgcalConcentrator_cff import *
    from L1Trigger.L1THGCal.hgcalBackEndLayer1_cff import *
    from L1Trigger.L1THGCal.hgcalBackEndLayer2_cff import *
    from L1Trigger.L1THGCal.hgcalTowerMap_cff import *
    from L1Trigger.L1THGCal.hgcalTower_cff import *
    
    hgcalTruthConcentratorProducer = hgcalConcentratorProducer.clone(
        InputTriggerCells = cms.InputTag('caloTruthCellsProducer:HGCalVFEProcessorSums')
    )
    
    hgcalTruthBackEndLayer1Producer = hgcalBackEndLayer1Producer.clone(
        InputTriggerCells = cms.InputTag('hgcalTruthConcentratorProducer:HGCalConcentratorProcessorSelection')
    )
    
    hgcalTruthBackEndLayer2Producer = hgcalBackEndLayer2Producer.clone(
        InputCluster = cms.InputTag('hgcalTruthBackEndLayer1Producer:HGCalBackendLayer1Processor2DClustering')
    )
    
    hgcalTruthTowerMapProducer = hgcalTowerMapProducer.clone(
        InputTriggerCells = cms.InputTag('caloTruthCellsProducer:HGCalVFEProcessorSums')
    )
    
    hgcalTruthTowerProducer = hgcalTowerProducer.clone(
        InputTowerMaps = cms.InputTag('hgcalTruthTowerMapProducer:HGCalTowerMapProcessor')
    )
    
    caloTruthCells += cms.Sequence(
        hgcalTruthConcentratorProducer *
        hgcalTruthBackEndLayer1Producer *
        hgcalTruthBackEndLayer2Producer *
        hgcalTruthTowerMapProducer *
        hgcalTruthTowerProducer
    )

    ## ntuplize clusters and towers

    ntuple_triggercells.caloParticlesToCells = cms.InputTag('caloTruthCellsProducer')
    ntuple_triggercells.FillTruthMap = cms.bool(True)
    
    ntuple_triggercells_truth = ntuple_triggercells.clone(
        TriggerCells = cms.InputTag('caloTruthCellsProducer:HGCalVFEProcessorSums'),
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


from L1Trigger.L1THGCal.hgcalTriggerPrimitives_cff import hgcalTriggerPrimitives
hgcalTriggerPrimitives += caloTruthCells
