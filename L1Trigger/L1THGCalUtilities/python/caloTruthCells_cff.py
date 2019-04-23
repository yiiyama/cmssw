import FWCore.ParameterSet.Config as cms

from L1Trigger.L1THGCalUtilities.caloTruthCellsProducer_cfi import *

caloTruthCells = cms.Sequence(caloTruthCellsProducer)

if caloTruthCellsProducer.makeCellsCollection:
    ## cluster and tower sequence

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
    
    caloTruthCells += cms.Sequence(
        hgcalTruthConcentratorProducer *
        hgcalTruthBackEndLayer1Producer *
        hgcalTruthBackEndLayer2Producer *
        hgcalTruthTowerMapProducer *
        hgcalTruthTowerProducer
    )

