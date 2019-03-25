import FWCore.ParameterSet.Config as cms

caloTruthCells = cms.EDProducer('CaloTruthCellsProducer',
    caloParticles = cms.InputTag('mix', 'MergedCaloTruth'),
    triggerCells = cms.InputTag('hgcalConcentratorProducer:HGCalConcentratorProcessorSelection'),
    simHitsEE = cms.InputTag('g4SimHits:HGCHitsEE'),
    simHitsHEfront = cms.InputTag('g4SimHits:HGCHitsHEfront'),
    simHitsHEback = cms.InputTag('g4SimHits:HcalHits')
)
