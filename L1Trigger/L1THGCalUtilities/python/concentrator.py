import FWCore.ParameterSet.Config as cms
import SimCalorimetry.HGCalSimProducers.hgcalDigitizer_cfi as digiparam
from L1Trigger.L1THGCal.hgcalConcentratorProducer_cfi import threshold_conc_proc, best_conc_proc, supertc_conc_proc


def create_supertriggercell(process, inputs,
                            stcSize=supertc_conc_proc.stcSize
                            ):
    producer = process.hgcalConcentratorProducer.clone()
    producer.ProcessorParameters = supertc_conc_proc.clone()
    producer.ProcessorParameters.stcSize = stcSize
    producer.InputTriggerCells = cms.InputTag('{}:HGCalVFEProcessorSums'.format(inputs))
    producer.InputTriggerSums = cms.InputTag('{}:HGCalVFEProcessorSums'.format(inputs))
    return producer


def create_threshold(process, inputs,
                     threshold_silicon=threshold_conc_proc.triggercell_threshold_silicon,  # in mipT
                     threshold_scintillator=threshold_conc_proc.triggercell_threshold_scintillator  # in mipT
                     ):
    producer = process.hgcalConcentratorProducer.clone()
    producer.ProcessorParameters = threshold_conc_proc.clone()
    producer.ProcessorParameters.triggercell_threshold_silicon = threshold_silicon  # MipT
    producer.ProcessorParameters.triggercell_threshold_scintillator = threshold_scintillator  # MipT
    producer.InputTriggerCells = cms.InputTag('{}:HGCalVFEProcessorSums'.format(inputs))
    producer.InputTriggerSums = cms.InputTag('{}:HGCalVFEProcessorSums'.format(inputs))
    return producer


def create_bestchoice(process, inputs,
                      triggercells=best_conc_proc.NData
                      ):
    producer = process.hgcalConcentratorProducer.clone()
    producer.ProcessorParameters = best_conc_proc.clone()
    producer.ProcessorParameters.NData = triggercells
    producer.InputTriggerCells = cms.InputTag('{}:HGCalVFEProcessorSums'.format(inputs))
    producer.InputTriggerSums = cms.InputTag('{}:HGCalVFEProcessorSums'.format(inputs))
    return producer
