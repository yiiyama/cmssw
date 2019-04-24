import FWCore.ParameterSet.Config as cms

import SimCalorimetry.HGCalSimProducers.hgcalDigitizer_cfi as digiparam
import RecoLocalCalo.HGCalRecProducers.HGCalUncalibRecHit_cfi as recoparam
import RecoLocalCalo.HGCalRecProducers.HGCalRecHit_cfi as recocalibparam

from L1Trigger.L1THGCal.egammaIdentification import egamma_identification_drnn_cone, \
                                                    egamma_identification_drnn_dbscan, \
                                                    egamma_identification_histomax

from Configuration.Eras.Modifier_phase2_hgcalV9_cff import phase2_hgcalV9


binSums = cms.vuint32(13,               # 0
                      11, 11, 11,       # 1 - 3
                      9, 9, 9,          # 4 - 6
                      7, 7, 7, 7, 7, 7,  # 7 - 12
                      5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,  # 13 - 27
                      3, 3, 3, 3, 3, 3, 3, 3  # 28 - 35
                      )

EE_DR_GROUP = 7
FH_DR_GROUP = 6
BH_DR_GROUP = 12
MAX_LAYERS = 52

dr_layerbylayer = ([0] + # no layer 0
        [0.015]*EE_DR_GROUP + [0.020]*EE_DR_GROUP + [0.030]*EE_DR_GROUP + [0.040]*EE_DR_GROUP + # EM
        [0.040]*FH_DR_GROUP + [0.050]*FH_DR_GROUP + # FH
        [0.050]*BH_DR_GROUP) # BH


dr_layerbylayer_Bcoefficient = ([0] + # no layer 0
        [0.020]*EE_DR_GROUP + [0.020]*EE_DR_GROUP + [0.02]*EE_DR_GROUP + [0.020]*EE_DR_GROUP + # EM
        [0.020]*FH_DR_GROUP + [0.020]*FH_DR_GROUP + # FH
        [0.020]*BH_DR_GROUP) # BH


neighbour_weights_1stOrder = cms.vdouble(0, 0.25, 0,
                                         0.25, 0, 0.25,
                                         0, 0.25, 0)

neighbour_weights_2ndOrder = cms.vdouble(-0.25, 0.5, -0.25,
                                         0.5, 0,  0.5,
                                         -0.25, 0.5, -0.25)


distance_C3d_params = cms.PSet(type_multicluster=cms.string('dRC3d'),
                               dR_multicluster=cms.double(0.01),
                               minPt_multicluster=cms.double(0.5),  # minimum pt of the multicluster (GeV)
                               dist_dbscan_multicluster=cms.double(0.),
                               minN_dbscan_multicluster=cms.uint32(0),
                               EGIdentification=egamma_identification_drnn_cone.clone(),
                               )


dbscan_C3d_params = cms.PSet(type_multicluster=cms.string('DBSCANC3d'),
                             dR_multicluster=cms.double(0.),
                             minPt_multicluster=cms.double(0.5),  # minimum pt of the multicluster (GeV)
                             dist_dbscan_multicluster=cms.double(0.005),
                             minN_dbscan_multicluster=cms.uint32(3),
                             EGIdentification=egamma_identification_drnn_dbscan.clone())


histoMax_C3d_params = cms.PSet(type_multicluster=cms.string('HistoMaxC3d'),
                               dR_multicluster=cms.double(0.03),
                               dR_multicluster_byLayer_coefficientA=cms.vdouble(),
                               dR_multicluster_byLayer_coefficientB=cms.vdouble(),
                               minPt_multicluster=cms.double(0.5),  # minimum pt of the multicluster (GeV)
                               nBins_R_histo_multicluster=cms.uint32(36),
                               nBins_Phi_histo_multicluster=cms.uint32(216),
                               binSumsHisto=binSums,
                               threshold_histo_multicluster=cms.double(10.),
                               cluster_association=cms.string("NearestNeighbour"),
                               EGIdentification=egamma_identification_histomax.clone(),
                               neighbour_weights=neighbour_weights_1stOrder
                               )
# V9 samples have a different defintiion of the dEdx calibrations. To account for it
# we reascale the thresholds of the clustering seeds
# (see https://indico.cern.ch/event/806845/contributions/3359859/attachments/1815187/2966402/19-03-20_EGPerf_HGCBE.pdf
# for more details)
phase2_hgcalV9.toModify(histoMax_C3d_params,
                        threshold_histo_multicluster=cms.double(7.5),  # MipT
                        )


histoMaxVariableDR_C3d_params = histoMax_C3d_params.clone()
histoMaxVariableDR_C3d_params.dR_multicluster = cms.double(0.)
histoMaxVariableDR_C3d_params.dR_multicluster_byLayer_coefficientA = cms.vdouble(dr_layerbylayer)
histoMaxVariableDR_C3d_params.dR_multicluster_byLayer_coefficientB = cms.vdouble([0]*53)


histoSecondaryMax_C3d_params = histoMax_C3d_params.clone()
histoSecondaryMax_C3d_params.type_multicluster = cms.string('HistoSecondaryMaxC3d')


histoInterpolatedMax_C3d_params = histoMax_C3d_params.clone()
histoInterpolatedMax_C3d_params.type_multicluster = cms.string('HistoInterpolatedMaxC3d')


histoThreshold_C3d_params = histoMax_C3d_params.clone()
histoThreshold_C3d_params.type_multicluster = cms.string('HistoThresholdC3d')


be_proc = cms.PSet(ProcessorName  = cms.string('HGCalBackendLayer2Processor3DClustering'),
                   C3d_parameters = histoMaxVariableDR_C3d_params.clone()
                   )

hgcalBackEndLayer2Producer = cms.EDProducer(
    "HGCalBackendLayer2Producer",
    InputCluster = cms.InputTag('hgcalBackEndLayer1Producer:HGCalBackendLayer1Processor2DClustering'),
    ProcessorParameters = be_proc.clone()
    )
