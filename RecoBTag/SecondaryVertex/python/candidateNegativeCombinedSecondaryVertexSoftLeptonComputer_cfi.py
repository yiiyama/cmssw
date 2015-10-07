import FWCore.ParameterSet.Config as cms

from RecoBTag.SecondaryVertex.candidateCombinedSecondaryVertexSoftLeptonComputer_cfi import *

candidateNegativeCombinedSecondaryVertexSoftLeptonComputer = candidateCombinedSecondaryVertexSoftLeptonComputer.clone()
candidateNegativeCombinedSecondaryVertexSoftLeptonComputer.vertexFlip = True
candidateNegativeCombinedSecondaryVertexSoftLeptonComputer.trackFlip = True
candidateNegativeCombinedSecondaryVertexSoftLeptonComputer.trackSelection.sip3dSigMax = 0
candidateNegativeCombinedSecondaryVertexSoftLeptonComputer.trackPseudoSelection.sip3dSigMax = 0
candidateNegativeCombinedSecondaryVertexSoftLeptonComputer.trackPseudoSelection.sip2dSigMin = -99999.9
candidateNegativeCombinedSecondaryVertexSoftLeptonComputer.trackPseudoSelection.sip2dSigMax = -2.0
