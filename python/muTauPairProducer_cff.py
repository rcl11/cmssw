import FWCore.ParameterSet.Config as cms
import copy

from TauAnalysis.CandidateTools.muTauPairSelector_cfi import *

allMuTauPairs = cms.EDProducer("PATMuTauPairProducer",
  useLeadingTausOnly = cms.bool(False),
  srcLeg1 = cms.InputTag('selectedLayer1MuonsTrkIPcumulative'),
  srcLeg2 = cms.InputTag('selectedLayer1TausForMuTauMuonVetoCumulative'),
  dRmin12 = cms.double(0.3),
  srcMET = cms.InputTag('allLayer1METs'),
  recoMode = cms.string(""),
  verbosity = cms.untracked.int32(0)
)

produceMuTauPairs = cms.Sequence( allMuTauPairs * selectMuTauPairs )


