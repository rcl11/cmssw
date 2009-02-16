import FWCore.ParameterSet.Config as cms
import copy

from TauAnalysis.CandidateTools.elecMuPairSelector_cfi import *

allElecMuPairs = cms.EDProducer("PATElecMuPairProducer",
  useLeadingTausOnly = cms.bool(False),
  srcLeg1 = cms.InputTag('selectedLayer1ElectronsTrkIPcumulative'),
  srcLeg2 = cms.InputTag('selectedLayer1MuonsTrkIPcumulative'),
  dRmin12 = cms.double(-1.),
  srcMET = cms.InputTag('allLayer1METs'),
  recoMode = cms.string(""),
  verbosity = cms.untracked.int32(0)
)

produceElecMuPairs = cms.Sequence( allElecMuPairs * selectElecMuPairs )


