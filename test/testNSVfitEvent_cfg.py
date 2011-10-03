import FWCore.ParameterSet.Config as cms

process = cms.Process("Test")

process.source = cms.Source("PoolSource",
  fileNames = cms.untracked.vstring(
    "file:patTuple.root"
  )
)
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(102) )

process.MessageLogger = cms.Service("MessageLogger")

process.TFileService = cms.Service("TFileService",
  fileName = cms.string('analyzeNSVfit.root')
)

## load nSVfitEventAnalyzer
process.load("TauAnalysis.CandidateTools.nSVfitEventAnalyzer_cfi")

## load nSVfitProducer (plugin version)
process.load("TauAnalysis.CandidateTools.nSVfitAlgorithmDiTau_cfi")

## Geometry and Detector Conditions (needed for a few patTuple production steps)
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
from Configuration.PyReleaseValidation.autoCond import autoCond
process.GlobalTag.globaltag = cms.string( autoCond[ 'startup' ] )
process.load("Configuration.StandardSequences.MagneticField_cff")

process.nSVfitProducerByLikelihoodMaximization.config.event.resonances.A.daughters.leg1.src = "cleanPatMuons"
process.nSVfitProducerByLikelihoodMaximization.config.event.resonances.A.daughters.leg2.src = "cleanPatElectrons"
process.nSVfitProducerByLikelihoodMaximization.config.event.resonances.A.daughters.leg2.likelihoodFunctions = cms.VPSet(process.nSVfitElectronLikelihoodPhaseSpace)
process.nSVfitProducerByLikelihoodMaximization.config.event.resonances.A.daughters.leg2.builder = process.nSVfitTauToElecBuilder
process.nSVfitProducerByLikelihoodMaximization.config.event.srcMEt = "patMETs"
process.nSVfitProducerByLikelihoodMaximization.algorithm.verbosity = 0

## check the event content
process.content = cms.EDAnalyzer("EventContentAnalyzer")

process.p = cms.Path(
    process.nSVfitProducerByLikelihoodMaximization
# * process.content
  * process.nSVfitEventAnalyzer
    )
