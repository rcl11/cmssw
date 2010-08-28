import FWCore.ParameterSet.Config as cms

svFitLikelihoodMEt = cms.PSet(
    pluginName = cms.string("svFitLikelihoodMEt"),
    pluginType = cms.string("SVfitLikelihoodMEtMuTau"),
    resolution = cms.PSet(
        parSigma = cms.string("2.85 + 0.02072*x"),
        parBias = cms.string("1.183"),
        perpSigma = cms.string("2.3 + 0.02284*x"),
        perpBias = cms.string("0.0"),
    )
)
