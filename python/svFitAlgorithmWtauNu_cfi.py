import FWCore.ParameterSet.Config as cms

svFitLikelihoodWtauNuKinematicsPhaseSpace = cms.PSet(
    pluginName = cms.string("svFitLikelihoodWtauNuKinematicsPhaseSpace"),
    pluginType = cms.string("SVfitLikelihoodTauNuPairKinematics"),
    # Always fit
    firstFitIteration = cms.uint32(0),
    leg1 = cms.PSet(
        pluginType = cms.string("SVfitTauLikelihoodPhaseSpace")
    )
)

svFitLikelihoodWtauNuMEt = cms.PSet(
    pluginName = cms.string("svFitLikelihoodWtauNuMEt"),
    pluginType = cms.string("SVfitLikelihoodTauNuPairMEt"),
    # Always fit
    firstFitIteration = cms.uint32(0),
    resolution = cms.PSet(
        parSigma = cms.string("7.54*(1 - 0.00542*x)"),
        parBias = cms.string("-0.96"),
        perpSigma = cms.string("6.85*(1 - 0.00547*x)"),
        perpBias = cms.string("0."),
    ),
    srcPFCandidates = cms.InputTag('particleFlow')
)

svFitLikelihoodWtauNuPtBalance = cms.PSet(
    pluginName = cms.string("svFitLikelihoodWtauNuPtBalance"),
    pluginType = cms.string("SVfitLikelihoodTauNuPairPtBalance"),
    # Always fit
    firstFitIteration = cms.uint32(0)
)
