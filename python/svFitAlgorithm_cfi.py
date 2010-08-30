import FWCore.ParameterSet.Config as cms

svFitLikelihoodDiTauKinematicsPhaseSpace = cms.PSet(
    pluginName = cms.string("svFitLikelihoodDiTauKinematicsPhaseSpace"),
    pluginType = cms.string("SVfitLikelihoodMuTauPairKinematics"),
    polarizationCoefficients = cms.PSet(
        LL = cms.double(0.5),
        LR = cms.double(0.5),
        RL = cms.double(0.5),
        RR = cms.double(0.5)
    ),
    leg = cms.PSet(
        pluginType = cms.string("SVfitMuonLikelihoodPhaseSpace")
    )
)

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

svFitLikelihoodDiTauPtBalance = cms.PSet(
    pluginName = cms.string("svFitLikelihoodDiTauPtBalance"),
    pluginType = cms.string("SVfitLikelihoodMuTauPairPtBalance")
)



