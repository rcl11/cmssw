import FWCore.ParameterSet.Config as cms

svFitLikelihoodDiTauKinematicsPhaseSpace = cms.PSet(
    pluginName = cms.string("svFitLikelihoodDiTauKinematicsPhaseSpace"),
    pluginType = cms.string("SVfitLikelihoodDiTauLegs"),
    # Always fit
    firstFitIteration = cms.uint32(0),
    leg1 = cms.PSet(
        pluginType = cms.string("SVfitLegLikelihoodPhaseSpace")
    ),
    leg2 = cms.PSet(
        pluginType = cms.string("SVfitLegLikelihoodPhaseSpace")
    )
)

svFitLikelihoodDiTauKinematicsPolarized = cms.PSet(
    pluginName = cms.string("svFitLikelihoodDiTauKinematicsPolarized"),
    pluginType = cms.string("SVfitLikelihoodDiTauLegs"),
    # Always fit
    firstFitIteration = cms.uint32(0),
    leg1 = cms.PSet(
        pluginType = cms.string("SVfitLegLikelihoodPolarizationBase"),
        usePolarization = cms.bool(True),
        useCollApproxFormulas = cms.bool(False)
        #useCollApproxFormulas = cms.bool(True)
    ),
    leg2 = cms.PSet(
        pluginType = cms.string("SVfitLegLikelihoodPolarizationBase"),
        ##mapRecToGenTauDecayModes = cms.PSet(
        ##    fileName = cms.string("/afs/cern.ch/user/v/veelken/public/plotsAHtoMuTau.root"),
        ##    meName = cms.string('DQMData/ahMuTauAnalyzer_woBtag/afterEvtSelNonCentralJetEt20bTag/TauQuantities/TauRecVsGenDecayMode')
        ##),
        decayModeParameters = cms.PSet(
            oneProngZeroPi0s = cms.PSet(
                pMin = cms.double(0.05)
            ),
            oneProngOnePi0 = cms.PSet(
                xSigma = cms.string("0.014"),
                xBias = cms.string("0.000"),
                pMin = cms.double(0.05)
            ),
            oneProngTwoPi0s = cms.PSet(
                xSigma = cms.string("0.013"),
                xBias = cms.string("0.000"),
                pMin = cms.double(0.05)
            ),
            threeProngZeroPi0s = cms.PSet(
                xSigma = cms.string("0.018"),
                xBias = cms.string("0.000"),
                pMin = cms.double(0.05)
            )
        ),
        usePolarization = cms.bool(True),
        useCollApproxFormulas = cms.bool(False)
        #useCollApproxFormulas = cms.bool(True)
    )
)

svFitLikelihoodMEt = cms.PSet(
    pluginName = cms.string("svFitLikelihoodMEt"),
    pluginType = cms.string("SVfitLikelihoodMEt"),
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

svFitLikelihoodTrackInfo = cms.PSet(
    pluginName = cms.string("svFitLikelihoodTrackInfo"),
    pluginType = cms.string("SVfitLikelihoodDiTauTrackInfo"),
    # Only fit on second iteration
    firstFitIteration = cms.uint32(1),
    leg1 = cms.PSet(
        pluginType = cms.string("SVfitLegLikelihoodTrackInfo"),
        useLinearApprox = cms.bool(True)
    ),
    leg2 = cms.PSet(
        pluginType = cms.string("SVfitLegLikelihoodTrackInfo"),
        useLinearApprox = cms.bool(True)
    ),

    useLifetimeConstraint = cms.bool(True)
)

svFitLikelihoodDiTauPtBalance = cms.PSet(
    pluginName = cms.string("svFitLikelihoodDiTauPtBalance"),
    pluginType = cms.string("SVfitLikelihoodDiTauPtBalance"),
    # Always fit
    firstFitIteration = cms.uint32(0),
)

svFitLikelihoodDiTauPtBalance2 = cms.PSet(
    pluginName = cms.string("svFitLikelihoodDiTauPtBalance2"),
    pluginType = cms.string("SVfitLikelihoodDiTauPtBalance2"),
    # Always fit
    firstFitIteration = cms.uint32(0),
    # Leptonic leg
    leg1 = cms.PSet(
        smear = cms.string("4.4 + 0.036*[0]"),
        gaussFrac = cms.string("0.93"),
        turnOnWidth = cms.string(
            "0.19 + (-0.0016*[0]) + (5.27e-6*[0]*[0]) + (-6.0e-9*[0]*[0]*[0])"),
        turnOnThreshold = cms.string("1.355 + 0.379*[0]"),
        gammaShape = cms.string("2"),
        gammaScale = cms.string("[0]/4"),
        overallNorm = cms.string("2"),
    ),
    # Tau leg
    leg2 = cms.PSet(
        smear = cms.string("6.3 + 0.019*[0]"),
        gaussFrac = cms.string("0.93"),
        turnOnWidth = cms.string(
            "0.23 + (-0.0022*[0]) + (7.91e-6*[0]*[0]) + (-9.4e-9*[0]*[0]*[0])"),
        turnOnThreshold = cms.string("2.2 + 0.365*[0]"),
        gammaShape = cms.string("2"),
        gammaScale = cms.string("[0]/4"),
        overallNorm = cms.string("2"),
    ),
)

svFitLikelihoodDiTauPt = cms.PSet(
    pluginName = cms.string("svFitLikelihoodDiTauPt"),
    pluginType = cms.string("SVfitLikelihoodDiTauPt"),
    # Always fit
    firstFitIteration = cms.uint32(0),
    pdf = cms.string("[0]*[1]/[2]*TMath::Sqrt(0.5*TMath::Pi())*TMath::Exp([1]*[1]/(2*[2]*[2]) - (x - [3])/[2])"
                    + "*(1. - TMath::Erf((([3] - x)/[1] + [1]/[2])/TMath::Sqrt(2.)))"),
    par0 = cms.string("5.38553e-1"),
    par1 = cms.string("1.39514e+0"),
    par2 = cms.string("1.08756e+1"),
    par3 = cms.string("8.85593e-1")
)

svFitLikelihoodDiTauProdZ0 = cms.PSet(
    pluginName = cms.string("svFitLikelihoodDiTauProdZ0"),
    pluginType = cms.string("SVfitLikelihoodDiTauProd"),
    # Always fit
    firstFitIteration = cms.uint32(0),
    process = cms.string("Z0"), # either 'Z0' or 'Higgs'
    pdfSet = cms.string("cteq65.LHgrid"),
    sqrtS = cms.double(7000.) # 7 TeV center-of-mass energy
)

