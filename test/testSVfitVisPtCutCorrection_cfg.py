import FWCore.ParameterSet.Config as cms

import copy
import os
import re

process = cms.Process("produceSVfitNtuple")

process.load('FWCore/MessageService/MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 100
process.MessageLogger.cerr.threshold = cms.untracked.string('INFO')
process.load('Configuration.StandardSequences.Geometry_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = cms.string('START42_V13::All')

#--------------------------------------------------------------------------------
# define configuration parameter default values

sample = 'Ztautau'
#sample = 'vbfHiggs120'
#sample = 'vbfHiggs130'
sample_type = None
#channel = 'eleTau'
channel = 'muTau'
#channel = 'eleMu'
#channel = 'diTau'
metResolution = None # take reconstructed PFMET
#metResolution = 5. # produce "toy" MET = generated MET plus 5 GeV Gaussian smearing in x/y direction
inputFileNames = None
maxEvents = 20
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
# define "hooks" for replacing configuration parameters
# in case running jobs on the CERN batch system/grid
#
#__sample = '#sample#'
##__sample_type = '#sample_type#'
##__channel = '#channel#'
#__maxEvents = #maxEvents#
##__inputFileNames = #inputFileNames#
##__outputFileName = '#outputFileName#'
#
#--------------------------------------------------------------------------------

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(maxEvents)
)

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring()
)

#--------------------------------------------------------------------------------
# read AOD/RECO input files from local directory/castor

if inputFileNames is None:    
    import TauAnalysis.Configuration.tools.castor as castor
    inputFilePath = None
    inputFile_regex = None
    if sample == 'Ztautau':
        inputFilePath = '/data1/veelken/CMSSW_4_2_x/skims/user/v/veelken/CMSSW_4_2_x/skims/ZtoMuTau/GenZtoMuTauWithinAcc/2012Jan11/'
        inputFile_regex = \
          r"[a-zA-Z0-9_/:.]*skimGenZtoMuTauWithinAcc_Ztautau_2012Jan11_(?P<gridJob>\d*)(_(?P<gridTry>\d*))*_(?P<hash>[a-zA-Z0-9]*).root"
        sample_type = 'Z'
    elif sample == 'vbfHiggs120':
        inputFilePath = '/data1/veelken/CMSSW_4_2_x/skims/user/v/veelken/CMSSW_4_2_x/skims/GenAHtoMuTauWithinAcc/vbfHiggs120/'
        inputFile_regex = \
          r"[a-zA-Z0-9_/:.]*skimGenAHtoMuTauWithinAcc_vbfHiggs120_2012Jan30_(?P<gridJob>\d*)(_(?P<gridTry>\d*))*_(?P<hash>[a-zA-Z0-9]*).root"
        sample_type = 'Higgs'
    elif sample == 'vbfHiggs130':
        inputFilePath = '/data1/veelken/CMSSW_4_2_x/skims/user/v/veelken/CMSSW_4_2_x/skims/GenAHtoMuTauWithinAcc/vbfHiggs130/'
        inputFile_regex = \
          r"[a-zA-Z0-9_/:.]*skimGenAHtoMuTauWithinAcc_vbfHiggs130_2012Jan11_(?P<gridJob>\d*)(_(?P<gridTry>\d*))*_(?P<hash>[a-zA-Z0-9]*).root"
        sample_type = 'Higgs'
    else:
        raise ValueError("Invalid sample = %s !!" % sample_type)

    inputFileNames = []
    if inputFilePath.find('/castor/') != -1:
        inputFileNames = [ 'rfio:%s' % file_info['path'] for file_info in castor.nslsl(inputFilePath) ]
    else:
        inputFileNames = [ 'file:%s' % os.path.join(inputFilePath, file_name) for file_name in os.listdir(inputFilePath) ]

    inputFile_matcher = re.compile(inputFile_regex)

    inputFileNames_matched = []
    for inputFileName in inputFileNames:
        if inputFile_matcher.match(inputFileName):
   	    inputFileNames_matched.append(inputFileName)

    #print "inputFileNames_matched = %s" % inputFileNames_matched

    setattr(process.source, "fileNames", cms.untracked.vstring(inputFileNames_matched))

    ##setattr(process.source, "eventsToProcess", cms.untracked.VEventRange(
    ##    '1:74:73837491',
    ##    '1:74:73389858'
    ##))
else:
    setattr(process.source, "fileNames", cms.untracked.vstring(inputFileNames))
#--------------------------------------------------------------------------------

process.testSVfitVisPtCutCorrSequence = cms.Sequence()

#--------------------------------------------------------------------------------
# select collections of electrons, muons and tau-jets
# matching genuine tau -> e, tau -> mu and tau -> hadronic decays on generator level

process.load("RecoTauTag/Configuration/RecoPFTauTag_cff")
process.testSVfitVisPtCutCorrSequence += process.PFTau

process.load("PhysicsTools/JetMCAlgos/TauGenJets_cfi")
process.testSVfitVisPtCutCorrSequence += process.tauGenJets

genElectronsFromTauDecays = None
genMuonsFromTauDecays = None
genTauJetsFromTauDecays = None
genTauPairs = None
genTaus = None
if sample_type == 'Z':
    process.load("TauAnalysis/GenSimTools/gen_decaysFromZs_cfi")
    process.testSVfitVisPtCutCorrSequence += process.produceGenDecayProductsFromZs
    genElectronsFromTauDecays = 'genElectronsFromZtautauDecays'
    genMuonsFromTauDecays = 'genMuonsFromZtautauDecays'
    genTauJetsFromTauDecays = 'genHadronsFromZtautauDecays'
    genTauPairs = 'genZdecayToTaus'
    genTaus = 'genTausFromZs'
elif sample_type == 'Higgs':
    process.load("TauAnalysis/GenSimTools/gen_decaysFromAHs_cfi")
    process.testSVfitVisPtCutCorrSequence += process.produceGenDecayProductsFromAHs
    genElectronsFromTauDecays = 'genElectronsFromAHtautauDecays'
    genMuonsFromTauDecays = 'genMuonsFromAHtautauDecays'
    genTauJetsFromTauDecays = 'genHadronsFromAHtautauDecays'
    genTauPairs = 'genAHdecayToTaus'
    genTaus = 'genTausFromAHs'
else:
    raise ValueError("Invalid sample type = %s !!" % sample_type)

process.genMatchedElectrons = cms.EDFilter("GsfElectronAntiOverlapSelector",
    src = cms.InputTag('gsfElectrons'),
    srcNotToBeFiltered = cms.VInputTag(genElectronsFromTauDecays),
    dRmin = cms.double(0.3),
    invert = cms.bool(True),
    filter = cms.bool(False)                                                          
)                                                  
process.testSVfitVisPtCutCorrSequence += process.genMatchedElectrons

process.genMatchedMuons = cms.EDFilter("MuonAntiOverlapSelector",
    src = cms.InputTag('muons'),
    srcNotToBeFiltered = cms.VInputTag(genMuonsFromTauDecays),
    dRmin = cms.double(0.3),
    invert = cms.bool(True),
    filter = cms.bool(False)                                                          
)                        
process.testSVfitVisPtCutCorrSequence += process.genMatchedMuons

process.selectedTauJets = cms.EDFilter("PFTauSelector",
    src = cms.InputTag('hpsPFTauProducer'),
    discriminators = cms.VPSet(
        cms.PSet(
            discriminator = cms.InputTag('hpsPFTauDiscriminationByDecayModeFinding'),
            selectionCut = cms.double(0.5)
        ),
        cms.PSet(
            discriminator = cms.InputTag('hpsPFTauDiscriminationByLooseCombinedIsolationDBSumPtCorr'),
            selectionCut = cms.double(0.5)
        )                        
    ),
    cut = cms.string("pt > 20. & abs(eta) < 2.3")                        
)
process.testSVfitVisPtCutCorrSequence += process.selectedTauJets

process.genMatchedTauJets = cms.EDFilter("PFTauAntiOverlapSelector",
    src = cms.InputTag('selectedTauJets'),
    srcNotToBeFiltered = cms.VInputTag(genTauJetsFromTauDecays),
    dRmin = cms.double(0.3),
    invert = cms.bool(True),
    filter = cms.bool(False)                                                          
)
process.testSVfitVisPtCutCorrSequence += process.genMatchedTauJets
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
# require event to contain generator level tau lepton pair,
# decaying in the sprecified channel
numElectrons = None
numMuons     = None
numTauJets   = None
if channel == 'muTau':
    numElectrons = 0
    numMuons     = 1
    numTauJets   = 1
elif channel == 'eleTau':
    numElectrons = 1
    numMuons     = 0
    numTauJets   = 1
elif channel == 'eleMu':
    numElectrons = 1
    numMuons     = 1
    numTauJets   = 0
elif channel == 'diTau':
    numElectrons = 0
    numMuons     = 0
    numTauJets   = 2  
else:
    raise ValueError("Invalid channel = %s !!" % channel)

if numElectrons > 0:
    process.electronFilter = cms.EDFilter("CandViewCountFilter",
        src = cms.InputTag(genElectronsFromTauDecays),
        minNumber = cms.uint32(1),
        maxNumber = cms.uint32(1)                      
    ) 
    process.testSVfitVisPtCutCorrSequence += process.electronFilter

if numMuons > 0:    
    process.muonFilter = cms.EDFilter("CandViewCountFilter",
        src = cms.InputTag(genMuonsFromTauDecays),
        minNumber = cms.uint32(1),
        maxNumber = cms.uint32(1)                      
    )
    process.testSVfitVisPtCutCorrSequence += process.muonFilter

if numTauJets > 0:
    process.tauFilter = cms.EDFilter("CandViewCountFilter",
        src = cms.InputTag(genTauJetsFromTauDecays),
        minNumber = cms.uint32(1),
        maxNumber = cms.uint32(1)                      
    )
    process.testSVfitVisPtCutCorrSequence += process.tauFilter
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
# produce PAT-tuple
process.load("PhysicsTools/PatAlgos/patSequences_cff")

# switch to ak5PFJets
from PhysicsTools.PatAlgos.tools.jetTools import *
switchJetCollection(
    process,
    cms.InputTag('ak5PFJets'),
    doJTA = True,
    doBTagging = False,
    jetCorrLabel = ( 'AK5PF', cms.vstring([ 'L1FastJet', 'L2Relative', 'L3Absolute' ]) ),
    doType1MET = False,
    doJetID = True,
    jetIdLabel = "ak5",
    outputModule = ''
)

# switch to PFMET and apply Type 1 MET corrections
process.load("PhysicsTools/PatUtils/patPFMETCorrections_cff")
process.testSVfitVisPtCutCorrSequence += process.kt6PFJets
process.testSVfitVisPtCutCorrSequence += process.ak5PFJets
process.patDefaultSequence.remove(process.patMETs)
process.patDefaultSequence.remove(process.patMETCorrections)
process.testSVfitVisPtCutCorrSequence += process.patDefaultSequence
process.testSVfitVisPtCutCorrSequence += process.producePatPFMETCorrections

# switch to HPS PFTaus (and disable all "cleaning" cuts)
from PhysicsTools.PatAlgos.tools.tauTools import *
switchToPFTauHPS(process)

# switch input for pat::Electrons, pat::Muons and pat::Taus to gen. matched objects
process.patElectrons.electronSource = cms.InputTag('genMatchedElectrons')
process.patMuons.muonSource = cms.InputTag('genMatchedMuons')
process.patTaus.tauSource = cms.InputTag('genMatchedTauJets')

# disable matching pat::Electrons, pat::Muons and pat::Taus to generator level quantities
import PhysicsTools.PatAlgos.tools.helpers as patutils
removeMCMatching(process, ["All"], outputInProcess = False)
process.patDefaultSequence.remove(process.patJetPartonMatch)

# disable all pat::Electron embedding
for objSelAttrName in dir(process.patElectrons):
    objSelAttr = getattr(process.patElectrons, objSelAttrName)
    if isinstance(objSelAttr, cms.bool) and (objSelAttrName.startswith("embed") or objSelAttrName.startswith("add")):        
        setattr(process.patElectrons, objSelAttrName, cms.bool(False))
        
# disable all pat::Muon embedding
for objSelAttrName in dir(process.patMuons):
    objSelAttr = getattr(process.patMuons, objSelAttrName)
    if isinstance(objSelAttr, cms.bool) and (objSelAttrName.startswith("embed") or objSelAttrName.startswith("add")):
        setattr(process.patMuons, objSelAttrName, cms.bool(False))

# disable all pat::Tau embedding
for objSelAttrName in dir(process.patTaus):
    objSelAttr = getattr(process.patTaus, objSelAttrName)
    if isinstance(objSelAttr, cms.bool) and (objSelAttrName.startswith("embed") or objSelAttrName.startswith("add")):
        setattr(process.patTaus, objSelAttrName, cms.bool(False))
process.patTaus.isoDeposits = cms.PSet()
process.patTaus.tauJetCorrFactorsSource = cms.VInputTag()
process.patTaus.tauIDSources = cms.PSet()
process.patTaus.userIsolation = cms.PSet()
process.cleanPatTaus.preselection = cms.string('')
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------    
# apply vertex selection
process.load("TauAnalysis.RecoTools.recoVertexSelection_cff")
process.testSVfitVisPtCutCorrSequence += process.selectPrimaryVertex
#--------------------------------------------------------------------------------  

#--------------------------------------------------------------------------------
# produce collections of jets not overlapping with leptons
# and with (corrected) Pt > 20 GeV
process.selectedPatJetsNotOverlappingWithLeptons = cms.EDFilter("PATJetAntiOverlapSelector",
    src = cms.InputTag('patJets'),
    # CV: set 'srcNotToBeFiltered' to collections of electrons, muons and taus passing analysis specific selection criteria
    srcNotToBeFiltered = cms.VInputTag(
        'genMatchedElectrons',
        'genMatchedMuons',
        'genMatchedTauJets'
    ),
    dRmin = cms.double(0.5),
    filter = cms.bool(False)      
)
process.testSVfitVisPtCutCorrSequence += process.selectedPatJetsNotOverlappingWithLeptons

process.selectedPatJetsPt20 = cms.EDFilter("PATJetSelector",
    src = cms.InputTag('selectedPatJetsNotOverlappingWithLeptons'),                                   
    cut = cms.string('pt > 20.'), 
    filter = cms.bool(False)
)
process.testSVfitVisPtCutCorrSequence += process.selectedPatJetsPt20
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
# select ak5PFJets corresponding to selected pat::Jets
process.ak5PFJetsNotOverlappingWithLeptonsCorrPtGt20 = cms.EDFilter("PFJetAntiOverlapSelector",
    src = cms.InputTag('ak5PFJets'),
    srcNotToBeFiltered = cms.VInputTag('selectedPatJetsPt20'),
    dRmin = cms.double(1.e-1),
    invert = cms.bool(True),
    filter = cms.bool(False)                                                          
)
process.testSVfitVisPtCutCorrSequence += process.ak5PFJetsNotOverlappingWithLeptonsCorrPtGt20
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
# select PFCandidates ("unclustered energy") not within jets
from CommonTools.ParticleFlow.TopProjectors.pfNoJet_cfi import pfNoJet
process.pfCandsNotInSelectedJets = pfNoJet.clone(
    topCollection = cms.InputTag('ak5PFJetsNotOverlappingWithLeptonsCorrPtGt20'),
    bottomCollection = cms.InputTag('particleFlow')
)
process.testSVfitVisPtCutCorrSequence += process.pfCandsNotInSelectedJets
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
# produce PFMET significance cov. matrix
from RecoMET.METProducers.METSigParams_cfi import *
process.pfMEtSignCovMatrix = cms.EDProducer("PFMEtSignCovMatrixProducer",
    METSignificance_params,                     
    src = cms.VInputTag(
        'genMatchedElectrons',
        'genMatchedMuons',
        'genMatchedTauJets',                                   
        'ak5PFJetsNotOverlappingWithLeptonsCorrPtGt20',
        'pfCandsNotInSelectedJets'
    ),
    addJERcorr = cms.PSet(
        inputFileName = cms.FileInPath('PhysicsTools/PatUtils/data/pfJetResolutionMCtoDataCorrLUT.root'),
        lutName = cms.string('pfJetResolutionMCtoDataCorrLUT')
    )
)
process.testSVfitVisPtCutCorrSequence += process.pfMEtSignCovMatrix
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
# produce genMET
process.load("RecoMET.Configuration.GenMETParticles_cff")
process.testSVfitVisPtCutCorrSequence += process.genParticlesForMETAllVisible

process.load("RecoMET.METProducers.genMetTrue_cfi")
process.genMetFromGenParticles = process.genMetTrue.clone(
    src = cms.InputTag('genParticlesForMETAllVisible'),
    alias = cms.string('genMetFromGenParticles')
)
process.testSVfitVisPtCutCorrSequence += process.genMetFromGenParticles
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
# run SVfit with different options, make plots

process.load("TauAnalysis.CandidateTools.nSVfitAlgorithmDiTau_cfi")
from TauAnalysis.CandidateTools.nSVfitAlgorithmVisPtCutCorrections_cfi import *

srcGenLeg1                        = None
srcRecLeg1                        = None
nSVfitLeg1LikelihoodPhaseSpace    = None
nSVfitLeg1LikelihoodMatrixElement = None
nSVfitLeg1visPtCutThreshold       = None
srcGenLeg2                        = None
srcRecLeg2                        = None
nSVfitLeg2LikelihoodPhaseSpace    = None
nSVfitLeg2LikelihoodMatrixElement = None
nSVfitLeg2visPtCutThreshold       = None
if channel == 'muTau':
    srcGenLeg1                        = genMuonsFromTauDecays
    srcRecLeg1                        = 'patMuons'
    nSVfitLeg1LikelihoodPhaseSpace    = process.nSVfitMuonLikelihoodPhaseSpace
    nSVfitLeg1LikelihoodMatrixElement = process.nSVfitMuonLikelihoodMatrixElement
    nSVfitLeg1visPtCutThreshold       = 15.
    srcGenLeg2                        = genTauJetsFromTauDecays
    srcRecLeg2                        = 'patTaus'
    nSVfitLeg2LikelihoodPhaseSpace    = process.nSVfitTauLikelihoodPhaseSpace
    nSVfitLeg2LikelihoodMatrixElement = process.nSVfitTauLikelihoodMatrixElement
    nSVfitLeg2visPtCutThreshold       = 20.
elif channel == 'eleTau':
    srcGenLeg1                        = genElectronsFromTauDecays
    srcRecLeg1                        = 'patElectrons'
    nSVfitLeg1LikelihoodPhaseSpace    = process.nSVfitElectronLikelihoodPhaseSpace
    nSVfitLeg1LikelihoodMatrixElement = process.nSVfitElectronLikelihoodMatrixElement
    nSVfitLeg1visPtCutThreshold       = 20.
    srcGenLeg2                        = genTauJetsFromTauDecays
    srcRecLeg2                        = 'patTaus'
    nSVfitLeg2LikelihoodPhaseSpace    = process.nSVfitTauLikelihoodPhaseSpace
    nSVfitLeg2LikelihoodMatrixElement = process.nSVfitTauLikelihoodMatrixElement
    nSVfitLeg2visPtCutThreshold       = 20.
elif channel == 'eleMu':
    srcGenLeg1                        = genElectronsFromTauDecays
    srcRecLeg1                        = 'patElectrons'
    nSVfitLeg1LikelihoodPhaseSpace    = process.nSVfitElectronLikelihoodPhaseSpace
    nSVfitLeg1LikelihoodMatrixElement = process.nSVfitElectronLikelihoodMatrixElement
    nSVfitLeg1visPtCutThreshold       = 15.
    srcGenLeg2                        = genMuonsFromTauDecays
    srcRecLeg2                        = 'patMuons'
    nSVfitLeg2LikelihoodPhaseSpace    = process.nSVfitMuonLikelihoodPhaseSpace 
    nSVfitLeg2LikelihoodMatrixElement = process.nSVfitMuonLikelihoodMatrixElement
    nSVfitLeg2visPtCutThreshold       = 15.
elif channel == 'diTau':
    srcGenLeg1                        = genTauJetsFromTauDecays
    srcRecLeg1                        = 'patTaus'
    nSVfitLeg1LikelihoodPhaseSpace    = process.nSVfitTauLikelihoodPhaseSpace
    nSVfitLeg1LikelihoodMatrixElement = process.nSVfitTauLikelihoodMatrixElement
    nSVfitLeg1visPtCutThreshold       = 35.
    srcGenLeg2                        = genTauJetsFromTauDecays
    srcRecLeg2                        = 'patTaus'
    nSVfitLeg2LikelihoodPhaseSpace    = process.nSVfitTauLikelihoodPhaseSpace
    nSVfitLeg2LikelihoodMatrixElement = process.nSVfitTauLikelihoodMatrixElement
    nSVfitLeg2visPtCutThreshold       = 35.
else:
    raise ValueError("Invalid channel = %s !!" % channel)

srcRecMEt           = None
srcRecMEtCovMatrix  = None
if metResolution is not None:
    process.toyMEt = cms.EDProducer("ToyPATMEtProducer",
        srcGenMEt = cms.InputTag('genMetFromGenParticles'),
        resolutionX = cms.double(metResolution),
        resolutionY = cms.double(metResolution)
    )
    process.testSVfitVisPtCutCorrSequence += process.toyMEt
    srcRecMEt = 'toyMEt'

    process.toyMEtCovMatrix = cms.EDProducer("ToyMEtSignCovMatrixProducer",
        resolutionX = cms.double(metResolution),
        resolutionY = cms.double(metResolution)
    )
    process.testSVfitVisPtCutCorrSequence += process.toyMEtCovMatrix
    srcRecMEtCovMatrix = 'toyMEtCovMatrix'
else:
    srcRecMEt = 'patType1CorrectedPFMet'

for idxSVfitOption in range(14):
    ##if not (idxSVfitOption == 3):
    ##    continue
    nSVfitProducer = None
    if idxSVfitOption == 3 or idxSVfitOption == 6 or idxSVfitOption == 11 or idxSVfitOption == 13:
        nSVfitProducer = copy.deepcopy(process.nSVfitProducerByIntegration)
    else:
        nSVfitProducer = copy.deepcopy(process.nSVfitProducerByLikelihoodMaximization)
    nSVfitProducer.config.event.resonances.A.daughters.leg1.src = cms.InputTag(srcRecLeg1)
    nSVfitProducer.config.event.resonances.A.daughters.leg2.src = cms.InputTag(srcRecLeg2)
    nSVfitProducer.config.event.resonances.A.builder.polStates = cms.vstring('undefined')
    if idxSVfitOption <= 1:
        nSVfitProducer.config.event.resonances.A.daughters.leg1.likelihoodFunctions = cms.VPSet(nSVfitLeg1LikelihoodPhaseSpace.clone())
        nSVfitProducer.config.event.resonances.A.daughters.leg2.likelihoodFunctions = cms.VPSet(nSVfitLeg2LikelihoodPhaseSpace.clone())
    elif idxSVfitOption <= 12:
        nSVfitProducer.config.event.resonances.A.daughters.leg1.likelihoodFunctions = cms.VPSet(nSVfitLeg1LikelihoodMatrixElement.clone())
        nSVfitProducer.config.event.resonances.A.daughters.leg2.likelihoodFunctions = cms.VPSet(nSVfitLeg2LikelihoodMatrixElement.clone())
    else:
        nSVfitProducer.config.event.resonances.A.daughters.leg1.likelihoodFunctions = cms.VPSet(
            nSVfitLeg1LikelihoodMatrixElement.clone(),
            cms.PSet(
                pluginName = cms.string("nSVfitLeg1LikelihoodTrackInfo"),
                pluginType = cms.string("NSVfitTauDecayLikelihoodTrackInfo"),
                useLifetimeConstraint = cms.bool(True),
                verbosity = cms.int32(0)  
            )
        )
        nSVfitProducer.config.event.resonances.A.daughters.leg2.likelihoodFunctions = cms.VPSet(
            nSVfitLeg2LikelihoodMatrixElement.clone(),
            cms.PSet(
                pluginName = cms.string("nSVfitLeg2LikelihoodTrackInfo"),
                pluginType = cms.string("NSVfitTauDecayLikelihoodTrackInfo"),
                useLifetimeConstraint = cms.bool(True),
                verbosity = cms.int32(0)  
            )
        )
    if idxSVfitOption >= 1:
        nSVfitProducer.config.event.resonances.A.likelihoodFunctions = cms.VPSet()        
        nSVfitProducer.config.event.resonances.A.daughters.leg1.likelihoodFunctions[0].applySinThetaFactor = cms.bool(False)
        nSVfitProducer.config.event.resonances.A.daughters.leg1.likelihoodFunctions[0].applyVisPtCutCorrection = cms.bool(False)
        nSVfitProducer.config.event.resonances.A.daughters.leg2.likelihoodFunctions[0].applySinThetaFactor = cms.bool(False)
        nSVfitProducer.config.event.resonances.A.daughters.leg2.likelihoodFunctions[0].applyVisPtCutCorrection = cms.bool(False)
    if idxSVfitOption == 10 or idxSVfitOption == 11:
        nSVfitProducer.config.event.resonances.A.daughters.leg1.likelihoodFunctions[0].applyVisPtCutCorrection = cms.bool(True)
        nSVfitProducer.config.event.resonances.A.daughters.leg1.likelihoodFunctions[0].visPtCutThreshold = \
          cms.double(nSVfitLeg2visPtCutThreshold)
        nSVfitProducer.config.event.resonances.A.daughters.leg2.likelihoodFunctions[0].applyVisPtCutCorrection = cms.bool(True)
        nSVfitProducer.config.event.resonances.A.daughters.leg2.likelihoodFunctions[0].visPtCutThreshold = \
          cms.double(nSVfitLeg2visPtCutThreshold)
    if idxSVfitOption == 4:
        nSVfitResonanceLikelihoodPrior_strength1 = copy.deepcopy(process.nSVfitResonanceLikelihoodPrior)
        nSVfitResonanceLikelihoodPrior_strength1.parameter.par0 = cms.double(1.)
        nSVfitProducer.config.event.resonances.A.likelihoodFunctions = cms.VPSet(nSVfitResonanceLikelihoodPrior_strength1)
    elif idxSVfitOption == 5 or idxSVfitOption == 6:
        nSVfitResonanceLikelihoodPrior_strength2 = copy.deepcopy(process.nSVfitResonanceLikelihoodPrior)
        nSVfitResonanceLikelihoodPrior_strength2.parameter.par0 = cms.double(2.)
        nSVfitProducer.config.event.resonances.A.likelihoodFunctions = cms.VPSet(nSVfitResonanceLikelihoodPrior_strength2)
    elif idxSVfitOption == 7:
        nSVfitResonanceLikelihoodPrior_strength3 = copy.deepcopy(process.nSVfitResonanceLikelihoodPrior)
        nSVfitResonanceLikelihoodPrior_strength3.parameter.par0 = cms.double(3.)
        nSVfitProducer.config.event.resonances.A.likelihoodFunctions = cms.VPSet(nSVfitResonanceLikelihoodPrior_strength3)
    elif idxSVfitOption == 8:
        nSVfitResonanceLikelihoodPrior_strength4 = copy.deepcopy(process.nSVfitResonanceLikelihoodPrior)
        nSVfitResonanceLikelihoodPrior_strength4.parameter.par0 = cms.double(4.)
        nSVfitProducer.config.event.resonances.A.likelihoodFunctions = cms.VPSet(nSVfitResonanceLikelihoodPrior_strength4)
    elif idxSVfitOption == 9:
        nSVfitResonanceLikelihoodPrior_strength5 = copy.deepcopy(process.nSVfitResonanceLikelihoodPrior)
        nSVfitResonanceLikelihoodPrior_strength5.parameter.par0 = cms.double(5.)
        nSVfitProducer.config.event.resonances.A.likelihoodFunctions = cms.VPSet(nSVfitResonanceLikelihoodPrior_strength5)
    nSVfitProducer.config.event.srcMEt = cms.InputTag(srcRecMEt)
    if srcRecMEtCovMatrix is not None:
        nSVfitProducer.config.event.likelihoodFunctions[0].srcMEtCovMatrix = cms.InputTag(srcRecMEtCovMatrix)
    nSVfitProducer.config.event.srcPrimaryVertex = cms.InputTag('selectedPrimaryVertexPosition')
    nSVfitProducerName = "nSVfitProducer%i" % idxSVfitOption
    setattr(process, nSVfitProducerName, nSVfitProducer)
    process.testSVfitVisPtCutCorrSequence += nSVfitProducer

    nSVfitAnalyzerType = None
    if idxSVfitOption == 3 or idxSVfitOption == 6 or idxSVfitOption == 11 or idxSVfitOption == 13:
        nSVfitAnalyzerType = "NSVfitEventHypothesisByIntegrationAnalyzer"
    else:  
        nSVfitAnalyzerType = "NSVfitEventHypothesisAnalyzer"
    nSVfitAnalyzer = cms.EDAnalyzer(nSVfitAnalyzerType,
        srcEventHypotheses = cms.InputTag(nSVfitProducerName),
        srcGenLeg1 = cms.InputTag(srcGenLeg1),
        srcGenLeg2 = cms.InputTag(srcGenLeg2),
        srcGenMEt = cms.InputTag('genMetFromGenParticles'),
        srcPFMEtCovMatrix = cms.InputTag('pfMEtSignCovMatrix'),
        srcWeights = cms.VInputTag(),                
        dqmDirectory = cms.string("nSVfitAnalyzerOption%i" % idxSVfitOption)
    )                                    
    nSVfitAnalyzerName = "nSVfitAnalyzer%i" % idxSVfitOption
    setattr(process, nSVfitAnalyzerName, nSVfitAnalyzer)
    process.testSVfitVisPtCutCorrSequence += nSVfitAnalyzer
#--------------------------------------------------------------------------------

process.DQMStore = cms.Service("DQMStore")

process.saveSVfitVisPtCutCorrectionPlots = cms.EDAnalyzer("DQMSimpleFileSaver",
    outputFileName = cms.string('testSVfitVisPtCutCorrection_%s_2012Mar20.root' % sample)
)

process.p = cms.Path(
    process.testSVfitVisPtCutCorrSequence
   + process.saveSVfitVisPtCutCorrectionPlots
)

processDumpFile = open('testSVfitVisPtCutCorrection.dump', 'w')
print >> processDumpFile, process.dumpPython()




