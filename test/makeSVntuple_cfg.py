import FWCore.ParameterSet.Config as cms
import copy

process = cms.Process('runAHtoMuTau')

# import of standard configurations for RECOnstruction
# of electrons, muons and tau-jets with non-standard isolation cones
process.load('Configuration/StandardSequences/Services_cff')
process.load('FWCore/MessageService/MessageLogger_cfi')
#process.MessageLogger.cerr.FwkReport.reportEvery = 100
process.MessageLogger.cerr.INFO.limit = cms.untracked.int32(100000)
process.MessageLogger.cerr.INFO.limit = cms.untracked.int32(0)
process.load('Configuration/StandardSequences/GeometryIdeal_cff')
process.load('Configuration/StandardSequences/MagneticField_cff')
process.load('Configuration/StandardSequences/Reconstruction_cff')
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = cms.string('MC_36Y_V7A::All')

# import particle data table
# needed for print-out of generator level information
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")

#--------------------------------------------------------------------------------
# import sequences for PAT-tuple production
process.load("TauAnalysis.Configuration.producePatTuple_cff")
process.load("TauAnalysis.Configuration.producePatTupleAHtoMuTauSpecific_cff")

# import sequence for event selection
process.load("TauAnalysis.Configuration.selectAHtoMuTau_cff")

# import sequence for filling of histograms, cut-flow table
# and of run + event number pairs for events passing event selection
process.load("TauAnalysis.Configuration.analyzeAHtoMuTau_cff")

# import configuration parameters for submission of jobs to CERN batch system
# (running over skimmed samples stored on CASTOR)
from TauAnalysis.Configuration.recoSampleDefinitionsAHtoMuTau_7TeV_cfi import *
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
# print memory consumed by cmsRun
# (for debugging memory leaks)
#process.SimpleMemoryCheck = cms.Service("SimpleMemoryCheck",
#    ignoreTotal = cms.untracked.int32(1) # default is one
#)

process.printGenParticleList = cms.EDAnalyzer("ParticleListDrawer",
    src = cms.InputTag("genParticles"),
    maxEventsToPrint = cms.untracked.int32(100)
)

# print event content 
process.printEventContent = cms.EDAnalyzer("EventContentAnalyzer")

# print debug information whenever plugins get loaded dynamically from libraries
# (for debugging problems with plugin related dynamic library loading)
#process.add_( cms.Service("PrintLoadingPlugins") )
#--------------------------------------------------------------------------------

process.DQMStore = cms.Service("DQMStore")

process.saveAHtoMuTauPlots = cms.EDAnalyzer("DQMSimpleFileSaver",
    outputFileName = cms.string('plotsAHtoMuTau.root')
)

process.maxEvents = cms.untracked.PSet(            
    input = cms.untracked.int32(-1)
)

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        #'/store/relval/CMSSW_3_1_2/RelValZTT/GEN-SIM-RECO/STARTUP31X_V2-v1/0007/A4DD1FAE-B178-DE11-B608-001D09F24EAC.root',
        #'/store/relval/CMSSW_3_1_2/RelValZTT/GEN-SIM-RECO/STARTUP31X_V2-v1/0007/9408B54D-CB78-DE11-9AEB-001D09F2503C.root'
        #'rfio:/castor/cern.ch/user/l/lusito/SkimOctober09/ZtautauSkimMT314_3/muTauSkim_1.root',
        #'rfio:/castor/cern.ch/user/l/lusito/SkimOctober09/ZtautauSkimMT314_3/muTauSkim_2.root'
        'file:/mnt/hadoop/store/user/cms1227/Ztautau/friis/Ztautau/ZtoMuTauSkimTest4/60dace99b1523ca69748e499060d8ab3/muTauSkim_1.root'
    )
    #skipBadFiles = cms.untracked.bool(True) 
)

#--------------------------------------------------------------------------------
# define "hooks" for replacing configuration parameters
# in case running jobs on the CERN batch system
#
#__process.source.fileNames = #inputFileNames#
#__process.maxEvents.input = cms.untracked.int32(#maxEvents#)
#__process.analyzeAHtoMuTauEvents_woBtag.filters[0] = copy.deepcopy(#genPhaseSpaceCut#)
#__process.analyzeAHtoMuTauEvents_wBtag.filters[0] = copy.deepcopy(#genPhaseSpaceCut#)
#__process.saveAHtoMuTauPlots.outputFileName = #plotsOutputFileName#
#__#isBatchMode#
#__#disableEventDump#
#
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
# import utility function for switching pat::Tau input
# to different reco::Tau collection stored on AOD
from PhysicsTools.PatAlgos.tools.tauTools import * 

# comment-out to take reco::CaloTaus instead of reco::PFTaus
# as input for pat::Tau production
#switchToCaloTau(process)

# comment-out to take shrinking dR = 5.0/Et(PFTau) signal cone
# instead of fixed dR = 0.07 signal cone reco::PFTaus
# as input for pat::Tau production
switchToPFTauShrinkingCone(process)
#switchToPFTauFixedCone(process)
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
# import utility function for managing pat::Jets
from PhysicsTools.PatAlgos.tools.jetTools import *

# uncomment to replace caloJets by pfJets
switchJetCollection(process, jetCollection = cms.InputTag("iterativeCone5PFJets"))
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
# import utility function for managing pat::METs
from TauAnalysis.Configuration.tools.metTools import *

# uncomment to add pfMET
# set Boolean swich to true in order to apply type-1 corrections
addPFMet(process, correct = False)

# uncomment to replace caloMET by pfMET in all di-tau objects
process.load("TauAnalysis.CandidateTools.diTauPairProductionAllKinds_cff")
replaceMETforDiTaus(process, cms.InputTag('patMETs'), cms.InputTag('patPFMETs'))
#--------------------------------------------------------------------------------

common_ntuple_quantities = [
    ('VisLeg1Pt', 'leg1().pt()'),
    ('VisLeg1Eta', 'leg1().eta()'),
    ('VisLeg1Phi', 'leg1().phi()'),

    ('GenLeg1Pt', 'p4Leg1gen().pt()'),
    ('GenLeg1Eta', 'p4Leg1gen().eta()'),
    ('GenLeg1Phi', 'p4Leg1gen().phi()'),

    ('GenVisLeg1Pt', 'p4VisLeg1gen().pt()'),
    ('GenVisLeg1Eta', 'p4VisLeg1gen().eta()'),
    ('GenVisLeg1Phi', 'p4VisLeg1gen().phi()'),

    ('VisLeg2Pt', 'leg2().pt()'),
    ('VisLeg2Eta', 'leg2().eta()'),
    ('VisLeg2Phi', 'leg2().phi()'),

    ('GenLeg2Pt', 'p4Leg2gen().pt()'),
    ('GenLeg2Eta', 'p4Leg2gen().eta()'),
    ('GenLeg2Phi', 'p4Leg2gen().phi()'),

    ('GenVisLeg2Pt', 'p4VisLeg2gen().pt()'),
    ('GenVisLeg2Eta', 'p4VisLeg2gen().eta()'),
    ('GenVisLeg2Phi', 'p4VisLeg2gen().phi()'),

    ('GenPVx', 'primaryVertexPosGen().x()'),
    ('GenPVy', 'primaryVertexPosGen().y()'),
    ('GenPVz', 'primaryVertexPosGen().z()'),

    ('GenLeg1SVx', 'decayVertexPosLeg1gen().x()'),
    ('GenLeg1SVy', 'decayVertexPosLeg1gen().y()'),
    ('GenLeg1SVz', 'decayVertexPosLeg1gen().z()'),

    ('GenLeg2SVx', 'decayVertexPosLeg2gen().x()'),
    ('GenLeg2SVy', 'decayVertexPosLeg2gen().y()'),
    ('GenLeg2SVz', 'decayVertexPosLeg2gen().z()'),

    ('CollinearApproxIsValid', 'collinearApproxIsValid()'),
    ('CollinearApproxMass', 'p4CollinearApprox().mass()'),

    ('Pt', 'pt()'),
    ('Eta', 'eta()'),
    ('Phi', 'phi()'),
    ('Mass', 'mass()'),

    ('GenPt', 'p4gen().pt()'),
    ('GenPz', 'p4gen().pz()'),
    ('GenEta', 'p4gen().eta()'),
    ('GenPhi', 'p4gen().phi()'),
    ('GenMass', 'p4gen().mass()'),

    ('GenMETx', 'p4InvisGen().px()'),
    ('GenMETy', 'p4InvisGen().py()'),

    ('GenMETPt', 'p4InvisGen().pt()'),
    ('GenMETPhi', 'p4InvisGen().phi()'),

    ('GenMETParMuon', '(p4InvisGen().px()*leg1().p4().px() + p4InvisGen().py()*leg1().p4().py())/leg1().pt()'),
    ('GenMETPerpMuon', '(p4InvisGen().px()*leg1().p4().py() - p4InvisGen().py()*leg1().p4().px())/leg1().pt()'),

    ('PZeta', 'pZeta()'),
    ('PZetaVis', 'pZetaVis()'),
    
    ('METx', 'met().px()'),
    ('METy', 'met().py()'),

    ('SumET', 'metSumEt()'),

    ('METParMuon', '(met().px()*leg1().p4().px() + met().py()*leg1().p4().py())/leg1().pt()'),
    ('METPerpMuon', '(met().px()*leg1().p4().py() - met().py()*leg1().p4().px())/leg1().pt()'),

    ('METPt', 'met().pt()'),
    ('METphi', 'met().phi()'),

    ('CAngle1', 'leg1().charge()*cos(svFitSolution(0).restFrameVisThetaLeg1())'),
    ('CAngle2', 'leg2().charge()*cos(svFitSolution(0).restFrameVisThetaLeg2())'),

    ('Angle1', 'cos(svFitSolution(0).restFrameVisThetaLeg1())'),
    ('Angle2', 'cos(svFitSolution(0).restFrameVisThetaLeg2())'),

    ('CPolarization', 'leg1().charge()*cos(svFitSolution(0).restFrameVisThetaLeg1()) + leg2().charge()*cos(svFitSolution(0).restFrameVisThetaLeg2())'),

    ('Polarization', 'cos(svFitSolution(0).restFrameVisThetaLeg1()) + cos(svFitSolution(0).restFrameVisThetaLeg2())'),
]

sv_soln_quantities = [
    ('Pt', 'p4().pt()'),
    ('Pz', 'p4().pz()'),
    ('Eta', 'p4().eta()'),
    ('Phi', 'p4().phi()'),
    ('Mass', 'p4().mass()'),

    ('Migrad', 'migradStatus()'),
    ('FitStatus', 'svFitStatus()'),
    ('MScale1', 'mScale1()'),
    ('MScale2', 'mScale2()'),

    ('PVx', 'primaryVertexPosSVrefitted().x()'),
    ('PVy', 'primaryVertexPosSVrefitted().y()'),
    ('PVz', 'primaryVertexPosSVrefitted().z()'),

    ('METx', 'p4Invis().px()'),
    ('METy', 'p4Invis().py()'),

    ('METPt', 'p4Invis().pt()'),
    ('METPhi', 'p4Invis().phi()'),

    ('Likelihood', 'logLikelihood()'),
    ('METLikelihood', 'logLikelihoodMEt()'),

    ('X1', 'x1()'),
    ('X2', 'x2()'),

    ('Leg1SVx', 'decayVertexPosLeg1().x()'),
    ('Leg1SVy', 'decayVertexPosLeg1().y()'),
    ('Leg1SVz', 'decayVertexPosLeg1().z()'),

    ('Leg2SVx', 'decayVertexPosLeg2().x()'),
    ('Leg2SVy', 'decayVertexPosLeg2().y()'),
    ('Leg2SVz', 'decayVertexPosLeg2().z()'),
]

# produce ntuple
process.ntupleProducer = cms.EDAnalyzer(
    "ObjValEDNtupleProducer",
    ntupleName = cms.string("exampleNtuple"),
    sources = cms.PSet(
        # Grouping of sources is for convenience of specifying pluginTypes, etc
        hadronicTaus = cms.PSet(
            # Select multiplicy of object(s) to store
            vector = cms.bool(True), # Store a value for all objects in this collection
            #indices = cms.vuint_32([0, 1, 2]) # Store values for first, second, third objects

            # Extractor plugin
            pluginType = cms.string("PATMuTauPairVectorValExtractor"),

            # Collection to extract from
            src = cms.InputTag("selectedMuTauPairsForAHtoMuTauAntiOverlapVetoCumulative"),

            # Variables to compute for this source
            columns = cms.PSet(
            )
        ),
    )
)

for name, cut_string in common_ntuple_quantities:
    my_columns = process.ntupleProducer.sources.hadronicTaus.columns
    setattr(my_columns, name, cms.string(cut_string))

for name, cut_string in sv_soln_quantities:
    my_columns = process.ntupleProducer.sources.hadronicTaus.columns
    for solution_index in range(4):
        new_name = "SV%i%s" % (solution_index, name)
        new_cut = "svFitSolution(%i).%s" % (solution_index, cut_string)
        setattr(my_columns, new_name, cms.string(new_cut))

#print process.ntupleProducer.sources.hadronicTaus.columns
# Save ntuple
process.out = cms.OutputModule(
    "PoolOutputModule",                                                                                                                                                        
    outputCommands = cms.untracked.vstring("drop *", "keep *_*ntupleProducer*_*_*" ),
    verbose = cms.untracked.bool(False),
    fileName = cms.untracked.string("sv_ntuple.root")      
)

process.p = cms.Path(
   process.producePatTupleAHtoMuTauSpecific
# + process.printGenParticleList # uncomment to enable print-out of generator level particles
# + process.printEventContent    # uncomment to enable dump of event content after PAT-tuple production
  + process.selectAHtoMuTauEvents
  + process.ntupleProducer
#  + process.analyzeAHtoMuTauEvents
#  + process.saveAHtoMuTauPlots
)
process.end = cms.EndPath(process.out)

#--------------------------------------------------------------------------------
# import utility function for switching HLT InputTags when processing
# RECO/AOD files produced by MCEmbeddingTool
from TauAnalysis.MCEmbeddingTools.tools.switchInputTags import switchInputTags
#
# comment-out to switch HLT InputTags
#switchInputTags(process)
#--------------------------------------------------------------------------------
 
#--------------------------------------------------------------------------------
# import utility function for factorization
from TauAnalysis.Configuration.tools.factorizationTools import enableFactorization_runAHtoMuTau
#
# define "hook" for enabling/disabling factorization
# in case running jobs on the CERN batch system
# (needs to be done after process.p has been defined)
#__#factorization#
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
# import utility function for disabling estimation of systematic uncertainties
#
# NOTE: per default, estimation of systematic uncertainties is **enabled** per default
#
from TauAnalysis.Configuration.tools.sysUncertaintyTools import disableSysUncertainties_runAHtoMuTau
#from TauAnalysis.Configuration.tools.sysUncertaintyTools import enableSysUncertainties_runAHtoMuTau
#
# define "hook" for keeping enabled/disabling estimation of systematic uncertainties
# in case running jobs on the CERN batch system
# (needs to be done after process.p has been defined)
#__#systematics#
if not hasattr(process, "isBatchMode"):
    disableSysUncertainties_runAHtoMuTau(process)
    #enableSysUncertainties_runAHtoMuTau(process)
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
# disable event-dump output
# in order to reduce size of log-files
if hasattr(process, "disableEventDump"):
    process.analyzeAHtoMuTauEvents_woBtag.eventDumps = cms.VPSet()
    process.analyzeAHtoMuTauEvents_wBtag.eventDumps = cms.VPSet()
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
#
process.producePatTupleAll = cms.Sequence(process.producePatTuple + process.producePatTupleAHtoMuTauSpecific)
#
# define "hook" for enabling/disabling production of PAT-tuple event content,
# depending on whether RECO/AOD or PAT-tuples are used as input for analysis
#
#__#patTupleProduction#
if not hasattr(process, "isBatchMode"):
    process.p.replace(process.producePatTupleAHtoMuTauSpecific, process.producePatTuple + process.producePatTupleAHtoMuTauSpecific)
#--------------------------------------------------------------------------------

# print-out all python configuration parameter information
#
# NOTE: need to delete empty sequence produced by call to "switchJetCollection"
#       in order to avoid error when calling "process.dumpPython"
#      ( cf. https://hypernews.cern.ch/HyperNews/CMS/get/physTools/1688/1/1/1/1/1.html )
#
#del process.patJetMETCorrections
#print process.dumpPython()
