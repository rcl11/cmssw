import FWCore.ParameterSet.Config as cms

#from SimGeneral.HepPDTESSource.pdt_cfi import *

#--------------------------------------------------------------------------------
# produce data-formats providing information 
# about compatibility of an elecron + tau-jet pair
# with the hypothesis of being a pair of electron,
# resulting from a Z --> e+ e- decay
#--------------------------------------------------------------------------------

elecTauPairZeeHypotheses = cms.EDProducer("PATElecTauPairZeeHypothesisProducer",
    elecTauPairSource = cms.InputTag('selectedElecTauPairsMt1METcumulative'),

    genElectronsFromZsSource = cms.InputTag('genElectronsFromZs'),
                                          
    caloJetSource = cms.InputTag('iterativeCone5CaloJets'),
    pfJetSource = cms.InputTag('iterativeCone5PFJets'),                                       
    trackSource = cms.InputTag('generalTracks'),
    gsfElectronSource = cms.InputTag('pixelMatchGsfElectrons'),
    gsfTrackSource = cms.InputTag('pixelMatchGsfFit'),
                                          
    tkminPixelHits = cms.int32(1),
    tkminTrackerHits = cms.int32(8),	
    tkmaxChi2 = cms.double(100.),

    dRmatch = cms.double(0.5),

    verbosity = cms.untracked.int32(1)                                      
)

produceElecTauPairZeeHypotheses = cms.Sequence( elecTauPairZeeHypotheses )
