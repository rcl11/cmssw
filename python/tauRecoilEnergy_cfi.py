import FWCore.ParameterSet.Config as cms

#--------------------------------------------------------------------------------
# produce data-formats providing information 
# about Et of jets (CaloTowers) opposite in transverse plane to tau-jet
#--------------------------------------------------------------------------------

tauRecoilEnergyFromJets = cms.EDProducer("PATTauRecoilEnergyFromJetsProducer",
    srcLeptons = cms.InputTag('selectedLayer1TausProngCumulative'),
    srcEnergyObjects = cms.InputTag('cleanLayer1Jets'),
    #srcEnergyObjects = cms.InputTag('towerMaker'),
    etaMin = cms.double(-2.5),
    etaMax = cms.double(+2.5),
    etMin = cms.double(0.),
    dPhiMin = cms.double(2.89),
    verbosity = cms.untracked.int32(0)
)

tauRecoilEnergyFromCaloTowers = cms.EDProducer("PATTauRecoilEnergyFromCaloTowersProducer",
    srcLeptons = cms.InputTag('selectedLayer1TausProngCumulative'),
    srcEnergyObjects = cms.InputTag('towerMaker'),
    etaMin = cms.double(-5.),
    etaMax = cms.double(+5.),
    etMin = cms.double(0.5),
    dPhiMin = cms.double(2.89),
    verbosity = cms.untracked.int32(0)
)

produceTauRecoilEnergy = cms.Sequence( tauRecoilEnergyFromJets * tauRecoilEnergyFromCaloTowers )
