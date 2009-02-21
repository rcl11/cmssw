import FWCore.ParameterSet.Config as cms
import copy

from TauAnalysis.RecoTools.pftauPatSelector_cfi import *

# require muon and tau-jet to be separated in eta-phi,
# in order to ensure that both do not refer to one and the same physical particle
# (NOTE: cut is already applied during skimming,
#        so should not reject any events)
selectedElecTauPairsAntiOverlapVeto = cms.EDFilter("PATElecTauPairSelector",
     src = cms.InputTag("allElecTauPairs"),
     cut = cms.string('dR12 > 0.7'),
     filter = cms.bool(False)
)

# require muon and tau not to be back-to-back
selectedElecTauPairsAcoplanarityIndividual = cms.EDFilter("PATElecTauPairSelector",
     src = selectedElecTauPairsAntiOverlapVeto.src,
     cut = cms.string('dPhi1MET < 2.4'),
     filter = cms.bool(False)
)

selectedElecTauPairsAcoplanarityCumulative = copy.deepcopy(selectedElecTauPairsAcoplanarityIndividual)
selectedElecTauPairsAcoplanarityCumulative.src = cms.InputTag("selectedElecTauPairsAntiOverlapVeto")

# require muon and tau to form a zero-charge pair
selectedElecTauPairsZeroChargeIndividual = cms.EDFilter("PATElecTauPairSelector",
     src = selectedElecTauPairsAntiOverlapVeto.src,
     cut = cms.string('charge = 0'),
     #cut = cms.string('(leg1.charge + leg2.leadTrack.charge) = 0'), # NOTE: to be used for background studies only !!                    
     filter = cms.bool(False)
)

selectedElecTauPairsZeroChargeCumulative = copy.deepcopy(selectedElecTauPairsZeroChargeIndividual)
selectedElecTauPairsZeroChargeCumulative.src = cms.InputTag("selectedElecTauPairsAcoplanarityCumulative")

selectElecTauPairs = cms.Sequence( selectedElecTauPairsAntiOverlapVeto
                                  *selectedElecTauPairsAcoplanarityIndividual
                                  *selectedElecTauPairsAcoplanarityCumulative
                                  *selectedElecTauPairsZeroChargeIndividual
                                  *selectedElecTauPairsZeroChargeCumulative )
