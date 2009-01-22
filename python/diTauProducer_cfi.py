import FWCore.ParameterSet.Config as cms


diTauProducer = cms.EDProducer(
    "DiTauProducer",

    # collection of Candidate based objects acting as the hadronic tau
    hadronicTaus = cms.InputTag(''),

    # collection of Candidate based objects acting as the leptonic tau
    leptonicTaus = cms.InputTag(''),

    # collection of Candidate based objects acting as MET
    # must contain only one object
    METs = cms.InputTag(''),

    #   enum MetModes {
    #    COLLINEAR_APPROX, 
    #    NO_MET, 
    #    TRANSVERSE_RECO
    #  };
    metMode = cms.int32(0),

    # if true, use the highest pT tau on both legs.
    # otherwise, try all combinations
    useLeadingTaus = cms.bool(True),
    
    verbose =  cms.untracked.bool(False)
    )
