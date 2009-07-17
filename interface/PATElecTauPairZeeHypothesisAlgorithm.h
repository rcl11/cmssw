#ifndef TauAnalysis_CandidateTools_PATElecTauPairZeeHypothesisAlgorithm_h
#define TauAnalysis_CandidateTools_PATElecTauPairZeeHypothesisAlgorithm_h

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/InputTag.h"

#include "DataFormats/Common/interface/Ptr.h"

#include "DataFormats/Candidate/interface/Particle.h" 

#include "AnalysisDataFormats/TauAnalysis/interface/PATElecTauPairZeeHypothesis.h"

class PATElecTauPairZeeHypothesisAlgorithm
{
 public:

  PATElecTauPairZeeHypothesisAlgorithm(const edm::ParameterSet&);
  ~PATElecTauPairZeeHypothesisAlgorithm();

  PATElecTauPairZeeHypothesis buildZeeHypothesis(edm::Ptr<PATElecTauPair>, const edm::Event&, const edm::EventSetup&);

 private:

  edm::InputTag srcGenElectrons_;
  edm::InputTag srcCaloJets_;
  edm::InputTag srcPFJets_;
  edm::InputTag srcTracks_;
  edm::InputTag srcGsfElectrons_;
  edm::InputTag srcGsfTracks_;

  int tkminPixelHits_;
  int tkminTrackerHits_;
  double tkmaxChi2_;

  double dRmatch_;

  int verbosity_;  
};

#endif 

