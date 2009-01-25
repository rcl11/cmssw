#ifndef TauAnalysis_CandidateTools_DiTauAntiOverlapSelector_h
#define TauAnalysis_CandidateTools_DiTauAntiOverlapSelector_h

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/InputTag.h"

#include "DataFormats/Common/interface/Handle.h"
#include "AnalysisDataFormats/TauAnalysis/interface/DiTauCandidateFwd.h"
#include "AnalysisDataFormats/TauAnalysis/interface/DiTauCandidate.h"

#include "PhysicsTools/Utilities/interface/deltaR.h"

class DiTauAntiOverlapSelectorImp
{
 public:
  typedef DiTauCandidateCollection collection;

  DiTauAntiOverlapSelectorImp(const edm::ParameterSet& cfg) 
    : dRmin_(cfg.getParameter<double>("dRmin")) {}

  std::vector<const DiTauCandidate*>::const_iterator begin() const { return selected_.begin(); }
  std::vector<const DiTauCandidate*>::const_iterator end() const { return selected_.end(); }

  void select(const edm::Handle<collection>& diTaus, const edm::Event& evt, const edm::EventSetup& es) 
  {
    selected_.clear();
  
    for ( collection::const_iterator diTau = diTaus->begin();
	  diTau != diTaus->end(); ++diTau ) {
      reco::CandidatePtr firstDecayProduct = diTau->hadrTauPtr();
      reco::CandidatePtr secondDecayProduct = diTau->leptTauPtr();

      double dR = reco::deltaR(firstDecayProduct->p4(), secondDecayProduct->p4());

      if ( dR > dRmin_ ) selected_.push_back(&(*diTau)); 
    }
  }

  size_t size() const { return selected_.size(); }

private:
  std::vector<const DiTauCandidate*> selected_;

  double dRmin_;
};

#endif
