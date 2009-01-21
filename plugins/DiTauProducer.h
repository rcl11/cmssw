#ifndef __HiggsAnalysis_CandidateTools_DiTauProducer__
#define __HiggsAnalysis_CandidateTools_DiTauProducer__

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "AnalysisDataFormats/TauAnalysis/interface/DiTauCandidate.h"

#include "DataFormats/Math/interface/LorentzVector.h"

#include <stdio.h>


class DiTauProducer : public edm::EDProducer {
 public:


  enum MetModes {
    COLLINEAR_APPROX,
    NO_MET, 
    TRANSVERSE_RECO
  };

  explicit DiTauProducer(const edm::ParameterSet&);
  ~DiTauProducer();

   
  
 private:
  virtual void beginJob(const edm::EventSetup&) ;
  virtual void produce(edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
    
  bool 
    computeDiTau4P( const math::XYZTLorentzVector& tau1, 
		    const math::XYZTLorentzVector& tau2, 
		    const math::XYZTLorentzVector& met,
		    double& a1, double &a2,
		    math::XYZTLorentzVector& diTau ) const;

  typedef reco::CandidateView::const_iterator CVIT; 
  typedef std::vector<DiTauCandidate> DiTauCollection;
  void fillDiTau( std::auto_ptr<DiTauCollection>& diTauCollection,
		  const reco::CandidatePtr& lTauIt, 
		  const reco::CandidatePtr& hTauIt, 
		  const reco::CandidatePtr& metIt) const;

  bool checkCombination(const reco::Candidate& tau1, 
			const reco::Candidate& tau2) const;

  // compute cosine of angle phi between tau1 and tau2 (in transversal plane)
  double cosPhiT( const math::XYZTLorentzVector& tau1,
		  const math::XYZTLorentzVector& tau2 ) const; 

  // ----------member data ---------------------------

  edm::InputTag   inputTagHadronicTaus_;
  edm::InputTag   inputTagLeptonicTaus_;
  edm::InputTag   inputTagMETs_;

  bool            verbose_;

  bool            useLeadingTaus_;
  int             metMode_;



};

#endif

