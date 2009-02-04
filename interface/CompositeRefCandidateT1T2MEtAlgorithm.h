#ifndef TauAnalysis_CandidateTools_CompositeRefCandidateT1T2MEtAlgorithm_h
#define TauAnalysis_CandidateTools_CompositeRefCandidateT1T2MEtAlgorithm_h

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "PhysicsTools/Utilities/interface/deltaR.h"
#include "DataFormats/Math/interface/normalizedPhi.h"

#include "AnalysisDataFormats/TauAnalysis/interface/CompositeRefCandidateT1T2MEt.h"

#include "DataFormats/Candidate/interface/CandidateFwd.h" 
#include "DataFormats/Candidate/interface/Candidate.h" 

#include <TMath.h>

template<typename T1, typename T2>
class CompositeRefCandidateT1T2MEtAlgorithm 
{
  typedef std::vector<T1> T1Collection;
  typedef edm::Ref<T1Collection> T1Ref;
  typedef std::vector<T2> T2Collection;
  typedef edm::Ref<T2Collection> T2Ref;

 public:

  CompositeRefCandidateT1T2MEtAlgorithm(const edm::ParameterSet& cfg)
  {
    recoMode_ = cfg.getParameter<std::string>("recoMode");
    verbosity_ = cfg.getUntrackedParameter<int>("verbosity", 0);
  }
  ~CompositeRefCandidateT1T2MEtAlgorithm() {}

  CompositeRefCandidateT1T2MEt<T1,T2> buildCompositeRefCandidate(const T1Ref leg1, 
								 const T2Ref leg2, 
								 const reco::CandidatePtr met)
  {
    CompositeRefCandidateT1T2MEt<T1,T2> compositeRefCandidate(leg1, leg2, met);
  
    if ( leg1.isNull() || 
         leg2.isNull() ) {
      edm::LogError ("CompositeRefCandidateT1T2MEtAlgorithm") << " References to visible Decay products invalid !!";
      return compositeRefCandidate;
    }

//--- compute quantities that are independent of MET	
    compositeRefCandidate.setCharge(leg1->charge() + leg2->charge());
    compositeRefCandidate.setP4Vis(leg1->p4() + leg2->p4());
    compositeRefCandidate.setDR12(reco::deltaR(leg1->p4(), leg2->p4()));
    compositeRefCandidate.setDPhi12(normalizedPhi(leg1->phi() - leg2->phi()));
    compositeRefCandidate.setVisEtaMin(TMath::Min(leg1->eta(), leg2->eta()));
    compositeRefCandidate.setVisEtaMax(TMath::Max(leg1->eta(), leg2->eta()));

//--- compute quantities that do dependent on MET
    if ( met.isNonnull() ) {
      compCollinearApprox(compositeRefCandidate, leg1->p4(), leg2->p4(), met->px(), met->py());
    
      compositeRefCandidate.setP4CDFmethod(compP4CDFmethod(leg1->p4(), leg2->p4(), met->px(), met->py()));
      compositeRefCandidate.setMt12MET(compMt(leg1->p4(), leg2->p4(), met->px(), met->py()));    
      compositeRefCandidate.setMt1MET(compMt(leg1->p4(), met->px(), met->py()));
      compositeRefCandidate.setMt2MET(compMt(leg2->p4(), met->px(), met->py()));
      compositeRefCandidate.setDPhi1MET(normalizedPhi(leg1->phi() - met->phi()));
      compositeRefCandidate.setDPhi2MET(normalizedPhi(leg2->phi() - met->phi()));
    } else {
      compositeRefCandidate.setCollinearApproxQuantities(reco::Candidate::LorentzVector(0,0,0,0), -1, -1, false);
    }
 
//--- set compositeRef four-momentum
//    (depending on recoMode configuration parameter)
    if ( recoMode_ == "collinearApprox" ) {
      if ( met.isNonnull() ) {
        compositeRefCandidate.setP4(compositeRefCandidate.p4CollinearApprox());
      } else {
        edm::LogError ("buildCompositeRefCandidate") << " Failed to set four-momentum:"
						     << " recoMode = " << recoMode_ << " requires MET pointer to be valid !!";
      }
    } else if ( recoMode_ == "cdfMethod" ) {
      if ( met.isNonnull() ) {
	compositeRefCandidate.setP4(compositeRefCandidate.p4CDFmethod());
      } else {
	edm::LogError ("buildCompositeRefCandidate") << " Failed to set four-momentum:"
						     << " recoMode = " << recoMode_ << " requires MET pointer to be valid !!";
      }
    } else if ( recoMode_ == "" ) {
      compositeRefCandidate.setP4(compositeRefCandidate.p4Vis());
    } else {
      edm::LogError ("buildCompositeRefCandidate") << " Failed to set four-momentum:"
						   << " recoMode = " << recoMode_ << " undefined !!";
    }  
    
    return compositeRefCandidate;
  }

 private: 

  void compCollinearApprox(CompositeRefCandidateT1T2MEt<T1,T2>& compositeRefCandidate,
			   const reco::Candidate::LorentzVector& leg1,
			   const reco::Candidate::LorentzVector& leg2,
			   double metPx, double metPy)
  {
    double x1_numerator = leg1.px()*leg2.py() - leg2.px()*leg1.py();
    double x1_denominator = leg2.py()*(leg1.px() + metPx) - leg2.px()*(leg1.py() + metPy);
    double x1 = ( x1_denominator != 0. ) ? x1_numerator/x1_denominator : -1.;

    double x2_numerator = x1_numerator;
    double x2_denominator = leg1.px()*(leg2.py() + metPy) - leg1.py()*(leg2.px() + metPx);
    double x2 = ( x2_denominator != 0. ) ? x2_numerator/x2_denominator : -1.;

    if ( (x1 > 0. && x1 < 1.) &&
         (x2 > 0. && x2 < 1.) ) {
      reco::Candidate::LorentzVector p4 = leg1/x1 + leg2/x2;
      compositeRefCandidate.setCollinearApproxQuantities(p4, x1, x2, true);
    } else {
      compositeRefCandidate.setCollinearApproxQuantities(reco::Candidate::LorentzVector(0,0,0,0), -1, -1, false);
    }
  }
  reco::Candidate::LorentzVector compP4CDFmethod(const reco::Candidate::LorentzVector& leg1, 
						 const reco::Candidate::LorentzVector& leg2, 
						 double metPx, double metPy)
  {
    double px = leg1.px() + leg2.px() + metPx;
    double py = leg1.py() + leg2.py() + metPy;
    double pz = leg1.pz() + leg2.pz();
    double e = leg1.energy() + leg2.energy() + TMath::Sqrt(metPx*metPx + metPy*metPy);
    reco::Candidate::LorentzVector p4(px, py, pz, e);
    return p4;
  }
  double compMt(const reco::Candidate::LorentzVector& leg1, 
		const reco::Candidate::LorentzVector& leg2, 
		double metPx, double metPy)
  {
    double px = leg1.px() + leg2.px() + metPx;
    double py = leg1.py() + leg2.py() + metPy;
    double et = leg1.Et() + leg2.Et() + TMath::Sqrt(metPx*metPx + metPy*metPy);
    double mt2 = et*et - (px*px + py*py);
    if ( mt2 < 0 ) {
      edm::LogError ("compMt") << " mt2 = " << mt2 << " must not be negative !!";
      return -1.;
    }
    return TMath::Sqrt(mt2);
  }
  double compMt(const reco::Candidate::LorentzVector& visParticle, 
		double metPx, double metPy)
  {
    double px = visParticle.px() + metPx;
    double py = visParticle.py() + metPy;
    double et = visParticle.Et() + TMath::Sqrt(metPx*metPx + metPy*metPy);
    double mt2 = et*et - (px*px + py*py);
    if ( mt2 < 0 ) {
      edm::LogError ("compMt") << " mt2 = " << mt2 << " must not be negative !!";
      return -1.;
    }
    return TMath::Sqrt(mt2);
  }

  std::string recoMode_;
  int verbosity_;
};
#endif 

