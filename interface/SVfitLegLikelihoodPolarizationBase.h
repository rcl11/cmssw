#ifndef TauAnalysis_CandidateTools_SVfitLegLikelihoodPolarizationBase_h
#define TauAnalysis_CandidateTools_SVfitLegLikelihoodPolarizationBase_h

/** \class SVfitLegLikelihoodPolarizationBase
 *
 * Base-class for plugins computing likelihood for tau lepton decay "leg"
 * to be compatible with decay of polarized tau lepton,
 * assuming  matrix element of V-A electroweak decay
 *
 * \author Christian Veelken, UC Davis
 *
 * \version $Revision: 1.2 $
 *
 * $Id: SVfitLegLikelihoodPolarizationBase.h,v 1.2 2010/09/08 13:28:09 veelken Exp $
 *
 */

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "TauAnalysis/CandidateTools/interface/SVfitLegLikelihoodBase.h"

#include "AnalysisDataFormats/TauAnalysis/interface/SVfitLegSolution.h"

#include <string>
#include <iostream>

template <typename T>
class SVfitLegLikelihoodPolarizationBase : public SVfitLegLikelihoodBase<T>
{
 public:
  SVfitLegLikelihoodPolarizationBase(const edm::ParameterSet& cfg) 
    : SVfitLegLikelihoodBase<T>(cfg)
  {
    usePolarization_ = cfg.getParameter<bool>("usePolarization");
  }
  virtual ~SVfitLegLikelihoodPolarizationBase() {}

  virtual void beginEvent(const edm::Event&, const edm::EventSetup&) {}

  virtual bool supportsPolarization() const
  {
    return usePolarization_;
  }

  virtual double operator()(const T& leg, const SVfitLegSolution& solution) const
  {
    if ( usePolarization_ ) {
//--- determine whether visible tau decay products come from a tau- or a tau+ decay;
//    depending on whether the tau lepton is a tau- or a tau+,
//    the following convention is used to relate the "handed-ness" of the tau to the polarization:
//   o P(tau-_{L}) = P(tau+_{R}) = -1
//   o P(tau-_{R}) = P(tau+_{L}) = +1
//
      double tauPolarizationSign = 0.;
      if      ( leg.charge() < 0. ) tauPolarizationSign = +1.; // tau- case
      else if ( leg.charge() > 0. ) tauPolarizationSign = -1.; // tau+ case
      else {
	edm::LogWarning ("SVfitLegLikelihoodPolarizationBase::operator()") 
	  << " Failed to identify whether visible tau decay products come from a tau- or a tau+ decay"
	  << " --> returning log-likelihood value for unpolarized case !!";
	return negLogLikelihoodUnpolarized(leg, solution);
      }

      if ( solution.polarizationHypothesis() == SVfitLegSolution::kLeftHanded  ) {
	return negLogLikelihoodPolarized(leg, solution, -1.*tauPolarizationSign);
      } else if ( solution.polarizationHypothesis() == SVfitLegSolution::kRightHanded ) {
	return negLogLikelihoodPolarized(leg, solution, +1.*tauPolarizationSign);
      } else {
	edm::LogWarning ("SVfitLegLikelihoodPolarizationBase::operator()") 
	  << " Unknown polarization of tau lepton" 
	  << " --> returning log-likelihood value for unpolarized case !!";
	return negLogLikelihoodUnpolarized(leg, solution);
      }
    } else {
      return negLogLikelihoodUnpolarized(leg, solution);
    }
  }
 protected:
  double negLogLikelihoodUnpolarized(const T& leg, const SVfitLegSolution& solution) const 
  {
    return 0.5*(negLogLikelihoodPolarized(leg, solution, +1.) +  negLogLikelihoodPolarized(leg, solution, -1.));
  }

  virtual double negLogLikelihoodPolarized(const T&, const SVfitLegSolution&, double) const = 0;

  bool usePolarization_;
};

#endif
