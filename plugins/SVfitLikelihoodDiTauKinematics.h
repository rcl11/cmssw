#ifndef TauAnalysis_CandidateTools_SVfitLikelihoodDiTauKinematics_h
#define TauAnalysis_CandidateTools_SVfitLikelihoodDiTauKinematics_h

/** \class SVfitLikelihoodDiTauKinematics
 *
 * Plugin for computing likelihood for decay kinematics of tau lepton pair;
 * different types of likelihoods (isotropic/unpolarized/polarized decay) 
 * can be specified via plugins for the two tau lepton decay "legs";
 * plugin is used by SVfit algorithm
 * 
 * \author Evan Friis, Christian Veelken; UC Davis
 *
 * \version $Revision: 1.5 $
 *
 * $Id: SVfitLikelihoodDiTauKinematics.h,v 1.5 2010/09/21 09:01:43 veelken Exp $
 *
 */

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "TauAnalysis/CandidateTools/interface/SVfitDiTauLikelihoodBase.h"
#include "TauAnalysis/CandidateTools/interface/SVfitLegLikelihoodBase.h"

#include "AnalysisDataFormats/TauAnalysis/interface/CompositePtrCandidateT1T2MEt.h"
#include "AnalysisDataFormats/TauAnalysis/interface/SVfitDiTauSolution.h"

#include <iomanip>

template <typename T1, typename T2>
  class SVfitLikelihoodDiTauKinematics : public SVfitDiTauLikelihoodBase<T1,T2>
{
 public:
  SVfitLikelihoodDiTauKinematics(const edm::ParameterSet&);
  ~SVfitLikelihoodDiTauKinematics();

  virtual void beginJob();
  virtual void beginEvent(edm::Event&, const edm::EventSetup&);
  virtual void beginCandidate(const CompositePtrCandidateT1T2MEt<T1,T2>&);

  void print(std::ostream&) const;

  bool isFittedParameter(int) const;
  bool supportsPolarization() const;

  double operator()(const CompositePtrCandidateT1T2MEt<T1,T2>&, const SVfitDiTauSolution&) const;

 private:
  SVfitLegLikelihoodBase<T1>* leg1Likelihood_;
  SVfitLegLikelihoodBase<T2>* leg2Likelihood_;
};

#endif
