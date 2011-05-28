#ifndef TauAnalysis_CandidateTools_plugins_NSVfitResonanceLikelihoodPtBalance2_h
#define TauAnalysis_CandidateTools_plugins_NSVfitResonanceLikelihoodPtBalance2_h

/** \class NSVfitResonanceLikelihoodPtBalance2
 *
 * New and improved plugin for computing likelihood for two tau leptons
 * to balance each other in transverse momentum
 *
 * \author Evan Friis, Christian Veelken; UC Davis
 *
 * \version $Revision: $
 *
 * $Id: $
 *
 */

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "TauAnalysis/CandidateTools/interface/NSVfitResonanceLikelihood.h"
#include "TauAnalysis/CandidateTools/interface/NSVfitPtBalanceCompositePdf.h"

#include "AnalysisDataFormats/TauAnalysis/interface/NSVfitResonanceHypothesis.h"


class NSVfitResonanceLikelihoodPtBalance2 : public NSVfitResonanceLikelihood
{
 public:
  NSVfitResonanceLikelihoodPtBalance2(const edm::ParameterSet&);
  ~NSVfitResonanceLikelihoodPtBalance2();

  void beginJob(NSVfitAlgorithmBase*) const {}

  void beginCandidate(const NSVfitResonanceHypothesis*) const;

  double operator()(const NSVfitResonanceHypothesis*) const;

 private:
  NSVfitPtBalanceCompositePdf leg1PDF_;
  NSVfitPtBalanceCompositePdf leg2PDF_;
  bool addMassFactor_;
  double maxMass_;
  double minMass_;
  mutable double highestLikelihood_;
};


#endif /* end of include guard: TauAnalysis_CandidateTools_plugins_NSVfitResonanceLikelihoodPtBalance22_h */
