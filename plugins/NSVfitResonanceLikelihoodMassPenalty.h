#ifndef TauAnalysis_CandidateTools_plugins_NSVfitResonanceLikelihoodMassPenalty_h
#define TauAnalysis_CandidateTools_plugins_NSVfitResonanceLikelihoodMassPenalty_h

/** \class NSVfitResonanceLikelihoodMassPenalty
 *
 * Adds a penalty term for high masses.  The return value of operator()
 * is
 *
 * -penaltyFactor_*Log(resonanceMass)
 *
 * \author Evan Friis, Christian Veelken; UC Davis
 *
 */

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "TauAnalysis/CandidateTools/interface/NSVfitResonanceLikelihood.h"

#include "AnalysisDataFormats/TauAnalysis/interface/NSVfitResonanceHypothesis.h"

class NSVfitResonanceLikelihoodMassPenalty : public NSVfitResonanceLikelihood
{
 public:
  NSVfitResonanceLikelihoodMassPenalty(const edm::ParameterSet&);
  ~NSVfitResonanceLikelihoodMassPenalty() {};

  void beginJob(NSVfitAlgorithmBase*) const {}

  void beginCandidate(const NSVfitResonanceHypothesis*) const {};

  double operator()(const NSVfitResonanceHypothesis*) const;

 private:
  double power_;
};


#endif /* end of include guard: TauAnalysis_CandidateTools_plugins_NSVfitResonanceLikelihoodMassPenalty2_h */
