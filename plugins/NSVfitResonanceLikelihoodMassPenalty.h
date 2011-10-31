#ifndef TauAnalysis_CandidateTools_plugins_NSVfitResonanceLikelihoodMassPenalty_h
#define TauAnalysis_CandidateTools_plugins_NSVfitResonanceLikelihoodMassPenalty_h

/** \class NSVfitResonanceLikelihoodMassPenalty
 *
 * Adds a penalty term for high masses.  
 * The return value of operator() is configurable via python.
 *
 * Options used in the past are:
 *  o TMath::Log(mass)
 *  o [0]*([1] - TMath::Erf((x - [2])*[3]))
 *    with p0 = 4.21e-2, p1 = 2.52e-2, p2 = 4.40e+1, p3 = -6.90e-3 
 *   (efficiency of gg --> Higgs --> mu + tau_had channel in 2010 analysis)
 *  o [0]*([1] - TMath::Erf((x - [2])*[3]))
 *    with p0 = 2.49e-2, p1 = 7.78e-2, p2 = 5.63e+1, p3 = -7.53e-3 
 *   (efficiency of gg --> Higgs --> e + tau_had channel in 2010 analysis)
 *
 * \author Evan Friis, Christian Veelken; UC Davis
 *
 */

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "TauAnalysis/CandidateTools/interface/NSVfitResonanceLikelihood.h"

#include "AnalysisDataFormats/TauAnalysis/interface/NSVfitResonanceHypothesis.h"

#include "TFormula.h"

class NSVfitResonanceLikelihoodMassPenalty : public NSVfitResonanceLikelihood
{
 public:

  NSVfitResonanceLikelihoodMassPenalty(const edm::ParameterSet&);
  ~NSVfitResonanceLikelihoodMassPenalty();

  void beginJob(NSVfitAlgorithmBase*) const {}

  void beginCandidate(const NSVfitResonanceHypothesis*) const {};

  double operator()(const NSVfitResonanceHypothesis*) const;

 private:

  TFormula* formula_;

  double power_;
};


#endif /* end of include guard: TauAnalysis_CandidateTools_plugins_NSVfitResonanceLikelihoodMassPenalty2_h */
