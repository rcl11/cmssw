#ifndef TauAnalysis_CandidateTools_plugins_NSVfitResonanceLikelihoodSculpting2_h
#define TauAnalysis_CandidateTools_plugins_NSVfitResonanceLikelihoodSculpting2_h

/** \class NSVfitResonanceLikelihoodSculpting2
 *
 * Correct tau decay likelihood for effect of visible Pt cuts
 *
 * \author Christian Veelken; LLR
 *
 */

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "TauAnalysis/CandidateTools/interface/NSVfitResonanceLikelihood.h"

#include "AnalysisDataFormats/TauAnalysis/interface/NSVfitResonanceHypothesis.h"

#include <TFile.h>
#include <TH2.h>

class NSVfitResonanceLikelihoodSculpting2 : public NSVfitResonanceLikelihood
{
 public:

  NSVfitResonanceLikelihoodSculpting2(const edm::ParameterSet&);
  ~NSVfitResonanceLikelihoodSculpting2();

  void beginJob(NSVfitAlgorithmBase*) const {}

  void beginCandidate(const NSVfitResonanceHypothesis*) const {};

  double operator()(const NSVfitResonanceHypothesis*) const;

 private:

  double getProbabilityX1andX2gt(double, double) const;

  std::string inputFileName_;
  TFile* inputFile_;
  std::string histogramNameX1_;
  std::string histogramNameX2_;

  TH2* histogramX1vsX2_;

  double minVisPt1_;
  double minVisPt2_;
  double minVisPt1alt_;
  double minVisPt2alt_;

  double power_;
};


#endif /* end of include guard: TauAnalysis_CandidateTools_plugins_NSVfitResonanceLikelihoodSculpting2_h */
