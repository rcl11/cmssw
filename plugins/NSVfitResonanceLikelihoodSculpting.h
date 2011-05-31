#ifndef TauAnalysis_CandidateTools_NSVfitResonanceLikelihoodSculpting_h
#define TauAnalysis_CandidateTools_NSVfitResonanceLikelihoodSculpting_h

/** \class NSVfitResonanceLikelihoodSculpting
 *
 * Plugin for computing likelihood of a resonance mass
 * given the reconstructed visible mass.
 *
 * \author Evan Friis, Christian Veelken; UC Davis
 *
 * \version $Revision: 1.2 $
 *
 * $Id: NSVfitResonanceLikelihoodSculpting.h,v 1.2 2011/03/03 13:04:47 veelken Exp $
 *
 */

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "TauAnalysis/CandidateTools/interface/NSVfitResonanceLikelihood.h"

#include "AnalysisDataFormats/TauAnalysis/interface/NSVfitResonanceHypothesis.h"

#include <TFormula.h>

class NSVfitResonanceLikelihoodSculpting : public NSVfitResonanceLikelihood
{
 public:
  NSVfitResonanceLikelihoodSculpting(const edm::ParameterSet&);
  ~NSVfitResonanceLikelihoodSculpting();

  void beginJob(NSVfitAlgorithmBase*) const;

  double operator()(const NSVfitResonanceHypothesis*) const;

 private:
  TFormula meanFunc_;
  TFormula rmsFunc_;
};

#endif
