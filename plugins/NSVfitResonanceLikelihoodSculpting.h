#ifndef TauAnalysis_CandidateTools_NSVfitResonanceLikelihoodSculpting_h
#define TauAnalysis_CandidateTools_NSVfitResonanceLikelihoodSculpting_h

/** \class NSVfitResonanceLikelihoodSculpting
 *
 * Plugin for computing likelihood of a resonance mass
 * given the reconstructed visible mass.
 *
 * \author Evan Friis, Christian Veelken; UC Davis
 *
 * \version $Revision: 1.1 $
 *
 * $Id: NSVfitResonanceLikelihoodSculpting.h,v 1.1 2011/05/31 09:00:47 friis Exp $
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
  double power_;
  bool normalize_;
};

#endif
