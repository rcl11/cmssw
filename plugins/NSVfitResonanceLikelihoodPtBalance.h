#ifndef TauAnalysis_CandidateTools_NSVfitResonanceLikelihoodPtBalance_h
#define TauAnalysis_CandidateTools_NSVfitResonanceLikelihoodPtBalance_h

/** \class NSVfitResonanceLikelihoodPtBalance
 *
 * Plugin for computing likelihood for two tau leptons
 * to balance each other in transverse momentum
 *
 * \author Evan Friis, Christian Veelken; UC Davis
 *
 * \version $Revision: 1.2 $
 *
 * $Id: NSVfitResonanceLikelihoodPtBalance.h,v 1.2 2011/03/03 13:04:47 veelken Exp $
 *
 */

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "TauAnalysis/CandidateTools/interface/NSVfitResonanceLikelihood.h"

#include "AnalysisDataFormats/TauAnalysis/interface/NSVfitResonanceHypothesis.h"

#include <TFormula.h>

class NSVfitResonanceLikelihoodPtBalance : public NSVfitResonanceLikelihood
{
 public:
  NSVfitResonanceLikelihoodPtBalance(const edm::ParameterSet&);
  ~NSVfitResonanceLikelihoodPtBalance();
  
  void beginJob(NSVfitAlgorithmBase*) const;

  double operator()(const NSVfitResonanceHypothesis*) const;
  
 private:  
  double power_;

  struct Parameterization {
    Parameterization(const edm::ParameterSet& cfg);
    std::auto_ptr<TFormula> smear_;
    std::auto_ptr<TFormula> gaussFrac_;
    std::auto_ptr<TFormula> turnOnThreshold_;
    std::auto_ptr<TFormula> turnOnWidth_;
    std::auto_ptr<TFormula> gammaShape_;
    std::auto_ptr<TFormula> gammaScale_;
    std::auto_ptr<TFormula> overallNorm_;
    double evaluate(double tauPt, double M) const;
  };
  
  Parameterization leg1Likelihood_;
  Parameterization leg2Likelihood_;
};

#endif
