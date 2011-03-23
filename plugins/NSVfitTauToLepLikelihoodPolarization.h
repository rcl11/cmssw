#ifndef TauAnalysis_CandidateTools_NSVfitTauToLepLikelihoodPolarization_h
#define TauAnalysis_CandidateTools_NSVfitTauToLepLikelihoodPolarization_h

/** \class NSVfitTauToLepLikelihoodPolarization
 *
 * Plugin for computing likelihood for tau lepton decay 
 *  tau- --> e- nu nu (tau- --> mu- nu nu)
 * to be compatible with matrix element of V-A electroweak decay
 *
 * \author Christian Veelken, UC Davis
 *
 * \version $Revision: 1.3 $
 *
 * $Id: NSVfitTauToLepLikelihoodPolarization.h,v 1.3 2011/01/18 16:47:16 friis Exp $
 *
 */

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "TauAnalysis/CandidateTools/interface/NSVfitSingleParticleLikelihood.h"

#include "AnalysisDataFormats/TauAnalysis/interface/NSVfitSingleParticleHypothesisBase.h"

template <typename T>
class NSVfitTauToLepLikelihoodPolarization : public NSVfitSingleParticleLikelihood
{
 public:
  NSVfitTauToLepLikelihoodPolarization(const edm::ParameterSet&);
  ~NSVfitTauToLepLikelihoodPolarization();
  
  void beginJob(NSVfitAlgorithmBase*);
  
  double operator()(const NSVfitSingleParticleHypothesisBase*) const;

 private:
  bool useCollApproxFormulas_;
};

#endif
