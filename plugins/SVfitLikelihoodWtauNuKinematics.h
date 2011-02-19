#ifndef TauAnalysis_CandidateTools_SVfitLikelihoodWtauNuKinematics_h
#define TauAnalysis_CandidateTools_SVfitLikelihoodWtauNuKinematics_h

/** \class SVfitLikelihoodWtauNuKinematics
 *
 * Plugin for computing likelihood for decay kinematics of tau lepton produced in W boson decay;
 * different types of likelihoods (isotropic/unpolarized/polarized decay) 
 * can be specified via plugins for the tau lepton decay "leg";
 * plugin is used by SVfit algorithm
 * 
 * \author Christian Veelken, UC Davis
 *
 * \version $Revision: 1.6 $
 *
 * $Id: SVfitLikelihoodWtauNuKinematics.h,v 1.6 2010/09/24 10:18:36 veelken Exp $
 *
 */

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "TauAnalysis/CandidateTools/interface/SVfitWtauNuLikelihoodBase.h"
#include "TauAnalysis/CandidateTools/interface/SVfitLegLikelihoodBase.h"

#include "AnalysisDataFormats/TauAnalysis/interface/CompositePtrCandidateTMEt.h"
#include "AnalysisDataFormats/TauAnalysis/interface/SVfitWtauNuSolution.h"

#include <iomanip>

template <typename T>
class SVfitLikelihoodWtauNuKinematics : public SVfitWtauNuLikelihoodBase<T>
{
 public:
  SVfitLikelihoodWtauNuKinematics(const edm::ParameterSet&);
  ~SVfitLikelihoodWtauNuKinematics();

  virtual void beginJob();
  virtual void beginEvent(edm::Event&, const edm::EventSetup&);
  virtual void beginCandidate(const CompositePtrCandidateTMEt<T>&);

  void print(std::ostream&) const;

  bool isFittedParameter(int) const;
  bool supportsPolarization() const;

  double operator()(const CompositePtrCandidateTMEt<T>&, const SVfitWtauNuSolution&) const;

 private:
  SVfitLegLikelihoodBase<T>* leg1Likelihood_;
};

#endif
