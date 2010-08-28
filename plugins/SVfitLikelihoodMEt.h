#ifndef TauAnalysis_CandidateTools_SVfitLikelihoodMEt_h
#define TauAnalysis_CandidateTools_SVfitLikelihoodMEt_h

/** \class SVfitLikelihoodMEt
 *
 * Plugin for computing likelihood for neutrinos produced in tau lepton decays
 * to match missing transverse momentum reconstructed in the event
 * 
 * \author Evan Friis, Christian Veelken; UC Davis
 *
 * \version $Revision: 1.2 $
 *
 * $Id: SVfitLikelihoodDiTauKinematics.h,v 1.2 2010/08/27 06:59:19 veelken Exp $
 *
 */

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "TauAnalysis/CandidateTools/interface/SVfitDiTauLikelihoodBase.h"
#include "TauAnalysis/CandidateTools/interface/SVfitLegLikelihoodBase.h"

#include "AnalysisDataFormats/TauAnalysis/interface/CompositePtrCandidateT1T2MEt.h"
#include "AnalysisDataFormats/TauAnalysis/interface/SVfitDiTauSolution.h"

#include <TFormula.h>

template <typename T1, typename T2>
class SVfitLikelihoodMEt : public SVfitDiTauLikelihoodBase<T1,T2>
{
 public:
  SVfitLikelihoodMEt(const edm::ParameterSet&);
  ~SVfitLikelihoodMEt();

  bool isFittedParameter(int) const;

  double operator()(const CompositePtrCandidateT1T2MEt<T1,T2>&, const SVfitDiTauSolution&) const;
 private:
  TFormula* parSigma_;
  TFormula* parBias_;

  TFormula* perpSigma_;
  TFormula* perpBias_;
};

#endif
