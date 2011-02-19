#ifndef TauAnalysis_CandidateTools_SVfitLikelihoodWtauNuPtPtBalance_h
#define TauAnalysis_CandidateTools_SVfitLikelihoodWtauNuPtBalance_h

/** \class SVfitLikelihoodWtauNuPtBalance
 *
 * Plugin for computing likelihood for tau leptons and neutrino produced in W boson decay
 * to balance each other in transverse momentum
 *
 * \author Christian Veelken, UC Davis
 *
 * \version $Revision: 1.5 $
 *
 * $Id: SVfitLikelihoodWtauNuPtBalance.h,v 1.5 2011/01/18 16:47:16 friis Exp $
 *
 */

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "TauAnalysis/CandidateTools/interface/SVfitWtauNuLikelihoodBase.h"
#include "TauAnalysis/CandidateTools/interface/SVfitLegLikelihoodBase.h"

#include "AnalysisDataFormats/TauAnalysis/interface/CompositePtrCandidateTMEt.h"
#include "AnalysisDataFormats/TauAnalysis/interface/SVfitWtauNuSolution.h"

template <typename T>
class SVfitLikelihoodWtauNuPtBalance : public SVfitWtauNuLikelihoodBase<T>
{
 public:
  SVfitLikelihoodWtauNuPtBalance(const edm::ParameterSet&);
  ~SVfitLikelihoodWtauNuPtBalance();

  double operator()(const CompositePtrCandidateTMEt<T>&, const SVfitWtauNuSolution&) const;
};

#endif
