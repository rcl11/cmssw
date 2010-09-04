#ifndef TauAnalysis_CandidateTools_SVfitLikelihoodDiTauPt_h
#define TauAnalysis_CandidateTools_SVfitLikelihoodDiTauPt_h

/** \class SVfitLikelihoodDiTauPt
 *
 * Plugin for computing likelihood for tau+ tau- pair to have a certain transverse momentum
 * 
 * \author Christian Veelken; UC Davis
 *
 * \version $Revision: 1.3 $
 *
 * $Id: SVfitLikelihoodMEt.h,v 1.3 2010/09/02 16:39:23 veelken Exp $
 *
 */

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "TauAnalysis/CandidateTools/interface/SVfitDiTauLikelihoodBase.h"
#include "TauAnalysis/CandidateTools/interface/SVfitLegLikelihoodBase.h"

#include "AnalysisDataFormats/TauAnalysis/interface/CompositePtrCandidateT1T2MEt.h"
#include "AnalysisDataFormats/TauAnalysis/interface/SVfitDiTauSolution.h"

#include <TFormula.h>

#include <iostream>

template <typename T1, typename T2>
class SVfitLikelihoodDiTauPt : public SVfitDiTauLikelihoodBase<T1,T2>
{
 public:
  SVfitLikelihoodDiTauPt(const edm::ParameterSet&);
  ~SVfitLikelihoodDiTauPt();

  double operator()(const CompositePtrCandidateT1T2MEt<T1,T2>&, const SVfitDiTauSolution&) const;
 private:
  void print(std::ostream&) const;

  TFormula* pdf_;
  std::map<int, TFormula*> pdfParameters_;
  int parseError_;
};

#endif
