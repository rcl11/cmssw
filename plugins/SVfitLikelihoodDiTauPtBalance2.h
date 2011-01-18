#ifndef TauAnalysis_CandidateTools_SVfitDiTauPtBalance2_h
#define TauAnalysis_CandidateTools_SVfitDiTauPtBalance2_h

/** \class SVfitDiTauPtBalance2
 *
 * Plugin for computing likelihood for two tau leptons
 * to balance each other in transverse momentum
 *
 * \author Evan Friis, Christian Veelken; UC Davis
 *
 * \version $Revision: 1.1 $
 *
 * $Id: SVfitLikelihoodDiTauPtBalance2.h,v 1.1 2010/12/13 15:56:42 friis Exp $
 *
 */

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "TauAnalysis/CandidateTools/interface/SVfitDiTauLikelihoodBase.h"
#include "TauAnalysis/CandidateTools/interface/SVfitLegLikelihoodBase.h"

#include "AnalysisDataFormats/TauAnalysis/interface/CompositePtrCandidateT1T2MEt.h"
#include "AnalysisDataFormats/TauAnalysis/interface/SVfitDiTauSolution.h"

#include <TFormula.h>

template <typename T1, typename T2>
class SVfitLikelihoodDiTauPtBalance2 : public SVfitDiTauLikelihoodBase<T1,T2>
{
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

  public:
    SVfitLikelihoodDiTauPtBalance2(const edm::ParameterSet&);
    ~SVfitLikelihoodDiTauPtBalance2();

    double operator()(const CompositePtrCandidateT1T2MEt<T1,T2>&, const SVfitDiTauSolution&) const;

  private:
    Parameterization leg1PDF_;
    Parameterization leg2PDF_;
};

#endif
