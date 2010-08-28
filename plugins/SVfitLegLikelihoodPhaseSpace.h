#ifndef TauAnalysis_CandidateTools_SVfitLegLikelihoodPhaseSpace_h
#define TauAnalysis_CandidateTools_SVfitLegLikelihoodPhaseSpace_h

/** \class SVfitLegLikelihoodPhaseSpace
 *
 * Plugin for computing likelihood for tau lepton decay "leg"
 * to be compatible with three-body decay,
 * assuming constant matrix element, 
 * so that energy and angular distribution of decay products is solely determined by phase-space
 * 
 * \author Evan Friis, Christian Veelken; UC Davis
 *
 * \version $Revision: 1.1 $
 *
 * $Id: SVfitLegLikelihoodPhaseSpace.h,v 1.1 2010/08/28 13:18:03 veelken Exp $
 *
 */

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "TauAnalysis/CandidateTools/interface/SVfitLegLikelihoodBase.h"

#include "AnalysisDataFormats/TauAnalysis/interface/SVfitLegSolution.h"

template <typename T>
class SVfitLegLikelihoodPhaseSpace : public SVfitLegLikelihoodBase<T>
{
 public:
  SVfitLegLikelihoodPhaseSpace(const edm::ParameterSet&);
  ~SVfitLegLikelihoodPhaseSpace();

  double operator()(const T&, const SVfitLegSolution&) const;
};

#endif
