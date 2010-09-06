#ifndef TauAnalysis_CandidateTools_SVfitTauLikelihoodPolarization_h
#define TauAnalysis_CandidateTools_SVfitTauLikelihoodPolarization_h

/** \class SVfitTauLikelihoodPolarization
 *
 * Plugin for computing likelihood for tau lepton decay "leg"
 * to be compatible with decay tau- --> X nu of polarized tau lepton into hadrons,
 * assuming  matrix element of V-A electroweak decay
 * 
 * NOTE: the system of hadrons X can either be pi-, rho- --> pi- pi0, 
 *       a1- --> pi- pi0 pi0 or a1- --> pi- pi+ pi-;
 *       tau decays into pi- pi+ pi- pi0 are **not** supported
 *
 * \author Christian Veelken, UC Davis
 *
 * \version $Revision: 1.1 $
 *
 * $Id: SVfitTauLikelihoodPolarization.h,v 1.1 2010/09/06 07:48:16 veelken Exp $
 *
 */

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/PatCandidates/interface/Tau.h"

#include "TauAnalysis/CandidateTools/interface/SVfitLegLikelihoodPolarizationBase.h"

class SVfitTauLikelihoodPolarization : public SVfitLegLikelihoodPolarizationBase<pat::Tau>
{
 public:
  SVfitTauLikelihoodPolarization(const edm::ParameterSet&);
  ~SVfitTauLikelihoodPolarization();

 private:
  double negLogLikelihoodPolarized(const pat::Tau&, const SVfitLegSolution&, double) const;

  double probOneProngZeroPi0(const pat::Tau&, const SVfitLegSolution&, double) const;
  double probOneProngOnePi0(const pat::Tau&, const SVfitLegSolution&, double) const;
  double probOneProngTwoPi0(const pat::Tau&, const SVfitLegSolution&, double) const;
  double probThreeProngZeroPi0(const pat::Tau&, const SVfitLegSolution&, double) const;
  double probOtherDecayMode(const pat::Tau&, const SVfitLegSolution&, double) const;

  std::vector<int> supportedTauDecayModes_;
  size_t numSupportedTauDecayModes_;

  mutable TVectorD vRec_;
  TMatrixD mapRecToGenTauDecayModes_;
  mutable TVectorD vGen_;
  mutable TVectorD vProb_;

  int tauDecayModeOther_index_;

  SVfitLegLikelihoodBase<pat::Tau>* likelihoodPhaseSpace_;

  bool useCollApproxFormulas_;
};

#endif
