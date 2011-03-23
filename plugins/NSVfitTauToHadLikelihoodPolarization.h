#ifndef TauAnalysis_CandidateTools_NSVfitTauToHadLikelihoodPolarization_h
#define TauAnalysis_CandidateTools_NSVfitTauToHadLikelihoodPolarization_h

/** \class NSVfitTauToHadLikelihoodPolarization
 *
 * Plugin for computing likelihood for tau lepton decay 
 *  tau- --> X nu
 * to be compatible with matrix element of V-A electroweak decay
 *
 * NOTE: the system of hadrons X can either be pi-, rho- --> pi- pi0,
 *       a1- --> pi- pi0 pi0 or a1- --> pi- pi+ pi-;
 *       tau decays into pi- pi+ pi- pi0 are **not** supported
 *
 * \author Christian Veelken, UC Davis
 *
 * \version $Revision: 1.7 $
 *
 * $Id: NSVfitTauToHadLikelihoodPolarization.h,v 1.7 2011/01/18 16:47:16 friis Exp $
 *
 */

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "TauAnalysis/CandidateTools/interface/NSVfitSingleParticleLikelihood.h"

#include "AnalysisDataFormats/TauAnalysis/interface/NSVfitSingleParticleHypothesisBase.h"
#include "TauAnalysis/CandidateTools/interface/SVfitVMlineShapeIntegral.h"

#include <TFormula.h>
#include <TVectorD.h>

class NSVfitTauToHadHypothesis;

class NSVfitTauToHadLikelihoodPolarization : public NSVfitSingleParticleLikelihood
{
 public:
  NSVfitTauToHadLikelihoodPolarization(const edm::ParameterSet&);
  ~NSVfitTauToHadLikelihoodPolarization();

  void beginJob(NSVfitAlgorithmBase*);
  void beginCandidate(const NSVfitSingleParticleHypothesisBase*);
  
  double operator()(const NSVfitSingleParticleHypothesisBase*) const;

 private:
  enum decayModes { kPion, kVMrho, kVMa1Neutral, kVMa1Charged, kOther };

  double probOneProngZeroPi0s(const NSVfitTauToHadHypothesis*) const;
  double probOneProngOnePi0(const NSVfitTauToHadHypothesis*) const; 
  double probOneProngTwoPi0s(const NSVfitTauToHadHypothesis*) const;
  double probThreeProngZeroPi0s(const NSVfitTauToHadHypothesis*) const;
  double probOtherDecayMode(const NSVfitTauToHadHypothesis*) const;

//--- auxiliary functions needed for computation of likelihood
//    for tau- --> a1- nu --> pi- pi0 pi0 nu, tau- --> a1- nu --> pi- pi+ pi- nu decays
  double compVMa1x(double, double, double, double, double) const;
  double compVMa1DecayProbL(double, double, double, double, double) const;
  double compVMa1DecayProbT(double, double, double, double, double, double) const;

  std::vector<int> supportedTauDecayModes_;
  size_t numSupportedTauDecayModes_;

  mutable TVectorD vRec_;
  TMatrixD mapRecToGenTauDecayModes_;
  mutable TVectorD vGen_;
  mutable TVectorD vProb_;

  struct decayModeEntryType
  {
    decayModeEntryType(const edm::ParameterSet&);
    ~decayModeEntryType();
    void print(std::ostream&) const;
    TFormula* xSigma_;
    TFormula* xBias_;
    double pMin_;
  };

  std::vector<decayModeEntryType*> decayModeParameters_;
  std::vector<bool> fitDecayMode_;

  SVfitVMlineShapeIntegral* rhoLpolLineShape_;
  SVfitVMlineShapeIntegral* rhoTpolLineShape_;
  SVfitVMlineShapeIntegral* a1LpolLineShape_;
  SVfitVMlineShapeIntegral* a1TpolLineShape_;

  NSVfitSingleParticleLikelihood* likelihoodPhaseSpace_;

  bool useCollApproxFormulas_;

//--- temporary variables to speed-up computations
//    (computed once in constructor)
  double a1posMassTerm_;
  double a1posMassTerm2_;
  double a1negMassTerm_;
  double a1q_;
  double a1_8piDiv9_;
  double a1_16piDiv9_;
};

#endif
