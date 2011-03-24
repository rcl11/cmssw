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
 * \version $Revision: 1.1 $
 *
 * $Id: NSVfitTauToHadLikelihoodPolarization.h,v 1.1 2011/03/23 17:46:39 veelken Exp $
 *
 */

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "TauAnalysis/CandidateTools/interface/NSVfitSingleParticleLikelihood.h"

#include "AnalysisDataFormats/TauAnalysis/interface/NSVfitSingleParticleHypothesisBase.h"
#include "TauAnalysis/CandidateTools/interface/NSVfitVMlineShape.h"

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

  struct decayModeEntryType
  {
    decayModeEntryType(const edm::ParameterSet& cfg)
      : vmLineShapeLpol_(0),
	vmLineShapeTpol_(0),
        xSigma_(0),
	xBias_(0)
    {
      if ( cfg.exists("xSigma") ) xSigma_ = new TFormula("xSigma", cfg.getParameter<std::string>("xSigma").data());
      if ( cfg.exists("xBias")  ) xBias_  = new TFormula("xBias",  cfg.getParameter<std::string>("xBias").data());
      pMin_ = cfg.getParameter<double>("pMin");
    }
    ~decayModeEntryType()
    {
      delete vmLineShapeLpol_;
      delete vmLineShapeTpol_;
      delete xSigma_;
      delete xBias_;
    }
    void print(std::ostream& stream) const
    {
      stream << "<decayModeEntryType::print>:" << std::endl;
      if ( xSigma_ != 0 ) std::cout << " xSigma = " << xSigma_->GetTitle() << std::endl;
      if ( xBias_  != 0 ) std::cout << " xBias = " << xBias_->GetTitle() << std::endl;
      std::cout << " pMin = " << pMin_ << std::endl;
    }
    NSVfitVMlineShape* vmLineShapeLpol_;
    NSVfitVMlineShape* vmLineShapeTpol_;
    TFormula* xSigma_;
    TFormula* xBias_;
    double pMin_;
  };

  double probChargedPionDecay(const NSVfitTauToHadHypothesis*) const;
  double probVMrhoDecay(const NSVfitTauToHadHypothesis*) const; 
  double probVMa1Decay(const NSVfitTauToHadHypothesis*, const decayModeEntryType*) const; 
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

  std::vector<decayModeEntryType*> decayModeParameters_;
  std::vector<bool> fitDecayMode_;

  NSVfitVMlineShape* rhoLpolLineShape_;
  NSVfitVMlineShape* rhoTpolLineShape_;
  NSVfitVMlineShape* a1NeutralLpolLineShape_;
  NSVfitVMlineShape* a1NeutralTpolLineShape_;
  NSVfitVMlineShape* a1ChargedLpolLineShape_;
  NSVfitVMlineShape* a1ChargedTpolLineShape_;

  NSVfitSingleParticleLikelihood* likelihoodPhaseSpace_;
  
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
