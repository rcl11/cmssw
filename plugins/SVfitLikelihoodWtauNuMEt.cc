#include "TauAnalysis/CandidateTools/plugins/SVfitLikelihoodWtauNuMEt.h"

#include "DataFormats/Candidate/interface/Candidate.h"

#include "TauAnalysis/CandidateTools/interface/SVfitAlgorithm.h"
#include "TauAnalysis/CandidateTools/interface/svFitAuxFunctions.h"
#include "TauAnalysis/CandidateTools/interface/candidateAuxFunctions.h"

#include <TMath.h>
#include <TVector2.h>

#include <string>

using namespace SVfit_namespace;

const double parSigmaMin = 5.0;
const double perpSigmaMin = 5.0;

template <typename T>
SVfitLikelihoodWtauNuMEt<T>::SVfitLikelihoodWtauNuMEt(const edm::ParameterSet& cfg)
  : SVfitWtauNuLikelihoodBase<T>(cfg),
    parSigma_(0),
    parBias_(0),
    perpSigma_(0),
    perpBias_(0)
{
  edm::ParameterSet cfgResolution = cfg.getParameter<edm::ParameterSet>("resolution");

  parSigma_ = new TFormula("parSigma", cfgResolution.getParameter<std::string>("parSigma").data());
  parBias_ = new TFormula("parBias", cfgResolution.getParameter<std::string>("parBias").data());

  perpSigma_ = new TFormula("perpSigma", cfgResolution.getParameter<std::string>("perpSigma").data());
  perpBias_ = new TFormula("perpBias", cfgResolution.getParameter<std::string>("perpBias").data());

  srcPFCandidates_ = cfg.getParameter<edm::InputTag>("srcPFCandidates");

  varyPhi_ = cfg.exists("varyPhi") ? cfg.getParameter<bool>("varyPhi") : true;
}

template <typename T>
SVfitLikelihoodWtauNuMEt<T>::~SVfitLikelihoodWtauNuMEt()
{
  delete parSigma_;
  delete parBias_;

  delete perpSigma_;
  delete perpBias_;
}

template <typename T>
void SVfitLikelihoodWtauNuMEt<T>::beginEvent(const edm::Event& evt, const edm::EventSetup& es)
{
  //std::cout << "<SVfitLikelihoodWtauNuMEt::beginEvent>:" << std::endl;

  evt.getByLabel(srcPFCandidates_, pfCandidates_);
}

template <typename T>
void SVfitLikelihoodWtauNuMEt<T>::beginCandidate(const CompositePtrCandidateTMEt<T>& candidate)
{
  qX_ = candidate.visDecayProducts()->px() + candidate.met()->px();
  qY_ = candidate.visDecayProducts()->py() + candidate.met()->py();
  qT_ = TMath::Sqrt(qX_*qX_ + qY_*qY_);
}

template <typename T>
bool SVfitLikelihoodWtauNuMEt<T>::isFittedParameter(int index) const
{
  if      (             index == SVfit_namespace::kLeg1thetaRest ) return true;
  else if ( varyPhi_ && index == SVfit_namespace::kLeg1phiLab    ) return true;
  else return false;
}

template <typename T>
double SVfitLikelihoodWtauNuMEt<T>::operator()(const CompositePtrCandidateTMEt<T>& candidate,
					       const SVfitWtauNuSolution& solution) const
{
//--- compute negative log-likelihood for neutrinos produced in W plus tau lepton decays
//    to match missing transverse momentum reconstructed in the event
//
//    NB: MET likelihood is split into perp/par components along (leptonic) leg1 of the W candidate
//
  if ( this->verbosity_ ) {
    std::cout << "SVfitLikelihoodWtauNuMEt::operator()>:" << std::endl;
    std::cout << " sumEt = " << candidate.met()->sumEt() << std::endl;
  }

  double parSigma = parSigma_->Eval(qT_);
  if ( parSigma < parSigmaMin ) parSigma = parSigmaMin;
  double parBias = parBias_->Eval(qT_);
  if ( this->verbosity_ ) std::cout << " parSigma = " << parSigma << ", parBias = " << parBias << std::endl;

  double perpSigma = perpSigma_->Eval(qT_);
  if ( perpSigma < perpSigmaMin ) perpSigma = perpSigmaMin;
  double perpBias = perpBias_->Eval(qT_);
  if ( this->verbosity_ ) std::cout << " perpSigma = " << perpSigma << ", perpBias = " << perpBias << std::endl;

  double projCosPhi,  projSinPhi;
  if ( qT_ > 0. ) {
    projCosPhi = (qX_/qT_);
    projSinPhi = (qY_/qT_);
  } else {
//--- make unit vector in tau direction
//    and project difference between "true" generated and reconstructed MET
//    in direction parallel and perpendicular to that vector
    TVector2 tauDirection(TMath::Cos(candidate.visDecayProducts()->phi()), TMath::Sin(candidate.visDecayProducts()->phi()));

    projCosPhi = tauDirection.X();
    projSinPhi = tauDirection.Y();
  }

  double metPx = candidate.met()->px();
  double metPy = candidate.met()->py();

  double recoMET_par = (metPx*projCosPhi + metPy*projSinPhi);
  double recoMET_perp = (metPx*projSinPhi - metPy*projCosPhi);
  if ( this->verbosity_ ) {
    std::cout << " recoMET_par = " << recoMET_par << std::endl;
    std::cout << " recoMET_perp = " << recoMET_perp << std::endl;
  }

  reco::Candidate::LorentzVector nuP4 = solution.p4Invis();
  double nuPx = nuP4.px();
  double nuPy = nuP4.py();

  double fittedMET_par = (nuPx*projCosPhi + nuPy*projSinPhi);
  double fittedMET_perp = (nuPx*projSinPhi - nuPy*projCosPhi);
  if ( this->verbosity_ ) {
    std::cout << " fittedMET_par = " << fittedMET_par << std::endl;
    std::cout << " fittedMET_perp = " << fittedMET_perp << std::endl;
  }

  double parResidual = (recoMET_par - fittedMET_par) - parBias;
  double perpResidual = (recoMET_perp - fittedMET_perp) - perpBias;
  if ( this->verbosity_ ) {
    std::cout << " parResidual = " << parResidual << std::endl;
    std::cout << " perpResidual = " << perpResidual << std::endl;
  }

  double negLogLikelihood = -(logGaussian(parResidual, parSigma) + logGaussian(perpResidual, perpSigma));
  if ( this->verbosity_ ) std::cout << "--> negLogLikelihood = " << negLogLikelihood << std::endl;

  return negLogLikelihood;
}

#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/Candidate/interface/Candidate.h"

typedef SVfitLikelihoodWtauNuMEt<pat::Tau> SVfitLikelihoodTauNuPairMEt;

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_EDM_PLUGIN(SVfitTauNuPairLikelihoodBasePluginFactory, SVfitLikelihoodTauNuPairMEt, "SVfitLikelihoodTauNuPairMEt");
