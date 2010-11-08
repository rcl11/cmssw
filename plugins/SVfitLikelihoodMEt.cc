#include "TauAnalysis/CandidateTools/plugins/SVfitLikelihoodMEt.h"

#include "DataFormats/Candidate/interface/Candidate.h"

#include "TauAnalysis/CandidateTools/interface/SVfitAlgorithm.h"
#include "TauAnalysis/CandidateTools/interface/svFitAuxFunctions.h"
#include "TauAnalysis/CandidateTools/interface/candidateAuxFunctions.h"

#include <TMath.h>
#include <TVector2.h>

#include <string>

using namespace SVfit_namespace;

template <typename T1, typename T2>
SVfitLikelihoodMEt<T1,T2>::SVfitLikelihoodMEt(const edm::ParameterSet& cfg)
  : SVfitDiTauLikelihoodBase<T1,T2>(cfg),
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
}

template <typename T1, typename T2>
SVfitLikelihoodMEt<T1,T2>::~SVfitLikelihoodMEt()
{
  delete parSigma_;
  delete parBias_;

  delete perpSigma_;
  delete perpBias_;
}

template <typename T1, typename T2>
void SVfitLikelihoodMEt<T1,T2>::beginEvent(const edm::Event& evt, const edm::EventSetup& es)
{
  //std::cout << "<SVfitLikelihoodMEt::beginEvent>:" << std::endl;

  evt.getByLabel(srcPFCandidates_, pfCandidates_);
}

template <typename T1, typename T2>
bool SVfitLikelihoodMEt<T1,T2>::isFittedParameter(int index) const
{
  if      ( index == SVfit_namespace::kLeg1thetaRest ) return true;
  //else if ( index == SVfit_namespace::kLeg1phiLab    ) return true;
  else if ( index == SVfit_namespace::kLeg2thetaRest ) return true;
  //else if ( index == SVfit_namespace::kLeg2phiLab    ) return true;
  else return false;
}

template <typename T1, typename T2>
double SVfitLikelihoodMEt<T1,T2>::operator()(const CompositePtrCandidateT1T2MEt<T1,T2>& diTau,
					     const SVfitDiTauSolution& solution) const
{
//--- compute negative log-likelihood for neutrinos produced in tau lepton decays
//    to match missing transverse momentum reconstructed in the event
//
//    NB: MET likelihood is split into perp/par components along (leptonic) leg1 of the diTau object
//
  if ( verbosity_ ) {
    std::cout << "SVfitLikelihoodMEt::operator()>:" << std::endl;
    std::cout << " sumEt = " << diTau.met()->sumEt() << std::endl;
  }

//--- make unit vector bisecting tau lepton "legs"
//    and project difference between "true" generated and reconstructed MET
//    in direction parallel and perpendicular to that vector
  TVector2 diTauDirection = getDiTauBisectorDirection(diTau.leg1()->p4(), diTau.leg2()->p4());

  double dummy, metSumP_par, metSumP_perp;
  computeMEtProjection(*pfCandidates_, diTauDirection, dummy, metSumP_par, metSumP_perp);
  if ( verbosity_ ) {
    std::cout << " metSumP_par = " << metSumP_par << std::endl;
    std::cout << " metSumP_perp = " << metSumP_perp << std::endl;
  }

  double parSigma = parSigma_->Eval(metSumP_par);
  double parBias = parBias_->Eval(metSumP_par);
  if ( verbosity_ ) std::cout << " parSigma = " << parSigma << ", parBias = " << parBias << std::endl;

  double perpSigma = perpSigma_->Eval(metSumP_perp);
  double perpBias = perpBias_->Eval(metSumP_perp);
  if ( verbosity_ ) std::cout << " perpSigma = " << perpSigma << ", perpBias = " << perpBias << std::endl;

  double projCosPhi = diTauDirection.X();
  double projSinPhi = diTauDirection.Y();

  double metPx = diTau.met()->px();
  double metPy = diTau.met()->py();

  double recoMET_par = (metPx*projCosPhi + metPy*projSinPhi);
  double recoMET_perp = (metPx*projSinPhi - metPy*projCosPhi);
  if ( verbosity_ ) {
    std::cout << " recoMET_par = " << recoMET_par << std::endl;
    std::cout << " recoMET_perp = " << recoMET_perp << std::endl;
  }

  reco::Candidate::LorentzVector nuP4 = solution.leg1().p4Invis() + solution.leg2().p4Invis();
  double nuPx = nuP4.px();
  double nuPy = nuP4.py();

  double fittedMET_par = (nuPx*projCosPhi + nuPy*projSinPhi);
  double fittedMET_perp = (nuPx*projSinPhi - nuPy*projCosPhi);
  if ( verbosity_ ) {
    std::cout << " fittedMET_par = " << fittedMET_par << std::endl;
    std::cout << " fittedMET_perp = " << fittedMET_perp << std::endl;
  }

  double parResidual = recoMET_par - fittedMET_par - parBias;
  double perpResidual = recoMET_perp - fittedMET_perp - perpBias;
  if ( verbosity_ ) {
    std::cout << " parResidual = " << parResidual << std::endl;
    std::cout << " perpResidual = " << perpResidual << std::endl;
  }

  double negLogLikelihood = -(logGaussian(parResidual, parSigma) + logGaussian(perpResidual, perpSigma));
  if ( verbosity_ ) std::cout << "--> negLogLikelihood = " << negLogLikelihood << std::endl;

  return negLogLikelihood;
}

#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/Candidate/interface/Candidate.h"

typedef SVfitLikelihoodMEt<pat::Electron, pat::Tau> SVfitLikelihoodMEtElecTau;
typedef SVfitLikelihoodMEt<pat::Muon, pat::Tau> SVfitLikelihoodMEtMuTau;
typedef SVfitLikelihoodMEt<pat::Tau, pat::Tau> SVfitLikelihoodMEtDiTau;
typedef SVfitLikelihoodMEt<pat::Electron, pat::Muon> SVfitLikelihoodMEtElecMu;
typedef SVfitLikelihoodMEt<reco::Candidate, reco::Candidate> SVfitLikelihoodMEtDiCandidate;

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_EDM_PLUGIN(SVfitElecTauPairLikelihoodBasePluginFactory, SVfitLikelihoodMEtElecTau, "SVfitLikelihoodMEtElecTau");
DEFINE_EDM_PLUGIN(SVfitMuTauPairLikelihoodBasePluginFactory, SVfitLikelihoodMEtMuTau, "SVfitLikelihoodMEtMuTau");
DEFINE_EDM_PLUGIN(SVfitDiTauPairLikelihoodBasePluginFactory, SVfitLikelihoodMEtDiTau, "SVfitLikelihoodMEtDiTau");
DEFINE_EDM_PLUGIN(SVfitElecMuPairLikelihoodBasePluginFactory, SVfitLikelihoodMEtElecMu, "SVfitLikelihoodMEtElecMu");
DEFINE_EDM_PLUGIN(SVfitDiCandidatePairLikelihoodBasePluginFactory, SVfitLikelihoodMEtDiCandidate, "SVfitLikelihoodMEtDiCandidate");
