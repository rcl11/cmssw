#include "TauAnalysis/CandidateTools/plugins/SVfitLikelihoodMEt.h"

#include "TauAnalysis/CandidateTools/interface/SVfitAlgorithm.h"
#include "TauAnalysis/CandidateTools/interface/svFitAuxFunctions.h"

#include "DataFormats/Candidate/interface/Candidate.h"

#include <TMath.h>

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
bool SVfitLikelihoodMEt<T1,T2>::isFittedParameter(int index) const
{
  if      ( index == SVfitAlgorithm<T1,T2>::kLeg1thetaRest ) return true;
  else if ( index == SVfitAlgorithm<T1,T2>::kLeg1phiLab    ) return true;
  else if ( index == SVfitAlgorithm<T1,T2>::kLeg2thetaRest ) return true;
  else if ( index == SVfitAlgorithm<T1,T2>::kLeg2phiLab    ) return true;
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
  //std::cout << "SVfitLikelihoodMEt::operator()>:" << std::endl;
  //std::cout << " sumEt = " << diTau.met()->sumEt() << std::endl;
  //std::cout << " metPx = " << diTau.met()->px() << std::endl;
  //std::cout << " metPy = " << diTau.met()->py() << std::endl;

  double parX = diTau.met()->sumEt();
  double perpX = diTau.met()->sumEt();

  double parSigma = parSigma_->Eval(parX);
  double parBias = parBias_->Eval(parX);
  //std::cout << " parSigma = " << parSigma << ", parBias = " << parBias << std::endl;

  double perpSigma = perpSigma_->Eval(perpX);
  double perpBias = perpBias_->Eval(perpX);
  //std::cout << " perpSigma = " << perpSigma << ", perpBias = " << perpBias << std::endl;
/*
  double parSigma = 2.85 + 0.02072*diTau.met()->sumEt();  
  double perpSigma = 2.3 + 0.02284*diTau.met()->sumEt();

  double parBias = 1.183; // ET sum reconstructed in (leptonic) leg1 direction overestimated
  double perpBias = 0.0;
*/
  reco::Candidate::LorentzVector projDirection = diTau.leg1()->p4();
  double projCosPhi = TMath::Cos(projDirection.phi());
  double projSinPhi = TMath::Sin(projDirection.phi());

  double metPx = diTau.met()->px();
  double metPy = diTau.met()->py();

  double recoMET_par = (metPx*projCosPhi + metPy*projSinPhi);
  //std::cout << " recoMET_par = " << recoMET_par << std::endl;
  double recoMET_perp = (metPx*projSinPhi - metPy*projCosPhi);
  //std::cout << " recoMET_perp = " << recoMET_perp << std::endl;
  
  reco::Candidate::LorentzVector nuP4 = solution.leg1().p4Invis() + solution.leg2().p4Invis();
  double nuPx = nuP4.px();
  double nuPy = nuP4.py();

  double fittedMET_par = (nuPx*projCosPhi + nuPy*projSinPhi);
  //std::cout << " fittedMET_par = " << fittedMET_par << std::endl;
  double fittedMET_perp = (nuPx*projSinPhi - nuPy*projCosPhi);
  //std::cout << " fittedMET_perp = " << fittedMET_perp << std::endl;
  
  double parResidual = recoMET_par - fittedMET_par - parBias;
  //std::cout << " parResidual = " << parResidual << std::endl;
  double perpResidual = recoMET_perp - fittedMET_perp - perpBias;
  //std::cout << " perpResidual = " << perpResidual << std::endl;

  double negLogLikelihood = -(logGaussian(parResidual, parSigma) + logGaussian(perpResidual, perpSigma));
  //std::cout << "--> negLogLikelihood = " << negLogLikelihood << std::endl;

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
