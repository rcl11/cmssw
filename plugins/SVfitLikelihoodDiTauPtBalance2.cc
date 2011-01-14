#include "TauAnalysis/CandidateTools/plugins/SVfitLikelihoodDiTauPtBalance2.h"

#include "TauAnalysis/CandidateTools/interface/svFitAuxFunctions.h"

#include "DataFormats/Candidate/interface/Candidate.h"

#include <TMath.h>

#include <string>

using namespace SVfit_namespace;

//-------------------------------------------------------------------------------
// Compute likelihood for tau leptons produced in decay of particle of mass M
// to have transverse momenta leg1Pt, leg2Pt
//
// The form of the likelihood function has been determined numerically by:
//  o making a Taylor expansion of Jacobian peak smeared by Gaussian distribution
//  o fitting the first five terms of the Taylor expansion plus a gamma distribution
//    to the tau lepton Pt distribution in simulated Z --> tau+ tau- and H/A --> tau+ tau- events
//  o parametrizing the fit coefficients as function of tau+ tau- mass
//    of the Z --> tau+ tau- and H/A --> tau+ tau- events
//
//-------------------------------------------------------------------------------

namespace {
double smearedKinematicDistribution(double x, double M, double s) {
  double num_first_term = TMath::Exp(-0.5*square(x)/square(s))
                         *8*s*(fourth(M) + 2*square(M)*(2*square(s) + square(x)) + 6*(square(s) + square(x))*(8*square(s) + square(x)));

  double num_second_term = TMath::Exp(-square(M - 2*x)/(8*square(s)))
                          *s*(15*fourth(M) + 14*cube(M) + 48*(square(s)+square(x))*(8*square(s)+square(x))
                             + 4*square(M)*(20*square(s) + 7*square(x)) + 24*M*(7*square(s)*x + cube(x)));

  double num_third_term = 4*TMath::Sqrt(TMath::TwoPi())
                         *x*(fourth(M) + 6*square(M)*square(s) + 90*fourth(s) + 2*(square(M) + 30*square(s))*square(x) + 6*fourth(x))
                         *(TMath::Erf((M - 2*x)/(2*TMath::Sqrt2()*s)) + TMath::Erf(x/(TMath::Sqrt2()*s)));
  double num_factor = 1/(2*TMath::Sqrt(TMath::TwoPi()));
  double numerator = num_factor*(num_first_term - num_second_term + num_third_term);

  // now compute normalization factor
  double den_first_term = (2*TMath::Sqrt(1.0/TMath::PiOver2())
                         *TMath::Exp(-square(M)/(8*square(s)))
                         *M*s*(11*fourth(M) + 44*square(M)*square(s) + 240*fourth(s)));
  double den_second_term = TMath::Erf(M/(2*TMath::Sqrt2()*s))
                          *(11*fourth(M)*square(M) - 32*fourth(M)*square(s) - 96*square(M)*fourth(s) - 960*fourth(s)*square(s));
  double denominator = (1./16)*(den_first_term + den_second_term);

  return numerator/denominator;
}
}

template <typename T1, typename T2>
SVfitLikelihoodDiTauPtBalance2<T1,T2>::Parameterization::Parameterization(
    const edm::ParameterSet& cfg) {
  smear_.reset(new TFormula(
          "smear", cfg.getParameter<std::string>("smear").data()));
  gaussFrac_.reset(new TFormula(
          "gaussFrac", cfg.getParameter<std::string>("gaussFrac").data()));
  turnOnWidth_.reset(new TFormula(
          "turnOnWidth", cfg.getParameter<std::string>("turnOnWidth").data()));
  turnOnThreshold_.reset(new TFormula(
          "turnOnThreshold", cfg.getParameter<std::string>("turnOnThreshold").data()));
  gammaScale_.reset(new TFormula(
          "gammaScale", cfg.getParameter<std::string>("gammaScale").data()));
  gammaShape_.reset(new TFormula(
          "gammaShape", cfg.getParameter<std::string>("gammaShape").data()));
  overallNorm_.reset(new TFormula(
          "overallNorm", cfg.getParameter<std::string>("overallNorm").data()));
}

template <typename T1, typename T2>
double
SVfitLikelihoodDiTauPtBalance2<T1,T2>::Parameterization::evaluate(
    double x, double M) const {
  // Get smeared kinematic contribution
  double smeared = smearedKinematicDistribution(x, M, smear_->Eval(M));
  // Get Gamma distribution (tails) contribution
  double gamma = TMath::GammaDist(x, gammaShape_->Eval(M), 0, gammaShape_->Eval(M));
  double gaussFrac = gaussFrac_->Eval(M);
  // Combine the two
  double additivePart = gaussFrac*smeared + (1-gaussFrac)*gamma;

  // Multiply by turnon and overal normalization factor
  double turnOn = 0.5*(1+TMath::Erf(
          turnOnWidth_->Eval(M)* (x - turnOnThreshold_->Eval(M))));
  double norm = overallNorm_->Eval(M);
  return norm*turnOn*additivePart;
}

template <typename T1, typename T2>
SVfitLikelihoodDiTauPtBalance2<T1,T2>::SVfitLikelihoodDiTauPtBalance2(const edm::ParameterSet& cfg)
  :SVfitDiTauLikelihoodBase<T1,T2>(cfg),
    leg1PDF_(cfg.getParameter<edm::ParameterSet>("leg1")),
    leg2PDF_(cfg.getParameter<edm::ParameterSet>("leg2")) {}

template <typename T1, typename T2>
SVfitLikelihoodDiTauPtBalance2<T1,T2>::~SVfitLikelihoodDiTauPtBalance2()
{
// nothing to be done yet...
}


template <typename T1, typename T2>
double SVfitLikelihoodDiTauPtBalance2<T1,T2>::operator()(const CompositePtrCandidateT1T2MEt<T1,T2>& diTau,
					     const SVfitDiTauSolution& solution) const
{
//--- compute negative log-likelihood for two tau leptons
//    to have transverse momenta leg1Pt, leg2Pt

  if ( verbosity_ ) std::cout << "<SVfitLikelihoodDiTauPtBalance2::operator()>:" << std::endl;

  double diTauMass = solution.p4().mass();
  if ( verbosity_ ) std::cout << " diTauMass = " << diTauMass << std::endl;

  double leg1Pt = solution.leg1().p4().pt();
  double leg2Pt = solution.leg2().p4().pt();
  if ( verbosity_ ) {
    std::cout << " leg1Pt = " << leg1Pt << std::endl;
    std::cout << " leg2Pt = " << leg2Pt << std::endl;
  }

  double leg1Prob = leg1PDF_.evaluate(leg1Pt, diTauMass);
  double leg2Prob = leg2PDF_.evaluate(leg2Pt, diTauMass);

  if ( verbosity_ ) std::cout << "--> leg1prob = " << leg1Prob << std::endl;
  if ( verbosity_ ) std::cout << "--> leg2prob = " << leg2Prob << std::endl;

  if ( !(leg1Prob > 0.) || !(leg2Prob > 0.) ) {
//    edm::LogWarning ("SVfitLikelihoodDiTauPtBalance2::operator()")
//      << " Unphysical solution --> returning very large negative number !!"
//      << " leg1Prob: " << leg1Prob << " leg2Prob: " << leg2Prob;
    return std::numeric_limits<float>::min();
  }

  double logLikelihood = TMath::Log(leg1Prob) + TMath::Log(leg2Prob);
  if ( verbosity_ ) std::cout << " -logLikelihood = " << -logLikelihood << std::endl;

  return -logLikelihood;
}

#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/Candidate/interface/Candidate.h"

typedef SVfitLikelihoodDiTauPtBalance2<pat::Electron, pat::Tau> SVfitLikelihoodElecTauPairPtBalance2;
typedef SVfitLikelihoodDiTauPtBalance2<pat::Muon, pat::Tau> SVfitLikelihoodMuTauPairPtBalance2;
typedef SVfitLikelihoodDiTauPtBalance2<pat::Tau, pat::Tau> SVfitLikelihoodDiTauPairPtBalance2;
typedef SVfitLikelihoodDiTauPtBalance2<pat::Electron, pat::Muon> SVfitLikelihoodElecMuPairPtBalance2;
typedef SVfitLikelihoodDiTauPtBalance2<pat::Electron, pat::Electron> SVfitLikelihoodDiElecPairPtBalance2;
typedef SVfitLikelihoodDiTauPtBalance2<pat::Muon, pat::Muon> SVfitLikelihoodDiMuPairPtBalance2;
typedef SVfitLikelihoodDiTauPtBalance2<reco::Candidate, reco::Candidate> SVfitLikelihoodDiCandidatePairPtBalance2;

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_EDM_PLUGIN(SVfitElecTauPairLikelihoodBasePluginFactory, SVfitLikelihoodElecTauPairPtBalance2, "SVfitLikelihoodElecTauPairPtBalance2");
DEFINE_EDM_PLUGIN(SVfitMuTauPairLikelihoodBasePluginFactory, SVfitLikelihoodMuTauPairPtBalance2, "SVfitLikelihoodMuTauPairPtBalance2");
DEFINE_EDM_PLUGIN(SVfitDiTauPairLikelihoodBasePluginFactory, SVfitLikelihoodDiTauPairPtBalance2, "SVfitLikelihoodDiTauPairPtBalance2");
DEFINE_EDM_PLUGIN(SVfitElecMuPairLikelihoodBasePluginFactory, SVfitLikelihoodElecMuPairPtBalance2, "SVfitLikelihoodElecMuPairPtBalance2");
DEFINE_EDM_PLUGIN(SVfitDiElecPairLikelihoodBasePluginFactory, SVfitLikelihoodDiElecPairPtBalance2, "SVfitLikelihoodDiElecPairPtBalance2");
DEFINE_EDM_PLUGIN(SVfitDiMuPairLikelihoodBasePluginFactory, SVfitLikelihoodDiMuPairPtBalance2, "SVfitLikelihoodDiMuPairPtBalance2");
DEFINE_EDM_PLUGIN(SVfitDiCandidatePairLikelihoodBasePluginFactory, SVfitLikelihoodDiCandidatePairPtBalance2, "SVfitLikelihoodDiCandidatePairPtBalance2");

