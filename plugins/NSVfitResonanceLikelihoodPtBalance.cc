#include "TauAnalysis/CandidateTools/plugins/NSVfitResonanceLikelihoodPtBalance.h"

#include "FWCore/Utilities/interface/Exception.h"

#include "TauAnalysis/CandidateTools/interface/NSVfitAlgorithmBase.h"
#include "TauAnalysis/CandidateTools/interface/svFitAuxFunctions.h"

#include <TMath.h>

#include <string>
#include <limits>

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

NSVfitResonanceLikelihoodPtBalance::Parameterization::Parameterization(const edm::ParameterSet& cfg) 
{
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

double
NSVfitResonanceLikelihoodPtBalance::Parameterization::evaluate(double x, double M) const 
{
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

NSVfitResonanceLikelihoodPtBalance::NSVfitResonanceLikelihoodPtBalance(const edm::ParameterSet& cfg)
  : NSVfitResonanceLikelihood(cfg),
    leg1Likelihood_(cfg.getParameter<edm::ParameterSet>("leg1")),
    leg2Likelihood_(cfg.getParameter<edm::ParameterSet>("leg2")) 
{
  power_ = ( cfg.exists("power") ) ?
    cfg.getParameter<double>("power") : 1.0;
}

NSVfitResonanceLikelihoodPtBalance::~NSVfitResonanceLikelihoodPtBalance()
{
// nothing to be done yet...
}

void NSVfitResonanceLikelihoodPtBalance::beginJob(NSVfitAlgorithmBase* algorithm) const 
{
  algorithm->requestFitParameter("allTauDecays", nSVfit_namespace::kTau_visEnFracX, pluginName_);
  algorithm->requestFitParameter("allTauDecays", nSVfit_namespace::kTau_phi_lab,    pluginName_);
  algorithm->requestFitParameter("allLeptons",   nSVfit_namespace::kLep_shiftEn,    pluginName_);
  algorithm->requestFitParameter("allNeutrinos", nSVfit_namespace::kNu_energy_lab,  pluginName_);
  algorithm->requestFitParameter("allNeutrinos", nSVfit_namespace::kNu_phi_lab,     pluginName_);
}

double 
NSVfitResonanceLikelihoodPtBalance::operator()(const NSVfitResonanceHypothesis* hypothesis) const
{
//--- compute negative log-likelihood for two tau leptons
//    to have transverse momenta leg1Pt, leg2Pt

  if ( this->verbosity_ ) std::cout << "<NSVfitLikelihoodDiTauPtBalance::operator()>:" << std::endl;

  if ( hypothesis->numDaughters() != 2 ) {
    throw cms::Exception("NSVfitResonanceLikelihoodPtBalance::operator()")
      << " Resonance hypothesis passed as function argument has " << hypothesis->numDaughters()
      << " daughter particles, exactly two expected !!\n";
  }

  const NSVfitSingleParticleHypothesis* daughter1 = hypothesis->daughter(0);
  const NSVfitSingleParticleHypothesis* daughter2 = hypothesis->daughter(1);

  double diTauMass = hypothesis->p4_fitted().mass();
  if ( this->verbosity_ ) std::cout << " diTauMass = " << diTauMass << std::endl;
  
  double leg1Pt = daughter1->p4_fitted().pt();
  double leg2Pt = daughter2->p4_fitted().pt();
  if ( this->verbosity_ ) {
    std::cout << " leg1Pt = " << leg1Pt << std::endl;
    std::cout << " leg2Pt = " << leg2Pt << std::endl;
  }

  double leg1Prob = leg1Likelihood_.evaluate(leg1Pt, diTauMass);
  double leg2Prob = leg2Likelihood_.evaluate(leg2Pt, diTauMass);
  if ( this->verbosity_ ) {
    std::cout << " leg1prob = " << leg1Prob << std::endl;
    std::cout << " leg2prob = " << leg2Prob << std::endl;
  }

  if ( !(leg1Prob > 0. && leg2Prob > 0.) ) {
    //edm::LogWarning ("SVfitLikelihoodDiTauPtBalance2::operator()")
    //  << " Unphysical solution --> returning very large number !!"
    //  << " leg1Prob: " << leg1Prob << " leg2Prob: " << leg2Prob;
    return std::numeric_limits<float>::max();
  }

  double nll = -(TMath::Log(leg1Prob) + TMath::Log(leg2Prob));
  if ( this->verbosity_ ) std::cout << "--> nll = " << nll << std::endl;

  return power_*nll;
}

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_EDM_PLUGIN(NSVfitResonanceLikelihoodPluginFactory, NSVfitResonanceLikelihoodPtBalance, "NSVfitResonanceLikelihoodPtBalance");
