#include "TauAnalysis/CandidateTools/plugins/SVfitLikelihoodWtauNuPtBalance.h"

#include "TauAnalysis/CandidateTools/interface/svFitAuxFunctions.h"

#include "DataFormats/Candidate/interface/Candidate.h"

#include <TMath.h>

#include <string>

using namespace SVfit_namespace;

template <typename T>
SVfitLikelihoodWtauNuPtBalance<T>::SVfitLikelihoodWtauNuPtBalance(const edm::ParameterSet& cfg)
  : SVfitWtauNuLikelihoodBase<T>(cfg)
{
// nothing to be done yet...
}

template <typename T>
SVfitLikelihoodWtauNuPtBalance<T>::~SVfitLikelihoodWtauNuPtBalance()
{
// nothing to be done yet...
}

//-------------------------------------------------------------------------------
// Compute likelihood for tau lepton and neutrino produced in decay of particle of mass M
// to have transverse momenta leg1Pt, nuPt
//-------------------------------------------------------------------------------
namespace {
  double smearedKinematicDistribution(double x, double Mt, double s)
  {
    double num_first_term = TMath::Exp(-0.5*square(x)/square(s))
                           *8*s*(fourth(Mt) + 2*square(Mt)*(2*square(s) + square(x)) + 6*(square(s) + square(x))*(8*square(s) + square(x)));

    double num_second_term = TMath::Exp(-square(Mt - 2*x)/(8*square(s)))
                            *s*(15*fourth(Mt) + 14*cube(Mt) + 48*(square(s)+square(x))*(8*square(s)+square(x))
                               + 4*square(Mt)*(20*square(s) + 7*square(x)) + 24*Mt*(7*square(s)*x + cube(x)));

    double num_third_term = 4*TMath::Sqrt(TMath::TwoPi())
                           *x*(fourth(Mt) + 6*square(Mt)*square(s) + 90*fourth(s) + 2*(square(Mt) + 30*square(s))*square(x) + 6*fourth(x))
                           *(TMath::Erf((Mt - 2*x)/(2*TMath::Sqrt2()*s)) + TMath::Erf(x/(TMath::Sqrt2()*s)));
    double num_factor = 1/(2*TMath::Sqrt(TMath::TwoPi()));
    double numerator = num_factor*(num_first_term - num_second_term + num_third_term);

    // now compute normalization factor
    double den_first_term = (2*TMath::Sqrt(1.0/TMath::PiOver2())
                           *TMath::Exp(-square(Mt)/(8*square(s)))
                           *Mt*s*(11*fourth(Mt) + 44*square(Mt)*square(s) + 240*fourth(s)));
    double den_second_term = TMath::Erf(Mt/(2*TMath::Sqrt2()*s))
                            *(11*fourth(Mt)*square(Mt) - 32*fourth(Mt)*square(s) - 96*square(Mt)*fourth(s) - 960*fourth(s)*square(s));
    double denominator = (1./16)*(den_first_term + den_second_term);

    return numerator/denominator;
  }

  double movingTauLeptonPtPDF(double daughterPt, double Mt)
  {
    double smearNorm = 0.52 + 0.000658*Mt;
    double smearWidth = 1.8 + 0.018*Mt;
    double smearedMt = 2.3 + 1.04*Mt;
    double gammaScale = 6.74 + 0.020*Mt;
    double gammaShape = 2.2 + 0.0364*Mt;

    return smearNorm*smearedKinematicDistribution(daughterPt, smearedMt, smearWidth)
          + (1 - smearNorm)*TMath::GammaDist(daughterPt, gammaShape, 0., gammaScale);
  }
}

template <typename T>
double SVfitLikelihoodWtauNuPtBalance<T>::operator()(const CompositePtrCandidateTMEt<T>& candidate,
						     const SVfitWtauNuSolution& solution) const
{
//--- compute negative log-likelihood for tau lepton and neutrino produced in W boson decay
//    to have transverse momenta leg1Pt, nu2Pt

  if ( this->verbosity_ ) std::cout << "<SVfitLikelihoodWtauNuPtBalance::operator()>:" << std::endl;

  double mt = solution.mt();
  if ( this->verbosity_ ) std::cout << " Mt = " << mt << std::endl;

  double leg1Pt = solution.leg1().p4().pt();
  double nuPt = solution.nu().pt();
  if ( this->verbosity_ ) {
    std::cout << " leg1Pt = " << leg1Pt << std::endl;
    std::cout << " nuPt = " << nuPt << std::endl;
  }

  double prob = movingTauLeptonPtPDF(leg1Pt, mt)*movingTauLeptonPtPDF(nuPt, mt);
  if ( this->verbosity_ ) std::cout << "--> prob = " << prob << std::endl;

  if ( !(prob > 0.) ) {
    //edm::LogWarning ("SVfitLikelihoodWtauNuPtBalance::operator()")
    //  << " Unphysical solution --> returning very large negative number !!";
    return std::numeric_limits<float>::min();
  }

  double logLikelihood = TMath::Log(prob);
  if ( this->verbosity_ ) std::cout << " -logLikelihood = " << -logLikelihood << std::endl;

  return -logLikelihood;
}

#include "DataFormats/PatCandidates/interface/Tau.h"

typedef SVfitLikelihoodWtauNuPtBalance<pat::Tau> SVfitLikelihoodTauNuPairPtBalance;

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_EDM_PLUGIN(SVfitTauNuPairLikelihoodBasePluginFactory, SVfitLikelihoodTauNuPairPtBalance, "SVfitLikelihoodTauNuPairPtBalance");
