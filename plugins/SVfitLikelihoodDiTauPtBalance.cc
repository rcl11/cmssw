#include "TauAnalysis/CandidateTools/plugins/SVfitLikelihoodDiTauPtBalance.h"

#include "TauAnalysis/CandidateTools/interface/svFitAuxFunctions.h"

#include "DataFormats/Candidate/interface/Candidate.h"

#include <TMath.h>

#include <string>

using namespace SVfit_namespace;

template <typename T1, typename T2>
SVfitLikelihoodDiTauPtBalance<T1,T2>::SVfitLikelihoodDiTauPtBalance(const edm::ParameterSet& cfg)
  : SVfitDiTauLikelihoodBase<T1,T2>(cfg)
{
// nothing to be done yet...
}

template <typename T1, typename T2>
SVfitLikelihoodDiTauPtBalance<T1,T2>::~SVfitLikelihoodDiTauPtBalance()
{
// nothing to be done yet...
}

double smearedKinematicDistribution(double x, double M, double s) 
{
//-- compute likelihood for tau lepton produced in decay of particle of mass M
//   (smeared by Gaussian distribution of width s)
//   to have transverse momentum x

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

double movingTauLeptonPtPDF(double tauPt, double diTauMass)
{
  double smearNorm = 0.52 + 0.000658*diTauMass;
  double smearWidth = 1.8 + 0.018*diTauMass;
  double M = 2.3 + 1.04*diTauMass;
  double gammaScale = 6.74 + 0.020*diTauMass;
  double gammaShape = 2.2 + 0.0364*diTauMass;
  
  return smearNorm*smearedKinematicDistribution(tauPt, M, smearWidth) +
    (1 - smearNorm)*TMath::GammaDist(tauPt, gammaShape, 0, gammaScale);
}

template <typename T1, typename T2>
double SVfitLikelihoodDiTauPtBalance<T1,T2>::operator()(const CompositePtrCandidateT1T2MEt<T1,T2>& diTau, 
					     const SVfitDiTauSolution& solution) const
{
//--- compute negative log-likelihood for two tau leptons 
//    to balance each other in transverse momentum
//
//   NOTE: the form of the likelihood function has been determined numerically
//         by fitting the tau lepton transverse momentum distribution in generated
//         Z --> tau+ tau- and H/A --> tau+ tau- decays
//         with a Taylor series approximating a Jacobian peak smeared by a Gaussian
//         plus a gamma distribution
//
  reco::Candidate::LorentzVector leg1P4 = solution.leg1().p4Vis() + solution.leg1().p4Invis();
  reco::Candidate::LorentzVector leg2P4 = solution.leg2().p4Vis() + solution.leg2().p4Invis();
  
  double diTauMass = (leg1P4 + leg2P4).mass();

  return -(TMath::Log(movingTauLeptonPtPDF(diTauMass, leg1P4.pt())) + TMath::Log(movingTauLeptonPtPDF(diTauMass, leg2P4.pt())));
}

#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/Candidate/interface/Candidate.h"

typedef SVfitLikelihoodDiTauPtBalance<pat::Electron, pat::Tau> SVfitLikelihoodElecTauPairPtBalance;
typedef SVfitLikelihoodDiTauPtBalance<pat::Muon, pat::Tau> SVfitLikelihoodMuTauPairPtBalance;
typedef SVfitLikelihoodDiTauPtBalance<pat::Tau, pat::Tau> SVfitLikelihoodDiTauPairPtBalance;
typedef SVfitLikelihoodDiTauPtBalance<pat::Electron, pat::Muon> SVfitLikelihoodElecMuPairPtBalance;
typedef SVfitLikelihoodDiTauPtBalance<reco::Candidate, reco::Candidate> SVfitLikelihoodDiCandidatePairPtBalance;

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_EDM_PLUGIN(SVfitElecTauPairLikelihoodBasePluginFactory, SVfitLikelihoodElecTauPairPtBalance, "SVfitLikelihoodElecTauPairPtBalance");
DEFINE_EDM_PLUGIN(SVfitMuTauPairLikelihoodBasePluginFactory, SVfitLikelihoodMuTauPairPtBalance, "SVfitLikelihoodMuTauPairPtBalance");
DEFINE_EDM_PLUGIN(SVfitDiTauPairLikelihoodBasePluginFactory, SVfitLikelihoodDiTauPairPtBalance, "SVfitLikelihoodDiTauPairPtBalance");
DEFINE_EDM_PLUGIN(SVfitElecMuPairLikelihoodBasePluginFactory, SVfitLikelihoodElecMuPairPtBalance, "SVfitLikelihoodElecMuPairPtBalance");
DEFINE_EDM_PLUGIN(SVfitDiCandidatePairLikelihoodBasePluginFactory, SVfitLikelihoodDiCandidatePairPtBalance, "SVfitLikelihoodDiCandidatePairPtBalance");

