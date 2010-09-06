#include "TauAnalysis/CandidateTools/plugins/SVfitLeptonLikelihoodPolarization.h"

#include "TauAnalysis/CandidateTools/interface/svFitAuxFunctions.h"

#include "DataFormats/Candidate/interface/Candidate.h"

#include <TMath.h>

#include <limits>

using namespace SVfit_namespace;

template <typename T>
SVfitLeptonLikelihoodPolarization<T>::SVfitLeptonLikelihoodPolarization(const edm::ParameterSet& cfg)
  : SVfitLegLikelihoodPolarizationBase<T>(cfg)
{
  useCollApproxFormulas_ = cfg.getParameter<bool>("useCollApproxFormulas");
}

template <typename T>
SVfitLeptonLikelihoodPolarization<T>::~SVfitLeptonLikelihoodPolarization()
{
// nothing to be done yet...
}

template <typename T>
double SVfitLeptonLikelihoodPolarization<T>::negLogLikelihoodPolarized(
  const T& leg, const SVfitLegSolution& solution, double tauLeptonPol) const
{
//--- compute negative log-likelihood for tau lepton decay "leg"
//    to be compatible with decay tau- --> e- nu nu (tau- --> mu- nu nu)
//    of polarized tau lepton into electron (muon),
//    assuming  matrix element of V-A electroweak decay
//
//    NOTE: The formulas taken from the paper
//           "Tau polarization and its correlations as a probe of new physics",
//           B.K. Bullock, K. Hagiwara and A.D. Martin,
//           Nucl. Phys. B395 (1993) 499.
//
  std::cout << "<SVfitLeptonLikelihoodPolarization::negLogLikelihoodPolarized>:" << std::endl;
          
  double prob = 0.;

  if ( !useCollApproxFormulas_ ) {
    double emuMass2 = square(leg.mass());                                            // electron/muon mass
    std::cout << " emuMass2 = " << emuMass2 << std::endl;
    double Emax = (tauLeptonMass2 + emuMass2)/(2*tauLeptonMass);                     // formula (2.6)    
    double E = leg.energy();                                                         // electron/muon energy
    double p = leg.p();                                                              // electron/muon momentum
    double cosTheta = solution.cosThetaRest();
    double theta = TMath::ACos(cosTheta);
    prob = p*E*(3*Emax - 2*E - emuMass2/E
               + tauLeptonPol*cosTheta*(p/E)*(Emax - 2*E + emuMass2/tauLeptonMass))
          *TMath::Sin(theta);                                                        // formula (2.5)
  } else {
    double z = solution.x();                                                         // tau lepton visible momentum fraction
    double z2 = square(z);
    prob = (1./3.)*(1 - z)*((5 + 5*z - 4*z2) + tauLeptonPol*(1 + z - 8*z2));         // formula (2.8)
  }

  if ( !(prob > 0.) ) {
    edm::LogWarning ("SVfitLeptonLikelihoodPolarization::operator()") 
      << " Unphysical solution --> returning very large negative number !!";
    return std::numeric_limits<float>::min();
  }
  
  double logLikelihood = TMath::Log(prob);
  std::cout << " -logLikelihood = " << -logLikelihood << std::endl;
  
  return -logLikelihood;
}

#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Muon.h"

typedef SVfitLeptonLikelihoodPolarization<pat::Electron> SVfitElectronLikelihoodPolarization;
typedef SVfitLeptonLikelihoodPolarization<pat::Muon> SVfitMuonLikelihoodPolarization;

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_EDM_PLUGIN(SVfitElectronLikelihoodBasePluginFactory, SVfitElectronLikelihoodPolarization, "SVfitElectronLikelihoodPolarization");
DEFINE_EDM_PLUGIN(SVfitMuonLikelihoodBasePluginFactory, SVfitMuonLikelihoodPolarization, "SVfitMuonLikelihoodPolarization");
