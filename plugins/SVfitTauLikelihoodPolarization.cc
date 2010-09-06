#include "TauAnalysis/CandidateTools/plugins/SVfitTauLikelihoodPolarization.h"

#include "TauAnalysis/CandidateTools/interface/svFitAuxFunctions.h"

#include <TMath.h>

#include <limits>

using namespace SVfit_namespace;

SVfitTauLikelihoodPolarization::SVfitTauLikelihoodPolarization(const edm::ParameterSet& cfg)
  : SVfitLegLikelihoodPolarizationBase<pat::Tau>(cfg)
{
  useCollApproxFormulas_ = cfg.getParameter<bool>("useCollApproxFormulas");
}

SVfitTauLikelihoodPolarization::~SVfitTauLikelihoodPolarization()
{
// nothing to be done yet...
}

double SVfitTauLikelihoodPolarization::negLogLikelihoodPolarized(
  const pat::Tau& leg, const SVfitLegSolution& solution, double tauLeptonPol) const
{
//--- compute negative log-likelihood for tau lepton decay "leg"
//    to be compatible with decay tau- --> X nu of polarized tau lepton into hadrons,
//    assuming  matrix element of V-A electroweak decay
//
//    NOTE: The formulas taken from the papers
//         [1] "Tau polarization and its correlations as a probe of new physics",
//             B.K. Bullock, K. Hagiwara and A.D. Martin,
//             Nucl. Phys. B395 (1993) 499.
//         [2] "Charged Higgs boson search at the TeVatron upgrade using tau polarization",
//             S. Raychaudhuri and D.P. Roy,
//             Phys. Rev.  D52 (1995) 1556.           
//
  std::cout << "<SVfitTauLikelihoodPolarization::negLogLikelihoodPolarized>:" << std::endl;
          
  double prob = 0.;
  
  double logLikelihood = TMath::Log(prob);
  std::cout << " -logLikelihood = " << -logLikelihood << std::endl;
  
  return -logLikelihood;
}

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_EDM_PLUGIN(SVfitTauLikelihoodBasePluginFactory, SVfitTauLikelihoodPolarization, "SVfitTauLikelihoodPolarization");
