#include "TauAnalysis/CandidateTools/plugins/NSVfitTauToLepLikelihoodPolarization.h"

#include "TauAnalysis/CandidateTools/interface/NSVfitAlgorithmBase.h"
#include "TauAnalysis/CandidateTools/interface/svFitAuxFunctions.h"

#include "AnalysisDataFormats/TauAnalysis/interface/NSVfitTauToLepHypothesis.h"

#include <TMath.h>

#include <limits>

using namespace SVfit_namespace;

template <typename T>
NSVfitTauToLepLikelihoodPolarization<T>::NSVfitTauToLepLikelihoodPolarization(const edm::ParameterSet& cfg)
  : NSVfitSingleParticleLikelihood(cfg)
{
  if ( this->verbosity_ ) std::cout << "<NSVfitTauToLepLikelihoodPolarization::NSVfitTauToLepLikelihoodPolarization>:" << std::endl;

  useCollApproxFormulas_ = cfg.exists("useCollApproxFormulas") ?
    cfg.getParameter<bool>("useCollApproxFormulas") : false;
}

template <typename T>
NSVfitTauToLepLikelihoodPolarization<T>::~NSVfitTauToLepLikelihoodPolarization()
{
// nothing to be done yet...
}

template <typename T>
void NSVfitTauToLepLikelihoodPolarization<T>::beginJob(NSVfitAlgorithmBase* algorithm)
{
  algorithm->requestFitParameter(prodParticleLabel_, kTau_visEnFracX, pluginName_);
  algorithm->requestFitParameter(prodParticleLabel_, kTau_phi_lab,    pluginName_);
  algorithm->requestFitParameter(prodParticleLabel_, kTau_nuInvMass,  pluginName_);
  algorithm->requestFitParameter(prodParticleLabel_, kTau_pol,        pluginName_);
}

template <typename T>
double NSVfitTauToLepLikelihoodPolarization<T>::operator()(const NSVfitSingleParticleHypothesisBase* hypothesis) const
{
//--- compute negative log-likelihood for tau lepton decay 
//    tau- --> e- nu nu (tau- --> mu- nu nu)
//    to be compatible with matrix element of V-A electroweak decay
//
//    NOTE: The formulas taken from the paper
//           "Tau polarization and its correlations as a probe of new physics",
//           B.K. Bullock, K. Hagiwara and A.D. Martin,
//           Nucl. Phys. B395 (1993) 499.
//
  const NSVfitTauToLepHypothesis<T>* hypothesis_T = dynamic_cast<const NSVfitTauToLepHypothesis<T>*>(hypothesis);
  assert(hypothesis_T != 0);

  if ( this->verbosity_ ) std::cout << "<NSVfitTauToLepLikelihoodPolarization::operator()>:" << std::endl;
  
  double prob = 0.;
  if ( !useCollApproxFormulas_ ) {
    double chargedLepMass2 = square(hypothesis_T->p4().mass());              // electron/muon mass
    double Emax = (tauLeptonMass2 + chargedLepMass2)/(2*tauLeptonMass);      // formula (2.6)
    double E = hypothesis_T->p4vis_rf().energy();                            // electron/muon energy (in tau lepton rest-frame)
    double p = hypothesis_T->p4vis_rf().P();                                 // electron/muon momentum (in tau lepton rest-frame)
    double theta = hypothesis_T->decay_angle_rf();
    double cosTheta = TMath::Cos(theta);
    double sinTheta = TMath::Sin(theta);
    double nuMass = hypothesis_T->p4invis_rf().mass();
    double tauLeptonPol = hypothesis_T->polarization();

    if ( this->verbosity_ ) {
      std::cout << " chargedLepMass2 = " << chargedLepMass2 << std::endl;
      std::cout << " Emax = " << Emax << std::endl;
      std::cout << " E = " << E << std::endl;
      std::cout << " p = " << p << std::endl;
      std::cout << " theta = " << theta << std::endl;
      std::cout << " nuMass = " << nuMass << std::endl;
    }

    prob = p*E*(3*Emax - 2*E - chargedLepMass2/E
               + tauLeptonPol*cosTheta*(p/E)*(Emax - 2*E + chargedLepMass2/tauLeptonMass))
          *sinTheta*(nuMass/tauLeptonMass);                                  // formula (2.5)
  } else {
    double z = hypothesis_T->visEnFracX();                                   // tau lepton visible momentum fraction (in laboratory frame)
    double z2 = square(z);
    double tauLeptonPol = hypothesis_T->polarization();
    prob = (1./3.)*(1 - z)*((5 + 5*z - 4*z2) + tauLeptonPol*(1 + z - 8*z2)); // formula (2.8)
  }

  double nll = 0.;
  if ( prob > 0. ) {
    nll = -TMath::Log(prob);
  } else {
    if ( prob < 0. ) 
      edm::LogWarning ("NSVfitTauToLepLikelihoodPolarization::operator()")
	<< " Unphysical solution: prob = " << prob << " --> returning very large negative number !!";
    nll = std::numeric_limits<float>::max();
  }

  if ( this->verbosity_ ) std::cout << "--> nll = " << nll << std::endl;

  return nll;
}

#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Muon.h"

typedef NSVfitTauToLepLikelihoodPolarization<pat::Electron> NSVfitTauToElecLikelihoodPolarization;
typedef NSVfitTauToLepLikelihoodPolarization<pat::Muon> NSVfitTauToMuLikelihoodPolarization;

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_EDM_PLUGIN(NSVfitSingleParticleLikelihoodPluginFactory, NSVfitTauToElecLikelihoodPolarization, "NSVfitTauToElecLikelihoodPolarization");
DEFINE_EDM_PLUGIN(NSVfitSingleParticleLikelihoodPluginFactory, NSVfitTauToMuLikelihoodPolarization, "NSVfitTauToMuLikelihoodPolarization");

