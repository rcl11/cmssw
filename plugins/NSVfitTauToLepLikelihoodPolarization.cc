#include "TauAnalysis/CandidateTools/plugins/NSVfitTauToLepLikelihoodPolarization.h"

#include "TauAnalysis/CandidateTools/interface/NSVfitAlgorithmBase.h"
#include "TauAnalysis/CandidateTools/interface/svFitAuxFunctions.h"

#include "AnalysisDataFormats/TauAnalysis/interface/NSVfitTauToDaughtersHypothesisBaseT1T2.h"

#include <TMath.h>

#include <limits>

using namespace SVfit_namespace;

template <typename T>
NSVfitTauToLepLikelihoodPolarization<T>::NSVfitTauToLepLikelihoodPolarization(const edm::ParameterSet& cfg)
  : NSVfitSingleParticleLikelihood(cfg)
{
  if ( this->verbosity_ ) std::cout << "<NSVfitTauToLepLikelihoodPolarization::NSVfitTauToLepLikelihoodPolarization>:" << std::endl;
}

template <typename T>
NSVfitTauToLepLikelihoodPolarization<T>::~NSVfitTauToLepLikelihoodPolarization()
{
// nothing to be done yet...
}

template <typename T>
void NSVfitTauToLepLikelihoodPolarization<T>::beginJob(NSVfitAlgorithmBase* algorithm)
{
  algorithm->requestFitParameter(prodParticleLabel_, nSVfit_namespace::kTau_visEnFracX, pluginName_);
  algorithm->requestFitParameter(prodParticleLabel_, nSVfit_namespace::kTau_phi_lab,    pluginName_);
  algorithm->requestFitParameter(prodParticleLabel_, nSVfit_namespace::kTau_nuInvMass,  pluginName_);
  algorithm->requestFitParameter(prodParticleLabel_, nSVfit_namespace::kTau_pol,        pluginName_);
}

template <typename T>
double NSVfitTauToLepLikelihoodPolarization<T>::operator()(const NSVfitSingleParticleHypothesis* hypothesis) const
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
  const NSVfitTauToDaughtersHypothesisBaseT1T2<NSVfitTauDecayHypothesis, T>* hypothesis_T = 
    dynamic_cast<const NSVfitTauToDaughtersHypothesisBaseT1T2<NSVfitTauDecayHypothesis, T>*>(hypothesis);
  assert(hypothesis_T != 0);

  if ( this->verbosity_ ) std::cout << "<NSVfitTauToLepLikelihoodPolarization::operator()>:" << std::endl;
  
  double chargedLepMass2 = square(hypothesis_T->p4().mass());         // electron/muon mass
  double Emax = (tauLeptonMass2 + chargedLepMass2)/(2*tauLeptonMass); // formula (2.6)
  double E = hypothesis_T->p4vis_rf().energy();                       // electron/muon energy (in tau lepton rest-frame)
  double p = hypothesis_T->p4vis_rf().P();                            // electron/muon momentum (in tau lepton rest-frame)
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

  double prob = p*E*(3*Emax - 2*E - chargedLepMass2/E + tauLeptonPol*cosTheta*(p/E)*(Emax - 2*E + chargedLepMass2/tauLeptonMass))
               *sinTheta*(nuMass/tauLeptonMass); // formula (2.5)

  if ( applyVisPtCutCorrection_ ) prob *= evaluateVisPtCutCorrection(hypothesis);

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

