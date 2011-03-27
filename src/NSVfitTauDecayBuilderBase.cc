#include "TauAnalysis/CandidateTools/interface/NSVfitTauDecayBuilderBase.h"
#include "AnalysisDataFormats/TauAnalysis/interface/NSVfitTauDecayHypothesis.h"
#include "TauAnalysis/CandidateTools/interface/NSVfitAlgorithmBase.h"
#include "TauAnalysis/CandidateTools/interface/svFitAuxFunctions.h"
#include "TauAnalysis/CandidateTools/interface/nSVfitParameter.h"

using namespace nSVfit_namespace;

// Map the fit parameters to indices.
void
NSVfitTauDecayBuilderBase::beginJob(NSVfitAlgorithmBase* algorithm) {
  algorithm_ = algorithm;
  idxFitParameter_visEnFracX_ = getFitParameterIdx(algorithm, prodParticleLabel_, nSVfit_namespace::kTau_visEnFracX);
  idxFitParameter_phi_lab_    = getFitParameterIdx(algorithm, prodParticleLabel_, nSVfit_namespace::kTau_phi_lab);
  // Nu mass is only a fit parameter if the nus can have mass.
  if (!nuSystemIsMassless()) {
    idxFitParameter_nuInvMass_  = getFitParameterIdx(algorithm, prodParticleLabel_, nSVfit_namespace::kTau_nuInvMass);
  }
  // The following parameters are only valid when the track likelihoods are
  // enabled.
  idxFitParameter_deltaR_ = getFitParameterIdx(algorithm, prodParticleLabel_, nSVfit_namespace::kTau_decayDistance_lab);
}

NSVfitSingleParticleHypothesisBase* NSVfitTauDecayBuilderBase::build(
    const NSVfitTauDecayBuilderBase::inputParticleMap& inputParticles) const {

  inputParticleMap::const_iterator particlePtr = inputParticles.find(prodParticleLabel_);
  assert(particlePtr != inputParticles.end());
  const reco::Candidate* reconstructedParticle = particlePtr->second.get();
  // Build our tau decay hypthosis
  NSVfitTauDecayHypothesis* hypothesis = buildSpecific(
      particlePtr->second, prodParticleLabel_, barcodeCounter_);
  ++barcodeCounter_;

  hypothesis->p3Vis_unit_ = reconstructedParticle->p4().Vect().Unit();
  hypothesis->visMass_ = reconstructedParticle->mass();
  // If this is a leptonic tau decay, we need to setup the limits on the
  // neutrino system invariant mass parameter.
  if (!nuSystemIsMassless()) {
    NSVfitAlgorithmBase::fitParameterType* fitParameter =
      algorithm_->getFitParameter(prodParticleLabel_, nSVfit_namespace::kTau_nuInvMass);
    assert(fitParameter);
    fitParameter->upperLimit_ =
      SVfit_namespace::tauLeptonMass - hypothesis->visMass_;
  }
  // Extract the associated tracks, and fit a vertex if possible.
  hypothesis->tracks_ = extractTracks(reconstructedParticle);
  return hypothesis;
}

void
NSVfitTauDecayBuilderBase::applyFitParameter(
    NSVfitSingleParticleHypothesisBase* hypothesis, double* param) const {
  using namespace SVfit_namespace;

  // Cast to the concrete tau decay hypothesis
  NSVfitTauDecayHypothesis* hypothesis_T =
    dynamic_cast<NSVfitTauDecayHypothesis*>(hypothesis);

  // Apply any decay-type specific fit parameters (like polarizaiton etc)
  applyFitParameterSpecific(hypothesis_T, param);

  double visEnFracX = param[idxFitParameter_visEnFracX_];
  double phi_lab    = param[idxFitParameter_phi_lab_];
  double pVis_lab   = hypothesis_T->p4().P();
  double enVis_lab  = hypothesis_T->p4().energy();
  double visMass    = hypothesis_T->visMass();
  double nuInvMass  = nuSystemIsMassless() ? 0.0 :
    param[idxFitParameter_nuInvMass_];

  const reco::Candidate::Vector& p3Vis_unit = hypothesis_T->p3Vis_unit();

//--- compute momentum of visible decay products in tau lepton rest frame
  double pVis_rf = pVisRestFrame(visMass, nuInvMass);

//--- decay angle in tau lepton rest frame as function of X
//    (= energy ratio of visible decay products/tau lepton energy)
  double gjAngle = SVfit_namespace::gjAngleFromX(visEnFracX, visMass, pVis_rf, enVis_lab);

//--- compute tau lepton decay angle in laboratory frame
  double angleVis_lab = SVfit_namespace::gjAngleToLabFrame(pVis_rf, gjAngle, pVis_lab);

//--- compute tau lepton momentum in laboratory frame
  double pTau_lab = SVfit_namespace::tauMomentumLabFrame(visMass, pVis_rf, gjAngle, pVis_lab);

//--- compute tau lepton direction in laboratory frame
  reco::Candidate::Vector tauFlight = SVfit_namespace::tauDirection(p3Vis_unit, angleVis_lab, phi_lab);

//--- compute tau lepton four-vector in laboratory frame
  reco::Candidate::LorentzVector p4Tau = SVfit_namespace::tauP4(tauFlight.Unit(), pTau_lab);

  hypothesis_T->p4_fitted_      = p4Tau;
  hypothesis_T->dp4_            = (p4Tau - hypothesis_T->p4_);

  hypothesis_T->p4invis_rf_     = boostToCOM(p4Tau, hypothesis_T->dp4_);
  hypothesis_T->p4vis_rf_       = boostToCOM(p4Tau, hypothesis_T->p4());

  hypothesis_T->visEnFracX_     = visEnFracX;
  hypothesis_T->decay_angle_rf_ = gjAngle;

  if ( verbosity_ ) {
    std::cout << "<NSVfitTauDecayBuilder::applyFitParameter>:" << std::endl;
    std::cout << " visEnFracX = " << param[idxFitParameter_visEnFracX_] << std::endl;
    std::cout << " phi_lab = " << param[idxFitParameter_phi_lab_] << std::endl;
    std::cout << " enVis_lab = " << enVis_lab << std::endl;
    std::cout << " visMass = " << visMass << std::endl;
    std::cout << " nuInvMass = " << param[idxFitParameter_nuInvMass_] << std::endl;
    std::cout << " gjAngle = " << gjAngle << std::endl;
    std::cout << " angleVis_lab = " << angleVis_lab << std::endl;
    std::cout << " pTau_lab = " << pTau_lab << std::endl;
    std::cout << "p4Vis: E = " << hypothesis_T->p4_.energy() << ","
	      << " px = " << hypothesis_T->p4_.px() << ", py = " << hypothesis_T->p4_.py() << ","
	      << " pz = " << hypothesis_T->p4_.pz() << std::endl;
    std::cout << "p4Tau: E = " << p4Tau.energy() << ","
	      << " px = " << p4Tau.px() << ", py = " << p4Tau.py() << ","
	      << " pz = " << p4Tau.pz() << std::endl;
  }
}

void NSVfitTauDecayBuilderBase::print(std::ostream& stream) const {
  stream << "<NSVfitTauDecayBuilder::print>:" << std::endl;
  stream << " pluginName = " << pluginName_ << std::endl;
  stream << " pluginType = " << pluginType_ << std::endl;
  stream << " prodParticleLabel = " << prodParticleLabel_ << std::endl;
  stream << " idxFitParameter_visEnFracX = " << idxFitParameter_visEnFracX_ << std::endl;
  stream << " idxFitParameter_phi_lab = " << idxFitParameter_phi_lab_ << std::endl;
  stream << " idxFitParameter_nuInvMass = " << idxFitParameter_nuInvMass_ << std::endl;
}
