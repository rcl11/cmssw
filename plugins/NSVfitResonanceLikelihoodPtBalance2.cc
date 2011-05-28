#include "TauAnalysis/CandidateTools/plugins/NSVfitResonanceLikelihoodPtBalance2.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "TMath.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include <math.h>
#include <limits.h>

NSVfitResonanceLikelihoodPtBalance2::~NSVfitResonanceLikelihoodPtBalance2(){}

// constructor
NSVfitResonanceLikelihoodPtBalance2::NSVfitResonanceLikelihoodPtBalance2(
    const edm::ParameterSet& pset):NSVfitResonanceLikelihood(pset),
  leg1PDF_(pset.getParameter<edm::ParameterSet>("leg1PDF")),
  leg2PDF_(pset.getParameter<edm::ParameterSet>("leg2PDF")) {

  // Check if we want to set an upper limit on the mass of the resonance.
  // This can ensure that we don't wander into regions where the PDFs are poorly
  // defined.
  // In general, this should be safe, as the PDFs are roughly mass independent.
  maxMass_ = pset.exists("maxMass") ?
    pset.getParameter<double>("maxMass") : -1.0;
  minMass_ = pset.exists("minMass") ?
    pset.getParameter<double>("minMass") : -1.0;
  addMassFactor_ = pset.getParameter<bool>("addMassFactor");
}

void NSVfitResonanceLikelihoodPtBalance2::beginCandidate(
    const NSVfitResonanceHypothesis* hyp) const {
  // Reset the highest likelihood found so far
  highestLikelihood_ = -std::numeric_limits<double>::max();
}

double NSVfitResonanceLikelihoodPtBalance2::operator()(
    const NSVfitResonanceHypothesis* resonance) const {

  const edm::OwnVector<NSVfitSingleParticleHypothesisBase>& daughters =
    resonance->daughters();
  assert(daughters.size() == 2);

  // Get the visible angles
  double deltaPhi = TMath::Abs(
      TMath::Pi() -
      TMath::Abs(reco::deltaPhi(daughters[0].p4().phi(),
        daughters[1].p4().phi())));

  double original_mass = resonance->p4_fitted().mass();
  double mass = original_mass;

  double leg1Pt = daughters[0].p4_fitted().pt();
  double leg2Pt = daughters[1].p4_fitted().pt();

  // Always use the reconstructed mass to do the scaling.
  double scaledLeg1Pt = 2.0*leg1Pt/mass;
  double scaledLeg2Pt = 2.0*leg2Pt/mass;

  // Check if we want to cut off the mass when parameterizing the PDFs.
  if (maxMass_ > 0 && mass > maxMass_) {
    mass = maxMass_;
  }
  if (minMass_ > 0 && mass < minMass_) {
    mass = minMass_;
  }

  if (isnan(leg1Pt)) {
    edm::LogWarning("BadPtBalanceInputValue") << "Leg1 pt is NaN!"
      << " Returning largest NLL (" << highestLikelihood_ << ") so far." << std::endl;
    return highestLikelihood_;
  }
  if (isnan(leg2Pt)) {
    edm::LogWarning("BadPtBalanceInputValue") << "Leg2 pt is NaN!"
      << " Returning largest NLL (" << highestLikelihood_ << ") so far." << std::endl;
    return highestLikelihood_;
  }
  if (isnan(mass)) {
    edm::LogWarning("BadPtBalanceInputValue") << "Resonance mass is NaN!"
      << " Returning largest NLL (" << highestLikelihood_ << ") so far." << std::endl;
    return highestLikelihood_;
  }

  double nll = -TMath::Log(leg1PDF_(scaledLeg1Pt, mass, deltaPhi))
    -TMath::Log(leg2PDF_(scaledLeg2Pt, mass, deltaPhi));

  if (nll > highestLikelihood_)
    highestLikelihood_ = nll;

  if (addMassFactor_) {
    nll += 2*TMath::Log(original_mass);
  }

  return nll;
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_EDM_PLUGIN(NSVfitResonanceLikelihoodPluginFactory, NSVfitResonanceLikelihoodPtBalance2, "NSVfitResonanceLikelihoodPtBalance2");
