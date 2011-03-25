#ifndef TauAnalysis_CandidateTools_nSVfitParameter_h
#define TauAnalysis_CandidateTools_nSVfitParameter_h

namespace nSVfit_namespace
{
  enum fitParameter { 
    // fit parameters related to shifts of primary event vertex
    kPV_shiftX, kPV_shiftY, kPV_shiftZ,
    // fit parameters specific to tau decays
    kTau_visEnFracX, kTau_phi_lab, kTau_decayDistance_lab, kTau_nuInvMass, kTau_pol, 
    kTauVM_theta_rho, kTauVM_mass2_rho, 
    kTauVM_theta_a1, kTauVM_theta_a1r, kTauVM_phi_a1r, kTauVM_mass2_a1, 
    // fit parameters specific to electrons, muons not originating from tau decays
    kLep_shiftEn,
    // fit parameters specific to neutrinos (not originating from tau decays)
    kNu_energy_lab, kNu_phi_lab
  };
}

#endif
