#include "TauAnalysis/CandidateTools/plugins/NSVfitResonanceLikelihoodSculpting.h"

#include "TauAnalysis/CandidateTools/interface/NSVfitAlgorithmBase.h"
#include "TauAnalysis/CandidateTools/interface/svFitAuxFunctions.h"

#include <string>

NSVfitResonanceLikelihoodSculpting::NSVfitResonanceLikelihoodSculpting(
    const edm::ParameterSet& cfg)
  : NSVfitResonanceLikelihood(cfg),
    meanFunc_("meanFunc", cfg.getParameter<std::string>("meanFunction").data()),
    rmsFunc_("rmsFunc", cfg.getParameter<std::string>("rmsFunction").data())
{}

NSVfitResonanceLikelihoodSculpting::~NSVfitResonanceLikelihoodSculpting() { }

void NSVfitResonanceLikelihoodSculpting::beginJob(NSVfitAlgorithmBase* algorithm) const
{
  algorithm->requestFitParameter("allTauDecays", nSVfit_namespace::kTau_visEnFracX, pluginName_);
  algorithm->requestFitParameter("allTauDecays", nSVfit_namespace::kTau_phi_lab,    pluginName_);
  algorithm->requestFitParameter("allLeptons",   nSVfit_namespace::kLep_shiftEn,    pluginName_);
  algorithm->requestFitParameter("allNeutrinos", nSVfit_namespace::kNu_energy_lab,  pluginName_);
  algorithm->requestFitParameter("allNeutrinos", nSVfit_namespace::kNu_phi_lab,     pluginName_);
}

double
NSVfitResonanceLikelihoodSculpting::operator()(const NSVfitResonanceHypothesis* hypothesis) const
{
//--- compute negative log-likelihood for two tau leptons
//    to have transverse momenta leg1Pt, leg2Pt

  if ( this->verbosity_ )
    std::cout << "<NSVfitLikelihoodDiTauSculpting::operator()>:" << std::endl;

  double diTauMass = hypothesis->p4_fitted().mass();
  double visMass = hypothesis->p4().mass();
  double scaledVisMass = visMass/diTauMass;
  double residual = scaledVisMass - meanFunc_.Eval(diTauMass);
  double sigma = rmsFunc_.Eval(diTauMass);
  double nll = -SVfit_namespace::logGaussian(residual, sigma);

  if ( this->verbosity_ ) {
    std::cout << " diTauMass = " << diTauMass
      << " visMass = " << visMass
      << " scaledVisMass = " << scaledVisMass
      << " residual = " << residual
      << " sigma = " << sigma << std::endl;
  }
  return nll;
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_EDM_PLUGIN(NSVfitResonanceLikelihoodPluginFactory,
    NSVfitResonanceLikelihoodSculpting, "NSVfitResonanceLikelihoodSculpting");
