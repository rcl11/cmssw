#include "TauAnalysis/CandidateTools/plugins/NSVfitResonanceLikelihoodMassPenalty.h"

#include <TMath.h>
#include <TString.h>

#include <string>

NSVfitResonanceLikelihoodMassPenalty::NSVfitResonanceLikelihoodMassPenalty(const edm::ParameterSet& cfg)
  : NSVfitResonanceLikelihood(cfg)
{
  TString nll_formula_string = cfg.getParameter<std::string>("nll").data();
  nll_formula_string.ReplaceAll("mass", "x");
  
  nll_formula_ = new TFormula(TString(pluginName_.data()).Append("_nll").Data(), nll_formula_string.Data());

  power_ = cfg.getParameter<double>("power");
}

NSVfitResonanceLikelihoodMassPenalty::~NSVfitResonanceLikelihoodMassPenalty()
{
  delete nll_formula_;
}

double NSVfitResonanceLikelihoodMassPenalty::operator()(const NSVfitResonanceHypothesis* resonance) const 
{
  assert(resonance);

  double mass = resonance->p4_fitted().mass();

  double nll = power_*nll_formula_->Eval(mass);
  //std::cout << "<NSVfitResonanceLikelihoodMassPenalty::operator()>:" << std::endl;
  //std::cout << " formula = " << nll_formula_->GetTitle() << ": mass = " << mass << " --> nll = " << nll << std::endl;
  
  return nll;
}

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_EDM_PLUGIN(NSVfitResonanceLikelihoodPluginFactory, NSVfitResonanceLikelihoodMassPenalty, "NSVfitResonanceLikelihoodMassPenalty");
