#include "TauAnalysis/CandidateTools/plugins/NSVfitResonanceLikelihoodMassPenalty.h"

#include <TMath.h>
#include <TString.h>

#include <string>

NSVfitResonanceLikelihoodMassPenalty::NSVfitResonanceLikelihoodMassPenalty(const edm::ParameterSet& cfg)
  : NSVfitResonanceLikelihood(cfg)
{
  TString formula_string = cfg.getParameter<std::string>("formula").data();
  formula_string.ReplaceAll("mass", "x");
  
  formula_ = new TFormula(TString(pluginName_.data()).Append("_formula").Data(), formula_string.Data());

  power_ = cfg.getParameter<double>("power");
}

NSVfitResonanceLikelihoodMassPenalty::~NSVfitResonanceLikelihoodMassPenalty()
{
  delete formula_;
}

double NSVfitResonanceLikelihoodMassPenalty::operator()(const NSVfitResonanceHypothesis* resonance) const 
{
  assert(resonance);

  double mass = resonance->p4_fitted().mass();

  double retVal = power_*formula_->Eval(mass);
  std::cout << "<NSVfitResonanceLikelihoodMassPenalty::operator()>:" << std::endl;
  std::cout << " formula = " << formula_->GetTitle() << ": mass = " << mass << " --> retVal = " << retVal << std::endl;
  
  return retVal;
}

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_EDM_PLUGIN(NSVfitResonanceLikelihoodPluginFactory, NSVfitResonanceLikelihoodMassPenalty, "NSVfitResonanceLikelihoodMassPenalty");
