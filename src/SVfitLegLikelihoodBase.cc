#include "TauAnalysis/CandidateTools/interface/SVfitLegLikelihoodBase.h"

template <typename T>
SVfitLegLikelihoodBase<T>::SVfitLegLikelihoodBase(const edm::ParameterSet& cfg)
{
// nothing to be done yet...
}

template <typename T>
SVfitLegLikelihoodBase<T>::~SVfitLegLikelihoodBase()
{
// nothing to be done yet...
}

template <typename T>
bool SVfitLegLikelihoodBase<T>::isFittedParameter(unsigned index)
{
  return false;
}

#include "FWCore/Framework/interface/MakerMacros.h"

EDM_REGISTER_PLUGINFACTORY(SVfitElectronLikelihoodBasePluginFactory, "SVfitElectronLikelihoodBasePluginFactory");
EDM_REGISTER_PLUGINFACTORY(SVfitMuonLikelihoodBasePluginFactory, "SVfitMuonLikelihoodBasePluginFactory");
EDM_REGISTER_PLUGINFACTORY(SVfitTauLikelihoodBasePluginFactory, "SVfitTauLikelihoodBasePluginFactory");
EDM_REGISTER_PLUGINFACTORY(SVfitCandidateLikelihoodBasePluginFactory, "SVfitCandidateLikelihoodBasePluginFactory");
