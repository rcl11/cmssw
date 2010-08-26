#include "TauAnalysis/CandidateTools/interface/SVfitDiTauLikelihoodBase.h"

template <typename T1, typename T2>
SVfitDiTauLikelihoodBase<T1,T2>::SVfitDiTauLikelihoodBase(const edm::ParameterSet& cfg)
{
// nothing to be done yet...
}

template <typename T1, typename T2>
SVfitDiTauLikelihoodBase<T1,T2>::~SVfitDiTauLikelihoodBase()
{
// nothing to be done yet...
}

template <typename T1, typename T2>
bool SVfitDiTauLikelihoodBase<T1,T2>::isFittedParameter(unsigned index)
{
  return false;
}

#include "FWCore/Framework/interface/MakerMacros.h"

EDM_REGISTER_PLUGINFACTORY(SVfitElecTauPairLikelihoodBasePluginFactory, "SVfitElecTauPairLikelihoodBasePluginFactory");
EDM_REGISTER_PLUGINFACTORY(SVfitMuTauPairLikelihoodBasePluginFactory, "SVfitMuTauPairLikelihoodBasePluginFactory");
EDM_REGISTER_PLUGINFACTORY(SVfitDiTauPairLikelihoodBasePluginFactory, "SVfitDiTauPairLikelihoodBasePluginFactory");
EDM_REGISTER_PLUGINFACTORY(SVfitElecMuPairLikelihoodBasePluginFactory, "SVfitElecMuPairLikelihoodBasePluginFactory");
EDM_REGISTER_PLUGINFACTORY(SVfitDiCandidatePairLikelihoodBasePluginFactory, "SVfitDiCandidatePairLikelihoodBasePluginFactory");
