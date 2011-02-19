#include "TauAnalysis/CandidateTools/plugins/SVfitLikelihoodWtauNuKinematics.h"

#include "TauAnalysis/CandidateTools/interface/SVfitAlgorithmWtauNu.h"
#include "TauAnalysis/CandidateTools/interface/svFitAuxFunctions.h"
#include "AnalysisDataFormats/TauAnalysis/interface/tauAnalysisAuxFunctions.h"

using namespace SVfit_namespace;

template <typename T>
SVfitLegLikelihoodBase<T>* createLikelihoodPlugin(const edm::ParameterSet& cfg)
{
  std::string pluginType = cfg.getParameter<std::string>("pluginType");

  typedef edmplugin::PluginFactory<SVfitLegLikelihoodBase<T>* (const edm::ParameterSet&)> SVfitLegLikelihoodPluginFactory;
  SVfitLegLikelihoodPluginFactory* pluginFactory = SVfitLegLikelihoodPluginFactory::get();
  
//--- print error message in case plugin of specified type cannot be created
  if ( !pluginFactory->tryToCreate(pluginType, cfg) ) {
    edm::LogError ("createLikelihoodPlugin") 
      << "Failed to create plugin of type = " << pluginType << " !!";
    std::cout << " category = " << pluginFactory->category() << std::endl;
    std::cout << " available plugins = { ";
    std::vector<edmplugin::PluginInfo> plugins = pluginFactory->available();
    unsigned numPlugins = plugins.size();
    for ( unsigned iPlugin = 0; iPlugin < numPlugins; ++iPlugin ) {
      std::cout << plugins[iPlugin].name_;
      if ( iPlugin < (numPlugins - 1) ) std::cout << ", ";
    }
    std::cout << " }" << std::endl;
  }

  return pluginFactory->create(pluginType, cfg);
}

template <typename T>
SVfitLikelihoodWtauNuKinematics<T>::SVfitLikelihoodWtauNuKinematics(const edm::ParameterSet& cfg)
  : SVfitWtauNuLikelihoodBase<T>(cfg)
{
//--- initialize plugins computing likelihoods for the tau lepton decay "leg"
  edm::ParameterSet cfgLeg1Likelihood;
  if ( cfg.exists("leg1") ) {
    cfgLeg1Likelihood = cfg.getParameter<edm::ParameterSet>("leg1");
  } else {
    cfgLeg1Likelihood = cfg.getParameter<edm::ParameterSet>("leg");
  }

  leg1Likelihood_ = createLikelihoodPlugin<T>(cfgLeg1Likelihood);
}

template <typename T>
SVfitLikelihoodWtauNuKinematics<T>::~SVfitLikelihoodWtauNuKinematics()
{
  delete leg1Likelihood_;
}

template <typename T>
void SVfitLikelihoodWtauNuKinematics<T>::beginJob()
{
  leg1Likelihood_->beginJob();
}

template <typename T>
void SVfitLikelihoodWtauNuKinematics<T>::beginEvent(edm::Event& evt, const edm::EventSetup& es)
{
  leg1Likelihood_->beginEvent(evt, es);
}

template <typename T>
void SVfitLikelihoodWtauNuKinematics<T>::beginCandidate(const CompositePtrCandidateTMEt<T>& candidate)
{
  leg1Likelihood_->beginCandidate(*candidate.visDecayProducts());
}

template <typename T>
void SVfitLikelihoodWtauNuKinematics<T>::print(std::ostream& stream) const
{
  SVfitWtauNuLikelihoodBase<T>::print(stream);
  leg1Likelihood_->print(stream);
}

template <typename T>
bool SVfitLikelihoodWtauNuKinematics<T>::isFittedParameter(int index) const
{
  if ( index == SVfit_namespace::kLeg1thetaRest        || 
       index == SVfit_namespace::kLeg1phiLab           ||
       index == SVfit_namespace::kLeg1decayDistanceLab ||
       index == SVfit_namespace::kLeg1nuInvMass        ||
       index == SVfit_namespace::kLeg1thetaVMrho       ||
       index == SVfit_namespace::kLeg1thetaVMa1        ||
       index == SVfit_namespace::kLeg1thetaVMa1r       ||
       index == SVfit_namespace::kLeg1phiVMa1r         ) 
    return leg1Likelihood_->isFittedParameter(SVfit_namespace::kLeg1, index);
  else return false;
}

template <typename T>
bool SVfitLikelihoodWtauNuKinematics<T>::supportsPolarization() const
{
  return ( leg1Likelihood_->supportsPolarization() );
}

template <typename T>
double SVfitLikelihoodWtauNuKinematics<T>::operator()(const CompositePtrCandidateTMEt<T>& candidate, 
						      const SVfitWtauNuSolution& solution) const
{
//--- compute negative log-likelihood for tau lepton hypothesis given as function argument
  double negativeLogLikelihood = (*leg1Likelihood_)(*candidate.visDecayProducts(), solution.leg1());
  
  return negativeLogLikelihood;
}

#include "DataFormats/PatCandidates/interface/Tau.h"

typedef SVfitLikelihoodWtauNuKinematics<pat::Tau> SVfitLikelihoodTauNuPairKinematics;

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_EDM_PLUGIN(SVfitTauNuPairLikelihoodBasePluginFactory, SVfitLikelihoodTauNuPairKinematics, "SVfitLikelihoodTauNuPairKinematics");
