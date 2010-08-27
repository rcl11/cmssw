#ifndef TauAnalysis_CandidateTools_SVfitLegLikelihoodBase_h
#define TauAnalysis_CandidateTools_SVfitLegLikelihoodBase_h

/** \class SVfitLegLikelihoodBase
 *
 * Abstract base-class for plugins computing likelihood for one tau lepton decay "leg";
 * used by SVfit algorithm
 * 
 * \author Evan Friis, Christian Veelken; UC Davis
 *
 * \version $Revision: 1.2 $
 *
 * $Id: SVfitLegLikelihoodBase.h,v 1.2 2010/08/27 06:59:19 veelken Exp $
 *
 */

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "AnalysisDataFormats/TauAnalysis/interface/SVfitLegSolution.h"

template <typename T>
class SVfitLegLikelihoodBase
{
 public:
  SVfitLegLikelihoodBase(const edm::ParameterSet&) {}
  virtual ~SVfitLegLikelihoodBase() {}

  virtual bool isFittedParameter(unsigned) const
  {
    return false;
  }

  virtual bool supportsPolarization() const
  {
    return false;
  }

  virtual double operator()(const T&, const SVfitLegSolution&) const = 0;
};

#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Tau.h"

typedef SVfitLegLikelihoodBase<pat::Electron> SVfitElectronLikelihoodBase;
typedef SVfitLegLikelihoodBase<pat::Muon> SVfitMuonLikelihoodBase;
typedef SVfitLegLikelihoodBase<pat::Tau> SVfitTauLikelihoodBase;

#include "DataFormats/Candidate/interface/Candidate.h"

typedef SVfitLegLikelihoodBase<reco::Candidate> SVfitCandidateLikelihoodBase;

#include "FWCore/PluginManager/interface/PluginFactory.h"

typedef edmplugin::PluginFactory<SVfitElectronLikelihoodBase* (const edm::ParameterSet&)> SVfitElectronLikelihoodBasePluginFactory;
typedef edmplugin::PluginFactory<SVfitMuonLikelihoodBase* (const edm::ParameterSet&)> SVfitMuonLikelihoodBasePluginFactory;
typedef edmplugin::PluginFactory<SVfitTauLikelihoodBase* (const edm::ParameterSet&)> SVfitTauLikelihoodBasePluginFactory;
typedef edmplugin::PluginFactory<SVfitCandidateLikelihoodBase* (const edm::ParameterSet&)> SVfitCandidateLikelihoodBasePluginFactory;

#endif
