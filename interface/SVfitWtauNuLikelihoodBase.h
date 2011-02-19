#ifndef TauAnalysis_CandidateTools_SVfitWtauNuLikelihoodBase_h
#define TauAnalysis_CandidateTools_SVfitWtauNuLikelihoodBase_h

/** \class SVfitWtauNuLikelihoodBase
 *
 * Abstract base-class for plugins computing likelihood for tau lepton plus neutrino pair produced in W boson decay;
 * used by SVfit algorithm
 *
 * \author Christian Veelken, UC Davis
 *
 * \version $Revision: 1.12 $
 *
 * $Id: SVfitWtauNuLikelihoodBase.h,v 1.12 2011/01/18 16:42:29 friis Exp $
 *
 */

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "AnalysisDataFormats/TauAnalysis/interface/CompositePtrCandidateTMEt.h"
#include "AnalysisDataFormats/TauAnalysis/interface/SVfitWtauNuSolution.h"

#include <string>
#include <iostream>

template <typename T>
class SVfitWtauNuLikelihoodBase
{
 public:
  SVfitWtauNuLikelihoodBase(const edm::ParameterSet& cfg)
  {
    verbosity_ = cfg.exists("verbosity") ?
      cfg.getParameter<int>("verbosity") : 0;
    pluginType_ = cfg.getParameter<std::string>("pluginType");
    pluginName_ = cfg.getParameter<std::string>("pluginName");
    firstFit_ = cfg.getParameter<unsigned int>("firstFitIteration");
  }
  virtual ~SVfitWtauNuLikelihoodBase() {}

  const std::string& name() const { return pluginName_; }

  virtual void beginJob() {}
  virtual void beginEvent(const edm::Event&, const edm::EventSetup&) {}
  virtual void beginCandidate(const CompositePtrCandidateTMEt<T>&) {}

  virtual void print(std::ostream& stream) const
  {
    stream << "<SVfitWtauNuLikelihoodBase::print>:" << std::endl;
    stream << " pluginType = " << pluginType_ << std::endl;
  }

  virtual bool isFittedParameter(int parNo) const {
    return false;
  }

  virtual bool supportsPolarization() const {
    return false;
  }

  // Check if this likelihood should be used in this fit iteration.
  bool isFitted(unsigned int fitIter) const {
    return fitIter >= firstFit_;
  }

  // Get the first fit iteration this is valid
  unsigned int firstFit() const {
    return firstFit_;
  }

  virtual double operator()(const CompositePtrCandidateTMEt<T>&, const SVfitWtauNuSolution&) const = 0;
 protected:
  std::string pluginType_;
  std::string pluginName_;
  unsigned int firstFit_;
  int verbosity_;
};

#include "DataFormats/PatCandidates/interface/Tau.h"

typedef SVfitWtauNuLikelihoodBase<pat::Tau> SVfitTauNuPairLikelihoodBase;

#include "FWCore/PluginManager/interface/PluginFactory.h"

typedef edmplugin::PluginFactory<SVfitTauNuPairLikelihoodBase* (const edm::ParameterSet&)> SVfitTauNuPairLikelihoodBasePluginFactory;

#endif
