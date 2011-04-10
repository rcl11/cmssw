#ifndef TauAnalysis_CandidateTools_NSVfitTauToLepBuilder_h
#define TauAnalysis_CandidateTools_NSVfitTauToLepBuilder_h

/** \class NSVfitSingleParticleBuilderBase
 *
 * Auxiliary class reconstructing tau --> had decays and
 * building NSVfitTauToLepHypothesis objects;
 * used by NSVfit algorithm
 *
 * \author Christian Veelken, UC Davis
 *
 * \version $Revision: 1.4 $
 *
 * $Id: NSVfitTauToLepBuilder.h,v 1.4 2011/03/27 14:22:35 friis Exp $
 *
 */

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "TauAnalysis/CandidateTools/interface/NSVfitSingleParticleBuilderBase.h"
#include "TauAnalysis/CandidateTools/interface/SVfitLegTrackExtractor.h"

#include <string>
#include <iostream>

template <typename T>
class NSVfitTauToLepBuilder : public NSVfitSingleParticleBuilderBase
{
 public:
  NSVfitTauToLepBuilder(const edm::ParameterSet&);
  ~NSVfitTauToLepBuilder();

  void beginJob(NSVfitAlgorithmBase*);

  virtual NSVfitSingleParticleHypothesisBase* build(const inputParticleMap&) const;

  void applyFitParameter(NSVfitSingleParticleHypothesisBase*, const double*) const;

  void print(std::ostream& stream) const;

 protected:
  int idxFitParameter_visEnFracX_;
  int idxFitParameter_phi_lab_;
  int idxFitParameter_nuInvMass_;

  SVfitLegTrackExtractor<T> trackExtractor_;

  NSVfitAlgorithmBase* algorithm_;
};

#include "FWCore/PluginManager/interface/PluginFactory.h"

typedef edmplugin::PluginFactory<NSVfitSingleParticleBuilderBase* (const edm::ParameterSet&)> NSVfitSingleParticleBuilderPluginFactory;

#endif



