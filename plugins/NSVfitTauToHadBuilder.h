#ifndef TauAnalysis_CandidateTools_NSVfitTauToHadBuilder_h
#define TauAnalysis_CandidateTools_NSVfitTauToHadBuilder_h

/** \class NSVfitSingleParticleBuilderBase
 *
 * Auxiliary class reconstructing tau --> had decays and
 * building NSVfitTauToHadHypothesis objects;
 * used by NSVfit algorithm
 *
 * \author Christian Veelken, UC Davis
 *
 * \version $Revision: 1.1 $
 *
 * $Id: NSVfitTauToHadBuilder.h,v 1.1 2011/02/27 16:45:16 veelken Exp $
 *
 */

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "TauAnalysis/CandidateTools/interface/NSVfitSingleParticleBuilderBase.h"
#include "TauAnalysis/CandidateTools/interface/SVfitLegTrackExtractor.h"

#include <string>
#include <iostream>

class NSVfitTauToHadBuilder : public NSVfitSingleParticleBuilderBase
{
 public:
  NSVfitTauToHadBuilder(const edm::ParameterSet&);
  ~NSVfitTauToHadBuilder();

  void initialize(NSVfitAlgorithmBase*);

  virtual NSVfitSingleParticleHypothesisBase* build(const inputParticleMap&) const;

  void applyFitParameter(NSVfitSingleParticleHypothesisBase*, double*) const;

  void print(std::ostream& stream) const;

protected:
  int idxFitParameter_theta_rf_;
  int idxFitParameter_phi_lab_;

  SVfitLegTrackExtractor<pat::Tau> trackExtractor_;
};

#include "FWCore/PluginManager/interface/PluginFactory.h"

typedef edmplugin::PluginFactory<NSVfitSingleParticleBuilderBase* (const edm::ParameterSet&)> NSVfitSingleParticleBuilderPluginFactory;

#endif



