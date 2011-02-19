#ifndef TauAnalysis_CandidateTools_SVfitLikelihoodWtauNuMEt_h
#define TauAnalysis_CandidateTools_SVfitLikelihoodWtauNuMEt_h

/** \class SVfitLikelihoodWtauNuMEt
 *
 * Plugin for computing likelihood for neutrinos produced in W boson and tau lepton decays
 * to match missing transverse momentum reconstructed in the event
 *
 * \author Christian Veelken, UC Davis
 *
 * \version $Revision: 1.8 $
 *
 * $Id: SVfitLikelihoodWtauNuMEt.h,v 1.8 2011/01/18 16:47:16 friis Exp $
 *
 */

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/Common/interface/Handle.h"

#include "TauAnalysis/CandidateTools/interface/SVfitWtauNuLikelihoodBase.h"
#include "TauAnalysis/CandidateTools/interface/SVfitLegLikelihoodBase.h"

#include "AnalysisDataFormats/TauAnalysis/interface/CompositePtrCandidateTMEt.h"
#include "AnalysisDataFormats/TauAnalysis/interface/SVfitDiTauSolution.h"

#include <TFormula.h>

template <typename T>
class SVfitLikelihoodWtauNuMEt : public SVfitWtauNuLikelihoodBase<T>
{
 public:
  SVfitLikelihoodWtauNuMEt(const edm::ParameterSet&);
  ~SVfitLikelihoodWtauNuMEt();

  void beginEvent(const edm::Event&, const edm::EventSetup&);
  void beginCandidate(const CompositePtrCandidateTMEt<T>&);

  bool isFittedParameter(int) const;

  double operator()(const CompositePtrCandidateTMEt<T>&, const SVfitWtauNuSolution&) const;
 private:
  TFormula* parSigma_;
  TFormula* parBias_;

  TFormula* perpSigma_;
  TFormula* perpBias_;

  double qX_;
  double qY_;
  double qT_;

  // Whether or not to allow the phi variable to vary
  bool varyPhi_;

  edm::InputTag srcPFCandidates_;
  edm::Handle<reco::PFCandidateCollection> pfCandidates_;
};

#endif
