#ifndef TauAnalysis_CandidateTools_SVfitLegLikelihoodTrackInfo_h
#define TauAnalysis_CandidateTools_SVfitLegLikelihoodTrackInfo_h

/** \class SVfitLegLikelihoodTrackInfo
 *
 * Plugin for computing likelihood for tracks of tau lepton decay "leg"
 * to be compatible with originating from hypothetic secondary (tau lepton decay) vertex
 * 
 * \author Evan Friis, Christian Veelken; UC Davis
 *
 * \version $Revision: 1.1 $
 *
 * $Id: SVfitLegLikelihoodTrackInfo.h,v 1.1 2010/08/28 14:12:06 veelken Exp $
 *
 */

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"

#include "TauAnalysis/CandidateTools/interface/SVfitLegLikelihoodBase.h"
#include "TauAnalysis/CandidateTools/interface/SVfitLegTrackExtractor.h"

#include "AnalysisDataFormats/TauAnalysis/interface/SVfitLegSolution.h"

template <typename T>
class SVfitLegLikelihoodTrackInfo : public SVfitLegLikelihoodBase<T>
{
 public:
  SVfitLegLikelihoodTrackInfo(const edm::ParameterSet&);
  ~SVfitLegLikelihoodTrackInfo();

  void beginEvent(const edm::Event&, const edm::EventSetup&);
  void beginCandidate(const T&);

  bool isFittedParameter(int, int) const;

  double operator()(const T&, const SVfitLegSolution&) const;
 private:
  SVfitLegTrackExtractor<T> trackExtractor_;

  const TransientTrackBuilder* trackBuilder_;

  mutable std::vector<reco::TransientTrack> selectedTracks_;

  unsigned minNumHits_;
  unsigned minNumPixelHits_;
  double maxChi2DoF_;
  double maxDeltaPoverP_;
  double minPt_;
};

#endif
