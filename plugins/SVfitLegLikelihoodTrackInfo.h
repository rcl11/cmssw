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
 * $Id: SVfitLegLikelihoodTrackInfo.h,v 1.1 2010/09/21 09:03:00 veelken Exp $
 *
 */

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"

#include "TauAnalysis/CandidateTools/interface/SVfitLegLikelihoodBase.h"
#include "TauAnalysis/CandidateTools/interface/SVfitLegTrackExtractor.h"

#include "AnalysisDataFormats/TauAnalysis/interface/SVfitLegSolution.h"

namespace SVfitLegLikelihoodTrackInfo_namespace 
{
  struct selTrackExtrapolation
  {
    selTrackExtrapolation(const reco::TransientTrack&, const AlgebraicVector3&);

    const AlgebraicVector3& tangent() const { return tangent_; }
    const AlgebraicVector3& dcaPosition() const { return dcaPosition_; }
    const AlgebraicVector3& refPoint() const { return dcaPosition_; }

    double logLikelihood(const AlgebraicVector3&) const;

    AlgebraicVector3 tangent_;
    AlgebraicVector3 dcaPosition_;
    AlgebraicMatrix33 invRotationMatrix_;
    AlgebraicMatrix33 rotCovMatrix_;
    AlgebraicMatrix22 rotCovMatrix2_;

    int errorFlag_;
  };
}

template <typename T>
class SVfitLegLikelihoodTrackInfo : public SVfitLegLikelihoodBase<T>
{
 public:
  SVfitLegLikelihoodTrackInfo(const edm::ParameterSet&);
  ~SVfitLegLikelihoodTrackInfo();

  void beginEvent(const edm::Event&, const edm::EventSetup&);
  void beginCandidate(const T&);

  bool isFittedParameter(int, int) const;

  void setEventVertexPos(const AlgebraicVector3& pvPosition)
  {
    // "original" (unshifted) position of primary event (tau production) vertex
    pvPosition_ = pvPosition;
  }

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

  bool useLinearApprox_;
  AlgebraicVector3 pvPosition_;
  mutable std::vector<SVfitLegLikelihoodTrackInfo_namespace::selTrackExtrapolation> selectedTrackInfo_;
  mutable bool isNewCandidate_;
};

#endif
