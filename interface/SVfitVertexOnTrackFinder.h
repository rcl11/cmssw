#ifndef TauAanlysis_CandidateTools_SVfitVertexOnTrackFinder_h
#define TauAanlysis_CandidateTools_SVfitVertexOnTrackFinder_h

/*
 * Given a set of tracks, a visible tau direction, the labframe opening angle of
 * the decay, and a primary vertex * compute the approximate decay vertex of a
 * tau lepton.
 *
 * Author: Evan K. Friis, Christian Veelken (UC Davis)
 *
 */

#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "AnalysisDataFormats/TauAnalysis/interface/SVfitLegSolution.h"

namespace SVfit { namespace track {

class VertexOnTrackFinder {
  public:
    VertexOnTrackFinder(const SVfitLegSolution& solution);

    GlobalPoint decayVertex(const GlobalPoint& pv,
                            double angleVisLabFrame);

    // Version with corrections applied.
    // The direction object is updated to point a long the reconstructed tau
    // direction. Note that this direction will still point in the intended
    // direction even if the flight correction is 'backwards' and larger than
    // the total flight distance of the tau.
    GlobalPoint decayVertex(const GlobalPoint& pv,
                            double angleVisLabFrame,
                            double phiCorrection,
                            double flightCorrection,
                            reco::Candidate::Vector *direction);

    // hate all these different types!
    GlobalPoint decayVertex(const AlgebraicVector3& pv,
                            double angleVisLabFrame,
                            double phiCorrection,
                            double flightCorrection,
                            reco::Candidate::Vector *direction);

    static void setEventSetup(const edm::EventSetup& es) {
      edm::ESHandle<TransientTrackBuilder> trackBuilderHandle;
      es.get<TransientTrackRecord>().get("TransientTrackBuilder",
                                         trackBuilderHandle);
      builder_ = trackBuilderHandle.product();
    }

  private:
    // We only need the track for one prongs.
    const SVfitLegSolution& soln_;
    reco::TransientTrack track_;
    static const TransientTrackBuilder* builder_;
};

}}  // end namespace SVfit::track
#endif
