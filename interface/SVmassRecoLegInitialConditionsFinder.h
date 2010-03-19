#ifndef TauAnalysisTools_CandidateTools_SVmassRecoLegInitialConditionsFinder_h
#define TauAnalysisTools_CandidateTools_SVmassRecoLegInitialConditionsFinder_h

/*
 * SVmassRecoLegInitialConditionsFinder
 *
 * Authors: Evan Friis, Christian Veelken (UC Davis)
 *
 * Function that provides an initial, reasonable guess of the location the
 * secondary vertex for a leg in a ditau candidate system, given an initial
 * condition for the primary vertex.  If the leg is a three prong tau, the
 * function will attempt to fit the vertex and return that is the initial SV.
 * Failing that, the method scans along the lead track of the object and finds
 * the point along the lead track that has the lowest negative log likelihood,
 * subject to the constraint that the secondary vertex selection results in a
 * physical solution.  The error returned is the spatial error of the track at
 * that point which may or may not be useful.
 *
 */

#include "TauAnalysis/CandidateTools/interface/SVmassRecoSingleLegLikelihood.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"

namespace svMassReco 
{
  // forward declaration
  template<typename T> class SVmassRecoSingleLegLikelihood;

  template<typename T> std::pair<GlobalPoint, GlobalError> 
  findInitialSecondaryVertex(const GlobalPoint& pv, SVmassRecoSingleLegLikelihood<T> *leg)
  {
    edm::LogInfo("SVMethod") << "Finding a valid SV with PV @ "
			     << "(" << pv.x() << ", " << pv.y() << ", " << pv.z() << ")";
    
    // Get the fourvector of the object
    FourVector visP4 = leg->uncorrectedP4();
    edm::LogInfo("SVMethod") << "Uncorrected vis p4 info: mass: " << visP4.M() << " P: " << visP4.P() 
			     << " eta: " << visP4.eta() << " phi: " << visP4.phi();
    
    // Start growing from PV
    GlobalVector basePoint = GlobalVector(pv.x(), pv.y(), pv.z());
    ThreeVector visDir = leg->uncorrectedP4().Vect().Unit();
    // Grow the candidate point along visible momentum
    GlobalVector svGrowthDireciton = GlobalVector(visDir.x(), visDir.y(), visDir.z());
    
    // Find lead track
    edm::LogInfo("SVMethod") << "Finding lead track out of " << leg->tracks().size() << " tracks";
    size_t highestIndex = 0;
    double highestPt = -1;
    for ( size_t itrk = 0; itrk < leg->tracks().size(); ++itrk ) {
      TrajectoryStateClosestToPoint tcsp = 
	leg->tracks()[itrk].trajectoryStateClosestToPoint(pv);
      double pt = tcsp.momentum().perp();
      if ( pt > highestPt ) {
	highestIndex = itrk;
	highestPt = pt;
      }
    }
    reco::TransientTrack leadTrack = leg->tracks()[highestIndex];
    edm::LogInfo("SVMethod") << "Found lead track with pt " << highestPt;
    assert(highestPt > 0);
    
    /// If we have enough tracks, do a vertex fit to find starting place (or starting R)
    if ( leg->tracks().size() > 2 ) {
      KalmanVertexFitter fitter(false);
      TransientVertex myVertex = fitter.vertex(leg->tracks()); 
      if ( myVertex.isValid() ) {
	GlobalPoint pos = myVertex.position();
	// Check if this is a valid point for the fitter
	int status = 0;
	leg->setPoints(pv, pos.x(), pos.y(), pos.z(), 0.5, status);
	
	edm::LogInfo("SVMethod") << "Found a vetex-type starting position at " 
				 << "(" << pos.x() << ", " << pos.y() << ", " << pos.z() << ")";
	
	// If vertex is physical return that
	if( status == 0 ) {
	  edm::LogInfo("SVMethod") << " the vertex is valid, using it!";
	  return std::make_pair(pos, myVertex.positionError());
	} 

	// Otherwise set start point as the fitted SV
	basePoint = GlobalVector(pos.x(), pos.y(), pos.z());
      }
    } 
   
    double stepsize = 5e-4; // 5 micron
    double currentScale = 0;
    double radius = 0;
    GlobalPoint correctedCandSV;
    
    // Walk out along the lead track until we are in a physical region
    int error = 1; 
    while ( error == 1 && radius < 3.0 ) {
      currentScale += stepsize;
      GlobalVector candSV = basePoint + svGrowthDireciton*currentScale;
      GlobalPoint candSVPoint(candSV.x(), candSV.y(), candSV.z()); // BS
      // Find the closest point on the lead track
      TrajectoryStateClosestToPoint tcsp = leadTrack.trajectoryStateClosestToPoint(candSVPoint);
      correctedCandSV = tcsp.position();
      // See if this point works (error will go to zero)
      leg->setPoints(pv, correctedCandSV.x(), correctedCandSV.y(), correctedCandSV.z(), 0.7, error);
      radius = correctedCandSV.perp();
    }
    
    // Now we shoudl be physical, check if we actually found a point or 
    // just reached the highest radius permitted (3cm)
    if ( error == 1 ) 
       edm::LogWarning("SVMethod") <<  "No valid starting point found!" << std::endl;
    
    // Now naively find the minimum NLL along the lead track
    double lowestNLL = leg->nllOfLeg();
    GlobalPoint bestCandSV = correctedCandSV;
    
    while ( error == 0 && radius < 3.0 ) {
      // Check if the nll of the last point is better than the current best
      double lastRoundNLL = leg->nllOfLeg();
      if ( lastRoundNLL < lowestNLL ) {
	bestCandSV = correctedCandSV;
	lowestNLL = lastRoundNLL;
      }
      
      currentScale += stepsize;
      GlobalVector candSV = basePoint + svGrowthDireciton*currentScale;
      GlobalPoint candSVPoint(candSV.x(), candSV.y(), candSV.z()); // BS
      // Find the closest point on the lead track
      TrajectoryStateClosestToPoint tcsp = leadTrack.trajectoryStateClosestToPoint(candSVPoint);
      // FIXME
      correctedCandSV = tcsp.position();
      // See if this point works (error will go to zero
      leg->setPoints(pv, correctedCandSV.x(), correctedCandSV.y(), correctedCandSV.z(), 0.7, error);
      radius = correctedCandSV.perp();
    }
    
    // Get the track error at our best point
    GlobalError errorAtBestPoint = 
      leadTrack.trajectoryStateClosestToPoint(bestCandSV).theState().cartesianError().position();
    
    edm::LogInfo("SVMethod") 
      << "Returning track-type starting position at " 
      << "(" << bestCandSV.x() << ", " << bestCandSV.y() << ", " 
      << bestCandSV.z() << ") NLL = " << lowestNLL;

    return std::make_pair(bestCandSV, errorAtBestPoint);
  }
}

#endif
