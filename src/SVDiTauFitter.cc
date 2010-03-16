#include "TauAnalysis/CandidateTools/interface/SVDiTauFitter.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"

namespace TauVertex {

/// Find a reasonable valid intial guess for the SV on a leg
std::pair<GlobalPoint, GlobalError> findValidSV(const GlobalPoint& pv, LegFitter* leg)
{
   edm::LogInfo("SVMethod") << "Finding a valid SV with PV @ "
      << "(" << pv.x() << ", " << pv.y() << ", " << pv.z() << ")";

   FourVector visP4 = leg->uncorrectedP4();
   edm::LogInfo("SVMethod") << "Uncorrected vis p4 info: mass: " << visP4.M() << " P: " << visP4.P() 
      << " eta: " << visP4.eta() << " phi: " << visP4.phi();

   // Start growing from PV
   GlobalVector basePoint = GlobalVector(pv.x(), pv.y(), pv.z());
   ThreeVector visDir = leg->uncorrectedP4().Vect().Unit();
   // Grow the candidate point along visible momentum
   GlobalVector svGrowthDireciton = GlobalVector(visDir.x(), visDir.y(), visDir.z());

   // Find lead track
   size_t highestIndex = 0;
   double highestPt = -1;
   size_t nTracks = leg->tracks().size();
   edm::LogInfo("SVMethod") << "Finding lead track out of " << nTracks << " tracks";
   for(size_t itrk = 0; itrk < nTracks; ++itrk)
   {
      TrajectoryStateClosestToPoint tcsp = 
         leg->tracks()[itrk].trajectoryStateClosestToPoint(pv);
      double pt = tcsp.momentum().perp();
      if(pt > highestPt)
      {
         highestIndex = itrk;
         highestPt = pt;
      }
   }
   reco::TransientTrack leadTrack = leg->tracks()[highestIndex];
   edm::LogInfo("SVMethod") << "Found lead track with pt " << highestPt;
   assert(highestPt > 0);

   /// If we have enough tracks, do a vertex fit to find starting place (or starting R)
   if(leg->tracks().size() > 2)
   {
      KalmanVertexFitter fitter(false);
      TransientVertex myVertex = fitter.vertex(leg->tracks()); 
      if(myVertex.isValid())
      {
         GlobalPoint pos = myVertex.position();
         // Check if this is a valid point for the fitter
         int status = 0;
         leg->setPoints(pv, pos.x(), pos.y(), pos.z(), 0.5, status);

         edm::LogInfo("SVMethod") << "Found a vetex-type starting position at " 
            << "(" << pos.x() << ", " << pos.y() << ", " << pos.z() << ")";

         // If vertex is physical return that
         if( status == 0 )
         {
            edm::LogInfo("SVMethod") << " the vertex is valid, using it!";
            return std::make_pair(pos, myVertex.positionError());
         } 

         // Otherwise set start point
         basePoint = GlobalVector(pos.x(), pos.y(), pos.z());
      }
   } 
   
   // Slowly expand the cylindar, find the lead track intersection with the cylinder,
   // and use that as the candidate SV point
   double stepsize = 5e-4; // 5 micron
   double currentScale = 0;
   double radius = 0;
   GlobalPoint correctedCandSV;

   int error = 1; // will stay 1 while in unphysical region
   while(error == 1 && radius < 3.0)
   {
      currentScale += stepsize;
      GlobalVector candSV = basePoint + svGrowthDireciton*currentScale;
      GlobalPoint candSVPoint(candSV.x(), candSV.y(), candSV.z()); // BS
      // Find the closest point on the lead track
      TrajectoryStateClosestToPoint tcsp = leadTrack.trajectoryStateClosestToPoint(candSVPoint);
      correctedCandSV = tcsp.position();
      // See if this point works (error will go to zero
      leg->setPoints(pv, correctedCandSV.x(), correctedCandSV.y(), correctedCandSV.z(), 0.7, error);
      radius = correctedCandSV.perp();
   }

   // Check if we actually found a point
   if(error == 1)
      edm::LogWarning("SVMethod") <<  "No valid starting point found!" << std::endl;

   // Now naively find the minimum NLL along the lead track
   double lowestNLL = leg->nllOfLeg();
   GlobalPoint bestCandSV = correctedCandSV;

   while(error == 0 && radius < 3.0)
   {
      // Check if the nll of the last point is better than the current best
      double lastRoundNLL = leg->nllOfLeg();
      if(lastRoundNLL < lowestNLL)
      {
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
      //correctedCandSV = candSVPoint;
      // See if this point works (error will go to zero
      leg->setPoints(pv, correctedCandSV.x(), correctedCandSV.y(), correctedCandSV.z(), 0.7, error);
      radius = correctedCandSV.perp();
      //edm::LogInfo("SVMethodScan") << "Radius " << radius << " NLL: " << lastRoundNLL << " LowestNLL: " << lowestNLL;
   }

   // Get the track error at our best point
   GlobalError errorAtBestPoint = 
      leadTrack.trajectoryStateClosestToPoint(bestCandSV).theState().cartesianError().position();

   edm::LogInfo("SVMethod") << "Returning track-type starting position at " 
      << "(" << bestCandSV.x() << ", " << bestCandSV.y() << ", " << bestCandSV.z() << ") NLL = " << lowestNLL;
   return std::make_pair(bestCandSV, errorAtBestPoint);
}
}
