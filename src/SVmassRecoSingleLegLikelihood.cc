#include "TauAnalysis/CandidateTools/interface/SVmassRecoSingleLegLikelihood.h"

#include "TauAnalysis/CandidateTools/interface/svMassRecoLikelihoodAuxFunctions.h"

#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"

using namespace svMassReco;

double compLandauValue(double rapidity, double mean, double sigma)
{
  double landau = TMath::Landau(rapidity, mean, sigma, true);
  if ( landau < 1.0e-10 ) landau = 1.0e-10; // sanity check
  return -1.*TMath::Log(landau);
}

double SVmassRecoSingleLegLikelihood::nllVisRapidityGivenMomentumElectronCase(double rapidity, double momentum) const
{
//--- pat::Electron specific customization

  double mean = 2.238*TMath::Power(momentum, 0.2036);
  double sigma = 0.0901 + 0.0006182*momentum;
  return compLandauValue(rapidity, mean, sigma);
}
  
double SVmassRecoSingleLegLikelihood::nllVisRapidityGivenMomentumMuonCase(double rapidity, double momentum) const
{
//--- pat::Muon specific customization  

  double mean = 2.251*TMath::Power(momentum, 0.2013);
  double sigma = 0.09365 + 0.0005191*momentum;
  return compLandauValue(rapidity, mean, sigma);
}
  
double SVmassRecoSingleLegLikelihood::nllVisRapidityGivenMomentumTauJetCase(int legTypeLabel, double rapidity, double momentum) const
{
//--- pat::Tau specific customization
//
//    NOTE: special case as distribution for tau-jet case depends on decay mode  
//

  int decayMode = legTypeLabel;

  double mean = 0;
  double sigma = 0;
  switch ( decayMode ) {
  case 0: // 1 prong 0 pi0
    mean = 1.987*TMath::Power(momentum, 0.2215);
    sigma = 0.2215 + 0.0007142*momentum;
    break;
  case 1: // 1 prong 1 pi0
    mean = 2.03*TMath::Power(momentum, 0.2073);
    sigma = 0.1162 + 0.0001587*momentum;
    break;
  case 2: // 1 prong 2 pi0
    mean = 1.864*TMath::Power(momentum, 0.2149);
    sigma = 0.0962 + 0.0001395*momentum;
    break;
  case 10: // 3 prong 0 pi0
    mean = 1.996*TMath::Power(momentum, 0.1987);
    sigma = 0.1429 + 4.695e-5*momentum;
    break;
  case 11: // 3 prong 1 pi0
    mean = 1.868*TMath::Power(momentum, 0.214);
    sigma = 0.1037 + 4.702e-5*momentum;
    break;
  default: // all others
    mean = 1.996*TMath::Power(momentum, 0.1987);
    sigma = 0.1429 + 4.695e-5*momentum;
    break;
  }
  // Single prong no pi0 case is landau distributed
  if ( decayMode == 0 ) {
    return compLandauValue(rapidity, mean, sigma);
  } else { // otherwise gaussian
    return (square(rapidity-mean)/(2*square(sigma)) + nlGaussianNorm(sigma));
  }
}

std::pair<GlobalPoint, GlobalError> SVmassRecoSingleLegLikelihood::findInitialSecondaryVertex(const GlobalPoint& pv)
{
  // Function that provides an initial, reasonable guess of the location the
  // secondary vertex for a leg in a ditau candidate system, given an initial
  // condition for the primary vertex.  If the leg is a three prong tau, the
  // function will attempt to fit the vertex and return that is the initial SV.
  // Failing that, the method scans along the lead track of the object and finds
  // the point along the lead track that has the lowest negative log likelihood,
  // subject to the constraint that the secondary vertex selection results in a
  // physical solution.  The error returned is the spatial error of the track at
  // that point which may or may not be useful.
  
  edm::LogInfo("findInitialSecondaryVertex") 
    << "Finding a valid SV with PV @ " << "(" << pv.x() << ", " << pv.y() << ", " << pv.z() << ")";
  
  // Get the fourvector of the object
  FourVector visP4 = this->uncorrectedP4();
  edm::LogInfo("findInitialSecondaryVertex") 
    << "Uncorrected vis p4 info: mass = " << visP4.M() << " P = " << visP4.P() 
    << " eta = " << visP4.eta() << " phi = " << visP4.phi();
  
  // Start growing from PV
  GlobalVector basePoint = GlobalVector(pv.x(), pv.y(), pv.z());
  ThreeVector visDir = this->uncorrectedP4().Vect().Unit();
  // Grow the candidate point along visible momentum
  GlobalVector svGrowthDireciton = GlobalVector(visDir.x(), visDir.y(), visDir.z());
  
  // Find lead track
  edm::LogInfo("findInitialSecondaryVertex") << "Finding lead track out of " << tracks_.size() << " tracks";
  size_t highestIndex = 0;
  double highestPt = -1;
  for ( size_t itrk = 0; itrk < tracks_.size(); ++itrk ) {
    TrajectoryStateClosestToPoint tcsp = 
      tracks_[itrk].trajectoryStateClosestToPoint(pv);
    double pt = tcsp.momentum().perp();
    if ( pt > highestPt ) {
      highestIndex = itrk;
      highestPt = pt;
    }
  }
  reco::TransientTrack leadTrack = tracks_[highestIndex];
  edm::LogInfo("SVMethod") << "Found lead track with pt " << highestPt;
  assert(highestPt > 0);
  
  // If we have enough tracks, do a vertex fit to find starting place (or starting R)
  if ( tracks_.size() > 2 ) {
    KalmanVertexFitter fitter(false);
    TransientVertex myVertex = fitter.vertex(tracks_); 
    if ( myVertex.isValid() ) {
      GlobalPoint pos = myVertex.position();
      // Check if this is a valid point for the fitter
      int status = 0;
      this->setPoints(pv, pos.x(), pos.y(), pos.z(), 0.5, status);
      
      edm::LogInfo("findInitialSecondaryVertex") 
	<< "Found a vetex-type starting position at " << "(" << pos.x() << ", " << pos.y() << ", " << pos.z() << ")";
      
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
    this->setPoints(pv, correctedCandSV.x(), correctedCandSV.y(), correctedCandSV.z(), 0.7, error);
    radius = correctedCandSV.perp();
  }
  
  // Now we shoudl be physical, check if we actually found a point or 
  // just reached the highest radius permitted (3cm)
  if ( error == 1 ) 
    edm::LogWarning("SVMethod") <<  "No valid starting point found!" << std::endl;
  
  // Now naively find the minimum NLL along the lead track
  double lowestNLL = this->nllOfLeg();
  GlobalPoint bestCandSV = correctedCandSV;
  
  while ( error == 0 && radius < 3.0 ) {
    // Check if the nll of the last point is better than the current best
    double lastRoundNLL = this->nllOfLeg();
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
    this->setPoints(pv, correctedCandSV.x(), correctedCandSV.y(), correctedCandSV.z(), 0.7, error);
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
