#include "TauAnalysis/CandidateTools/plugins/SVfitLegLikelihoodTrackInfo.h"

#include "TauAnalysis/CandidateTools/interface/svFitAuxFunctions.h"

#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateClosestToPoint.h"

#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/CLHEP/interface/AlgebraicObjects.h"

#include <TRotation.h>

#include <limits>

using namespace SVfit_namespace;
using namespace SVfitLegLikelihoodTrackInfo_namespace;

unsigned defaultMinNumHits = 5;
unsigned defaultMinNumPixelHits = 1;
double defaultMaxChi2DoF = 10.;
unsigned defaultMaxDeltaPoverP = 1.e+3;
double defaultMinPt = 5.;

selTrackExtrapolation::selTrackExtrapolation(const reco::TransientTrack& transientTrack, const AlgebraicVector3& refPoint)
  : errorFlag_(0)
{
  //std::cout << "<selTrackExtrapolation::selTrackExtrapolation>:" << std::endl;

  //std::cout << "refPoint:" << std::endl;
  //std::cout << " x = " << refPoint.At(0) << ", y = " << refPoint.At(1) << ", z = " << refPoint.At(2) << std::endl;

//--- compute point of closest approach of track to reference point
  GlobalPoint refPoint_global(refPoint.At(0), refPoint.At(1), refPoint.At(2));

  TrajectoryStateClosestToPoint dcaPosition = transientTrack.trajectoryStateClosestToPoint(refPoint_global);
  if ( TMath::IsNaN(dcaPosition.position().x()) ||
       TMath::IsNaN(dcaPosition.position().y()) || 
       TMath::IsNaN(dcaPosition.position().z()) ) {
    edm::LogWarning ("selTrackExtrapolation") 
      << " Failed to extrapolate track: Pt = " << transientTrack.track().pt() << "," 
      << " eta = " << transientTrack.track().eta() << ", phi = " << transientTrack.track().phi()*180./TMath::Pi() 
      << " --> skipping !!";
    errorFlag_ = 1;
  }

  dcaPosition_(0) = dcaPosition.position().x();
  dcaPosition_(1) = dcaPosition.position().y();
  dcaPosition_(2) = dcaPosition.position().z();
  //std::cout << "dcaPosition:" << std::endl;
  //std::cout << " x = " << dcaPosition_.At(0) << ", y = " << dcaPosition_.At(1) << ", z = " << dcaPosition_.At(2) << std::endl;

  tangent_(0) = dcaPosition.momentum().x();
  tangent_(1) = dcaPosition.momentum().y();
  tangent_(2) = dcaPosition.momentum().z();
  //std::cout << "tangent:" << std::endl;
  //std::cout << " x = " << tangent_.At(0) << ", y = " << tangent_.At(1) << ", z = " << tangent_.At(2) << std::endl;
  
  TVector3 tangent_vector(tangent_.At(0), tangent_.At(1), tangent_.At(2));

  TRotation rotation;
  rotation.SetZAxis(tangent_vector);
  TRotation invRotation = rotation.Inverse();

  TVector3 test = invRotation*tangent_vector;
  //std::cout << "test:" << std::endl;
  //std::cout << " x = " << test.x() << ", y = " << test.y() << ", z = " << test.z() << std::endl;
  const double epsilon = 1.e-5;
  assert(TMath::Abs(test.x()) < epsilon && TMath::Abs(test.y()) < epsilon);
    
  invRotationMatrix_(0, 0) = invRotation.XX();
  invRotationMatrix_(0, 1) = invRotation.XY();
  invRotationMatrix_(0, 2) = invRotation.XZ();
  invRotationMatrix_(1, 0) = invRotation.YX();
  invRotationMatrix_(1, 1) = invRotation.YY();
  invRotationMatrix_(1, 2) = invRotation.YZ();
  invRotationMatrix_(2, 0) = invRotation.ZX();
  invRotationMatrix_(2, 1) = invRotation.ZY();
  invRotationMatrix_(2, 2) = invRotation.ZZ();
  //std::cout << "invRotationMatrix:" << std::endl;
  //invRotationMatrix_.Print(std::cout);
  //std::cout << std::endl;

//--- compute covariance matrix in rotated coordinates:
//     return x^T Vxx^-1 x = x^T R^-1 R Vxx^-1 R^-1   R x  // R^-1 = R^T (rotation matrices are orthogonal)
//                         = x^T R^T  R Vxx^-1 R^-1   R x  // (R x)^T = x^T R^T
//                         = (R x)^T (R^-1 Vxx R)^-1 (R x) // y := R x 
//                         = y^T (R^T Vxx R)^-1 y          // Vyy := R^T Vxx R
//                         = y^T Vyy^-1 y 
  AlgebraicMatrix33 covMatrix = dcaPosition.theState().cartesianError().position().matrix_new();
  //std::cout << "covMatrix:" << std::endl;
  //covMatrix.Print(std::cout);
  //std::cout << std::endl;
  
  rotCovMatrix_ = ROOT::Math::Transpose(invRotationMatrix_)*covMatrix*invRotationMatrix_;  
  //std::cout << "covMatrix in rotated coordinates:" << std::endl;
  //rotCovMatrix_.Print(std::cout);
  //std::cout << std::endl;

  for ( unsigned iRow = 0; iRow < 2; ++iRow ) {
    for ( unsigned iColumn = 0; iColumn < 2; ++iColumn ) {
      rotCovMatrix2_(iRow, iColumn) = rotCovMatrix_(iRow, iColumn);
    }
  }
}

double selTrackExtrapolation::logLikelihood(const AlgebraicVector3& displacement) const
{
  //std::cout << "<selTrackExtrapolation::logLikelihood>:" << std::endl;

  if ( errorFlag_ ) return 0.;

  AlgebraicVector3 rotDisplacement = invRotationMatrix_*displacement;
  //std::cout << "displacement in rotated coordinates:" << std::endl;
  //rotDisplacement.Print(std::cout);
  //std::cout << std::endl;

  AlgebraicVector2 rotDisplacement2;
  rotDisplacement2(0) = rotDisplacement(0); 
  rotDisplacement2(1) = rotDisplacement(1); 

  double logLikelihood = logGaussianNd(rotDisplacement2, rotCovMatrix2_);

//--- add "penalty" term in case displacement has component opposite to direction of track momentum 
  if ( rotDisplacement(2) < 0. ) {
    //double penaltyTerm = -0.5*square(rotDisplacement(2)/TMath::Sqrt(rotCovMatrix_(2, 2)));
    double penaltyTerm = -square(10.*rotDisplacement(2));
    //std::cout << "displacement in track direction is negative" 
    //	        << " --> adding 'penalty' term  = " << penaltyTerm << " to likelihood !!" << std::endl;
    logLikelihood += penaltyTerm;
  }
 
  return logLikelihood;
}

//
//-------------------------------------------------------------------------------
//

template <typename T>
SVfitLegLikelihoodTrackInfo<T>::SVfitLegLikelihoodTrackInfo(const edm::ParameterSet& cfg)
  : SVfitLegLikelihoodBase<T>(cfg),
    trackBuilder_(0)
{
  //std::cout << "<SVfitLegLikelihoodTrackInfo::SVfitLegLikelihoodTrackInfo>:" << std::endl;

  minNumHits_ = ( cfg.exists("minNumHits") ) ? cfg.getParameter<unsigned>("minNumHits") : defaultMinNumHits;
  minNumPixelHits_ = ( cfg.exists("minNumPixelHits") ) ? cfg.getParameter<unsigned>("minNumPixelHits") : defaultMinNumPixelHits;
  maxChi2DoF_ = ( cfg.exists("maxChi2DoF") ) ? cfg.getParameter<double>("maxChi2DoF") : defaultMaxChi2DoF;
  maxDeltaPoverP_ = ( cfg.exists("maxDeltaPoverP") ) ? cfg.getParameter<double>("maxDeltaPoverP") : defaultMaxDeltaPoverP;
  minPt_ = ( cfg.exists("minPt") ) ? cfg.getParameter<double>("minPt") : defaultMinPt;

  useLinearApprox_ = ( cfg.exists("useLinearApprox") ) ? cfg.getParameter<bool>("useLinearApprox") : true;
  //std::cout << " useLinearApprox = " << useLinearApprox_ << std::endl;
}

template <typename T>
SVfitLegLikelihoodTrackInfo<T>::~SVfitLegLikelihoodTrackInfo()
{
// nothing to be done yet...
}

template <typename T>
void SVfitLegLikelihoodTrackInfo<T>::beginEvent(const edm::Event& evt, const edm::EventSetup& es)
{
//--- get pointer to TransientTrackBuilder
  edm::ESHandle<TransientTrackBuilder> trackBuilderHandle;
  es.get<TransientTrackRecord>().get("TransientTrackBuilder", trackBuilderHandle);
  trackBuilder_ = trackBuilderHandle.product();
  if ( !trackBuilder_ ) {
    edm::LogError ("SVfitLegLikelihoodTrackInfo::beginEvent") 
      << " Failed to access TransientTrackBuilder !!";
  }
}

template <typename T>
void SVfitLegLikelihoodTrackInfo<T>::beginCandidate(const T& leg)
{
  //std::cout << "<SVfitLegLikelihoodTrackInfo::beginCandidate>:" << std::endl;
  //std::cout << " trackBuilder = " << trackBuilder_ << std::endl;

  selectedTracks_.clear();

  std::vector<reco::TrackBaseRef> legTracks = trackExtractor_(leg);
  for ( std::vector<reco::TrackBaseRef>::const_iterator track = legTracks.begin();
	track != legTracks.end(); ++track ) {
    const reco::HitPattern& trackHitPattern = (*track)->hitPattern();
    if ( trackHitPattern.numberOfValidTrackerHits() >= (int)minNumHits_ &&
	 trackHitPattern.numberOfValidPixelHits() >= (int)minNumPixelHits_ &&
	 (*track)->normalizedChi2() < maxChi2DoF_ &&
	 ((*track)->ptError()/(*track)->pt()) < maxDeltaPoverP_ &&
	 (*track)->pt() > minPt_ ) {
      reco::TransientTrack transientTrack = trackBuilder_->build(track->get());
      selectedTracks_.push_back(transientTrack);
    }
  }
  
  isNewCandidate_ = true;
}

template <typename T>
bool SVfitLegLikelihoodTrackInfo<T>::isFittedParameter(int legIndex, int parIndex) const
{
  if ( selectedTracks_.size() > 0 && ((legIndex == SVfit_namespace::kLeg1 && parIndex == SVfit_namespace::kLeg1phiLab              ) ||
				      (legIndex == SVfit_namespace::kLeg1 && parIndex == SVfit_namespace::kLeg1sqrtDecayDistanceLab) ||
				      (legIndex == SVfit_namespace::kLeg2 && parIndex == SVfit_namespace::kLeg2phiLab              ) ||
				      (legIndex == SVfit_namespace::kLeg2 && parIndex == SVfit_namespace::kLeg2sqrtDecayDistanceLab)) )
    return true;
  else 
    return SVfitLegLikelihoodBase<T>::isFittedParameter(legIndex, parIndex);
}

double compScalarProduct(const AlgebraicVector3& v1, const AlgebraicVector3& v2)
{
  return v1.At(0)*v2.At(0) + v1.At(1)*v2.At(1) + v1.At(2)*v2.At(2);
}

double compNorm(const AlgebraicVector3& v)
{
  return TMath::Sqrt(square(v.At(0)) + square(v.At(1)) + square(v.At(2)));
}

template <typename T>
double SVfitLegLikelihoodTrackInfo<T>::operator()(const T& leg, const SVfitLegSolution& solution) const
{
//--- compute negative log-likelihood for tracks of tau lepton decay "leg"
//    to be compatible with originating from hypothetic secondary (tau lepton decay) vertex
//
//    The likelihood is computed as the product of probabilities for the tracks 
//    to be compatible with the hypothetic secondary vertex of the tau lepton decay
//   (distance of closest approach of track to secondary vertex divided by estimated uncertainties of track extrapolation)

  //std::cout << "<SVfitLegLikelihoodTrackInfo::operator()>:" << std::endl;

  double logLikelihood = 0.;

  if ( useLinearApprox_ ) {

//--- "linearize" computation of point of closest approach of track to hypothetic seconday (tau decay) vertex,
//    in order to improve numerical convergence of Minuit fit;
//    compute track direction and covariance matrix for one reference point
//   ("original" primary event (tau production) vertex position 
//    plus expected "mean" tau lepton flight path, given by (visibleEnergy/tauMass)*c*tauLifetime)
//    instead of doing a "full" helix track extrapolation for each hypothetic seconday (tau decay) vertex position;
//    given the "linearized" track direction vector, approximate point of closest approach by:
//
//     primaryVertex + trackVector*scalarProduct(trackVector, secondaryVertex)/(|trackVector|*|secondaryVertex|)
//
    if ( isNewCandidate_ ) {
      //std::cout << "--> computing linear approximation of helix track extrapolation..." << std::endl;
      
      double decayDistance0 = (leg.energy()/SVfit_namespace::tauLeptonMass)*SVfit_namespace::cTauLifetime;

      AlgebraicVector3 direction;
      direction(0) = leg.px()/leg.p();
      direction(1) = leg.py()/leg.p();
      direction(2) = leg.pz()/leg.p();
      
      AlgebraicVector3 refPoint = pvPosition_ + decayDistance0*direction;
      
      selectedTrackInfo_.clear();
      for ( std::vector<reco::TransientTrack>::const_iterator selectedTrack = selectedTracks_.begin();
	    selectedTrack != selectedTracks_.end(); ++selectedTrack ) {
	selectedTrackInfo_.push_back(selTrackExtrapolation(*selectedTrack, refPoint));
      }

      isNewCandidate_ = false;
    }

    //std::cout << "--> computing distance between extrapolated track and (hypothetic) tau decay vertex..." << std::endl;
    for ( std::vector<selTrackExtrapolation>::const_iterator selectedTrackInfo = selectedTrackInfo_.begin();
	  selectedTrackInfo != selectedTrackInfo_.end(); ++selectedTrackInfo ) {
      AlgebraicVector3 relDecayVertexPos = solution.decayVertexPos() - selectedTrackInfo->refPoint();

      double projection = compScalarProduct(relDecayVertexPos, selectedTrackInfo->tangent());
      projection /= compNorm(relDecayVertexPos);
      projection /= compNorm(selectedTrackInfo->tangent());
      
      AlgebraicVector3 displacement = relDecayVertexPos - projection*selectedTrackInfo->tangent();
      
      logLikelihood += selectedTrackInfo->logLikelihood(displacement);
    }
  } else {
    //std::cout << "--> computing distance to (hypothetic) tau decay vertex using full helix extrapolation of track..." << std::endl;
    for ( std::vector<reco::TransientTrack>::const_iterator selectedTrack = selectedTracks_.begin();
	  selectedTrack != selectedTracks_.end(); ++selectedTrack ) {
      selTrackExtrapolation selectedTrackInfo(*selectedTrack, solution.decayVertexPos());
      
      AlgebraicVector3 displacement = solution.decayVertexPos() - selectedTrackInfo.dcaPosition();
      
      logLikelihood += selectedTrackInfo.logLikelihood(displacement);
    }
  }
  
  //std::cout << "--> -log(likelihood) = " << -logLikelihood << std::endl;
  
  return -logLikelihood;
}

#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/Candidate/interface/Candidate.h"

typedef SVfitLegLikelihoodTrackInfo<pat::Electron> SVfitElectronLikelihoodTrackInfo;
typedef SVfitLegLikelihoodTrackInfo<pat::Muon> SVfitMuonLikelihoodTrackInfo;
typedef SVfitLegLikelihoodTrackInfo<pat::Tau> SVfitTauLikelihoodTrackInfo;
typedef SVfitLegLikelihoodTrackInfo<reco::Candidate> SVfitCandidateLikelihoodTrackInfo;

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_EDM_PLUGIN(SVfitElectronLikelihoodBasePluginFactory, SVfitElectronLikelihoodTrackInfo, "SVfitElectronLikelihoodTrackInfo");
DEFINE_EDM_PLUGIN(SVfitMuonLikelihoodBasePluginFactory, SVfitMuonLikelihoodTrackInfo, "SVfitMuonLikelihoodTrackInfo");
DEFINE_EDM_PLUGIN(SVfitTauLikelihoodBasePluginFactory, SVfitTauLikelihoodTrackInfo, "SVfitTauLikelihoodTrackInfo");
DEFINE_EDM_PLUGIN(SVfitCandidateLikelihoodBasePluginFactory, SVfitCandidateLikelihoodTrackInfo, "SVfitCandidateLikelihoodTrackInfo");
