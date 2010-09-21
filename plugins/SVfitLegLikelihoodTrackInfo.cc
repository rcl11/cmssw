#include "TauAnalysis/CandidateTools/plugins/SVfitLegLikelihoodTrackInfo.h"

#include "TauAnalysis/CandidateTools/interface/svFitAuxFunctions.h"

#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateClosestToPoint.h"

#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/CLHEP/interface/AlgebraicObjects.h"

#include <TRotation.h>

#include <limits>

using namespace SVfit_namespace;

unsigned defaultMinNumHits = 5;
unsigned defaultMinNumPixelHits = 1;
double defaultMaxChi2DoF = 10.;
unsigned defaultMaxDeltaPoverP = 1.e+3;
double defaultMinPt = 5.;

template <typename T>
SVfitLegLikelihoodTrackInfo<T>::SVfitLegLikelihoodTrackInfo(const edm::ParameterSet& cfg)
  : SVfitLegLikelihoodBase<T>(cfg),
    trackBuilder_(0)
{
  minNumHits_ = ( cfg.exists("minNumHits") ) ? cfg.getParameter<unsigned>("minNumHits") : defaultMinNumHits;
  minNumPixelHits_ = ( cfg.exists("minNumPixelHits") ) ? cfg.getParameter<unsigned>("minNumPixelHits") : defaultMinNumPixelHits;
  maxChi2DoF_ = ( cfg.exists("maxChi2DoF") ) ? cfg.getParameter<double>("maxChi2DoF") : defaultMaxChi2DoF;
  maxDeltaPoverP_ = ( cfg.exists("maxDeltaPoverP") ) ? cfg.getParameter<double>("maxDeltaPoverP") : defaultMaxDeltaPoverP;
  minPt_ = ( cfg.exists("minPt") ) ? cfg.getParameter<double>("minPt") : defaultMinPt;
}

template <typename T>
SVfitLegLikelihoodTrackInfo<T>::~SVfitLegLikelihoodTrackInfo()
{
// nothing to be done yet...
}

template <typename T>
void SVfitLegLikelihoodTrackInfo<T>::beginEvent(const edm::Event& evt, const edm::EventSetup& es)
{
  std::cout << "<SVfitLegLikelihoodTrackInfo::beginEvent>:" << std::endl;

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
  std::cout << "<SVfitLegLikelihoodTrackInfo::beginCandidate>:" << std::endl;
  std::cout << " trackBuilder = " << trackBuilder_ << std::endl;

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
}

template <typename T>
bool SVfitLegLikelihoodTrackInfo<T>::isFittedParameter(int legIndex, int parIndex) const
{
  if ( selectedTracks_.size() > 0 && ((legIndex == SVfit_namespace::kLeg1 && parIndex == SVfit_namespace::kLeg1phiLab       ) ||
				      (legIndex == SVfit_namespace::kLeg1 && parIndex == SVfit_namespace::kLeg1flightPathLab) ||
				      (legIndex == SVfit_namespace::kLeg2 && parIndex == SVfit_namespace::kLeg2phiLab       ) ||
				      (legIndex == SVfit_namespace::kLeg2 && parIndex == SVfit_namespace::kLeg2flightPathLab)) )
    return true;
  else 
    return SVfitLegLikelihoodBase<T>::isFittedParameter(legIndex, parIndex);
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

  GlobalPoint svPosition(solution.decayVertexPos().At(0), solution.decayVertexPos().At(1), solution.decayVertexPos().At(2));

  double logLikelihood = 0.;

  for ( std::vector<reco::TransientTrack>::const_iterator track = selectedTracks_.begin();
	track != selectedTracks_.end(); ++track ) {

//--- compute point of closest approach of track to secondary (tau lepton dercay) vertex
    TrajectoryStateClosestToPoint dcaPosition = track->trajectoryStateClosestToPoint(svPosition);

//--- compute distance between point of closest approach and secondary (tau lepton dercay) vertex position
    AlgebraicVector3 displacement3(dcaPosition.position().x() - svPosition.x(), 
				   dcaPosition.position().y() - svPosition.y(), 
				   dcaPosition.position().z() - svPosition.z());
    std::cout << "displacement3:" << std::endl;
    displacement3.Print(std::cout);

//--- compute uncertainty on extrapolated track position
//    projected onto plane normal to track direction 
//   (at point of closest approach)
    AlgebraicMatrix33 displacement3Err(dcaPosition.theState().cartesianError().position().matrix_new());
    std::cout << "displacement3Err:" << std::endl;
    displacement3Err.Print(std::cout);

    TVector3 tangent(dcaPosition.momentum().x(), dcaPosition.momentum().y(), dcaPosition.momentum().z());
    TRotation rot;
    rot.SetZAxis(tangent);
    TRotation rotInverse = rot.Inverse();
    
    AlgebraicMatrix33 rotInverseMatrix;
    rotInverseMatrix(0, 0) = rotInverse.XX();
    rotInverseMatrix(0, 1) = rotInverse.XY();
    rotInverseMatrix(0, 2) = rotInverse.XZ();
    rotInverseMatrix(1, 0) = rotInverse.YX();
    rotInverseMatrix(1, 1) = rotInverse.YY();
    rotInverseMatrix(1, 2) = rotInverse.YZ();
    rotInverseMatrix(2, 0) = rotInverse.ZX();
    rotInverseMatrix(2, 1) = rotInverse.ZY();
    rotInverseMatrix(2, 2) = rotInverse.ZZ();

    displacement3 = rotInverseMatrix*displacement3;
    std::cout << "displacement in rotated coordinates:" << std::endl;
    displacement3.Print(std::cout);

    AlgebraicMatrix33 rotInverseMatrixT(rotInverseMatrix); // rotation matrices are orthogonal, i.e. rotInverse^T = rotInverse

    displacement3Err = rotInverseMatrix*displacement3Err*rotInverseMatrixT;

    std::cout << "displacement3Err in rotated coordinates:" << std::endl;
    displacement3Err.Print(std::cout);

    AlgebraicVector2 displacement2;
    displacement2(0) = displacement3(0); 
    displacement2(1) = displacement3(1); 

    AlgebraicMatrix22 displacement2Err;
    for ( unsigned iRow = 0; iRow < 2; ++iRow ) {
      for ( unsigned iColumn = 0; iColumn < 2; ++iColumn ) {
	displacement2Err(iRow, iColumn) = displacement3Err(iRow, iColumn);
      }
    }

    logLikelihood += logGaussianNd(displacement2, displacement2Err);
  }

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
