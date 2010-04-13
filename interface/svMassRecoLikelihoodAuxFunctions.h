#ifndef TauAnalysis_CandidateTools_svMassRecoLikelihoodAuxFunctions_h
#define TauAnalysis_CandidateTools_svMassRecoLikelihoodAuxFunctions_h

#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
#include "DataFormats/METReco/interface/MET.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateClosestToPoint.h"
#include "DataFormats/Candidate/interface/Candidate.h"

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Tau.h"

namespace svMassReco {

  /// Should use PDT tables...
  const double tauMass = 1.77685;
  const double tauLifetime = 290.6e-15;

  /// Tired of all these namespaces
  typedef reco::Candidate::LorentzVector FourVector;
  typedef std::pair<FourVector, FourVector> FourVectorPair;
  typedef reco::Candidate::Vector ThreeVector;

  /// inline definitions of auxiliary functions
  inline double square(double x) { return x*x; }
  inline double cube(double x) { return x*x*x; }

  inline double nlGaussianNorm(double sigma, int dimension=1) { 
    // Norm of Gaussian = 
    //  1/(sigma*(2*pi)^(k/2))
    return (log(sigma) + (dimension/2.0)*log(TMath::TwoPi()));
  }

  /**************************************************************
   * Likelihoods for the various constraints used in the fitter *
   **************************************************************/

  /// Negative log likelihood for a given lepton energy in a tau rest frame
  double nllLeptonTauRestFrameEnergy(double m12Squared, double leptonMass);

  /// Negative log likelihood for a point to be compatabile with a measured vertex
  double nllPointGivenVertex(const GlobalPoint& point, const TransientVertex& vertex);

  /// Negative log likelihood for a point to be compatabile with a measured track
  double nllPointGivenTrack(const TrajectoryStateClosestToPoint& tcsp);

  /// Negative log likelihood for a tau of given energy to have reconstructed decay length
  double nllTauDecayLengthGivenMomentum(double length, double momentum);

  /// Returns negative log likelihood function for the visible decay products to
  /// have given rapidity (defined w.r.t. the reconstructed tau direction) given
  /// the measured momentum. Determined from MC, for each decay type.  
  template<typename T>
  double nllVisRapidityGivenMomentum(const T& obj, double rapidity, double momentum);

  /// Compute the negative log likelihood of measured MET given fitted invisible four
  /// vector sum
  double nllNuSystemGivenMET(const FourVector& nus, const reco::MET* met); 

  /**************************************************************
   * Utilities to separate leg tracks from the PV tracks        *
   **************************************************************/
  /// Predicate to key tracks
  template<typename T>
  struct RefToBaseLess : public std::binary_function<edm::RefToBase<T>, edm::RefToBase<T>, bool> 
  {
    inline bool operator()(const edm::RefToBase<T> &r1, const edm::RefToBase<T> &r2) const
    {
      return r1.id() < r2.id() || (r1.id() == r2.id() && r1.key() < r2.key());
    }
  };

  /// auxiliary Function to compare two tracks
  bool tracksAreMatched(const reco::TrackBaseRef&, const reco::TrackBaseRef&); 

  /**************************************************************
   * Functions to geometrically solve for invisible momenta     *
   **************************************************************/

  /// Compute the two solutions for the missing four vector of a tau decay, given the measured
  /// tau direction, visible four vector, and estimated (or exact) mass of the missing energy system.
  /// Error returns as 1 if unphysical result
  FourVectorPair compInvisibleLeg(const ThreeVector&, const FourVector&, const double, const double, int&);

  /// Compute the upper limit on the neutrino system mass^2 given the tau direction
  /// and visible momentum
  double m12SquaredUpperBound(const FourVector&, const ThreeVector&);
}

#endif

