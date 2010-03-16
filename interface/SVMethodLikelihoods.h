#ifndef TauAnalysis_RecoTools_SVMethodLikelihoods_h
#define TauAnalysis_RecoTools_SVMethodLikelihoods_h

#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
#include "DataFormats/METReco/interface/MET.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateClosestToPoint.h"
#include "DataFormats/Candidate/interface/Candidate.h"

namespace TauVertex {

/// Should use PDT tables...
const double tauMass = 1.77685;
const double tauLifetime = 290.6e-15;
const double chargedPionMass = 0.13957;
const double muonMass = 0.105658;
const double electronMass = 5.109989e-4;

/// Tired of all these namespaces
typedef reco::Candidate::LorentzVector FourVector;
typedef std::pair<FourVector, FourVector> FourVectorPair;
typedef reco::Candidate::Vector ThreeVector;

/// Negative log likelihood for a point to be compatabile with a measured vertex
double nllPointGivenVertex(const GlobalPoint& point, const TransientVertex& vertex);

/// Negative log likelihood for a point to be compatabile with a measured track
double nllPointGivenTrack(const TrajectoryStateClosestToPoint& tcsp);

/// Negative log likelihood for a tau of given energy to have reconstructed decay length
double nllTauDecayLengthGivenMomentum(double length, double momentum);

/// Returns negative log likelihood function for the visible decay products to
// have given rapidity (defined w.r.t. the reconstructed tau direction) given
// the measured momentum.  Determined from MC, for each decay type.  Requires PFTauDecayMode 
// be embedded in the pat::Tau
//template<int decayindex>
//double nllVisRapidityGivenMomentum(double rapidity, double momentum);
template<typename T>
double nllVisRapidityGivenMomentum(const T& obj, double rapidity, double momentum);

/// Compute the two solutions for the missing four vector of a tau decay, given the measured
///  tau direction, visible four vector, and estimated (or exact) mass of the missing energy system.
///  Error returns as 1 if unphysical result
FourVectorPair compInvisibleLeg( const ThreeVector& motherDirection,
         const FourVector& visibleP4, const double massMother, const double m12Squared, int& error);

double nllNuSystemGivenMET(const FourVector& nus, const reco::MET* met); 

}

#endif

