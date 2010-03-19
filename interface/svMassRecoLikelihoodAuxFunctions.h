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
  const double chargedPionMass = 0.13957;
  const double muonMass = 0.105658;
  const double electronMass = 5.109989e-4;

  /// Tired of all these namespaces
  typedef reco::Candidate::LorentzVector FourVector;
  typedef std::pair<FourVector, FourVector> FourVectorPair;
  typedef reco::Candidate::Vector ThreeVector;

  /**************************************************************
   * Likelihoods for the various constraints used in the fitter *
   **************************************************************/

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
  FourVectorPair compInvisibleLeg(const ThreeVector& motherDirection,
                                  const FourVector& visibleP4, const double massMother, 
                                  const double m12Squared, int& error);

  /// Compute the upper limit on the neutrino system mass^2 given the tau direction
  /// and visible momentum
  double m12SquaredUpperBound(const FourVector& visP4, const ThreeVector& tauDir);

  /**************************************************************
   * Helper functions for type dependent functionality          *
   **************************************************************/

///--- Index the different types used
  template<typename T> int legTypeLabel(const T& leg) { return -30; }
  template<> inline int legTypeLabel<pat::Electron>(const pat::Electron& leg) { return -2; }
  template<> inline int legTypeLabel<pat::Muon>(const pat::Muon& leg) { return -1; }
  template<> inline int legTypeLabel<pat::Tau>(const pat::Tau& leg) { return leg.decayMode(); }

///--- Get the tracks associated to the leg
  template<typename T> std::vector<reco::TrackBaseRef> getTracks(const T& leg) {
     // default case
     return std::vector<reco::TrackBaseRef>();
  }
  // Get tracks from tau
  template<> inline std::vector<reco::TrackBaseRef> getTracks<pat::Tau>(const pat::Tau& tau) {
     std::vector<reco::TrackBaseRef> output;
    const reco::PFCandidateRefVector& signalCHs = tau.signalPFChargedHadrCands();
    for ( size_t itrk = 0; itrk < signalCHs.size(); itrk++ ) {
      output.push_back(reco::TrackBaseRef(signalCHs.at(itrk)->trackRef()));
    }
    return output;
  }
  // Get track from muon
  template<> inline std::vector<reco::TrackBaseRef> getTracks<pat::Muon>(const pat::Muon& muon) {
     std::vector<reco::TrackBaseRef> output;
    output.push_back(reco::TrackBaseRef(muon.track())); //track() returns only inner track
    return output;
  }
  // Get track from electron
  template<> inline std::vector<reco::TrackBaseRef> getTracks<pat::Electron>(const pat::Electron& elec) {
     std::vector<reco::TrackBaseRef> output;
    output.push_back(reco::TrackBaseRef(elec.gsfTrack())); // Track or GSF track?
    return output;
  }

///--- Get the hypothesis for the mass of a charged constituent
  template<typename T> double chargedMass2ByType() { return -1; } //default
  // Hadronic case (pion)
  template<> inline double chargedMass2ByType<pat::Tau>() { return chargedPionMass*chargedPionMass; }
  // Muonic case
  template<> inline double chargedMass2ByType<pat::Muon>() { return muonMass*muonMass; }
  // Electronic case
  template<> inline double chargedMass2ByType<pat::Electron>() { return electronMass*electronMass; }

  
///--- Get the visible neutral components, if they exist
  template<typename T> FourVector getNeutralP4(const T& object) { return FourVector(); }
  /// Only the tau has associated neutral objects
  template<> inline FourVector getNeutralP4<pat::Tau>(const pat::Tau& tau) {
    FourVector output;
    const reco::PFCandidateRefVector& signalGammas = tau.signalPFGammaCands();
    for ( size_t gamma = 0; gamma < signalGammas.size(); ++gamma ) {
      output += signalGammas[gamma]->p4();
    }
    return output;
  }

///--- Determine is the leg's invisible component can have mass (i.e. one or two neutrinos)
  template <typename T> bool nuSystemIsMassless() { return true; }
  template<> inline bool nuSystemIsMassless<pat::Tau>() { return true;}
  template<> inline bool nuSystemIsMassless<pat::Electron>() { return false;}
  template<> inline bool nuSystemIsMassless<pat::Muon>() { return false;}

///--- Helper function to determine if the given type is supported
  template<typename T> bool typeIsSupportedBySVFitter()  { return false; }
  template<> inline bool typeIsSupportedBySVFitter<pat::Tau>() { return true; }
  template<> inline bool typeIsSupportedBySVFitter<pat::Muon>() { return true; }
  template<> inline bool typeIsSupportedBySVFitter<pat::Electron>() { return true; }
}

#endif

