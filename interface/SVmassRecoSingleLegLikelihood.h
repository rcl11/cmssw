#ifndef TauAnalysis_CandidateTools_SVmassRecoSingleLegLikelihood_h
#define TauAnalysis_CandidateTools_SVmassRecoSingleLegLikelihood_h

/* 
 * svMassReco::SVmassRecoSingleLegLikelihood
 *
 * Author: Evan K. Friis, UC Davis
 *
 * Class that computes the negative log likelihood for a 'leg' in a ditau candidate.
 * All pat::Electron, pat::Muon and pat::Tau specific code to access variables used in the fit
 * is encapsulated in an object of type SVmassRecoSingleLegExtractorBase passed
 * to SVmassRecoSingleLegLikelihood object in the constructor.
 *
 */

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

#include "TauAnalysis/CandidateTools/interface/svMassRecoLikelihoodAuxFunctions.h"
#include "TauAnalysis/CandidateTools/interface/SVmassRecoSingleLegExtractorBase.h"
#include "TauAnalysis/CandidateTools/interface/SVFitTrackPositionFinder.h"

#include "TrackingTools/GeomPropagators/interface/AnalyticalTrajectoryExtrapolatorToLine.h"

#include <memory>
#include <TVector3.h>

namespace svMassReco 
{
  
  class SVmassRecoSingleLegLikelihood
  {
   public:
    SVmassRecoSingleLegLikelihood(const SVmassRecoSingleLegExtractorBase& extractor, 
				  const std::vector<reco::TransientTrack>& tracks);

    virtual ~SVmassRecoSingleLegLikelihood() {};

    /// Set points to determine NLL 
    void setPoints(const GlobalPoint& pv, double thetaRest, 
          double phiLab, double radiusLab, double m12); 

    /// Do a naive opimization
    void findIntialValues(const GlobalPoint& pv, 
          double &thetaGuess, double &thetaError, 
          double &phiGuess, double &phiError,
          double &radiusGuess, double &radiusError,
          double &m12Guess, double &m12Error);

    /// Build the direction of tau given the fit parameters 
    GlobalVector buildTauDirection(double thetaRest, double phiLab, double m12);

    /// Given a value for the neutrino system mass, determine the tau rest frame
    /// visible momentum
    double findPGivenM12(double M12) const;

    /// Get the current fitted mass of neutrinos
    double m12() const { return M12_; };

    /// Get the current visible decay angle in the rest frame
    double thetaRest() const { return thetaRest_; };

    /// Get the maximum possible neutrino system mass
    double maxM12() const;

    /// Total NLL for this leg
    double nllOfLeg() const; 

    /// NLL of the kinematics in the rest frame of the tau
    double nllRestFrameKinematics() const; 

    /// NLL for this leg from the decay length constraint.
    double nllDecayLength() const;

    /// NLL for SV given tracker measurements
    double nllTopological() const; 

    /// Secondary vertex associated with this leg
    const GlobalPoint& sv() const { return sv_; };

    /// Inferred tau diretion and decay length
    const ThreeVector& dir() const { return legDir_; };

    /// P4 of this visible part of this leg
    const FourVector& visP4() const { return visP4_; };

    /// P4 of this invisible part of this leg
    const FourVector& nuP4() const { return nuP4_; };

    /// P4 of this leg
    const FourVector& fittedP4() const { return p4_; };

    /// Method to get the type of Leg 
    int legType() const { return extractor_.legTypeLabel(); }

    bool nuSystemIsMassless() const { return extractor_.nuSystemIsMassless(); }

    friend std::ostream& operator<< (std::ostream &out, const SVmassRecoSingleLegLikelihood& fit) 
    { 
      fit.printTo(out); 
      return out; 
    };

    /// Pretty print leg information
    void printTo(std::ostream &out) const;

  private:
    // Leg object
    const SVmassRecoSingleLegExtractorBase& extractor_;
    
    /// The associated tracks
    const std::vector<reco::TransientTrack>& tracks_;

    AnalyticalTrajectoryExtrapolatorToLine propagator_;

    // Visible quantities that are used a lot
    const FourVector visP4_;
    const ThreeVector visP3_;
    const GlobalVector visP3GlobalVector_; // let's only do this once
    const TVector3 visPDirection_;
    const double visibleMass_;

    /// The track with the best fit
    reco::TransientTrack bestTrack_;
    /// Fitted vertex, if a three prong tau.  Null if unfilled
    TransientVertex vertex_;


    /* Fit parameters */
    GlobalPoint pv_;
    // Mass of two neutrino system
    double M12_;
    // Current rest frame theta value in fit
    double thetaRest_;
    // Current phi value in fit
    double phiLab_;
    // Current radius value in fit
    double radiusLab_;

    // The secondary vertex
    GlobalPoint sv_;
    // Inferred direction of tau lepton
    ThreeVector legDir_;
    // Total fitted p4
    FourVector p4_;
    // Nu p4
    FourVector nuP4_;     


  };

}

#endif
