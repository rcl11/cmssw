#ifndef TauAnalysis_CandidateTools_SVLegFitter_h
#define TauAnalysis_CandidateTools_SVLegFitter_h
/* 
 * TauVertex::LegFitter
 *
 * Author: Evan K. Friis, UC Davis
 *
 * Class that computes the negative log likelihood for a 'leg' in a ditau candidate.
 * Main functionality is held in the abstract base class LegFitter, with specific implementations
 * for type T = pat::Muons, Electrons, and Taus defined in LegFitterSpecific<T>
 *
 */

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

#include "TMinuit.h"

// Load the actual NLL distributions
#include "TauAnalysis/CandidateTools/interface/SVMethodLikelihoods.h"

namespace TauVertex {

class LegFitter
{
   public:
      LegFitter(const std::vector<reco::TransientTrack>& tracks):
         tracks_(tracks),tscps_(std::vector<TrajectoryStateClosestToPoint>(tracks.size())),nTracks_(tracks.size()){};
      virtual ~LegFitter() {};

      /// Set points to determine NLL 
      void setPoints(const GlobalPoint& pv, double x, double y, double z, double m12scale, int& error);

      /// Total NLL for this leg
      double nllOfLeg() const { return 
            nllTopological() + 
            nllRapidity() + 
            nllDecayLength() +
            nllM12Penalty(); 
      };

      /// NLL to keep the fit physical
      double nllM12Penalty() const;

      /// nll for this leg from the decay length constraint.
      double nllDecayLength() const; 

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

      /// NLL for SV given tracker measurements
      double nllTopological() const;

      /// Rapidity of visible w.r.t tau direction
      double visRapidity() const;

      /// Uncorrected visible P4 (i.e. straight from the pat::Tau, etc)
      virtual FourVector uncorrectedP4() const = 0;

      /* Abstract functions specified for each leg type */
      /// NLL for rapidity given momentum
      virtual double nllRapidity() const = 0;

      /// Get the total four momentum of the tracks at the point closes to the SV
      FourVector visChargedP4() const;

      /// Get the neutral visible p4, specific to each data type
      virtual FourVector visNeutralP4() const = 0;

      /// Abstract method to assign correct mass hypothesis to tracks
      virtual FourVector chargedP4FromMomentum(const GlobalVector& p3) const = 0;

      /// Method to get the type of Leg 
      virtual int legType() const = 0;

      /// Access to tracks
      const std::vector<reco::TransientTrack>& tracks() const {return tracks_;};

      friend std::ostream& operator<< (std::ostream &out, LegFitter &fit) { fit.printTo(out); return out; };

      virtual void printTo(std::ostream &out) const;

   protected:

      /// Total visible p4
      FourVector fitVisP4() const { return visNeutralP4() + visChargedP4(); };

      /// Fit neutrino momentum 
      virtual TauVertex::FourVector fitNuP4(double m12scale, int& error) const = 0;

      /// The associated tracks
      const std::vector<reco::TransientTrack>& tracks_;
      /// The trajectory states closes to the secondary vertex
      std::vector<TrajectoryStateClosestToPoint> tscps_;
      size_t nTracks_;

   private:
      // Fit parameters
      GlobalPoint sv_;
      ThreeVector legDir_;
      // Total fitted p4
      TauVertex::FourVector p4_;
      // Visible p4
      TauVertex::FourVector visP4_;
      // Nu p4
      TauVertex::FourVector nuP4_;
         
};

/// HELPER FUNCTIONS

/// Get the tracks from a given type
template<typename T> int legTypeLabel(const T& leg);

/// Get the tracks from a given type
template<typename T> std::vector<reco::TrackBaseRef> getTracks(const T& leg);

/// Template to determine the mass hypothesis for tracks for the various decay types (muon/electon/pion)
template<typename T> double chargedMass2ByType();

/* 
 * Helper template used to distinguish between leptonic decays and hadronic
 * decays.  In leptonic decays, the 2 neutrino system can have a mass, which
 * must be fitted.  For 1-nu hadronic decays, this always zero.  For e/mu
 * decays, it is defined by m12^2 = mTau^2 + mLepton^2 - 2 mTau E_rest, where
 * E_rest is the energy of the lepton in tau rest frame.  E_rest is bounded
 * from below by sqrt(mLepton^2 + pLeptonPerp^2)
*/
double 
m12SquaredUpperBound(const FourVector& visP4, const ThreeVector& tauDir);

/*** Utility function to determine whether Minuit should allow the scale parameter to float */
// Default case is false - allow it to float. 
template <typename T> bool fixM12ScaleParameter();

/// Template to get the neutral components of a compound object
template<typename T> FourVector getNeutralP4(const T& object);

// Type specific implementation
template<typename T>
class LegFitterSpecific : public LegFitter
{
   public: 
      LegFitterSpecific(const T& object, const std::vector<reco::TransientTrack>& tracks, bool ansatzForward):
         LegFitter(tracks),object_(object),ansatzForward_(ansatzForward){};
      virtual ~LegFitterSpecific(){};

      /// Return the appropriate NLL for the vis rapidity
      double nllRapidity() const { 
         double nll = nllVisRapidityGivenMomentum<T>(object_, this->visRapidity(), this->visP4().P()); 
         return nll;
      }

      /// Helper function to build fourvector from momentum
      FourVector chargedP4FromMomentum(const GlobalVector& p) const {
         return FourVector(p.x(), p.y(), p.z(), sqrt(p.mag2() + chargedMass2ByType<T>()));
      }

      FourVector visNeutralP4() const { return getNeutralP4<T>(object_); };

      /// Raw p4 of the object
      FourVector uncorrectedP4() const { return object_.p4(); };

      /// Log helper function
      int legType() const { return legTypeLabel<T>(object_); };

   protected:
      /// Total fit of P4 (including neutrino)
      FourVector fitNuP4(double m12scale, int& error) const 
      { 
         double m12SquaredUpperBoundVal = m12SquaredUpperBound(this->visP4(), this->dir());
         double m12Squared = m12scale*m12scale*m12SquaredUpperBoundVal;
         FourVectorPair solutions = 
            compInvisibleLeg(this->dir(), this->visP4(), tauMass, m12Squared, error);
         // Determine which solution to take
         if(ansatzForward_)
            return solutions.first;
         return solutions.second;
      }
   private:
      const T& object_;
      bool ansatzForward_;
};
}
#endif
