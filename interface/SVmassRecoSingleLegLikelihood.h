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
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

#include "TauAnalysis/CandidateTools/interface/svMassRecoLikelihoodAuxFunctions.h"
#include "TauAnalysis/CandidateTools/interface/SVmassRecoSingleLegExtractorBase.h"

namespace svMassReco 
{

class SVmassRecoSingleLegLikelihood
{
   public:
      SVmassRecoSingleLegLikelihood(const SVmassRecoSingleLegExtractorBase& extractor, 
				    const std::vector<reco::TransientTrack>& tracks, bool ansatzForward)
	: extractor_(extractor),
	  tracks_(tracks),
	  tscps_(std::vector<TrajectoryStateClosestToPoint>(tracks.size())), nTracks_(tracks.size()),
  	  ansatzForward_(ansatzForward)
      {};

      virtual ~SVmassRecoSingleLegLikelihood() {};

      /// Get a valid initial condition for the SV associated to this leg
      std::pair<GlobalPoint, GlobalError> findInitialSecondaryVertex(const GlobalPoint&);

      /// Set points to determine NLL 
      void setPoints(const GlobalPoint& pv, double x, double y, double z, 
            double m12scale, int& error) {
         // secondary vertex for this leg
         sv_ = GlobalPoint(x,y,z); 
         legDir_ = ThreeVector(x-pv.x(), y-pv.y(), z-pv.z());
         // Update the trajectory states closest to the SV for each track
         for ( size_t itrk = 0; itrk < nTracks_; ++itrk ) {
            tscps_[itrk] = tracks_[itrk].trajectoryStateClosestToPoint(sv_);
         }
         /// Update all the kinematic quantities
         visP4_ = fitVisP4();
         nuP4_ = fitNuP4(m12scale, error);
         p4_ = visP4_ + nuP4_;
      }

      /// Total NLL for this leg
      double nllOfLeg() const 
      { 
         return nllTopological() + nllRapidity() + nllDecayLength() + nllM12Penalty(); 
      };

      /// NLL to keep the fit physical
      double nllM12Penalty() const {
         double m12SquaredUpperBoundVal = m12SquaredUpperBound(this->visP4(), this->dir());
         if ( m12SquaredUpperBoundVal < 0 ) {
            return (m12SquaredUpperBoundVal*m12SquaredUpperBoundVal/1e-6);
         } else return 0;
      }

      /// nll for this leg from the decay length constraint.
      double nllDecayLength() const {
         return nllTauDecayLengthGivenMomentum(legDir_.r(), p4_.P());
      }

      /// NLL for SV given tracker measurements
      double nllTopological() const {
         /// TODO does this work in three prong case, instead of vertex fit?
         double output = 0.0;
         for ( size_t itrk = 0; itrk < nTracks_; ++itrk ) {
            output += nllPointGivenTrack(tscps_[itrk]);
         }
         return output;
      }

      /// Return the appropriate NLL for the vis rapidity
      double nllVisRapidityGivenMomentumElectronCase(double, double) const;
      double nllVisRapidityGivenMomentumMuonCase(double, double) const;
      double nllVisRapidityGivenMomentumTauJetCase(int, double, double) const;
      double nllRapidity() const { 
	int legTypeLabel = extractor_.legTypeLabel();
	double rapidity = this->visRapidity();
	double momentum = this->visP4().P(); 
	switch ( legTypeLabel ) { // legTypeLabel values defined in TauAnalysis/CandidateTools/interface/SVmassRecoSingleLegExtractorT.h
	case -30:
	  //std::cout << "<SVmassRecoSingleLegLikelihood::nllRapidity(reco::Candidate)>:" << std::endl;
	  return 0.;
	case -2:
	  //std::cout << "<SVmassRecoSingleLegLikelihood::nllRapidity(pat::Electron)>:" << std::endl;
	  return nllVisRapidityGivenMomentumElectronCase(rapidity, momentum);
	case -1:
	  //std::cout << "<SVmassRecoSingleLegLikelihood::nllRapidity(pat::Muon)>:" << std::endl;
	  return nllVisRapidityGivenMomentumMuonCase(rapidity, momentum);
	default:
	  //std::cout << "<SVmassRecoSingleLegLikelihood::nllRapidity(pat::Tau)>:" << std::endl;
	  return nllVisRapidityGivenMomentumTauJetCase(legTypeLabel, rapidity, momentum);
	}
      }

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

      /// Rapidity of visible w.r.t tau direction
      double visRapidity() const {
         return atanh(visP4_.Vect().Dot(legDir_.unit())/visP4_.e()); 
      }

      /// Uncorrected visible P4 (i.e. straight from the pat::Tau, etc)
      virtual FourVector uncorrectedP4() const { return extractor_.p4(); };

      /// Get the total four momentum of the tracks at the point closes to the SV
      FourVector visChargedP4() const {
         FourVector output;
         for ( size_t itrk = 0; itrk < nTracks_; ++itrk ) {
            GlobalVector trkMomentumAtCP = tscps_[itrk].momentum();
            output += chargedP4FromMomentum(trkMomentumAtCP);
         }
         return output;
      }

      /// Get the neutral visible p4, specific to each data type
      FourVector visNeutralP4() const { 
	return extractor_.getNeutralP4(); 
      }

      /// Helper function to build fourvector from momentum, with the correct mass
      FourVector chargedP4FromMomentum(const GlobalVector& p) const 
      {
         return FourVector(p.x(), p.y(), p.z(), sqrt(p.mag2() + extractor_.chargedMass2()));
      }

      /// Method to get the type of Leg 
      int legType() const { 
	return extractor_.legTypeLabel();
      }

      bool nuSystemIsMassless() const {
	return extractor_.nuSystemIsMassless();
      }

      /// Access to tracks
      //const std::vector<reco::TransientTrack>& tracks() const { return tracks_; };

      friend std::ostream& operator<< (std::ostream &out, const SVmassRecoSingleLegLikelihood& fit) 
      { 
	fit.printTo(out); 
	return out; 
      };

      void printTo(std::ostream &out) const
      {
         using namespace std;
         out << setw(10) << "Type: " << legType() << endl;
         out << setw(10) << "NLL" << setw(10) << nllOfLeg() << endl;
         out << setw(10) << "- NLLTopo" << setw(10) << nllTopological() << endl;
         out << setw(10) << "- NLLRapidity" << setw(10) << nllRapidity() << setw(10) << "y:" 
            << setw(10) << visRapidity() << setw(10) << "p:" << setw(10) << visP4().P() << endl;
         out << setw(10) << "- NLLDecay" << setw(10) << nllDecayLength() << setw(10) << "r:" 
            << setw(10) << legDir_.r() << setw(10) << "p:" << setw(10) << p4_.P() << endl;
         out << setw(10) << "- NLLPenalty" << setw(10) << nllM12Penalty() << endl;
         out << setw(10) << "-- SV" << setw(30) << sv_ << endl;
         out << setw(10) << "-- Dir" << setw(10) << legDir_ << endl;
         out << setw(10) << "-- VisP4" << setw(30) << visP4_ << " Mass: " << visP4_.mass() << endl;
         out << setw(10) << "-- NuP4" << setw(30) << nuP4_ << endl;
         out << setw(10) << "-- M12Up" << setw(10) << m12SquaredUpperBound(this->visP4(), this->dir()) << endl;
      }

   protected:
      /// Total visible p4
      FourVector fitVisP4() const { return visNeutralP4() + visChargedP4(); };

      /// Fit neutrino momentum 
      FourVector fitNuP4(double m12scale, int& error) const 
      { 
         double m12SquaredUpperBoundVal = m12SquaredUpperBound(this->visP4(), this->dir());
         double m12Squared = m12scale*m12scale*m12SquaredUpperBoundVal;
         FourVectorPair solutions = compInvisibleLeg(this->dir(), this->visP4(), tauMass, m12Squared, error);
         // Determine which solution to take
         if(ansatzForward_) return solutions.first;
         return solutions.second;
      }


   private:
      // Leg object
      const SVmassRecoSingleLegExtractorBase& extractor_;

      /// The associated tracks
      const std::vector<reco::TransientTrack>& tracks_;
      /// The trajectory states closes to the secondary vertex
      std::vector<TrajectoryStateClosestToPoint> tscps_;
      size_t nTracks_;

      /* Fit parameters */
      // Which neutrino solution to use
      bool ansatzForward_;
      // The secondary vertex
      GlobalPoint sv_;
      // Inferred direction of tau lepton
      ThreeVector legDir_;
      // Total fitted p4
      FourVector p4_;
      // Visible p4
      FourVector visP4_;
      // Nu p4
      FourVector nuP4_;     
};

}

#endif
