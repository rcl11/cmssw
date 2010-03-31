#ifndef TauAnalysis_CandidateTools_SVmassRecoFitter_h
#define TauAnalysis_CandidateTools_SVmassRecoFitter_h

/*
 * SVmassRecoFitter
 *
 * Authors: Evan K. Friis, Christian Veelken (UC Davis)
 *
 * Class used to fit the composite mass of diTau candidates.
 *
 * The method fits the three coordinates of the PV and the three coordinates of
 * both secondary vertices. 
 *
 * There are four possible solutions for any fit, corresponding to the 2
 * possible alignments of the neutrino direction in the rest rame of each leg.
 * All four solutions are returned, in the format
 * std::vector<svMassReco::Solution>.  The solutions are ordered by
 * increasing negative log likelihood, so the first solution shoudl generally
 * be the best.
 *
 * The computation of the likelihood for a given set of fit parameters is done
 * by the template class SVmassRecoDiTauLikelihood.  This class computes
 * likelihoods that are specific to the di-tau system, such as the
 * compatability with the measured PV and the measured MET.  The
 * SVmassRecoDiTauLikelihood additionally owns two
 * SVmassRecoSingleLegLikelihood instances.  Each of these instances computes
 * the likelihoods that are leg specific - the compatability of the associated
 * SV with the leg tracks, the likelihood of the visible decay angle, and the
 * likelihood of the PV-SV leg decay length given the total energy of the leg.
 * 
 */

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"

#include "DataFormats/METReco/interface/MET.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"

#include "TauAnalysis/CandidateTools/interface/SVmassRecoDiTauLikelihood.h"
#include "TauAnalysis/CandidateTools/interface/SVmassRecoSingleLegExtractorT.h"
#include "AnalysisDataFormats/TauAnalysis/interface/SVmassRecoSolution.h"

using namespace reco;
using namespace std;

namespace svMassReco {

  /// Function for Minuit to minimize
  template<typename T1, typename T2>
  void objectiveFcn(Int_t& npars, Double_t* grad, Double_t& fcn, Double_t* pars, Int_t flag)
  {
    // Get the fit object
    SVmassRecoDiTauLikelihood* fit = dynamic_cast<SVmassRecoDiTauLikelihood*>(gMinuit->GetObjectFit());
    if ( !fit ) {
      edm::LogError("objectiveFcn") << "Call of gMinuit::GetObjectFit returned NULL pointer !!";
      return;
    }

    Int_t error = 0;
    fcn = fit->nll(pars, error);
    if ( error != 0 ) {
      //std::cout << "UNPHYSICAL FIT" << std::endl;
    }
  }

  template<typename T1, typename T2>
  class SVmassRecoFitter 
  {
   public:
    SVmassRecoFitter()
      : minuit_(11)
    {
      std::cout << "<SVmassRecoFitter::SVmassRecoFitter>:" << std::endl;
      std::cout << " disabling MINUIT output..." << std::endl;
      minuit_.SetPrintLevel(-1);
      minuit_.SetErrorDef(0.5);
    }
    ~SVmassRecoFitter() {}

    /// Class to fit a ditau candidate.  T1 and T2 must be either pat::Muon, pat::Electron, or pat::Tau type
    std::vector<SVmassRecoSolution> fitVertices(const edm::Ptr<T1>& leg1Ptr, const edm::Ptr<T2>& leg2Ptr, const CandidatePtr metCandPtr, 
						const Vertex& pv, const BeamSpot& bs, const TransientTrackBuilder* trackBuilder)
    {
      if ( !(leg1extractor_.typeIsSupportedBySVFitter() && leg2extractor_.typeIsSupportedBySVFitter()) ) {
	return std::vector<SVmassRecoSolution>();
      }

      const T1& leg1 = *leg1Ptr;
      const T2& leg2 = *leg2Ptr;

      edm::Ptr<reco::MET> metPtr = edm::Ptr<reco::MET>(metCandPtr);
      // TODO make const ref.
      const reco::MET* met = metPtr.get();
      
      // initialize extractor objects with leg1, leg2 objects passed as function arguments
      leg1extractor_.setLeg(leg1);
      leg2extractor_.setLeg(leg2);

      // Get the tracks associated on each leg
      vector<TrackBaseRef> leg1Tracks = leg1extractor_.getTracks();
      vector<TrackBaseRef> leg2Tracks = leg2extractor_.getTracks();
      
      // Containers for the transient tracks on each leg
      vector<TransientTrack> leg1TransTracks;
      vector<TransientTrack> leg2TransTracks;
      
      TransientVertex cleanPV;
      // Get the clean PV and trans tracks.
      cleanupPV(pv, bs, trackBuilder, leg1Tracks, leg2Tracks, leg1TransTracks, leg2TransTracks, cleanPV);
      
      // There are four fits, for all comibination of forward/backward ansatzs.
      // We fit all four, then take the best as the solution
      std::vector<SVmassRecoSolution> solutions;
      
      // Do each fit
      for ( size_t iSolution = 0; iSolution < 4; ++iSolution ) {
	gMinuit = &minuit_;
	minuit_.SetFCN(objectiveFcn<T1,T2>);

        // Solutions not yet sorted by log-likelihood; set rank to -1
	SVmassRecoSolution solution(true, -1);

	// Build the fitter
	boost::shared_ptr<SVmassRecoDiTauLikelihood> fitter(
	  new SVmassRecoDiTauLikelihood(leg1extractor_, leg1TransTracks, leg2extractor_, leg2TransTracks, cleanPV, met, iSolution));

	// Set as active fit
	minuit_.SetObjectFit(fitter.get());
	// Prepare Minuit for this fit, including choosing initial parameters
	fitter->setupParametersFromScratch(minuit_, true);
	fitter->enableMET(true);

	edm::LogInfo("SVMethod") << "Locking Nu Mass Scales";
	fitter->lockNuMassScalers(minuit_);

	// Minimize!
	edm::LogInfo("SVMethod") << "Beginning minimization";
	minuit_.Migrad();
	//int migradStatus = minuit_.Migrad();

	// Hopefully we are near the minimum now. Release the mass scalers
	edm::LogInfo("SVMethod") << "Release Nu Mass Scales";
	fitter->releaseNuMassScalers(minuit_);
	
	// Now allow MET in the fit
	edm::LogInfo("SVMethod") << "Reminimize";
	int migradStatus = minuit_.Migrad();

	// Get results
	Double_t pv_x, pv_y, pv_z, dummyError;
	minuit_.GetParameter(0, pv_x, dummyError);
	minuit_.GetParameter(1, pv_y, dummyError);
	minuit_.GetParameter(2, pv_z, dummyError);
	
	Double_t sv1_x, sv1_y, sv1_z, m1scale;
	minuit_.GetParameter(3, sv1_x, dummyError);
	minuit_.GetParameter(4, sv1_y, dummyError);
	minuit_.GetParameter(5, sv1_z, dummyError);
	minuit_.GetParameter(6, m1scale, dummyError);
	
	Double_t sv2_x, sv2_y, sv2_z, m2scale;
	minuit_.GetParameter(7, sv2_x, dummyError);
	minuit_.GetParameter(8, sv2_y, dummyError);
	minuit_.GetParameter(9, sv2_z, dummyError);
	minuit_.GetParameter(10, m2scale, dummyError);

	solution.setSVrefittedPrimaryVertexPos(reco::Candidate::Point(pv_x, pv_y, pv_z));
	solution.setDecayVertexPosLeg1(reco::Candidate::Point(sv1_x, sv1_y, sv1_z));
	solution.setDecayVertexPosLeg2(reco::Candidate::Point(sv2_x, sv2_y, sv2_z));
        solution.setMscale1(m1scale);
	solution.setMscale2(m2scale);
	
	// Make sure this is a legal solution
	int error = 0;
	Double_t pars[11];
	pars[0] = pv_x;
	pars[1] = pv_y;
	pars[2] = pv_z;
	pars[3] = sv1_x;
	pars[4] = sv1_y;
	pars[5] = sv1_z;
	pars[6] = m1scale;
	pars[7] = sv2_x;
	pars[8] = sv2_y;
	pars[9] = sv2_z;
	pars[10] = m2scale;

	solution.setMigradStatus(migradStatus);
	solution.setLogLikelihood(fitter->nll(pars, error));
	solution.setLogLikelihoodMEt(fitter->metNLL());
	solution.setSVfitStatus(error); // will be 1,2 or 3 in case solution is non-physical, 4 in case any parameter is "not-a-number"
	solution.setP4VisLeg1(fitter->leg1Likelihood()->visP4());
	solution.setP4VisLeg2(fitter->leg2Likelihood()->visP4());
	reco::Candidate::LorentzVector p4Leg1 = fitter->leg1Likelihood()->visP4() + fitter->leg1Likelihood()->nuP4();
	double x1 = ( p4Leg1.P() > 0. ) ? (fitter->leg1Likelihood()->visP4().P()/p4Leg1.P()) : -1.;
	solution.setX1(x1);
	reco::Candidate::LorentzVector p4Leg2 = fitter->leg2Likelihood()->visP4() + fitter->leg2Likelihood()->nuP4();
	double x2 = ( p4Leg2.P() > 0. ) ? (fitter->leg2Likelihood()->visP4().P()/p4Leg2.P()) : -1.;
	solution.setX2(x2);
	solution.setP4(p4Leg1 + p4Leg2);

	solutions.push_back(solution);
      }

      // Sort the solution by how good the fit is
      sort(solutions.begin(), solutions.end());

      size_t numSolutions = solutions.size();
      for ( size_t iSolution = 0; iSolution < numSolutions; ++iSolution ) {
	solutions[iSolution].setRank(iSolution);
	edm::LogInfo("fitResults") << solutions[iSolution];
      }

      return solutions;
    }

  protected:    
    /// Function to remove all tracks belonging to the legs from a PV and refit it
    void cleanupPV(const reco::Vertex& pv, const reco::BeamSpot& bs, const TransientTrackBuilder* trackBuilder,
		   const vector<TrackBaseRef>& leg1Tracks, const vector<TrackBaseRef>& leg2Tracks, 
		   vector<TransientTrack>& leg1TransTracks, vector<TransientTrack>& leg2TransTracks, 
		   TransientVertex& cleanPV)
    {
      if ( !trackBuilder ) {
	edm::LogError("cleanupPV") << "Invalid pointer to trackBuilder object --> skipping !!";
	return;
      }
      
      typedef std::map<TrackBaseRef, TransientTrack, RefToBaseLess<Track> > TransientTrackMap;
      // Only used if there are not enough unique tracks to make a fit
      vector<TransientTrack> pvOriginalTracks; 
      
      // Get the tracks associated to the PV
      TransientTrackMap pvTracks;
      for ( Vertex::trackRef_iterator iter = pv.tracks_begin(); iter != pv.tracks_end(); ++iter ) {
	// Build the transient track
	TransientTrack track = trackBuilder->build(iter->castTo<TrackRef>());
	pvTracks.insert(make_pair(*iter, track));
	pvOriginalTracks.push_back(track);
      }
      
      edm::LogInfo("cleanupPV")  << "Before cleaning PV has " << pvOriginalTracks.size() << " tracks";
      
      // Remove the leg1 and leg2 tracks from the PV
      for ( vector<TrackBaseRef>::const_iterator track = leg1Tracks.begin(); track != leg1Tracks.end(); ++track ) {
	TransientTrackMap::iterator pos = pvTracks.find(*track);
	// If this track is in the PV, delete it!
	if ( pos != pvTracks.end() ) {
	  edm::LogInfo("cleanupPV") << "Matched a leg1 track by ref";
	  pvTracks.erase(pos);
	  leg1TransTracks.push_back(pos->second);
	} else {
	  // Delete it if we can match it
	  for ( TransientTrackMap::iterator pvTrack = pvTracks.begin(); pvTrack != pvTracks.end(); ++pvTrack ) {
	    if ( tracksAreMatched(pvTrack->first, *track) ) {
	      edm::LogInfo("cleanupPV")  << "Matched a leg1 track by comparison";
	      pvTracks.erase(pvTrack);
	      break;
	    }
	  }
	  TrackRef trackRef = track->castTo<TrackRef>();
	  TransientTrack newTrack = trackBuilder->build(trackRef);
	  leg1TransTracks.push_back(newTrack);
	}
      }
      
      for ( vector<TrackBaseRef>::const_iterator track = leg2Tracks.begin(); track != leg2Tracks.end(); ++track ) {
	TransientTrackMap::iterator pos = pvTracks.find(*track);
	if ( pos != pvTracks.end() ) {
	  edm::LogInfo("cleanupPV")  << "Matched a leg2 track by ref";
	  pvTracks.erase(pos);
	  leg2TransTracks.push_back(pos->second);
	} else {
	  // Delete it if we can match it
	  for ( TransientTrackMap::iterator pvTrack = pvTracks.begin(); pvTrack != pvTracks.end(); ++pvTrack ) {
	    if ( tracksAreMatched(pvTrack->first, *track) ) {
	      edm::LogInfo("cleanupPV")  << "Matched a leg2 track by comparison";
	      pvTracks.erase(pvTrack);
	      break;
	    }
	  }
	  TrackRef trackRef = track->castTo<TrackRef>();
	  TransientTrack newTrack = trackBuilder->build(trackRef);
	  leg2TransTracks.push_back(newTrack);
	}
      }
      
      vector<TransientTrack> pvTransTracks;
      for ( TransientTrackMap::iterator pvTrack = pvTracks.begin(); pvTrack != pvTracks.end(); ++pvTrack ) {
	pvTransTracks.push_back(pvTrack->second);
      }
      
      // Refit the cleaned PV
      KalmanVertexFitter vertexFitter(true);
      if ( pvTransTracks.size() > 1 ) {
	cleanPV = vertexFitter.vertex(pvTransTracks, bs);
      } else {
	// If there aren't enough tracks, just refit the PV 
	edm::LogWarning("cleanupPV") << "Not enough remaining tracks, PV not refitted!" << endl;
	cleanPV = vertexFitter.vertex(pvOriginalTracks, bs);
      }
      
      edm::LogInfo("cleanupPV")  
	<< "After cleaning:  PV (" << cleanPV.refittedTracks().size() << ")"
	<< " Leg1 (" << leg1TransTracks.size() << ") Leg2 (" << leg2TransTracks.size() << ")";
    }
    
  private:
    SVmassRecoSingleLegExtractorT<T1> leg1extractor_;
    SVmassRecoSingleLegExtractorT<T2> leg2extractor_;

    TMinuit minuit_;
  };

}

#endif
