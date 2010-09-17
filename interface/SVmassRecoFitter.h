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
#include "FWCore/ParameterSet/interface/ParameterSet.h"

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

  const int verbosity = 0;

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
    //for(int i = 0; i < npars; i++) { cout << "  " << pars[i]; } cout << endl;
    fcn = fit->nll(pars, error);
    //std::cout << "FCN: " << fcn << endl;
    if ( error != 0 ) {
      //std::cout << "UNPHYSICAL FIT" << std::endl;
    }
  }

  template<typename T1, typename T2>
  class SVmassRecoFitter 
  {
   public:
    SVmassRecoFitter(const edm::ParameterSet& options)
      : options_(options),minuit_(11)
    {
      //std::cout << "<SVmassRecoFitter::SVmassRecoFitter>:" << std::endl;
      //std::cout << " disabling MINUIT output..." << std::endl;
      minuit_.SetPrintLevel(-1);
      minuit_.SetErrorDef(0.5);
    }
    ~SVmassRecoFitter() {}

    /// Class to fit a ditau candidate.  T1 and T2 must be either pat::Muon, pat::Electron, or pat::Tau type
    std::vector<SVmassRecoSolution> fitVertices(const edm::Ptr<T1>& leg1Ptr, 
          const edm::Ptr<T2>& leg2Ptr, const CandidatePtr metCandPtr, const Vertex& pv, const BeamSpot& bs, 
          const TransientTrackBuilder* trackBuilder)
    {

      if ( verbosity ) std::cout << "<SVmassRecoFitter::fitVertices>:" << std::endl;

       // Ensure both legs are supported by the fitter
      if ( !(leg1extractor_.typeIsSupportedBySVFitter() && leg2extractor_.typeIsSupportedBySVFitter()) ) 
      {
         edm::LogWarning("SVmassRecoFitter") << " One or both of the input collection types are unsupported," 
            << " SV fit will not be run!";
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
      
      // Get the clean PV and trans tracks.
      TransientVertex cleanPV;
      edm::LogInfo("SVRecoFitter") << "Cleaning primary vertex";
      cleanupPV(pv, bs, trackBuilder, leg1Tracks, leg2Tracks, leg1TransTracks, leg2TransTracks, cleanPV);
      
      // There are four fits, for all comibination of forward/backward ansatzs.
      // We fit all four, then take the best as the solution
      std::vector<SVmassRecoSolution> solutions;

      gMinuit = &minuit_;
      minuit_.SetFCN(objectiveFcn<T1,T2>);

      // build fitter
      edm::LogInfo("SVRecoFitter") << "Building di tau likelihood";

      SVmassRecoDiTauLikelihood fitter(options_, leg1extractor_, leg1TransTracks, leg2extractor_, leg2TransTracks, cleanPV, met);

      // Set as active fit
      minuit_.SetObjectFit(&fitter);
      fitter.setupParametersFromScratch(minuit_);

      // Set precision
      //Double_t precision[1];
      //precision[0] = 1e-6;
      //Int_t precErrFlag;
      //minuit_.mnexcm("SET EPS", precision, 1, precErrFlag);
      //edm::LogInfo("SVRecoFitter") << "Set precision to " << precision[0] << " got status" << precErrFlag;

      // First find the Hessian error matrix
      //Int_t hesseStatus;
      //minuit_.mncomd("HES", hesseStatus);
      //edm::LogInfo("SVRecoFitter") << "Built initial error matrix using HESSE.  Status: " << hesseStatus;

      int migradStatus = minuit_.Command("MIN");
      edm::LogInfo("SVRecoFitter") << "Final fit result: " << migradStatus << "\n" << fitter;

      // Retrive the information
      
      Double_t pv_x, pv_y, pv_z, dummyError;
      minuit_.GetParameter(0, pv_x, dummyError);
      minuit_.GetParameter(1, pv_y, dummyError);
      minuit_.GetParameter(2, pv_z, dummyError);

      Double_t sv1_theta, 
               sv1_phiLab,
               sv1_radiusLab,
               sv1_m12;

      minuit_.GetParameter(3, sv1_theta, dummyError);
      minuit_.GetParameter(4, sv1_phiLab, dummyError);
      minuit_.GetParameter(5, sv1_radiusLab, dummyError);
      minuit_.GetParameter(6, sv1_m12, dummyError);

      Double_t sv2_theta, 
               sv2_phiLab,
               sv2_radiusLab,
               sv2_m12;

      minuit_.GetParameter(7, sv2_theta, dummyError);
      minuit_.GetParameter(8, sv2_phiLab, dummyError);
      minuit_.GetParameter(9, sv2_radiusLab, dummyError);
      minuit_.GetParameter(10, sv2_m12, dummyError);

      // must be an easier way...
      Double_t pars[11];
      pars[0] = pv_x;
      pars[1] = pv_y;
      pars[2] = pv_z;
      pars[3] = sv1_theta;
      pars[4] = sv1_phiLab;
      pars[5] = sv1_radiusLab;
      pars[6] = sv1_m12;
      pars[7] = sv2_theta;
      pars[8] = sv2_phiLab;
      pars[9] = sv2_radiusLab;
      pars[10] = sv2_m12;
      if ( verbosity ) {
	for ( int iPar = 0; iPar < 11; ++iPar ) {
	  std::cout << " Parameter #" << iPar << " = " << pars[iPar] << std::endl;
	}
      }

      Int_t dummy_status;
      // Make sure our working point is fresh
      //fitter.nll(pars, dummy_status, true); // debug print-out enabled
      fitter.nll(pars, dummy_status, false); // debug print-out disabled

      SVmassRecoSolution solution(true, -1);
      solution.setP4VisLeg1(fitter.leg1Likelihood()->visP4());
      solution.setP4VisLeg2(fitter.leg2Likelihood()->visP4());
      solution.setP4Leg1(fitter.leg1Likelihood()->fittedP4());
      solution.setP4Leg2(fitter.leg2Likelihood()->fittedP4());
      solution.setP4(fitter.p4());
      solution.setSVrefittedPrimaryVertexPos(fitter.fittedPV());
      solution.setDecayVertexPosLeg1(fitter.leg1Likelihood()->sv());
      solution.setDecayVertexPosLeg2(fitter.leg2Likelihood()->sv());
      solution.setMscale1(fitter.leg1Likelihood()->m12());
      solution.setMscale2(fitter.leg2Likelihood()->m12());
      solution.setLogLikelihood(fitter.getNLL());
      //std::cout << " fitter.getNLL() = " << fitter.getNLL() << std::endl;
      solution.setLogLikelihoodMEt(fitter.metNLL());
      solution.setMigradStatus(migradStatus);

      // compute X1 and X2
      FourVector p4Leg1 = fitter.leg1Likelihood()->fittedP4();
      double x1 = ( p4Leg1.P() > 0. ) ? (fitter.leg1Likelihood()->visP4().P()/p4Leg1.P()) : -1.;
      solution.setX1(x1);
      FourVector p4Leg2 = fitter.leg2Likelihood()->fittedP4();
      double x2 = ( p4Leg2.P() > 0. ) ? (fitter.leg2Likelihood()->visP4().P()/p4Leg2.P()) : -1.;
      solution.setX2(x2);

      solution.setRestFrameVisThetaLeg1(fitter.leg1Likelihood()->thetaRest());
      solution.setRestFrameVisThetaLeg2(fitter.leg2Likelihood()->thetaRest());

      edm::LogInfo("SVRecoFitter") << "Final solution: " << solution;

      // ek dub
      for(size_t i = 0; i < 4; i++)
         solutions.push_back(solution);
      
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
    edm::ParameterSet options_;
    SVmassRecoSingleLegExtractorT<T1> leg1extractor_;
    SVmassRecoSingleLegExtractorT<T2> leg2extractor_;

    TMinuit minuit_;
  };

}

#endif
