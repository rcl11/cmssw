#ifndef TauAnalysis_CandidateTools_SVmassRecoFitter_h
#define TauAnalysis_CandidateTools_SVmassRecoFitter_h

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

using namespace reco;
using namespace std;

namespace svMassReco {

  /// Function for Minuit to minimize
  void objectiveFcn(Int_t& npars, Double_t* grad, Double_t& fcn, Double_t* pars, Int_t flag)
  {
    // Get the fit object
    SVmassRecoDiTauLikelihoodBase* fit = dynamic_cast<SVmassRecoDiTauLikelihoodBase*>(gMinuit->GetObjectFit());
    if ( !fit ) {
      edm::LogError("objectiveFcn") << "Call of gMinuit::GetObjectFit returned NULL pointer !!";
      return;
    }

    Int_t error = 0;
    fcn = fit->nll(pars, error);
    //std::cout << *fit << std::endl;
    if ( error != 0 ) {
      //std::cout << "UNPHYSICAL FIT" << std::endl;
    }
  }

  struct Solution 
  {
    boost::shared_ptr<SVmassRecoDiTauLikelihoodBase> fitter;
    int solnType;
    double nllOfFit;
    double metNLL;
    int migradResult;
    int svFitResult;
    reco::Candidate::Point pv;
    reco::Candidate::Point sv1;
    reco::Candidate::Point sv2;
    double leg1MScale;
    double leg2MScale;
    int leg1Type;
    int leg2Type;
    FourVector leg1VisP4;
    FourVector leg1NuP4;
    FourVector leg2VisP4;
    FourVector leg2NuP4;
    // Sort operator
    bool operator<(const Solution& rhs) const { 
      // Prefer physical solutions first
      //if(svFitResult != rhs.svFitResult)
      //   return (svFitResult < rhs.svFitResult);
      //else // take the one with a better NLL
      return (nllOfFit < rhs.nllOfFit);
    }
    friend std::ostream& operator<< (std::ostream& o, const Solution& soln);
  };

  /// Pretty print a solution
  std::ostream& operator<< (std::ostream& o, const Solution& soln);

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
  bool isMatched(const TrackBaseRef&, const TrackBaseRef&); 

  /// Function to determine whether the given type is valid for the SV fitter
  /// Default case - specific implemntations in SVMethodT1T2Algorithm.cc
  template<typename T> bool typeIsSupportedBySVFitter() 
  { 
    return false; 
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
      minuit_.SetErrorDef(0.5);
      minuit_.SetPrintLevel(-1);
    }
    ~SVmassRecoFitter() {}

    /// Class to fit a ditau candidate.  T1 and T2 must be either pat::Muon, pat::Electron, or pat::Tau type
    std::vector<Solution> fitVertices(const edm::Ptr<T1>& leg1Ptr, const edm::Ptr<T2>& leg2Ptr, const CandidatePtr metCandPtr, 
				      const Vertex& pv, const BeamSpot& bs, const TransientTrackBuilder* trackBuilder)
    {
      if ( !(typeIsSupportedBySVFitter<T1>() && typeIsSupportedBySVFitter<T2>()) ) return std::vector<Solution>();

      const T1& leg1 = *leg1Ptr;
      const T2& leg2 = *leg2Ptr;

      edm::Ptr<reco::MET> metPtr = edm::Ptr<reco::MET>(metCandPtr);
      // TODO make const ref.
      const reco::MET* met = metPtr.get();
      
      // Get the tracks associated on each leg
      vector<TrackBaseRef> leg1Tracks = getTracks<T1>(leg1);
      vector<TrackBaseRef> leg2Tracks = getTracks<T2>(leg2);
      
      // Containers for the transient tracks on each leg
      vector<TransientTrack> leg1TransTracks;
      vector<TransientTrack> leg2TransTracks;
      
      TransientVertex cleanPV;
      // Get the clean PV and trans tracks.
      cleanupPV(pv, bs, trackBuilder, leg1Tracks, leg2Tracks, leg1TransTracks, leg2TransTracks, cleanPV);
      
      // There are four fits, for all comibination of forward/backward ansatzs.
      // We fit all four, then take the best as the solution
      std::vector<Solution> fits;
      
      // Do each fit
      for ( size_t soln = 0; soln < 4; soln++ ) {
	// Our fitter of 11 parameters
	//TMinuit minuit(11);
	gMinuit = &minuit_;
	minuit_.SetFCN(objectiveFcn);
	// Neg. log likelihood
	//minuit_.SetErrorDef(0.5);
	//minuit_.SetPrintLevel(-1);

	Solution mySoln;
	mySoln.solnType = soln;

	// Build the fitter
	mySoln.fitter = boost::shared_ptr<SVmassRecoDiTauLikelihoodBase>(
          new SVmassRecoDiTauLikelihood<T1,T2>(leg1, leg1TransTracks, leg2, leg2TransTracks, cleanPV, met, soln));

	// Set as active fit
	minuit_.SetObjectFit(mySoln.fitter.get());
	// Prepare Minuit for this fit, including choosing initial parameters
	mySoln.fitter->setupParametersFromScratch(minuit_, true);
	mySoln.fitter->enableMET(false);

	edm::LogInfo("SVMethod") << "Locking Nu Mass Scales";
	mySoln.fitter->lockNuMassScalers(minuit_);

	// Minimize!
	edm::LogInfo("SVMethod") << "Beginning minimization";
	minuit_.Migrad();

	// Hopefully we are near the minimum now. Release the mass scalers
	edm::LogInfo("SVMethod") << "Release Nu Mass Scales";
	mySoln.fitter->releaseNuMassScalers(minuit_);

	// Minimize again
	//mySoln.migradResult = minuit_.Migrad();
	
	// Now finally fit the PV simultaneously
	//mySoln.fitter->releasePV(minuit_);
	
	// Minimize
	//mySoln.migradResult = minuit_.Migrad();
	
	// Now allow MET in the fit
	edm::LogInfo("SVMethod") << "Reminimize";
	mySoln.migradResult = minuit_.Migrad();

	//edm::LogInfo("SVMethod") << "Removing nu scale limits...";
	// Now remove all limits
	//mySoln.fitter->setupParametersFromPrevious(minuit_, minuit_, false, true);
	//edm::LogInfo("SVMethod") << "Reminimizing";
	//mySoln.migradResult = minuit_.Migrad();
	
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
	
	mySoln.pv = reco::Candidate::Point(pv_x, pv_y, pv_z);
	mySoln.sv1 = reco::Candidate::Point(sv1_x, sv1_y, sv1_z);
	mySoln.sv2 = reco::Candidate::Point(sv2_x, sv2_y, sv2_z);
	mySoln.leg1MScale = m1scale;
	mySoln.leg2MScale = m2scale;
	
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
	
	mySoln.nllOfFit = mySoln.fitter->nll(pars, error);
	mySoln.metNLL = mySoln.fitter->metNLL();
	mySoln.svFitResult = error; // will be 1,2 or 3 if non-physical
	mySoln.leg1VisP4 = mySoln.fitter->leg1()->visP4();
	mySoln.leg2VisP4 = mySoln.fitter->leg2()->visP4();
	mySoln.leg1NuP4 = mySoln.fitter->leg1()->nuP4();
	mySoln.leg2NuP4 = mySoln.fitter->leg2()->nuP4();
	mySoln.leg1Type = mySoln.fitter->leg1()->legType();
	mySoln.leg2Type = mySoln.fitter->leg2()->legType();
	fits.push_back(mySoln);
	//edm::LogInfo("SVMethodFitResults") << mySoln;
      }

      // Sort the solution by how good the fit is
      sort(fits.begin(), fits.end());

      // Print
      for ( size_t iFit = 0; iFit < 4; iFit++ ) edm::LogInfo("fitResults") << fits[iFit];

      return fits;
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
	    if ( isMatched(pvTrack->first, *track) ) {
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
	    if ( isMatched(pvTrack->first, *track) ) {
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
    TMinuit minuit_;
  };

}

#endif
