#ifndef TauAnalysis_RecoTools_SVMethodT1T2Algorithm_h
#define TauAnalysis_RecoTools_SVMethodT1T2Algorithm_h

#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"

#include "DataFormats/METReco/interface/MET.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"

#include "TauAnalysis/CandidateTools/interface/SVDiTauFitter.h"

using namespace reco;
using namespace std;

namespace TauVertex {

/// Function for Minuit to minimize
void objectiveFcn(Int_t& npars, Double_t* grad, Double_t& fcn, Double_t* pars, Int_t flag)
{
   // Get the fit object
   DiTauFitterInterface* fit = dynamic_cast<DiTauFitterInterface*>(gMinuit->GetObjectFit());
   Int_t error = 0;
   fcn = fit->nll(pars, error);
   //std::cout << *fit << std::endl;
   if(error != 0)
   {
      //std::cout << "UNPHYSICAL FIT" << std::endl;
   }
}

struct Solution {
   boost::shared_ptr<DiTauFitterInterface> fitter;
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
   TauVertex::FourVector leg1VisP4;
   TauVertex::FourVector leg1NuP4;
   TauVertex::FourVector leg2VisP4;
   TauVertex::FourVector leg2NuP4;
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

/// Function to remove all tracks belonging to the legs from a PV and refit it
void cleanupPV(const reco::Vertex& pv, const reco::BeamSpot& bs, const TransientTrackBuilder* trackBuilder,
      const vector<TrackBaseRef>& leg1Tracks, const vector<TrackBaseRef>& leg2Tracks,
      vector<TransientTrack>& leg1TransTracks, vector<TransientTrack>& leg2TransTracks,
      TransientVertex& cleanPV);

// Function to determine whether the given type is valid for the SV fitter
// Default case - specific implemntations in SVMethodT1T2Algorithm.cc
template<typename T> bool typeIsSupportedBySVFitter() { return false; }

// Foward declaration
template<typename T1, typename T2> std::vector<Solution>
fitVertices(const edm::Ptr<T1>& leg1Ptr, const edm::Ptr<T2>& leg2Ptr, const CandidatePtr metCandPtr, 
            const Vertex& pv, const BeamSpot& bs, const TransientTrackBuilder* trackBuilder);

// No candidates allowed!
template<> std::vector<Solution> 
fitVertices<reco::Candidate,reco::Candidate>(const edm::Ptr<reco::Candidate>& leg1Ptr, const edm::Ptr<reco::Candidate>& leg2Ptr, const CandidatePtr metCandPtr, 
            const Vertex& pv, const BeamSpot& bs, const TransientTrackBuilder* trackBuilder) { return std::vector<Solution>(); }

/// Class to fit a ditau candidate.  T1 and T2 must be in pat::Muon, pat::Electron, or pat::Tau
template<typename T1, typename T2> std::vector<Solution>
fitVertices(const edm::Ptr<T1>& leg1Ptr, const edm::Ptr<T2>& leg2Ptr, const CandidatePtr metCandPtr, 
            const Vertex& pv, const BeamSpot& bs, const TransientTrackBuilder* trackBuilder)
{
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
   cleanupPV(pv, bs, trackBuilder, leg1Tracks, leg2Tracks, leg1TransTracks, leg2TransTracks, cleanPV)

   // There are four fits, for all comibination of forward/backward ansatzs.
   // We fit all four, then take the best as the solution
   std::vector<Solution> fits;

   // Do each fit
   for(size_t soln = 0; soln < 4; soln++) {
      // Our fitter of 11 parameters
      TMinuit minuit(11);
      minuit.SetFCN(TauVertex::objectiveFcn);
      // Neg. log likelihood
      minuit.SetErrorDef(0.5);
      minuit.SetPrintLevel(-1);

      Solution mySoln;
      mySoln.solnType = soln;

      // Build the fitter
      mySoln.fitter = 
         boost::shared_ptr<DiTauFitterInterface>(new DiTauFitter<T1,T2>(leg1, leg1TransTracks, leg2, leg2TransTracks, 
									cleanPV, met, soln));

      // Set as active fit
      minuit.SetObjectFit(mySoln.fitter.get());
      // Prepare Minuit for this fit, including choosing initial parameters
      mySoln.fitter->setupParametersFromScratch(minuit, true);
      mySoln.fitter->enableMET(false);

      edm::LogInfo("SVMethod") << "Locking Nu Mass Scales";
      mySoln.fitter->lockNuMassScalers(minuit);

      // Minimize!
      edm::LogInfo("SVMethod") << "Beginning minimization";
      minuit.Migrad();

      // Hopefully we are near the minimum now. Release the mass scalers
      edm::LogInfo("SVMethod") << "Release Nu Mass Scales";
      mySoln.fitter->releaseNuMassScalers(minuit);

      // Minimize again
      //mySoln.migradResult = minuit.Migrad();

      // Now finally fit the PV simultaneously
      //mySoln.fitter->releasePV(minuit);

      // Minimize
      //mySoln.migradResult = minuit.Migrad();

      // Now allow MET in the fit
      edm::LogInfo("SVMethod") << "Reminimize";
      mySoln.migradResult = minuit.Migrad();

      //edm::LogInfo("SVMethod") << "Removing nu scale limits...";
      // Now remove all limits
      //mySoln.fitter->setupParametersFromPrevious(minuit, minuit, false, true);
      //edm::LogInfo("SVMethod") << "Reminimizing";
      //mySoln.migradResult = minuit.Migrad();

      // Get results
      Double_t pv_x, pv_y, pv_z, dummyError;
      minuit.GetParameter(0, pv_x, dummyError);
      minuit.GetParameter(1, pv_y, dummyError);
      minuit.GetParameter(2, pv_z, dummyError);

      Double_t sv1_x, sv1_y, sv1_z, m1scale;
      minuit.GetParameter(3, sv1_x, dummyError);
      minuit.GetParameter(4, sv1_y, dummyError);
      minuit.GetParameter(5, sv1_z, dummyError);
      minuit.GetParameter(6, m1scale, dummyError);

      Double_t sv2_x, sv2_y, sv2_z, m2scale;
      minuit.GetParameter(7, sv2_x, dummyError);
      minuit.GetParameter(8, sv2_y, dummyError);
      minuit.GetParameter(9, sv2_z, dummyError);
      minuit.GetParameter(10, m2scale, dummyError);

      mySoln.pv = reco::Candidate::Point(pv_x, pv_y, pv_z);
      mySoln.sv1 = reco::Candidate::Point(sv1_x, sv1_y, sv1_z);
      mySoln.sv2 = reco::Candidate::Point(sv2_x, sv2_y, sv2_z);
      mySoln.leg1MScale = m1scale;
      mySoln.leg2MScale = m2scale;

      /// Make sure this is a legal solution
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
   for(size_t iFit = 0; iFit < 4; iFit++)
      edm::LogInfo("fitResults") << fits[iFit];

   return fits;
}

}
#endif
