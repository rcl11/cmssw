#include "TauAnalysis/CandidateTools/interface/SVMethodFitter.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Tau.h"

/// Predicate to key tracks
namespace {
   template<typename T>
      struct RefToBaseLess : public std::binary_function<edm::RefToBase<T>,
      edm::RefToBase<T>,
      bool> {
         inline bool operator()(const edm::RefToBase<T> &r1,
               const edm::RefToBase<T> &r2) const
         {
            return r1.id() < r2.id() || (r1.id() == r2.id() && r1.key() < r2.key());
         }
      };
}

namespace TauVertex
{ 
// Give supported types
template<> bool typeIsSupportedBySVFitter<pat::Tau>() { return true; }
template<> bool typeIsSupportedBySVFitter<pat::Muon>() { return true; }
template<> bool typeIsSupportedBySVFitter<pat::Electron>() { return true; }

// Dump a solution
std::ostream& operator<< (std::ostream& o, const Solution& soln)
{
   FourVector totalP4 = soln.leg1VisP4 + soln.leg2VisP4 + soln.leg1NuP4 + soln.leg2NuP4;
   return o 
      << " Fit Results:  Solution index (" << soln.solnType << ")" << std::endl
      << " * NLL = " << soln.nllOfFit  << " MET NLL " << soln.metNLL << std::endl
      << " Migrad(" << soln.migradResult << ") Physicality(" << soln.svFitResult << ")" << std::endl
      << " * Leg 1: Vis (Pt,Eta,Phi) = (" << soln.leg1VisP4.pt() << ", " << soln.leg1VisP4.eta() << ", " << soln.leg1VisP4.phi() << ") Type: " << soln.leg1Type << std::endl
      << " Nu (Pt,Eta,Phi) = (" << soln.leg1NuP4.pt() << ", " << soln.leg1NuP4.eta() << ", " << soln.leg1NuP4.phi() << ")" << std::endl
      << " * Leg 2: Vis (Pt,Eta,Phi) = (" << soln.leg2VisP4.pt() << ", " << soln.leg2VisP4.eta() << ", " << soln.leg2VisP4.phi() << ") Type: " << soln.leg2Type << std::endl
      << " Nu (Pt,Eta,Phi) = (" << soln.leg2NuP4.pt() << ", " << soln.leg2NuP4.eta() << ", " << soln.leg2NuP4.phi() << ")" << std::endl
      << " Total (Pt,Eta,Phi,M) = (" << totalP4.pt() << ", " << totalP4.eta() << ", " << totalP4.phi() << "," << totalP4.mass() << ")" << std::endl;
}

// Adapted from RecoBTag/SecondaryVertex/plugins/SecondaryVertexProducer.cc
// key comparison for tracking map


bool isMatched(const TrackBaseRef& trk1, const TrackBaseRef& trk2)
{
   double dEta = trk1->eta() - trk2->eta();
   double dPhi = trk1->phi() - trk2->phi();
   double dR2 = dEta*dEta + dPhi*dPhi;
   if(dR2 < (0.01*0.01) && trk1->charge() == trk2->charge())
   {
      if(abs(trk1->pt() - trk2->pt())/(trk2->pt()+trk1->pt()) < 0.05 )
         return true;
   }
   return false;
}


// Remove tracks associated w/ the legs from the PV and refit it
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
   for(Vertex::trackRef_iterator iter = pv.tracks_begin(); iter != pv.tracks_end(); ++iter) 
   {
     // Build the transient track
     TransientTrack track = trackBuilder->build(iter->castTo<TrackRef>());
     pvTracks.insert(make_pair(*iter, track));
     pvOriginalTracks.push_back(track);
   }

   edm::LogInfo("cleanupPV")  << "Before cleaning PV has " << pvOriginalTracks.size() << " tracks";

   // Remove the leg1 and leg2 tracks from the PV
   for(vector<TrackBaseRef>::const_iterator track = leg1Tracks.begin(); track != leg1Tracks.end(); ++track)
   {
      TransientTrackMap::iterator pos = pvTracks.find(*track);
      // If this track is in the PV, delete it!
      if(pos != pvTracks.end())
      {
         edm::LogInfo("cleanupPV")  << "Matched a leg1 track by ref";
         pvTracks.erase(pos);
         leg1TransTracks.push_back(pos->second);
      } else 
      {
         // Delete it if we can match it
         for(TransientTrackMap::iterator pvTrack = pvTracks.begin(); pvTrack != pvTracks.end(); ++pvTrack)
         {
            if(isMatched(pvTrack->first, *track))
            {
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

   for(vector<TrackBaseRef>::const_iterator track = leg2Tracks.begin(); track != leg2Tracks.end(); ++track)
   {
      TransientTrackMap::iterator pos = pvTracks.find(*track);
      if(pos != pvTracks.end())
      {
         edm::LogInfo("cleanupPV")  << "Matched a leg2 track by ref";
         pvTracks.erase(pos);
         leg2TransTracks.push_back(pos->second);
      } else
      {
         // Delete it if we can match it
         for(TransientTrackMap::iterator pvTrack = pvTracks.begin(); pvTrack != pvTracks.end(); ++pvTrack)
         {
            if(isMatched(pvTrack->first, *track))
            {
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
   for(TransientTrackMap::iterator pvTrack = pvTracks.begin(); pvTrack != pvTracks.end(); ++pvTrack)
   {
      pvTransTracks.push_back(pvTrack->second);
   }

   // Refit the cleaned PV
   KalmanVertexFitter vertexFitter(true);
   if(pvTransTracks.size() > 1)
   {
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

} //end namespace
