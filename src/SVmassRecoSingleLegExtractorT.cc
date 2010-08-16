#include "TauAnalysis/CandidateTools/interface/SVmassRecoSingleLegExtractorT.h"

#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Tau.h"

namespace svMassRecoSingleLegExtractorTImpl {

   //-------------------------------------------------------------------------------
   // Unsupported Candidate customization
   //-------------------------------------------------------------------------------

   template<> bool typeIsSupportedBySVFitter<reco::Candidate>()  
   { 
      return false; 
   }

   template<> bool nuSystemIsMassless<reco::Candidate>() 
   {
      return false;
   }

   template<> std::vector<reco::TrackBaseRef> getTracks<reco::Candidate>(const reco::Candidate* leg) 
   {
      return std::vector<reco::TrackBaseRef>();
   }

   template<> int legTypeLabel<reco::Candidate>(const reco::Candidate* leg) 
   {
      return -30; 
   }

   //-------------------------------------------------------------------------------
   // pat::Electron specific customization 
   //-------------------------------------------------------------------------------

   template<> bool typeIsSupportedBySVFitter<pat::Electron>()  
   { 
      return true; 
   }

   template<> bool nuSystemIsMassless<pat::Electron>() 
   {
      return false;
   }

   template<> std::vector<reco::TrackBaseRef> getTracks<pat::Electron>(const pat::Electron* leg) 
   {
      std::vector<reco::TrackBaseRef> tracks;
      tracks.push_back(reco::TrackBaseRef(leg->gsfTrack()));
      return tracks;
   }

   template<> int legTypeLabel<pat::Electron>(const pat::Electron* leg) 
   {
      return -2; 
   }

   //-------------------------------------------------------------------------------
   // pat::Muon specific customization
   //-------------------------------------------------------------------------------

   template<> bool typeIsSupportedBySVFitter<pat::Muon>()  
   { 
      return true; 
   }

   template<> bool  nuSystemIsMassless<pat::Muon>() 
   {
      return false;
   }

   template<> std::vector<reco::TrackBaseRef> getTracks<pat::Muon>(const pat::Muon* leg) 
   {
      std::vector<reco::TrackBaseRef> tracks;
      tracks.push_back(reco::TrackBaseRef(leg->innerTrack()));
      return tracks;
   }

   template<> int legTypeLabel<pat::Muon>(const pat::Muon* leg) 
   {
      return -1; 
   }

   //-------------------------------------------------------------------------------
   // pat::Tau specific customization
   //-------------------------------------------------------------------------------

   template<> bool typeIsSupportedBySVFitter<pat::Tau>()  
   { 
      return true; 
   }

   template<> bool nuSystemIsMassless<pat::Tau>() 
   {
      return true;
   }

   template<> std::vector<reco::TrackBaseRef> getTracks<pat::Tau>(const pat::Tau* leg) 
   {
      std::vector<reco::TrackBaseRef> tracks;
      const reco::PFCandidateRefVector& signalChargedHadrons = leg->signalPFChargedHadrCands();
      unsigned numChargedHadrons = signalChargedHadrons.size();
      for ( unsigned iChargedHadron = 0; iChargedHadron < numChargedHadrons; ++iChargedHadron ) {
         tracks.push_back(reco::TrackBaseRef(signalChargedHadrons.at(iChargedHadron)->trackRef()));
      }
      return tracks;
   }

   template<> int legTypeLabel<pat::Tau>(const pat::Tau* leg) 
   {
      return leg->decayMode(); 
   }


   // Declare our objects to the linker
   template bool typeIsSupportedBySVFitter<pat::Muon>();
   template bool nuSystemIsMassless<pat::Muon>();
   template std::vector<reco::TrackBaseRef> getTracks<pat::Muon>(const pat::Muon*);
   template int legTypeLabel<pat::Muon>(const pat::Muon*);

   template bool typeIsSupportedBySVFitter<pat::Electron>();
   template bool nuSystemIsMassless<pat::Electron>();
   template std::vector<reco::TrackBaseRef> getTracks<pat::Electron>(const pat::Electron*);
   template int legTypeLabel<pat::Electron>(const pat::Electron*);

   template bool typeIsSupportedBySVFitter<pat::Tau>();
   template bool nuSystemIsMassless<pat::Tau>();
   template std::vector<reco::TrackBaseRef> getTracks<pat::Tau>(const pat::Tau*);
   template int legTypeLabel<pat::Tau>(const pat::Tau*);

} // end namespace
