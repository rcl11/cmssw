/*
 * SVLegFitter.cc
 *
 * Implements TauAnalysis/CandidateTools/interface/SVLegFitter.h
 * for Muon, Electron, and Tau type legs.
 *
 */

#include "TauAnalysis/CandidateTools/interface/SVLegFitter.h"
#include "TauAnalysis/CandidateTools/interface/SVMethodLikelihoods.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"

using namespace std;
using namespace reco;

namespace TauVertex {

/** Common leg fitter functions **/
void 
LegFitter::setPoints(const GlobalPoint& pv, double x, double y, double z, double m12scale, int& error)
{
   sv_ = GlobalPoint(x,y,z); 
   legDir_ = ThreeVector(x-pv.x(), y-pv.y(), z-pv.z());
   // Update the trajectory states closest to the SV
   for(size_t itrk = 0; itrk < nTracks_; ++itrk)
   {
      //std::cout << "Finding tcsp " << itrk << " @ " << sv_ << std::endl;
      tscps_[itrk] = tracks_[itrk].trajectoryStateClosestToPoint(sv_);
   }
   /// Update all the kinematic quantities
   visP4_ = fitVisP4();
   nuP4_ = fitNuP4(m12scale, error);
   p4_ = visP4_ + nuP4_;
}

double 
LegFitter::nllTopological() const
{
   /// TODO does this work in three prong case, instead of vertex fit?
   double output = 0.0;
   for(size_t itrk = 0; itrk < nTracks_; ++itrk)
      output += nllPointGivenTrack(tscps_[itrk]);
   return output;
}

double LegFitter::nllM12Penalty() const 
{
   double m12SquaredUpperBoundVal = m12SquaredUpperBound(this->visP4(), this->dir());
   //return -1.0*log(0.5 * (TMath::Erf(1000*m12SquaredUpperBoundVal)+1));
   if(m12SquaredUpperBoundVal < 0)
      return (m12SquaredUpperBoundVal*m12SquaredUpperBoundVal/1e-6);
   else 
      return 0;
}

double 
LegFitter::visRapidity() const
{ 
   return atanh(visP4_.Vect().Dot(legDir_.unit())/visP4_.e()); 
}

double 
LegFitter::nllDecayLength() const 
{
   double nll = nllTauDecayLengthGivenMomentum(legDir_.r(), p4_.P());
   return nll;
}

FourVector 
LegFitter::visChargedP4() const {
   TauVertex::FourVector output;
   for(size_t itrk = 0; itrk < nTracks_; ++itrk)
   {
      GlobalVector trkMomentumAtCP = tscps_[itrk].momentum();
      output += chargedP4FromMomentum(trkMomentumAtCP);
   }
   return output;
}

void
LegFitter::printTo(std::ostream &out) const {
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

template<> int legTypeLabel<pat::Electron>(const pat::Electron& leg) { return -2; }
template<> int legTypeLabel<pat::Muon>(const pat::Muon& leg) { return -1; }
template<> int legTypeLabel<pat::Tau>(const pat::Tau& leg) { return leg.decayMode(); }
template<> int legTypeLabel<reco::Candidate>(const reco::Candidate& leg) { return -3; }


/**** Functions to get tracks ****/

// Get tracks from tau
template<> vector<TrackBaseRef> getTracks<pat::Tau>(const pat::Tau& tau) {
   vector<TrackBaseRef> output;
   const PFCandidateRefVector& signalCHs = tau.signalPFChargedHadrCands();
   size_t nTracks = signalCHs.size();
   for(size_t itrk = 0; itrk < nTracks; itrk++)
   {
      output.push_back(TrackBaseRef(signalCHs.at(itrk)->trackRef()));
   }
   return output;
}

// Get track from muon
template<> vector<TrackBaseRef> getTracks<pat::Muon>(const pat::Muon& muon) {
   vector<TrackBaseRef> output;
   output.push_back(TrackBaseRef(muon.track())); //track() returns only inner track
   return output;
}

// Get track from electron
template<> vector<TrackBaseRef> getTracks<pat::Electron>(const pat::Electron& elec) {
   vector<TrackBaseRef> output;
   output.push_back(TrackBaseRef(elec.gsfTrack())); // Track or GSF track?
   return output;
}

// Get track from Candidate - not supported
template<> vector<TrackBaseRef> getTracks<reco::Candidate>(const reco::Candidate& elec) {
   vector<TrackBaseRef> output;
   return output;
}

/****    Mass template specializations ****/
// Hadronic case (pion)
template<> double chargedMass2ByType<pat::Tau>() { return chargedPionMass*chargedPionMass; }
// Muonic case
template<> double chargedMass2ByType<pat::Muon>() { return muonMass*muonMass; }
// Electronic case
template<> double chargedMass2ByType<pat::Electron>() { return electronMass*electronMass; }
// Candidate - not supported
template<> double chargedMass2ByType<reco::Candidate>() { return -1; }

// Get the upper limit on the neutrino system mass
double
m12SquaredUpperBound(const FourVector& visP4, const ThreeVector& tauDir) 
{
   double visPerpSquared = tauDir.unit().Cross(visP4.Vect()).mag2();
   double visMassSquared = visP4.mass()*visP4.mass();
   return (tauMass*tauMass + 
           visMassSquared - 2*tauMass*
           sqrt(visMassSquared + visPerpSquared)
          );
}


/****    Neutral component accessor templates ****/
// Candidate case - not supported
template<> FourVector getNeutralP4<reco::Candidate>(const reco::Candidate& object) { return FourVector(); }
/// The leptonic cases have no neutral component
template<> FourVector getNeutralP4<pat::Muon>(const pat::Muon& object) { return FourVector(); }
template<> FourVector getNeutralP4<pat::Electron>(const pat::Electron& object) { return FourVector(); }
template<> FourVector getNeutralP4<pat::Tau>(const pat::Tau& tau)
{
   FourVector output;
   const reco::PFCandidateRefVector& signalGammas = tau.signalPFGammaCands();
   // FIXME?
   //const reco::PFCandidateRefVector& signalNeutrHadrons = tau.signalPFGammaCands();
   for(size_t gamma = 0; gamma < signalGammas.size(); ++gamma)
   {
      output += signalGammas[gamma]->p4();
   }
   return output;
}

// M12 is always zero for hadronic decays, constrain it
template<> bool fixM12ScaleParameter<pat::Tau>() { return true;}
template<> bool fixM12ScaleParameter<pat::Electron>() { return false;}
template<> bool fixM12ScaleParameter<pat::Muon>() { return false;}
template<> bool fixM12ScaleParameter<reco::Candidate>() { return false;} // not supported

// Make available to linker
template class LegFitterSpecific<pat::Muon>;
template class LegFitterSpecific<pat::Electron>;
template class LegFitterSpecific<pat::Tau>;
template class LegFitterSpecific<reco::Candidate>;
}

