#ifndef TauAnalysis_CandidateTools_SVmassRecoSingleLegExtractorT_h
#define TauAnalysis_CandidateTools_SVmassRecoSingleLegExtractorT_h

/** \class SVmassRecoSingleLegExtractorT
 *
 * Auxiliary class for extracting variables used by secondary vertex based mass reconstruction algorithm
 * from pat::Electron, pat::Muon and pat::Tau objects
 * 
 * \author Christian Veelken, UC Davis
 *
 * \version $Revision: 1.1 $
 *
 * $Id: SVmassRecoSingleLegExtractorT.h,v 1.1 2009/06/11 07:23:28 veelken Exp $
 *
 */

#include "DataFormats/TauReco/interface/PFTauDecayMode.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

template<typename T>
class SVmassRecoSingleLegExtractorT : public SVmassRecoSingleLegExtractorBase
{
 public:
  // constructor 
  explicit SVmassRecoSingleLegExtractorT() 
    : leg_(0)
  {}
  
  // destructor
  virtual ~SVmassRecoSingleLegExtractorT() {}

  void setLeg(const T& leg) { leg_ = &leg; }

  // define "default" return values
  // for generic reco::Candidate case
  const reco::Candidate::LorentzVector& p4() const { return leg_->p4(); }
  bool typeIsSupportedBySVFitter() const { return false; }
  reco::Candidate::LorentzVector getNeutralP4() const { return reco::Candidate::LorentzVector(0,0,0,0); }
  double chargedMass2() const { return 0.; }
  bool nuSystemIsMassless() const { return true; }
  std::vector<reco::TrackBaseRef> getTracks() const { return std::vector<reco::TrackBaseRef>(); }

  int legTypeLabel() const { return -30; }

 private:
  const T* leg_;
};

//-------------------------------------------------------------------------------
// pat::Electron specific customization 
//-------------------------------------------------------------------------------

template<>
bool SVmassRecoSingleLegExtractorT<pat::Electron>::typeIsSupportedBySVFitter() const 
{ 
  return true; 
}

template<>
double SVmassRecoSingleLegExtractorT<pat::Electron>::chargedMass2() const
{
  const double electronMass = 5.109989e-4;
  return electronMass*electronMass;
}

template<>
bool SVmassRecoSingleLegExtractorT<pat::Electron>::nuSystemIsMassless() const
{
  return false;
}

template<>
std::vector<reco::TrackBaseRef> SVmassRecoSingleLegExtractorT<pat::Electron>::getTracks() const
{
  std::vector<reco::TrackBaseRef> tracks;
  tracks.push_back(reco::TrackBaseRef(leg_->gsfTrack()));
  return tracks;
}

template<>
int SVmassRecoSingleLegExtractorT<pat::Electron>::legTypeLabel() const
{
  return -2; 
}

//-------------------------------------------------------------------------------
// pat::Muon specific customization
//-------------------------------------------------------------------------------

template<>
bool SVmassRecoSingleLegExtractorT<pat::Muon>::typeIsSupportedBySVFitter() const 
{ 
  return true; 
}

template<>
double SVmassRecoSingleLegExtractorT<pat::Muon>::chargedMass2() const
{
  const double muonMass = 0.105658;
  return muonMass*muonMass;
}

template<>
bool SVmassRecoSingleLegExtractorT<pat::Muon>::nuSystemIsMassless() const
{
  return false;
}

template<>
std::vector<reco::TrackBaseRef> SVmassRecoSingleLegExtractorT<pat::Muon>::getTracks() const
{
  std::vector<reco::TrackBaseRef> tracks;
  tracks.push_back(reco::TrackBaseRef(leg_->track()));
  return tracks;
}

template<>
int SVmassRecoSingleLegExtractorT<pat::Muon>::legTypeLabel() const
{
  return -1; 
}

//-------------------------------------------------------------------------------
// pat::Tau specific customization
//-------------------------------------------------------------------------------

template<>
bool SVmassRecoSingleLegExtractorT<pat::Tau>::typeIsSupportedBySVFitter() const 
{ 
  return true; 
}

template<>
reco::Candidate::LorentzVector SVmassRecoSingleLegExtractorT<pat::Tau>::getNeutralP4() const 
{ 
  reco::Candidate::LorentzVector p4(0,0,0,0);
  const reco::PFCandidateRefVector& signalGammas = leg_->signalPFGammaCands();
  for ( reco::PFCandidateRefVector::const_iterator signalGamma = signalGammas.begin();
	signalGamma != signalGammas.end(); ++signalGamma ) {
    p4 += (*signalGamma)->p4();
  }
  return p4;
}

template<>
double SVmassRecoSingleLegExtractorT<pat::Tau>::chargedMass2() const
{
  const double chargedPionMass = 0.13957;
  return chargedPionMass*chargedPionMass;
}

template<>
bool SVmassRecoSingleLegExtractorT<pat::Tau>::nuSystemIsMassless() const
{
  return true;
}

template<>
std::vector<reco::TrackBaseRef> SVmassRecoSingleLegExtractorT<pat::Tau>::getTracks() const
{
  std::vector<reco::TrackBaseRef> tracks;
  const reco::PFCandidateRefVector& signalChargedHadrons = leg_->signalPFChargedHadrCands();
  unsigned numChargedHadrons = signalChargedHadrons.size();
  for ( unsigned iChargedHadron = 0; iChargedHadron < numChargedHadrons; ++iChargedHadron ) {
    tracks.push_back(reco::TrackBaseRef(signalChargedHadrons.at(iChargedHadron)->trackRef()));
  }
  return tracks;
}

template<>
int SVmassRecoSingleLegExtractorT<pat::Tau>::legTypeLabel() const
{
  return leg_->decayMode(); 
}

#endif  

