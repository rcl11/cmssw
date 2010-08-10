#ifndef TauAnalysis_CandidateTools_SVmassRecoSingleLegExtractorBase_h
#define TauAnalysis_CandidateTools_SVmassRecoSingleLegExtractorBase_h

/** \class SVmassRecoSingleLegExtractorBase
 *
 * Base-class for extracting variables used by secondary vertex based mass reconstruction algorithm
 * from pat::Electron, pat::Muon and pat::Tau objects
 * (derived classes are templated and implement different behavior 
 *  for pat::Electron, pat::Muon and pat::Tau and generic reco::Candidate case)
 * 
 * \author Christian Veelken, UC Davis
 *
 * \version $Revision: 1.2 $
 *
 * $Id: SVmassRecoSingleLegExtractorBase.h,v 1.2 2010/06/06 21:55:20 friis Exp $
 *
 */

#include "FWCore/ParameterSet/interface/ParameterSet.h"

class SVmassRecoSingleLegExtractorBase
{
 public:
  // constructor 
  explicit SVmassRecoSingleLegExtractorBase() {}
  
  // destructor
  virtual ~SVmassRecoSingleLegExtractorBase() {}

  virtual reco::Candidate::LorentzVector p4() const = 0;
  virtual int charge() const = 0;
  virtual bool typeIsSupportedBySVFitter() const = 0;
  virtual reco::Candidate::LorentzVector getNeutralP4() const = 0;
  virtual double chargedMass2() const = 0;
  virtual bool nuSystemIsMassless() const = 0;
  virtual std::vector<reco::TrackBaseRef> getTracks() const = 0;

  virtual int legTypeLabel() const = 0;
};

#endif  

