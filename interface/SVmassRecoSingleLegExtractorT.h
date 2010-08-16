#ifndef TauAnalysis_CandidateTools_SVmassRecoSingleLegExtractorT_h
#define TauAnalysis_CandidateTools_SVmassRecoSingleLegExtractorT_h

/** \class SVmassRecoSingleLegExtractorT
 *
 * Auxiliary class for extracting variables used by secondary vertex based mass reconstruction algorithm
 * from pat::Electron, pat::Muon and pat::Tau objects
 * 
 * \author Christian Veelken, UC Davis
 *
 * \version $Revision: 1.3 $
 *
 * $Id: SVmassRecoSingleLegExtractorT.h,v 1.3 2010/08/10 13:43:25 friis Exp $
 *
 */

#include "TauAnalysis/CandidateTools/interface/SVmassRecoSingleLegExtractorBase.h"
#include <assert.h>

// Type specific template specializatiosn found in corresponding .cc file
/*
namespace svMassRecoSingleLegExtractorTImpl {
   template<typename T> bool typeIsSupportedBySVFitter() { 
      assert(false); 
      std::cout << "What the hell" << std::endl;
      return false; 
   }
   template<typename T> bool nuSystemIsMassless() { return false; }
   template<typename T> std::vector<reco::TrackBaseRef> getTracks(const T* leg) { 
      assert(false);
      return std::vector<reco::TrackBaseRef>(); 
   }
   template<typename T> int legTypeLabel(const T* leg) { return -30; }
}
*/

namespace svMassRecoSingleLegExtractorTImpl {
   template<typename T> bool typeIsSupportedBySVFitter();
   template<typename T> bool nuSystemIsMassless();
   template<typename T> std::vector<reco::TrackBaseRef> getTracks(const T* leg);
   template<typename T> int legTypeLabel(const T* leg);
}

template<typename T>
class SVmassRecoSingleLegExtractorT : public SVmassRecoSingleLegExtractorBase
{
 public:
  // constructor 
  explicit SVmassRecoSingleLegExtractorT():leg_(0) {}
  
  // destructor
  virtual ~SVmassRecoSingleLegExtractorT() {}

  void setLeg(const T& leg) { leg_ = &leg; }
  reco::Candidate::LorentzVector p4() const { return leg_->p4(); }
  int charge() const { return leg_->charge(); }

  // Type specific functions are factored out into helper functions above
  bool typeIsSupportedBySVFitter() const { 
     return svMassRecoSingleLegExtractorTImpl::typeIsSupportedBySVFitter<T>(); 
  }
  bool nuSystemIsMassless() const { 
     return svMassRecoSingleLegExtractorTImpl::nuSystemIsMassless<T>(); 
  }
  std::vector<reco::TrackBaseRef> getTracks() const { 
     return svMassRecoSingleLegExtractorTImpl::getTracks<T>(leg_); 
  }
  int legTypeLabel() const { 
     return svMassRecoSingleLegExtractorTImpl::legTypeLabel<T>(leg_); 
  }

 private:
  const T* leg_;
};

#endif  

