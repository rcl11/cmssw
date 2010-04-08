#ifndef TauAnalysis_CandidateTools_CollinearApproxFitter_h
#define TauAnalysis_CandidateTools_CollinearApproxFitter_h

/** \class CollinearApproxFitter
 *
 * Fit x and y components of missing transverse momentum
 * so as to best match visible tau+ tau- decay products
 * to a specific invariant mass (specified via configuration parameters).
 * The idea is to quantify the compatibility of observed tau decay products
 * with M --> tau+ tau- mass hypothesis of decay of resonance of mass M into tau-pairs,
 * reconstructed via collinear approximation.
 * 
 * \authors Christian Veelken
 *
 * \version $Revision: 1.2 $
 *
 * $Id: CollinearApproxFitter.h,v 1.2 2009/10/25 12:38:15 veelken Exp $
 *
 */

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/METReco/interface/MET.h"

#include "AnalysisDataFormats/TauAnalysis/interface/CollinearApproxCompatibility.h"

#include <TObject.h>
#include <TMinuit.h>
#include <TF1.h>

class CollinearApproxFitter : public TObject
{
  typedef edm::Ptr<reco::MET> MEtPtr;

 public:
  CollinearApproxFitter(const edm::ParameterSet&);
  ~CollinearApproxFitter();
     
//--- actual "public" fitting function
  CollinearApproxCompatibility fit(const reco::Candidate::LorentzVector&, const reco::Candidate::LorentzVector&, edm::Ptr<reco::MET>);

//--- auxiliary function needed to pass fit parameters to Minuit
  double compChi2(Double_t* parameter) const;

 protected:  

//--- configuration parameter
  double resonanceMass_;
  double resonanceWidth_;
  
//--- pointer to Minuit fitting algorithm
  TMinuit minuit_;

//--- pointers to MEt resolution functions
  TF1* metResolutionPx_;
  TF1* metResolutionPy_;

//--- temporary copies of parameters 
//    passed to Minuit during fit
  reco::Candidate::LorentzVector p4leg1_;
  reco::Candidate::LorentzVector p4leg2_;
  MEtPtr met_;
  double metPx_;
  double metPxErr_;
  double metPy_;
  double metPyErr_;
};

#endif
