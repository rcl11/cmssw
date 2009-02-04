#ifndef TauAnalysis_CandidateTools_CompositeRefCandidateT1T2MEtProducer_h
#define TauAnalysis_CandidateTools_CompositeRefCandidateT1T2MEtProducer_h

/** \class CompositeRefCandidateT1T2MEtProducer
 *
 * Produce combinations of leptonic and hadronic decay products 
 * of a pair of tau leptons plus missing transverse momentum 
 * (representing the undetected momentum carried away by the neutrinos 
 *  produced in the two tau decays) 
 * 
 * \authors Colin Bernet,
 *          Michal Bluj,
 *          Christian Veelken
 *
 * \version $Revision: 1.1 $
 *
 * $Id: CompositeRefCandidateT1T2MEtProducer.h,v 1.1 2009/01/29 13:22:18 veelken Exp $
 *
 */

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/InputTag.h"

#include "TauAnalysis/CandidateTools/interface/CompositeRefCandidateT1T2MEtAlgorithm.h"

#include <string>

template<typename T1, typename T2>
class CompositeRefCandidateT1T2MEtProducer : public edm::EDProducer 
{
  typedef std::vector<T1> T1Collection;
  typedef edm::Ref<T1Collection> T1Ref;
  typedef std::vector<T2> T2Collection;
  typedef edm::Ref<T2Collection> T2Ref;

  typedef std::vector<CompositeRefCandidateT1T2MEt<T1,T2> > CompositeRefCandidateCollection;
  
 public:

  explicit CompositeRefCandidateT1T2MEtProducer(const edm::ParameterSet&);
  ~CompositeRefCandidateT1T2MEtProducer();

  void produce(edm::Event&, const edm::EventSetup&);

 private:

  CompositeRefCandidateT1T2MEtAlgorithm<T1,T2> algorithm_;
  
  bool useLeadingTausOnly_;
  edm::InputTag srcLeg1_;
  edm::InputTag srcLeg2_;
  double dRmin12_;
  edm::InputTag srcMET_;
  std::string recoMode_;
  int verbosity_;

  int cfgError_;
};

#endif

