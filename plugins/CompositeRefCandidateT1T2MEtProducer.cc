#include "TauAnalysis/CandidateTools/plugins/CompositeRefCandidateT1T2MEtProducer.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "PhysicsTools/Utilities/interface/deltaR.h"

#include "TauAnalysis/CandidateTools/interface/FetchCollection.h"

#include "AnalysisDataFormats/TauAnalysis/interface/CompositeRefCandidateT1T2MEt.h"
#include "AnalysisDataFormats/TauAnalysis/interface/CompositeRefCandidateT1T2MEtFwd.h"

template<typename T1, typename T2>
CompositeRefCandidateT1T2MEtProducer<T1,T2>::CompositeRefCandidateT1T2MEtProducer(const edm::ParameterSet& cfg)
  : algorithm_(cfg), cfgError_(0)
{
  useLeadingTausOnly_ = cfg.getParameter<bool>("useLeadingTausOnly");
  srcLeg1_ = cfg.getParameter<edm::InputTag>("srcLeg1");
  srcLeg2_ = cfg.getParameter<edm::InputTag>("srcLeg2");
  dRmin12_ = cfg.getParameter<double>("dRmin12");
  srcMET_ = ( cfg.exists("srcMET") ) ? cfg.getParameter<edm::InputTag>("srcMET") : edm::InputTag();
  recoMode_ = cfg.getParameter<std::string>("recoMode");
  verbosity_ = cfg.getUntrackedParameter<int>("verbosity", 0);

//--- check that InputTag for MET collection has been defined,
//    in case it is needed for the reconstruction mode 
//    specified in the configuration parameter set
  if ( srcMET_.label() == "" && recoMode_ != "" ) {
    edm::LogError ("CompositeRefCandidateT1T2MEtProducer") << " Configuration Parameter srcMET undefined," 
							   << " needed for recoMode = " << recoMode_ << " !!";
    cfgError_ = 1;
  }

  produces<CompositeRefCandidateCollection>("");
}
  
template<typename T1, typename T2>
CompositeRefCandidateT1T2MEtProducer<T1,T2>::~CompositeRefCandidateT1T2MEtProducer()
{
// nothing to be done yet...
}

template<typename T1, typename T2>
void CompositeRefCandidateT1T2MEtProducer<T1,T2>::produce(edm::Event& evt, const edm::EventSetup& es)
{
//--- print-out an error message and add an empty collection to the event 
//    in case of erroneous configuration parameters
  if ( cfgError_ ) {
    edm::LogError ("produce") << " Error in Configuration ParameterSet" 
			      << " --> CompositeRefCandidateT1T2MEt collection will NOT be produced !!";
    std::auto_ptr<CompositeRefCandidateCollection> emptyCompositeRefCandidateCollection(new CompositeRefCandidateCollection());
    evt.put(emptyCompositeRefCandidateCollection);
    return;
  }

  edm::Handle<std::vector<T1> > leg1Collection;
  pf::fetchCollection(leg1Collection, srcLeg1_, evt);
  edm::Handle<std::vector<T2> > leg2Collection;
  pf::fetchCollection(leg2Collection, srcLeg2_, evt);
  
  reco::CandidatePtr metPtr;
  if ( srcMET_.label() != "" ) {
    edm::Handle<reco::CandidateView> metCollection;
    pf::fetchCollection(metCollection, srcMET_, evt);

//--- check that there is exactly one MET object in the event
//    (missing transverse momentum is an **event level** quantity)
    if ( metCollection->size() == 1 ) {
      metPtr = metCollection->ptrAt(0);
    } else {
      edm::LogError ("produce") << " Found " << metCollection->size() << " MET objects in collection = " << srcMET_ << ","
				<< " --> CompositeRefCandidateT1T2MEt collection will NOT be produced !!";
      std::auto_ptr<CompositeRefCandidateCollection> emptyCompositeRefCandidateCollection(new CompositeRefCandidateCollection());
      evt.put(emptyCompositeRefCandidateCollection);
      return;
    }
  } 

//--- check if only one combination of tau decay products 
//    (the combination of highest Pt object in leg1 collection + highest Pt object in leg2 collection)
//    shall be produced, or all possible combinations of leg1 and leg2 objects
  std::auto_ptr<CompositeRefCandidateCollection> compositeRefCandidateCollection(new CompositeRefCandidateCollection());
  if ( useLeadingTausOnly_ ) {

//--- find highest Pt particles in leg1 and leg2 collections
    int idxLeadingLeg1 = -1;
    double leg1PtMax = 0.;
    for ( unsigned idxLeg1 = 0, numLeg1 = leg1Collection->size(); 
	  idxLeg1 < numLeg1; ++idxLeg1 ) {
      const T1& leg1 = leg1Collection->at(idxLeg1);
      if ( idxLeadingLeg1 == -1 || leg1.pt() > leg1PtMax ) {
	idxLeadingLeg1 = idxLeg1;
	leg1PtMax = leg1.pt();
      }
    }

    int idxLeadingLeg2 = -1;
    double leg2PtMax = 0.;
    for ( unsigned idxLeg2 = 0, numLeg2 = leg2Collection->size(); 
	  idxLeg2 < numLeg2; ++idxLeg2 ) {
      const T2& leg2 = leg2Collection->at(idxLeg2);

//--- do not create CompositeRefCandidateT1T2MEt object 
//    for combination of particle with itself
      if ( idxLeadingLeg1 != -1 ) {
	const T1& leadingLeg1 = leg1Collection->at(idxLeadingLeg1);
	double dR = reco::deltaR(leadingLeg1.p4(), leg2.p4());
	if ( dR < dRmin12_ ) continue;
      }

      if ( idxLeadingLeg2 == -1 || leg2.pt() > leg2PtMax ) {
	idxLeadingLeg2 = idxLeg2;
	leg2PtMax = leg2.pt();
      }
    }

    if ( idxLeadingLeg1 != -1 &&
	 idxLeadingLeg2 != -1 ) {
      T1Ref leadingLeg1Ref(leg1Collection, idxLeadingLeg1);
      T2Ref leadingLeg2Ref(leg2Collection, idxLeadingLeg2);

      CompositeRefCandidateT1T2MEt<T1,T2> compositeRefCandidate = 
	algorithm_.buildCompositeRefCandidate(leadingLeg1Ref, leadingLeg2Ref, metPtr);
      compositeRefCandidateCollection->push_back(compositeRefCandidate);
    } else {
      if ( verbosity_ >= 1 ) {
	edm::LogInfo ("produce") << " Found no combination of particles in Collections" 
				 << " leg1 = " << srcLeg1_ << " and leg2 = " << srcLeg2_ << ".";
      }
    }
  } else {
    for ( unsigned idxLeg1 = 0, numLeg1 = leg1Collection->size(); 
	  idxLeg1 < numLeg1; ++idxLeg1 ) {
      T1Ref leg1Ref(leg1Collection, idxLeg1);
      for ( unsigned idxLeg2 = 0, numLeg2 = leg2Collection->size(); 
	    idxLeg2 < numLeg2; ++idxLeg2 ) {
	T2Ref leg2Ref(leg2Collection, idxLeg2);

//--- do not create CompositeRefCandidateT1T2MEt object 
//    for combination of particle with itself
	double dR = reco::deltaR(leg1Ref->p4(), leg2Ref->p4());
	if ( dR < dRmin12_ ) continue;

	CompositeRefCandidateT1T2MEt<T1,T2> compositeRefCandidate = 
	  algorithm_.buildCompositeRefCandidate(leg1Ref, leg2Ref, metPtr);
	compositeRefCandidateCollection->push_back(compositeRefCandidate);
      }
    }
  }

//--- add the collection of reconstructed CompositeRefCandidateT1T2MEts to the event
  evt.put(compositeRefCandidateCollection);
}

#include "DataFormats/Candidate/interface/Candidate.h" 
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Tau.h"

typedef CompositeRefCandidateT1T2MEtProducer<reco::Candidate, reco::Candidate> DiCandidatePairProducer;
typedef CompositeRefCandidateT1T2MEtProducer<pat::Electron, pat::Tau> PATElecTauPairProducer;
typedef CompositeRefCandidateT1T2MEtProducer<pat::Muon, pat::Tau> PATMuTauPairProducer;
typedef CompositeRefCandidateT1T2MEtProducer<pat::Tau, pat::Tau> PATDiTauPairProducer;
typedef CompositeRefCandidateT1T2MEtProducer<pat::Electron, pat::Muon> PATElecMuPairProducer;

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_ANOTHER_FWK_MODULE(DiCandidatePairProducer);
DEFINE_ANOTHER_FWK_MODULE(PATElecTauPairProducer);
DEFINE_ANOTHER_FWK_MODULE(PATMuTauPairProducer);
DEFINE_ANOTHER_FWK_MODULE(PATDiTauPairProducer);
DEFINE_ANOTHER_FWK_MODULE(PATElecMuPairProducer);
