#include "TauAnalysis/CandidateTools/plugins/PATElecTauPairZeeHypothesisProducer.h"

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/View.h"
#include "DataFormats/Common/interface/Ptr.h"

#include "AnalysisDataFormats/TauAnalysis/interface/CompositePtrCandidateT1T2MEt.h"
#include "AnalysisDataFormats/TauAnalysis/interface/CompositePtrCandidateT1T2MEtFwd.h"

#include "AnalysisDataFormats/TauAnalysis/interface/PATElecTauPairZeeHypothesis.h"
#include "AnalysisDataFormats/TauAnalysis/interface/PATElecTauPairZeeHypothesisFwd.h"

#include "TauAnalysis/CandidateTools/interface/FetchCollection.h"

PATElecTauPairZeeHypothesisProducer::PATElecTauPairZeeHypothesisProducer(const edm::ParameterSet& cfg)
  : algorithm_(cfg)
{
  srcElecTauPairs_ = cfg.getParameter<edm::InputTag>("elecTauPairSource");

  produces<PATElecTauPairZeeHypothesisCollection>("");
}

PATElecTauPairZeeHypothesisProducer::~PATElecTauPairZeeHypothesisProducer()
{
//--- nothing to be done yet...
}

void PATElecTauPairZeeHypothesisProducer::produce(edm::Event& evt, const edm::EventSetup& es)
{
  edm::Handle<edm::View<PATElecTauPair> > elecTauPairs;
  evt.getByLabel(srcElecTauPairs_, elecTauPairs);

  std::auto_ptr<PATElecTauPairZeeHypothesisCollection> elecTauPairZeeHypotheses(new PATElecTauPairZeeHypothesisCollection());

  for ( unsigned idxElecTauPair = 0, numElecTauPairs = elecTauPairs->size(); 
	idxElecTauPair < numElecTauPairs; ++idxElecTauPair ) {
    edm::Ptr<PATElecTauPair> elecTauPair = elecTauPairs->ptrAt(idxElecTauPair);

    PATElecTauPairZeeHypothesis elecTauPairZeeHypothesis = algorithm_.buildZeeHypothesis(elecTauPair, evt, es);
    elecTauPairZeeHypotheses->push_back(elecTauPairZeeHypothesis);
  }

  evt.put(elecTauPairZeeHypotheses);
}

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_ANOTHER_FWK_MODULE(PATElecTauPairZeeHypothesisProducer);



