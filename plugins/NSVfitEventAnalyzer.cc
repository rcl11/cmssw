#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/EventSetup.h"

#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"

#include "TauAnalysis/CandidateTools/plugins/NSVfitEventAnalyzer.h"
#include "TauAnalysis/CandidateTools/interface/PFMEtSignInterface.h"
#include "TauAnalysis/CandidateTools/interface/NSVfitStandaloneAlgorithm.h"


/// default constructor
NSVfitEventAnalyzer::NSVfitEventAnalyzer(const edm::ParameterSet& cfg): 
  met_  (cfg.getParameter<edm::InputTag>("met"  )),
  leps1_(cfg.getParameter<edm::InputTag>("muons")),
  leps2_(cfg.getParameter<edm::InputTag>("elecs"))
{
  // initialize MET significance calculation 
  metSign_ = new PFMEtSignInterface(cfg.getParameter<edm::ParameterSet>("metSignificance"));
}

/// default destructor
NSVfitEventAnalyzer::~NSVfitEventAnalyzer()
{
  delete metSign_;
}

/// everything that needs to be done during the event loop
void 
NSVfitEventAnalyzer::analyze(const edm::Event& event, const edm::EventSetup& eventSetup)
{
  // fetch MET
  edm::Handle<std::vector<pat::MET> > met;
  event.getByLabel(met_, met);
  // fetch selected muon
  edm::Handle<std::vector<pat::Muon> > leps1;
  event.getByLabel(leps1_, leps1);
  // fetch selected electron
  edm::Handle<std::vector<pat::Electron> > leps2;
  event.getByLabel(leps2_, leps2);
  // ...
  metSign_->beginEvent(event, eventSetup);
  // list of reco::Candidates to be subtracted from the list
  // of particle flow candidates, when calculating the MET
  // significance 
  std::list<const reco::Candidate*> candList;

  // make sure to prevent ambiguities
  if(leps1->size()>1 || leps2->size()>1){
    std::cout << "sorry got mixed up, too many e and mu candidates" << std::endl;
    return;
  }

  // loop leptons
  for(std::vector<pat::Muon>::const_iterator lep1=leps1->begin(); lep1!=leps1->end(); ++lep1){
    for(std::vector<pat::Electron>::const_iterator lep2=leps2->begin(); lep2!=leps2->end(); ++lep2){
      // determine met significance
      candList.push_back(dynamic_cast<const reco::Candidate*>(&*lep1));
      candList.push_back(dynamic_cast<const reco::Candidate*>(&*lep2));

      TMatrixD covMET(2,2);
      covMET = (*metSign_)(candList);

      // setup measure tau lepton vectors 
      std::vector<NSVfitStandalone::MeasuredTauLepton> measuredTauLeptons;
      measuredTauLeptons.push_back(NSVfitStandalone::MeasuredTauLepton(NSVfitStandalone::kLepDecay, lep1->p4()));
      measuredTauLeptons.push_back(NSVfitStandalone::MeasuredTauLepton(NSVfitStandalone::kLepDecay, lep2->p4()));
      
      // construct the class object from the minimal necesarry information
      NSVfitStandaloneAlgorithm algo(measuredTauLeptons, met->front().momentum(), covMET, 0);
      // apply customized configurations if wanted (examples are given below)
      algo.maxObjFunctionCalls(5000);
      //algo.addLogM(false);
      //algo.metPower(0.5)
      // run the fit
      algo.fit();
      // retrieve the results upon success
      //if(algo.isValidSolution()){
      std::cout << "--> getting mass = " << algo.fittedDiTauSystem().mass() << " + " << algo.massUncert() << " - " << algo.massUncert() << " [" << algo.fitStatus() << "]" << std::endl;
      //}
    }
  }
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(NSVfitEventAnalyzer);
