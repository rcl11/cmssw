#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "AnalysisDataFormats/TauAnalysis/interface/CompositePtrCandidateT1T2MEt.h"

#include "TauAnalysis/CandidateTools/interface/svMassRecoLikelihoodAuxFunctions.h"

#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TauAnalysis/CandidateTools/interface/SVmassRecoSingleLegExtractorT.h"

#include "DQMServices/Core/interface/DQMStore.h"
#include "DQMServices/Core/interface/MonitorElement.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include <TMath.h>

using namespace reco;
using namespace edm;
using namespace std;
using namespace svMassReco;

class SVmassLikelihoodValidator : public EDAnalyzer {
   public:
      explicit SVmassLikelihoodValidator(const ParameterSet& pset);
      ~SVmassLikelihoodValidator(){};
      virtual void analyze(const Event&, const EventSetup&);
   private:
      InputTag diTauCandidateSrc_;
      MonitorElement* hNLLTrueVertexLeg1Track_;
      MonitorElement* hNLLTrueMET_;
      MonitorElement* hProbChi2Vertex_;
      MonitorElement* hProbChi2MEt_ ;
      DQMStore* dqmStore_;
};

SVmassLikelihoodValidator::SVmassLikelihoodValidator(const ParameterSet& pset)
{
   diTauCandidateSrc_ = pset.getParameter<InputTag>("diTauCandidateSrc");
   dqmStore_ = &(*edm::Service<DQMStore>());
   dqmStore_->setCurrentFolder("SVmassValidation");
   // book histograms
   hNLLTrueVertexLeg1Track_ = dqmStore_->book1D("NLLTrueVertexLeg1", 
         "NLL of true SV given reco track", 80, -20, 20);
   hProbChi2Vertex_ = dqmStore_->book1D("ProbChi2Vertex", 
         "Prob of Vertex Chi2 for NDF=2", 102, -0.01, 1.01);
   hProbChi2MEt_ = dqmStore_->book1D("ProbChi2MEt", 
         "Prob of MEt Chi2 for NDF=2", 102, -0.01, 1.01);
   hNLLTrueMET_ = dqmStore_->book1D("NLLTrueMET", 
         "NLL of the true MET given reco MET", 80, -20, 20);
}

void SVmassLikelihoodValidator::analyze(const Event& evt, const EventSetup& es)
{
   const TransientTrackBuilder* trackBuilder = NULL;
   edm::ESHandle<TransientTrackBuilder> myTransientTrackBuilder;
   es.get<TransientTrackRecord>().get("TransientTrackBuilder", myTransientTrackBuilder);
   trackBuilder = myTransientTrackBuilder.product();

   typedef pat::Muon T1;
   typedef pat::Tau T2;
   // Get the diTau candidates
   typedef CompositePtrCandidateT1T2MEt<T1,T2> DiTauCandidate; 
   typedef std::vector<DiTauCandidate> DiTauCandidateCollection;

   edm::Handle<DiTauCandidateCollection> diTauCandidates;
   evt.getByLabel(diTauCandidateSrc_, diTauCandidates);

   SVmassRecoSingleLegExtractorT<T1> leg1Extractor;
   SVmassRecoSingleLegExtractorT<T2> leg2Extractor;

   for(size_t iDiTau = 0; iDiTau < diTauCandidates->size(); ++iDiTau)
   {
      const DiTauCandidate& ditau = diTauCandidates->at(iDiTau);

      ///  Compute the NLL of the LEG 1 track w.r.t true secondary vertex
      // Get the vertex of Leg1
      const reco::Candidate::Point& leg1TrueSVPoint = ditau.decayVertexPosLeg1gen();
      GlobalPoint leg1TrueSV = GlobalPoint(leg1TrueSVPoint.x(), leg1TrueSVPoint.y(), leg1TrueSVPoint.z());
      //cout << "Muon SV: " << leg1TrueSV << endl;
      //cout << "True PV: " << ditau.primaryVertexPosGen() << endl;

      // get leg 1 trans track
      leg1Extractor.setLeg(*ditau.leg1());
      leg2Extractor.setLeg(*ditau.leg2());

      vector<TrackBaseRef> leg1Tracks = leg1Extractor.getTracks();
      vector<TrackBaseRef> leg2Tracks = leg2Extractor.getTracks();

      cout << "Leg 1 reco eta: " << ditau.leg1()->p4().eta() << " leg 2 gen eta: " << ditau.p4VisLeg1gen().eta() << endl;
      // Build transient tracks
      if(leg1Tracks.size())
      {
         TransientTrack leg1Track = trackBuilder->build(leg1Tracks[0].castTo<TrackRef>());

         // Get the TSCP at the 
         TrajectoryStateClosestToPoint tscp = leg1Track.trajectoryStateClosestToPoint(leg1TrueSV);

         // Get the NLL
         double nll = nllPointGivenTrack(tscp);
         hNLLTrueVertexLeg1Track_->Fill(nll);
         hProbChi2Vertex_->Fill(TMath::Prob(nll, 2));

         //cout << "Got NLL: " << nll << endl;
      } 

      /// Compute the NLL of the reco MET w.r.t. true MET
      double metNLL = nllNuSystemGivenMET(ditau.p4InvisGen(), ditau.leg1()->p4(), ditau.met().get());
      hNLLTrueMET_->Fill(metNLL);
      hProbChi2MEt_->Fill(TMath::Prob(metNLL, 2));
   }
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(SVmassLikelihoodValidator);

