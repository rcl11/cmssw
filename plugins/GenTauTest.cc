#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "PhysicsTools/JetMCUtils/interface/JetMCTag.h"
#include "DataFormats/JetReco/interface/GenJet.h"

#include "TauAnalysis/CandidateTools/interface/SVMethodLikelihoods.h"

using namespace reco;
using namespace edm;
using namespace std;

class GenTauTestProducer : public EDProducer {
   public:
      explicit GenTauTestProducer(const ParameterSet& pset);
      ~GenTauTestProducer(){};
      virtual void produce(Event&, const EventSetup&);
   private:
      InputTag genTauJetSrc_;
};

GenTauTestProducer::GenTauTestProducer(const ParameterSet& pset)
{
   genTauJetSrc_ = pset.getParameter<InputTag>("genTauJetSrc");

   produces<vector<int> >("decayMode").setBranchAlias("decayMode");

   produces<std::vector<float> >("tauP4#px").setBranchAlias("tauP4#px");
   produces<std::vector<float> >("tauP4#py").setBranchAlias("tauP4#py");
   produces<std::vector<float> >("tauP4#pz").setBranchAlias("tauP4#pz");
   produces<std::vector<float> >("tauP4#e").setBranchAlias("tauP4#e");

   produces<std::vector<float> >("visP4#px").setBranchAlias("visP4#px");
   produces<std::vector<float> >("visP4#py").setBranchAlias("visP4#py");
   produces<std::vector<float> >("visP4#pz").setBranchAlias("visP4#pz");
   produces<std::vector<float> >("visP4#e").setBranchAlias("visP4#e");

   produces<std::vector<float> >("nuP4#px").setBranchAlias("nuP4#px");
   produces<std::vector<float> >("nuP4#py").setBranchAlias("nuP4#py");
   produces<std::vector<float> >("nuP4#pz").setBranchAlias("nuP4#pz");
   produces<std::vector<float> >("nuP4#e").setBranchAlias("nuP4#e");

   produces<std::vector<float> >("nuCompP41#px").setBranchAlias("nuCompP41#px");
   produces<std::vector<float> >("nuCompP41#py").setBranchAlias("nuCompP41#py");
   produces<std::vector<float> >("nuCompP41#pz").setBranchAlias("nuCompP41#pz");
   produces<std::vector<float> >("nuCompP41#e").setBranchAlias("nuCompP41#e");

   produces<std::vector<float> >("nuCompP42#px").setBranchAlias("nuCompP42#px");
   produces<std::vector<float> >("nuCompP42#py").setBranchAlias("nuCompP42#py");
   produces<std::vector<float> >("nuCompP42#pz").setBranchAlias("nuCompP42#pz");
   produces<std::vector<float> >("nuCompP42#e").setBranchAlias("nuCompP42#e");

   produces<std::vector<float> >("visP").setBranchAlias("visP");
   produces<std::vector<float> >("visCosTheta").setBranchAlias("visCosTheta");
   produces<std::vector<float> >("visRelRapidity").setBranchAlias("visRelRapidity");

   produces<std::vector<float> >("visEta").setBranchAlias("visEta");
   produces<std::vector<float> >("visPt").setBranchAlias("visPt");
   produces<std::vector<float> >("visMass").setBranchAlias("visMass");
   produces<std::vector<float> >("nuMass").setBranchAlias("nuMass");

   produces<std::vector<float> >("visPerp").setBranchAlias("visPerp");
   produces<std::vector<float> >("visPar").setBranchAlias("visPar");

   produces<std::vector<float> >("nuPerp").setBranchAlias("nuPerp");
   produces<std::vector<float> >("nuPar").setBranchAlias("nuPar");
}

void GenTauTestProducer::produce(Event& evt, const EventSetup& setup)
{
   Handle<GenJetCollection> tauGenJets;
   evt.getByLabel(genTauJetSrc_, tauGenJets);

   // our products
   auto_ptr<vector<int> > decayMode(new vector<int>);

   // SO LAME
   auto_ptr<std::vector<float> > tauP4_px(new std::vector<float>);
   auto_ptr<std::vector<float> > tauP4_py(new std::vector<float>);
   auto_ptr<std::vector<float> > tauP4_pz(new std::vector<float>);
   auto_ptr<std::vector<float> > tauP4_e(new std::vector<float>);

   auto_ptr<std::vector<float> > visP4_px(new std::vector<float>);
   auto_ptr<std::vector<float> > visP4_py(new std::vector<float>);
   auto_ptr<std::vector<float> > visP4_pz(new std::vector<float>);
   auto_ptr<std::vector<float> > visP4_e(new std::vector<float>);

   auto_ptr<std::vector<float> > nuP4_px(new std::vector<float>);
   auto_ptr<std::vector<float> > nuP4_py(new std::vector<float>);
   auto_ptr<std::vector<float> > nuP4_pz(new std::vector<float>);
   auto_ptr<std::vector<float> > nuP4_e(new std::vector<float>);

   auto_ptr<std::vector<float> > nuCompP41_px(new std::vector<float>);
   auto_ptr<std::vector<float> > nuCompP41_py(new std::vector<float>);
   auto_ptr<std::vector<float> > nuCompP41_pz(new std::vector<float>);
   auto_ptr<std::vector<float> > nuCompP41_e(new std::vector<float>);

   auto_ptr<std::vector<float> > nuCompP42_px(new std::vector<float>);
   auto_ptr<std::vector<float> > nuCompP42_py(new std::vector<float>);
   auto_ptr<std::vector<float> > nuCompP42_pz(new std::vector<float>);
   auto_ptr<std::vector<float> > nuCompP42_e(new std::vector<float>);

   auto_ptr<std::vector<float> > visP(new std::vector<float>);
   auto_ptr<std::vector<float> > visCosTheta(new std::vector<float>);
   auto_ptr<std::vector<float> > visRelRapidity(new std::vector<float>);

   auto_ptr<std::vector<float> > visPt(new std::vector<float>);
   auto_ptr<std::vector<float> > visEta(new std::vector<float>);
   auto_ptr<std::vector<float> > visMass(new std::vector<float>);
   // nu refers to Sum of all neutrinos, so this can be nonzero
   auto_ptr<std::vector<float> > nuMass(new std::vector<float>);

   auto_ptr<std::vector<float> > visPerp(new std::vector<float>);
   auto_ptr<std::vector<float> > visPar(new std::vector<float>);

   auto_ptr<std::vector<float> > nuPerp(new std::vector<float>);
   auto_ptr<std::vector<float> > nuPar(new std::vector<float>);

   // loop over taus
   for(GenJetCollection::const_iterator tau = tauGenJets->begin();
         tau != tauGenJets->end(); ++tau)
   {
      // Split into vis/invisible parts
      math::XYZTLorentzVector theVisP4;
      math::XYZTLorentzVector theNuP4;

      const CompositePtrCandidate::daughters& daughters = tau->daughterPtrVector();
      for(CompositePtrCandidate::daughters::const_iterator daughter = daughters.begin(); 
            daughter != daughters.end(); ++daughter ) 
      {
         int pdg_id = abs((*daughter)->pdgId());
         if(pdg_id == 12 || pdg_id == 14 || pdg_id == 16)
         {
            theNuP4 += (*daughter)->p4();
         } else
         {
            theVisP4 += (*daughter)->p4();
         }
      }
      if (theVisP4.pt() < 5 || abs(theVisP4.Eta()) > 2.1)
         continue;
      // Get tau level quantities 
      
      // NO reflex dictionariy for std::string??
      string namedDM = JetMCTagUtils::genTauDecayMode(*tau);
      int dm = -20;
      if (namedDM == "electron")
         dm = -1;
      else if (namedDM == "muon")
         dm = -2;
      else if (namedDM == "oneProng0Pi0") 
         dm = 0;
      else if (namedDM == "oneProng1Pi0")
         dm = 1;
      else if (namedDM == "oneProng2Pi0")
         dm = 2;
      else if (namedDM == "threeProng0Pi0")
         dm = 10;
      else if (namedDM == "threeProng1Pi0")
         dm = 11;

      decayMode->push_back(dm);;
      tauP4_px->push_back(tau->p4().px());
      tauP4_py->push_back(tau->p4().py());
      tauP4_pz->push_back(tau->p4().pz());
      tauP4_e->push_back(tau->p4().e());

      // Compute some helper quantities
      visP->push_back(theVisP4.Vect().r());
      double cosTheta = theVisP4.Vect().Unit().Dot(tau->p4().Vect().Unit());
      visCosTheta->push_back(cosTheta);
      double visRapid = atanh(theVisP4.Vect().Dot(tau->p4().Vect().Unit())/theVisP4.energy());
      visRelRapidity->push_back(visRapid);

      math::XYZVector visPar3 = theVisP4.Vect().Dot(tau->p4().Vect().Unit())*tau->p4().Vect().Unit();
      math::XYZVector visPerp3 = theVisP4.Vect() - visPar3;

      math::XYZVector nuPar3 = theNuP4.Vect().Dot(tau->p4().Vect().Unit())*tau->p4().Vect().Unit();
      math::XYZVector nuPerp3 = theNuP4.Vect() - nuPar3;

      visPerp->push_back(visPerp3.r());
      visPar->push_back(visPar3.r());

      nuPerp->push_back(nuPerp3.r());
      nuPar->push_back(nuPar3.r());

      double estNuMass = sqrt(
            (1.78)*(1.78) + theVisP4.mass()*theVisP4.mass() - 
            2*(1.78)*sqrt(theVisP4.mass()*theVisP4.mass() + visPerp3.r()*visPerp3.r()));

      /*
      std::cout << "********************** DM: " << dm <<  std::endl;
      std::cout << " Tau energy: " << tau->p4().e() << " Tau eta: " << tau->p4().eta() << " Tau phi: " << tau->p4().phi() << " mass: " << tau->p4().mass() << std::endl ;
      std::cout << " Vis energy: " << theVisP4.e() << " Vis eta: " << theVisP4.eta() <<  " Vis phi: " << theVisP4.phi() << std::endl;
      std::cout << " Nu energy: " << theNuP4.e() << " Nu mass: " << theNuP4.mass() << "(" << estNuMass << ")" << " Nu eta: " << theNuP4.eta() <<  " Nu phi: " << theNuP4.phi() << std::endl;
      math::XYZTLorentzVector  testicles = theVisP4 + theNuP4;
      std::cout << "TEST: " << testicles.mass() << std::endl;
      */

      // not for hadronic modes
      if(dm >= 0 || dm == -20)
         estNuMass = 0;

      int error;
      std::pair<reco::Candidate::LorentzVector,reco::Candidate::LorentzVector> theCompNuP4s = 
         TauVertex::compInvisibleLeg(tau->p4().Vect().Unit(), theVisP4, 1.78, estNuMass, error);

      visP4_px->push_back(theVisP4.px());
      visP4_py->push_back(theVisP4.py());
      visP4_pz->push_back(theVisP4.pz());
      visP4_e->push_back(theVisP4.e());

      nuP4_px->push_back(theNuP4.px());
      nuP4_py->push_back(theNuP4.py());
      nuP4_pz->push_back(theNuP4.pz());
      nuP4_e->push_back(theNuP4.e());

      nuCompP41_px->push_back(theCompNuP4s.first.px());
      nuCompP41_py->push_back(theCompNuP4s.first.py());
      nuCompP41_pz->push_back(theCompNuP4s.first.pz());
      nuCompP41_e->push_back(theCompNuP4s.first.e());

      nuCompP42_px->push_back(theCompNuP4s.second.px());
      nuCompP42_py->push_back(theCompNuP4s.second.py());
      nuCompP42_pz->push_back(theCompNuP4s.second.pz());
      nuCompP42_e->push_back(theCompNuP4s.second.e());

      visMass->push_back(theVisP4.mass());
      visEta->push_back(theVisP4.eta());
      visPt->push_back(theVisP4.pt());
      nuMass->push_back(theNuP4.mass());

   }

   evt.put(decayMode, "decayMode");

   evt.put(tauP4_px, "tauP4#px");
   evt.put(tauP4_py, "tauP4#py");
   evt.put(tauP4_pz, "tauP4#pz");
   evt.put(tauP4_e, "tauP4#e");

   evt.put(visP4_px, "visP4#px");
   evt.put(visP4_py, "visP4#py");
   evt.put(visP4_pz, "visP4#pz");
   evt.put(visP4_e, "visP4#e");

   evt.put(nuP4_px, "nuP4#px");
   evt.put(nuP4_py, "nuP4#py");
   evt.put(nuP4_pz, "nuP4#pz");
   evt.put(nuP4_e, "nuP4#e");

   evt.put(nuCompP41_px, "nuCompP41#px");
   evt.put(nuCompP41_py, "nuCompP41#py");
   evt.put(nuCompP41_pz, "nuCompP41#pz");
   evt.put(nuCompP41_e, "nuCompP41#e");

   evt.put(nuCompP42_px, "nuCompP42#px");
   evt.put(nuCompP42_py, "nuCompP42#py");
   evt.put(nuCompP42_pz, "nuCompP42#pz");
   evt.put(nuCompP42_e, "nuCompP42#e");

   evt.put(visP, "visP");
   evt.put(visCosTheta, "visCosTheta");
   evt.put(visRelRapidity, "visRelRapidity");
   evt.put(visMass, "visMass");
   evt.put(visPt, "visPt");
   evt.put(visEta, "visEta");
   evt.put(nuMass, "nuMass");
   evt.put(visPerp, "visPerp");
   evt.put(visPar, "visPar");
   evt.put(nuPerp, "nuPerp");
   evt.put(nuPar, "nuPar");
}
   
DEFINE_ANOTHER_FWK_MODULE(GenTauTestProducer);

