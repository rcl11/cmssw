#include "TauAnalysis/CandidateTools/interface/SVmassRecoFitter.h"

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Tau.h"

namespace svMassReco
{ 

  // Give supported types
  template<> bool typeIsSupportedBySVFitter<pat::Tau>() { return true; }
  template<> bool typeIsSupportedBySVFitter<pat::Muon>() { return true; }
  template<> bool typeIsSupportedBySVFitter<pat::Electron>() { return true; }

  // Dump a solution
  std::ostream& operator<< (std::ostream& o, const Solution& soln)
  {
    FourVector totalP4 = soln.leg1VisP4 + soln.leg2VisP4 + soln.leg1NuP4 + soln.leg2NuP4;
    return o 
      << " Fit Results:  Solution index (" << soln.solnType << ")" << std::endl
      << " * NLL = " << soln.nllOfFit  << " MET NLL " << soln.metNLL << std::endl
      << " Migrad(" << soln.migradResult << ") Physicality(" << soln.svFitResult << ")" << std::endl
      << " * Leg 1: Vis (Pt,Eta,Phi) = (" << soln.leg1VisP4.pt() << ", " << soln.leg1VisP4.eta() << ", " << soln.leg1VisP4.phi() << ") Type: " << soln.leg1Type << std::endl
      << " Nu (Pt,Eta,Phi) = (" << soln.leg1NuP4.pt() << ", " << soln.leg1NuP4.eta() << ", " << soln.leg1NuP4.phi() << ")" << std::endl
      << " * Leg 2: Vis (Pt,Eta,Phi) = (" << soln.leg2VisP4.pt() << ", " << soln.leg2VisP4.eta() << ", " << soln.leg2VisP4.phi() << ") Type: " << soln.leg2Type << std::endl
      << " Nu (Pt,Eta,Phi) = (" << soln.leg2NuP4.pt() << ", " << soln.leg2NuP4.eta() << ", " << soln.leg2NuP4.phi() << ")" << std::endl
      << " Total (Pt,Eta,Phi,M) = (" << totalP4.pt() << ", " << totalP4.eta() << ", " << totalP4.phi() << "," << totalP4.mass() << ")" << std::endl;
  }

  // Adapted from RecoBTag/SecondaryVertex/plugins/SecondaryVertexProducer.cc
  // key comparison for tracking map
  bool isMatched(const TrackBaseRef& trk1, const TrackBaseRef& trk2) 
  {
    double dEta = trk1->eta() - trk2->eta();
    double dPhi = trk1->phi() - trk2->phi();
    double dR2 = dEta*dEta + dPhi*dPhi;
    if ( dR2 < (0.01*0.01) && trk1->charge() == trk2->charge() ) {
      if ( abs(trk1->pt() - trk2->pt())/(trk2->pt()+trk1->pt()) < 0.05 ) return true;
    }
    return false;
  }

}
