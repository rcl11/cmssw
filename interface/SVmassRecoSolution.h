#ifndef TauAnalysis_CandidateTools_SVmassSolution_h
#define TauAnalysis_CandidateTools_SVmassSolution_h

/*
 * Container class to hold information about a solution from the SV mass
 * reconstruction algorithm
 * 
 * Author: Evan K. Friis (UC Davis)
 *
 */

#include "TauAnalysis/CandidateTools/interface/svMassRecoLikelihoodAuxFunctions.h"

namespace svMassReco {

   template<typename T1, typename T2>
   struct Solution 
   {
      boost::shared_ptr<SVmassRecoDiTauLikelihood<T1,T2> > fitter;
      int solnType;
      double nllOfFit;
      double metNLL;
      int migradResult;
      int svFitResult;
      reco::Candidate::Point pv;
      reco::Candidate::Point sv1;
      reco::Candidate::Point sv2;
      double leg1MScale;
      double leg2MScale;
      int leg1Type;
      int leg2Type;
      FourVector leg1VisP4;
      FourVector leg1NuP4;
      FourVector leg2VisP4;
      FourVector leg2NuP4;
      // Sort operator
      bool operator<(const Solution& rhs) const { 
         // Prefer physical solutions first
         //if(svFitResult != rhs.svFitResult)
         //   return (svFitResult < rhs.svFitResult);
         //else // take the one with a better NLL
         return (nllOfFit < rhs.nllOfFit);
      }

      friend std::ostream& operator<< (std::ostream& o, const Solution<T1,T2>& soln)
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

   };
}
#endif
