#ifndef TauAnalysis_CandidateTools_SVmassRecoDiTauLikelihood_h
#define TauAnalysis_CandidateTools_SVmassRecoDiTauLikelihood_h

/*
* classes SVmassRecoDiTauLikelihood*
*
 * Authors: Evan K. Friis, Christian Veelken, UC Davis
 *
* A composite object to compute the likelihood of a diTau system.  Each diTau
* system has likelihood components that below to both 'legs', and then each
* leg has individual likelihood components.
*
* Common likelihoods:
*  o Likelihood of candidate primary vertex given measured primary vertex fit 
*  o Likelihood of MET given the missing energy on each leg
*
* Descriptions of the likelihood contributions from the legs are defined in
* TauAnalysis/CandidateTools/interface/SVmassRecoSingleLegLikelihood.h
*
* The negative log likelihood of the entire system is given by the function
*   double nll(Double_t* pars, Int_t &status, bool verbose=true);
*
* which computes the negative log likelihood of the entire system given
* the vector of paramters defined in Double_t *pars.
*
* The numbering of parameters is defined as follows
* -  0	   pv - x component relative to initial clean PV
* -  1	   pv - y component relative to initial clean PV
* -  2	   pv - z component relative to initial clean PV
* -  3	   sv - r component for leg 1 relative to initial PV
* -  4	   sv - theta component for leg 1 relative to fitted PV
* -  5	   sv - phi component for leg 1 relative to fitted PV
* -  6	   neutrino system mass scaler for leg 1 
* -  7	   sv - r component for leg 2 relative to fitted PV
* -  8	   sv - theta component for leg 2 relative to fitted PV
* -  9	   sv - z component for leg 2 relative to fitted PV
* -  10	   neutrino system mass scaler for leg 2 
*
* The neutrino system mass scaler only applies for leptonic tau decays. For
* hadronic decays, it is fixed at zero.
*/

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"

#include "TauAnalysis/CandidateTools/interface/SVmassRecoSingleLegLikelihood.h"
#include "TauAnalysis/CandidateTools/interface/svMassRecoLikelihoodAuxFunctions.h"
#include "TauAnalysis/CandidateTools/interface/SVmassRecoSingleLegExtractorBase.h"

#include "DataFormats/METReco/interface/MET.h"

#include <TObject.h>
#include <TMinuit.h>

namespace svMassReco 
{
  /// Pure abstract class to interface with Minuit
  class SVmassRecoDiTauLikelihood : public TObject
  {
   public:
     SVmassRecoDiTauLikelihood(const SVmassRecoSingleLegExtractorBase& leg1extractor, 
			       const std::vector<reco::TransientTrack>& leg1Tracks,
			       const SVmassRecoSingleLegExtractorBase& leg2extractor, 
			       const std::vector<reco::TransientTrack>& leg2Tracks, 
			       const TransientVertex& pv, const reco::MET* met)
       : pv_(pv),
         met_(met)
     {
        if(!pv_.isValid())
        {
           edm::LogError("SVDiTauLikelihood") << " PV is invalid!";
        }
        // build each leg fitter
        leg1Likelihood_ = 
           std::auto_ptr<SVmassRecoSingleLegLikelihood>(new SVmassRecoSingleLegLikelihood(leg1extractor, leg1Tracks));
        leg2Likelihood_ = 
           std::auto_ptr<SVmassRecoSingleLegLikelihood>(new SVmassRecoSingleLegLikelihood(leg2extractor, leg2Tracks));
        enableMET_ = true;
     }
     virtual ~SVmassRecoDiTauLikelihood(){};

     /// Get NLL given parameters
     double nll(Double_t* pars, Int_t& status, bool verbose=false)
     {

       double nllOutput = 0;
       /* npars layout:  (PV dx,dy,dz), (Leg1 x,y,z,m12) (Leg2 x,y,z,m12) */

       // Fitted PV is original PV position plus a fitted correction factor
       fitPV_ = pv_.position() + GlobalVector(pars[0], pars[1], pars[2]);

       // Get the NLL of the PV given the vertex fit
       pvNLL_ = nllPointGivenVertex(fitPV_, pv_);
       nllOutput += pvNLL_;
	 
       // Get the NLL for each leg
       leg1Likelihood_->setPoints(fitPV_, pars[3], pars[4], pars[5], pars[6]);

       nllOutput += leg1Likelihood_->nllOfLeg();

       leg2Likelihood_->setPoints(fitPV_, pars[7], pars[8], pars[9], pars[10]);

       nllOutput += leg2Likelihood_->nllOfLeg();

       // Constrain by MET, if desired.
       // Compute the four vector sum of invisible objects
       nus_ = leg2Likelihood_->nuP4() + leg1Likelihood_->nuP4();
       // MET is parameterized w.r.t muon direction
       metNLL_ = nllNuSystemGivenMET(nus_, leg1Likelihood_->visP4(), met_);
       if ( enableMET_ ) nllOutput += metNLL_;

       // Cache NLL result
       lastNLL_ = nllOutput;
       
       return nllOutput;
     }

     // P4 of both di tau system
     FourVector p4() const { return leg1Likelihood_->fittedP4() + leg2Likelihood_->fittedP4(); }

     // P4 of all fitted invisible objects 
     const FourVector& invisibleP4() const { return nus_; }

     const GlobalPoint& fittedPV() const { return fitPV_; };

     /// Get MET NLL for the last points fitted
     double metNLL() const { return metNLL_; };

     /// Get eotal NLL for the last points fitted
     double getNLL() const { return lastNLL_; };

     /// Setup the parameters in Minuit and find reasonable physical initial
     /// values
     void setupParametersFromScratch(TMinuit& minuit)
     {

       edm::LogInfo("SVMethod") << "Setting up initial PV";
       // ParNo, Name, InitValue, InitError, Limits
       minuit.DefineParameter(0, "pv_x", 0, pv_.positionError().cxx(), 0, 0);
       minuit.DefineParameter(1, "pv_y", 0, pv_.positionError().cyy(), 0, 0);
       minuit.DefineParameter(2, "pv_z", 0, pv_.positionError().czz(), 0, 0);
       
       edm::LogInfo("SVMethod") << "Finding initial conditions for leg 1";
       // Get initial guesses for each leg
       double thetaGuess, thetaError, 
              phiGuess, phiError,
              radiusGuess, radiusError,
              m12Guess, m12Error;

       // Find initial values for this leg
       leg1Likelihood_->findIntialValues(pv_.position(),
             thetaGuess, thetaError, 
             phiGuess, phiError,
             radiusGuess, radiusError,
             m12Guess, m12Error);

       minuit.DefineParameter(3, "sv1_thetaRest", thetaGuess, thetaError, 0, TMath::Pi());
       minuit.DefineParameter(4, "sv1_phiLab", phiGuess, phiError, 0, 0);
       minuit.DefineParameter(5, "sv1_radiusLab", radiusGuess, radiusError, 0, 0); //limits?
       minuit.DefineParameter(6, "sv1_m12", m12Guess, m12Error, 0, leg1Likelihood_->maxM12()); 

       edm::LogInfo("SVMethod") << "Finding initial SV for leg 2";
       leg2Likelihood_->findIntialValues(pv_.position(),
             thetaGuess, thetaError, 
             phiGuess, phiError,
             radiusGuess, radiusError,
             m12Guess, m12Error);

       minuit.DefineParameter(7, "sv2_thetaRest", thetaGuess, thetaError, 0, TMath::Pi());
       minuit.DefineParameter(8, "sv2_phiLab", phiGuess, phiError, 0, 0);
       minuit.DefineParameter(9, "sv2_radiusLab", radiusGuess, radiusError, 0, 0); //limits?
       minuit.DefineParameter(10, "sv2_m12", m12Guess, m12Error, 0, leg2Likelihood_->maxM12()); 

       // Determine whether or not m12 is a free parameter and needs to be scaled.
       // For hadronic taus, there is only 1 nu, so it's mass is always zero and we 
       // can ignore m12
       if ( leg1Likelihood_->nuSystemIsMassless() ) minuit.FixParameter(6);
       if ( leg2Likelihood_->nuSystemIsMassless() ) minuit.FixParameter(10);
     }

     /// Lock the PV variables in the fit
     void lockPV(TMinuit& minuit) {
       minuit.FixParameter(0);
       minuit.FixParameter(1);
       minuit.FixParameter(2);
     }

     /// Allow PV variables to float
     void releasePV(TMinuit& minuit) {
       minuit.Release(0);
       minuit.Release(1);
       minuit.Release(2);
     }

     /// Fix the Nu mass scale variables
     void lockNuMassScalers(TMinuit& minuit) {
       minuit.FixParameter(6);
       minuit.FixParameter(10);
     }

     /// Allow the Nu mass scale variables to float, if applicable
     void releaseNuMassScalers(TMinuit& minuit) {
       // Only let the ones that can actually be variable vary
       if ( !leg1Likelihood_->nuSystemIsMassless() ) minuit.Release(6);
       if ( !leg2Likelihood_->nuSystemIsMassless() ) minuit.Release(10);
     }

     /// Enable/disable MET in the fit
     void enableMET(bool enable) {
        enableMET_ = enable;
     }

    /// Access to the legs
    const SVmassRecoSingleLegLikelihood* leg1Likelihood() const { return leg1Likelihood_.get(); };
    const SVmassRecoSingleLegLikelihood* leg2Likelihood() const { return leg2Likelihood_.get(); };

    friend std::ostream& operator<< (std::ostream &out, const SVmassRecoDiTauLikelihood& fit) 
    { 
      fit.printTo(out); 
      return out; 
    };

    void printTo(std::ostream &out) const
    {
       using namespace std;
       FourVector total = leg1Likelihood_->fittedP4() + leg2Likelihood_->fittedP4();
       out << "=========== SVmassRecoDiTauLikelihood =========== " << endl;
       out << setw(10) <<"Total NLL:" << setw(10) << lastNLL_ << setw(10) << "Status:" << setw(10) << lastStatus_ << endl;
       out << setw(10) <<"PV" << setw(10) << pvNLL_ << setw(10) << "Fit:" << setw(10) << fitPV_ << "Meas:" << setw(10) << pv_.position() << endl;
       if ( enableMET_ )
          out << setw(10) <<"MET" << setw(10) << metNLL_ << setw(10) << "Fit:" << setw(30) << nus_ << setw(10) << "Meas:" << setw(30) << met_->p4() << endl;
       else
          out << "MET disabled." << endl;
       out << "Total mass: " << total.mass() << endl;
       out << " *** Leg 1 *** " << endl;
       leg1Likelihood_->printTo(out);
       out << " *** Leg 2 *** " << endl;
       leg2Likelihood_->printTo(out);
    }
     
   private:
     const TransientVertex& pv_;
     const reco::MET* met_;
     
     // Logging stuff
     GlobalPoint fitPV_;
     double pvNLL_;
     double metNLL_;
     FourVector nus_;
     
     double lastNLL_;
     int lastStatus_;
     std::auto_ptr<SVmassRecoSingleLegLikelihood> leg1Likelihood_;
     std::auto_ptr<SVmassRecoSingleLegLikelihood> leg2Likelihood_;
     
     bool enableMET_;
  };
}

#endif
