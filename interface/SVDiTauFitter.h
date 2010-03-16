#ifndef TauAnalysis_CandidateTools_SVDiTauFitter_h
#define TauAnalysis_CandidateTools_SVDiTauFitter_h
/* 
 * TauVertex::DiTauFitter
 *
 * Author: Evan K. Friis, UC Davis
 *
 * Class to compute the NLL of a di-tau system
 *
 */

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "TauAnalysis/CandidateTools/interface/SVLegFitter.h"
#include <TObject.h>
#include <TMinuit.h>

#include "TauAnalysis/CandidateTools/interface/SVMethodLikelihoods.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"

#include "DataFormats/METReco/interface/MET.h"
//#include "DataFormats/PatCandidates/interface/Muon.h"
//#include "DataFormats/PatCandidates/interface/Electron.h"
//#include "DataFormats/PatCandidates/interface/Tau.h"

namespace TauVertex {

/// Find a valid initial position for a leg
std::pair<GlobalPoint, GlobalError> findValidSV(const GlobalPoint& pv, LegFitter* leg);

/// Pure abstract class to interface with Minuit
class DiTauFitterInterface : public TObject
{
   public:
      DiTauFitterInterface(){};
      virtual ~DiTauFitterInterface(){};
      /// Get nll given parameters
      virtual double nll(Double_t* pars, Int_t &status, bool verbose=true) = 0;
      /// Get MET NLL
      virtual double metNLL() const = 0;
      /// To set up parameters for the fit
      virtual void setupParametersFromScratch(TMinuit& minuit, bool useGeoLimits) = 0;
      virtual void setupParametersFromPrevious(const TMinuit& from, TMinuit& to, 
            bool useM12Limits, bool useGeoLimits) = 0;

      /// Fix the PV in the fit
      virtual void lockPV(TMinuit& minuit) = 0;
      /// Release the PV in the fit
      virtual void releasePV(TMinuit& minuit) = 0;
      /// Fix the Nu mass scales in the fit
      virtual void lockNuMassScalers(TMinuit& minuit) = 0;
      /// Release the nu mass scales in the fit
      virtual void releaseNuMassScalers(TMinuit& minuit) = 0;
      
      /// Enable/disable MET in the fit
      virtual void enableMET(bool enable=true)=0;

      // Access to the legs
      virtual const LegFitter* leg1() const = 0;
      virtual const LegFitter* leg2() const = 0;

      friend std::ostream& operator<< (std::ostream &out, DiTauFitterInterface &fit) { fit.printTo(out); return out; };
      virtual void printTo(std::ostream &out) const = 0;
};

template<typename T1, typename T2>
class DiTauFitter : public DiTauFitterInterface
{
   public:
      // Each nu-fit has two solutions, corresponding to the nu direction in the rest frame
      DiTauFitter(const T1& leg1, const std::vector<reco::TransientTrack>& leg1Tracks, 
            const T2& leg2, const std::vector<reco::TransientTrack>& leg2Tracks, 
            const TransientVertex& pv, const reco::MET* met, unsigned char soln):pv_(pv),met_(met)
      {
         // build each leg fitter
         leg1_ = std::auto_ptr<LegFitter>(new LegFitterSpecific<T1>(leg1, leg1Tracks, soln & 0x01));
         leg2_ = std::auto_ptr<LegFitter>(new LegFitterSpecific<T2>(leg2, leg2Tracks, soln & 0x02));
         enableMET_ = true;
      }

      const LegFitter* leg1() const {return leg1_.get(); };
      const LegFitter* leg2() const {return leg2_.get(); };

      /// Interface to Minuit
      double nll(Double_t *pars, Int_t &status, bool verbose)
      {
         // Update our status
         status = 0;

         // Make sure no parameters have blown up
         for(size_t ipar=0; ipar<11; ++ipar)
         {
            if(isnan(pars[ipar]) || isinf(pars[ipar]))
            {
               status |= 0x4;
               return 1e6;
            }
         }

         double nllOutput = 0;
         /* npars layout:  (PV x,y,z), (Leg1 x,y,z,m12) (Leg2 x,y,z,m12) */
         fitPV_ = GlobalPoint(pars[0], pars[1], pars[2]);

         /// Get the NLL of the PV
         pvNLL_ = nllPointGivenVertex(fitPV_, pv_);
         nllOutput += pvNLL_;

         int leg1Status=0; 
         int leg2Status=0;
         /// Get the NLL for each leg
         leg1_->setPoints(fitPV_, pars[3], pars[4], pars[5], pars[6], leg1Status);
         nllOutput += leg1_->nllOfLeg();
         leg2_->setPoints(fitPV_, pars[7], pars[8], pars[9], pars[10], leg2Status);
         nllOutput += leg2_->nllOfLeg();

         // Check if either leg failed
         if(leg1Status)
            status |= 0x1;
         if(leg2Status)
            status |= 0x2;

         lastStatus_ = status;

         // Constrain by MET
         nus_ = leg2_->nuP4() + leg1_->nuP4();
         metNLL_ = nllNuSystemGivenMET(nus_, met_);
         if(enableMET_)
            nllOutput += metNLL_;

         // Log NLL result
         lastNLL_ = nllOutput;

         return nllOutput;
      }

      double metNLL() const { return metNLL_; };

      /// Pretty print the current status
      void printTo(std::ostream &out) const
      {
         using namespace std;
         out << "=========== DiTauFitter =========== " << endl;
         out << setw(10) <<"Total NLL:" << setw(10) << lastNLL_ << setw(10) << "Status:" << setw(10) << lastStatus_ << endl;
         out << setw(10) <<"PV" << setw(10) << pvNLL_ << setw(10) << "Fit:" << setw(10) << fitPV_ << "Meas:" << setw(10) << pv_.position() << endl;
         if(enableMET_)
            out << setw(10) <<"MET" << setw(10) << metNLL_ << setw(10) << "Fit:" << setw(30) << nus_ << setw(10) << "Meas:" << setw(30) << met_->p4() << endl;
         else
            out << "MET disabled." << endl;
         out << " *** Leg 1 *** " << endl;
         leg1_->printTo(out);
         out << " *** Leg 2 *** " << endl;
         leg2_->printTo(out);
      }

      /// Setup minuit parameters using reasonable names
      void setupParametersFromPrevious(const TMinuit& from, TMinuit& to, bool useM12Limits, bool useGeoLimits)
      {
         double xyLimit = 0.;
         double zLimit = 0.;
         double m12Limits = 0.;
         if(useGeoLimits)
         {
            xyLimit = 3.0; //cm
            zLimit = 30.0; //cm
         }
         if(useM12Limits) 
         {
            m12Limits = 1.0;
         }

         Double_t par, parError;
         // ParNo, Name, InitValue, InitError, Limits
         from.GetParameter(0, par, parError);
         to.DefineParameter(0, "pv_x", par, parError, -1*xyLimit, xyLimit);
         from.GetParameter(1, par, parError);
         to.DefineParameter(1, "pv_y", par, parError, -1*xyLimit, xyLimit);
         from.GetParameter(2, par, parError);
         to.DefineParameter(2, "pv_z", par, parError, -1*zLimit, zLimit);
         
         from.GetParameter(3, par, parError);
         to.DefineParameter(3, "sv1_x", par, parError, -1*xyLimit, xyLimit);
         from.GetParameter(4, par, parError);
         to.DefineParameter(4, "sv1_y", par, parError, -1*xyLimit, xyLimit);
         from.GetParameter(5, par, parError);
         to.DefineParameter(5, "sv1_z", par, parError, -1*zLimit, zLimit);
         from.GetParameter(6, par, parError);
         to.DefineParameter(6, "mNuScale1", 
               (!fixM12ScaleParameter<T1>()) ? par : 0.0, 
               parError, -1., 1.); 

         from.GetParameter(7, par, parError);
         to.DefineParameter(7, "sv2_x", par, parError, -1*xyLimit, xyLimit);
         from.GetParameter(8, par, parError);
         to.DefineParameter(8, "sv2_y", par, parError, -1*xyLimit, xyLimit);
         from.GetParameter(9, par, parError);
         to.DefineParameter(9, "sv2_z", par, parError, -1*zLimit, zLimit);
         from.GetParameter(10, par, parError);
         to.DefineParameter(10, "mNuScale2", 
               (!fixM12ScaleParameter<T2>()) ? par : 0.0, 
               parError, -1., 1.);

         // Determine whether or not m12 is a free parameter and needs to be scaled.
         // For hadronic taus, there is only 1 nu, so it's mass is always zero and we 
         // can ignore m12
         if(fixM12ScaleParameter<T1>())
            to.FixParameter(6);
         if(fixM12ScaleParameter<T2>())
            to.FixParameter(10);
      }


      /// Setup minuit parameters using reasonable names
      void setupParametersFromScratch(TMinuit& minuit, bool useGeoLimits)
      {
         double xyLimit = 0.;
         double zLimit = 0.;
         if(useGeoLimits)
         {
            xyLimit = 3.0; //cm
            zLimit = 30.0; //cm
         }

         edm::LogInfo("SVMethod") << "Setting up initial PV";
         // ParNo, Name, InitValue, InitError, Limits
         minuit.DefineParameter(0, "pv_x", pv_.position().x(), pv_.positionError().cxx(), -1*xyLimit, xyLimit);
         minuit.DefineParameter(1, "pv_y", pv_.position().y(), pv_.positionError().cyy(), -1*xyLimit, xyLimit);
         minuit.DefineParameter(2, "pv_z", pv_.position().z(), pv_.positionError().czz(), -1*zLimit, zLimit);
         
         edm::LogInfo("SVMethod") << "Finding initial SV for leg 1";
         // Get initial guesses for each leg
         std::pair<GlobalPoint, GlobalError> leg1Init = findValidSV(pv_.position(), leg1_.get());
         GlobalPoint leg1pos = leg1Init.first;
         GlobalError leg1err = leg1Init.second;
         minuit.DefineParameter(3, "sv1_x", leg1pos.x(), leg1err.cxx(), -1*xyLimit, xyLimit);
         minuit.DefineParameter(4, "sv1_y", leg1pos.y(), leg1err.cyy(), -1*xyLimit, xyLimit);
         minuit.DefineParameter(5, "sv1_z", leg1pos.z(), leg1err.czz(), -1*zLimit, zLimit);
         minuit.DefineParameter(6, "mNuScale1", 
               (!fixM12ScaleParameter<T1>()) ? 0.7 : 0.0, 
               0.1, 0., 0.); 

         edm::LogInfo("SVMethod") << "Finding initial SV for leg 2";
         std::pair<GlobalPoint, GlobalError> leg2Init = findValidSV(pv_.position(), leg2_.get());
         GlobalPoint leg2pos = leg2Init.first;
         GlobalError leg2err = leg2Init.second;
         minuit.DefineParameter(7, "sv2_x", leg2pos.x(), leg2err.cxx(), -1*xyLimit, xyLimit);
         minuit.DefineParameter(8, "sv2_y", leg2pos.y(), leg2err.cyy(), -1*xyLimit, xyLimit);
         minuit.DefineParameter(9, "sv2_z", leg2pos.z(), leg2err.czz(), -1*zLimit, zLimit);
         minuit.DefineParameter(10, "mNuScale2", 
               (!fixM12ScaleParameter<T2>()) ? 0.7 : 0.0, 
               0.1, 0., 0.);

         // Determine whether or not m12 is a free parameter and needs to be scaled.
         // For hadronic taus, there is only 1 nu, so it's mass is always zero and we 
         // can ignore m12
         if(fixM12ScaleParameter<T1>())
            minuit.FixParameter(6);
         if(fixM12ScaleParameter<T2>())
            minuit.FixParameter(10);
      }

      /// Lock the PV variables in the fit
      void lockPV(TMinuit& minuit)
      {
         minuit.FixParameter(0);
         minuit.FixParameter(1);
         minuit.FixParameter(2);
      }

      /// Allow PV variables to float
      void releasePV(TMinuit& minuit)
      {
         minuit.Release(0);
         minuit.Release(1);
         minuit.Release(2);
      }

      /// Fix the Nu mass scale variables
      void lockNuMassScalers(TMinuit& minuit)
      {
         minuit.FixParameter(6);
         minuit.FixParameter(10);
      }

      /// Allow the Nu mass scale variables to float, if applicable
      void releaseNuMassScalers(TMinuit& minuit)
      {
         // Only let the ones that can actually be variable vary
         if(!fixM12ScaleParameter<T1>())
            minuit.Release(6);
         if(!fixM12ScaleParameter<T2>())
            minuit.Release(10);
      }

      void enableMET(bool enable)
      {
         enableMET_ = enable;
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
      std::auto_ptr<LegFitter> leg1_;
      std::auto_ptr<LegFitter> leg2_;

      bool enableMET_;

};
}
#endif
