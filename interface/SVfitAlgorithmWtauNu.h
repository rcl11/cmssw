#ifndef TauAnalysis_CandidateTools_SVfitAlgorithmWtauNu_h
#define TauAnalysis_CandidateTools_SVfitAlgorithmWtauNu_h

/** \class SVfitAlgorithmWtauNu
 *
 * Class used to reconstruct via fitting/likelihood techniques 
 * the transverse mass of tau lepton and neutrino produced in W boson decays
 *
 * The fit is performed by passing log-likelihood functions to be minimized to Minuit;
 * the actual likelihood functions are implememted as plugins,
 * providing flexibility to reconstruct multiple transverse mass solutions
 * using different combinations of
 *  kinematic, missing transverse momentum, tau lepton lifetime,...
 * information.
 *
 * \author Christian Veelken, UC Davis
 *
 * \version $Revision: 1.1 $
 *
 * $Id: SVfitAlgorithmWtauNu.h,v 1.1 2011/02/19 13:36:27 veelken Exp $
 *
 */

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "TauAnalysis/CandidateTools/interface/SVfitAlgorithm.h"
#include "TauAnalysis/CandidateTools/interface/SVfitWtauNuLikelihoodBase.h"
#include "TauAnalysis/CandidateTools/interface/SVfitEventVertexRefitter.h"
#include "TauAnalysis/CandidateTools/interface/SVfitVertexOnTrackFinder.h"
#include "TauAnalysis/CandidateTools/interface/SVfitLegTrackExtractor.h"
#include "TauAnalysis/CandidateTools/interface/svFitAuxFunctions.h"

#include "AnalysisDataFormats/TauAnalysis/interface/SVfitWtauNuSolution.h"
#include "AnalysisDataFormats/TauAnalysis/interface/SVfitLegSolution.h"

#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

#include <TFitterMinuit.h>
#include <Minuit2/FCNBase.h>
#include <TMath.h>
#include <TMatrixDSym.h>
#include <TDecompChol.h>
#include <TRandom3.h>

#include <vector>
#include <algorithm>
#include <iostream>

#include <boost/foreach.hpp>

// forward declaration of SVfitAlgorithmWtauNu class
template<typename T> class SVfitAlgorithmWtauNu;

namespace SVfitWtauNu_namespace
{
  enum fitParameter { 
    kNuPtLab = SVfit_namespace::kLeg2thetaRest, kNuPhiLab
  };
}

template<typename T>
class SVfitMinuitFCNadapterWtauNu : public ROOT::Minuit2::FCNBase
{
  public:
    SVfitMinuitFCNadapterWtauNu()
      : svFitAlgorithm_(0)
    {}
    ~SVfitMinuitFCNadapterWtauNu() {}

    void setSVfitAlgorithmWtauNu(const SVfitAlgorithmWtauNu<T>* svFitAlgorithm) { svFitAlgorithm_ = svFitAlgorithm; }

    /// define "objective" function called by Minuit
    double operator()(const std::vector<double>& x) const {
      return svFitAlgorithm_->negLogLikelihood(x);
    }

    /// define increase in "objective" function used by Minuit to determine one-sigma error contours;
    /// in case of (negative) log-likelihood functions, the value needs to be 0.5
    double Up() const { return 0.5; }
  private:
    const SVfitAlgorithmWtauNu<T>* svFitAlgorithm_;
};

template<typename T>
class SVfitAlgorithmWtauNu
{
  public:
    SVfitAlgorithmWtauNu(const edm::ParameterSet& cfg) 
      : currentWtauNuCandidate_(0)
    {
      verbosity_ = cfg.exists("verbosity") ? cfg.getParameter<int>("verbosity") : 0;

      name_ = cfg.getParameter<std::string>("name");
      if ( verbosity_ ) std::cout << "initializing SVfit algorithm with name = " << name_ << std::endl;

      eventVertexRefitAlgorithm_ = new SVfitEventVertexRefitter(cfg);

      typedef std::vector<edm::ParameterSet> vParameterSet;
      vParameterSet cfgLikelihoodFunctions = cfg.getParameter<vParameterSet>("likelihoodFunctions");
      for ( vParameterSet::const_iterator cfgLikelihoodFunction = cfgLikelihoodFunctions.begin();
	    cfgLikelihoodFunction != cfgLikelihoodFunctions.end(); ++cfgLikelihoodFunction ) {
        std::string pluginType = cfgLikelihoodFunction->getParameter<std::string>("pluginType");
        typedef edmplugin::PluginFactory<SVfitWtauNuLikelihoodBase<T>* (const edm::ParameterSet&)> SVfitWtauNuLikelihoodPluginFactory;
        SVfitWtauNuLikelihoodBase<T>* likelihoodFunction
	  = SVfitWtauNuLikelihoodPluginFactory::get()->create(pluginType, *cfgLikelihoodFunction);
        likelihoodFunctions_.push_back(likelihoodFunction);
      }

      minuitMaxInterations_ = cfg.exists("minuitMaxInterations") ?
	cfg.getParameter<unsigned>("minuitMaxInterations") : SVfit_namespace::defaultMinuitMaxInterations;

//--- initialize Minuit
      minuitFCNadapter_.setSVfitAlgorithmWtauNu(this);

      minuit_ = new TFitterMinuit();
      minuit_->SetMinuitFCN(&minuitFCNadapter_);
//--- set Minuit strategy = 2,
//    in order to enable reliable error estimates
//    ( cf. http://www-cdf.fnal.gov/physics/statistics/recommendations/minuit.html )
      minuit_->SetStrategy(2);
      minuit_->SetMaxIterations(minuitMaxInterations_);

      //std::cout << "<SVfitAlgorithmWtauNu::SVfitAlgorithmWtauNu>:" << std::endl;
      //std::cout << " disabling MINUIT output..." << std::endl;
      //minuit_->SetPrintLevel(1);
      minuit_->SetPrintLevel(-1);
      if ( verbosity_ ) minuit_->SetPrintLevel(0);
      minuit_->SetErrorDef(0.5);

      minuit_->CreateMinimizer();

      minuitFittedParameterValues_.resize(minuitNumParameters_);

      if ( cfg.exists("parameterizeVertexAlongTrackLeg1") ) {
	if ( verbosity_ ) std::cout << "--> enabling parameterizeVertexAlongTrackLeg1" << std::endl;
        parameterizeVertexAlongTrackLeg1_ = cfg.getParameter<bool>("parameterizeVertexAlongTrackLeg1");
      } else {
	if ( verbosity_ ) std::cout << "--> disabling parameterizeVertexAlongTrackLeg1" << std::endl;
        parameterizeVertexAlongTrackLeg1_ = false;
      }

      //print(std::cout);
    }

    ~SVfitAlgorithmWtauNu() 
    {
      delete eventVertexRefitAlgorithm_;

      for ( typename std::vector<SVfitWtauNuLikelihoodBase<T>*>::iterator it = likelihoodFunctions_.begin();
           it != likelihoodFunctions_.end(); ++it ) {
        delete (*it);
      }
    }

    void beginJob() {
      for ( typename std::vector<SVfitWtauNuLikelihoodBase<T>*>::const_iterator likelihoodFunction = likelihoodFunctions_.begin();
           likelihoodFunction != likelihoodFunctions_.end(); ++likelihoodFunction ) {
        (*likelihoodFunction)->beginJob();
      }
    }

    void beginEvent(edm::Event& evt, const edm::EventSetup& es) 
    {
      //std::cout << "<SVfitAlgorithmWtauNu::beginEvent>:" << std::endl;
      SVfit::track::VertexOnTrackFinder::setEventSetup(es);

      eventVertexRefitAlgorithm_->beginEvent(evt, es);

      for ( typename std::vector<SVfitWtauNuLikelihoodBase<T>*>::const_iterator likelihoodFunction = likelihoodFunctions_.begin();
           likelihoodFunction != likelihoodFunctions_.end(); ++likelihoodFunction ) {
        (*likelihoodFunction)->beginEvent(evt, es);
      }
    }

    void print(std::ostream& stream) const 
    {
      stream << "<SVfitAlgorithmWtauNu::print>" << std::endl;
      stream << " name = " << name_ << std::endl;
      for ( typename std::vector<SVfitWtauNuLikelihoodBase<T>*>::const_iterator likelihoodFunction = likelihoodFunctions_.begin();
           likelihoodFunction != likelihoodFunctions_.end(); ++likelihoodFunction ) {
        (*likelihoodFunction)->print(stream);
      }
      stream << " minuitNumParameters = " << minuitNumParameters_ << std::endl;
      stream << " numSamplings = " << numSamplings_ << std::endl;
      stream << std::endl;
    }

    std::vector<SVfitWtauNuSolution> fit(const CompositePtrCandidateTMEt<T>& candidate)
    {
      if ( verbosity_ ) {
        std::cout << "<SVfitAlgorithmWtauNu::fit>:" << std::endl;
        std::cout << " name = " << name_ << std::endl;
      }

      std::vector<SVfitWtauNuSolution> solutions;

//--- refit primary event vertex
//    excluding tracks associated to tau decay products
      std::vector<reco::TrackBaseRef> leg1Tracks = leg1TrackExtractor_(*candidate.visDecayProducts());
      TransientVertex pv = eventVertexRefitAlgorithm_->refit(&leg1Tracks);
      if ( verbosity_ ) {
        std::cout << " refitted event vertex (#tracks = " << pv.originalTracks().size() << "):"
		  << " x = " << pv.position().x() << " +/- " << TMath::Sqrt(pv.positionError().cxx()) << ","
		  << " y = " << pv.position().y() << " +/- " << TMath::Sqrt(pv.positionError().cyy()) << ","
		  << " z = " << pv.position().z() << " +/- " << TMath::Sqrt(pv.positionError().czz()) << std::endl;
      }

      currentWtauNuSolution_ = SVfitWtauNuSolution();
      // Associate the tracks
      currentWtauNuSolution_.leg1_.tracks_ = leg1Tracks;
      // Refit the vertices if possible
      currentWtauNuSolution_.leg1_.recoVertex_ = eventVertexRefitAlgorithm_->fitSecondaryVertex(leg1Tracks);

      fitHypothesis(candidate, currentWtauNuSolution_, pv);
      solutions.push_back(currentWtauNuSolution_);

      return solutions;
    }

    double negLogLikelihood(const std::vector<double>& x) const 
    {
      ++indexFitFunctionCall_;

      if ( verbosity_ ) {
        std::cout << "<SVfitAlgorithmWtauNu::negLogLikelihood>:" << std::endl;
        std::cout << " indexFitFunctionCall = " << indexFitFunctionCall_ << std::endl;
        for ( unsigned iParameter = 0; iParameter < minuitNumParameters_; ++iParameter ) {
          // only print unfixed parameters
          if ( !minuit_->IsFixed(iParameter) )
	    std::cout << " Parameter #" << iParameter << " "
              << minuit_->GetParName(iParameter) << " = "
              << x[iParameter] << std::endl;
        }
      }

      if ( !currentWtauNuCandidate_ ) {
        edm::LogError("SVfitAlgorithmWtauNu::logLikelihood")
	  << " Pointer to currentWtauNuCandidate has not been initialized --> skipping !!";
        return 0.;
      }


      applyParameters(currentWtauNuSolution_, x);

      double negLogLikelihood = 0.;
      for ( typename std::vector<SVfitWtauNuLikelihoodBase<T>*>::const_iterator likelihoodFunction = likelihoodFunctions_.begin();
           likelihoodFunction != likelihoodFunctions_.end(); ++likelihoodFunction ) {
        // Only use this likelihood if it is enabled in this iteration.
        if ( (*likelihoodFunction)->isFitted(fitIteration_) ) {
          negLogLikelihood += (**likelihoodFunction)(*currentWtauNuCandidate_, currentWtauNuSolution_);
        }
      }

      if ( verbosity_ ) std::cout << "--> negLogLikelihood = " << negLogLikelihood << ","
				  << " SVfit transverse mass = " << currentWtauNuSolution_.mt() << std::endl;

//-- in order to resolve ambiguities and improve convergence of the fit,
//   add "penalty" terms in case leg1phiLab, nuPhiLab parameters
//   are outside of the "physical" interval -pi..+pi
      double penalty = 0.;
      double leg1phiLab = x[SVfit_namespace::kLeg1phiLab];
      if ( TMath::Abs(leg1phiLab) > TMath::Pi() ) penalty += SVfit_namespace::square(TMath::Abs(leg1phiLab) - TMath::Pi());
      double nuPhiLab = x[SVfitWtauNu_namespace::kNuPhiLab];
      if ( TMath::Abs(nuPhiLab)   > TMath::Pi() ) penalty += SVfit_namespace::square(TMath::Abs(nuPhiLab)   - TMath::Pi());

      return negLogLikelihood + penalty;
    }

  private:

    unsigned int fitIteration_;

    void fitHypothesis(const CompositePtrCandidateTMEt<T>& candidate,
		       SVfitWtauNuSolution& solution,
		       const TransientVertex& pv)
    {
      if ( verbosity_ ) std::cout << "<SVfitAlgorithmWtauNu::fitHypothesis>:" << std::endl;

//--- initialize pointer to current candidate object
      currentWtauNuCandidate_ = &candidate;

//--- initialize data-members of solution object
      if ( pv.isValid() ) {
        currentWtauNuSolution_.eventVertexPosition_(0) = pv.position().x();
        currentWtauNuSolution_.eventVertexPosition_(1) = pv.position().y();
        currentWtauNuSolution_.eventVertexPosition_(2) = pv.position().z();
        currentWtauNuSolution_.eventVertexPositionErr_ = pv.positionError().matrix_new();
      }
      // CV: need to add protection against case that primary event vertex is not valid <-- FIXME ?
      currentWtauNuSolution_.eventVertexIsValid_ = pv.isValid();
      currentWtauNuSolution_.leg1_.p4Vis_ = candidate.visDecayProducts()->p4();
      currentWtauNuSolution_.nu_ = candidate.met()->p4();

      // Turn off vertex parameterization in case neutrals are present - leads
      // to situations with non-existent solutions.
      if ( leg1NeutralActivity_(*candidate.visDecayProducts()) ) {
        if ( verbosity_ )
          std::cout << "Disabling vertex parameterization by leg 1, it has neutrals."
		    << std::endl;
        parameterizeVertexAlongTrackLeg1_ = false;
      }

//--- initialize start-values of Minuit fit parameters
//
//    CV: how to deal with measurement errors in the visible momenta of the tau decay "leg"
//        when setting the maximum mass of the neutrino system ?
//

      double pvPositionXerr, pvPositionYerr, pvPositionZerr;
      if ( pv.isValid() ) {
        pvPositionXerr = TMath::Sqrt(pv.positionError().cxx());
        pvPositionYerr = TMath::Sqrt(pv.positionError().cyy());
        pvPositionZerr = TMath::Sqrt(pv.positionError().czz());
      } else {
        pvPositionXerr = 0.01;
        pvPositionYerr = 0.01;
        pvPositionZerr = 0.01;
      }

      minuit_->SetParameter(SVfit_namespace::kPrimaryVertexShiftX, "pv_dx", 0., pvPositionXerr,  -1.,  +1.);
      minuit_->SetParameter(SVfit_namespace::kPrimaryVertexShiftY, "pv_dy", 0., pvPositionYerr,  -1.,  +1.);
      minuit_->SetParameter(SVfit_namespace::kPrimaryVertexShiftZ, "pv_dz", 0., pvPositionZerr, -50., +50.);
      minuit_->SetParameter(SVfit_namespace::kLeg1thetaRest, "sv1_thetaRest", 0.25*TMath::Pi(), 0.5*TMath::Pi(), 0., TMath::Pi());
      double leg1PhiLabStepSize = ( parameterizeVertexAlongTrackLeg1_ ) ?
	TMath::Pi()/100. : TMath::Pi();
      minuit_->SetParameter(SVfit_namespace::kLeg1phiLab, "sv1_phiLab", 0., leg1PhiLabStepSize, 0., 0.); // do not set limits for phiLab

      double leg1decayDistanceLab0, leg1decayDistanceLabStepSize;
      if ( parameterizeVertexAlongTrackLeg1_ ) {
	leg1decayDistanceLab0 = 0.;
	leg1decayDistanceLabStepSize = 0.0100; // 100 microns
      } else {
	double gamma_cTauLifetime = (candidate.visDecayProducts()->energy()/SVfit_namespace::tauLeptonMass)*SVfit_namespace::cTauLifetime;
	if ( verbosity_ ) {
	  std::cout << "leg1: energy = " << candidate.visDecayProducts()->energy()
		    << " --> gamma = " << (candidate.visDecayProducts()->energy()/SVfit_namespace::tauLeptonMass) << std::endl;
	  std::cout << " gamma_cTauLifetime = " << gamma_cTauLifetime << std::endl;
	}

	// do not set limits for decayDistanceLab; guarantee that decayDistanceLab is positive by taking square-root as fit parameter
	leg1decayDistanceLab0 = TMath::Sqrt(gamma_cTauLifetime);
	leg1decayDistanceLabStepSize = 0.1*leg1decayDistanceLab0;
      }
      minuit_->SetParameter(SVfit_namespace::kLeg1decayDistanceLab, "sv1_decayDistanceLab",
			    leg1decayDistanceLab0, leg1decayDistanceLabStepSize, 0., 0.);

      double leg1NuMass0, leg1NuMassErr, leg1NuMassMax;
      if ( !SVfit_namespace::isMasslessNuSystem<T>() ) {
        leg1NuMass0 = 0.8;
        leg1NuMassErr = 0.4;
        leg1NuMassMax = SVfit_namespace::tauLeptonMass - candidate.visDecayProducts()->mass();
      } else {
        leg1NuMass0 = 0.;
        leg1NuMassErr = 1.;
        leg1NuMassMax = 0.;
      }
      if ( verbosity_ ) std::cout << " leg1NuMass0 = " << leg1NuMass0 << ", leg1NuMassMax = " << leg1NuMassMax << std::endl;
      minuit_->SetParameter(SVfit_namespace::kLeg1nuInvMass, "sv1_m12", leg1NuMass0, leg1NuMassErr, 0., leg1NuMassMax);
      minuit_->SetParameter(SVfit_namespace::kLeg1thetaVMrho, "sv1_thetaVMrho", 0.25*TMath::Pi(), 0.5*TMath::Pi(), 0., TMath::Pi());
      minuit_->SetParameter(SVfit_namespace::kLeg1thetaVMa1, "sv1_thetaVMa1", 0.25*TMath::Pi(), 0.5*TMath::Pi(), 0., TMath::Pi());
      minuit_->SetParameter(SVfit_namespace::kLeg1thetaVMa1r, "sv1_thetaVMa1r", 0.25*TMath::Pi(), 0.5*TMath::Pi(), 0., TMath::Pi());
      minuit_->SetParameter(SVfit_namespace::kLeg1phiVMa1r, "sv1_phiVMa1r", 0., TMath::Pi(), 0., 0.); // do not set limits for phiVMa1r

      minuit_->SetParameter(SVfitWtauNu_namespace::kNuPtLab, "nu_ptLab", candidate.met()->pt(), 5., 0., 0.); 
      minuit_->SetParameter(SVfitWtauNu_namespace::kNuPhiLab, "nu_phiLab", 0., 0.1*TMath::Pi(), 0., 0.); // do not set limits for phiLab

      unsigned int fitIterations = 0;

      // The list of all parameters that will be fitted.
      std::set<std::string> fittedParameterNames;
      for ( typename std::vector<SVfitWtauNuLikelihoodBase<T>*>::const_iterator likelihoodFunction = likelihoodFunctions_.begin();
	    likelihoodFunction != likelihoodFunctions_.end(); ++likelihoodFunction ) {
        // Load information about the candidate
        (*likelihoodFunction)->beginCandidate(candidate);
        // Figure out how many fit iterations our selected likelihoods need
        if ((*likelihoodFunction)->firstFit() > fitIterations) {
          fitIterations = (*likelihoodFunction)->firstFit();
        }
        // Get the list of all fitted parameters names (for all iterations)
        // We use this to book our monitoring histograms.
        for ( unsigned iParameter = 0; iParameter < minuitNumParameters_; ++iParameter ) {
          if ((*likelihoodFunction)->isFittedParameter(iParameter)) {
            fittedParameterNames.insert(minuit_->GetParName(iParameter));
          }
        }
      }

      minuitStatus_ = -20;
      minuitEdm_ = 0;
      // Loop over the desired number of fit iterations.
      indexFitFunctionCall_ = 0;
      for ( fitIteration_ = 0; fitIteration_ < (fitIterations + 1); ++fitIteration_ ) {
//--- lock (i.e. set to fixed values) Minuit parameters
//    which are not constrained by any likelihood function (in this iteration)
        for ( unsigned iParameter = 0; iParameter < minuitNumParameters_; ++iParameter ) {
          bool minuitLockParameter = true;
          for ( typename std::vector<SVfitWtauNuLikelihoodBase<T>*>::const_iterator likelihoodFunction = likelihoodFunctions_.begin();
		likelihoodFunction != likelihoodFunctions_.end(); ++likelihoodFunction ) {
            if ( (*likelihoodFunction)->isFittedParameter(iParameter) &&
		 (*likelihoodFunction)->isFitted(fitIteration_) ) {
              minuitLockParameter = false;
              break;
            }
          }

          if (  minuitLockParameter && !minuit_->IsFixed(iParameter) ) minuit_->FixParameter(iParameter);
          if ( !minuitLockParameter &&  minuit_->IsFixed(iParameter) ) minuit_->ReleaseParameter(iParameter);

	  if ( verbosity_ ) {
            std::cout << " Parameter #" << iParameter << " (" << minuit_->GetParName(iParameter) << "): ";
            if ( minuit_->IsFixed(iParameter) ) std::cout << "LOCKED";
            else std::cout << "FITTED";
            std::cout << std::endl;
          }
        }

        minuitNumFreeParameters_ = minuit_->GetNumberFreeParameters();
        minuitNumFixedParameters_ = minuit_->GetNumberTotalParameters() - minuitNumFreeParameters_;

        if ( verbosity_ ) {
          std::cout << " minuitNumParameters = " << minuit_->GetNumberTotalParameters()
		    << " (free = " << minuitNumFreeParameters_ << ", fixed = " << minuitNumFixedParameters_ << ")" << std::endl;
        }
        assert((minuitNumFreeParameters_ + minuitNumFixedParameters_) == minuitNumParameters_);

        if ( verbosity_ ) std::cout << " BEGINNING FIT ITERATION #" << fitIteration_ << std::endl;

        minuitStatus_ = minuit_->Minimize();

        // Get stats about the minimizer;
        Double_t dummyErrDef;
        Int_t dummyNpvar, dummyNparx;
        minuit_->GetStats(minuitNll_, minuitEdm_, dummyErrDef, dummyNpvar, dummyNparx);

        if ( verbosity_ ) {
          std::cout << " DONE WITH FIT ITERATION #" << fitIteration_
		    << " FIT RESULT: " << minuitStatus_
		    << " NLL: " << minuitNll_ << " EDM: " << minuitEdm_ << std::endl;
          std::cout << " COVARIANCE MATRIX:" << std::endl;
        }
      }

      edm::LogInfo("SVfitAlgorithmWtauNu::fit") << " Minuit fit Status = " << minuitStatus_ << std::endl;

      for ( unsigned iParameter = 0; iParameter < minuitNumParameters_; ++iParameter ) {
        minuitFittedParameterValues_[iParameter] = minuit_->GetParameter(iParameter);
        if ( verbosity_ ) {
          std::cout << " Parameter #" << iParameter <<  ":"
		    << " " << minuit_->GetParName(iParameter) << " = "
		    << minuitFittedParameterValues_[iParameter] << std::endl;
        }
      }

      applyParameters(currentWtauNuSolution_, minuitFittedParameterValues_);

      for ( typename std::vector<SVfitWtauNuLikelihoodBase<T>*>::const_iterator likelihoodFunction = likelihoodFunctions_.begin();
           likelihoodFunction != likelihoodFunctions_.end(); ++likelihoodFunction ) {
        double likelihoodFunctionValue = (**likelihoodFunction)(*currentWtauNuCandidate_, currentWtauNuSolution_);
        currentWtauNuSolution_.negLogLikelihoods_.insert(std::make_pair((*likelihoodFunction)->name(), likelihoodFunctionValue));
      }

      currentWtauNuSolution_.minuitStatus_ = minuitStatus_;
    }

    void applyParameters(SVfitWtauNuSolution& solution, const std::vector<double>& x) const
    {
//--- set primary event vertex position (tau lepton production vertex)
      solution.eventVertexPositionShift_(0) = x[SVfit_namespace::kPrimaryVertexShiftX];
      solution.eventVertexPositionShift_(1) = x[SVfit_namespace::kPrimaryVertexShiftY];
      solution.eventVertexPositionShift_(2) = x[SVfit_namespace::kPrimaryVertexShiftZ];

//--- build tau decay "leg"
      applyParametersToLeg(solution.leg1_,
			   x, solution.eventVertexPosSVrefitted(),
			   parameterizeVertexAlongTrackLeg1_);

//--- build neutrino
      applyParametersToNeutrino(solution.nu_, x);
    }

    void applyParametersToLeg(SVfitLegSolution& legSolution,
        const std::vector<double>& x,
        const AlgebraicVector3& eventVertexPos,
        bool parameterizeVertexAlongTrack) const
    {
      double gjAngle  = x[SVfit_namespace::kLeg1thetaRest];
      double phiLab   = x[SVfit_namespace::kLeg1phiLab];
      double massNuNu = x[SVfit_namespace::kLeg1nuInvMass];

      const reco::Candidate::LorentzVector& p4Vis = legSolution.p4Vis();

      // Compute the tau momentum in the rest frame
      double pVisRestFrame = SVfit_namespace::pVisRestFrame(p4Vis.mass(), massNuNu);
      // Get the opening angle in the lab frame
      double angleVisLabFrame = SVfit_namespace::gjAngleToLabFrame(pVisRestFrame, gjAngle, p4Vis.P());
      // Compute the tau momentum in the lab frame
      double momentumLabFrame = SVfit_namespace::tauMomentumLabFrame(p4Vis.mass(), pVisRestFrame, gjAngle, p4Vis.P());

      reco::Candidate::Vector direction;

//--- if parameterizeVertexAlongTrack is enabled, compute approximate position
//    of tau lepton decay vertex as intersection of tau momentum vector with track
//    and take phiLab and decayDistanceToTrack parameters as corrections relative to that position,
//    else compute tau lepton decay vertex from thetaRest, phiLab and decayDistanceLab parameters directly
      if ( !parameterizeVertexAlongTrack ) {
	double flightDistance = SVfit_namespace::square(x[SVfit_namespace::kLeg1decayDistanceLab]);
	//std::cout << " flightDistance = " << flightDistance << std::endl;

        // Determine the direction of the tau
        direction = SVfit_namespace::tauDirection(p4Vis.Vect().Unit(), angleVisLabFrame, phiLab);
        //std::cout << " direction: x = " << direction.x() << ", y = " << direction.y() << ", z = " << direction.z() << std::endl;

        // Set the flight path
        legSolution.decayVertexPos_(0) = eventVertexPos(0) + flightDistance*direction.x();
        legSolution.decayVertexPos_(1) = eventVertexPos(1) + flightDistance*direction.y();
        legSolution.decayVertexPos_(2) = eventVertexPos(2) + flightDistance*direction.z();
      } else {
        // Find the approximate decay vertex given the track information.
	double flightDistanceCorr = x[SVfit_namespace::kLeg1decayDistanceLab];

        SVfit::track::VertexOnTrackFinder svFinder(legSolution);
        // Get the reconstructed vertex and the direction
        GlobalPoint secondaryVertexOnTrack = svFinder.decayVertex(
            eventVertexPos, angleVisLabFrame,
            phiLab, flightDistanceCorr, &direction);

        legSolution.decayVertexPos_(0) = secondaryVertexOnTrack.x();
        legSolution.decayVertexPos_(1) = secondaryVertexOnTrack.y();
        legSolution.decayVertexPos_(2) = secondaryVertexOnTrack.z();

        if ( verbosity_ ) std::cout << "Found vertex via track @ " << secondaryVertexOnTrack << std::endl;
      }

      // Throw a warning if the tau direction is opposite to that of the visible
      // momentum.
      if (direction.Unit().Dot(p4Vis.Vect().Unit()) < 0) {
        edm::LogWarning("BadTauDirection")
          << "The computed tau direction points opposite to the vis momentum!"
          << " tau: " << direction.Unit() << " vis: " << p4Vis.Vect().Unit();
      }

      // Compute the tau P4
      reco::Candidate::LorentzVector tauP4 = SVfit_namespace::tauP4(direction.Unit(), momentumLabFrame);
      if ( verbosity_ ) std::cout << "tauP4: p = " << tauP4.P() << ","
				  << " eta = " << tauP4.eta() << ", phi = " << tauP4.phi() << std::endl;

      // Build the tau four vector. By construction, the neutrino is tauP4 - visP4
      legSolution.p4Invis_ = tauP4 - p4Vis;

      // Build boost vector and compute the rest frame quanitites
      legSolution.p4VisRestFrame_ = SVfit_namespace::boostToCOM(tauP4, legSolution.p4Vis_);
      legSolution.p4InvisRestFrame_ = SVfit_namespace::boostToCOM(tauP4, legSolution.p4Invis_);

      // Set meson decay angles for tau- --> rho- nu --> pi- pi0 nu
      // and tau- --> a1- nu --> pi- pi0 pi0 nu, tau- --> a1- nu --> pi- pi+ pi- nu decay modes
      // (needed in case likelihood functions for decays of polarized tau leptons are used)
      legSolution.thetaVMrho_ = x[SVfit_namespace::kLeg1thetaVMrho];
      legSolution.thetaVMa1_  = x[SVfit_namespace::kLeg1thetaVMa1];
      legSolution.thetaVMa1r_ = x[SVfit_namespace::kLeg1thetaVMa1r];
      legSolution.phiVMa1r_   = x[SVfit_namespace::kLeg1phiVMa1r];
    }

    void applyParametersToNeutrino(
        reco::Candidate::LorentzVector& nuSolution,
	const std::vector<double>& x) const
    {
      double pt  = x[SVfitWtauNu_namespace::kNuPtLab];
      double phi = x[SVfitWtauNu_namespace::kNuPhiLab];

      double px  = pt*TMath::Cos(phi); 
      double py  = pt*TMath::Sin(phi); 

      nuSolution = reco::Candidate::LorentzVector(px, py, 0., pt);
    }

    std::string name_;
    bool parameterizeVertexAlongTrackLeg1_;

    SVfitEventVertexRefitter* eventVertexRefitAlgorithm_;
    SVfitLegTrackExtractor<T> leg1TrackExtractor_;
    SVfitLegHasNeutralsExtractor<T> leg1NeutralActivity_;

    std::vector<SVfitWtauNuLikelihoodBase<T>*> likelihoodFunctions_;

    mutable const CompositePtrCandidateTMEt<T>* currentWtauNuCandidate_;
    mutable SVfitWtauNuSolution currentWtauNuSolution_;

    mutable TFitterMinuit* minuit_;
    SVfitMinuitFCNadapterWtauNu<T> minuitFCNadapter_;
    unsigned minuitMaxInterations_;
    const static unsigned minuitNumParameters_ = 13;
    mutable unsigned minuitNumFreeParameters_;
    mutable unsigned minuitNumFixedParameters_;
    mutable std::vector<double> minuitFittedParameterValues_;

    int numSamplings_;
    mutable TRandom3 rnd_;

    mutable long indexFitFunctionCall_;

    // Class level variables for extracting monitoring from minuit.
    int minuitStatus_;
    Double_t minuitEdm_, minuitNll_;

    int verbosity_;
};

#endif
