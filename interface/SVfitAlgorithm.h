#ifndef TauAnalysis_CandidateTools_SVfitAlgorithm_h
#define TauAnalysis_CandidateTools_SVfitAlgorithm_h

/** \class SVfitAlgorithm
 *
 * Class used to reconstruct di-tau invariant mass via fitting/likelihood techniques.
 *
 * The fit is performed by passing log-likelihood functions to be minimized to Minuit;
 * the actual likelihood functions are implememted as plugins,
 * providing flexibility to reconstruct multiple di-tau invariant mass solutions
 * using different combinations of
 *  kinematic, missing transverse momentum, tau lepton lifetime,...
 * information.
 * 
 * \author Evan Friis, Christian Veelken; UC Davis
 *
 * \version $Revision: 1.6 $
 *
 * $Id: SVfitAlgorithm.h,v 1.6 2010/08/30 10:11:54 friis Exp $
 *
 */

#include "TauAnalysis/CandidateTools/interface/SVfitDiTauLikelihoodBase.h"
#include "TauAnalysis/CandidateTools/interface/SVfitEventVertexRefitter.h"
#include "TauAnalysis/CandidateTools/interface/SVfitLegTrackExtractor.h"
#include "TauAnalysis/CandidateTools/interface/svFitAuxFunctions.h"

#include "AnalysisDataFormats/TauAnalysis/interface/SVfitDiTauSolution.h"
#include "AnalysisDataFormats/TauAnalysis/interface/SVfitLegSolution.h"

#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

#include <TMinuit.h>
#include <TObject.h>

#include <Math/VectorUtil.h>

#include <vector>
#include <iostream>

// forward declaration of SVfitAlgorithm class
template<typename T1, typename T2> class SVfitAlgorithm;

namespace SVfitAlgorithm_namespace 
{
  // define function to be minimized by Minuit
  template <typename T1, typename T2>
  void objectiveFcn(Int_t& nParameter, Double_t*, Double_t& fcn, Double_t* parameters, Int_t errorFlag)
  {
    // retrieve the fit object
    SVfitAlgorithm<T1,T2>* algorithm = dynamic_cast<SVfitAlgorithm<T1,T2>*>(gMinuit->GetObjectFit());
    if ( !algorithm ) {
      edm::LogError("objectiveFcn") 
	<< "Call to gMinuit::GetObjectFit returned NULL pointer !!";
      return;
    }

    fcn = algorithm->negLogLikelihood(nParameter, parameters);    
  }
}

template<typename TT1, typename TT2>
class SVfitAlgorithm : public TObject
{
 public:
  enum fitParameter { kPrimaryVertexX, kPrimaryVertexY, kPrimaryVertexZ,
                      kLeg1thetaRest, kLeg1phiLab, kLeg1flightPathLab, kLeg1nuInvMass, 
                      kLeg2thetaRest, kLeg2phiLab, kLeg2flightPathLab, kLeg2nuInvMass };
  
  SVfitAlgorithm(const edm::ParameterSet& cfg)
    : currentDiTau_(0),
      minuit_(11)
  {
    name_ = cfg.getParameter<std::string>("name");

    eventVertexRefitAlgorithm_ = new SVfitEventVertexRefitter(cfg);

    likelihoodsSupportPolarization_ = false;

    typedef std::vector<edm::ParameterSet> vParameterSet;
    vParameterSet cfgLikelihoodFunctions = cfg.getParameter<vParameterSet>("likelihoodFunctions");
    for ( vParameterSet::const_iterator cfgLikelihoodFunction = cfgLikelihoodFunctions.begin();
	  cfgLikelihoodFunction != cfgLikelihoodFunctions.end(); ++cfgLikelihoodFunction ) {
      std::string pluginType = cfgLikelihoodFunction->getParameter<std::string>("pluginType");
      typedef edmplugin::PluginFactory<SVfitDiTauLikelihoodBase<TT1,TT2>* (const edm::ParameterSet&)> SVfitDiTauLikelihoodPluginFactory;
      SVfitDiTauLikelihoodBase<TT1,TT2>* likelihoodFunction 
	= SVfitDiTauLikelihoodPluginFactory::get()->create(pluginType, *cfgLikelihoodFunction);
      likelihoodsSupportPolarization_ |= likelihoodFunction->supportsPolarization();
      likelihoodFunctions_.push_back(likelihoodFunction);
    }
    
    std::cout << "<SVfitAlgorithm::SVfitAlgorithm>:" << std::endl;
    std::cout << " disabling MINUIT output..." << std::endl;
    minuit_.SetPrintLevel(-1);
    minuit_.SetErrorDef(0.5);
    
    minuitNumParameters_ = minuit_.GetNumPars();
    minuitParameterValues_ = new Double_t[minuitNumParameters_];
    minuitLockParameters_ = new bool[minuitNumParameters_];
    
//--- lock (i.e. set to fixed values) Minuit parameters
//    which are not constrained by any likelihood function
    for ( int iParameter = 0; iParameter < minuitNumParameters_; ++iParameter ) {
      minuitLockParameters_[iParameter] = true;
      
      for ( typename std::vector<SVfitDiTauLikelihoodBase<TT1,TT2>*>::const_iterator likelihoodFunction = likelihoodFunctions_.begin();
	    likelihoodFunction != likelihoodFunctions_.end(); ++likelihoodFunction ) {
	if ( (*likelihoodFunction)->isFittedParameter(iParameter) ) minuitLockParameters_[iParameter] = false;
      }
      
      if ( minuitLockParameters_[iParameter] ) minuit_.FixParameter(iParameter);
    }

    print(std::cout);
  }

  ~SVfitAlgorithm()
  {
    delete eventVertexRefitAlgorithm_;

    for ( typename std::vector<SVfitDiTauLikelihoodBase<TT1,TT2>*>::iterator it = likelihoodFunctions_.begin();
	  it != likelihoodFunctions_.end(); ++it ) {
      delete (*it);
    }

    delete [] minuitParameterValues_;
    delete [] minuitLockParameters_;
  }

  void beginEvent(edm::Event& evt, const edm::EventSetup& es)
  {
    for ( typename std::vector<SVfitDiTauLikelihoodBase<TT1,TT2>*>::const_iterator likelihoodFunction = likelihoodFunctions_.begin();
	  likelihoodFunction != likelihoodFunctions_.end(); ++likelihoodFunction ) {
      (*likelihoodFunction)->beginEvent(evt, es);
    }
  }

  void print(std::ostream& stream) const
  {
    stream << "<SVfitAlgorithm::print>" << std::endl;    
    stream << " name = " << name_ << std::endl;
    for ( typename std::vector<SVfitDiTauLikelihoodBase<TT1,TT2>*>::const_iterator likelihoodFunction = likelihoodFunctions_.begin();
	  likelihoodFunction != likelihoodFunctions_.end(); ++likelihoodFunction ) {
      (*likelihoodFunction)->print(stream);
      stream << std::endl;
    }
    for ( int iParameter = 0; iParameter < minuitNumParameters_; ++iParameter ) {
      stream << " Parameter #" << iParameter << ": ";
      if ( minuitLockParameters_[iParameter] ) stream << " LOCKED";
      else stream << " FITTED";
      stream << std::endl;
    }
    stream << std::endl;
  }

  std::vector<SVfitDiTauSolution> fit(const CompositePtrCandidateT1T2MEt<TT1,TT2>& diTauCandidate)
  {
    std::vector<SVfitDiTauSolution> solutions;
    
//--- refit primary event vertex
//    excluding tracks associated to tau decay products
    std::vector<reco::TrackBaseRef> leg1Tracks = leg1TrackExtractor_(*diTauCandidate.leg1());
    std::vector<reco::TrackBaseRef> leg2Tracks = leg2TrackExtractor_(*diTauCandidate.leg2());
    TransientVertex pv = eventVertexRefitAlgorithm_->refit(leg1Tracks, leg2Tracks);
    
    currentDiTau_ = &diTauCandidate;
    
    minuit_.SetFCN(SVfitAlgorithm_namespace::objectiveFcn<TT1,TT2>);
    minuit_.SetObjectFit(this);
    minuit_.SetMaxIterations(1000);
    gMinuit = &minuit_;
    
//--- check if at least one likelihood function supports polarization 
//    (i.e. returns a polarization depent likelihood value);
//    in case none of the likelihood functions support polarization,
//    run the fit only once (for "unknown" polarization), 
//    in order to save computing time
    if ( likelihoodsSupportPolarization_ ) {
      for ( int leg1PolarizationHypothesis = SVfitLegSolution::kLeftHanded; 
	    leg1PolarizationHypothesis <= SVfitLegSolution::kRightHanded; ++leg1PolarizationHypothesis ) {
	for ( int leg2PolarizationHypothesis = SVfitLegSolution::kLeftHanded; 
	      leg2PolarizationHypothesis <= SVfitLegSolution::kRightHanded; ++leg2PolarizationHypothesis ) {
	  currentDiTauSolution_ = SVfitDiTauSolution((SVfitLegSolution::polarizationHypothesisType)leg1PolarizationHypothesis, 
						     (SVfitLegSolution::polarizationHypothesisType)leg2PolarizationHypothesis);
	  currentDiTauSolution_.eventVertexPosition_.SetXYZ(pv.position().x(), pv.position().y(), pv.position().z());
	  currentDiTauSolution_.eventVertexPositionErr_ = pv.positionError().matrix();
	  fitPolarizationHypothesis(currentDiTauSolution_);
	  solutions.push_back(currentDiTauSolution_);
	} 
      }
    } else {
      currentDiTauSolution_ = SVfitDiTauSolution(SVfitLegSolution::kUnknown, SVfitLegSolution::kUnknown);
      currentDiTauSolution_.eventVertexPosition_.SetXYZ(pv.position().x(), pv.position().y(), pv.position().z());
      currentDiTauSolution_.eventVertexPositionErr_ = pv.positionError().matrix();
      fitPolarizationHypothesis(currentDiTauSolution_);
      solutions.push_back(currentDiTauSolution_);
    }
    
    return solutions;
  }
  
  double negLogLikelihood(Int_t nParameter, Double_t* parameters) const
  {
    if ( !currentDiTau_ ) {
      edm::LogError("SVfitAlgorithm::logLikelihood") 
	<< " Pointer to currentDiTau has not been initialized --> skipping !!";
    }
    
    readMinuitParameters();
    applyMinuitParameters(currentDiTauSolution_);
    
    double negLogLikelihood = 0.;
    
    for ( typename std::vector<SVfitDiTauLikelihoodBase<TT1,TT2>*>::const_iterator likelihoodFunction = likelihoodFunctions_.begin();
	  likelihoodFunction != likelihoodFunctions_.end(); ++likelihoodFunction ) {
      negLogLikelihood += (**likelihoodFunction)(*currentDiTau_, currentDiTauSolution_);
    }
    
    return negLogLikelihood;
  }
  
 private:

  void fitPolarizationHypothesis(SVfitDiTauSolution& solution)
  {
    int minuitStatus = minuit_.Command("MIN"); 
    edm::LogInfo("SVfitAlgorithm::fit") 
      << " Minuit fit Status = " << minuitStatus << std::endl;
	
    readMinuitParameters();
    applyMinuitParameters(currentDiTauSolution_);
    
    currentDiTauSolution_.minuitStatus_ = minuitStatus;
  }

  void readMinuitParameters() const
  {
    Double_t dummy;
    for ( Int_t iParameter = 0; iParameter < minuitNumParameters_; ++iParameter ) {
      minuit_.GetParameter(iParameter, minuitParameterValues_[iParameter], dummy);
    }
  }
  
  void applyMinuitParameters(SVfitDiTauSolution& diTauSolution) const 
  {
//--- set primary event vertex position (tau lepton production vertex)
    diTauSolution.eventVertexPositionCorr_.SetX(minuitParameterValues_[kPrimaryVertexX]);
    diTauSolution.eventVertexPositionCorr_.SetY(minuitParameterValues_[kPrimaryVertexY]);
    diTauSolution.eventVertexPositionCorr_.SetZ(minuitParameterValues_[kPrimaryVertexZ]);

    // Build leg 1
    applyMinuitParametersToLeg(kLeg1thetaRest, diTauSolution.leg1_);
    // Build leg 2
    applyMinuitParametersToLeg(kLeg2thetaRest, diTauSolution.leg2_);
  }

  void applyMinuitParametersToLeg(int parameterStart, SVfitLegSolution& legSolution) const
  {
    double gjAngle         = minuitParameterValues_[parameterStart + 0];
    double phiLab          = minuitParameterValues_[parameterStart + 1];
    double flightDistance_ = minuitParameterValues_[parameterStart + 2];
    double massNuNu        = minuitParameterValues_[parameterStart + 3];

    const reco::Candidate::LorentzVector& p4Vis = legSolution.p4Vis();

    // Compute the tau momentum in the rest frame
    double pVisRestFrame = SVfit_namespace::pVisRestFrame(p4Vis.mass(), massNuNu);
    // Get the opening angle in the lab frame
    double angleVisLabFrame = SVfit_namespace::gjAngleToLabFrame(pVisRestFrame, gjAngle, p4Vis.P());
    // Compute the tau momentum in the lab frame
    double momentumLabFrame = SVfit_namespace::tauMomentumLabFrame(p4Vis.mass(), pVisRestFrame, gjAngle, p4Vis.P());
    // Determine the direction of the tau
    reco::Candidate::Vector direction = SVfit_namespace::tauDirection(p4Vis.Vect().Unit(), angleVisLabFrame, phiLab);

    reco::Candidate::LorentzVector tauP4 = SVfit_namespace::tauP4(direction, momentumLabFrame);

    // Buid the tau four vector.  By construction, the neutrino is tauP4 - visP4
    legSolution.p4Invis_ =  tauP4 - p4Vis;

    // Build boost vector and compute the rest frame quanitites
    reco::Candidate::Vector boost = tauP4.BoostToCM();
    legSolution.p4VisRestFrame_ = ROOT::Math::VectorUtil::boost(legSolution.p4Vis_, boost);
    legSolution.p4InvisRestFrame_ = ROOT::Math::VectorUtil::boost(legSolution.p4Invis_, boost);

    // Set the flight path
    legSolution.tauFlightPath_ = direction*flightDistance_;
  }

  std::string name_;
  
  SVfitEventVertexRefitter* eventVertexRefitAlgorithm_;
  SVfitLegTrackExtractor<TT1> leg1TrackExtractor_;
  SVfitLegTrackExtractor<TT2> leg2TrackExtractor_;

  std::vector<SVfitDiTauLikelihoodBase<TT1,TT2>*> likelihoodFunctions_;
  bool likelihoodsSupportPolarization_;
  
  mutable const CompositePtrCandidateT1T2MEt<TT1,TT2>* currentDiTau_;
  mutable SVfitDiTauSolution currentDiTauSolution_;
  
  mutable TMinuit minuit_;
  Int_t minuitNumParameters_;
  mutable Double_t* minuitParameterValues_;
  bool* minuitLockParameters_;
};

#endif
