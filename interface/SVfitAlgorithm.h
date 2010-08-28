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
 * \version $Revision: 1.4 $
 *
 * $Id: SVfitAlgorithm.h,v 1.4 2010/08/27 17:21:50 veelken Exp $
 *
 */

#include "TauAnalysis/CandidateTools/interface/SVfitDiTauLikelihoodBase.h"

#include "AnalysisDataFormats/TauAnalysis/interface/SVfitDiTauSolution.h"
#include "AnalysisDataFormats/TauAnalysis/interface/SVfitLegSolution.h"

#include <TMinuit.h>
#include <TObject.h>

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

template<typename T1, typename T2>
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

    likelihoodsSupportPolarization_ = false;

    typedef std::vector<edm::ParameterSet> vParameterSet;
    vParameterSet cfgNegLogLikelihoodFunctions = cfg.getParameter<vParameterSet>("negLogLikelihoodFunctions");
    for ( vParameterSet::const_iterator cfgNegLogLikelihoodFunction = cfgNegLogLikelihoodFunctions.begin();
	  cfgNegLogLikelihoodFunction != cfgNegLogLikelihoodFunctions.end(); ++cfgNegLogLikelihoodFunction ) {
      std::string pluginType = cfgNegLogLikelihoodFunction->getParameter<std::string>("pluginType");
      typedef edmplugin::PluginFactory<SVfitDiTauLikelihoodBase<T1,T2>* (const edm::ParameterSet&)> SVfitDiTauLikelihoodPluginFactory;
      SVfitDiTauLikelihoodBase<T1,T2>* negLogLikelihoodFunction 
	= SVfitDiTauLikelihoodPluginFactory::get()->create(pluginType, *cfgNegLogLikelihoodFunction);
      likelihoodsSupportPolarization_ |= negLogLikelihoodFunction->supportsPolarization();
      negLogLikelihoodFunctions_.push_back(negLogLikelihoodFunction);
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
      
      for ( typename std::vector<SVfitDiTauLikelihoodBase<T1,T2>*>::const_iterator negLogLikelihoodFunction = negLogLikelihoodFunctions_.begin();
	    negLogLikelihoodFunction != negLogLikelihoodFunctions_.end(); ++negLogLikelihoodFunction ) {
	if ( (*negLogLikelihoodFunction)->isFittedParameter(iParameter) ) minuitLockParameters_[iParameter] = false;
      }
      
      if ( minuitLockParameters_[iParameter] ) minuit_.FixParameter(iParameter);
    }

    print(std::cout);
  }

  ~SVfitAlgorithm()
  {
    for ( typename std::vector<SVfitDiTauLikelihoodBase<T1,T2>*>::iterator it = negLogLikelihoodFunctions_.begin();
	  it != negLogLikelihoodFunctions_.end(); ++it ) {
      delete (*it);
    }

    delete [] minuitParameterValues_;
    delete [] minuitLockParameters_;
  }

  void print(std::ostream& stream) const
  {
    stream << "<SVfitAlgorithm::print>" << std::endl;    
    stream << " name = " << name_ << std::endl;
    for ( typename std::vector<SVfitDiTauLikelihoodBase<T1,T2>*>::const_iterator negLogLikelihoodFunction = negLogLikelihoodFunctions_.begin();
	  negLogLikelihoodFunction != negLogLikelihoodFunctions_.end(); ++negLogLikelihoodFunction ) {
      (*negLogLikelihoodFunction)->print(stream);
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

  std::vector<SVfitDiTauSolution> fit(const CompositePtrCandidateT1T2MEt<T1,T2>& diTau)
  {
    std::vector<SVfitDiTauSolution> solutions;
    
    currentDiTau_ = &diTau;
    
    minuit_.SetFCN(SVfitAlgorithm_namespace::objectiveFcn<T1,T2>);
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
	  fitPolarizationHypothesis(currentDiTauSolution_);
	  solutions.push_back(currentDiTauSolution_);
	} 
      }
    } else {
      currentDiTauSolution_ = SVfitDiTauSolution(SVfitLegSolution::kUnknown, SVfitLegSolution::kUnknown);
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
    
    for ( typename std::vector<SVfitDiTauLikelihoodBase<T1,T2>*>::const_iterator negLogLikelihoodFunction = negLogLikelihoodFunctions_.begin();
	  negLogLikelihoodFunction != negLogLikelihoodFunctions_.end(); ++negLogLikelihoodFunction ) {
      negLogLikelihood += (**negLogLikelihoodFunction)(*currentDiTau_, currentDiTauSolution_);
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

//--- set secondary vertex position (tau lepton decay vertex)

//--- set momenta of visible and invisible tau decay products

  }

  std::string name_;
  
  std::vector<SVfitDiTauLikelihoodBase<T1,T2>*> negLogLikelihoodFunctions_;
  bool likelihoodsSupportPolarization_;
  
  mutable const CompositePtrCandidateT1T2MEt<T1,T2>* currentDiTau_;
  mutable SVfitDiTauSolution currentDiTauSolution_;
  
  mutable TMinuit minuit_;
  Int_t minuitNumParameters_;
  mutable Double_t* minuitParameterValues_;
  bool* minuitLockParameters_;
};

#endif
