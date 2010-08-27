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
 * \version $Revision: 1.2 $
 *
 * $Id: SVfitAlgorithm.h,v 1.2 2009/05/26 12:36:29 veelken Exp $
 *
 */

#include "TauAnalysis/CandidateTools/interface/SVfitDiTauLikelihoodBase.h"

#include "AnalysisDataFormats/TauAnalysis/interface/SVfitDiTauSolution.h"
#include "AnalysisDataFormats/TauAnalysis/interface/SVfitLegSolution.h"

#include <TMinuit.h>
#include <TObject.h>

#include <vector>

namespace SVfitAlgorithm_namespace 
{
  // forward declaration of SVfitAlgorithm class
  template<typename T1, typename T2> class SVfitAlgorithm;

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

    fcn = algorithm->logLikelihood(nParameter, parameters);    
  }
}

template<typename T1, typename T2>
class SVfitAlgorithm : public TObject
{
 public:
  SVfitAlgorithm(const edm::ParameterSet& cfg)
    : currentDiTau_(0),
      currentDiTauSolution_(0),
      minuit_(11)
  {
    typedef std::vector<edm::ParameterSet> vParameterSet;
    vParameterSet cfgLogLikelihoodFunctions = cfg.getParameter<vParameterSet>("logLikelihoodFunctions");
    for ( vParameterSet::const_iterator cfgLogLikelihoodFunction = cfgLogLikelihoodFunctions.begin();
	  cfgLogLikelihoodFunction != cfgLogLikelihoodFunctions.end(); ++cfgLogLikelihoodFunction ) {
      std::string pluginType = cfgLogLikelihoodFunction->getParameter<std::string>("pluginType");
      typedef edmplugin::PluginFactory<SVfitDiTauLikelihoodBase<T1,T2>* (const edm::ParameterSet&)> SVfitDiTauLikelihoodPluginFactory;
      SVfitDiTauLikelihoodBase<T1,T2>* logLikelihoodFunction 
	= SVfitDiTauLikelihoodPluginFactory::get()->create(pluginType, *cfgLogLikelihoodFunction);
      logLikelihoodFunctions_.push_back(logLikelihoodFunction);
    }
    
    std::cout << "<SVfitAlgorithm::SVfitAlgorithm>:" << std::endl;
    std::cout << " disabling MINUIT output..." << std::endl;
    minuit_.SetPrintLevel(-1);
    minuit_.SetErrorDef(0.5);
    
    minuitNumParameters_ = minuit_.GetNumPars();
    minuitParameterValues_ = new Double_t[minuitNumParameters_];
  }

  ~SVfitAlgorithm()
  {
    for ( typename std::vector<SVfitDiTauLikelihoodBase<T1,T2>*>::iterator it = logLikelihoodFunctions_.begin();
	  it != logLikelihoodFunctions_.end(); ++it ) {
      delete (*it);
    }

    delete [] minuitParameterValues_;
  }

  std::vector<SVfitDiTauSolution> fit(const CompositePtrCandidateT1T2MEt<T1,T2>& diTau)
  {
    std::vector<SVfitDiTauSolution> solutions;
    
    currentDiTau_ = &diTau;
    
    minuit_.SetFCN(SVfitAlgorithm_namespace::objectiveFcn<T1,T2>);
    minuit_.SetObjectFit(this);
    minuit_.SetMaxIterations(1000);
    gMinuit = &minuit_;
    
    for ( int leg1PolarizationHypothesis = SVfitLegSolution::kLeftHanded; 
	  leg1PolarizationHypothesis <= SVfitLegSolution::kRightHanded; ++leg1PolarizationHypothesis ) {
      for ( int leg2PolarizationHypothesis = SVfitLegSolution::kLeftHanded; 
	    leg2PolarizationHypothesis <= SVfitLegSolution::kRightHanded; ++leg2PolarizationHypothesis ) {
	currentDiTauSolution_ = SVfitDiTauSolution(leg1PolarizationHypothesis, leg2PolarizationHypothesis);
	
	int minuitStatus = minuit_.Command("MIN"); 
	edm::LogInfo("SVfitAlgorithm::fit") 
	  << " Minuit fit Status = " << minuitStatus << std::endl;
	
	readMinuitParameters();
	applyMinuitParameters(currentDiTauSolution_);
	currentDiTauSolution_.minuitStatus_ = minuitStatus;
	
	solutions.push_back(currentDiTauSolution_);
      }
    }

    return solutions;
  }
  
  double logLikelihood(Int_t nParameter, Double_t* parameters) const 
  {
    if ( !currentDiTau_ ) {
      edm::LogError("SVfitAlgorithm::logLikelihood") 
	<< " Pointer to currentDiTau has not been initialized --> skipping !!";
    }
    
    readMinuitParameters();
    applyMinuitParameters(currentDiTauSolution_);
    
    double logLikelihood = 0.;
    
    for ( typename std::vector<SVfitDiTauLikelihoodBase<T1,T2>*>::const_iterator logLikelihoodFunction = logLikelihoodFunctions_.begin();
	  logLikelihoodFunction != logLikelihoodFunctions_.end(); ++logLikelihoodFunction ) {
      logLikelihood += (**logLikelihoodFunction)(*currentDiTau_, currentDiTauSolution_);
    }
    
    return logLikelihood;
  }
  
 private:
  void readMinuitParameters()
  {
    Double_t dummy;
    for ( Int_t iParameter = 0; iParameter < minuitNumParameters_; ++iParameter ) {
      minuit_.GetParameter(iParameter, minuitParameterValues_[iParameter], dummy);
    }
  }
  
  void applyMinuitParameters(SVfitDiTauSolution& diTauSolution)
  {
    std::cout << "blah" << std::endl;
  }
  
  std::vector<SVfitDiTauLikelihoodBase<T1,T2>*> logLikelihoodFunctions_;
  
  const CompositePtrCandidateT1T2MEt<T1,T2>* currentDiTau_;
  SVfitDiTauSolution currentDiTauSolution_;
  
  TMinuit minuit_;
  Int_t minuitNumParameters_;
  Double_t* minuitParameterValues_;
};

#endif
