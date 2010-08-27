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
 * \version $Revision: 1.1 $
 *
 * $Id: SVfitAlgorithm.h,v 1.1 2010/08/27 07:00:09 veelken Exp $
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
  enum fitParameter { kPrimaryVertexX, kPrimaryVertexY, kPrimaryVertexZ };

  void readMinuitParameters()
  {
    Double_t dummy;
    for ( Int_t iParameter = 0; iParameter < minuitNumParameters_; ++iParameter ) {
      minuit_.GetParameter(iParameter, minuitParameterValues_[iParameter], dummy);
    }
  }
  
  void applyMinuitParameters(SVfitDiTauSolution& diTauSolution)
  {
//--- set primary event vertex position (tau lepton production vertex)
    diTauSolution.eventVertexPositionCorr_.SetX(minuitParameterValues_[kPrimaryVertexX]);
    diTauSolution.eventVertexPositionCorr_.SetY(minuitParameterValues_[kPrimaryVertexY]);
    diTauSolution.eventVertexPositionCorr_.SetZ(minuitParameterValues_[kPrimaryVertexZ]);

//--- set secondary vertex position (tau lepton decay vertex)
/*
   // Determine lab frame opening angle between tau direction
   // and vis. momentum
   //
   // pl_perp = pr(m12)*sin(theta_r) ==>
   // pl*sin(thetal) = pr(m12)*sin(theta_r)
   // thetal = asin(pr(m12)*sin(theta_r)/pl)

   double thetaLab = TMath::ASin(restFrameVisMomentum()*
               TMath::Sin(thetaRest_)/visible->p4().P());

   // Build our displacement vector assuming visP parallel to Z axis
   reco::Candidate::Vector secondaryVertexDirection(
         TMath::Sin(thetaLab)*TMath::Cos(phiLab_),
         TMath::Sin(thetaLab)*TMath::Sin(phiLab_),
         TMath::Cos(thetaLab));

   // Rotate such that Z is along visible momentum
   reco::Candidate::Vector flight = 
      rotateUz(secondaryVertexDirection, visible->p4().P().Unit()) *
      flightDistance_;

   // Set the decay vertices of the two objects
   visible()->setVertex(this->vertex() + flight);
   invisible()->setVertex(this->vertex() + flight);

   // Determine tau direction
   reco::Particle::Vector tauFlight = visible()->vertex() - this->vertex();
   reco::Particle::Vector tauDir = tauFlight.Unit();

   // Visible lab frame momentum par/perp to tau
   double labFramePPerp = visible()->p4().Vect().Cross(tauDir).R();
   double labFramePParallel = visible()->p4().Vect().Dot(tauDir);

   // Find the rest frame momentum of the visible stuff
   // pLabPerp = pRestPerp = pRest*sin(theta) ==> pRest = pLabPerp/sin(theta)
   double restFrameP = labFramePPerp/TMath::Sin(thetaRest_);
   // The parallel component
   double restFramePParallel = restFrameP*TMath::Cos(thetaRest_);
   double restFrameE = energy(visible()->mass(), restFrameP);

   double gamma = (
         restFrameE * TMath::Sqrt(square(restFrameE) + square(labFramePParallel) - square(restFramePParallel)) -
         restFramePParallel*labFramePParallel ) /
         (square(restFrameE) - square(restFramePParallel));

   double tauEnergy = gamma*tauMass;
   double tauMomentum = momentum(tauMass, tauEnergy);

   reco::Candidate::LorentzVector fourVector(
         tauDir.X()*tauMomentum,
         tauDir.Y()*tauMomentum,
         tauDir.Z()*tauMomentum,
         tauEnergy);

   // Set the total p4 
   setP4(fourVector);
   // Set the neutrino p4.  Difference between tau and visible.
   invisible()->setP4(fourVector - this->visible()->p4());
 */
  }
  
  std::vector<SVfitDiTauLikelihoodBase<T1,T2>*> logLikelihoodFunctions_;
  
  const CompositePtrCandidateT1T2MEt<T1,T2>* currentDiTau_;
  SVfitDiTauSolution currentDiTauSolution_;
  
  TMinuit minuit_;
  Int_t minuitNumParameters_;
  Double_t* minuitParameterValues_;
};

#endif
