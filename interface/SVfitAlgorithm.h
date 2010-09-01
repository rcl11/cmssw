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
 * \version $Revision: 1.7 $
 *
 * $Id: SVfitAlgorithm.h,v 1.7 2010/08/30 13:29:45 veelken Exp $
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
#include <TMath.h>

#include <Math/VectorUtil.h>

#include <vector>
#include <iostream>

// forward declaration of SVfitAlgorithm class
template<typename T1, typename T2> class SVfitAlgorithm;

//
// WARNING: numbering of TMinuit parameters starts at 1 (Fortran convention);
//          an offset of 1 needs hence to be added to all calls to the TMinuit functions
//         o DefineParameter
//         o FixParameter
//         o GetParameter
//
const int minuitParameterOffset = 1;

namespace SVfitAlgorithm_namespace 
{
  // define function to be minimized by Minuit
  template <typename T1, typename T2>
  void objectiveFcn(Int_t& numParameters, Double_t*, Double_t& fcn, Double_t* parameter, Int_t errorFlag)
  {
    // retrieve the fit object
    SVfitAlgorithm<T1,T2>* algorithm = dynamic_cast<SVfitAlgorithm<T1,T2>*>(gMinuit->GetObjectFit());
    if ( !algorithm ) {
      edm::LogError("objectiveFcn") 
	<< "Call to gMinuit::GetObjectFit returned NULL pointer !!";
      return;
    }

    fcn = algorithm->negLogLikelihood(numParameters, parameter);    
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
      minuit_(minuitNumParameters_)
  {
    name_ = cfg.getParameter<std::string>("name");

    eventVertexRefitAlgorithm_ = new SVfitEventVertexRefitter(cfg);

    likelihoodsSupportPolarization_ = false;

    typedef std::vector<edm::ParameterSet> vParameterSet;
    vParameterSet cfgLikelihoodFunctions = cfg.getParameter<vParameterSet>("likelihoodFunctions");
    for ( vParameterSet::const_iterator cfgLikelihoodFunction = cfgLikelihoodFunctions.begin();
	  cfgLikelihoodFunction != cfgLikelihoodFunctions.end(); ++cfgLikelihoodFunction ) {
      std::string pluginType = cfgLikelihoodFunction->getParameter<std::string>("pluginType");
      typedef edmplugin::PluginFactory<SVfitDiTauLikelihoodBase<T1,T2>* (const edm::ParameterSet&)> SVfitDiTauLikelihoodPluginFactory;
      SVfitDiTauLikelihoodBase<T1,T2>* likelihoodFunction 
	= SVfitDiTauLikelihoodPluginFactory::get()->create(pluginType, *cfgLikelihoodFunction);
      likelihoodsSupportPolarization_ |= likelihoodFunction->supportsPolarization();
      likelihoodFunctions_.push_back(likelihoodFunction);
    }
    
    std::cout << "<SVfitAlgorithm::SVfitAlgorithm>:" << std::endl;
    std::cout << " disabling MINUIT output..." << std::endl;
    minuit_.SetPrintLevel(-1);
    minuit_.SetErrorDef(0.5);
    
    minuitParameterValues_ = new Double_t[minuitNumParameters_];
    minuitLockParameters_ = new bool[minuitNumParameters_];
    
//--- lock (i.e. set to fixed values) Minuit parameters
//    which are not constrained by any likelihood function
    for ( int iParameter = 0; iParameter < minuitNumParameters_; ++iParameter ) {
      minuitLockParameters_[iParameter] = true;
      
      for ( typename std::vector<SVfitDiTauLikelihoodBase<T1,T2>*>::const_iterator likelihoodFunction = likelihoodFunctions_.begin();
	    likelihoodFunction != likelihoodFunctions_.end(); ++likelihoodFunction ) {
	if ( (*likelihoodFunction)->isFittedParameter(iParameter) ) minuitLockParameters_[iParameter] = false;
      }
      
      if ( minuitLockParameters_[iParameter] ) minuit_.FixParameter(iParameter + minuitParameterOffset);
    }

    print(std::cout);
  }

  ~SVfitAlgorithm()
  {
    delete eventVertexRefitAlgorithm_;

    for ( typename std::vector<SVfitDiTauLikelihoodBase<T1,T2>*>::iterator it = likelihoodFunctions_.begin();
	  it != likelihoodFunctions_.end(); ++it ) {
      delete (*it);
    }

    delete [] minuitParameterValues_;
    delete [] minuitLockParameters_;
  }

  void beginEvent(edm::Event& evt, const edm::EventSetup& es)
  {
    eventVertexRefitAlgorithm_->beginEvent(evt, es);

    for ( typename std::vector<SVfitDiTauLikelihoodBase<T1,T2>*>::const_iterator likelihoodFunction = likelihoodFunctions_.begin();
	  likelihoodFunction != likelihoodFunctions_.end(); ++likelihoodFunction ) {
      (*likelihoodFunction)->beginEvent(evt, es);
    }
  }

  void print(std::ostream& stream) const
  {
    stream << "<SVfitAlgorithm::print>" << std::endl;    
    stream << " name = " << name_ << std::endl;
    for ( typename std::vector<SVfitDiTauLikelihoodBase<T1,T2>*>::const_iterator likelihoodFunction = likelihoodFunctions_.begin();
	  likelihoodFunction != likelihoodFunctions_.end(); ++likelihoodFunction ) {
      (*likelihoodFunction)->print(stream);
    }
    stream << " minuitNumParameters = " << minuitNumParameters_ << std::endl;
    for ( int iParameter = 0; iParameter < minuitNumParameters_; ++iParameter ) {
      stream << " Parameter #" << iParameter << ": ";
      if ( minuitLockParameters_[iParameter] ) stream << "LOCKED";
      else stream << "FITTED";
      stream << std::endl;
    }
    stream << std::endl;
  }

  std::vector<SVfitDiTauSolution> fit(const CompositePtrCandidateT1T2MEt<T1,T2>& diTauCandidate)
  {
    //std::cout << "<SVfitAlgorithm::fit>:" << std::endl;

    std::vector<SVfitDiTauSolution> solutions;
    
//--- refit primary event vertex
//    excluding tracks associated to tau decay products
    std::vector<reco::TrackBaseRef> leg1Tracks = leg1TrackExtractor_(*diTauCandidate.leg1());
    std::vector<reco::TrackBaseRef> leg2Tracks = leg2TrackExtractor_(*diTauCandidate.leg2());
    TransientVertex pv = eventVertexRefitAlgorithm_->refit(leg1Tracks, leg2Tracks);

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
	  fitPolarizationHypothesis(diTauCandidate, currentDiTauSolution_, pv);
	  solutions.push_back(currentDiTauSolution_);
	} 
      }
    } else {
      currentDiTauSolution_ = SVfitDiTauSolution(SVfitLegSolution::kUnknown, SVfitLegSolution::kUnknown);
      fitPolarizationHypothesis(diTauCandidate, currentDiTauSolution_, pv);
      solutions.push_back(currentDiTauSolution_);
    }

    gMinuit = 0;
    
    return solutions;
  }
  
  double negLogLikelihood(Int_t numParameters, Double_t* parameter) const
  {
    //std::cout << "<SVfitAlgorithm::negLogLikelihood>:" << std::endl;
    //std::cout << " numParameters = " << numParameters << std::endl;

    if ( !currentDiTau_ ) {
      edm::LogError("SVfitAlgorithm::logLikelihood") 
	<< " Pointer to currentDiTau has not been initialized --> skipping !!";
      return 0.;
    }
    
    readMinuitParameters();
/*
    for ( Int_t iParameter = 0; iParameter <= numParameters; ++iParameter ) {
      minuitParameterValues_[iParameter] = parameters[iParameter];
      //std::cout << " Parameter #" << iParameter << " = " << minuitParameterValues_[iParameter] << std::endl;
    }
 */
    applyMinuitParameters(currentDiTauSolution_);
    
    double negLogLikelihood = 0.;    
    for ( typename std::vector<SVfitDiTauLikelihoodBase<T1,T2>*>::const_iterator likelihoodFunction = likelihoodFunctions_.begin();
	  likelihoodFunction != likelihoodFunctions_.end(); ++likelihoodFunction ) {
      negLogLikelihood += (**likelihoodFunction)(*currentDiTau_, currentDiTauSolution_);
    }
    
    return negLogLikelihood;
  }
  
 private:

  void fitPolarizationHypothesis(const CompositePtrCandidateT1T2MEt<T1,T2>& diTauCandidate,
				 SVfitDiTauSolution& solution,
				 const TransientVertex& pv)
  {
    std::cout << "<SVfitAlgorithm::fitPolarizationHypothesis>:" << std::endl;
    
//--- initialize pointer to current diTau object
    currentDiTau_ = &diTauCandidate;

//--- initialize data-members of diTauSolution object
    if ( pv.isValid() ) {
      currentDiTauSolution_.eventVertexPosition_.SetXYZ(pv.position().x(), pv.position().y(), pv.position().z());
      currentDiTauSolution_.eventVertexPositionErr_ = pv.positionError().matrix();
    }
    currentDiTauSolution_.eventVertexIsValid_ = pv.isValid();
    currentDiTauSolution_.leg1_.p4Vis_ = diTauCandidate.leg1()->p4();
    currentDiTauSolution_.leg2_.p4Vis_ = diTauCandidate.leg2()->p4();
    
//--- initialize start-values of Minuit fit parameters
//
//    CV: how to deal with measurement errors in the visible momenta of the two tau decay "legs"
//        when setting the maximum mass of the neutrino system ?
//
    double pvPositionX, pvPositionXerr, pvPositionY, pvPositionYerr, pvPositionZ, pvPositionZerr;
    if ( pv.isValid() ) {
      pvPositionX = pv.position().x();
      pvPositionXerr = pv.positionError().cxx();
      pvPositionY = pv.position().y();
      pvPositionYerr = pv.positionError().cyy();
      pvPositionZ = pv.position().z();
      pvPositionZerr = pv.positionError().czz();
    } else {
      pvPositionX = 0.;
      pvPositionXerr = 0.1;
      pvPositionY = 0.;
      pvPositionYerr = 0.1;
      pvPositionZ = 0.;
      pvPositionZerr = 10.;
    }

    minuit_.DefineParameter(0 + minuitParameterOffset, "pv_x", pvPositionX, pvPositionXerr,  -1.,  +1.);
    minuit_.DefineParameter(1 + minuitParameterOffset, "pv_y", pvPositionY, pvPositionYerr,  -1.,  +1.);
    minuit_.DefineParameter(2 + minuitParameterOffset, "pv_z", pvPositionZ, pvPositionZerr, -50., +50.);
    minuit_.DefineParameter(3 + minuitParameterOffset, "sv1_thetaRest", 0.5*TMath::Pi(), 0.5*TMath::Pi(), 0., TMath::Pi());
    minuit_.DefineParameter(4 + minuitParameterOffset, "sv1_phiLab", 0., TMath::Pi(), 0., 0.); // do not set limits for phiLab
    double parameter4StartValue, parameter4Error;
    minuit_.GetParameter(4 + minuitParameterOffset, parameter4StartValue, parameter4Error);
    std::cout << " parameter4StartValue = " << parameter4StartValue << " +/- " << parameter4Error << std::endl;
    double leg1Radius0 = diTauCandidate.leg1()->energy()*SVfit_namespace::cTauLifetime/SVfit_namespace::tauLeptonMass;
    minuit_.DefineParameter(5 + minuitParameterOffset, "sv1_radiusLab", leg1Radius0, leg1Radius0, 0., 100.*leg1Radius0); 
    double leg1NuMass0, leg1NuMassErr, leg1NuMassMax;
    if ( !SVfit_namespace::isMasslessNuSystem<T1>() ) {
      leg1NuMass0 = 0.5;
      leg1NuMassErr = 0.5;
      leg1NuMassMax = SVfit_namespace::tauLeptonMass - diTauCandidate.leg1()->mass();
    } else {
      leg1NuMass0 = 0.;
      leg1NuMassErr = 0.;
      leg1NuMassMax = 0.;
    }
    //std::cout << " leg1NuMassMax = " << leg1NuMassMax << std::endl;
    minuit_.DefineParameter(6 + minuitParameterOffset, "sv1_m12", leg1NuMass0, leg1NuMassErr, 0., leg1NuMassMax);
    minuit_.DefineParameter(7 + minuitParameterOffset, "sv2_thetaRest", 0.5*TMath::Pi(), 0.5*TMath::Pi(), 0., TMath::Pi());
    minuit_.DefineParameter(8 + minuitParameterOffset, "sv2_phiLab", 0., TMath::Pi(), 0., 0.); // do not set limits for phiLab
    double leg2Radius0 = diTauCandidate.leg2()->energy()*SVfit_namespace::cTauLifetime/SVfit_namespace::tauLeptonMass;
    minuit_.DefineParameter(9 + minuitParameterOffset, "sv2_radiusLab", leg2Radius0, leg2Radius0, 0., 100.*leg2Radius0); 
    double leg2NuMass0, leg2NuMassErr, leg2NuMassMax;
    if ( !SVfit_namespace::isMasslessNuSystem<T2>() ) {
      leg2NuMass0 = 0.5;
      leg2NuMassErr = 0.5;
      leg2NuMassMax = SVfit_namespace::tauLeptonMass - diTauCandidate.leg2()->mass();
    } else {
      leg2NuMass0 = 0.;
      leg2NuMassErr = 0.;
      leg2NuMassMax = 0.;
    }
    //std::cout << " leg2NuMassMax = " << leg2NuMassMax << std::endl;
    minuit_.DefineParameter(10 + minuitParameterOffset, "sv2_m12", leg2NuMass0, leg2NuMassErr, 0., leg2NuMassMax);

    std::cout << " minuitNumParameters = " << minuit_.GetNumPars()
	      << " (free = " << minuit_.GetNumFreePars() << ", fixed = " << minuit_.GetNumFixedPars() << ")" << std::endl;
    assert(minuit_.GetNumPars() == minuitNumParameters_);
   
    int minuitStatus = minuit_.Command("MIN"); 
    edm::LogInfo("SVfitAlgorithm::fit") 
      << " Minuit fit Status = " << minuitStatus << std::endl;
	
    readMinuitParameters();
    for ( Int_t iParameter = 0; iParameter < minuitNumParameters_; ++iParameter ) {
      std::cout << " Parameter #" << iParameter << " = " << minuitParameterValues_[iParameter] << std::endl;
    }
    applyMinuitParameters(currentDiTauSolution_);
    
    for ( typename std::vector<SVfitDiTauLikelihoodBase<T1,T2>*>::const_iterator likelihoodFunction = likelihoodFunctions_.begin();
	  likelihoodFunction != likelihoodFunctions_.end(); ++likelihoodFunction ) {
      double likelihoodFunctionValue = (**likelihoodFunction)(*currentDiTau_, currentDiTauSolution_);
      currentDiTauSolution_.logLikelihoods_.insert(std::make_pair((*likelihoodFunction)->name(), likelihoodFunctionValue));
    }

    currentDiTauSolution_.minuitStatus_ = minuitStatus;
  }

  void readMinuitParameters() const
  {
    std::cout << "<readMinuitParameters>:" << std::endl;

    Double_t dummy;
    for ( Int_t iParameter = 0; iParameter < minuitNumParameters_; ++iParameter ) {
      minuit_.GetParameter(iParameter + minuitParameterOffset, minuitParameterValues_[iParameter], dummy);
      //std::cout << " Parameter #" << iParameter << " = " << minuitParameterValues_[iParameter] << std::endl;
    }

    if ( TMath::Abs(minuitParameterValues_[4]) > 0.01 ) {
      double dummyValue, dummyError;
      minuit_.GetParameter(4 + minuitParameterOffset, dummyValue, dummyError);
      std::cout << " Parameter #4 = " << minuitParameterValues_[4] << " (" << dummyValue << ")" << std::endl;
    }
  }
  
  void applyMinuitParameters(SVfitDiTauSolution& diTauSolution) const 
  {
//--- set primary event vertex position (tau lepton production vertex)
    diTauSolution.eventVertexPositionCorr_.SetX(minuitParameterValues_[kPrimaryVertexX]);
    diTauSolution.eventVertexPositionCorr_.SetY(minuitParameterValues_[kPrimaryVertexY]);
    diTauSolution.eventVertexPositionCorr_.SetZ(minuitParameterValues_[kPrimaryVertexZ]);

//--- build first tau decay "leg"
    applyMinuitParametersToLeg(kLeg1thetaRest, kLeg1phiLab, kLeg1flightPathLab, kLeg1nuInvMass, diTauSolution.leg1_);

//--- build second tau decay "leg"
    applyMinuitParametersToLeg(kLeg2thetaRest, kLeg2phiLab, kLeg2flightPathLab, kLeg2nuInvMass, diTauSolution.leg2_);
  }

  void applyMinuitParametersToLeg(int gjAngleIndex, int phiLabIndex, int flightDistanceIndex, int massNuNuIndex,
				  SVfitLegSolution& legSolution) const
  {
    double gjAngle        = minuitParameterValues_[gjAngleIndex];
    double phiLab         = minuitParameterValues_[phiLabIndex];
    double flightDistance = minuitParameterValues_[flightDistanceIndex];
    double massNuNu       = minuitParameterValues_[massNuNuIndex];

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

    // Build the tau four vector. By construction, the neutrino is tauP4 - visP4
    legSolution.p4Invis_ =  tauP4 - p4Vis;

    // Build boost vector and compute the rest frame quanitites
    reco::Candidate::Vector boost = tauP4.BoostToCM();
    legSolution.p4VisRestFrame_ = ROOT::Math::VectorUtil::boost(legSolution.p4Vis_, boost);
    legSolution.p4InvisRestFrame_ = ROOT::Math::VectorUtil::boost(legSolution.p4Invis_, boost);

    // Set the flight path
    legSolution.tauFlightPath_ = direction*flightDistance;
  }

  std::string name_;
  
  SVfitEventVertexRefitter* eventVertexRefitAlgorithm_;
  SVfitLegTrackExtractor<T1> leg1TrackExtractor_;
  SVfitLegTrackExtractor<T2> leg2TrackExtractor_;

  std::vector<SVfitDiTauLikelihoodBase<T1,T2>*> likelihoodFunctions_;
  bool likelihoodsSupportPolarization_;
  
  mutable const CompositePtrCandidateT1T2MEt<T1,T2>* currentDiTau_;
  mutable SVfitDiTauSolution currentDiTauSolution_;
  
  mutable TMinuit minuit_;
  const static Int_t minuitNumParameters_ = 11;
  mutable Double_t* minuitParameterValues_;
  bool* minuitLockParameters_;
};

#endif
