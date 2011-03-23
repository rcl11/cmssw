#include "TauAnalysis/CandidateTools/plugins/NSVfitTauToHadLikelihoodPolarization.h"

#include "DataFormats/TauReco/interface/PFTauDecayMode.h"

#include "TauAnalysis/CandidateTools/interface/SVfitVMlineShapeIntegrand.h"
#include "TauAnalysis/CandidateTools/interface/candidateAuxFunctions.h"
#include "TauAnalysis/CandidateTools/interface/NSVfitAlgorithmBase.h"
#include "TauAnalysis/CandidateTools/interface/svFitAuxFunctions.h"

#include "AnalysisDataFormats/TauAnalysis/interface/NSVfitTauToHadHypothesis.h"

#include <TFile.h>
#include <TH2.h>
#include <TMath.h>

#include <algorithm>
#include <limits>

using namespace SVfit_namespace;

NSVfitTauToHadLikelihoodPolarization::decayModeEntryType::decayModeEntryType(const edm::ParameterSet& cfg)
  : xSigma_(0),
    xBias_(0)
{
  if ( cfg.exists("xSigma") ) xSigma_ = new TFormula("xSigma", cfg.getParameter<std::string>("xSigma").data());
  if ( cfg.exists("xBias")  ) xBias_  = new TFormula("xBias",  cfg.getParameter<std::string>("xBias").data());
  pMin_ = cfg.getParameter<double>("pMin");
}

NSVfitTauToHadLikelihoodPolarization::decayModeEntryType::~decayModeEntryType()
{
  delete xSigma_;
  delete xBias_;
}

void NSVfitTauToHadLikelihoodPolarization::decayModeEntryType::print(std::ostream& stream) const
{
  stream << "<decayModeEntryType::print>:" << std::endl;
  if ( xSigma_ != 0 ) std::cout << " xSigma = " << xSigma_->GetTitle() << std::endl;
  if ( xBias_  != 0 ) std::cout << " xBias = " << xBias_->GetTitle() << std::endl;
  std::cout << " pMin = " << pMin_ << std::endl;
}

//
//-------------------------------------------------------------------------------
//

size_t getSupportedTauDecayModeIndex(
         const std::vector<int>& supportedTauDecayModes, const std::string& tauDecayModeName, int tauDecayModeOther)
{
//--- check if tau decay mode name given as function argument
//    matches one of the decay modes for which dedicated likelihood functions 
//    are implemented in the NSVfitTauToHadLikelihoodPolarization class
  size_t numSupportedTauDecayModes = supportedTauDecayModes.size();
  for ( size_t iDecayMode = 0; iDecayMode < numSupportedTauDecayModes; ++iDecayMode ) {
    if ( getTauDecayModeName(supportedTauDecayModes[iDecayMode]) == tauDecayModeName ) return iDecayMode;
  }
  
//--- tau decay mode name given as function argument is
//    not within the list of supported decay modes
  return tauDecayModeOther;
}

void normalizeMatrixColumns(TMatrixD& matrix)
{
  unsigned numRows = matrix.GetNrows();
  unsigned numColumns = matrix.GetNcols();
  
  for ( unsigned iColumn = 0; iColumn < numColumns; ++iColumn ) {
    double sum = 0.;
    
    for ( unsigned iRow = 0; iRow < numRows; ++iRow ) {
      sum += matrix(iRow, iColumn);
    }
    
    if ( sum > 0. ) {
      for ( unsigned iRow = 0; iRow < numRows; ++iRow ) {
	matrix(iRow, iColumn) /= sum;
      }
    } else {
      edm::LogError("normalizeMatrixRows") 
	<< " Sum of elements = " << sum << " for column = " << iColumn << " --> matrix will **not** be normalized !!";
    }
  }
}

NSVfitTauToHadLikelihoodPolarization::NSVfitTauToHadLikelihoodPolarization(const edm::ParameterSet& cfg)
  : NSVfitSingleParticleLikelihood(cfg),
    rhoLpolLineShape_(0),
    rhoTpolLineShape_(0),
    a1LpolLineShape_(0),
    a1TpolLineShape_(0),
    likelihoodPhaseSpace_(0)
{
  //std::cout << "<NSVfitTauToHadLikelihoodPolarization::NSVfitTauToHadLikelihoodPolarization>:" << std::endl;

//--- create list of supported decay modes
//   o tau- --> pi- nu
//   o tau- --> rho- nu --> pi- pi0 nu
//   o tau- --> a1- nu --> pi- pi0 pi0 nu or pi- pi+ pi- nu
  supportedTauDecayModes_.resize(kOther + 1);
  supportedTauDecayModes_[kPion]        = reco::PFTauDecayMode::tauDecay1ChargedPion0PiZero;
  supportedTauDecayModes_[kVMrho]       = reco::PFTauDecayMode::tauDecay1ChargedPion1PiZero;
  supportedTauDecayModes_[kVMa1Neutral] = reco::PFTauDecayMode::tauDecay1ChargedPion2PiZero;
  supportedTauDecayModes_[kVMa1Charged] = reco::PFTauDecayMode::tauDecay3ChargedPion0PiZero;
  supportedTauDecayModes_[kOther]       = reco::PFTauDecayMode::tauDecayOther;
  numSupportedTauDecayModes_ = supportedTauDecayModes_.size();

  vRec_.ResizeTo(numSupportedTauDecayModes_);
  mapRecToGenTauDecayModes_.ResizeTo(numSupportedTauDecayModes_, numSupportedTauDecayModes_);
  mapRecToGenTauDecayModes_.Zero();
  vGen_.ResizeTo(numSupportedTauDecayModes_);
  vProb_.ResizeTo(numSupportedTauDecayModes_);

//--- create "transfer matrix"
//
//    CV: In the histogram, the generated (reconstructed) tau decay modes are on the x-axis (y-axis);
//        the x-axis (y-axis) needs to be mapped to rows (columns) of the "transfer matrix".
//
//        The aim of the "transfer matrix" M is to map from reconstructed to generated ("true") tau lepton hadronic decay modes.
//        The mapping is implemented by matrix multiplication:
//         vGen = M * vRec,
//        where vRec is a vector encoding the reconstructed hadronic tau decay mode
//        (exactly one entry in this vector is 1, all other entries are 0)
//        and vGen gives the probabilities for the "true" decay mode of the tau lepton.
//
//        The required format of M is:
//         | p(gen=1Prong0Pi0|rec=1Prong0Pi0) .. p(gen=1Prong0Pi0|rec=3Prong0Pi0) p(gen=1Prong0Pi0|rec=other) |
//         | p(gen=1Prong1Pi0|rec=1Prong0Pi0) .. p(gen=1Prong1Pi0|rec=3Prong0Pi0) p(gen=1Prong1Pi0|rec=other) |
//         | p(gen=1Prong2Pi0|rec=1Prong0Pi0) .. p(gen=1Prong2Pi0|rec=3Prong0Pi0) p(gen=1Prong2Pi0|rec=other) |
//         | p(gen=3Prong0Pi0|rec=1Prong0Pi0) .. p(gen=3Prong0Pi0|rec=3Prong0Pi0) p(gen=3Prong0Pi0|rec=other) |
//         | p(gen=other     |rec=1Prong0Pi0) .. p(gen=other     |rec=3Prong0Pi0) p(gen=other     |rec=other) |
//
//        Note that **column** of the "transfer matrix" need to be normalized to unit probability.
//

//--- load histogram correlating reconstructed to generated hadronic decay modes
//    if it exists, else take "transfer matrix" to be diagonal
  if ( cfg.exists("mapRecToGenTauDecayModes") ) {
    edm::ParameterSet cfgMapRecToGenTauDecayModes = cfg.getParameter<edm::ParameterSet>("mapRecToGenTauDecayModes");
    
    std::string fileName_mapRecToGenTauDecayModes = cfgMapRecToGenTauDecayModes.getParameter<std::string>("fileName");
    std::string meName_mapRecToGenTauDecayModes = cfgMapRecToGenTauDecayModes.getParameter<std::string>("meName");
    
    TFile* file_mapRecToGenTauDecayModes = 0;
    try {
      file_mapRecToGenTauDecayModes = TFile::Open(fileName_mapRecToGenTauDecayModes.data());
    } catch (...) { 
      edm::LogError("NSVfitTauToHadLikelihoodPolarization") 
	<< " Failed to open file = " << fileName_mapRecToGenTauDecayModes << " !!";
    }
    
    
    TH2* histogram_mapRecToGenTauDecayModes 
      = dynamic_cast<TH2*>(file_mapRecToGenTauDecayModes->Get(meName_mapRecToGenTauDecayModes.data()));
    if ( histogram_mapRecToGenTauDecayModes ) {
      int numBinsX = histogram_mapRecToGenTauDecayModes->GetNbinsX();
      for ( int binIndex_x = 1; binIndex_x <= numBinsX; ++binIndex_x ) {
	std::string binLabel_x = histogram_mapRecToGenTauDecayModes->GetXaxis()->GetBinLabel(binIndex_x);
	
	size_t tauDecayMode_row = getSupportedTauDecayModeIndex(supportedTauDecayModes_, binLabel_x, kOther);
	assert(tauDecayMode_row >= 0 && tauDecayMode_row < numSupportedTauDecayModes_);
	
	int numBinsY = histogram_mapRecToGenTauDecayModes->GetNbinsY();
	for ( int binIndex_y = 1; binIndex_y <= numBinsY; ++binIndex_y ) {
	  std::string binLabel_y = histogram_mapRecToGenTauDecayModes->GetYaxis()->GetBinLabel(binIndex_y);
	  
	  size_t tauDecayMode_column = getSupportedTauDecayModeIndex(supportedTauDecayModes_, binLabel_y, kOther);
	  assert(tauDecayMode_column >= 0 && tauDecayMode_column < numSupportedTauDecayModes_);
	  
	  mapRecToGenTauDecayModes_(tauDecayMode_row, tauDecayMode_column) 
	    += histogram_mapRecToGenTauDecayModes->GetBinContent(binIndex_x, binIndex_y);
	}
      }
    }

    normalizeMatrixColumns(mapRecToGenTauDecayModes_);

//--- close file from which histogram was loaded
    delete file_mapRecToGenTauDecayModes;
  } else {
    for ( size_t iDecayMode = 0; iDecayMode < numSupportedTauDecayModes_; ++iDecayMode ) {
      mapRecToGenTauDecayModes_(iDecayMode, iDecayMode) = 1.;
    }
  }

  //std::cout << " mapRecToGenTauDecayModes:" << std::endl;
  //mapRecToGenTauDecayModes_.Print();

//--- initialize resolution and bias values for reconstructing
//    momentum fraction x = E(dist. pion)/E(vector meson)
//    of rho/a1 vector meson carried by "distinguishable" pion
//   (computed in laboratory frame)
  edm::ParameterSet cfgDecayModes = cfg.getParameter<edm::ParameterSet>("decayModeParameters");
  decayModeParameters_.resize(numSupportedTauDecayModes_);
  decayModeParameters_[kPion]        = new decayModeEntryType(cfgDecayModes.getParameter<edm::ParameterSet>("oneProngZeroPi0s")); 
  decayModeParameters_[kVMrho]       = new decayModeEntryType(cfgDecayModes.getParameter<edm::ParameterSet>("oneProngOnePi0"));
  decayModeParameters_[kVMa1Neutral] = new decayModeEntryType(cfgDecayModes.getParameter<edm::ParameterSet>("oneProngTwoPi0s"));
  decayModeParameters_[kVMa1Charged] = new decayModeEntryType(cfgDecayModes.getParameter<edm::ParameterSet>("threeProngZeroPi0s"));

//--- create auxiliary classes for computation of vector meson line-shape integrals
  useCollApproxFormulas_ = cfg.exists("useCollApproxFormulas") ?
    cfg.getParameter<bool>("useCollApproxFormulas") : false;

  rhoLpolLineShape_ = new SVfitVMlineShapeIntegral(SVfitVMlineShapeIntegrand::kVMrho, 
						   SVfitVMlineShapeIntegrand::kVMlongitudinalPol, useCollApproxFormulas_);
  rhoTpolLineShape_ = new SVfitVMlineShapeIntegral(SVfitVMlineShapeIntegrand::kVMrho, 
						   SVfitVMlineShapeIntegrand::kVMtransversePol, useCollApproxFormulas_);
  a1LpolLineShape_  = new SVfitVMlineShapeIntegral(SVfitVMlineShapeIntegrand::kVMa1, 
						   SVfitVMlineShapeIntegrand::kVMlongitudinalPol, useCollApproxFormulas_);
  a1TpolLineShape_  = new SVfitVMlineShapeIntegral(SVfitVMlineShapeIntegrand::kVMa1, 
						   SVfitVMlineShapeIntegrand::kVMtransversePol, useCollApproxFormulas_);

//--- generic "phase-space" plugin to be used for computing likelihood 
//    for hadronic tau decays not in the list of decay modes 
//    supported by the NSVfitTauToHadLikelihoodPolarization plugin
  edm::ParameterSet cfgLikelihoodPhaseSpace;
  std::string pluginTypeLikelihoodPhaseSpace = "NSVfitTauToHadLikelihoodPhaseSpace";
  cfgLikelihoodPhaseSpace.addParameter<std::string>("pluginType", pluginTypeLikelihoodPhaseSpace);
  cfgLikelihoodPhaseSpace.addParameter<std::string>("pluginName", "nSVfitTauToHadLikelihoodPhaseSpace");
  cfgLikelihoodPhaseSpace.addParameter<std::string>("prodParticleLabel", prodParticleLabel_);
  cfgLikelihoodPhaseSpace.addParameter<int>("verbosity", verbosity_);
  likelihoodPhaseSpace_ = 
    NSVfitSingleParticleLikelihoodPluginFactory::get()->create(pluginTypeLikelihoodPhaseSpace, cfgLikelihoodPhaseSpace);

//--- initialize temporary variables,
//    defined to speed-up computations
  a1posMassTerm_ = (a1MesonMass2 + rhoMesonMass2)/(a1MesonMass*rhoMesonMass);
  a1posMassTerm2_ = square(a1posMassTerm_);
  a1negMassTerm_ = (a1MesonMass2 - rhoMesonMass2)/(a1MesonMass*rhoMesonMass);
  a1q_ = 0.5*TMath::Sqrt(a1MesonMass2 - 4.*chargedPionMass2);
  a1_8piDiv9_ = 8.*TMath::Pi()/9.;
  a1_16piDiv9_ = 2.*a1_8piDiv9_;
}

NSVfitTauToHadLikelihoodPolarization::~NSVfitTauToHadLikelihoodPolarization()
{
  for ( std::vector<decayModeEntryType*>::iterator it = decayModeParameters_.begin();
	it != decayModeParameters_.end(); ++it ) {
    delete (*it);
  }

  delete rhoLpolLineShape_;
  delete rhoTpolLineShape_;
  delete a1LpolLineShape_;
  delete a1TpolLineShape_;

  delete likelihoodPhaseSpace_;
}

void NSVfitTauToHadLikelihoodPolarization::beginJob(NSVfitAlgorithmBase* algorithm)
{
  algorithm->requestFitParameter(prodParticleLabel_, kTau_visEnFracX,  pluginName_);
  algorithm->requestFitParameter(prodParticleLabel_, kTau_phi_lab,     pluginName_);
  algorithm->requestFitParameter(prodParticleLabel_, kTau_pol,         pluginName_);
  algorithm->requestFitParameter(prodParticleLabel_, kTauVM_theta_rho, pluginName_);
  algorithm->requestFitParameter(prodParticleLabel_, kTauVM_theta_a1,  pluginName_);
  algorithm->requestFitParameter(prodParticleLabel_, kTauVM_theta_a1r, pluginName_);
  algorithm->requestFitParameter(prodParticleLabel_, kTauVM_phi_a1r,   pluginName_);
}

void NSVfitTauToHadLikelihoodPolarization::beginCandidate(const NSVfitSingleParticleHypothesisBase* hypothesis)
{
  const NSVfitTauToHadHypothesis* hypothesis_T = dynamic_cast<const NSVfitTauToHadHypothesis*>(hypothesis);
  assert(hypothesis_T != 0);
  
  if ( this->verbosity_ ) std::cout << "<NSVfitTauToHadLikelihoodPolarization::beginCandidate>:" << std::endl;

  int recTauDecayMode = hypothesis_T->decayMode();
  //std::cout << " recTauDecayMode = " << recTauDecayMode << std::endl;
  
  vRec_.Zero();
  std::vector<int>::const_iterator tauDecayMode_index 
    = std::find(supportedTauDecayModes_.begin(), supportedTauDecayModes_.end(), recTauDecayMode);
  if ( tauDecayMode_index != supportedTauDecayModes_.end() ) {
    vRec_(tauDecayMode_index - supportedTauDecayModes_.begin()) = 1.;
  } else {
    vRec_(kOther) = 1.;
  }
  
  //std::cout << " vRec:" << std::endl;
  //vRec_.Print();

  vGen_ = vRec_;
  vGen_ *= mapRecToGenTauDecayModes_;

  //std::cout << " vGen:" << std::endl;
  //vGen_.Print();
}

double NSVfitTauToHadLikelihoodPolarization::operator()(const NSVfitSingleParticleHypothesisBase* hypothesis) const
{
//--- compute negative log-likelihood for tau lepton decay "leg"
//    to be compatible with decay tau- --> X nu of polarized tau lepton into hadrons,
//    assuming  matrix element of V-A electroweak decay
//
//    NOTE: The formulas taken from the papers
//         [1] "Tau polarization and its correlations as a probe of new physics",
//             B.K. Bullock, K. Hagiwara and A.D. Martin,
//             Nucl. Phys. B395 (1993) 499.
//         [2] "Charged Higgs boson search at the TeVatron upgrade using tau polarization",
//             S. Raychaudhuri and D.P. Roy,
//             Phys. Rev.  D52 (1995) 1556.           
//
  const NSVfitTauToHadHypothesis* hypothesis_T = dynamic_cast<const NSVfitTauToHadHypothesis*>(hypothesis);
  assert(hypothesis_T != 0);

  if ( this->verbosity_ ) std::cout << "<NSVfitTauToHadLikelihoodPolarization::negLogLikelihoodPolarized>:" << std::endl;

  double normProb = 0.;
  for ( size_t iDecayMode = 0; iDecayMode < numSupportedTauDecayModes_; ++iDecayMode ) {
    double probDecayMode = 0.;
    if      ( supportedTauDecayModes_[iDecayMode] == reco::PFTauDecayMode::tauDecay1ChargedPion0PiZero ) 
      probDecayMode = probOneProngZeroPi0s(hypothesis_T);
    else if ( supportedTauDecayModes_[iDecayMode] == reco::PFTauDecayMode::tauDecay1ChargedPion1PiZero )
      probDecayMode = probOneProngOnePi0(hypothesis_T);
    else if ( supportedTauDecayModes_[iDecayMode] == reco::PFTauDecayMode::tauDecay1ChargedPion2PiZero )
      probDecayMode = probOneProngTwoPi0s(hypothesis_T);
    else if ( supportedTauDecayModes_[iDecayMode] == reco::PFTauDecayMode::tauDecay3ChargedPion0PiZero )
      probDecayMode = probThreeProngZeroPi0s(hypothesis_T);
    else 
      probDecayMode = probOtherDecayMode(hypothesis_T);
    vProb_(iDecayMode) = probDecayMode;
    normProb += vGen_(iDecayMode);
  }

  if ( this->verbosity_ ) {
    std::cout << " vProb:" << std::endl;
    vProb_.Print();
  }

  double prob = (vGen_*vProb_)/normProb;

  double nll = 0.;
  if ( prob > 0. ) {
    nll = -TMath::Log(prob);
  } else {
    if ( prob < 0. ) 
      edm::LogWarning ("NSVfitTauToHadLikelihoodPolarization::operator()")
	<< " Unphysical solution: prob = " << prob << " --> returning very large negative number !!";
    nll = std::numeric_limits<float>::max();
  }

  if ( this->verbosity_ ) std::cout << "--> nll = " << nll << std::endl;

  return nll;
}

//
//-------------------------------------------------------------------------------
//

double NSVfitTauToHadLikelihoodPolarization::probOneProngZeroPi0s(const NSVfitTauToHadHypothesis* hypothesis_T) const
{
  if ( this->verbosity_ ) std::cout << "<NSVfitTauToHadLikelihoodPolarization::probOneProngZeroPi0s>:" << std::endl;
          
  double prob = 0.;
  if ( !useCollApproxFormulas_ ) {
    double theta = hypothesis_T->decay_angle_rf();
    double cosTheta = TMath::Cos(theta);
    double sinTheta = TMath::Sin(theta);
    double tauLeptonPol = hypothesis_T->polarization();
    prob = 0.5*(1. + tauLeptonPol*cosTheta)*sinTheta; // [1], formula (2.1)
  } else {    
    double z = hypothesis_T->visEnFracX();            // tau lepton visible momentum fraction
    double tauLeptonPol = hypothesis_T->polarization();
    prob = (1. + tauLeptonPol*(2*z - 1.));            // [1], formula (2.4)
  }

//--- multiply tau- --> pi- nu decay probability by
//    average over angles { thetaVMa1, thetaVMa1r, phiVMa1r }
//   (of tau- --> a1- nu --> pi- pi0 pi0 nu, tau- --> a1- nu --> pi- pi+ pi- nu decay),
//    in order to compare probabilities fo different hadronic decay modes on "equal footing"
//   (the computed prob values are in fact probability densities, the integral of which is normalized to one;
//    as a consequence, decay modes with more prob factors would get "penalized" otherwise)
  double decayModeCorrFactor = 1./(2*cube(TMath::Pi()));

  return decayModeCorrFactor*prob;
}

double compProbSmear(double xResidual, double xSigma)
{
//--- add protection against sigma equals zero case
  const double epsilon = 1.e-6;
  if ( xSigma < epsilon ) xSigma = epsilon;

//--- add protection against case of very large residuals (in terms of sigma),
//    in order to avoid problems with numerical stability 
//   (increase sigma by sqrt(num. standard deviations),
//    so that protection has no effect once the fit converges to small residuals)
  double numSigma = TMath::Abs(xResidual/xSigma);
  const double numSigma_cutoff = 2.;
  if ( numSigma > numSigma_cutoff ) xSigma *= TMath::Sqrt(numSigma - numSigma_cutoff + 1.);

  return TMath::Gaus(xResidual, 0., xSigma);
}

double NSVfitTauToHadLikelihoodPolarization::probOneProngOnePi0(const NSVfitTauToHadHypothesis* hypothesis_T) const
{
//--- compute likelihood for tau- --> rho- nu --> pi- pi0 nu decay

  if ( this->verbosity_ ) std::cout << "<NSVfitTauToHadLikelihoodPolarization::probOneProngOnePi0>:" << std::endl;

  const pat::Tau* tauPtr = dynamic_cast<const pat::Tau*>(hypothesis_T->particle().get());
  assert(tauPtr);

  double theta = hypothesis_T->decay_angle_rf();
  double z = hypothesis_T->visEnFracX();
  double tauLeptonPol = hypothesis_T->polarization();
  
  double probTauDecayL = (*rhoLpolLineShape_)(theta, tauLeptonPol, z);
  double probTauDecayT = (*rhoTpolLineShape_)(theta, tauLeptonPol, z);
  //std::cout << " probTauDecayL = " << probTauDecayL << ", probTauDecayT = " << probTauDecayT << std::endl;

//--- find "distinguishable" pion in tau-jet;
//    in case "distinguishable" pion cannot be found 
//   (e.g. in case "wrong" tau decay mode is reconstructed),
//    assume the "distinguishable" pion to be very soft
  const reco::Candidate* distPion = getDistPion(*tauPtr);
  double xMeasured = ( distPion != 0 ) ? (distPion->energy()/tauPtr->energy()) : 0.;
  //std::cout << " xMeasured = " << xMeasured << std::endl;
  double thetaVMrho = hypothesis_T->decay_angle_VMrho();
  double cosThetaVMrho = TMath::Cos(thetaVMrho);
  double sinThetaVMrho = TMath::Sin(thetaVMrho);
  double xFitted = 0.5*(1. + TMath::Sqrt(1 - 4.*(chargedPionMass2/rhoMesonMass2))*cosThetaVMrho); // [2], formula (41)
  //std::cout << " xFitted = " << xFitted << std::endl;

  double probVMrhoDecayL, probVMrhoDecayT;
  if ( !useCollApproxFormulas_ ) {
    probVMrhoDecayL = 1.5*square(cosThetaVMrho)*sinThetaVMrho; // [2], formula (39)
    probVMrhoDecayT = 0.75*cube(sinThetaVMrho);                // [2], formula (40)
  } else {
    probVMrhoDecayL = 1.5*(2*xFitted - 1.);                    // [2], formula (39)
    probVMrhoDecayT = 3*xFitted*(1. - xFitted) ;               // [2], formula (40)
  }

  //std::cout << " probVMrhoDecayL = " << probVMrhoDecayL << ", probVMrhoDecayT = " << probVMrhoDecayT << std::endl;

  double xSigma = decayModeParameters_[kVMrho]->xSigma_->Eval(tauPtr->pt());
  double xBias = decayModeParameters_[kVMrho]->xBias_->Eval(tauPtr->pt());
  //std::cout << " xSigma = " << xSigma << ", xBias = " << xBias << std::endl;

  double xResidual = xMeasured - xFitted - xBias;
  //std::cout << " xResidual = " << xResidual << std::endl;

  double probSmear = compProbSmear(xResidual, xSigma);
  //std::cout << " probSmear = " << probSmear << std::endl;

  double probL = probTauDecayL*probVMrhoDecayL*probSmear;
  double probT = probTauDecayT*probVMrhoDecayT*probSmear;
  //std::cout << "--> probL = " << probL << ", probT = " << probT << std::endl;

//--- multiply tau- --> rho- nu --> pi- pi0 nu decay probability by
//    average over angles { thetaVMa1r, phiVMa1r } 
//   (of tau- --> a1- nu --> pi- pi0 pi0 nu, tau- --> a1- nu --> pi- pi+ pi- nu decay),
//    in order to compare probabilities fo different hadronic decay modes on "equal footing"
//   (the computed prob values are in fact probability densities, the integral of which is normalized to one;
//    as a consequence, decay modes with more prob factors would get "penalized" otherwise)
  double decayModeCorrFactor = 1./(2*square(TMath::Pi()));

  return decayModeCorrFactor*(probL + probT);
}

double NSVfitTauToHadLikelihoodPolarization::compVMa1x(
         double cosThetaVMa1, double sinThetaVMa1, double cosThetaVMa1r, double sinThetaVMa1r, double cosPhiVMa1r) const
{
  return (1./a1MesonMass)*(0.25*rhoMesonMass*(a1posMassTerm_ + a1negMassTerm_*cosThetaVMa1) 
                          + 0.5*a1q_*(a1posMassTerm_*cosThetaVMa1*cosThetaVMa1r + a1negMassTerm_*cosThetaVMa1r 
				     - 2.*sinThetaVMa1*sinThetaVMa1r*cosPhiVMa1r)); // [2], formula (48)
}

double NSVfitTauToHadLikelihoodPolarization::compVMa1DecayProbL(
	 double cosThetaVMa1, double sinThetaVMa1, double cosThetaVMa1r, double sinThetaVMa1r, double cosPhiVMa1r) const
{
  return (square(a1posMassTerm_*cosThetaVMa1*cosThetaVMa1r - 2.*sinThetaVMa1*sinThetaVMa1r*cosPhiVMa1r)/
	  (a1_8piDiv9_*(a1posMassTerm2_ + 8.)))*sinThetaVMa1*sinThetaVMa1r; // [2], formula (46)
}

double NSVfitTauToHadLikelihoodPolarization::compVMa1DecayProbT(
	 double cosThetaVMa1, double sinThetaVMa1, 
	 double cosThetaVMa1r, double sinThetaVMa1r, double cosPhiVMa1r, double sinPhiVMa1r) const
{
  return ((square(a1posMassTerm_*sinThetaVMa1*cosThetaVMa1r + 2.*cosThetaVMa1*sinThetaVMa1r*cosPhiVMa1r)
          + 4.*square(sinThetaVMa1r)*square(sinPhiVMa1r))/
	  (a1_16piDiv9_*(a1posMassTerm2_ + 8.)))*sinThetaVMa1*sinThetaVMa1r; // [2], formula (47)
}

double NSVfitTauToHadLikelihoodPolarization::probOneProngTwoPi0s(const NSVfitTauToHadHypothesis* hypothesis_T) const
{
//--- compute likelihood for tau- --> a1- nu --> pi- pi0 pi0 nu decay

  if ( this->verbosity_ ) std::cout << "<NSVfitTauToHadLikelihoodPolarization::probOneProngTwoPi0s>:" << std::endl;

  const pat::Tau* tauPtr = dynamic_cast<const pat::Tau*>(hypothesis_T->particle().get());
  assert(tauPtr);

  double theta = hypothesis_T->decay_angle_rf();
  double z = hypothesis_T->visEnFracX();
  double tauLeptonPol = hypothesis_T->polarization();
  
  double probTauDecayL = (*a1LpolLineShape_)(theta, tauLeptonPol, z);
  double probTauDecayT = (*a1TpolLineShape_)(theta, tauLeptonPol, z);
  //std::cout << " probTauDecayL = " << probTauDecayL << ", probTauDecayT = " << probTauDecayT << std::endl;

//--- find "distinguishable" pion in tau-jet;
//    in case "distinguishable" pion cannot be found 
//   (e.g. in case "wrong" tau decay mode is reconstructed),
//    assume the "distinguishable" pion to be very soft
  const reco::Candidate* distPion = getDistPion(*tauPtr);
  double xMeasured = ( distPion != 0 ) ? (distPion->energy()/tauPtr->energy()) : 0.;
  //std::cout << " xMeasured = " << xMeasured << std::endl;
  double thetaVMa1 = hypothesis_T->decay_angle_VMa1();
  double cosThetaVMa1  = TMath::Cos(thetaVMa1);
  double sinThetaVMa1  = TMath::Sin(thetaVMa1);
  double thetaVMa1r = hypothesis_T->decay_angle_VMa1r_theta();
  double cosThetaVMa1r = TMath::Cos(thetaVMa1r);
  double sinThetaVMa1r = TMath::Sin(thetaVMa1r);
  double phiVMa1r = hypothesis_T->decay_angle_VMa1r_phi();
  double cosPhiVMa1r   = TMath::Cos(phiVMa1r);
  double sinPhiVMa1r   = TMath::Sin(phiVMa1r);
  double xFitted = compVMa1x(cosThetaVMa1, sinThetaVMa1, cosThetaVMa1r, sinThetaVMa1r, cosPhiVMa1r);
  //std::cout << " xFitted = " << xFitted << std::endl;

//--- CV: only non-collinear approximation type formulas available in literature
//        for tau- --> a1- nu --> pi- pi0 pi0 nu decay
  double probVMa1DecayL = compVMa1DecayProbL(cosThetaVMa1, sinThetaVMa1, cosThetaVMa1r, sinThetaVMa1r, cosPhiVMa1r);
  double probVMa1DecayT = compVMa1DecayProbT(cosThetaVMa1, sinThetaVMa1, cosThetaVMa1r, sinThetaVMa1r, cosPhiVMa1r, sinPhiVMa1r);
  //std::cout << " probVMa1DecayL = " << probVMa1DecayL << ", probVMa1DecayT = " << probVMa1DecayT << std::endl;

  double xSigma = decayModeParameters_[kVMa1Neutral]->xSigma_->Eval(tauPtr->pt());
  double xBias = decayModeParameters_[kVMa1Neutral]->xBias_->Eval(tauPtr->pt());
  //std::cout << " xSigma = " << xSigma << ", xBias = " << xBias << std::endl;

  double xResidual = xMeasured - xFitted - xBias;
  //std::cout << " xResidual = " << xResidual << std::endl;

  double probSmear = compProbSmear(xResidual, xSigma);
  //std::cout << " probSmear = " << probSmear << std::endl;

  double probL = probTauDecayL*probVMa1DecayL*probSmear;
  double probT = probTauDecayT*probVMa1DecayT*probSmear;
  //std::cout << "--> probL = " << probL << ", probT = " << probT << std::endl;

  return (probL + probT);
}

double NSVfitTauToHadLikelihoodPolarization::probThreeProngZeroPi0s(const NSVfitTauToHadHypothesis* hypothesis_T) const
{
//--- compute likelihood for tau- --> a1- nu --> pi- pi+ pi- nu decay

  if ( this->verbosity_ ) std::cout << "<NSVfitTauToHadLikelihoodPolarization::probThreeProngZeroPi0s>:" << std::endl;

  const pat::Tau* tauPtr = dynamic_cast<const pat::Tau*>(hypothesis_T->particle().get());
  assert(tauPtr);

  double theta = hypothesis_T->decay_angle_rf();
  double z = hypothesis_T->visEnFracX();
  double tauLeptonPol = hypothesis_T->polarization();
  
  double probTauDecayL = (*a1LpolLineShape_)(theta, tauLeptonPol, z);
  double probTauDecayT = (*a1TpolLineShape_)(theta, tauLeptonPol, z);
  //std::cout << " probTauDecayL = " << probTauDecayL << ", probTauDecayT = " << probTauDecayT << std::endl;

//--- find "distinguishable" pion in tau-jet;
//    in case "distinguishable" pion cannot be found 
//   (e.g. in case "wrong" tau decay mode is reconstructed),
//    assume the "distinguishable" pion to be very soft
  const reco::Candidate* distPion = getDistPion(*tauPtr);
  double xMeasured = ( distPion != 0 ) ? (distPion->energy()/tauPtr->energy()) : 0.;
  //std::cout << " xMeasured = " << xMeasured << std::endl;
  double thetaVMa1 = hypothesis_T->decay_angle_VMa1();
  double cosThetaVMa1  = TMath::Cos(thetaVMa1);
  double sinThetaVMa1  = TMath::Sin(thetaVMa1);
  double thetaVMa1r = hypothesis_T->decay_angle_VMa1r_theta();
  double cosThetaVMa1r = TMath::Cos(thetaVMa1r);
  double sinThetaVMa1r = TMath::Sin(thetaVMa1r);
  double phiVMa1r = hypothesis_T->decay_angle_VMa1r_phi();
  double cosPhiVMa1r   = TMath::Cos(phiVMa1r);
  double sinPhiVMa1r   = TMath::Sin(phiVMa1r);
  double xFitted = compVMa1x(cosThetaVMa1, sinThetaVMa1, cosThetaVMa1r, sinThetaVMa1r, cosPhiVMa1r);
  //std::cout << " xFitted = " << xFitted << std::endl;

//--- CV: only non-collinear approximation type formulas available in literature
//        for tau- --> a1- nu --> pi- pi+ pi- nu decay
  double probVMa1DecayL = compVMa1DecayProbL(cosThetaVMa1, sinThetaVMa1, cosThetaVMa1r, sinThetaVMa1r, cosPhiVMa1r);
  double probVMa1DecayT = compVMa1DecayProbT(cosThetaVMa1, sinThetaVMa1, cosThetaVMa1r, sinThetaVMa1r, cosPhiVMa1r, sinPhiVMa1r);
  //std::cout << " probVMa1DecayL = " << probVMa1DecayL << ", probVMa1DecayT = " << probVMa1DecayT << std::endl;

  double xSigma = decayModeParameters_[kVMa1Charged]->xSigma_->Eval(tauPtr->pt());
  double xBias = decayModeParameters_[kVMa1Charged]->xBias_->Eval(tauPtr->pt());
  //std::cout << " xSigma = " << xSigma << ", xBias = " << xBias << std::endl;

  double xResidual = xMeasured - xFitted - xBias;
  //std::cout << " xResidual = " << xResidual << std::endl;

  double probSmear = compProbSmear(xResidual, xSigma);
  //std::cout << " probSmear = " << probSmear << std::endl;

  double probL = probTauDecayL*probVMa1DecayL*probSmear;
  double probT = probTauDecayT*probVMa1DecayT*probSmear;
  //std::cout << "--> probL = " << probL << ", probT = " << probT << std::endl;

  return (probL + probT);
}

double NSVfitTauToHadLikelihoodPolarization::probOtherDecayMode(const NSVfitTauToHadHypothesis* hypothesis_T) const
{
  return (*likelihoodPhaseSpace_)(hypothesis_T);
}

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_EDM_PLUGIN(NSVfitSingleParticleLikelihoodPluginFactory, NSVfitTauToHadLikelihoodPolarization, "NSVfitTauToHadLikelihoodPolarization");
