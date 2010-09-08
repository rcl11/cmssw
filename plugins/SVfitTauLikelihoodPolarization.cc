#include "TauAnalysis/CandidateTools/plugins/SVfitTauLikelihoodPolarization.h"

#include "DataFormats/TauReco/interface/PFTauDecayMode.h"

#include "TauAnalysis/CandidateTools/interface/SVfitVMlineShapeIntegrand.h"
#include "TauAnalysis/CandidateTools/interface/candidateAuxFunctions.h"
#include "TauAnalysis/CandidateTools/interface/svFitAuxFunctions.h"

#include <TFile.h>
#include <TH2.h>
#include <TMath.h>

#include <algorithm>
#include <limits>

using namespace SVfit_namespace;

size_t getSupportedTauDecayModeIndex(const std::vector<int>& supportedTauDecayModes, const std::string& tauDecayModeName)
{
  size_t numSupportedTauDecayModes = supportedTauDecayModes.size();
  for ( size_t iDecayMode = 0; iDecayMode < numSupportedTauDecayModes; ++iDecayMode ) {
    if ( getTauDecayModeName(iDecayMode) == tauDecayModeName ) return iDecayMode;
  }
  
  edm::LogError("getSupportedTauDecayModeIndex") 
    << " Invalid tau decay mode = " << tauDecayModeName << " !!";
  return -1;
}

void normalizeMatrixRows(TMatrixD& matrix)
{
  unsigned numRows = matrix.GetNrows();
  unsigned numColumns = matrix.GetNcols();
  
  for ( unsigned iRow = 0; iRow < numRows; ++iRow ) {
    double rowSum = 0.;
    
    for ( unsigned iColumn = 0; iColumn < numColumns; ++iColumn ) {
      rowSum += matrix(iRow, iColumn);
    }
    
    if ( rowSum > 0. ) {
      for ( unsigned iColumn = 0; iColumn < numColumns; ++iColumn ) {
	matrix(iRow, iColumn) /= rowSum;
      }
    } else {
      edm::LogError("normalizeMatrixRows") 
	<< " Sum of elements = " << rowSum << " for row = " << iRow << " --> matrix will **not** be normalized !!";
    }
  }

  matrix.Print();
}

SVfitTauLikelihoodPolarization::SVfitTauLikelihoodPolarization(const edm::ParameterSet& cfg)
  : SVfitLegLikelihoodPolarizationBase<pat::Tau>(cfg),
    xVMrhoSigma_(0),
    xVMrhoBias_(0),
    xVMa1NeutralSigma_(0),
    xVMa1NeutralBias_(0),
    xVMa1ChargedSigma_(0),
    xVMa1ChargedBias_(0),
    rhoLpolLineShape_(0),
    rhoTpolLineShape_(0),
    a1LpolLineShape_(0),
    a1TpolLineShape_(0),
    likelihoodPhaseSpace_(0)
{
//--- create list of supported decay modes
//   o tau- --> pi- nu
//   o tau- --> rho- nu --> pi- pi0 nu
//   o tau- --> a1- nu --> pi- pi0 pi0 nu or pi- pi+ pi- nu
  supportedTauDecayModes_.resize(5);
  supportedTauDecayModes_[0] = reco::PFTauDecayMode::tauDecay1ChargedPion0PiZero;
  supportedTauDecayModes_[1] = reco::PFTauDecayMode::tauDecay1ChargedPion1PiZero;
  supportedTauDecayModes_[2] = reco::PFTauDecayMode::tauDecay1ChargedPion1PiZero;
  supportedTauDecayModes_[3] = reco::PFTauDecayMode::tauDecay3ChargedPion0PiZero;
  supportedTauDecayModes_[4] = reco::PFTauDecayMode::tauDecayOther;
  
  numSupportedTauDecayModes_ = supportedTauDecayModes_.size();

  vRec_.ResizeTo(numSupportedTauDecayModes_);
  mapRecToGenTauDecayModes_.ResizeTo(numSupportedTauDecayModes_, numSupportedTauDecayModes_);
  mapRecToGenTauDecayModes_.Zero();
  vGen_.ResizeTo(numSupportedTauDecayModes_);
  vProb_.ResizeTo(numSupportedTauDecayModes_);

  std::vector<int>::const_iterator tauDecayModeOther_index 
    = std::find(supportedTauDecayModes_.begin(), supportedTauDecayModes_.end(), reco::PFTauDecayMode::tauDecayOther);
  assert(tauDecayModeOther_index != supportedTauDecayModes_.end());
  tauDecayModeOther_index_ = (*tauDecayModeOther_index);

//--- load histogram correlating reconstructed to generated hadronic decay modes
  edm::ParameterSet cfgMapRecToGenTauDecayModes = cfg.getParameter<edm::ParameterSet>("mapRecToGenTauDecayModes");

  std::string fileName_mapRecToGenTauDecayModes = cfgMapRecToGenTauDecayModes.getParameter<std::string>("fileName");
  std::string meName_mapRecToGenTauDecayModes = cfgMapRecToGenTauDecayModes.getParameter<std::string>("meName");

  TFile* file_mapRecToGenTauDecayModes = 0;
  try {
    file_mapRecToGenTauDecayModes = TFile::Open(fileName_mapRecToGenTauDecayModes.data());
  } catch (...) { 
    edm::LogError("SVfitTauLikelihoodPolarization") 
      << " Failed to open file = " << fileName_mapRecToGenTauDecayModes << " !!";
  }
  
//--- create "transfer matrix"
//
//    CV: In the histogram, the generated (reconstructed) tau decay modes are on the x-axis (y-axis);
//        the x-axis (y-axis) needs to be mapped to columns (rows) of the "transfer matrix".
//
//        The aim of the "transfer matrix" M is to map from reconstructed to generated ("true") tau lepton hadronic decay modes.
//        The mapping is implemented by matrix multiplication:
//         vGen = M * vRec,
//        where vRec is a vector encoding the reconstructed hadronic tau decay mode
//        (exactly one entry in this vector is 1, all other entries are 0)
//        and vGen gives the probabilities for the "true" decay mode of the tau lepton.
//
//        The required format of M is:
//         | p(gen=oneProngZeroPi0|rec=oneProngZeroPi0) .. p(gen=thrProgZeroPi0|rec=oneProngZeroPi0) p(gen=other|rec=oneProngZeroPi0) |
//         | p(gen=oneProngZeroPi0|rec=oneProngOnePi0 ) .. p(gen=thrProgZeroPi0|rec=oneProngOnePi0 ) p(gen=other|rec=oneProngOnePi0 ) |
//         | p(gen=oneProngZeroPi0|rec=oneProngTwoPi0 ) .. p(gen=thrProgZeroPi0|rec=oneProngTwoPi0 ) p(gen=other|rec=oneProngTwoPi0 ) |
//         | p(gen=oneProngZeroPi0|rec=thrProngZeroPi0) .. p(gen=thrProgZeroPi0|rec=thrProngZeroPi0) p(gen=other|rec=thrProngZeroPi0) |
//         | p(gen=oneProngZeroPi0|rec=other          ) .. p(gen=thrProgZeroPi0|rec=other          ) p(gen=other|rec=other          ) |
//
//        Note that **rows** of the "transfer matrix" need to be normalized to unit probability.
//
  TH2* histogram_mapRecToGenTauDecayModes 
    = dynamic_cast<TH2*>(file_mapRecToGenTauDecayModes->Get(meName_mapRecToGenTauDecayModes.data()));
  if ( histogram_mapRecToGenTauDecayModes ) {
    int numBinsX = histogram_mapRecToGenTauDecayModes->GetNbinsX();
    for ( int binIndex_x = 1; binIndex_x <= numBinsX; ++binIndex_x ) {
      std::string binLabel_x = histogram_mapRecToGenTauDecayModes->GetXaxis()->GetBinLabel(binIndex_x);

      size_t tauDecayMode_column = getSupportedTauDecayModeIndex(supportedTauDecayModes_, binLabel_x);
      assert(tauDecayMode_column >= 0 && tauDecayMode_column < numSupportedTauDecayModes_);

      int numBinsY = histogram_mapRecToGenTauDecayModes->GetNbinsY();
      for ( int binIndex_y = 1; binIndex_y <= numBinsY; ++binIndex_y ) {
	std::string binLabel_y = histogram_mapRecToGenTauDecayModes->GetYaxis()->GetBinLabel(binIndex_y);

	size_t tauDecayMode_row = getSupportedTauDecayModeIndex(supportedTauDecayModes_, binLabel_x);
	assert(tauDecayMode_row >= 0 && tauDecayMode_row < numSupportedTauDecayModes_);

	mapRecToGenTauDecayModes_(tauDecayMode_row, tauDecayMode_column) 
	  += histogram_mapRecToGenTauDecayModes->GetBinContent(binIndex_x, binIndex_y);
      }
    }
  }

  normalizeMatrixRows(mapRecToGenTauDecayModes_);
   
//--- close file from which histogram was loaded
  delete file_mapRecToGenTauDecayModes;

//--- initialize resolution and bias values for reconstructing
//    momentum fraction x = E(dist. pion)/E(vector meson)
//    of rho/a1 vector meson carried by "distinguishable" pion
//   (computed in laboratory frame)
  edm::ParameterSet cfgResolution = cfg.getParameter<edm::ParameterSet>("resolution");

  xVMrhoSigma_ = new TFormula("xVMrhoSigma", cfgResolution.getParameter<std::string>("xVMrhoSigma").data());
  xVMrhoBias_ = new TFormula("xVMrhoBias", cfgResolution.getParameter<std::string>("xVMrhoBias").data());
  xVMa1NeutralSigma_ = new TFormula("xVMa1NeutralSigma", cfgResolution.getParameter<std::string>("xVMa1NeutralSigma").data());
  xVMa1NeutralBias_ = new TFormula("xVMa1NeutralBias", cfgResolution.getParameter<std::string>("xVMa1NeutralBias").data());
  xVMa1ChargedSigma_ = new TFormula("xVMa1ChargedSigma", cfgResolution.getParameter<std::string>("xVMa1ChargedSigma").data());
  xVMa1ChargedBias_ = new TFormula("xVMa1ChargedBias", cfgResolution.getParameter<std::string>("xVMa1ChargedBias").data());
 
//--- create auxiliary classes for computation of vector meson line-shape integrals
  useCollApproxFormulas_ = cfg.getParameter<bool>("useCollApproxFormulas");

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
//    supported by the SVfitTauLikelihoodPolarization plugin
  edm::ParameterSet cfgLikelihoodPhaseSpace;
  typedef edmplugin::PluginFactory<SVfitLegLikelihoodBase<pat::Tau>* (const edm::ParameterSet&)> SVfitTauLikelihoodPluginFactory;
  likelihoodPhaseSpace_ = SVfitTauLikelihoodPluginFactory::get()->create("SVfitTauLikelihoodPhaseSpace", cfgLikelihoodPhaseSpace);

//--- initialize temporary variables,
//    defined to speed-up computations
  a1posMassTerm_ = (a1MesonMass2 + rhoMesonMass2)/(a1MesonMass*rhoMesonMass);
  a1posMassTerm2_ = square(a1posMassTerm_);
  a1negMassTerm_ = (a1MesonMass2 - rhoMesonMass2)/(a1MesonMass*rhoMesonMass);
  a1q_ = 0.5*TMath::Sqrt(a1MesonMass2 - 4.*chargedPionMass2);
  a1_8piDiv9_ = 8.*TMath::Pi()/9.;
  a1_16piDiv9_ = 2.*a1_8piDiv9_;
}

SVfitTauLikelihoodPolarization::~SVfitTauLikelihoodPolarization()
{
  delete xVMrhoSigma_;
  delete xVMrhoBias_;
  delete xVMa1NeutralSigma_;
  delete xVMa1NeutralBias_;
  delete xVMa1ChargedSigma_;
  delete xVMa1ChargedBias_;

  delete rhoLpolLineShape_;
  delete rhoTpolLineShape_;
  delete a1LpolLineShape_;
  delete a1TpolLineShape_;

  delete likelihoodPhaseSpace_;
}

double SVfitTauLikelihoodPolarization::negLogLikelihoodPolarized(
         const pat::Tau& tau, const SVfitLegSolution& solution, double tauLeptonPol) const
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
  std::cout << "<SVfitTauLikelihoodPolarization::negLogLikelihoodPolarized>:" << std::endl;
          
  int recTauDecayMode = tau.decayMode();

  vRec_.Zero();
  std::vector<int>::const_iterator tauDecayMode_index 
    = std::find(supportedTauDecayModes_.begin(), supportedTauDecayModes_.end(), recTauDecayMode);
  if ( tauDecayMode_index != supportedTauDecayModes_.end() ) {
    vRec_(*tauDecayMode_index) = 1.;
  } else {
    vRec_(tauDecayModeOther_index_) = 1.;
  }
  
  vGen_ = vRec_;
  vGen_ *= mapRecToGenTauDecayModes_;

  for ( size_t iDecayMode = 0; iDecayMode < numSupportedTauDecayModes_; ++iDecayMode ) {
    double probDecayMode = 0.;
    if      ( supportedTauDecayModes_[iDecayMode] == reco::PFTauDecayMode::tauDecay1ChargedPion0PiZero ) 
      probDecayMode = probOneProngZeroPi0s(tau, solution, tauLeptonPol);
    else if ( supportedTauDecayModes_[iDecayMode] == reco::PFTauDecayMode::tauDecay1ChargedPion1PiZero )
      probDecayMode = probOneProngOnePi0(tau, solution, tauLeptonPol);
    else if ( supportedTauDecayModes_[iDecayMode] == reco::PFTauDecayMode::tauDecay1ChargedPion1PiZero )
      probDecayMode = probOneProngTwoPi0s(tau, solution, tauLeptonPol);
    else if ( supportedTauDecayModes_[iDecayMode] == reco::PFTauDecayMode::tauDecay3ChargedPion0PiZero )
      probDecayMode = probThreeProngZeroPi0s(tau, solution, tauLeptonPol);
    else 
      probDecayMode = probOtherDecayMode(tau, solution, tauLeptonPol);
    vProb_(iDecayMode) = probDecayMode;
  }

  double prob = vGen_*vProb_;
  std::cout << "--> prob = " << prob << std::endl;

  if ( !(prob > 0.) ) {
    edm::LogWarning ("SVfitTauLikelihoodPolarization::operator()") 
      << " Unphysical solution --> returning very large negative number !!";
    return std::numeric_limits<float>::min();
  }
  
  double logLikelihood = TMath::Log(prob);
  std::cout << " -logLikelihood = " << -logLikelihood << std::endl;
  
  return -logLikelihood;
}

//
//-------------------------------------------------------------------------------
//

double SVfitTauLikelihoodPolarization::probOneProngZeroPi0s(
         const pat::Tau& tau, const SVfitLegSolution& solution, double tauLeptonPol) const
{
  std::cout << "<SVfitTauLikelihoodPolarization::probOneProngZeroPi0s>:" << std::endl;
          
  double prob = 0.;

  if ( !useCollApproxFormulas_ ) {
    double cosTheta = solution.cosThetaRest();
    double theta = TMath::ACos(cosTheta);
    double sinTheta = TMath::Sin(theta);
    prob = 0.5*(1. + tauLeptonPol*cosTheta)*sinTheta; // [1], formula (2.1)
  } else {    
    double z = solution.x();                          // tau lepton visible momentum fraction
    prob = (1. + tauLeptonPol*(2*z - 1.));            // [1], formula (2.4)
  }

  std::cout << "--> prob = " << prob << std::endl;

  return prob;
}

double SVfitTauLikelihoodPolarization::probOneProngOnePi0(
         const pat::Tau& tau, const SVfitLegSolution& solution, double tauLeptonPol) const
{
//--- compute likelihood for tau- --> rho- nu --> pi- pi0 nu decay

  std::cout << "<SVfitTauLikelihoodPolarization::probOneProngOnePi0>:" << std::endl;

  double cosTheta = solution.cosThetaRest();
  double theta = TMath::ACos(cosTheta);
  double z = tau.energy()/solution.p4().energy();
  
  double probTauDecayL = (*rhoLpolLineShape_)(theta, tauLeptonPol, z);
  double probTauDecayT = (*rhoTpolLineShape_)(theta, tauLeptonPol, z);
  std::cout << " probTauDecayL = " << probTauDecayL << ", probTauDecayT = " << probTauDecayT << std::endl;

//--- find "distinguishable" pion in tau-jet;
//    in case "distinguishable" pion cannot be found 
//   (e.g. in case "wrong" tau decay mode is reconstructed),
//    assume the "distinguishable" pion to be very soft
  const reco::Candidate* distPion = getDistPion(tau);
  double xMeasured = ( distPion != 0 ) ? distPion->energy()/tau.energy() : 0.;
  double thetaVMrho = solution.thetaVMrho();
  double cosThetaVMrho = TMath::Cos(thetaVMrho);
  double sinThetaVMrho = TMath::Sin(thetaVMrho);
  double xFitted = 0.5*(1. + TMath::Sqrt(1 - 4.*(chargedPionMass2/rhoMesonMass2))*cosThetaVMrho); // [2], formula (41)

  double probVMrhoDecayL, probVMrhoDecayT;
  if ( !useCollApproxFormulas_ ) {
    probVMrhoDecayL = 1.5*square(cosThetaVMrho)*sinThetaVMrho; // [2], formula (39)
    probVMrhoDecayT = 0.75*cube(sinThetaVMrho);                // [2], formula (40)
  } else {
    probVMrhoDecayL = 1.5*(2*xFitted - 1.);                    // [2], formula (39)
    probVMrhoDecayT = 3*xFitted*(1. - xFitted) ;               // [2], formula (40)
  }

  std::cout << " probVMrhoDecayL = " << probVMrhoDecayL << ", probVMrhoDecayT = " << probVMrhoDecayT << std::endl;

  double xSigma = xVMrhoSigma_->Eval(tau.pt());
  double xBias = xVMrhoBias_->Eval(tau.pt());
  std::cout << " xSigma = " << xSigma << ", xBias = " << xBias << std::endl;
  double probSmear = TMath::Gaus(xMeasured, xFitted + xBias, xSigma);
  std::cout << " probSmear = " << probSmear << std::endl;

  double probL = probTauDecayL*probVMrhoDecayL*probSmear;
  double probT = probTauDecayT*probVMrhoDecayT*probSmear;
  std::cout << "--> probL = " << probL << ", probT = " << probT << std::endl;

  return (probL + probT);
}

double SVfitTauLikelihoodPolarization::compVMa1x(
         double cosThetaVMa1, double sinThetaVMa1, double cosThetaVMa1r, double sinThetaVMa1r, double cosPhiVMa1r) const
{
  return (1./a1MesonMass)*(0.25*rhoMesonMass*(a1posMassTerm_ + a1negMassTerm_*cosThetaVMa1) 
                          + 0.5*a1q_*(a1posMassTerm_*cosThetaVMa1*cosThetaVMa1r + a1negMassTerm_*cosThetaVMa1r 
				     - 2.*sinThetaVMa1*sinThetaVMa1r*cosPhiVMa1r)); // [2], formula (48)
}

double SVfitTauLikelihoodPolarization::compVMa1DecayProbL(
	 double cosThetaVMa1, double sinThetaVMa1, double cosThetaVMa1r, double sinThetaVMa1r, double cosPhiVMa1r) const
{
  return (square(a1posMassTerm_*cosThetaVMa1*cosThetaVMa1r - 2.*sinThetaVMa1*sinThetaVMa1r*cosPhiVMa1r)/
	  (a1_8piDiv9_*(a1posMassTerm2_ + 8.)))*sinThetaVMa1*sinThetaVMa1r; // [2], formula (46)
}

double SVfitTauLikelihoodPolarization::compVMa1DecayProbT(
	 double cosThetaVMa1, double sinThetaVMa1, 
	 double cosThetaVMa1r, double sinThetaVMa1r, double cosPhiVMa1r, double sinPhiVMa1r) const
{
  return ((square(a1posMassTerm_*sinThetaVMa1*cosThetaVMa1r + 2.*cosThetaVMa1*sinThetaVMa1r*cosPhiVMa1r)
          + 4.*square(sinThetaVMa1r)*square(sinPhiVMa1r))/
	  (a1_16piDiv9_*(a1posMassTerm2_ + 8.)))*sinThetaVMa1*sinThetaVMa1r; // [2], formula (47)
}

double SVfitTauLikelihoodPolarization::probOneProngTwoPi0s(
         const pat::Tau& tau, const SVfitLegSolution& solution, double tauLeptonPol) const
{
//--- compute likelihood for tau- --> a1- nu --> pi- pi0 pi0 nu decay

  std::cout << "<SVfitTauLikelihoodPolarization::probOneProngTwoPi0s>:" << std::endl;

  double cosTheta = solution.cosThetaRest();
  double theta = TMath::ACos(cosTheta);
  double z = tau.energy()/solution.p4().energy();
  
  double probTauDecayL = (*a1LpolLineShape_)(theta, tauLeptonPol, z);
  double probTauDecayT = (*a1TpolLineShape_)(theta, tauLeptonPol, z);
  std::cout << " probTauDecayL = " << probTauDecayL << ", probTauDecayT = " << probTauDecayT << std::endl;

//--- find "distinguishable" pion in tau-jet;
//    in case "distinguishable" pion cannot be found 
//   (e.g. in case "wrong" tau decay mode is reconstructed),
//    assume the "distinguishable" pion to be very soft
  const reco::Candidate* distPion = getDistPion(tau);
  double xMeasured = ( distPion != 0 ) ? distPion->energy()/tau.energy() : 0.;
  double thetaVMa1 = solution.thetaVMa1();
  double cosThetaVMa1  = TMath::Cos(thetaVMa1);
  double sinThetaVMa1  = TMath::Sin(thetaVMa1);
  double thetaVMa1r = solution.thetaVMa1r();
  double cosThetaVMa1r = TMath::Cos(thetaVMa1r);
  double sinThetaVMa1r = TMath::Sin(thetaVMa1r);
  double phiVMa1r = solution.phiVMa1r();
  double cosPhiVMa1r   = TMath::Cos(phiVMa1r);
  double sinPhiVMa1r   = TMath::Sin(phiVMa1r);
  double xFitted = compVMa1x(cosThetaVMa1, sinThetaVMa1, cosThetaVMa1r, sinThetaVMa1r, cosPhiVMa1r);

//--- CV: only non-collinear approximation type formulas available in literature
//        for tau- --> a1- nu --> pi- pi0 pi0 nu decay
  double probVMa1DecayL = compVMa1DecayProbL(cosThetaVMa1, sinThetaVMa1, cosThetaVMa1r, sinThetaVMa1r, cosPhiVMa1r);
  double probVMa1DecayT = compVMa1DecayProbT(cosThetaVMa1, sinThetaVMa1, cosThetaVMa1r, sinThetaVMa1r, cosPhiVMa1r, sinPhiVMa1r);
  std::cout << " probVMa1DecayL = " << probVMa1DecayL << ", probVMa1DecayT = " << probVMa1DecayT << std::endl;

  double xSigma = xVMa1NeutralSigma_->Eval(tau.pt());
  double xBias = xVMa1NeutralBias_->Eval(tau.pt());
  std::cout << " xSigma = " << xSigma << ", xBias = " << xBias << std::endl;
  double probSmear = TMath::Gaus(xMeasured, xFitted + xBias, xSigma);
  std::cout << " probSmear = " << probSmear << std::endl;

  double probL = probTauDecayL*probVMa1DecayL*probSmear;
  double probT = probTauDecayT*probVMa1DecayT*probSmear;
  std::cout << "--> probL = " << probL << ", probT = " << probT << std::endl;

  return (probL + probT);
}

double SVfitTauLikelihoodPolarization::probThreeProngZeroPi0s(
         const pat::Tau& tau, const SVfitLegSolution& solution, double tauLeptonPol) const
{
//--- compute likelihood for tau- --> a1- nu --> pi- pi+ pi- nu decay

  std::cout << "<SVfitTauLikelihoodPolarization::probThreeProngZeroPi0s>:" << std::endl;

  double cosTheta = solution.cosThetaRest();
  double theta = TMath::ACos(cosTheta);
  double z = tau.energy()/solution.p4().energy();
  
  double probTauDecayL = (*a1LpolLineShape_)(theta, tauLeptonPol, z);
  double probTauDecayT = (*a1TpolLineShape_)(theta, tauLeptonPol, z);
  std::cout << " probTauDecayL = " << probTauDecayL << ", probTauDecayT = " << probTauDecayT << std::endl;

//--- find "distinguishable" pion in tau-jet;
//    in case "distinguishable" pion cannot be found 
//   (e.g. in case "wrong" tau decay mode is reconstructed),
//    assume the "distinguishable" pion to be very soft
  const reco::Candidate* distPion = getDistPion(tau);
  double xMeasured = ( distPion != 0 ) ? distPion->energy()/tau.energy() : 0.;
  double thetaVMa1 = solution.thetaVMa1();
  double cosThetaVMa1  = TMath::Cos(thetaVMa1);
  double sinThetaVMa1  = TMath::Sin(thetaVMa1);
  double thetaVMa1r = solution.thetaVMa1r();
  double cosThetaVMa1r = TMath::Cos(thetaVMa1r);
  double sinThetaVMa1r = TMath::Sin(thetaVMa1r);
  double phiVMa1r = solution.phiVMa1r();
  double cosPhiVMa1r   = TMath::Cos(phiVMa1r);
  double sinPhiVMa1r   = TMath::Sin(phiVMa1r);
  double xFitted = compVMa1x(cosThetaVMa1, sinThetaVMa1, cosThetaVMa1r, sinThetaVMa1r, cosPhiVMa1r);

//--- CV: only non-collinear approximation type formulas available in literature
//        for tau- --> a1- nu --> pi- pi+ pi- nu decay
  double probVMa1DecayL = compVMa1DecayProbL(cosThetaVMa1, sinThetaVMa1, cosThetaVMa1r, sinThetaVMa1r, cosPhiVMa1r);
  double probVMa1DecayT = compVMa1DecayProbT(cosThetaVMa1, sinThetaVMa1, cosThetaVMa1r, sinThetaVMa1r, cosPhiVMa1r, sinPhiVMa1r);
  std::cout << " probVMa1DecayL = " << probVMa1DecayL << ", probVMa1DecayT = " << probVMa1DecayT << std::endl;

  double xSigma = xVMa1ChargedSigma_->Eval(tau.pt());
  double xBias = xVMa1ChargedBias_->Eval(tau.pt());
  std::cout << " xSigma = " << xSigma << ", xBias = " << xBias << std::endl;
  double probSmear = TMath::Gaus(xMeasured, xFitted + xBias, xSigma);
  std::cout << " probSmear = " << probSmear << std::endl;

  double probL = probTauDecayL*probVMa1DecayL*probSmear;
  double probT = probTauDecayT*probVMa1DecayT*probSmear;
  std::cout << "--> probL = " << probL << ", probT = " << probT << std::endl;

  return (probL + probT);
}

double SVfitTauLikelihoodPolarization::probOtherDecayMode(
         const pat::Tau& tau, const SVfitLegSolution& solution, double tauLeptonPol) const
{
  return (*likelihoodPhaseSpace_)(tau, solution);
}

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_EDM_PLUGIN(SVfitTauLikelihoodBasePluginFactory, SVfitTauLikelihoodPolarization, "SVfitTauLikelihoodPolarization");
