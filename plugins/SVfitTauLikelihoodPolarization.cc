#include "TauAnalysis/CandidateTools/plugins/SVfitTauLikelihoodPolarization.h"

#include "DataFormats/TauReco/interface/PFTauDecayMode.h"

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
    if ( SVfit_namespace::getTauDecayModeName(iDecayMode) == tauDecayModeName ) return iDecayMode;
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
 
//--- generic "phase-space" plugin to be used for computing likelihood 
//    for hadronic tau decays not in the list of decay modes 
//    supported by the SVfitTauLikelihoodPolarization plugin
  edm::ParameterSet cfgLikelihoodPhaseSpace;
  typedef edmplugin::PluginFactory<SVfitLegLikelihoodBase<pat::Tau>* (const edm::ParameterSet&)> SVfitTauLikelihoodPluginFactory;
  likelihoodPhaseSpace_ = SVfitTauLikelihoodPluginFactory::get()->create("SVfitTauLikelihoodPhaseSpace", cfgLikelihoodPhaseSpace);
  
  useCollApproxFormulas_ = cfg.getParameter<bool>("useCollApproxFormulas");
}

SVfitTauLikelihoodPolarization::~SVfitTauLikelihoodPolarization()
{
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
      probDecayMode = probOneProngZeroPi0(tau, solution, tauLeptonPol);
    else if ( supportedTauDecayModes_[iDecayMode] == reco::PFTauDecayMode::tauDecay1ChargedPion1PiZero )
      probDecayMode = probOneProngOnePi0(tau, solution, tauLeptonPol);
    else if ( supportedTauDecayModes_[iDecayMode] == reco::PFTauDecayMode::tauDecay1ChargedPion1PiZero )
      probDecayMode = probOneProngTwoPi0(tau, solution, tauLeptonPol);
    else if ( supportedTauDecayModes_[iDecayMode] == reco::PFTauDecayMode::tauDecay3ChargedPion0PiZero )
      probDecayMode = probThreeProngZeroPi0(tau, solution, tauLeptonPol);
    else 
      probDecayMode = probOtherDecayMode(tau, solution, tauLeptonPol);
    vProb_(iDecayMode) = probDecayMode;
  }

  double prob = vGen_*vProb_;
  std::cout << " prob = " << prob << std::endl;

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

double SVfitTauLikelihoodPolarization::probOneProngZeroPi0(
         const pat::Tau& tau, const SVfitLegSolution& solution, double tauLeptonPol) const
{
  // dummy implementation
  return 0.;
}

double SVfitTauLikelihoodPolarization::probOneProngOnePi0(
         const pat::Tau& tau, const SVfitLegSolution& solution, double tauLeptonPol) const
{
  // dummy implementation
  return 0.;
}

double SVfitTauLikelihoodPolarization::probOneProngTwoPi0(
         const pat::Tau& tau, const SVfitLegSolution& solution, double tauLeptonPol) const
{
  // dummy implementation
  return 0.;
}

double SVfitTauLikelihoodPolarization::probThreeProngZeroPi0(
         const pat::Tau& tau, const SVfitLegSolution& solution, double tauLeptonPol) const
{
  // dummy implementation
  return 0.;
}

double SVfitTauLikelihoodPolarization::probOtherDecayMode(
         const pat::Tau& tau, const SVfitLegSolution& solution, double tauLeptonPol) const
{
  return (*likelihoodPhaseSpace_)(tau, solution);
}

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_EDM_PLUGIN(SVfitTauLikelihoodBasePluginFactory, SVfitTauLikelihoodPolarization, "SVfitTauLikelihoodPolarization");
