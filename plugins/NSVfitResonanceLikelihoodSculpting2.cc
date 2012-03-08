#include "TauAnalysis/CandidateTools/plugins/NSVfitResonanceLikelihoodSculpting2.h"

#include "TauAnalysis/CandidateTools/interface/svFitAuxFunctions.h"

#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TAxis.h>
#include <TArrayD.h>
#include <TMath.h>

#include <string>

using namespace SVfit_namespace;

namespace
{
  TH1* getHistogram(TFile* inputFile, const std::string& histogramName)
  {
    TH1* histogram = (TH1*)inputFile->Get(histogramName.data());
    if ( !histogram) 
      throw cms::Exception("getHistogram") 
	<< " Failed to load histogram = " << histogramName.data() << " from file = " << inputFile->GetName() << " !!\n";
    return histogram;
  }

  TArrayD getBinning(const TH1* histogram)
  {
    TAxis* xAxis = histogram->GetXaxis();
    int numBins = xAxis->GetNbins();
    TArrayD binning(numBins + 1);
    for ( int iBin = 0; iBin < numBins; ++iBin ) {
      binning[iBin] = xAxis->GetBinLowEdge(iBin + 1); // CV: TAxis bins start @ index = 1
      binning[iBin + 1] = xAxis->GetBinUpEdge(iBin + 1);
    }
    //for ( int iBin = 0; iBin <= numBins; ++iBin ) {
    //  std::cout << "binning[" << iBin << "] = " << binning[iBin] << std::endl;
    //}
    return binning;
  }

  double getMinX(double tauLeptonEnergy, double visTheta, double visMass, double minVisPt)
  { 
    double minVisP  = minVisPt/TMath::Sin(visTheta);
    double minVisEn = TMath::Sqrt(square(minVisP) + square(visMass));
    double minX     = minVisEn/tauLeptonEnergy;
    return minX;
  }
}

NSVfitResonanceLikelihoodSculpting2::NSVfitResonanceLikelihoodSculpting2(const edm::ParameterSet& cfg)
  : NSVfitResonanceLikelihood(cfg),
    inputFile_(0),
    histogramX1vsX2_(0)
{
//--- load probability distributions of X1 = visPt1/tauLeptonPt1 and X2 = visPt2/tauLeptonPt2 
  edm::FileInPath inputFileName = cfg.getParameter<edm::FileInPath>("inputFileName");
  if ( !inputFileName.isLocal()) throw cms::Exception("NSVfitResonanceLikelihoodSculpting2") 
    << " Failed to find File = " << inputFileName << " !!\n";
  inputFile_ = new TFile(inputFileName.fullPath().data());
  
  std::string histogramNameX1 = cfg.getParameter<std::string>("histogramNameX1");
  TH1* histogramX1 = getHistogram(inputFile_, histogramNameX1);
  TArrayD binningX1(getBinning(histogramX1));
  int numBinsX1 = binningX1.GetSize() - 1;
  //std::cout << "numBinsX1 = " << numBinsX1 << std::endl;

  std::string histogramNameX2 = cfg.getParameter<std::string>("histogramNameX2");
  TH1* histogramX2 = getHistogram(inputFile_, histogramNameX2);
  TArrayD binningX2 = getBinning(histogramX2);
  int numBinsX2 = binningX2.GetSize() - 1;
  //std::cout << "numBinsX2 = " << numBinsX2 << std::endl;

  std::string histogramNameX1vsX2_tmp = Form("%s_histogramX1vsX2_tmp", pluginName_.data());
  TH2* histogramX1vsX2_tmp = 
    new TH2D(histogramNameX1vsX2_tmp.data(), 
	     histogramNameX1vsX2_tmp.data(), numBinsX1, binningX1.GetArray(), numBinsX2, binningX2.GetArray());
//--- compute product of probabilities p(x1 = X1 && x2 = X2) = p(x1 = X1) * p(x2 = X2)
//   (assuming that X1 and X2 are uncorrelated, which in good 
//    approximation they actually are; checked with Monte Carlo)
  for ( int iBinX1 = 1; iBinX1 <= numBinsX1; ++iBinX1 ) {
    double binContentX1 = histogramX1->GetBinContent(iBinX1);
    for ( int iBinX2 = 1; iBinX2 <= numBinsX2; ++iBinX2 ) {
      double binContentX2 = histogramX2->GetBinContent(iBinX2);
      double binContentX1timesX2 = binContentX1*binContentX2;      
      histogramX1vsX2_tmp->SetBinContent(iBinX1, iBinX2, binContentX1timesX2);
    }
  }
  if ( histogramX1vsX2_tmp->Integral() > 0. ) histogramX1vsX2_tmp->Scale(1./histogramX1vsX2_tmp->Integral());

  //TFile* outputFile_a = new TFile("histogramX1vsX2a.root", "RECREATE");
  //histogramX1vsX2_tmp->Write();
  //delete outputFile_a;

//--- sum bin-contents such that bin-content(X1,X2) corresponds to the probability
//    p(x1 >= X1 && x2 > = X2) for X1 and X2 to be above some threshold
//   (the thresholds will be X1 = visPtCut1/tauLeptonPt1 and X2 = visPtCut2/tauLeptonPt2 later,
//    when computing the likelihood correction for the effect of visible Pt cuts...)
  std::string histogramNameX1vsX2 = Form("%s_histogramX1vsX2", pluginName_.data());
  histogramX1vsX2_ = 
    new TH2D(histogramNameX1vsX2.data(), 
	     histogramNameX1vsX2.data(), numBinsX1, binningX1.GetArray(), numBinsX2, binningX2.GetArray());
  for ( int iBinX1 = 1; iBinX1 <= numBinsX1; ++iBinX1 ) {
    for ( int iBinX2 = 1; iBinX2 <= numBinsX2; ++iBinX2 ) {
      double binContentSum = histogramX1vsX2_tmp->Integral(iBinX1, numBinsX1, iBinX2, numBinsX2);
      histogramX1vsX2_->SetBinContent(iBinX1, iBinX2, binContentSum);
    }
  }

  delete histogramX1vsX2_tmp;

  //TFile* outputFile_b = new TFile("histogramX1vsX2b.root", "RECREATE");
  //histogramX1vsX2_->Write();
  //delete outputFile_b;

  minVisPt1_    = cfg.getParameter<double>("minVisPt1");
  minVisPt2_    = cfg.getParameter<double>("minVisPt2");
  minVisPt1alt_ = cfg.getParameter<double>("minVisPt1alt");
  minVisPt2alt_ = cfg.getParameter<double>("minVisPt2alt");

  power_ = cfg.getParameter<double>("power");
}

NSVfitResonanceLikelihoodSculpting2::~NSVfitResonanceLikelihoodSculpting2()
{
  delete histogramX1vsX2_;
  delete inputFile_;
}

double NSVfitResonanceLikelihoodSculpting2::operator()(const NSVfitResonanceHypothesis* resonance) const 
{
  if ( verbosity_ ) {
    std::cout << "<NSVfitResonanceLikelihoodSculpting2::operator()>:" << std::endl;
    std::cout << " mass = " << resonance->p4_fitted().mass() << std::endl;
  }

  assert(resonance);

  if ( resonance->numDaughters() != 2 ) {
    throw cms::Exception("NSVfitResonanceLikelihoodSculpting2::operator()")
      << " Resonance hypothesis passed as function argument has " << resonance->numDaughters()
      << " daughter particles, exactly two expected !!\n";
  }

  const NSVfitSingleParticleHypothesis* daughter1 = resonance->daughter(0);
  const NSVfitSingleParticleHypothesis* daughter2 = resonance->daughter(1);

  // CV: protection against "unphysical" directions of visible decay products
  //    (|sin(theta)| < 1.e-3 should never happen, because it is far outside of detector acceptance)
  if ( TMath::Abs(TMath::Sin(daughter1->p4().theta())) < 1.e-3 ||
       daughter1->p4_fitted().energy()                 < 1.e-3 ||
       TMath::Abs(TMath::Sin(daughter2->p4().theta())) < 1.e-3 ||
       daughter2->p4_fitted().energy()                 < 1.e-3 ) return 0.;

  double leg1En       = daughter1->p4_fitted().energy();
  double leg1VisTheta = daughter1->p4().theta();
  double leg1VisMass  = daughter1->p4().mass();
  double minX1        = getMinX(leg1En, leg1VisTheta, leg1VisMass, minVisPt1_);
  if ( verbosity_ ) std::cout << "x1 = " << (daughter1->p4().energy()/leg1En) << " (min = " << minX1 << ")" << std::endl;

  double leg2En       = daughter2->p4_fitted().energy();
  double leg2VisTheta = daughter2->p4().theta();
  double leg2VisMass  = daughter2->p4().mass();
  double minX2        = getMinX(leg2En, leg2VisTheta, leg2VisMass, minVisPt2_);
  if ( verbosity_ ) std::cout << "x2 = " << (daughter2->p4().energy()/leg2En) << " (min = " << minX2 << ")" << std::endl;

  double prob = 1.;
  if ( minVisPt1alt_ <= 0. || minVisPt2alt_ <= 0. ) {
    prob = getProbabilityX1andX2gt(minX1, minX2);
  } else {
    double minX1alt = getMinX(leg1En, leg1VisTheta, leg1VisMass, minVisPt1alt_);
    double minX2alt = getMinX(leg2En, leg2VisTheta, leg2VisMass, minVisPt2alt_);
    //
    //          leg2Pt
    //             ^
    //             |        |            |
    //             |        |     R1     |     R2 
    //             |        |            |
    // mivPisPt2alt+        +------------+-------------
    //             |        |            |
    //             |        |     R3     |     R4
    //             |        |            |
    //    minVisPt2+        +------------+-------------
    //             |
    //             |
    //             +--------+------------+-------------> leg1Pt
    //                minVisPt1alt   minVisPt1
    //
    //
    // CV: probability for visible tau decay products to pass cuts,
    //       p((leg1Pt > minVisPt1 && leg2Pt > minVisPt2) || (leg1Pt > minVisPt1alt && leg2Pt > minVisPt2alt))
    //      = R1 + R2 + R4 = (R1 + R2) + (R2 + R4) - R2
    //      =  p((leg1Pt > min(minVisPt1,minVisPt1alt) && leg2Pt > max(minVisPt2,minVisPt2alt))
    //       + p((leg1Pt > max(minVisPt1,minVisPt1alt) && leg2Pt > min(minVisPt2,minVisPt2alt))
    //       - p((leg1Pt > max(minVisPt1,minVisPt1alt) && leg2Pt > max(minVisPt2,minVisPt2alt))
    //
    prob = getProbabilityX1andX2gt(TMath::Min(minX1, minX1alt), TMath::Max(minX2, minX2alt))
          + getProbabilityX1andX2gt(TMath::Max(minX1, minX1alt), TMath::Min(minX2, minX2alt)) 
          - getProbabilityX1andX2gt(TMath::Max(minX1, minX1alt), TMath::Max(minX2, minX2alt));
  }
  
  double nllCorr = 0.;
  if ( prob > 0. ) {
    double probCorr = 1./prob;
    if ( verbosity_ ) std::cout << "probCorr = " << probCorr << std::endl;
    nllCorr = -power_*TMath::Log(probCorr);
  } else {
    nllCorr = std::numeric_limits<float>::min();
  }
  
  if ( this->verbosity_ ) std::cout << "--> nllCorr = " << nllCorr << std::endl;
  
  return nllCorr;
}

double interpolateLinearly(double x, double x1, double f1, double x2, double f2)
{
  double df_by_dx = (f2 - f1)/(x2 - x1);
  double f = f1 + (x - x1)*df_by_dx;
  return f;
}

double NSVfitResonanceLikelihoodSculpting2::getProbabilityX1andX2gt(double minX1, double minX2) const
{
  if ( minX1 < 0. ) minX1 = 0.;
  if ( minX1 > 1. ) minX1 = 1.;
  if ( minX2 < 0. ) minX2 = 0.;
  if ( minX2 > 1. ) minX2 = 1.;

  TAxis* axisX1 = histogramX1vsX2_->GetXaxis();
  int binX1_center = axisX1->FindBin(minX1);
  double binCenterX1 = axisX1->GetBinCenter(binX1_center);  
  int numBinsX1 = axisX1->GetNbins();

  int binX1_left  = ( minX1      >  binCenterX1  ) ? binX1_center : (binX1_center - 1);
  int binX1_right = ( binX1_left != binX1_center ) ? binX1_center : (binX1_center + 1);
  if ( binX1_left < 1 ) {
    binX1_left  = 1;
    binX1_right = binX1_left + 1;
  }
  if ( binX1_right > numBinsX1 ) {
    binX1_right = numBinsX1;
    binX1_left  = binX1_right - 1;
  }
  double binCenterX1_left  = axisX1->GetBinCenter(binX1_left);
  double binCenterX1_right = axisX1->GetBinCenter(binX1_right);

  TAxis* axisX2 = histogramX1vsX2_->GetYaxis();
  int binX2_center = axisX2->FindBin(minX2);
  double binCenterX2 = axisX2->GetBinCenter(binX2_center);
  int numBinsX2 = axisX2->GetNbins();

  int binX2_bottom = ( minX2        > binCenterX2   ) ? binX2_center : (binX2_center - 1);
  int binX2_top    = ( binX2_bottom != binX2_center ) ? binX2_center : (binX2_center + 1);
  if ( binX2_bottom < 1 ) {
    binX2_bottom = 1;
    binX2_top = binX2_bottom + 1;
  }
  if ( binX2_top > numBinsX2 ) {
    binX2_top = numBinsX2;
    binX2_bottom  = binX2_top - 1;
  }
  double binCenterX2_bottom = axisX2->GetBinCenter(binX2_bottom);
  double binCenterX2_top    = axisX2->GetBinCenter(binX2_top);

  double binContent_CC = histogramX1vsX2_->GetBinContent(binX1_center, binX2_center);
  double binContent_BL = histogramX1vsX2_->GetBinContent(binX1_left,   binX2_bottom);
  double binContent_BR = histogramX1vsX2_->GetBinContent(binX1_right,  binX2_bottom);
  double binContent_TL = histogramX1vsX2_->GetBinContent(binX1_left,   binX2_top);
  double binContent_TR = histogramX1vsX2_->GetBinContent(binX1_right,  binX2_top);

  double binContent_BC = interpolateLinearly(minX1, binCenterX1_left, binContent_BL, binCenterX1_right, binContent_BR);
  double binContent_TC = interpolateLinearly(minX1, binCenterX1_left, binContent_TL, binCenterX1_right, binContent_TR);

  double prob = interpolateLinearly(minX2, binCenterX2_bottom, binContent_BC, binCenterX2_top, binContent_TC); 

  if ( verbosity_ ) {
    std::cout << "binContent (CC) = " << binContent_CC << std::endl;
    std::cout << "binContent (BL) = " << binContent_BL << std::endl;
    std::cout << "binContent (BR) = " << binContent_BR << std::endl;
    std::cout << "binContent (TL) = " << binContent_TL << std::endl;
    std::cout << "binContent (TR) = " << binContent_TR << std::endl;
    std::cout << "binContent (BC) = " << binContent_BC << std::endl;
    std::cout << "binContent (TC) = " << binContent_TC << std::endl;
    std::cout << "--> prob = " << prob << std::endl;
  }

  return prob;
}

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_EDM_PLUGIN(NSVfitResonanceLikelihoodPluginFactory, NSVfitResonanceLikelihoodSculpting2, "NSVfitResonanceLikelihoodSculpting2");
