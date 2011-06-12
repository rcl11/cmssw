#include "TauAnalysis/CandidateTools/plugins/NSVfitResonanceLikelihoodSculpting.h"

#include "TauAnalysis/CandidateTools/interface/NSVfitAlgorithmBase.h"
#include "TauAnalysis/CandidateTools/interface/svFitAuxFunctions.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/FileInPath.h"
#include "FWCore/Utilities/interface/Exception.h"

#include <string>

#include "TFile.h"
#include "RooAbsPdf.h"
#include "RooAbsReal.h"
#include "RooArgSet.h"
#include "RooWorkspace.h"
#include "RooRealVar.h"

NSVfitResonanceLikelihoodSculpting::NSVfitResonanceLikelihoodSculpting(
    const edm::ParameterSet& cfg)
  : NSVfitResonanceLikelihood(cfg) {
  power_ = cfg.getParameter<double>("power");

  std::string workspaceName = cfg.getParameter<std::string>("workspaceName");
  std::string modelName = cfg.getParameter<std::string>("modelName");
  std::string xVarName = cfg.getParameter<std::string>("xVarName");
  std::string massVarName = cfg.getParameter<std::string>("massVarName");
  edm::FileInPath pdfFilePath = cfg.getParameter<edm::FileInPath>("pdfFile");

  if (!pdfFilePath.isLocal()) {
    throw cms::Exception("NSVfitPtBalanceCompositePdf")
      << " Failed to find file = " << pdfFilePath.fullPath() << " !!\n";
  }

  TFile* file = TFile::Open(pdfFilePath.fullPath().data(), "READ");
  RooWorkspace* ws = dynamic_cast<RooWorkspace*>(
      file->Get(workspaceName.c_str()));

  RooAbsPdf* pdf = ws->pdf(modelName.c_str());
  if (!pdf) {
    throw cms::Exception("NSVfitPtBalanceCompositePdf")
      << " Couldn't load pdf " << modelName << " from "
      << pdfFilePath.fullPath().data() << std::endl;
  }

  RooRealVar* xVar = ws->var(xVarName.c_str());
  if (!xVar) {
    throw cms::Exception("NSVfitPtBalanceCompositePdf")
      << " Couldn't load dependent var " << xVarName << " from "
      << pdfFilePath.fullPath().data() << std::endl;
  }

  RooRealVar* massVar = ws->var(massVarName.c_str());
  if (!massVar) {
    throw cms::Exception("NSVfitPtBalanceCompositePdf")
      << " Couldn't load mass var " << massVarName << " from "
      << pdfFilePath.fullPath().data() << std::endl;
  }

  pdf_ = NSVfitCachingPdfWrapper(pdf, xVar, massVar, 400, 0, 1, 300, 0, 600);

  delete file;
}

NSVfitResonanceLikelihoodSculpting::~NSVfitResonanceLikelihoodSculpting() { }

void NSVfitResonanceLikelihoodSculpting::beginJob(NSVfitAlgorithmBase* algorithm) const
{
  algorithm->requestFitParameter("allTauDecays", nSVfit_namespace::kTau_visEnFracX, pluginName_);
  algorithm->requestFitParameter("allTauDecays", nSVfit_namespace::kTau_phi_lab,    pluginName_);
  algorithm->requestFitParameter("allLeptons",   nSVfit_namespace::kLep_shiftEn,    pluginName_);
  algorithm->requestFitParameter("allNeutrinos", nSVfit_namespace::kNu_energy_lab,  pluginName_);
  algorithm->requestFitParameter("allNeutrinos", nSVfit_namespace::kNu_phi_lab,     pluginName_);
}

double
NSVfitResonanceLikelihoodSculpting::operator()(const NSVfitResonanceHypothesis* hypothesis) const
{
//--- compute negative log-likelihood for two tau leptons
//    to have transverse momenta leg1Pt, leg2Pt

  if ( this->verbosity_ )
    std::cout << "<NSVfitLikelihoodDiTauSculpting::operator()>:" << std::endl;

  double diTauMass = hypothesis->p4_fitted().mass();
  double visMass = hypothesis->p4().mass();
  double scaledVisMass = visMass/diTauMass;

  double nll = -TMath::Log(pdf_.getVal(scaledVisMass, diTauMass));

  if ( this->verbosity_ ) {
    std::cout << " diTauMass = " << diTauMass
      << " visMass = " << visMass
      << " scaledVisMass = " << scaledVisMass
      << std::endl;
  }
  return power_*nll;
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_EDM_PLUGIN(NSVfitResonanceLikelihoodPluginFactory,
    NSVfitResonanceLikelihoodSculpting, "NSVfitResonanceLikelihoodSculpting");
