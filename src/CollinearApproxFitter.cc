#include "TauAnalysis/CandidateTools/interface/CollinearApproxFitter.h"

#include "TauAnalysis/CandidateTools/interface/candidateAuxFunctions.h"

#include <TMath.h>

const double epsilon = 1.e-3;

double compX(double par)
{
//--- "smooth" mapping of Minuit fit parameter
//    into interval = 0..1 representing "physical" solutions
//    for quantity x defined as: energy of visible decay products / tau lepton energy

  return 0.5*(TMath::Sin(par) + 1.);
}

double compX_Err(double par, double parErr)
{
  double x = compX(par);

  double xUpShift = compX(par + parErr);
  double xDownShift = compX(par - parErr);

  return TMath::Sqrt(0.5*(TMath::Power(xUpShift - x, 2) + TMath::Power(xDownShift - x, 2)));
}

//
//-----------------------------------------------------------------------------------------------------------------------
//

CollinearApproxFitter::CollinearApproxFitter(const edm::ParameterSet& cfg)
  : resonanceMass_(cfg.getParameter<double>("resonanceMass")),
    resonanceWidth_(cfg.getParameter<double>("resonanceWidth")),
    minuit_(2)
{
  std::cout << "<CollinearApproxFitter::CollinearApproxFitter>:" << std::endl;
  std::cout << " disabling MINUIT output..." << std::endl;
  minuit_.SetPrintLevel(-1);
  minuit_.SetErrorDef(0.5);

  metResolutionPx_ = new TF1("metResolutionPx", cfg.getParameter<std::string>("metResolutionPx").data(), 0., 1.e+6);
  metResolutionPy_ = new TF1("metResolutionPy", cfg.getParameter<std::string>("metResolutionPy").data(), 0., 1.e+6);
}

CollinearApproxFitter::~CollinearApproxFitter()
{
  delete metResolutionPx_;
  delete metResolutionPy_;
}

void objectiveFcn(Int_t& nParameter, Double_t*, Double_t& fcn, Double_t* parameter, Int_t)
{
  CollinearApproxFitter* fitter = dynamic_cast<CollinearApproxFitter*>(gMinuit->GetObjectFit());
  if ( !fitter ) {
    edm::LogError("objectiveFcn") << "Call of gMinuit::GetObjectFit returned NULL pointer !!";
    return;
  }

  fcn = fitter->compChi2(parameter);
}

double CollinearApproxFitter::compChi2(Double_t* parameter) const
{
  //std::cout << "<CollinearApproxFitter::compChi2>:" << std::endl;
  //std::cout << " parameter[0] = " << parameter[0] << std::endl;
  //std::cout << " parameter[1] = " << parameter[1] << std::endl;

  double x1_fit = compX(parameter[0]);
  double x2_fit = compX(parameter[1]);

  if ( x1_fit == 0. ) x1_fit = epsilon;
  if ( x2_fit == 0. ) x2_fit = epsilon;

  double metPx_fit = (1./x1_fit - 1.)*p4leg1_.px() + (1./x2_fit - 1.)*p4leg2_.px();
  double metPy_fit = (1./x1_fit - 1.)*p4leg1_.py() + (1./x2_fit - 1.)*p4leg2_.py();

  reco::Candidate::LorentzVector p4Resonance = p4leg1_/x1_fit + p4leg2_/x2_fit;
  
  double chi2 = TMath::Power((metPx_fit - metPx_)/metPxErr_, 2) + TMath::Power((metPy_fit - metPy_)/metPyErr_, 2);
  chi2 += TMath::Power((p4Resonance.mass() - resonanceMass_)/resonanceWidth_, 2);
  chi2 /= 2.;
  //std::cout << " chi2 = " << chi2 << std::endl;

  return chi2;
}

CollinearApproxCompatibility CollinearApproxFitter::fit(const reco::Candidate::LorentzVector& p4leg1, 
							const reco::Candidate::LorentzVector& p4leg2, edm::Ptr<reco::MET> met)
{
  //std::cout << "<CollinearApproxFitter::fit>:" << std::endl;

  p4leg1_ = p4leg1;
  p4leg2_ = p4leg2;
  met_ = met;

//--- compute initial values for solutions x1, x2 of collinear approximation
  double x1, x2;
  compX1X2byCollinearApprox(x1, x2, p4leg1.px(), p4leg1.py(), p4leg2.px(), p4leg2.py(), met->px(), met->py());
  //std::cout << " x1 (initial) = " << x1 << std::endl;
  //std::cout << " x2 (initial) = " << x2 << std::endl;

  bool isX1withinPhysRange, isX2withinPhysRange;
  double x1phys = getPhysX(x1, isX1withinPhysRange);
  double x2phys = getPhysX(x2, isX2withinPhysRange);

  minuit_.DefineParameter(0, "par1", TMath::ASin(x1phys), 0.1, -1.e+3, +1.e+3);
  minuit_.DefineParameter(1, "par2", TMath::ASin(x2phys), 0.1, -1.e+3, +1.e+3);

  metPx_ = met_->px();  
  metPxErr_ = metResolutionPx_->Eval(met->sumEt());
  metPy_ = met_->py();
  metPyErr_ = metResolutionPy_->Eval(met->sumEt());
  //std::cout << "metPx = " << metPx_ << std::endl;
  //std::cout << "metPxErr = " << metPxErr_ << std::endl;
  //std::cout << "metPy = " << metPy_ << std::endl;
  //std::cout << "metPyErr = " << metPyErr_ << std::endl;

  minuit_.SetObjectFit(this);
  minuit_.SetFCN(objectiveFcn);
  gMinuit = &minuit_;

  std::cout << " starting Migrad minimization..." << std::endl;
  minuit_.Migrad();

  int minuitStatus = minuit_.GetStatus();
  //std::cout << " Migrad status = " << minuitStatus << std::endl;

  Double_t gradients[2];
  Double_t parameter[2];
  Double_t dummy;
  minuit_.GetParameter(0, parameter[0], dummy);
  minuit_.GetParameter(1, parameter[1], dummy);
  double minuitChi2;
  minuit_.Eval(2, gradients, minuitChi2, parameter, 0);
  //std::cout << " Chi2 = " << minuitChi2 << std::endl;

  double par1, par1Err, par2, par2Err;
  minuit_.GetParameter(0, par1, par1Err);
  minuit_.GetParameter(1, par2, par2Err);

  x1 = compX(par1);
  double x1Err = compX_Err(par1, par1Err);
  x2 = compX(par2);
  double x2Err = compX_Err(par2, par2Err);
  //std::cout << " x1 = " << x1 << " +/- " << x1Err << std::endl;
  //std::cout << " x2 = " << x2 << " +/- " << x2Err << std::endl;

  CollinearApproxCompatibility collinearApproxCompatibility;
  collinearApproxCompatibility.setMinuitStatus(minuitStatus);
  collinearApproxCompatibility.setMinuitChi2(minuitChi2);
  collinearApproxCompatibility.setX1(x1);
  collinearApproxCompatibility.setX1Err(x1Err);
  collinearApproxCompatibility.setX2(x2);
  collinearApproxCompatibility.setX2Err(x2Err);

  return collinearApproxCompatibility;
}
