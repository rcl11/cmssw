#include "TauAnalysis/CandidateTools/interface/SVMethodLikelihoods.h"
#include "math.h"
#include <TMath.h>
#include <TMatrixD.h>
#include <TVectorD.h>
#include <TMatrixDSym.h>
#include <TRotation.h>
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Muon.h"

namespace TauVertex {

inline double square(double x) { return x*x; }
inline double cube(double x) { return x*x*x; }

inline double nlGaussianNorm(double sigma, int dimension=1) { 
   // Norm of Gaussian = 
   //  1/(sigma*(2*pi)^(k/2))
   return (log(sigma) + (dimension/2.0)*log(TMath::TwoPi()));
}

double nllPointGiven2DError(const GlobalPoint& point, const GlobalPoint& central, const GlobalError& error)
{
   AlgebraicVector3 displacement(point.x()-central.x(), point.y()-central.y(), point.z()-central.z());
   AlgebraicSymMatrix33 errorMatrix(error.matrix_new());
   // Load into TMatrix
   TMatrixDSym matrix(3);
   matrix[0][0] = errorMatrix(0,0);
   matrix[0][1] = errorMatrix(0,1);
   matrix[0][2] = errorMatrix(0,2);
   matrix[1][0] = errorMatrix(1,0);
   matrix[1][1] = errorMatrix(1,1);
   matrix[1][2] = errorMatrix(1,2);
   matrix[2][0] = errorMatrix(2,0);
   matrix[2][1] = errorMatrix(2,1);
   matrix[2][2] = errorMatrix(2,2);

   //std::cout << "Error matrix is " << std::endl;
   //matrix.Print();

   // Get Eigenstuff
   TVectorD eigenvalues(3);

   //std::cout << "Eigenvalues " << eigenvalues[0] << ", " << eigenvalues[1] << ", " << eigenvalues[2] << std::endl;

   TMatrixD eigenvectors = matrix.EigenVectors(eigenvalues);

   //std::cout << "Eigenvectors " << std::endl;
   //eigenvectors.Print();

   // Get the lowest eigenvector
   TVector3 lowestEigenvector(eigenvectors[0][2], eigenvectors[1][2], eigenvectors[2][2]);

   //std::cout << "Lowest eigenvector " << std::endl;
   //lowestEigenvector.Print();

   // Make a rotation such that this eigenvector is aligned along Z
   TRotation rot;
   rot.SetZAxis(lowestEigenvector);

   // Okay, now let's move back to our OTHER set of matrix types :(
   AlgebraicMatrix33 rotSMatrix;
   rotSMatrix(0,0) = rot[0][0];
   rotSMatrix(0,1) = rot[0][1];
   rotSMatrix(0,2) = rot[0][2];
   rotSMatrix(1,0) = rot[1][0];
   rotSMatrix(1,1) = rot[1][1];
   rotSMatrix(1,2) = rot[1][2];
   rotSMatrix(2,0) = rot[2][0];
   rotSMatrix(2,1) = rot[2][1];
   rotSMatrix(2,2) = rot[2][2];

   // Transform our displacement vector
   AlgebraicVector3 rotDisp = rotSMatrix*displacement;

   //std::cout << "New disp. vector " << rotDisp(0) << ", " << rotDisp(1) << ", " << rotDisp(2) << std::endl;

   // Transform our error matrix
   AlgebraicMatrix33 rotError = ROOT::Math::Transpose(rotSMatrix)*errorMatrix*rotSMatrix;
   //AlgebraicMatrix33 rotError = rotSMatrix*errorMatrix*ROOT::Math::Transpose(rotSMatrix);

   //std::cout << "New error matrix " << rotError << std::endl;

   // Now throw away the Z components
   AlgebraicMatrix22 rotError2;
   rotError2(0,0) = rotError(0,0);
   rotError2(0,1) = rotError(0,1);
   rotError2(1,0) = rotError(1,0);
   rotError2(1,1) = rotError(1,1);

   AlgebraicVector2 rotDisp2;
   rotDisp2(0) = rotDisp(0);
   rotDisp2(1) = rotDisp(1);

   double determinant = 0;
   bool detOK = rotError2.Det2(determinant);

   //std::cout << "Determinant " << determinant << std::endl;

   int inversionResult = rotError2.Invert();
   if(!inversionResult || !detOK)
   {
      edm::LogWarning("SVMethodLikelihoods") << "2D MATRIX ERROR INVERSION FAILED\n" << errorMatrix;
   }

   double expResult = ROOT::Math::Dot(rotDisp2, rotError2*rotDisp2)/2.0;
   double normResult = nlGaussianNorm(sqrt(determinant), 2);
   return expResult + normResult;
}

double nllPointGiven3DError(const GlobalPoint& point, const GlobalPoint& central, const GlobalError& error)
{
   AlgebraicVector3 displacement(point.x()-central.x(), point.y()-central.y(), point.z()-central.z());
   AlgebraicSymMatrix33 errorMatrix(error.matrix_new());
   double determinant = 0;
   bool detOK = errorMatrix.Det2(determinant);
   // Invert the matrix in place
   int inversionResult = errorMatrix.Invert();
   if(!inversionResult || !detOK)
   {
      edm::LogWarning("SVMethodLikelihoods") << "3D MATRIX ERROR INVERSION FAILED\n" << errorMatrix;
   }
   // Now do d^T*e^-1*d
   double result = ROOT::Math::Dot(displacement, errorMatrix*displacement);
   return result/2.0 + nlGaussianNorm(sqrt(determinant), 3);
}

double nllPointGivenTrack(const TrajectoryStateClosestToPoint& tcsp)
{
   return nllPointGiven2DError(tcsp.referencePoint(), tcsp.position(), 
         tcsp.theState().cartesianError().position());
}

double nllPointGivenVertex(const GlobalPoint& point, const TransientVertex& vertex)
{
   return nllPointGiven3DError(point, vertex.position(), vertex.positionError());
}

double nllTauDecayLengthGivenMomentum(double length, double momentum)
{
   // see PDG kinematics section (http://www-pdg.lbl.gov/2005/reviews/kinemarpp.pdf)
   // Probability to decay at length x0 or longer = exp(-mass*x0*width/momentum)
   // normalized PDF = mass*width/momentum*exp(-mass*x0*width/momentum)
   // Taking negative log yields -log(constant) + -log(exp(-mass*x0*width/momentum))
   //  = mass*x*width/momentum = mass*x/momentum*lifetime
   //  tauMass = 1.77684, mean lifetime = 290.6e-15
   //  Length must be passed in centimeters!
   //return tauMass*(length / 100.0)/(tauLifetime*momentum);
   double ctau = 8.711e-3; //centimeters
   double factor = tauMass/(ctau*momentum);
   // NLL = (length/100)*factor + (-1)*(log(ctau*momentum) - log(tauMass))
   //return length*factor + log(tauMass) - log(ctau) - log(momentum);
   return length*factor - log(tauMass) + log(ctau) + log(momentum);
}

double nllNuSystemGivenMET(const FourVector& nus, const reco::MET* met)
{
   // ROUGHLY parameterized BY EYE from CMS AN-2009/039 Fig. 25a (For PF Jets)
   double sigma = 6.56 + (8.47/500)*met->sumEt();
   // FIXME FUDGE factor
   sigma = sigma/4.0;
   double sigmaSquared = square(sigma);
   double nllX = square(nus.px()-met->px())/(2*sigmaSquared);
   double nllY = square(nus.py()-met->py())/(2*sigmaSquared);
   return (nllX + nllY + 2*nlGaussianNorm(sigma));
}

// Unsupported reco::Candidate case
template<> double nllVisRapidityGivenMomentum<reco::Candidate>(const reco::Candidate& obj, double rapidity, double momentum)
{
   double mean = 0;
   double sigma = 0;
   mean = 2.238*TMath::Power(momentum, 0.2036);
   sigma = 0.0901 + 0.0006182*momentum;
   double landau = TMath::Landau(rapidity, mean, sigma, true);
   if(landau < 1.0e-10) // sanity check
      landau = 1.0e-10;
   return -1*log(landau);
}

// Different vis rapidity distributions
template<> double nllVisRapidityGivenMomentum<pat::Electron>(const pat::Electron& obj, double rapidity, double momentum)
{
   double mean = 0;
   double sigma = 0;
   mean = 2.238*TMath::Power(momentum, 0.2036);
   sigma = 0.0901 + 0.0006182*momentum;
   double landau = TMath::Landau(rapidity, mean, sigma, true);
   if(landau < 1.0e-10) // sanity check
      landau = 1.0e-10;
   return -1*log(landau);
}

template<> double nllVisRapidityGivenMomentum<pat::Muon>(const pat::Muon& obj, double rapidity, double momentum)
{
   double mean = 0;
   double sigma = 0;
   mean = 2.251*TMath::Power(momentum, 0.2013);
   sigma = 0.09365 + 0.0005191*momentum;
   double landau = TMath::Landau(rapidity, mean, sigma, true);
   if(landau < 1.0e-10) // sanity check
      landau = 1.0e-10;
   return -1*log(landau);
}

// Tau case is special as distribution depends on decay mode
template<> double nllVisRapidityGivenMomentum<pat::Tau>(const pat::Tau& tau, double rapidity, double momentum)
{
   double mean = 0;
   double sigma = 0;
   int decayMode = tau.decayMode();
   switch(decayMode)
   {
      case 0: // 1 prong 0 pi0
         mean = 1.987*TMath::Power(momentum, 0.2215);
         sigma = 0.2215 + 0.0007142*momentum;
         break;
      case 1: // 1 prong 1 pi0
         mean = 2.03*TMath::Power(momentum, 0.2073);
         sigma = 0.1162 + 0.0001587*momentum;
         break;
      case 2: // 1 prong 2 pi0
         mean = 1.864*TMath::Power(momentum, 0.2149);
         sigma = 0.0962 + 0.0001395*momentum;
         break;
      case 10: // 3 prong 0 pi0
         mean = 1.996*TMath::Power(momentum, 0.1987);
         sigma = 0.1429 + 4.695e-5*momentum;
         break;
      case 11: // 3 prong 1 pi0
         mean = 1.868*TMath::Power(momentum, 0.214);
         sigma = 0.1037 + 4.702e-5*momentum;
         break;
      default: // all others
         mean = 1.996*TMath::Power(momentum, 0.1987);
         sigma = 0.1429 + 4.695e-5*momentum;
         break;
   }
   // Single prong no pi0 case is landau distributed
   if (decayMode == 0) {
      return -1*log(TMath::Landau(rapidity, mean, sigma, true));
   }
   else { // otherwise gaussian
      return (square(rapidity-mean)/(2*square(sigma)) + nlGaussianNorm(sigma));
   }
}

/*
template<int decayindex>
double nllVisRapidityGivenMomentum(double rapidity, double momentum)
{
   double mean = 0;
   double sigma = 0;
   switch(decayindex)
   {
      case -2: //muon
         mean = 2.251*TMath::Power(momentum, 0.2013);
         sigma = 0.09365 + 0.0005191*momentum;
         break;
      case -1: //electron
         mean = 2.238*TMath::Power(momentum, 0.2036);
         sigma = 0.0901 + 0.0006182*momentum;
         break;
      case 0: // 1 prong 0 pi0
         mean = 1.987*TMath::Power(momentum, 0.2215);
         sigma = 0.2215 + 0.0007142*momentum;
         break;
      case 1: // 1 prong 1 pi0
         mean = 2.03*TMath::Power(momentum, 0.2073);
         sigma = 0.1162 + 0.0001587*momentum;
         break;
      case 2: // 1 prong 2 pi0
         mean = 1.864*TMath::Power(momentum, 0.2149);
         sigma = 0.0962 + 0.0001395*momentum;
         break;
      case 10: // 3 prong 0 pi0
         mean = 1.996*TMath::Power(momentum, 0.1987);
         sigma = 0.1429 + 4.695e-5*momentum;
         break;
      case 11: // 3 prong 1 pi0
         mean = 1.868*TMath::Power(momentum, 0.214);
         sigma = 0.1037 + 4.702e-5*momentum;
         break;
      case -3: // All others
         mean = 1.996*TMath::Power(momentum, 0.1987);
         sigma = 0.1429 + 4.695e-5*momentum;
         break;
   }

   double result = 0;
   switch (decayindex)
   {
      case -2:
      case -1:
      case 0:
         result = -1*log(TMath::Landau(rapidity, mean, sigma, true));
         break;
      case 1:
      case 2:
      case 10:
      case 11:
      case -3: 
         result = -1*log(TMath::Gaus(rapidity, mean, sigma, true));
         break;
   }

   return result;
}

template double nllVisRapidityGivenMomentum<-3>(double, double);
template double nllVisRapidityGivenMomentum<-2>(double, double);
template double nllVisRapidityGivenMomentum<-1>(double, double);
template double nllVisRapidityGivenMomentum<0>(double, double);
template double nllVisRapidityGivenMomentum<1>(double, double);
template double nllVisRapidityGivenMomentum<2>(double, double);
template double nllVisRapidityGivenMomentum<10>(double, double);
template double nllVisRapidityGivenMomentum<11>(double, double);
*/

FourVectorPair
compInvisibleLeg( const ThreeVector& motherDirection, 
      const FourVector& visibleP4, 
      const double massMother, const double m12Squared, int& error)
{
   /* Compute the two solutions for the missing four vector of a tau decay, given the measured
    * tau direction, visible four vector, and estimated (or exact) mass of the missing energy system.
    *
    * Based on http://arxiv.org/pdf/hep-ph/0607294
    *
    * The following defintions are used:
    *  o p3  : visible momentum
    *  o p12 : neutrino system (with mass m12, leptonic decays m12 > 0, hadronic decays m12 = 0)
    */

   // Visible momentum along tau direction
   const double p3par = motherDirection.Unit().Dot(visibleP4.Vect());
   reco::Candidate::Vector visParallelVect = p3par*motherDirection.Unit();
   // Get the perpindicular compenent
   reco::Candidate::Vector visPerpVect = visibleP4.Vect()-visParallelVect;
   // Visible momentum perpindicular to tau direction
   const double p3perp = sqrt(visPerpVect.Mag2());

   // By construction, the invisible perpindicular momentum must be equal in magnitude 
   // of the visibible perpendicular momentum
   const double p12perp = -1.0*p3perp;
   // Visible energy
   const double e3 = visibleP4.E();

   // p12par is given by (-a +- sqrt(r))/d
   double a = p3par*(square(e3) + m12Squared - square(massMother) + square(p12perp));
   double r = square(e3)*(
           square(square(e3)) 
         + square(m12Squared - square(massMother) + square(p12perp)) 
         + 2*square(p3par)*(m12Squared + square(massMother) + square(p12perp)) 
         + square(square(p3par)) 
         - 2*square(e3)*(m12Squared + square(massMother) + square(p12perp) + square(p3par)) 
         );
   double d = 2*(e3-p3par)*(e3+p3par);

   // Check for unphysical solutions
   if(r < 0) {
      double relativeR = sqrt(fabs(r))/fabs((cube(p3par)-a));
      if(relativeR < 1e-4)
         error = 1;
      r = 0;
   } else { error = 0;}

   // The two solutions
   double factor1 = (cube(p3par) - a + sqrt(r))/d;
   double factor2 = (cube(p3par) - a - sqrt(r))/d;

   reco::Candidate::Vector nuPerpVect = (-1)*visPerpVect;
   reco::Candidate::Vector nuVector1 = factor1*motherDirection.Unit() + nuPerpVect;
   reco::Candidate::Vector nuVector2 = factor2*motherDirection.Unit() + nuPerpVect;

   reco::Candidate::LorentzVector nuFourVector1(nuVector1.x(), nuVector1.y(), nuVector1.z(), sqrt(nuVector1.Mag2() + m12Squared));
   reco::Candidate::LorentzVector nuFourVector2(nuVector2.x(), nuVector2.y(), nuVector2.z(), sqrt(nuVector2.Mag2() + m12Squared));

   //reco::Candidate::LorentzVector total1 = nuFourVector1 + visibleP4;
   //reco::Candidate::LorentzVector total2 = nuFourVector2 + visibleP4;
   //std::cout << "Vis energy: " << visibleP4.energy() << std::endl;
   //std::cout << "Nu1 energy: " << nuFourVector1.energy() << " tot: " << total1.mass() << " N1 eta: " << nuFourVector1.eta() << " N1 phi: " << nuFourVector1.phi()  << std::endl;
   //std::cout << "Nu2 energy: " << nuFourVector2.energy() << " tot: " << total2.mass() << " N2 eta: " << nuFourVector2.eta() << " N2 phi: " << nuFourVector2.phi()  << std::endl;
   return std::make_pair(nuFourVector1, nuFourVector2);
}


}

