#include "TauAnalysis/CandidateTools/interface/svMassRecoLikelihoodAuxFunctions.h"

#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Tau.h"

#include <TMath.h>
#include <TMatrixD.h>
#include <TVectorD.h>
#include <TMatrixDSym.h>
#include <TRotation.h>

#include "math.h"

namespace svMassReco {

   GlobalPoint computeLabSVPosition(const GlobalPoint& fittedPV, double r, double theta, double phi) 
   {
      double x = TMath::Abs(r)*TMath::Cos(phi)*TMath::Sin(theta);
      double y = TMath::Abs(r)*TMath::Sin(phi)*TMath::Sin(theta);
      double z = TMath::Abs(r)*TMath::Cos(theta);
      return fittedPV + GlobalVector(x, y, z);
   }

   double nllPointGiven2DError(const GlobalPoint& point, const GlobalPoint& central, const GlobalError& error)
   {
      AlgebraicVector3 displacement(point.x()-central.x(), point.y()-central.y(), point.z()-central.z());
      //TVector3 displacement2(point.x()-central.x(), point.y()-central.y(), point.z()-central.z());

      //std::cout << "Displacement @ DCA: " << displacement << std::endl;

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

      TMatrixD eigenvectors = matrix.EigenVectors(eigenvalues);
      //td::cout << "Eigenvalues (" << eigenvalues[0] << ", " << eigenvalues[1] << ", " << eigenvalues[2] << ") " << std::endl;

      //std::cout << "Eigenvectors " << std::endl;
      //eigenvectors.Print();

      // Get the lowest eigenvector
      TVector3 lowestEigenvector(eigenvectors[0][2], eigenvectors[1][2], eigenvectors[2][2]);

      /*
      TVector3 middleEigenvector(eigenvectors[0][1], eigenvectors[1][1], eigenvectors[2][1]);
      TVector3 highestEigenvector(eigenvectors[0][0], eigenvectors[1][0], eigenvectors[2][0]);

      TVectorD lowestEigenvectorD(3);
      TVectorD middleEigenvectorD(3);
      TVectorD highestEigenvectorD(3);
      lowestEigenvectorD[0] = lowestEigenvector[0];
      middleEigenvectorD[0] = middleEigenvector[0];
      highestEigenvectorD[0] = highestEigenvector[0];
      lowestEigenvectorD[1] = lowestEigenvector[1];
      middleEigenvectorD[1] = middleEigenvector[1];
      highestEigenvectorD[1] = highestEigenvector[1];
      lowestEigenvectorD[2] = lowestEigenvector[2];
      middleEigenvectorD[2] = middleEigenvector[2];
      highestEigenvectorD[2] = highestEigenvector[2];

      lowestEigenvectorD *= matrix;
      middleEigenvectorD *= matrix;
      highestEigenvectorD *= matrix;
      std::cout << "x M.lowest " << (lowestEigenvectorD)[0] << std::endl;
      std::cout << "x M.middle " << (middleEigenvectorD)[0] << std::endl;
      std::cout << "x M.highest " << (highestEigenvectorD)[0] << std::endl;

      std::cout << "Lowest eigenvector dot displacement " << lowestEigenvector.Dot(displacement2) <<  std::endl;
      std::cout << "Middle eigenvector dot displacement " << middleEigenvector.Dot(displacement2) << std::endl;
      std::cout << "Highest eigenvector dot displacement " << highestEigenvector.Dot(displacement2) << std::endl;
      */

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

      // Transform our error matrix
      AlgebraicMatrix33 rotError = ROOT::Math::Transpose(rotSMatrix)*errorMatrix*rotSMatrix;

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
      if ( !inversionResult || !detOK ) {
         edm::LogWarning("SVMethodLikelihoods") << "2D MATRIX ERROR INVERSION FAILED\n" << errorMatrix;
      }

      double expResult = ROOT::Math::Dot(rotDisp2, rotError2*rotDisp2)/2.0;
      double normResult = nlGaussianNorm(sqrt(determinant), 2);
      //std::cout << " Chi2: " << expResult << " log(Norm): " << normResult << std::endl;
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
      if ( !inversionResult || !detOK ) {
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
      // return tauMass*(length / 100.0)/(tauLifetime*momentum);
      double ctau = 8.711e-3; //centimeters
      double factor = tauMass/(ctau*momentum);
      return length*factor - log(tauMass) + log(ctau) + log(momentum);
   }

   double nllNuSystemGivenMET(const FourVector& nus, const FourVector& direction, const reco::MET* met)
   {
      //std::cout << "<nllNuSystemGivenMET>:" << std::endl;
      //std::cout << " sumEt = " << met->sumEt() << std::endl;
      //std::cout << " metPx = " << met->px() << std::endl;
      //std::cout << " metPy = " << met->py() << std::endl;

      // MET likelihood split into perp/par components along the leptonic leg1.
      //double parSigma =  1.14770e+00 + 3.62242e-02*met->sumEt();
      //double perpSigma = 2.75926e-01 + 3.70582e-02*met->sumEt();
      //double parSigma =  8.068;
      //double perpSigma = 6.905;
      double parSigma = 2.85 + 0.02072*met->sumEt();
      double perpSigma = 2.3 + 0.02284*met->sumEt();

      double parBias = 1.183; // RECO is overestimated
      double perpBias = 0.0;

      //std::cout << " parSigma = " << parSigma << ", parBias = " << parBias << std::endl;
      //std::cout << " perpSigma = " << perpSigma << ", perpBias = " << perpBias << std::endl;

      double output = 0.0;

      double recoMETparToDir = (met->px()*direction.px() + met->py()*direction.py())/direction.pt();
      //std::cout << " recoMETparToDir = " << recoMETparToDir << std::endl;
      double recoMETperpToDir = (met->px()*direction.py() - met->py()*direction.px())/direction.pt();
      //std::cout << " recoMETperpToDir = " << recoMETperpToDir << std::endl;

      double fitMETparToDir = (nus.px()*direction.px() + nus.py()*direction.py())/direction.pt();
      //std::cout << " fitMETparToDir = " << fitMETparToDir << std::endl;
      double fitMETperpToDir = (nus.px()*direction.py() - nus.py()*direction.px())/direction.pt();
      //std::cout << " fitMETperpToDir = " << fitMETperpToDir << std::endl;

      double parResidual = recoMETparToDir - fitMETparToDir - parBias;
      //std::cout << " parResidual = " << parResidual << std::endl;
      double perpResidual = recoMETperpToDir - fitMETperpToDir - perpBias;
      //std::cout << " perpResidual = " << perpResidual << std::endl;

      output += 0.5*square(parResidual/parSigma);
      output += nlGaussianNorm(parSigma);

      output += 0.5*square(perpResidual/perpSigma);
      output += nlGaussianNorm(perpSigma);
      //std::cout << "--> output = " << output << std::endl;
      return output;
   }

   FourVector computeTauMomentum(const ThreeVector& tauDir, const FourVector& labFrameVisible, double theta)
   {
      // Split the lab frame momentum into two components
      double labFramePPerp = labFrameVisible.Vect().Cross(tauDir.Unit()).R();
      double labFramePParallel = labFrameVisible.Vect().Dot(tauDir.Unit());

      // Find the rest frame momentum of the visible stuff
      // pLabPerp = pRestPerp = pRest*sin(theta) ==> pRest = pLabPerp/sin(theta)
      double restFrameP = labFramePPerp/TMath::Sin(theta);
      // The parallel component
      double restFramePParallel = restFrameP*TMath::Cos(theta);
      double restFrameE = TMath::Sqrt(square(labFrameVisible.mass()) + square(restFrameP));

      // Now we can solve for the boost by comparing the visible parallel momentum
      // in the rest and lab frames
      // par_lab = gamma*par_rest + beta*gamma*e_rest
      //         = gamma*par_rest + sqrt(gamma^2 - 1)*e_rest
      // solving for gamma, we get
      // gamma = (e_rest * sqrt( e_rest^2 + par_lab^2 - par_rest^2 ) - par_rest*par_lab)/(e_rest^2 - par_rest^2);

      double gamma = (
            restFrameE * TMath::Sqrt(square(restFrameE) + square(labFramePParallel) - square(restFramePParallel)) -
            restFramePParallel*labFramePParallel ) /
         (square(restFrameE) - square(restFramePParallel));

      // now get the tau energy and momentum
      double tauEnergy = gamma*tauMass;
      double tauMomentum = TMath::Sqrt(square(gamma) - 1)*tauMass;

      // build the four vector
      ThreeVector tauThreeVector  = tauDir.Unit()*tauMomentum;
      FourVector output(tauThreeVector.X(), tauThreeVector.Y(), tauThreeVector.Z(), tauEnergy);

      return output;
   }

   // Determine if two tracks are identically (but from different collections)
   bool tracksAreMatched(const reco::TrackBaseRef& trk1, const reco::TrackBaseRef& trk2) 
   {
      double dEta = trk1->eta() - trk2->eta();
      double dPhi = trk1->phi() - trk2->phi();
      double dR2 = dEta*dEta + dPhi*dPhi;
      if ( dR2 < (0.01*0.01) && trk1->charge() == trk2->charge() ) {
         if ( fabs(trk1->pt() - trk2->pt())/(trk2->pt()+trk1->pt()) < 0.05 ) return true;
      }
      return false;
   }

}

