#include "TauAnalysis/CandidateTools/interface/SVmassRecoSingleLegLikelihood.h"
#include "TauAnalysis/CandidateTools/interface/svMassRecoLikelihoodAuxFunctions.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "DataFormats/GeometrySurface/interface/Line.h"

namespace svMassReco { 

SVmassRecoSingleLegLikelihood::SVmassRecoSingleLegLikelihood(const SVmassRecoSingleLegExtractorBase& extractor, 
      const std::vector<reco::TransientTrack>& tracks)
   : extractor_(extractor),
   tracks_(tracks),
   propagator_(tracks[0].isValid() ? tracks[0].field() : NULL),
   visP4_(extractor.p4()),
   visP3_(extractor.p4().Vect()),
   visP3GlobalVector_(extractor.p4().px(), extractor.p4().py(), extractor.p4().pz()),
   visPDirection_(TVector3(extractor.p4().px(), extractor.p4().py(), extractor.p4().pz()).Unit()),
   visibleMass_(extractor.p4().mass())
{
   edm::LogInfo("SingleLegLikelihood") << "Initializing...";

   // Check if we have enough tracks to fit a vertex
   if(tracks_.size() > 2) {
      KalmanVertexFitter vertexFitter(true);
      edm::LogInfo("SingleLegLikelihood") << "Fitting multi-track leg vertex";
      // fit the vertex
      vertex_ = vertexFitter.vertex(tracks_);
   }

   // Now select the best track to use for finding the initial point
   double bestNormedChi2 = -1;
   for(size_t iTrack = 0; iTrack < tracks_.size(); ++iTrack)
   {
      if(!tracks_[iTrack].isValid())
      {
         edm::LogError("SVSingleLegLikelihood") << "The " << iTrack << " track is invalid";
      }
      double chi2 = tracks_[iTrack].track().normalizedChi2();
      double pt = tracks_[iTrack].track().pt();
      edm::LogInfo("SingleLegLikelihood") << "Track " << iTrack << " has chi2 " << chi2 << " and pt " << pt;
      if(bestNormedChi2 < 0 || chi2 < bestNormedChi2)
      {
         bestTrack_ = tracks_[iTrack];
         bestNormedChi2 = chi2;
      }
   }
   edm::LogInfo("SingleLegLikelihood") << "Found best track with chi2 " << bestNormedChi2;
   if(!bestTrack_.isValid()) 
   {
      edm::LogError("SVSingleLegLikelihood") << "No valid track found";
   }
}

double SVmassRecoSingleLegLikelihood::findPGivenM12(double M12) const
{
   const double tauMassSq = tauMass*tauMass;
   return TMath::Sqrt(
         (tauMassSq - square(M12 + visibleMass_))*
         (tauMassSq - square(M12 - visibleMass_)))/(2*tauMass);
}

double SVmassRecoSingleLegLikelihood::maxM12() const 
{
   return (tauMass - visibleMass_);
}

void SVmassRecoSingleLegLikelihood::findIntialValues(const GlobalPoint& pv,
      double &thetaGuess, double &thetaError, 
      double &phiGuess, double &phiError,
      double &radiusGuess, double &radiusError,
      double &m12Guess, double &m12Error)
{
   size_t thetaSteps = 10;
   size_t phiSteps = 10;

   // only need to step m12 if it is non zero
   size_t m12Steps = (nuSystemIsMassless()) ? 1 : 5;

   // get the valid range (away from corners)
   const double thetaMin = 1e-5;
   const double thetaMax = TMath::Pi() - 1e-5;
   const double phiMin = -TMath::Pi();
   const double phiMax = TMath::Pi();
   const double m12Min = 1e-5;
   const double m12Max = maxM12() - 1e-5;

   // Set errors for the parameters we are optimizing
   thetaError = fabs(thetaMax - thetaMin)/thetaSteps;
   phiError = fabs(phiMax - phiMin)/phiSteps;
   m12Error = (m12Max-m12Min)/m12Steps;

   double bestLikelihood = 1e30;
   double bestTheta = 0;
   double bestPhi = 0;
   double bestRadius = 0;
   double bestRadiusError = 0;
   double bestM12 = 0;
     
   // loop over theta steps
   for(size_t iTheta = 0; iTheta < thetaSteps; ++iTheta)
   {
      double theta = thetaMin + iTheta*(thetaMax - thetaMin)*1.0/(thetaSteps-1);
      for(size_t iPhi = 0; iPhi < phiSteps; ++iPhi)
      {
         // don't need full (steps - 1) range as it's periodic on this interval
         double phi = phiMin + iPhi*(phiMax-phiMin)*1.0/(phiSteps);
         // loop over m12 steps (trivial in massless case)
         for(size_t iM12 = 0; iM12 < m12Steps; ++iM12)
         {
            // keep away from corners
            double m12 = (nuSystemIsMassless()) ? 0.0 : m12Min + (m12Max-m12Min)*iM12*1.0/(m12Steps-1);
            // Now we need to figure out what the radius of our test point
            // is.  We're looping over the directions, but let's find the radius
            // that gets us closest to the best track
            // find theta in lab frame
            GlobalVector direction = buildTauDirection(theta, phi, m12);
            GlobalPoint nonConstOrigin = pv; // wtf
            Line line(nonConstOrigin, direction);
            TrajectoryStateOnSurface tsos = propagator_.extrapolate(bestTrack_.initialFreeState(), line);
            double radius = 1;
            if(!tsos.isValid()) {
               edm::LogWarning("SVSingleLegLikelihood") << "Can't find a valid tsos for intial point";
            } else 
            {
               GlobalPoint pca = tsos.globalPosition();
               GlobalPoint displacement(pca.x()-pv.x(), pca.y()-pv.y(), pca.z()-pv.z());
               radius = displacement.mag();
            }
            //double err_radius = pca.mag(); // not quite right..
            double err_radius = 0.1;
            // Setup our guess
            setPoints(pv, theta, phi, radius, m12);
            // See if it is the best so far
            if(nllOfLeg() < bestLikelihood)
            {
               bestTheta = theta;
               bestPhi = phi;
               bestRadius = radius;
               bestRadiusError = err_radius;
               bestM12 = m12;
               bestLikelihood = nllOfLeg();
            }
         }
      }
   }

   // return output
   thetaGuess = bestTheta;
   m12Guess = bestM12;
   radiusGuess = bestRadius;
   phiGuess = bestPhi;
   radiusError = bestRadiusError;
   setPoints(pv, thetaGuess, phiGuess, radiusGuess, m12Guess);
   edm::LogInfo("SingleLegLikelihoodInitial") << "Found initial condition @" << *this;
}

GlobalVector SVmassRecoSingleLegLikelihood::buildTauDirection(double thetaRest,
      double phiLab, double m12)
{
   // Determine lab frame opening angle between tau direction
   // and vis. momentum
   //
   // pl_perp = pr(m12)*sin(theta_r) ==>
   // pl*sin(thetal) = pr(m12)*sin(theta_r)
   // thetal = asin(pr(m12)*sin(theta_r)/pl)

   double thetaLab = 
         TMath::ASin(findPGivenM12(m12)*TMath::Sin(thetaRest)/visP4_.P());

   // Build our displacement vector assuming visP parallel to Z axis
   TVector3 secondaryVertex(
         TMath::Sin(thetaLab)*TMath::Cos(phiLab),
         TMath::Sin(thetaLab)*TMath::Sin(phiLab),
         TMath::Cos(thetaLab));
   // Rotate such that Z is along visible momentum
   secondaryVertex.RotateUz(visPDirection_);
   return GlobalVector(secondaryVertex.x(), secondaryVertex.y(), secondaryVertex.z());
}

void SVmassRecoSingleLegLikelihood::setPoints(const GlobalPoint& pv, double thetaRest, 
      double phiLab, double radiusLab, double m12)
{
   pv_ = pv;
   thetaRest_ = thetaRest;
   phiLab_ = phiLab;
   radiusLab_ = radiusLab;
   M12_ = m12;

   // Find the tau direction, given the fit parameters
   GlobalVector secondaryVertex = buildTauDirection(thetaRest, phiLab, m12);

   // Scale 
   secondaryVertex *= radiusLab;

   // Set SV point
   sv_ = pv;
   sv_ += secondaryVertex;

   // Set tau direction
   legDir_ = ThreeVector(secondaryVertex.x(), secondaryVertex.y(), secondaryVertex.z());
   // Compute tau momentum
   p4_ =  computeTauMomentum(legDir_, visP4_, thetaRest);
   // By defintion, nu p4 is tauP4 - visP4
   nuP4_ = p4_ - visP4_;
}

double SVmassRecoSingleLegLikelihood::nllOfLeg() const
{ 
   return nllTopological() + nllDecayLength() + nllRestFrameKinematics(); 
}

double SVmassRecoSingleLegLikelihood::nllRestFrameKinematics() const 
{
   // in all cases we get a constraint from the solid angle differential term
   double output = -1.0*log(TMath::Sin(thetaRest_));
   if(!nuSystemIsMassless())
   {
      // In a three body system, there is an additional constraint on
      // the mass of the two neutrino sysemts
      // see PDG kinematics summary, equations 38.20(a,b) 
      double logOfM12Terms = 
         log(M12_) - log(2) + // p1 term
         0.5*(log(square(tauMass) - square(M12_ + visibleMass_)) +   //p3 term numerator
               log(square(tauMass) - square(M12_ - visibleMass_))) +  //p3 term numerator
         -1*(log(tauMass) + log(2));

      // return sum of m12 terms and sine log term
      output -= logOfM12Terms;
   }

   if(isnan(output) || isinf(output)) {
      //edm::LogWarning("SVSingleLeg") << " Got nan/inf for nllRestFrame!  nuSystemIsMassless: " << nuSystemIsMassless() << " thetaRest: " << thetaRest_ << " M12: " << M12_ << " Returning 30 instead.";
      return 30;
   }

   return output;
}

double SVmassRecoSingleLegLikelihood::nllDecayLength() const {
   return nllTauDecayLengthGivenMomentum(legDir_.r(), p4_.P());
}

double SVmassRecoSingleLegLikelihood::nllTopological() const
{
   // If we fit a SV, see how compatabile our fitted SV is with it
   if(vertex_.isValid()) {
      return nllPointGivenVertex(sv_, vertex_);
   } else 
   {
      TrajectoryStateClosestToPoint tcsp = bestTrack_.trajectoryStateClosestToPoint(sv_);
      if(!tcsp.isValid())
      {
         edm::LogWarning("SVSingleLeg") << "Can't find DCA for point " << sv_ << " returning 30."; 
         return 30;
      }
      // we only have one track, see how compatable our point is to the track
      return nllPointGivenTrack(tcsp);
   }

}

void SVmassRecoSingleLegLikelihood::printTo(std::ostream &out) const
{
   using namespace std;
   out << "Type: " << legType() << endl;
   out << "NLL" << setw(10) << nllOfLeg() << endl;
   out << "- NLLTopo" << setw(10) << nllTopological() << endl;
   out << "- NLLDecay" << setw(10) << nllDecayLength() << setw(10) << endl;
   out << "- NLLRestFrame" << setw(10) << nllRestFrameKinematics()  << endl;
   out << "-- SV" << setw(30) << sv_ << endl;
   out << "-- PV" << setw(30) << pv_ << endl;
   out << "-- ThetaRest" << setw(30) << thetaRest_ << endl;
   out << "-- PhiLab" << setw(30) << phiLab_ << endl;
   out << "-- RadiusLab" << setw(30) << radiusLab_ << endl;
   out << "-- M12 " << setw(30) << M12_ << endl;
   out << "-- Dir" << setw(10) << legDir_ << endl;
   out << "-- TauP4" << setw(30) << p4_ << " Mass: " << p4_.mass() << endl;
   out << "-- VisP4" << setw(30) << visP4_ << " Mass: " << visP4_.mass() << endl;
   out << "-- NuP4" << setw(30) << nuP4_ << endl;
}

}
