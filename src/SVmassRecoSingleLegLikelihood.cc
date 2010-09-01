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
   // All values computed assuming isotropic decay
   thetaGuess = TMath::PiOver2();
   thetaError = 0.25*(TMath::Pi()*TMath::Pi() - 8);
   phiGuess = TMath::Pi();
   phiError = square(TMath::Pi())/3.0;

   m12Guess = 0;
   m12Error = 0;

   // Check if we need to guess for M12
   if(legType() == -1) // muon
   {
      m12Guess = 0.921254;
      m12Error = 0.142983;
   }
   if(legType() == -2) // electron
   {
      m12Guess = 0.943999;
      m12Error = 0.15368;
   }

   // Figure what the average decay length is given our guess.  Note that
   // the reconstructed tau energy does not depend on the radius
   
   setPoints(pv, thetaGuess, phiGuess, 1.0, m12Guess);

   double tauMomentum = fittedP4().P();
   const double ctau = 8.711e-3; //centimeters
   radiusGuess = (ctau*tauMomentum)/tauMass;
   radiusError = square(radiusGuess);
}

GlobalVector SVmassRecoSingleLegLikelihood::buildTauDirection(double thetaRest, double phiLab, double m12)
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
   if(isnan(thetaRest) || isnan(phiLab) || isnan(radiusLab) || isnan(m12)) {
      edm::LogWarning("SVSingleLegLikelihood") << "Got NaN value in input! " 
         << "Input    -- M12: " << m12 << " theta: " << thetaRest << " phi: " << phiLab << " radius: " << radiusLab 
         << "\n"
         << "Previous -- M12: " << M12_ << " theta: " << thetaRest_ << " phi: " << phiLab_ << " radius: " << radiusLab_
	 << " NLLkin: " << nllRestFrameKinematics(90, true)
         << " NLLdl: " << nllDecayLength() 
         << " NLLtopo: " << nllTopological();
   }

   // Update the fit parameters
   pv_ = pv;
   thetaRest_ = thetaRest;
   phiLab_ = phiLab;
   radiusLab_ = radiusLab;
   M12_ = (nuSystemIsMassless() ? 0 : m12);


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
/*
double SVmassRecoSingleLegLikelihood::nllOfLeg() const
{ 
   return nllTopological() + nllDecayLength() + nllRestFrameKinematics(); 
}
 */
// Helper functions for the rest frame distributions
namespace restFrameDistributions {
   /// Compute theta from pseudorapidity
   double thetaFromEta(double eta) { return TMath::ACos(TMath::TanH(eta)); }

   /// Compute the rest frame energy and momentum from the neutrino mass
   double visRestEnergyFromM12(double m12, double visMass) {
      return (square(tauMass) + square(visMass) - square(m12))/(2*tauMass);
   }

   double momentumFromEnergy(double energy, double mass) {
      return TMath::Sqrt(square(energy) - square(mass));
   }

   /// Return visible PT given tau PT and rest frame decay parameters
   double F(double taupt, double theta, double m12, double eta, double visMass) {
      using namespace TMath;
      double thetaLab = thetaFromEta(eta);
      double Er = visRestEnergyFromM12(m12, visMass);
      double pr = momentumFromEnergy(Er, visMass);
      return (taupt*(Er + pr*Cos(theta)*Sqrt(1 - (Power(tauMass,2)*Power(Sin(thetaLab),2))/Power(taupt,2))))/tauMass;
   }

   /// Return derivate of F(...) w.r.t. tau PT
   double Fprime(double taupt, double theta, double m12, double eta, double visMass) {
      using namespace TMath;
      double thetaLab = thetaFromEta(eta);
      double Er = visRestEnergyFromM12(m12, visMass);
      double pr = momentumFromEnergy(Er, visMass);
      return (Er + (pr*Cos(theta))/Sqrt(1 - (Power(tauMass,2)*Power(Sin(thetaLab),2))/Power(taupt,2)))/tauMass;
   }

   /// Compute tau PT given visible PT and rest frame decay parameters
   double G(double vispt, double theta, double m12, double eta, double visMass) {
      using namespace TMath;
      double thetaLab = thetaFromEta(eta);
      double Er = visRestEnergyFromM12(m12, visMass);
      double pr = momentumFromEnergy(Er, visMass);
      return (Er*tauMass*vispt - Cos(theta)*Sqrt(Power(tauMass,2)*Power(pr,2)*
               (Power(vispt,2) + (-Power(Er,2) + Power(pr,2)*Power(Cos(theta),2))*
                Power(Sin(thetaLab),2))))/(Power(Er,2) - Power(pr,2)*Power(Cos(theta),2));
   }

   /// Expected distribution of underlying tau lepton PDF.  Currently fitted
   /// from MC
   double tauLeptonPtPDF(double tau_pt) {
      double first_gauss = 2.44e+03*TMath::Gaus(tau_pt, 3.29e+01, 1.37e+01);
      double second_gauss = 1.45e+03*TMath::Gaus(tau_pt, 4.27e+01, 3.32);
      // Add some extra stuff
      second_gauss += 1.45e+03*TMath::Exp(tau_pt/-45.);
      return (first_gauss+second_gauss)/95858.5;
   }

   double smearedKinematicDistribution(double x, double M, double s) {
      // Compute kinematic distribution of decay product PT
      // from a particle of mass <M> at rest, smeared by a Saussian 
      // of width <s>.
      double num_first_term = TMath::Exp(-0.5*square(x)/square(s))*8*s*(
            fourth(M) + 2*square(M)*(2*square(s)+square(x))+
            6*(square(s)+square(x))*(8*square(s)+square(x)));

      double num_second_term = TMath::Exp(-square(M - 2*x)/(8*square(s)))*s*(
            15*fourth(M) + 14*cube(M) + 
            48*(square(s)+square(x))*(8*square(s)+square(x)) +
            4*square(M)*(20*square(s) + 7*square(x)) + 24*M*(7*square(s)*x + cube(x)));

      double num_third_term = 4*TMath::Sqrt(TMath::TwoPi())*x*(
            fourth(M) + 6*square(M)*square(s) + 90*fourth(s) + 
            2*(square(M) + 30*square(s))*square(x) +
            6*fourth(x))*(
            TMath::Erf((M- 2*x)/(2*TMath::Sqrt2()*s)) +
            TMath::Erf(x/(TMath::Sqrt2()*s)));

      double num_factor = 1/(2*TMath::Sqrt(TMath::TwoPi()));
      double numerator = num_factor*(num_first_term - num_second_term + num_third_term);
      // Now make normalization factor
      double den_first_term = (2*TMath::Sqrt(1.0/TMath::PiOver2()) * 
            TMath::Exp(-square(M)/(8*square(s)))*M*s*(
               11*fourth(M) + 44*square(M)*square(s) + 
               240*fourth(s)));
      double den_second_term = TMath::Erf(M/(2*TMath::Sqrt2()*s))*(
            11*fourth(M)*square(M) - 32*fourth(M)*square(s) - 
            96*square(M)*fourth(s) - 960*fourth(s)*square(s));

      double denominator = (1./16)*(den_first_term + den_second_term);
      return numerator/denominator;
   }

   double movingTauLeptonPtPDF(double tauPt, double diTauMass)
   {
      double smearNorm = 0.52 + 0.000658*diTauMass;
      double smearWidth = 1.8 + 0.018*diTauMass;
      double Mfit = 2.3 + 1.04*diTauMass;
      double gammaScale = 6.74 + 0.020*diTauMass;
      double gammaShape = 2.2 + 0.0364*diTauMass;

      return smearNorm*smearedKinematicDistribution(tauPt, Mfit, smearWidth) +
         (1 - smearNorm)*TMath::GammaDist(tauPt, gammaShape, 0, gammaScale);
   }

   double negativeLogMovingVisiblePtPDF(double vispt, double theta, double m12, double eta, double visMass, double diTauMass) {
      double numerator = movingTauLeptonPtPDF(G(vispt, theta, m12, eta, visMass), diTauMass);
      double denominator = 
         Fprime(G(vispt, theta, m12, eta, visMass), theta, m12, eta, visMass);
      return -1.0*TMath::Log(numerator) + TMath::Log(denominator);
   }

   double negativeLogVisiblePtPDF(double vispt, double theta, double m12, double eta, double visMass) {
      double numerator = tauLeptonPtPDF(G(vispt, theta, m12, eta, visMass));
      double denominator = 
         Fprime(G(vispt, theta, m12, eta, visMass), theta, m12, eta, visMass);
      return -1.0*TMath::Log(numerator) + TMath::Log(denominator);
   }
}

double SVmassRecoSingleLegLikelihood::nllRestFrameKinematics(double diTauMass, bool usePtBalanceInFit) const {
   /*
    * The likelihood for a given rest frame decay angle and nu system mass
    * depends on lab frame PT of the visible objects AND the unknown PT spectrum
    * of the underlying tau leptons.  For example, taus from Zs are softer than
    * a tau from a 200 GeV Higgs, so the decays from Zs THAT PASS THE ANALYSIS
    * PT CUTS will be oriented such that the opening angle between the tau and
    * decay products is smaller.  The formulation below assumes that the
    * underlying spectrum comes from Z decays.  For a given visible PT, the
    * likelihood for theta and Mnunu to be compatible with he underlying
    * distributions is computed.
    *
    * By Bayes theorem:
    *
    * P(theta, m12|pt_vis, eta_vis) = P(pt_vis|theta, m12, eta_vis)*P(theta,m12)
    *
    * P(pt_vis|theta, m12) is computed by relating pt_vis to pt_tau, theta, m12,
    * and eta_vis by the function F.  Eta_vis and Eta_tau are assumed to be
    * equivalent.
    *
    * pt_vis = F(pt_tau, theta, m12, eta_vis)
    *
    * The inverse of this function is G.
    * pt_tau = G(pt_vis, theta, m12, eta_vis)
    *
    * The derivative of F with respect to pt_tau is written as F'
    *
    * We can transform the ansatz tau lepton PT sepctrum P_tau(pt_tau) into the
    * visible spectrum according to:
    *
    * P(pt_vis|theta, m12, eta_vis) = 
    *     P_tau( G(pt_vis, ...) ) * [F'( G(pt_vis, ...) )]^(-1)
    */

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

   //output += restFrameDistributions::negativeLogMovingVisiblePtPDF(
   //      visP4_.pt(), thetaRest_, M12_, visP4_.eta(), visibleMass_, diTauMass);
   if ( usePtBalanceInFit )
     output += -log(restFrameDistributions::movingTauLeptonPtPDF(p4_.pt(), diTauMass));

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
   // If we fit a reco. SV, see how compatible our fitted SV is with it
   double output = 0;
   if(vertex_.isValid()) {
      output = nllPointGivenVertex(sv_, vertex_);
      if(isnan(output)) {
         edm::LogWarning("SVSingleLeg") << " Got nan/inf for nllTopo vertex constraint! SV:" << sv_;
      }
   } else 
   {
      TrajectoryStateClosestToPoint tcsp = bestTrack_.trajectoryStateClosestToPoint(sv_);
      if(!tcsp.isValid())
      {
         edm::LogWarning("SVSingleLeg") << "Can't find DCA for point " << sv_ << " returning 30."; 
         return 30;
      } else {
         // we only have one track, see how compatable our point is to the track
         output = nllPointGivenTrack(tcsp);
         if(isnan(output)) {
            edm::LogWarning("SVSingleLeg") << " Got nan for nllTopo track constraint! SV:" << sv_ 
               << " PV: " << pv_ << " radius: " << radiusLab_;
         }
      }
   }
   return output;
}

void SVmassRecoSingleLegLikelihood::printTo(std::ostream &out) const
{
   using namespace std;
   out << "Type: " << legType() << endl;
   //out << "NLL" << setw(10) << nllOfLeg() << endl;
   out << "NLL" << setw(10) << (nllTopological() + nllDecayLength() + nllRestFrameKinematics(90, true)) << endl;
   out << "- NLLTopo" << setw(10) << nllTopological() << endl;
   out << "- NLLDecay" << setw(10) << nllDecayLength() << setw(10) << endl;
   out << "- NLLRestFrame" << setw(10) << nllRestFrameKinematics(90, true)  << endl;
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
