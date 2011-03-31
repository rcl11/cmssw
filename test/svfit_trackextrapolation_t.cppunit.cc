/*
 * Tests for the linearized track facilities used in the SVfit
 *
 * To run: $CMSSW_BASE/test/$SCRAM_ARCH/TestSVfitTrackExtrapolation
 *
 * Author: Evan K. Friis (UC Davis)
 *
 */

#include <cppunit/extensions/HelperMacros.h>
#include <sstream>
#include <Utilities/Testing/interface/CppUnit_testdriver.icpp>

#include "TMath.h"
#include "TauAnalysis/CandidateTools/interface/SVfitTrackExtrapolation.h"

using namespace SVfit::track;

class testSVFitTrackExtrapolation : public CppUnit::TestFixture {
  CPPUNIT_TEST_SUITE(testSVFitTrackExtrapolation);
  CPPUNIT_TEST(testLogLikelihoodFromDisplacementNoRotation);
  CPPUNIT_TEST(testLogLikelihoodFromDisplacement);
  CPPUNIT_TEST(testLogLikelihoodNoDisplacement);
  CPPUNIT_TEST(testSVDisplacementEquivalence);
  CPPUNIT_TEST(testApproximateTrackError);
  CPPUNIT_TEST(testDisplacementCalculation);
  CPPUNIT_TEST_SUITE_END();

  public:
    void setUp() {
      origin_ = GlobalPoint(0,0,0);
      // The fake point about which the DCA is "computed"
      refPoint_ = GlobalPoint(-1, 0, 0);
      zAxis_ = GlobalVector(0,0,1);
      yAxis_ = GlobalVector(0,1,0);
      xAxis_ = GlobalVector(1,0,0);
      // A standard offset
      stdOffset_ = GlobalPoint(1, 1, 1);

      // Make a fake error matrix about z.  No correlations.
      // Sigma_x = 1, Sigma_y = 2, Sigma_z = 0
      double fakeErrorAboutZVals[9] = {
          1, 0, 0,
          0, 2, 0,
          0, 0, 0};
      AlgebraicMatrix33 fakeErrorMatrixAboutZ(fakeErrorAboutZVals, 9);

      // Make a fake error matrix about x.  No correlations.
      // Sigma_x = 0, Sigma_y = 2, Sigma_z = 1
      double fakeErrorAboutXVals[9] = {
          0, 0, 0,
          0, 2, 0,
          0, 0, 1};
      AlgebraicMatrix33 fakeErrorMatrixAboutX(fakeErrorAboutXVals, 9);

      atOriginAlongXAxis_ = TrackExtrapolation(
          refPoint_, origin_, xAxis_, fakeErrorMatrixAboutX);
      atOriginAlongZAxis_ = TrackExtrapolation(
          refPoint_, origin_, zAxis_, fakeErrorMatrixAboutZ);

      atOffsetAlongXAxis_ = TrackExtrapolation(
          refPoint_, stdOffset_, xAxis_, fakeErrorMatrixAboutX);
      atOffsetAlongZAxis_ = TrackExtrapolation(
          refPoint_, stdOffset_, zAxis_, fakeErrorMatrixAboutZ);
    }

    void testLogLikelihoodFromDisplacementNoRotation() {
      // Test in non-rotated frame
      AlgebraicVector3 testDisplacement(1, 0, 0);
      double result = atOffsetAlongZAxis_.logLikelihoodFromDisplacement(testDisplacement);
      // We expect to be 1 sigma away in X, zero sigma away in Y
      double determinant = 2;
      // From http://en.wikipedia.org/wiki/Multivariate_normal_distribution
      // The matrix inverse is ((1,0), (0, 0.5))
      // So the exponent term is just 1, since we only go along x
      double logLikelihood =
        -TMath::Log(TMath::TwoPi())
        - (0.5)*TMath::Log(determinant) +
        -(0.5)*1;
      CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE(
          "Testing loglikelihood from displacement along X",
          logLikelihood, result, 1e-8);

      // Now try against displacement in Y which has sigma = 2
      testDisplacement = AlgebraicVector3(0,1,0);
      result = atOffsetAlongZAxis_.logLikelihoodFromDisplacement(testDisplacement);

      logLikelihood =
        -TMath::Log(TMath::TwoPi())
        - (0.5)*TMath::Log(determinant) +
        - (0.5)*0.5;

      CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE(
          "Testing loglikelihood from displacement along Y",
          logLikelihood, result, 1e-8);
    }

    void testLogLikelihoodNoDisplacement() {
      AlgebraicVector3 nullDisplacement(0, 0, 0);
      double result = atOffsetAlongZAxis_.logLikelihoodFromDisplacement(nullDisplacement);
      // We expect to be 1 sigma away in X, zero sigma away in Y
      double determinant = 2;
      double logLikelihood =
        -TMath::Log(TMath::TwoPi())
        - (0.5)*TMath::Log(determinant);
      CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE(
          "Testing loglikelihood from null displacement",
          logLikelihood, result, 1e-8);
    }

    // Test the case where the extraplator has to rotate the coordinates
    // internally. (it always assumes track lies along z)
    void testLogLikelihoodFromDisplacement() {
      // Test in non-rotated frame
      AlgebraicVector3 testDisplacement(0, 0, 1);
      double result = atOffsetAlongXAxis_.logLikelihoodFromDisplacement(testDisplacement);
      // We expect to be 1 sigma away in X, zero sigma away in Y
      double determinant = 2;
      // From http://en.wikipedia.org/wiki/Multivariate_normal_distribution
      // The matrix inverse is ((1,0), (0, 0.5))
      // So the exponent term is just 1, since we only go along x
      double logLikelihood =
        -TMath::Log(TMath::TwoPi())
        - (0.5)*TMath::Log(determinant) +
        - (0.5)*1;
      CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE(
          "Testing loglikelihood from displacement along Z",
          logLikelihood, result, 1e-8);

      // Now try against displacement in Y which has sigma = 2
      testDisplacement = AlgebraicVector3(0,1,0);
      result = atOffsetAlongZAxis_.logLikelihoodFromDisplacement(testDisplacement);
      logLikelihood =
        -TMath::Log(TMath::TwoPi())
        - (0.5)*TMath::Log(determinant) +
        - (0.5)*0.5;
      CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE(
          "Testing loglikelihood from displacement along Y",
          logLikelihood, result, 1e-8);
    }

    void testSVDisplacementEquivalence() {
      AlgebraicVector3 testDisplacement(1, 0, 0);
      // should have same displacment, as track is along Z
      AlgebraicVector3 testSV(1, 0, 20);
      double displacementResult =
        atOriginAlongZAxis_.logLikelihoodFromDisplacement(testDisplacement);
      double svResult =
        atOriginAlongZAxis_.logLikelihood(testSV);
      CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE(
          "Testing discplacment-SV consistency",
          displacementResult, svResult, 1e-8);

      // Make sure displacement function does not depend on origin of track.
      double displacementResultAtOffset =
        atOffsetAlongZAxis_.logLikelihoodFromDisplacement(testDisplacement);
      CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE(
          "Testing discplacment-displacement(offset) consistency",
          displacementResult, displacementResultAtOffset, 1e-8);
    }

    void testApproximateTrackError() {
      // Test along X
      double fakeErrorAboutZVals[9] = {
          3, 2, 0,
          2, 5, 0,
          0, 0, 0};
      AlgebraicMatrix33 fakeErrorMatrixAboutZ(fakeErrorAboutZVals, 9);
      TrackExtrapolation testTrack(
          GlobalPoint(-1, 0, 0), // ref point
          origin_, // track dca point
          zAxis_,
          fakeErrorMatrixAboutZ);
      // Sigma^2 is 3 in x direction, the displacement should be pure X
      CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE(
          "Testing approximate error on prinicipal axis",
          testTrack.approximateTrackError(),
          TMath::Sqrt(3), 1e-8);

      TrackExtrapolation testTrackDiagonalOffset(
          GlobalPoint(-1, -1, 0), // ref point
          origin_, // track dca point
          zAxis_,
          fakeErrorMatrixAboutZ);
      CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE(
          "Testing approximate error on non-principal axis",
          testTrackDiagonalOffset.approximateTrackError(),
          TMath::Sqrt(12)/TMath::Sqrt(2), 1e-8);
    }

    void testDisplacementCalculation() {
      AlgebraicVector3 testSV(1, 1, 200);
      AlgebraicVector3 displacement =
        atOriginAlongZAxis_.displacementFromTrack(testSV);
      // We expect displacement to be @ (1, 1, 0);
      CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE(
          "Testing x displacment", displacement(0), 1, 1e-6);
      CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE(
          "Testing y displacment", displacement(1), 1, 1e-6);
      CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE(
          "Testing z displacment", displacement(2), 0, 1e-6);
    }

  private:
    GlobalPoint origin_;
    GlobalPoint refPoint_;
    GlobalPoint stdOffset_;
    GlobalVector zAxis_;
    GlobalVector yAxis_;
    GlobalVector xAxis_;
    TrackExtrapolation atOriginAlongXAxis_;
    //TrackExtrapolation atOriginAlongYAxis_;
    TrackExtrapolation atOriginAlongZAxis_;
    TrackExtrapolation atOffsetAlongXAxis_;
    //TrackExtrapolation atOffsetAlongYAxis_;
    TrackExtrapolation atOffsetAlongZAxis_;
};

CPPUNIT_TEST_SUITE_REGISTRATION(testSVFitTrackExtrapolation);
