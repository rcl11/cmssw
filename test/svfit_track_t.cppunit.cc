/*
 * Closure tests track geometery extractions used in the
 * SV fit algorithm.
 *
 * To run: $CMSSW_BASE/test/$SCRAM_ARCH/TestSVTrack
 *
 * Author: Evan K. Friis (UC Davis)
 *
 */

#include <Utilities/Testing/interface/CppUnit_testdriver.icpp>
#include <cppunit/extensions/HelperMacros.h>

#include "TauAnalysis/CandidateTools/interface/SVfitTrackTools.h"
#include "TMath.h"
#include <sstream>

using namespace SVfit::track;

namespace {
  template<class T>
  std::string dump(const T& pt) {
    std::stringstream out;
    out << pt << std::endl;
    return out.str();
  }
}

class testSVFitTrack : public CppUnit::TestFixture {
  CPPUNIT_TEST_SUITE(testSVFitTrack);
  CPPUNIT_TEST(testIntersection);
  CPPUNIT_TEST(testIntersectionWithLineOriginInsideCone);
  CPPUNIT_TEST(testMissingPCA);
  CPPUNIT_TEST(testPCAtoPoint);
  CPPUNIT_TEST(testTransformClosure);
  CPPUNIT_TEST(testTranslationTransform);
  CPPUNIT_TEST(testRotationTransform);
  CPPUNIT_TEST(testPCA);
  CPPUNIT_TEST(testDoubleIntersection);
  CPPUNIT_TEST(testPCAPointBehindPV);
  CPPUNIT_TEST(testLineOriginInsideConeCheck);
  CPPUNIT_TEST(testPropagateLineOriginToConeVertexPlane);
  CPPUNIT_TEST(testLineConeCollinear);
  CPPUNIT_TEST(testLinePropagation);
  CPPUNIT_TEST(testTrackCorrections);

  CPPUNIT_TEST_SUITE_END();
  public:
     void setUp() {
       origin_ = GlobalPoint(0,0,0);
       zAxis_ = GlobalVector(0,0,1);
       yAxis_ = GlobalVector(0,1,0);
       xAxis_ = GlobalVector(1,0,0);
       pi_ = 3.14159265;
     }
     // Check for a valid intersection
     void testIntersection() {
       double angle = pi_/5;
       // Cone is at origin and points in Z
       // Line is offset 1 unit in X and points in Z
       GlobalPoint coneOrigin = origin_;
       GlobalVector coneDirection = zAxis_;
       GlobalPoint lineOrigin = GlobalPoint(1, 0, 0);
       GlobalVector lineDirection = zAxis_;
       int status = -999;
       GlobalPoint result = intersectionOfLineAndCone(
           lineOrigin, lineDirection,
           coneOrigin, coneDirection, angle, status);
       CPPUNIT_ASSERT(status == 1);
       CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE(
           "Checking x position", result.x(), 1.0, 1e-6);
       CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE(
           "Checking y position", result.y(), 0.0, 1e-6);
       CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE(
           "Checking z position", TMath::Tan(angle)*result.z(),
           1.0, 1e-6);
     }

     // Check the case where the line origin point is inside the cone.
     // (The intersection finder should internally move it outside)
     void testIntersectionWithLineOriginInsideCone() {
       double angle = pi_/5;
       // Cone is at origin and points in Z
       // Line is offset 1 unit in X and points in Z
       GlobalPoint coneOrigin = origin_;
       GlobalVector coneDirection = zAxis_;
       GlobalPoint lineOrigin = GlobalPoint(1, 0, 20); // <-- here is the diff.
       GlobalVector lineDirection = zAxis_;
       int status = -999;
       CPPUNIT_ASSERT_MESSAGE(
           "Sanity check: Is inside cone",
           pointIsInsideCone(lineOrigin, coneOrigin, coneDirection, angle));
       GlobalPoint result = intersectionOfLineAndCone(
           lineOrigin, lineDirection,
           coneOrigin, coneDirection, angle, status);
       CPPUNIT_ASSERT(status == 1);
       CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE(
           "Checking x position", result.x(), 1.0, 1e-6);
       CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE(
           "Checking y position", result.y(), 0.0, 1e-6);
       CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE(
           "Checking z position", TMath::Tan(angle)*result.z(),
           1.0, 1e-6);
     }

     void testPCA() {
       double angle = pi_/5;
       GlobalPoint coneOrigin = origin_;
       GlobalVector coneDirection = zAxis_;
       // Line is offset outside the cone and orthogonal to cone axis
       GlobalPoint lineOrigin = GlobalPoint(3, 0, 3);
       // The line has to have some Z component for the equations to find
       // solutions.
       GlobalVector lineDirection = GlobalVector(0, 1, 0.000000001);
       int status = -999;
       // Make sure we don't find an intersection (shouldn't exist)
       GlobalPoint badResult = intersectionOfLineAndCone(
           lineOrigin, lineDirection,
           coneOrigin, coneDirection, angle, status);
       CPPUNIT_ASSERT_MESSAGE(
           "Intersection exists when it shouldn't!", status == 0);
       status = -999;
       // Find the point *one the line* closest to the cone.
       GlobalPoint goodResult = pcaOfLineToCone(
           lineOrigin, lineDirection,
           coneOrigin, coneDirection, angle, status);
       CPPUNIT_ASSERT_MESSAGE(
           "Can't find the line PCA to the cone!", status == 1);
       // By symmetry the PCA of line to cone should be y = 0
       CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE(
           "Checking line PCA x " + dump(goodResult), goodResult.x(), 3, 1e-6);
       CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE(
           "Checking line PCA y", goodResult.y(), 0, 1e-6);
       CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE(
           "Checking line PCA z", goodResult.z(), 3, 1e-6);

       // Check the point we find on the cone.
       GlobalPoint coneResult = pcaOfConeToPoint(
           goodResult, coneOrigin, coneDirection, angle, status);
       CPPUNIT_ASSERT_MESSAGE(
           "Can't find the cone PCA to the line PCA point", status == 1);
       CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE(
           "Checking line PCA y", goodResult.y(), 0, 1e-6);
       // Expected radius of cone PCA point
       double conePCARadius = lineOrigin.mag()*TMath::Cos(
           lineOrigin.theta() - angle);
       double conePCA_x = conePCARadius*TMath::Sin(angle);
       double conePCA_z = conePCARadius*TMath::Cos(angle);
       CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE(
           "Checking line PCA y", coneResult.y(), 0, 1e-6);
       CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE(
           "Checking line PCA z", coneResult.z(), conePCA_z, 1e-6);
       CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE(
           "Checking line PCA x", coneResult.x(), conePCA_x, 1e-6);
     }

     /// The the ability to correctly find the point of closest approach
     /// on the cone to an arbritrary point
     void testPCAtoPoint() {
       double angle = pi_/5;
       GlobalPoint coneOrigin = origin_;
       GlobalVector coneDirection = zAxis_;
       // Line is offset outside the cone and orthogonal to cone axis
       GlobalPoint point = GlobalPoint(3, 0, 3);
       int status = -999;
       GlobalPoint coneResult = pcaOfConeToPoint(
           point, coneOrigin, coneDirection, angle, status);
       CPPUNIT_ASSERT(status == 1);
       // Expected radius of cone PCA point
       double conePCARadius = point.mag()*TMath::Cos(
           point.theta() - angle);
       double conePCA_x = conePCARadius*TMath::Sin(angle);
       double conePCA_z = conePCARadius*TMath::Cos(angle);
       CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE(
           "Checking cone PCA y", coneResult.y(), 0, 1e-6);
       CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE(
           "Checking cone PCA z", coneResult.z(), conePCA_z, 1e-6);
       CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE(
           "Checking cone PCA x", coneResult.x(), conePCA_x, 1e-6);

     }

     void testMissingPCA() {
       double angle = pi_/5;
       GlobalPoint coneOrigin = origin_;
       GlobalVector coneDirection = zAxis_;
       // Line is offset outside the cone and it's angle is larger than the
       // opening angle of the cone.  The only intersection can be behind
       // the primary vertex, which is not allowed.
       GlobalPoint lineOrigin = GlobalPoint(3, 0, 0);
       double lineAngle = angle + 0.1;
       GlobalVector lineDirection = GlobalVector(
           TMath::Cos(lineAngle),
           0,
           TMath::Sin(lineAngle));
       int status = -999;
       // Make sure we don't find an intersection (shouldn't exist)
       GlobalPoint badResult = intersectionOfLineAndCone(
           lineOrigin, lineDirection,
           coneOrigin, coneDirection, angle, status);
       CPPUNIT_ASSERT_MESSAGE(
           "Intersection exists when it shouldn't!", status == 0);
       status = -999;
       // Find the point *one the line* closest to the cone.
       GlobalPoint badResultPCA = pcaOfLineToCone(
           lineOrigin, lineDirection,
           coneOrigin, coneDirection, angle, status);
       CPPUNIT_ASSERT_MESSAGE(
           "PCA exists when it shouldn't!", status == 0);
     }

     // Test internal transformation functions used by the track tools
     void testTransformClosure() {
       GlobalPoint crazyOrigin(56, 100, 2);
       GlobalVector crazyUz(200, 0, 150);
       GlobalPoint ptToTransform(1, 0, 1);
       GlobalVector vecToTransform(1, 0, 1);

       // Test inversion
       GlobalPoint crazyPt = transform(crazyOrigin, crazyUz, ptToTransform);
       GlobalVector crazyVec = transform(crazyOrigin, crazyUz, vecToTransform);

       GlobalPoint sanePt = untransform(crazyOrigin, crazyUz, crazyPt);
       GlobalVector saneVec = untransform(crazyOrigin, crazyUz, crazyVec);

       // Test we are back where we started
       CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE(
           "TestPtTransform x", sanePt.x(), 1, 1e-6);
       CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE(
           "TestPtTransform y", sanePt.y(), 0, 1e-6);
       CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE(
           "TestPtTransform z", sanePt.z(), 1, 1e-5);

       // Test we are back where we started
       CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE(
           "TestVecTransform x", saneVec.x(), 1, 1e-6);
       CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE(
           "TestVecTransform y", saneVec.y(), 0, 1e-6);
       CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE(
           "TestVecTransform z", saneVec.z(), 1, 1e-6);
     }

     void testTranslationTransform() {
       GlobalPoint newOrigin(10, 0, 0);
       // Don't change
       GlobalVector newUz(0, 0, 1);

       GlobalPoint ptToTransform(1, 0, 1);
       GlobalVector vecToTransform(1, 0, 1);

       GlobalPoint newPt = transform(newOrigin, newUz, ptToTransform);
       GlobalVector newVec = transform(newOrigin, newUz, vecToTransform);

       CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE(
           "TestPtTranslate_x", newPt.x(), -9, 1e-6);
       CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE(
           "TestPtTranslate_y", newPt.y(), ptToTransform.y(), 1e-6);
       CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE(
           "TestPtTranslate_z", newPt.z(), ptToTransform.z(), 1e-6);

       // the vector shouldn't change under translation
       CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE(
           "TestVecTranslate", (newVec - vecToTransform).mag(), 0, 1e-6);
     }

     void testRotationTransform() {
       // Z now points along y
       GlobalVector newUz(0, 1, 0);
       GlobalPoint ptToTransform(1, 0, 1);
       GlobalVector vecToTransform(1, 0, 1);
       GlobalPoint newPt = transform(origin_, newUz, ptToTransform);
       GlobalVector newVec = transform(origin_, newUz, vecToTransform);

       CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE(
           "TestPtRotate_x", newPt.x(), 1, 1e-6);
       // The component that pointed in z should now be along y
       CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE(
           "TestPtRotate_y", newPt.y(), -1, 1e-6);
       CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE(
           "TestPtRotate_z", newPt.z(), 0, 1e-6);

       CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE(
           "TestVecRotate_x", newVec.x(), 1, 1e-6);
       // The component that pointed in z should now be along y
       CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE(
           "TestVecRotate_y", newVec.y(), -1, 1e-6);
       CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE(
           "TestVecRotate_z", newVec.z(), 0, 1e-6);

       // As y is the new Z axis, the new Y should lie along Z
       GlobalVector newYAxis = transform(origin_, newUz, yAxis_);
       CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE(
           "TestNewYAxis", newYAxis.dot(zAxis_), 1, 1e-6);
     }

     // Test what happens in the case that the line intersects with the cone
     // in two places.  In this case, the one closer to the origin
     // should be taken.
     void testDoubleIntersection() {
       double angle = pi_/5;
       GlobalPoint lineOrigin(3, 0, 0);
       GlobalVector lineDirection(-1, 0, 1); // 45 degrees

       int status = -999;
       GlobalPoint intersection =
         intersectionOfLineAndCone(
             lineOrigin, lineDirection,
             origin_, zAxis_, angle, status);
       CPPUNIT_ASSERT(status == 1);

       // Make sure we haven't chosen the -X solution that has a larger radius
       CPPUNIT_ASSERT(intersection.x() > 0);

       // Using the law of sines to double check the expected radius
       CPPUNIT_ASSERT_DOUBLES_EQUAL(
           intersection.mag()*TMath::Sqrt(2), // intersection radius/sin(45)
           3/TMath::Sin(angle + pi_/4),
           1e-6);
     }

     // Test the functionality that determines if the line origin is inside
     // the cone.
     void testLineOriginInsideConeCheck() {
       double angle = pi_/4; // 45 degrees
       CPPUNIT_ASSERT_MESSAGE(
           "Checking inside returns true",
           pointIsInsideCone(
             GlobalPoint(2, 2, 3), origin_, zAxis_, angle) == true);
       CPPUNIT_ASSERT_MESSAGE(
           "Checking outside returns false",
           pointIsInsideCone(
             GlobalPoint(4, 4, 3), origin_, zAxis_, angle) == false);
     }

     // Test the functionality to move the line origin to the plane
     // perpindicular to the cone direction that contains the cone vertex
     void testPropagateLineOriginToConeVertexPlane() {
       double angle = pi_/5;
       GlobalPoint coneOrigin(2, 2, 2);

       GlobalPoint lineOriginInside(2, 2, 3);
       GlobalVector lineDirection(2, 0, 2); // 45 degrees

       // Making sure we have setup the test correction
       CPPUNIT_ASSERT_MESSAGE(
           "Test Sanity check: line origin is inside cone",
           pointIsInsideCone(
             lineOriginInside, coneOrigin, zAxis_, angle));

       int status = -999;
       GlobalPoint newOrigin = originAtConeVertexPlane(
           lineOriginInside, lineDirection,
           coneOrigin, zAxis_, status);
       CPPUNIT_ASSERT_MESSAGE("Checking success", status == 1);

       CPPUNIT_ASSERT_MESSAGE(
           "Testing to make sure point is out of cone",
           pointIsInsideCone(
             newOrigin, coneOrigin, zAxis_, angle) == false);

       // We expect the Z of the new origin to be the same as the cone apex
       // since the cone points along z
       CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE(
           "Check z coordinate", newOrigin.z(), 2, 1e-6);
     }

     // Make sure no valid point is return when trying to find the PCA
     // of a point behind the apex.
     void testPCAPointBehindPV() {
       double angle = pi_/5;
       GlobalPoint ptToTest(3, 0, -1);
       int status = -999;
       GlobalPoint badPoint = pcaOfConeToPoint(
           ptToTest, origin_, zAxis_, angle, status);
       std::stringstream error;
       error <<  "Shoudl fail but doesnt! testPt: " << ptToTest << std::endl;
       error <<  "coneApex: " << origin_ << std::endl;;
       error <<  "returned point: " << badPoint << std::endl;;
       CPPUNIT_ASSERT_MESSAGE(error.str(), status == 0);
     }

     // See what happens if the track passes through the origin
     void testLineConeCollinear() {
       double angle = pi_/5;
       GlobalPoint lineOriginOnZ(0, 0, 3);
       int status = -999;
       GlobalPoint intersection = intersectionOfLineAndCone(
           lineOriginOnZ, zAxis_, origin_, zAxis_, angle, status);
       CPPUNIT_ASSERT_MESSAGE(dump(intersection), status == 1);
       CPPUNIT_ASSERT_DOUBLES_EQUAL(intersection.mag(), 0, 1e-8);

       // Check if we are behind the origin
       lineOriginOnZ = GlobalPoint(0, 0, -3);
       status = -999;
       intersection = intersectionOfLineAndCone(
           lineOriginOnZ, zAxis_, origin_, zAxis_, angle, status);
       CPPUNIT_ASSERT_MESSAGE(dump(intersection), status == 1);
       CPPUNIT_ASSERT_DOUBLES_EQUAL(intersection.mag(), 0, 1e-8);
     }

     void testLinePropagation() {
       GlobalPoint lineOrigin(3, 0, 3);
       GlobalVector lineDirection(0, 0, 1);
       int status = -999;
       GlobalPoint propagatedPoint = propagateTrackToDistanceWrtConeVertex(
           lineOrigin, lineDirection, origin_, zAxis_, 10, status);
       CPPUNIT_ASSERT(status);
       CPPUNIT_ASSERT_DOUBLES_EQUAL(propagatedPoint.x(), 3, 1e-9);
       CPPUNIT_ASSERT_DOUBLES_EQUAL(propagatedPoint.y(), 0, 1e-9);
       CPPUNIT_ASSERT_DOUBLES_EQUAL(propagatedPoint.z(), 10, 1e-9);
     }

     void testTrackCorrections() {
       // Try and rotate X into Y
       GlobalVector newY = applyPhiAndRadiusCorrections(
           zAxis_,
           xAxis_, // to correct,
           pi_/2, // angle to rotate
           0 // radial correction
           );
       CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE(
           "Testing rotation of X onto Y",
           newY.dot(yAxis_), 1, 1e-8);
       CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE(
           "Testing magnitude doesn't change on rotation",
           newY.mag(), 1, 1e-8);

       GlobalVector longer = applyPhiAndRadiusCorrections(
           zAxis_,
           xAxis_,
           0, 2);
       CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE(
           "Testing magnitude positive shift",
           longer.x(), 3, 1e-8);
       GlobalVector shorter = applyPhiAndRadiusCorrections(
           zAxis_,
           xAxis_,
           0, -2);
       CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE(
           "Testing magnitude negative shift",
           longer.x(), -1, 1e-8);
     }

  private:
     double pi_;
     GlobalPoint origin_;
     GlobalVector zAxis_;
     GlobalVector yAxis_;
     GlobalVector xAxis_;

};

CPPUNIT_TEST_SUITE_REGISTRATION(testSVFitTrack);
