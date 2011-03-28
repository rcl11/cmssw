#include "TauAnalysis/CandidateTools/interface/SVfitTrackTools.h"
#include "TMath.h"
#include <TRotation.h>
#include <iostream>
#include "DataFormats/Math/interface/angle.h"

using namespace TMath;

namespace {
  inline double square(double x) { return x*x; }
}

namespace SVfit { namespace track {

GlobalPoint propagateLine(const GlobalPoint& origin,
    const GlobalVector& tangent, double pathLength) {
  return origin + tangent*pathLength;
}

// Get the correct path length solution of the quadratic equation.
// We always prefer the closest solution that is not negative - i.e. below the
// the cone vertex.
double solveQuadraticForPathLength(double a, double b, double c, int& status) {
  // Check if solutions are real and finite.
  if ( b*b - 4*a*c < 0 || a == 0 ) {
    status = 0;
    return 0;
  }
  // We always want the smaller one, if it is above zero
  double firstSolution =  (-b + Sqrt(b*b - 4*a*c))/(2*a);
  double secondSolution =  (-b - Sqrt(b*b - 4*a*c))/(2*a);
  bool oppositeSignSolutions = ((firstSolution * secondSolution) < 0);

  //std::cout << "first soln: " << firstSolution <<
    //"second soln: " << secondSolution << std::endl;
  // If one is negative, and one is positive, return the positive one.
  if (oppositeSignSolutions) {
    status = 1;
    double output = std::max(firstSolution, secondSolution);
    //    std::cout << "opp sign-returning: " << output << std::endl;
    return output;
  }
  // Check if the solutions are both negative, that's an error.
  if (firstSolution < 0) {
    status = 0;
    //    std::cout << "both negative!" << std::endl;
    return 0;
  }
  // The solutions are both positive.  Take the one that is closer to the origin
  status = 1;
  double output = std::min(firstSolution, secondSolution);
  //  std::cout << "returning: " << output << std::endl;
  return output;
}


GlobalPoint intersectionOfLineAndCone(
    const GlobalPoint &lineOffsetIn, const GlobalVector &lineDirectionIn,
    const GlobalPoint &coneVertexIn, const GlobalVector &coneDirectionIn,
    double alpha /* cone angle */, int &status) {

  // Check if line origin is inside the cone.  If so, move it out.
  GlobalPoint lineOffset = lineOffsetIn;

  GlobalPoint originAtVertexPlane = originAtConeVertexPlane(
      lineOffsetIn, lineDirectionIn, coneVertexIn, coneDirectionIn, status);

  // Check if the line passes through the cone vertex - if so take it as it
  // is a special case for the equations used to find the equations.
  if (status &&
      vectorSubtract(originAtVertexPlane, coneVertexIn).mag2() < 1e-12) {
    return coneVertexIn;
  }

  // We only do this in the case that the direction of the line is in the same
  // direction as the cone.  We use a backward facing line, from the axis
  // in the PCA computations.
  if (lineDirectionIn.dot(coneDirectionIn) > 0 &&
      pointIsInsideCone(lineOffsetIn, coneVertexIn, coneDirectionIn, alpha)) {
    lineOffset = originAtVertexPlane;
    // This can only happen in the case that the track is perpendicular to
    // the cone direction.
    if (!status) {
      return GlobalPoint();
    }
  }

  // Shift coordinates such that the cone vertex is at the origin.
  lineOffset = transform(coneVertexIn, coneDirectionIn, lineOffset);

  // Now rotate so the cone is aligned along Z
  GlobalVector lineDirection = transform(coneVertexIn, coneDirectionIn,
                                         lineDirectionIn);

//  std::cout << "after transform: " << lineDirection << std::endl;

  // Okay, now we are ready to go.  Define our variables.  These match up to the
  // variable names in the corresponding Mathematic derivation.
  double trackDirX = lineDirection.x();
  double trackDirY = lineDirection.y();
  double trackDirZ = lineDirection.z();
  double trackOffsetX = lineOffset.x();
  double trackOffsetY = lineOffset.y();
  double trackOffsetZ = lineOffset.z();

  // Get quadratic coefficients
  double a = square(trackDirX) + square(trackDirY) - square(trackDirZ)*square(Tan(alpha));

  double b = 2*trackDirX*trackOffsetX + 2*trackDirY*trackOffsetY - 2*trackDirZ*trackOffsetZ*square(Tan(alpha));

  double c = square(trackOffsetX) + square(trackOffsetY) - square(trackOffsetZ)*square(Tan(alpha));

  // Find the correct solution to our quadratic equation.
  double solution = solveQuadraticForPathLength(a, b, c, status);
  // If it is a good solution, find the corresponding point in space and then
  // transform back to the original coordinate system.
  if (status) {
    return untransform(coneVertexIn, coneDirectionIn,
                       propagateLine(lineOffset, lineDirection, solution));
  } else
    return GlobalPoint();
}

GlobalPoint pcaOfLineToCone(
    const GlobalPoint &lineOffsetIn, const GlobalVector &lineDirectionIn,
    const GlobalPoint &coneVertexIn, const GlobalVector &coneDirectionIn,
    double alpha /* cone angle */, int &status) {

  // Shift coordinates such that the cone vertex is at the origin.
  GlobalPoint lineOffset = transform(coneVertexIn, coneDirectionIn,
                                     lineOffsetIn);

  // Now rotate so the cone is aligned along Z
  GlobalVector lineDirection = transform(coneVertexIn, coneDirectionIn,
                                         lineDirectionIn);

  double trackDirX = lineDirection.x();
  double trackDirY = lineDirection.y();
  double trackDirZ = lineDirection.z();
  double trackOffsetX = lineOffset.x();
  double trackOffsetY = lineOffset.y();

  // If the track is perpindicular to the cone direction (which should never
  // happen in practice), it's not possible to find a solution with the current
  // parameterization.
  if (trackDirZ == 0.) {
    status = 0;
    return GlobalPoint();
  }

  // Get quadratic coefficients for the equation to solve.
  double a = -((square(trackDirX) + square(trackDirY))*
               (square(trackDirX) + square(trackDirY) - square(trackDirZ) +
                (square(trackDirX) + square(trackDirY) +
                 square(trackDirZ))*Cos(2*alpha))) /
      (2.*square(trackDirZ));

  double b = -(((trackDirX*trackOffsetX + trackDirY*trackOffsetY)*
                (square(trackDirX) + square(trackDirY) - square(trackDirZ) +
                 (square(trackDirX) + square(trackDirY) + square(trackDirZ))
                 *Cos(2*alpha)))/square(trackDirZ));

  double c = (-square(trackDirX*trackOffsetX + trackDirY*trackOffsetY) + (square(trackDirX)*square(trackOffsetX) + 2*trackDirX*trackDirY*trackOffsetX*trackOffsetY +
        square(trackDirY)*square(trackOffsetY) + square(trackDirZ)*(square(trackOffsetX) + square(trackOffsetY)))*square(Sin(alpha)))/square(trackDirZ);

  // Find the correct solution to our quadratic equation.  This will give us the
  // t parameter that gives the point on the line closest to the
  // cone.
  double solution = solveQuadraticForPathLength(a, b, c, status);

  // Check if a solution exists
  if (!status) {
    return GlobalPoint();
  }

  GlobalPoint linePCA = propagateLine(lineOffset, lineDirection, solution);

  // Now find the point on the cone closest to this point.  Note that we've
  // already transformed our coordinates, we don't need to do it again.
  //GlobalPoint conePCA = pcaOfConeToPoint(linePCA, GlobalPoint(0,0,0),
  //                                       GlobalVector(0, 0, 1), alpha, status);
  //return untransform(coneVertexIn, coneDirectionIn, conePCA);
  return untransform(coneVertexIn, coneDirectionIn, linePCA);
}

GlobalPoint pcaOfConeToPoint(
    const GlobalPoint& pointIn,
    const GlobalPoint &coneVertex, const GlobalVector &coneDirection,
    double alpha, int &status) {

  // If the point is "behind" the cone apex, there is no valid PCA (other
  // than the cone tip).
  if (vectorSubtract(pointIn, coneVertex).dot(coneDirection) < 0)  {
    status = 0;
    return coneVertex;
  }

  // Shift coordinates such that the cone vertex is at the origin.
  GlobalPoint point = transform(coneVertex, coneDirection, pointIn);

  // Now we need to find the point on the cone.  We do this by first finding the
  // z0 position of the apex of the inverted cone (with angle 90-alpha) that
  // contains our point.
  double z0 = point.z() +
      Sqrt(point.x()*point.x() + point.y()*point.y())*Tan(alpha);

  GlobalPoint z0point(0, 0, z0);

  // Construct a vector from the z0 point to the point.  The
  // intersection of this line to the cone is the PCA for the CONE.  This line
  // is by construction perpindicular to the cone at the intersection.
  GlobalVector displacementDirection = point - z0point;

  // Note that when we call this function, we have already transformed into a
  // coordinate system where the cone apex is at the origin and the cone
  // direction is along the z axis.
  GlobalPoint conePCA = intersectionOfLineAndCone(z0point,
                                                  displacementDirection,
                                                  GlobalPoint(0,0,0),
                                                  GlobalVector(0,0,1),
                                                  alpha, status);
  if (!status) {
    return GlobalPoint();
  }
  // Get back to our original coordinate system.
  return untransform(coneVertex, coneDirection,conePCA);
}

/// Check whether a point is inside or outside a cone.
bool pointIsInsideCone(const GlobalPoint& point,
    const GlobalPoint &coneVertex, const GlobalVector &coneDirection,
    double coneAngle) {
  GlobalVector displacementFromVertex = point - coneVertex;
  return angle(displacementFromVertex, coneDirection) < coneAngle;
}

GlobalPoint originAtConeVertexPlane(
    const GlobalPoint &lineOffset, const GlobalVector &lineDirection,
    const GlobalPoint &coneVertex, const GlobalVector &coneDirection,
    int& status) {
  // In the case that the line direction and the cone direction are orthogonal,
  // its impossible to find a valid point.
  double projLineOnCone = lineDirection.dot(coneDirection);
  if (projLineOnCone == 0) {
    status = 0;
    return lineOffset;
  }
  status = 1;
  double projLineOffsetOnCone =
    vectorSubtract(lineOffset,coneVertex).dot(coneDirection);
  double pathLength = -projLineOffsetOnCone/projLineOnCone;
  return propagateLine(lineOffset, lineDirection, pathLength);
}

GlobalPoint propagateTrackToDistanceWrtConeVertex(
    const GlobalPoint &lineOffset, const GlobalVector &lineDirection,
    const GlobalPoint &coneVertex, const GlobalVector &coneDirection,
    double distance, int &status) {
  GlobalPoint correctedLineOrigin = originAtConeVertexPlane(
      lineOffset, lineDirection, coneVertex, coneDirection, status);
  return propagateLine(correctedLineOrigin, lineDirection.unit(), distance);
}

GlobalVector applyPhiAndRadiusCorrections(
    const GlobalVector& axis,
    const GlobalVector& toCorrect,
    double phiCorrection, double radiusCorrection) {
  TRotation rotation;
  rotation.Rotate(phiCorrection, convert<TVector3>(axis));
  TVector3 rotated = rotation*convert<TVector3>(toCorrect);
  double length = rotated.Mag();
  double scaleFactor = (length > 0) ?
    (length + radiusCorrection)/length : 0;
  rotated *= scaleFactor;
  return convert<GlobalVector>(rotated);
}

GlobalPoint transform(const GlobalPoint& newOrigin, const GlobalVector &newUz,
                      const GlobalPoint& toTransform) {
  // Shift origin
  GlobalPoint shifted = toTransform;
  shifted -= convert<GlobalVector>(newOrigin);
  // Rotate
  TRotation rotation;
  rotation.SetZAxis(convert<TVector3>(newUz).Unit());
  TVector3 rotated = rotation.Inverse()*convert<TVector3>(shifted);
  return convert<GlobalPoint>(rotated);
}

GlobalVector transform(const GlobalPoint& newOrigin, const GlobalVector &newUz,
                      const GlobalVector& toTransform) {
  // We can ignore the origin shift for vectors.
  TRotation rotation;
  rotation.SetZAxis(convert<TVector3>(newUz).Unit());
  TVector3 rotated = rotation.Inverse()*convert<TVector3>(toTransform);
  return convert<GlobalVector>(rotated);
}

GlobalPoint untransform(const GlobalPoint& newOrigin, const GlobalVector &newUz,
                      const GlobalPoint& toTransform) {
  // Rotate
  TRotation rotation;
  rotation.SetZAxis(convert<TVector3>(newUz).Unit());
  // Invert rotation.
  TVector3 unrotated = rotation*convert<TVector3>(toTransform);
  // Shift origin back
  GlobalPoint output = convert<GlobalPoint>(unrotated);
  output += convert<GlobalVector>(newOrigin);
  return output;
}

GlobalVector untransform(const GlobalPoint& newOrigin, const GlobalVector &newUz,
                         const GlobalVector& toTransform) {
  TRotation rotation;
  rotation.SetZAxis(convert<TVector3>(newUz).Unit());
  // Invert rotation.
  return convert<GlobalVector>(
      rotation*convert<TVector3>(toTransform));
}

}}  // end namespace SVfit::track
