#ifndef TauAanlysis_CandidateTools_SVfitTrackTools_h
#define TauAanlysis_CandidateTools_SVfitTrackTools_h

/*
 * Geometrical tools used in the SV fit method.  Specifically, functions to find
 * the intersection or point of closest approach between lines and cones.
 *
 * Authors: Evan K. Friis, Christian Veelken (UC Davis)
 *
 */

#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/GeometryVector/interface/GlobalVector.h"

namespace SVfit { namespace track {

// Stupid type conversion
template<typename Out, typename In>
inline Out convert(const In& in) { return Out(in.x(), in.y(), in.z()); }

template<typename T1, typename T2>
inline GlobalVector vectorSubtract(const T1& t1, const T2& t2) {
  return GlobalVector(t1.x() - t2.x(), t1.y() - t2.y(), t1.z() - t2.z());
}

/// Find the point of intersection between a line and a cone.  status = 1 if a
/// solution found, 0 if not.
GlobalPoint intersectionOfLineAndCone(
    const GlobalPoint &lineOffset, const GlobalVector &lineDirection,
    const GlobalPoint &coneVertex, const GlobalVector &coneDirection,
    double coneAngle, int &status);


/// Find the point of closest approach of a line to a cone.  status = 1 if
/// solution found, 0 if not.
GlobalPoint pcaOfLineToCone(
    const GlobalPoint &lineOffset, const GlobalVector &lineDirection,
    const GlobalPoint &coneVertex, const GlobalVector &coneDirection,
    double coneAngle, int &status);


/// Find the point *on the cone* closest to the input point.
GlobalPoint pcaOfConeToPoint(
    const GlobalPoint& point,
    const GlobalPoint &coneVertex, const GlobalVector &coneDirection,
    double coneAngle, int &status);

/// Check whether a point is inside or outside a cone.
bool pointIsInsideCone(const GlobalPoint& point,
    const GlobalPoint &coneVertex, const GlobalVector &coneDirection,
    double coneAngle);

///  Find a new line offset such that the lineOffset point lies in the plane
/// that is perpindicualr to the cone direction and contains the cone apex.
/// Sets status to 1 if successfull
GlobalPoint originAtConeVertexPlane(
    const GlobalPoint &lineOffset, const GlobalVector &lineDirection,
    const GlobalPoint &coneVertex, const GlobalVector &coneDirection,
    int &status);

/// Propagate the track to a point such that the path length between the
/// distance between the propagated point and the point on the line that lies
/// on the cone vertex plane (see above function) is d
GlobalPoint propagateTrackToDistanceWrtConeVertex(
    const GlobalPoint &lineOffset, const GlobalVector &lineDirection,
    const GlobalPoint &coneVertex, const GlobalVector &coneDirection,
    double distance, int &status);

GlobalVector applyPhiAndRadiusCorrections(
    const GlobalVector& axis,
    const GlobalVector& toCorrect,
    double phiCorrection, double radiusCorrection);

GlobalPoint transform(const GlobalPoint& newOrigin, const GlobalVector &newUz,
                      const GlobalPoint& toTransform);
GlobalPoint untransform(const GlobalPoint& newOrigin, const GlobalVector &newUz,
                      const GlobalPoint& toTransform);

GlobalVector transform(const GlobalPoint& newOrigin, const GlobalVector &newUz,
                      const GlobalVector& toTransform);
GlobalVector untransform(const GlobalPoint& newOrigin, const GlobalVector &newUz,
                      const GlobalVector& toTransform);

}}  // end namespace SVfit::track

#endif
