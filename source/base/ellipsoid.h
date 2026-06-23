/***********************************************/
/**
* @file ellipsoid.h
*
* @brief An ellipsoidal model.
*
* @author Torsten Mayer-Guerr
* @date 2004-10-25
*
*/
/***********************************************/

#ifndef __GROOPS_ELLIPSOID__
#define __GROOPS_ELLIPSOID__

#include "base/importStd.h"
#include "base/angle.h"
#include "base/vector3d.h"

/***** CLASS ***********************************/

/** @brief An ellipsoidal model.
* @ingroup vector3dGroup */
class Ellipsoid
{
  Double _a,_b;

public:
  /**
   * @brief Initializes an ellipsoid with given semi-major axis and inverse flattening.
   * @param a Semi-major axis [m]
   * @param f Inverse flattening. If f=0, a sphere is assumed.
   */
  Ellipsoid(Double a=DEFAULT_GRS80_a, Double f=DEFAULT_GRS80_f) : _a(a), _b((f != 0.) ? (a*(1-1/f)) : a) {}

  /**
   * @brief Computes geodetic/ellipsoidal coordinates from Cartesian coordinates.
   * @param[in] point Cartesian coordinates
   * @param[out] L longitude (-PI,PI]
   * @param[out] B latitude [-PI,PI]
   * @param[out] h height [m]
   */
  void operator()(const Vector3d &point, Angle &L, Angle &B, Double &h) const;

  /**
   * @brief Computes Cartesian coordinates from geodetic/ellipsoidal coordinates.
   * @param L longitude (-PI,PI]
   * @param B latitude [-PI/2,PI/2]
   * @param h height [m]
   */
  const Vector3d operator()(Angle L, Angle B, Double h) const;

  /** @brief Returns the semi-major axis. */
  Double a() const {return _a;}

  /** @brief Returns the semi-minor axis. */
  Double b() const {return _b;}

  /** @brief Returns the first numerical excentricity. */
  Double e() const {return sqrt(_a*_a-_b*_b)/_a;}

  /** @brief Returns the flattening. */
  Double f() const {return (_a-_b)/_a;}

  /**
   * @brief Calculates the ellipsoidal velocity from the local NEU velocity.
   * @param B Latitude [-PI/2,PI/2]
   * @param h Height [m]
   * @param vneu Velocity in local NEU system [m/s]
   * @return Latitude rate, longitude rate, and height rate [rad/s, rad/s, m/s]
   */
  const Vector3d ellipsoidalVelocity(Angle B, Double h, const Vector3d &vneu) const;
};

/*************************************************/

#endif
