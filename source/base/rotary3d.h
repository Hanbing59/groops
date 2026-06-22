/***********************************************/
/**
* @file rotary3d.h
*
* @brief Rotations in 3d space.
*
* @author Torsten Mayer-Guerr
* @date 2001-05-31
*
*/
/***********************************************/

#ifndef __GROOPS_ROTARY3D__
#define __GROOPS_ROTARY3D__

#include "base/importStd.h"
#include "base/matrix.h"

/** @addtogroup vector3dGroup */
/// @{

class Angle;
class Vector3d;
class Tensor3d;
class Transform3d;
class Ellipsoid;
class SphericalHarmonics;

/***** CLASS ***********************************/

/** @brief Rotation of a @a Vector3d or @a Tensor3d object in 3d space. */
class Rotary3d
{
  /// Rotation matrix (3x3) in the form of a 2d std::array.
  std::array<std::array<Double,3>,3> field;

public:
  /** @brief Default constructor (Unitary matrix). */
  Rotary3d() : field{{{1.,0.,0.}, {0.,1.,0.}, {0.,0.,1.}}} {}

  /** @brief Constructor with a 2d std::array. */
  Rotary3d(const std::array<std::array<Double,3>,3> &x) : field(x) {}

  /** @brief Constructor with a quaternions vector (q0, qx, qy, qz) or a 3x3 rotary matrix. */
  explicit Rotary3d(const_MatrixSliceRef q);

  /**
   * @brief Constructor with the x- and y-axis of a local system in the global system.
   * It fills @a field with the rotation matrix from the local system to the global system.
   * The z-axis will be defined to be orthogonal to the x- and y-axis according to
   * the right-hand rule. And the y-axis will be re-defined to be orthogonal to
   * the new z-axis and the given x-axis according to the right-hand rule.
   * @param x x-axis of the local system in the global system.
   * @param y y-axis of the local system in the global system.
   */
  Rotary3d(Vector3d x, Vector3d y);

  /** @brief Returns the quaternion representation (q0, qx, qy, qz) of the rotation. */
  Vector quaternion() const;

  /** @brief Returns @a field as a 3x3 matrix. */
  Matrix matrix() const;

  /** @brief Applies the rotation to a @a Vector3d object.
  * \f[ y = D \cdot v \f]. */
  Vector3d rotate(const Vector3d &v) const;
  /** @copydoc rotate */
  Vector3d transform(const Vector3d &v) const;

  /** @brief Applies the inverse rotation to a @a Vector3d object.
  * \f[ y = D^T \cdot v \f]. */
  Vector3d inverseRotate(const Vector3d &v) const;
  /** @copydoc inverseRotate */
  Vector3d inverseTransform(const Vector3d &v) const;

  /** @brief Applies the rotation to a @a Tensor3d object.
  * Both sides of the dyadic product are rotated.
  * \f[ y = D \cdot t \cdot D^T \f]. */
  Tensor3d rotate(const Tensor3d &t) const;
  /** @copydoc rotate */
  Tensor3d transform(const Tensor3d &t) const;

  /** @brief Applies the inverse rotation to a @a Tensor3d object.
  * Both sides of the dyadic product are rotated.
  * \f[ y = D^T \cdot t \cdot D \f]. */
  Tensor3d inverseRotate(const Tensor3d &t) const;
  /** @copydoc inverseRotate */
  Tensor3d inverseTransform(const Tensor3d &t) const;

  /** @brief Applies the rotation to a @a SphericalHarmonics object. */
  SphericalHarmonics rotate(const SphericalHarmonics &harm) const;

  /** @brief Applies the inverse rotation to a @a SphericalHarmonics object. */
  SphericalHarmonics inverseRotate(const SphericalHarmonics &harm) const;

  /** @brief Returns proper/classic Euler angles representation of the rotation.
  * @code rotary = rotaryZ(gamma)*rotaryX(beta)*rotaryZ(alpha); @endcode */
  void euler(Angle &alpha, Angle &beta, Angle &gamma) const;

  /** @brief Returns Cardan angles representation of the rotation.
  * @code rotary = rotaryZ(yaw)*rotaryY(pitch)*rotaryX(roll); @endcode */
  void cardan(Angle &roll, Angle &pitch, Angle &yaw) const;

  /** @brief Multiplies the current rotation with a given rotation. */
  Rotary3d    &operator*=(Rotary3d const &b);
  /** @brief Multiplies the current rotation with a given rotation. */
  Rotary3d     operator* (Rotary3d const &b) const;

  /** @brief Multiplies the current rotation with a given transformation and
   * returns the result as a @a Transform3d object. */
  Transform3d  operator* (Transform3d const &b) const;

  friend class Transform3d;
  friend Rotary3d rotaryX(Angle angle);
  friend Rotary3d rotaryY(Angle angle);
  friend Rotary3d rotaryZ(Angle angle);
  friend Rotary3d inverse(const Rotary3d &b);
  friend Rotary3d localNorthEastDown(const Vector3d &point);
};

/***********************************************/

/** @brief Rotation about x-axis.
* @param angle in [rad].
* \f[ D=\left(\begin{array}{ccc}
      1 & 0          & 0          \\
      0 & \cos\alpha & \sin\alpha \\
      0 &-\sin\alpha & \cos\alpha
      \end{array}\right) \f] */
Rotary3d rotaryX(Angle angle);

/** @brief Rotation about y-axis.
* @param angle in [rad].
* \f[ D=\left(\begin{array}{ccc}
      \cos\alpha & 0 &-\sin\alpha \\
      0          & 1 & 0          \\
      \sin\alpha & 0 & \cos\alpha
      \end{array}\right) \f] */
Rotary3d rotaryY(Angle angle);

/** @brief Rotation about z-axis.
* @param angle in [rad].
* \f[ D=\left(\begin{array}{ccc}
      \cos\alpha & \sin\alpha & 0 \\
     -\sin\alpha & \cos\alpha & 0 \\
      0          & 0          & 1
      \end{array}\right) \f] */
Rotary3d rotaryZ(Angle angle);

/** @brief Returns the inverse rotation of a given rotation by transposing the rotation matrix. */
Rotary3d inverse(const Rotary3d &b);

/**
 * @brief Returns a rotation matrix from the local right-handed NED system
 * at a given point to the global system.
 * @note The @a down direction of the local NED system is defined to be along the geocentric radius vector.
 */
//[[deprecated]]
Rotary3d localNorthEastDown(const Vector3d &point);

/**
 * @brief Returns a rotation matrix from the local right-handed NED system
 * at a given point to the global system.
 * @note The @a down direction of the local NED system is defined to be the ellipsoidal normal vector.
 */
//[[deprecated]]
Rotary3d localNorthEastDown(const Vector3d &point, const Ellipsoid &ellipsoid);

/// @}

/***********************************************/
/***** INLINES *********************************/
/***********************************************/

inline Matrix Rotary3d::matrix() const
{
  Matrix R(3, 3, Matrix::NOFILL);
  for(UInt i=0; i<3; i++)
    for(UInt k=0; k<3; k++)
      R(i,k) = field[i][k];
  return R;
}

/*************************************************/

#endif /* __GROOPS_VECTOR3D__ */
