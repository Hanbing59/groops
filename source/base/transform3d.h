/***********************************************/
/**
* @file transform3d.h
*
* @brief Coordinate transformations, i.e. rotations and reflections in 3d space.
*
* @author Torsten Mayer-Guerr
* @date 2019-03-03
*
*/
/***********************************************/

#ifndef __GROOPS_TRANSFORM3D__
#define __GROOPS_TRANSFORM3D__

#include "base/importStd.h"
#include "base/matrix.h"

/** @addtogroup vector3dGroup */
/// @{

class Angle;
class Vector3d;
class Tensor3d;
class Rotary3d;
class Ellipsoid;

/***** CLASS ***********************************/

/** @brief Coordinate transformations, i.e. rotations and reflections in 3d space. */
class Transform3d
{
  /// Transformation matrix (3x3) in the form of a 2d std::array.
  std::array<std::array<Double,3>,3> field;

public:
  /** @brief Default constructor (Unitary matrix). */
  Transform3d() : field{{{1.,0.,0.}, {0.,1.,0.}, {0.,0.,1.}}} {}

  /** @brief Constructor with a 2d std::array. */
  Transform3d(const std::array<std::array<Double,3>,3> &x) : field(x) {}

  /** @brief Constructor with a @a Rotary3d object. */
  Transform3d(const Rotary3d &rot);

  /** @brief Constructor with a 3x3 transformation matrix. */
  explicit Transform3d(const_MatrixSliceRef A);

  /**
   * @brief Constructor with the x- and y-axis of a local system in the global system.
   * It fills @a field with the rotation matrix from the local system to the global system.
   * The z-axis will be defined to be orthogonal to the x- and y-axis according to
   * the right-hand rule. And the y-axis will be re-defined to be orthogonal to
   * the new z-axis and the given x-axis according to the right-hand rule.
   * @param x x-axis of the local system in the global system.
   * @param y y-axis of the local system in the global system.
   */
  Transform3d(Vector3d x, Vector3d y);

  /** @brief Returns @a field as a 3x3 matrix. */
  Matrix matrix() const;

  /** @brief Applies the transformation to a given vector.
  * \f[ y = T \cdot v \f]. */
  Vector3d transform(const Vector3d &v) const;

  /** @brief Applies the inverse transformation to a given vector.
  * \f[ y = T^T \cdot v \f]. */
  Vector3d inverseTransform(const Vector3d &v) const;

  /** @brief Applies the transformation to a given tensor.
  * Both sides of the dyadic product are rotated.
  * \f[ y = T \cdot t \cdot T^T \f]. */
  Tensor3d transform(const Tensor3d &t) const;

  /** @brief Applies the inverse transformation to a given tensor.
  * Both sides of the dyadic product are rotated.
  * \f[ y = T^T \cdot t \cdot T \f]. */
  Tensor3d inverseTransform(const Tensor3d &t) const;

  Transform3d &operator*=(const Transform3d &b);
  Transform3d &operator*=(const Rotary3d &b);
  /** @brief Multiplication of two transformations both represented as @a Transform3d objects. */
  Transform3d  operator* (const Transform3d &b) const;
  /** @brief Multiplication of two transformations represented as @a Transform3d and @a Rotary3d objects, respectively. */
  Transform3d  operator* (const Rotary3d &b) const;

  friend class Rotary3d;
  friend Transform3d inverse(const Transform3d &b);
  friend Transform3d localNorthEastUp(const Vector3d &point);
  friend Transform3d localNorthEastUp(const Vector3d &point, const Ellipsoid &ellipsoid);
};

/***********************************************/

/** @brief Returns the transformation for flipping the x-axis. */
inline Transform3d flipX() {return Transform3d(std::array<std::array<Double,3>,3>{{{-1.,0.,0.}, {0.,1.,0.}, {0.,0.,1.}}});}

/** @brief Returns the transformation for flipping the y-axis. */
inline Transform3d flipY() {return Transform3d(std::array<std::array<Double,3>,3>{{{1.,0.,0.}, {0.,-1.,0.}, {0.,0.,1.}}});}

/** @brief Returns the transformation for flipping the z-axis. */
inline Transform3d flipZ() {return Transform3d(std::array<std::array<Double,3>,3>{{{1.,0.,0.}, {0.,1.,0.}, {0.,0.,-1.}}});}

/**
 * @brief Returns the inverse transformation, i.e. the transposed transformation matrix.
 * @ingroup vector3dGroup
 */
Transform3d inverse(const Transform3d &b);

/**
 * @brief Returns a rotational transformation from the local left-handed NEU system
 * at a given point to the global system.
 * @note The @a up direction of the local NEU system is defined to be along the geocentric radius vector.
 * @param point Cartesian coordinates of the point.
 * @ingroup vector3dGroup
 */
Transform3d localNorthEastUp(const Vector3d &point);

/**
 * @brief Returns a rotational transformation from the local left-handed NEU system
 * at a given point to the global system.
 * @note The @a up direction of the local NEU system is defined to be the ellipsoidal normal vector
 * @param point Cartesian coordinates of the point.
 * @param ellipsoid Reference ellipsoid.
 * @ingroup vector3dGroup
* Converts @p point to ellipsoidal coordinates beforehand. */
Transform3d localNorthEastUp(const Vector3d &point, const Ellipsoid &ellipsoid);

/// @}

/***********************************************/
/***** INLINES *********************************/
/***********************************************/

inline Matrix Transform3d::matrix() const
{
  Matrix R(3, 3, Matrix::NOFILL);
  for(UInt i=0; i<3; i++)
    for(UInt k=0; k<3; k++)
      R(i,k) = field[i][k];
  return R;
}

/*************************************************/

#endif /* __GROOPS_TRANSFORM3D__ */
