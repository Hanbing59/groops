/***********************************************/
/**
* @file griddedData.h
*
* @brief Gridded values.
*
* @author Torsten Mayer-Guerr
* @date 2005-01-14
*
*/
/***********************************************/

#ifndef __GROOPS_GRIDDEDDATA__
#define __GROOPS_GRIDDEDDATA__

#include "base/importStd.h"
#include "base/vector3d.h"
#include "base/ellipsoid.h"

/** @addtogroup base */
/// @{

class GriddedDataRectangular;

/***** CLASS ***********************************/

/** @brief Point list with (multiple) data. */
class GriddedData
{
public:
  /// The reference ellipsoid for conversion of geodetic/ellipsoidal coordinates
  Ellipsoid                        ellipsoid;
  /// Cartesian coordinates of grid points.
  std::vector<Vector3d>            points;
  /// Area elements (projected to unit sphere) at each point.
  std::vector<Double>              areas;
  /// Data values at each point.
  std::vector<std::vector<Double>> values;

  /** @brief Default constructor. */
  GriddedData() = default;

  /** @brief Constructor with points, area elements, and multiple values for each point. */
  GriddedData(const Ellipsoid &ellip, const std::vector<Vector3d> &_points, const std::vector<Double> &_areas, const std::vector<std::vector<Double>> &_values) : ellipsoid(ellip), points(_points), areas(_areas), values(_values) {}

  /** @brief Constructor with a @a GriddedDataRectangular object. */
  GriddedData(const GriddedDataRectangular &grid) {init(grid);}

  /** @brief Initialization from a @a GriddedDataRectangular object. */
  void init(const GriddedDataRectangular &grid);

  /** @brief Sort points geographically (North/West->South/East). */
  void sort();

  /** @brief Define points a rectangular grid?
  * if function returns FALSE, @a lambda, @a phi, @a radius contain garbage.
  * if function returns TRUE, points are in same order as:
  @code
  UInt idx = 0;
  for(UInt i=0; i<phi.size(); i++)
    for(UInt k=0; k<lambda.size(); k++)
      points.at(idx++) == polar(lambda.at(k), phi.at(i), r.at(i));
  @endcode */
  Bool isRectangle(std::vector<Angle> &lambda, std::vector<Angle> &phi, std::vector<Double> &radius) const;

  /** @brief Automatically area computation of rectangular grids (overwrite areas). */
  Bool computeArea();

  /** @brief Checks if the number of points, areas, and values of this grid data are consistent.*/
  Bool isValid() const;
};

/***** CLASS ***********************************/

/** @brief Rectangular grid with (multiple) data. */
class GriddedDataRectangular
{
public:
  /// The reference ellipsoid for conversion of geodetic/ellipsoidal coordinates
  Ellipsoid           ellipsoid;
  /// Longitudes (columns) of the grid points.
  std::vector<Angle>  longitudes;
  /// Latitudes (rows) of the grid points.
  std::vector<Angle>  latitudes;
  /// Elliposoidal heights (rows) of the grid points.
  std::vector<Double> heights;
  /// Data values at each point.
  std::vector<Matrix> values;

  /** @brief Initialization from GriddedData. */
  Bool init(const GriddedData &grid);

  /** @brief Gets the geocentric polar coordinates of the grid points. */
  void geocentric(std::vector<Angle> &lambda, std::vector<Angle> &phi, std::vector<Double> &radius) const;

  /** @brief borders of grid cell (i,k): (lat(i) - lat(i+1)) x (lon(k) - lon(k+1)). */
  void cellBorders(std::vector<Double> &longitudes, std::vector<Double> &latitudes) const;

  /** @brief Area elements projected on the unit sphere.
  * area(i,k) = dPhi(i)*dLambda(k).
  * @return total area (4pi for global grids). */
  Double areaElements(std::vector<Double> &dLambda, std::vector<Double> &dPhi) const;

  /** @brief Checks if the number of latitudes, heights, and values of this grid data are consistent. */
  Bool isValid() const;
};

/// @}

/***********************************************/

#endif /* __GROOPS__ */
