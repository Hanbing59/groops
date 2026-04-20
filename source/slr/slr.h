/***********************************************/
/**
* @file slr.h
*
* @brief Satellite Laser Ranging classes.
*
* @author Torsten Mayer-Guerr
* @date 2022-04-28
*
*/
/***********************************************/

#ifndef __GROOPS_SLR__
#define __GROOPS_SLR__

#include "parallel/parallel.h"
#include "base/parameterName.h"
#include "classes/noiseGenerator/noiseGenerator.h"
#include "slr/slrObservation.h"
#include "slr/slrDesignMatrix.h"
#include "slr/slrSatellite.h"
#include "slr/slrStation.h"
#include "slr/slrNormalEquationInfo.h"

/** @addtogroup slrGroup */
/// @{

/***** TYPES ***********************************/

class Slr;
typedef std::shared_ptr<Slr> SlrPtr;

class MatrixDistributed;
class SlrSatelliteGenerator;
class SlrStationGenerator;
class SlrParametrization;
class EarthRotation;
class PlatformSelector;
typedef std::shared_ptr<SlrSatelliteGenerator> SlrSatelliteGeneratorPtr;
typedef std::shared_ptr<SlrStationGenerator>   SlrStationGeneratorPtr;
typedef std::shared_ptr<SlrParametrization>    SlrParametrizationPtr;
typedef std::shared_ptr<EarthRotation>         EarthRotationPtr;
typedef std::shared_ptr<PlatformSelector>      PlatformSelectorPtr;

/***** CLASS ***********************************/

/** @brief Satellite Laser Ranging Class */
class Slr
{
public:
  /// The list of satellites
  std::vector<SlrSatellitePtr>                     satellites;
  /// The list of stations
  std::vector<SlrStationPtr>                       stations;
  SlrParametrizationPtr                            parametrization;
  std::function<void(SlrObservationEquation &eqn)> funcReduceModels;
  std::function<Rotary3d(const Time &time)>        funcRotationCrf2Trf;

  Polynomial        polynomialEop;
  std::vector<Time> times;
  /// Matrix eop columns: xp, yp, sp, deltaUT, LOD, X, Y, S.
  Matrix            eop;

  /** @brief Do some initializations.
   * @param times
   * @param satelliteGenerator
   * @param stationGenerator
   * @param earthRotation
   * @param parametrization
   */
  void init(const std::vector<Time> &times, SlrSatelliteGeneratorPtr satelliteGenerator, SlrStationGeneratorPtr stationGenerator,
            EarthRotationPtr earthRotation, SlrParametrizationPtr parametrization);

  /** @brief Inertial system (CRF) -> earth fixed system (TRF). */
  Rotary3d rotationCrf2Trf(const Time &time) const;

  /** @brief Initialize the normal equation. */
  void   initParameter            (SlrNormalEquationInfo &normalEquationInfo);
  /** @brief Returns a-priori values of parameters */
  Vector aprioriParameter         (const SlrNormalEquationInfo &normalEquationInfo) const;

  /** @brief Form the basic observation equations of one pass of SLR observations between a station and a satellite. */
  Bool   basicObservationEquations(const SlrNormalEquationInfo &normalEquationInfo, UInt idStat, UInt idSat, UInt idPass, SlrObservationEquation &eqn) const;
  /** @brief Generate the design matrix */
  void   designMatrix             (const SlrNormalEquationInfo &normalEquationInfo, const SlrObservationEquation &eqn, SlrDesignMatrix &A) const;
  /** @brief Apply the constraints */
  void   constraints              (const SlrNormalEquationInfo &normalEquationInfo, MatrixDistributed &normals, std::vector<Matrix> &n, Double &lPl, UInt &obsCount) const;
  /** @brief Update parameters */
  Double updateParameter          (const SlrNormalEquationInfo &normalEquationInfo, const_MatrixSliceRef x, const_MatrixSliceRef Wz);
  /** @brief Update the covariance */
  void   updateCovariance         (const SlrNormalEquationInfo &normalEquationInfo, const MatrixDistributed &covariance);
  /** @brief Write results out */
  void   writeResults             (const SlrNormalEquationInfo &normalEquationInfo, const std::string &suffix="");

  /** @brief Select only usable stations */
  std::vector<Byte> selectStations(PlatformSelectorPtr selector);
  /** @brief Select only usable satellites */
  std::vector<Byte> selectSatellites(PlatformSelectorPtr selector);

  /** @brief  Internal class, for parameters change */
  class InfoParameterChange
  {
  public:
    std::string unit;
    UInt        count;
    Double      rms;
    Double      maxChange;
    std::string info;

    /** @brief Constructor. */
    InfoParameterChange(const std::string &unit) : unit(unit), count(0), rms(0), maxChange(0) {}
    /** @brief Update the RMS and maximum change */
    Bool update(Double change);
    /** @brief Print the RMS and maximum change */
    void print(Double convertToMeter, Double &maxChangeTotal);
  };
};

/// @}

/***********************************************/

#endif /* __GROOPS___ */
