/***********************************************/
/**
* @file slrStation.h
*
* @brief SLR station.
*
* @author Torsten Mayer-Guerr
* @date 2022-04-28
*
*/
/***********************************************/

#ifndef __GROOPS_SLRSTATION__
#define __GROOPS_SLRSTATION__

#include "parser/expressionParser.h"
#include "files/fileInstrument.h"
#include "slr/slrObservation.h"
#include "slr/slrPlatform.h"
#include "slr/slrSatellite.h"

/** @addtogroup slrGroup */
/// @{

/***** TYPES ***********************************/

class SlrStation;
typedef std::shared_ptr<SlrStation> SlrStationPtr;

/***** CLASS ***********************************/

/** @brief Class for a SLR station. */
class SlrStation : public SlrPlatform
{
  Polynomial            polynomial;
  /// Rotation matrix from the global system (TRF) to the local system.
  Transform3d           global2local;
  mutable VariableList  accuracyVariableList;
  ExpressionVariablePtr accuracyExpr;

public:
  /// of offset and timeBiases
  std::vector<Time>         times;
  /// Regularized marker position in the global system (TRF)
  Vector3d                  pos;
  /// Offsets of the marker position to the System Reference Point (SRP) in the local system
  Matrix                    offset;
  /// Observed time - corrected time [seconds]
  Vector                    timeBiases;
  std::vector<std::string>  preprocessingInfos;
  std::string               disableReason;
  /// Observations at this station (for each satellite, for each pass)
  std::vector<std::vector<SlrObservationPtr>> observations;

  /** @brief Constructor. */
  SlrStation(const Platform &platform, const std::vector<Time> &times, const Vector3d &pos, const_MatrixSliceRef offset,
             ExpressionVariablePtr accuracyExpr, UInt interpolationDegree);

  /// Destructor.
  virtual ~SlrStation() {}

  /** @brief Gets the ID of the station in the SLR system. */
  UInt idStat() const {return id_;}

  /** @brief Disables the station completely. */
  void disable(const std::string &reason) override;

  /** @brief Gets the position of the System Reference Point (SRP) of the station in TRF. */
  Vector3d position(const Time &time) const;

  /** @brief Gets the rotation matrix from Terrestrial Reference Frame (TRF) to the local frame (north, east, up). */
  Transform3d global2localFrame(const Time &time) const;

  /** @brief Gets the @a times corrected for system time biases. */
  std::vector<Time> correctedTimes(const std::vector<Time> &times) const;

  /** @brief Gets the direction- (and other parameters) dependent standard deviation.
  * @param time Not used yet
  * @param residual 
  * @param accuracy
  * @param redundancy
  * @param laserWavelength
  * @param azimut Azimuth of the observation in the antenna frame (left-handed).
  * @param elevation Elevation of the observation in the antenna frame (left-handed).
  */
  Double accuracy(const Time &time, Double residual, Double accuracy, Double redundancy,
                  Double laserWavelength, Angle azimut, Angle elevation) const;
};

/***********************************************/

/// @}

#endif /* __GROOPS___ */
