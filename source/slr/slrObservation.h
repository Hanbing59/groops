/***********************************************/
/**
* @file slrObservation.h
*
* @brief SLR observations and observation equations.
*
* @author Torsten Mayer-Guerr
* @date 2022-04-28
*
*/
/***********************************************/

#ifndef __GROOPS_SLROBSERVATION__
#define __GROOPS_SLROBSERVATION__

/** @addtogroup slrGroup */
/// @{

/***** TYPES ***********************************/

class SlrStation;
class SlrSatellite;
class SlrObservation;
typedef std::shared_ptr<SlrObservation> SlrObservationPtr;

/***** CLASS ***********************************/

/** @brief SLR observations for one pass of a station-satellite pair. */
class SlrObservation
{
public:
  /// transmitted times of observations
  std::vector<Time> timesTrans;
  Vector            observations;    ///< original observations
  Vector            residuals;       ///< estimated postfit residuals
  Vector            redundancies;    ///< partial redundancies of the least squares adjustment
  Vector            sigmas0;         ///< expected (apriori) accuracies
  Vector            sigmas;          ///< modfied accuracies (downweighted outliers)
  Vector            laserWavelength; ///< laser wavelength

  SlrObservation() {}

  /** @brief Set standard deviations of the observations, and remove invalid ones.
   * Invalid observations include those have low elevations.
   * @param[in] station
   * @param[in] satellite
   * @param[in] rotationCrf2Trf
   * @param[in] elevationCutOff Cut-off elevation for valid observations
   */
  Bool init(const SlrStation &station, const SlrSatellite &satellite,
            const std::function<Rotary3d(const Time &time)> &rotationCrf2Trf, Angle elevationCutOff);

  /** @brief Calculates homogenized residuals */
  void setHomogenizedResiduals(const_MatrixSliceRef residuals, const_MatrixSliceRef redundancy);
};

/***** CLASS ***********************************/

/** @brief Reduced SLR observations (obs - computed) and design matrix for one pass of a station-satellite pair. */
class SlrObservationEquation
{
public:
  /// Indexes of different (group of) parameters in the design matrix
  enum {idxPosStat = 0,  // x,y,z (CRF)
        idxPosSat  = 3,  // x,y,z (CRF)
        idxTime    = 6,  // transmit time offset
        idxRange   = 7}; // One way range

  /// Supported SLR observation types
  enum Type {RANGE, DIRECTIONS};

  /// The type of SLR observations
  Type  type;
  /// The affiliated station
  const SlrStation   *station;
  /// The affiliated satellite
  const SlrSatellite *satellite;

  /// Weighted (reduced) observations (with 1/sigma)
  Vector l;
  /// The design matrix
  Matrix A;
  /// A-priori accuracies of observations
  Vector sigmas;
  Vector sigmas0;

  /// Indexes of observations, starting from 0
  std::vector<UInt>     index;
  /// Number of observations per epoch, always 1 for the observation type of @a RANGE 
  std::vector<UInt>     count;
  /// Positions of the station at transmit times in CRF
  std::vector<Vector3d> posStat;
  /// Positions of the satellite at bounce times in CRF
  std::vector<Vector3d> posSat;
  /// Transmit times at the station of the observations, corrected for system time biases
  std::vector<Time>     timesStat;
  /// Bounce times at the satellite of the observations
  std::vector<Time>     timesSat;
  /// Azimuths of the observations at the station
  std::vector<Angle>    azimutStat;
  /// Elevations of the observations at the station
  std::vector<Angle>    elevationStat;
  Vector                laserWavelength;

  /** @brief Constructor */
  SlrObservationEquation() : station(nullptr), satellite(nullptr) {}

  /** @brief Computes the reduced observations (OmCs) and forms the design matrix for one pass of SLR observations between a station and a satellite.
   * @param[in] observation One pass of SLR observations of the station-satellite pair.
   * @param[in] station the SLR station
   * @param[in] satellite the SLR satellite
   * @param[in] rotationCrf2Trf
   * @param[in] reduceModels
   * @param[in] homogenize
  */
  void compute(const SlrObservation &observation, const SlrStation &station, const SlrSatellite &satellite,
               const std::function<Rotary3d(const Time &time)> &rotationCrf2Trf,
               const std::function<void(SlrObservationEquation &eqn)> &reduceModels, Bool homogenize);
};

/***********************************************/

/// @}

#endif /* __GROOPS___ */
