/***********************************************/
/**
* @file gnssObservation.h
*
* @brief GNSS code & phase observations.
*
* @author Torsten Mayer-Guerr
* @author Sebastian Strasser
* @date 2012-04-18
*
*/
/***********************************************/

#ifndef __GROOPS_GNSSOBSERVATION__
#define __GROOPS_GNSSOBSERVATION__

#include "base/gnssType.h"

/** @addtogroup gnssGroup */
/// @{

/***** TYPES ***********************************/

class GnssTrack;
class GnssReceiver;
class GnssTransmitter;

/***** CLASS ***********************************/

/** @brief Class for one single GNSS observation. */
class GnssSingleObservation
{
public:
  /// GNSS measurement types (phases, pseudo ranges, ...)
  GnssType type;
  /// original observations
  Double   observation;
  /// estimated postfit residuals
  Double   residuals;
  /// partial redundancies of the least squares adjustment
  Double   redundancy;
  /// a priori accuracy, including signal direction-dependent STD
  Double   sigma0;
  /// modified accuracy (especially for downweighted observations), initially set to the same value as sigma0
  Double   sigma;

  /** @brief Constructor. */
  GnssSingleObservation() {}

  /** @brief Constructor. */
  GnssSingleObservation(GnssType _type, Double _observation=0., Double _residuals=0., Double _redundancy=0., Double _sigma0=0., Double _sigma=0.)
    : type(_type), observation(_observation), residuals(_residuals), redundancy(_redundancy), sigma0(_sigma0), sigma(_sigma) {}
};

/***** CLASS ***********************************/

/** @brief Class for GNSS observations between one receiver and one transmitter at one epoch. */
class GnssObservation
{
  /// List of GNSS observations
  std::vector<GnssSingleObservation> obs;

public:
  /// Bit mask for used GNSS observation groups
  typedef UInt Group;
  /// use code observations only
  constexpr static Group RANGE = 1<<0;
  /// use phase observations only
  constexpr static Group PHASE = 1<<1;
  /// use ionosphere observations only
  constexpr static Group IONO  = 1<<2;
  /// use all available observations
  constexpr static Group ALL   = ~0;

  /// the phase track to which this observation belongs
  GnssTrack *track;
  /// the total STEC
  Double     STEC;
  /// estimated part of the STEC
  Double     dSTEC;
  /// standard deviation of the STEC
  Double     sigmaSTEC;


  /** @brief Constructor. */
  GnssObservation() : track(nullptr), STEC(0.), dSTEC(0.), sigmaSTEC(0.) {}

  /** @brief Returns the size of the observations list. */
  UInt size() const                                      {return obs.size();}

  /** @brief Returns the reference to an GNSS observation indexed by its position in the observations list. */
  GnssSingleObservation       &at(UInt idType)           {return obs.at(idType);}

  /** @brief Returns the const reference to an GNSS observation indexed by its position in the observations list. */
  const GnssSingleObservation &at(UInt idType) const     {return obs.at(idType);}

  /** @brief Returns the reference to an GNSS observation indexed by its GNSS observation type. */
  GnssSingleObservation       &at(GnssType type)         {return obs.at(index(type));}

  /** @brief Returns the const reference to an GNSS observation indexed by its GNSS observation type. */
  const GnssSingleObservation &at(GnssType type) const   {return obs.at(index(type));}

  /** @brief Returns the position of the GNSS observation of a given type in the observations list. */
  UInt index(GnssType type) const;

  /** @brief Adds one GNSS observation to the observations list. */
  void push_back(const GnssSingleObservation &singleObs) {obs.push_back(singleObs);}

  /** @brief Removes one GNSS observation indexed by its position in the observations list. */
  void erase(UInt idType)                                {obs.erase(obs.begin()+idType); obs.shrink_to_fit();}
  void shrink_to_fit()                                   {obs.shrink_to_fit();}

  /** @brief Sorts the GNSS observations in the observations list by their GNSS observation type. */
  void sort()                                            {std::sort(obs.begin(), obs.end(), [](auto &a, auto &b){return a.type < b.type;});}

  /**
   * @brief Initializes observations for a given receiver and transmitter at one epoch,
   * performs checks to remove invalid observations and applies phase wind-up correction.
   * @param receiver The GNSS receiver.
   * @param transmitter The GNSS transmitter.
   * @param rotationCrf2Trf A function to compute the rotation from CRF to TRF.
   * @param idEpoch The epoch index for which to initialize the observation.
   * @param elevationCutOff The elevation cut-off angle at the receiver for valid observations.
   * @param phaseWindupOld The previous phase windup for input and the updated phase windup for output.
   * @return True if the initialization is successful, false otherwise.
   */
  Bool init(const GnssReceiver &receiver, const GnssTransmitter &transmitter, const std::function<Rotary3d(const Time &time)> &rotationCrf2Trf,
            UInt idEpoch, Angle elevationCutOff, Double &phaseWindupOld);

  /**
   * @brief Returns true if required GNSS observation types specified by @a group are available.
   * @param group The GNSS observation group to check for availability.
   * @param types List of available GNSS observation types.
   * @return True if the required GNSS observation types are available, false otherwise.
   * @note If @a group contains CODE or PHASE, code or phase observations from dual-frequency signals are required, respectively.
   */
  Bool observationList        (Group group, std::vector<GnssType> &types) const;
  void setHomogenizedResiduals(const std::vector<GnssType> &types, const_MatrixSliceRef residuals, const_MatrixSliceRef redundancy);
  void updateParameter        (const_MatrixSliceRef x, const_MatrixSliceRef covariance=Matrix());
};

/***** CLASS ***********************************/

/** @brief Observation equations (reduced observations) for one GNSS receiver to one GNSS transmitter at one epoch. */
class GnssObservationEquation
{
public:
  enum {idxPosRecv    = 0, // x,y,z (CRF)
        idxClockRecv  = 3,
        idxPosTrans   = 4, // x,y,z (CRF)
        idxClockTrans = 7,
        idxRange      = 8,
        idxSTEC       = 9,
        idxUnit       = 10};

  /// index of the epoch
  UInt  idEpoch;
  // the associated track
  const GnssTrack       *track;
  const GnssReceiver    *receiver;
  const GnssTransmitter *transmitter;

  /// observed types, may include composite signal types
  std::vector<GnssType> types;
  /// list of basic transmitted signals constituting the composite signal types in @a types
  std::vector<GnssType> typesTransmitted;
  /// from eliminated group parameters
  UInt   rankDeficit;
  /// weighted reduced observations, i.e. scaled by 1/sigma
  Vector l;
  /// modified accuracies, initially set to the same value as sigma0
  Vector sigma;
  /// a priori accuracies
  Vector sigma0;

  // design matrix
  Matrix A;      ///< columns: dl/dx, dl/dy, dl/dz, dl/dClock, unit matrix, transformation matrix (typesTransmitter->types)
  Matrix B;      ///< ionosphere, ...

  // approximate values (Taylor point)
  Time     timeRecv, timeTrans;
  Vector3d posRecv,  posTrans;
  Vector3d velocityRecv, velocityTrans;
  Angle    azimutRecvLocal, elevationRecvLocal;
  Angle    azimutRecvAnt,   elevationRecvAnt;
  Angle    azimutTrans,     elevationTrans;
  Double   STEC, dSTEC;
  Double   sigmaSTEC;

  GnssObservationEquation() : idEpoch(NULLINDEX), track(nullptr), receiver(nullptr), transmitter(nullptr), rankDeficit(0), STEC(0), dSTEC(0), sigmaSTEC(0) {}

  GnssObservationEquation(const GnssObservation &observation, const GnssReceiver &receiver, const GnssTransmitter &transmitter,
                          const std::function<Rotary3d(const Time &time)> &rotationCrf2Trf, const std::function<void(GnssObservationEquation &eqn)> &reduceModels,
                          UInt idEpoch, Bool homogenize, const std::vector<GnssType> &types)
    {compute(observation, receiver, transmitter, rotationCrf2Trf, reduceModels, idEpoch, homogenize, types);}

  /**
   * @brief Computes the observation equations (reduced observations) for a given receiver and transmitter at one epoch.
   */
  void compute(const GnssObservation &observation, const GnssReceiver &receiver, const GnssTransmitter &transmitter,
               const std::function<Rotary3d(const Time &time)> &rotationCrf2Trf, const std::function<void(GnssObservationEquation &eqn)> &reduceModels,
               UInt idEpoch, Bool homogenize, const std::vector<GnssType> &types);

  void eliminateGroupParameters(Bool removeRows=TRUE);
};

/***********************************************/

/// @}

#endif /* __GROOPS___ */
