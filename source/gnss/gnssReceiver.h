/***********************************************/
/**
* @file gnssReceiver.h
*
* @brief GNSS receiver.
*
* @author Torsten Mayer-Guerr
* @author Norbert Zehentner
* @author Sebastian Strasser
* @date 2013-06-28
*
*/
/***********************************************/

#ifndef __GROOPS_GNSSRECEIVER__
#define __GROOPS_GNSSRECEIVER__

#include "files/fileInstrument.h"
#include "classes/noiseGenerator/noiseGenerator.h"
#include "gnss/gnssObservation.h"
#include "gnss/gnssTransceiver.h"
#include "gnss/gnssTransmitter.h"

/** @addtogroup gnssGroup */
/// @{

/***** TYPES ***********************************/

class GnssReceiver;
class GnssTrack;
typedef std::shared_ptr<GnssReceiver> GnssReceiverPtr;
typedef std::shared_ptr<GnssTrack>    GnssTrackPtr;

/***** CLASS ***********************************/

/** @brief Abstract class for GNSS receivers. */
class GnssReceiver : public GnssTransceiver
{
  /// Observations of the receiver organized in a continuous memory block
  std::vector<GnssObservation> obsMem;
  /// Observations of the receiver organized in nested vectors (for each epoch, for each transmitter)
  std::vector<std::vector<GnssObservation*>> observations_;

  /** @brief Copy observations organized in nested vectors to a continuous memory block,
   * which seems to improve the performance.
  */
  void copyObservations2ContinuousMemoryBlock();

public:
  Bool                      isMyRank_;
  Bool                      isEarthFixed_;
  /// Processing epochs
  std::vector<Time>         times;
  /// Clock errors
  std::vector<Double>       clk;
  /// Regularized marker position in global system
  std::vector<Vector3d>     pos, vel;
  /// pos to ARF in local system
  std::vector<Vector3d>     offset;
  std::vector<Transform3d>  global2local;
  std::vector<Transform3d>  global2antenna;
  /// GNSS tracks
  std::vector<GnssTrackPtr> tracks;
  /// Median of the observation sampling interval in seconds (for all epochs)
  Double                    observationSampling;
  Bool                      integerAmbiguities;
  /// Factor to account for half-wavelength observations, i.e., collected by codeless squaring techniques.
  Double                    wavelengthFactor;
  /// Preprocessing informations.
  std::vector<std::string>  preprocessingInfos;
  /// Reason for disabling the receiver, which is empty if the receiver is enabled.
  std::string               disableReason;

  /** @brief Constructor. */
  GnssReceiver(Bool isMyRank, Bool isEarthFixed, const Platform &platform,
               GnssAntennaDefinition::NoPatternFoundAction noPatternFoundAction, const Vector &useableEpochs,
               Bool integerAmbiguities, Double wavelengthFactor);

  /** @brief Destructor. */
 ~GnssReceiver() {}

  /** @brief Gets ID of the receiver. */
  UInt  idRecv() const {return id_;}

  /** @brief Disables a given epoch for a given reason. */
  void disable(UInt idEpoch, const std::string &reason) override;

  /** @brief Disable the receiver for a given reason, which will be recorded in the member variable @a disableReason. */
  void disable(const std::string &reason) override;

  Bool isMyRank() const {return isMyRank_;}

  /** @brief Gets the clock error at a given epoch.
  * error = clock time - system time [s] */
  Double clockError(UInt idEpoch) const {return clk.at(idEpoch);}

  /** @brief Accumulates the clock error by a given amount.
  * error = clock time - system time [s] */
  void updateClockError(UInt idEpoch, Double deltaClock) {clk.at(idEpoch) += deltaClock;}

  /** @brief Applies the clock error correction of the time tag for a given epoch.
  * (receiver clock time - clock error) */
  Time timeCorrected(UInt idEpoch) const {return times.at(idEpoch) - seconds2time(clockError(idEpoch));}

  /** @brief Whether this is an earth-fixed receiver, e.g. on a ground GNSS station. */
  Bool isEarthFixed() const {return isEarthFixed_;}

  /** @brief Gets the position of the antenna reference point in TRF or CRF. */
  Vector3d position(UInt idEpoch) const {return pos.at(idEpoch) + global2local.at(idEpoch).inverseTransform(offset.at(idEpoch));}

  /** @brief Gets the velocity in TRF or CRF [m/s]. */
  Vector3d velocity(UInt idEpoch) const {return vel.at(idEpoch);}

  /** @brief Rotation from terrestrial reference frame (TRF) or celestial reference frame (CRF) to local horizont system (north, east, up). */
  Transform3d global2localFrame(UInt idEpoch) const {return global2local.at(idEpoch);}

  /** @brief Rotation from terrestrial reference frame (TRF) or celestial reference frame (CRF) to left-handed antenna system (usually north, east, up). */
  Transform3d global2antennaFrame(UInt idEpoch) const {return global2antenna.at(idEpoch);}

  /** @brief Transformation matrix for observed (composed) signal types from original transmitted signal types.
   * E.g. C2DG = C1CG - C1WG + C2WG.
   * Returns the @a typesTrans and the transformation matrix @a A (dimension: types.size() times typesTrans.size()). */
  static void signalCompositionDefault(const std::vector<GnssType> &types, std::vector<GnssType> &typesTrans, Matrix &A);

  /** @brief Transformation matrix for observed (composed) signal types from original transmitted signal types.
  * E.g. C2DG = C1CG - C1WG + C2WG.
  * Returns the @a typesTrans and the transformation matrix @a A (dimension: types.size() times typesTrans.size()). */
  virtual void signalComposition(UInt /*idEpoch*/, const std::vector<GnssType> &types, std::vector<GnssType> &typesTrans, Matrix &A) const;

  /** @brief Gets all observations between the receiver and a given transmitter at a given epoch. */
  GnssObservation *observation(UInt idTrans, UInt idEpoch) const;

  /** @brief Returns the number of observed epochs. */
  UInt idEpochSize() const {return observations_.size();}

  /** @brief Returns the number of observed transmitters at a given epoch. */
  UInt idTransmitterSize(UInt idEpoch) const {return (idEpoch < observations_.size()) ? observations_.at(idEpoch).size() : 0;}

  /** @brief Deletes the observations for a given transmitter and epoch.
   * If no transmitters remain for the epoch after the deletion, the epoch will be disabled.
  */
  void deleteObservation(UInt idTrans, UInt idEpoch);

  /** @brief Deletes all empty tracks. */
  void deleteEmptyTracks();

  // Preprocessing
  // -------------
  /** @brief Compute observation equations (reduced observations). */
  class ObservationEquationList
  {
    std::vector<std::vector<std::unique_ptr<GnssObservationEquation>>> eqn;

  public:
    ObservationEquationList() {}
    ObservationEquationList(const GnssReceiver &receiver, const std::vector<GnssTransmitterPtr> &transmitters,
                            const std::function<Rotary3d(const Time &time)> &rotationCrf2Trf,
                            const std::function<void(GnssObservationEquation &eqn)> &reduceModels, GnssObservation::Group group);

    GnssObservationEquation *operator()(UInt idTrans, UInt idEpoch) const;
  };

  /**
   * @brief Adds a preprocessing information message into the preprocessing informations list @a preprocessingInfos.
   * @param info The information message to add.
   * @param countEpochs The number of epochs.
   * @param countObservations The number of observations.
   * @param countTracks The number of tracks.
   */
  void preprocessingInfo(const std::string &info, UInt countEpochs=NULLINDEX, UInt countObservations=NULLINDEX, UInt countTracks=NULLINDEX);

  /**
   * @brief Reads observations from a file.
   * @param fileName The name of the file to read observations from.
   * @param transmitters A vector of shared pointers to GNSS transmitters.
   * @param rotationCrf2Trf A function that returns the rotation from celestial reference frame (CRF) to terrestrial reference frame (TRF) at a given time.
   * @param timeMargin The time margin for reading observations.
   * @param elevationCutOff The elevation cut-off angle for observations.
   * @param useType A vector of GNSS signal types to use.
   * @param ignoreType A vector of GNSS signal types to ignore.
   * @param group The group of observations to read.
   * @note The member variable @a times must be set before calling this function.
   * @note Receiver and Transmitter positions, orientations, etc. must be initialized beforehand.
   * @note Observations that don't match the types from receiver and transmitter definition will be deleted.
   */
  void readObservations(const FileName &fileName, const std::vector<GnssTransmitterPtr> &transmitters, const std::function<Rotary3d(const Time &time)> &rotationCrf2Trf,
                        const Time &timeMargin, Angle elevationCutOff, const std::vector<GnssType> &useType, const std::vector<GnssType> &ignoreType, GnssObservation::Group group);

  /**
   * @brief Simulates observations.
   * @note The member variable @a times must be set.
   * @note Receiver and Transmitter positions, orientations, ... must be initialized beforehand.
   * @note Observations that don't match the types from receiver and transmitter definition will be deleted.
   */
  void simulateZeroObservations(const std::vector<GnssType> &types,
                                const std::vector<GnssTransmitterPtr> &transmitters,
                                const std::function<Rotary3d(const Time &time)> &rotationCrf2Trf, Angle elevationCutOff,
                                const std::vector<GnssType> &useType, const std::vector<GnssType> &ignoreType, GnssObservation::Group group);

  /**
   * @brief Simulates observations.
   * @note The member variable @a times must be set.
   * @note Receiver and Transmitter positions, orientations, ... must be initialized beforehand.
   * @note Observations that don't match the types from receiver and transmitter definition will be deleted.
   */
  void simulateObservations(NoiseGeneratorPtr noiseClock, NoiseGeneratorPtr noiseObs,
                            const std::vector<GnssTransmitterPtr> &transmitters,
                            const std::function<Rotary3d(const Time &time)> &rotationCrf2Trf,
                            const std::function<void(GnssObservationEquation &eqn)> &reduceModels,
                            UInt minObsCountPerTrack, Angle elevationTrackMinimum, GnssObservation::Group group);

  /**
   * @brief Estimates coarse receiver clock errors from a code-based Precise Point Positioning (PPP) solution.
   * @param transmitters A vector of shared pointers to GNSS transmitters.
   * @param rotationCrf2Trf A function that returns the rotation from celestial reference frame (CRF) to terrestrial reference frame (TRF) at a given time.
   * @param reduceModels A function that reduces the observation equations.
   * @param huber The Huber loss parameter.
   * @param huberPower The power for the Huber loss function.
   * @param estimateKinematicPosition If TRUE, the receiver position is estimated at each epoch, otherwise it is estimated once for all epochs.
   * @return the estimated positions
   */
  std::vector<Vector3d> estimateInitialClockErrorFromCodeObservations(const std::vector<GnssTransmitterPtr> &transmitters,
                                                                      const std::function<Rotary3d(const Time &time)> &rotationCrf2Trf,
                                                                      const std::function<void(GnssObservationEquation &eqn)> &reduceModels,
                                                                      Double huber, Double huberPower, Bool estimateKinematicPosition);

  /**
   * @brief Disables epochs if reduced code observations @p eqn exceed @p threshold (e.g. 100+ km for code cycle slips).
   * @param eqn The observation equations containing the reduced observations.
   * @param threshold The threshold for detecting gross code outliers.
   * @param outlierRatio The ratio of observed satellites that must exceed the threshold to disable the epoch.
   * @note Epochs for which a ratio of @p outlierRatio or more of the observed satellites have gross code outliers will be disabled.
   */
  void disableEpochsWithGrossCodeObservationOutliers(ObservationEquationList &eqn, Double threshold, Double outlierRatio=0.5);

  /**
   * @brief Creates tracks with continuously identical phase observations.
   * @param transmitters A vector of shared pointers to GNSS transmitters.
   * @param minObsCountPerTrack The minimum number of observations per track.
   * @param extraTypes A vector of extra GNSS signal types to include in the tracks.
   *
   * @note Tracks may contain gaps but must contain observations of at least @p minObsCountPerTrack epochs.
   * @note Extra types are included (e.g. L5*G), but tracks must have at least two others phases at different frequencies.
   */
  void createTracks(const std::vector<GnssTransmitterPtr> &transmitters, UInt minObsCountPerTrack, const std::vector<GnssType> &extraTypes={});

  /** @brief Deletes a track and all its related observations. */
  void deleteTrack(UInt idTrack);

  /** @brief Removes tracks that never exceed @p minElevation (in radian). */
  void removeLowElevationTracks(ObservationEquationList &eqn, Angle minElevation);

  /** @brief Split a @p track at @p idEpochSplit into two new tracks.
  * Shortens the old track and returns the new track. Updates observation and observation equation @p eqn track assignments. */
  GnssTrackPtr splitTrack(ObservationEquationList &eqn, GnssTrackPtr track, UInt idEpochSplit);

  /**
   * @brief Computes Melbourne-Wuebbena-like linear combinations for a given track.
   * @param eqnList The list of observation equations.
   * @param track The track for which to compute combinations.
   * @param extraTypes The extra GNSS signal types to include.
   * @param typesPhase The phase types found in the track.
   * @param idEpochs The epochs in the track.
   * @param combinations The resulting combinations.
   * @param cycles2tecu The conversion factor from cycles to TECU.
   */
  void linearCombinations(ObservationEquationList &eqnList, GnssTrackPtr track, const std::vector<GnssType> &extraTypes,
                          std::vector<GnssType> &typesPhase, std::vector<UInt> &idEpochs, Matrix &combinations, Double &cycles2tecu) const;

  /** @brief Computes range and TEC from phase observations.
   * @param eqnList The list of observation equations.
   * @param idTrans The ID of the transmitter.
   * @param idEpochs The epochs for which to compute range and TEC.
   * @param typesPhase The phase types found in the track.
   * @param range The resulting range.
   * @param tec The resulting TEC.
   */
  void rangeAndTec(ObservationEquationList &eqnList, UInt idTrans, const std::vector<UInt> &idEpochs,
                   const std::vector<GnssType> &typesPhase, Vector &range, Vector &tec) const;

  /**
   * @brief Writes tracks to a file.
   * @param fileName The name of the file to write tracks to.
   * @param eqnList The list of observation equations.
   * @param extraTypes The extra GNSS signal types to include.
   *
   * @note The tracks are written in a format that can be read by the @a readTracks function.
   */
  void writeTracks(const FileName &fileName, ObservationEquationList &eqnList, const std::vector<GnssType> &extraTypes) const;

  /**
   * @brief Splits tracks at detected cycle slips based on all Melbourne-Wuebbena like combinations.
   * Solves the total variation regularized least-squares problem (total variation denoising).
   * In a second step, smoothness of TEC is evaluated using a moving window peak/outlier detection based on AR model residuals.
   *
   * @param eqnList Observation equations (reduced observations).
   * @param minObsCountPerTrack Removes tracks that don't have enough observations.
   * @param lambda Regularization parameter (@see @a totalVariationDenoising) (e.g. @p lambda = 5 for GPS ground stations).
   * @param windowSize Size of the moving window used for the TEC smoothness evaluation. If 0, TEC is not analyzed.
   * @param tecSigmaFactor Factor applied to moving standard deviation of AR model residuals to determine threshold for peak/outlier detection.
   * @param extraTypes GPS L5 observations are handled separately due to temporal changing bias.*/
  void cycleSlipsDetection(ObservationEquationList &eqnList, UInt minObsCountPerTrack, Double lambda, UInt windowSize, Double tecSigmaFactor, const std::vector<GnssType> &extraTypes={});

  /**
   * @brief Splits tracks at detected cycle slips based on all Melbourne-Wuebbena like combinations.
   */
  void cycleSlipsDetection(ObservationEquationList &eqnList, GnssTrackPtr track, Double lambda, UInt windowSize, Double tecSigmaFactor, const std::vector<GnssType> &extraTypes);

  /**
   * @brief Repairs cycle slip differences at same frequencies (e.g. between L1CG and L1WG), which allows to reduce the number of integer ambiguities.
   * @param eqnList Observation equations (reduced observations). */
  void cycleSlipsRepairAtSameFrequency(ObservationEquationList &eqnList);

  /**
   * @brief Track outlier detection based on robust least squares estimation.
   * @param eqn Observation equations (reduced observations).
   * @param ignoreTypes Types of observations to ignore.
   * @param huber Huber loss parameter.
   * @param huberPower Power for the Huber loss function.
   */
  void trackOutlierDetection(const ObservationEquationList &eqn, const std::vector<GnssType> &ignoreTypes, Double huber, Double huberPower);

  /**
   * @brief Total variation denoising, which solves the total variation regularized least-squares problem.
   * Laurent Condat. A Direct Algorithm for 1D Total Variation Denoising. IEEE Signal Processing Letters,
   * Institute of Electrical and Electronics Engineers, 2013, 20 (11), pp.1054-1057. DOI: 10.1109/LSP.2013.2278339.
   *
   * @param y noisy input time series.
   * @param lambda Regularization parameter.
   * @return denoised time series. */
  static Matrix totalVariationDenoising(const_MatrixSliceRef y, Double lambda);
}; //class GnssReceiver

/***** CLASS ***********************************/

class GnssAmbiguity;

/** @brief Class for GNSS tracking phases and ranges along epochs. */
class GnssTrack
{
public:
  GnssReceiver         *receiver;
  GnssTransmitter      *transmitter;
  UInt                  idEpochStart, idEpochEnd;
  std::vector<GnssType> types;
  GnssAmbiguity        *ambiguity;

  /** @brief Constructor. */
  GnssTrack(GnssReceiver *_receiver, GnssTransmitter *_transmitter, UInt _idEpochStart, UInt _idEpochEnd, const std::vector<GnssType> &_types);

  /** @brief Destructor. */
 ~GnssTrack();

  /** @brief Returns the number of observations in the track. */
  UInt countObservations() const;

  /**
   * @brief Removes ambiguities from observations in the track.
   * @param type Types of observations for which to remove ambiguities.
   * @param value Values of the ambiguities to remove. */
  void removeAmbiguitiesFromObservations(const std::vector<GnssType> &type, const std::vector<Double> &value);
};

/***** CLASS ***********************************/

/** @brief Abstract class for GNSS ambiguities. */
class GnssAmbiguity
{
public:
  /// Pointer to the track that owns this ambiguity. The ambiguity will be deleted when the track is destructed.
  GnssTrack *track;

  explicit GnssAmbiguity(GnssTrack *track) : track(track) {track->ambiguity = this;}
  virtual ~GnssAmbiguity() {}
  virtual Vector ambiguities(const std::vector<GnssType> &type) const = 0;
};

/***********************************************/

/// @}

#endif /* __GROOPS___ */
