/***********************************************/
/**
* @file gnss.h
*
* @brief global navigation satellite system classes.
*
* @author Torsten Mayer-Guerr
* @date 2010-08-03
*
*/
/***********************************************/

#ifndef __GROOPS_GNSS__
#define __GROOPS_GNSS__

#include "parallel/parallel.h"
#include "base/gnssType.h"
#include "base/parameterName.h"
#include "classes/noiseGenerator/noiseGenerator.h"
#include "gnss/gnssObservation.h"
#include "gnss/gnssDesignMatrix.h"
#include "gnss/gnssTransmitter.h"
#include "gnss/gnssReceiver.h"
#include "gnss/gnssNormalEquationInfo.h"

/** @addtogroup gnssGroup */
/// @{

/***** TYPES ***********************************/

class Gnss;
typedef std::shared_ptr<Gnss> GnssPtr;

class MatrixDistributed;
class GnssTransmitterGenerator;
class GnssReceiverGenerator;
class GnssParametrization;
class EarthRotation;
class PlatformSelector;
typedef std::shared_ptr<GnssTransmitterGenerator> GnssTransmitterGeneratorPtr;
typedef std::shared_ptr<GnssReceiverGenerator>    GnssReceiverGeneratorPtr;
typedef std::shared_ptr<GnssParametrization>      GnssParametrizationPtr;
typedef std::shared_ptr<EarthRotation>            EarthRotationPtr;
typedef std::shared_ptr<PlatformSelector>         PlatformSelectorPtr;

/***** CLASS ***********************************/

class Gnss
{
public:
  /// Epochs
  std::vector<Time>               times;
  /// GNSS transmitters (satellites)
  std::vector<GnssTransmitterPtr> transmitters;
  /// GNSS receivers (stations & LEOs)
  std::vector<GnssReceiverPtr>    receivers;
  /// Parameterization of GNSS observations
  GnssParametrizationPtr          parametrization;
  std::function<void(GnssObservationEquation &eqn)> funcReduceModels;
  std::function<Rotary3d(const Time &time)>         funcRotationCrf2Trf;
  /// Earth Orientation Parameters (EOP): xp, yp, sp, deltaUT, LOD, X, Y, S
  Matrix                          eop;
  /// List of GNSS types used for each receiver and transmitter
  std::vector<std::vector<std::vector<GnssType>>> typesRecvTrans;

  /**
   * @brief Initializes the GNSS object.
   * @param simulationTypes List of GNSS types to be simulated.
   * @param times List of epochs.
   * @param timeMargin Time margin for the epochs.
   * @param transmitterGenerator Generator of GNSS transmitters.
   * @param receiverGenerator Generator of GNSS receivers.
   * @param earthRotation Earth rotation model.
   * @param parametrization GNSS parametrization.
   * @param comm Parallel communicator.
   */
  void init(std::vector<GnssType> simulationTypes, const std::vector<Time> &times, const Time &timeMargin,
            GnssTransmitterGeneratorPtr transmitterGenerator, GnssReceiverGeneratorPtr receiverGenerator,
            EarthRotationPtr earthRotation, GnssParametrizationPtr parametrization, Parallel::CommunicatorPtr comm);

  /**
   * @brief Computes the rotation matrix from the inertial system (CRF) to the earth fixed system (TRF) at a given time.
   * @param time Time in GPS time system.
   * @return Rotary matrix (3x3).
   */
  Rotary3d rotationCrf2Trf(const Time &time) const;

  /**
   * @brief Synchronizes the transceivers.
   * @param comm Parallel communicator.
   */
  void     synchronizeTransceivers(Parallel::CommunicatorPtr comm);

  /**
   * @brief Initializes the parameter vector.
   * @param normalEquationInfo Normal equation information.
   */
  void   initParameter            (GnssNormalEquationInfo &normalEquationInfo);

  /**
   * @brief Computes the apriori parameter vector.
   * @param normalEquationInfo Normal equation information.
   * @return Apriori parameter vector.
   */
  Vector aprioriParameter         (const GnssNormalEquationInfo &normalEquationInfo) const;

  /**
   * @brief Sets up the basic GNSS observation equations for a given receiver, transmitter, and epoch.
   * @param normalEquationInfo Normal equation information.
   * @param idRecv Receiver ID.
   * @param idTrans Transmitter ID.
   * @param idEpoch Epoch ID.
   * @param eqn Output GNSS observation equation.
   * @return TRUE if the observation equations were successfully set up, FALSE otherwise.
   */
  Bool   basicObservationEquations(const GnssNormalEquationInfo &normalEquationInfo, UInt idRecv, UInt idTrans, UInt idEpoch, GnssObservationEquation &eqn) const;

  /**
   * @brief Builds the design matrix for a given GNSS observation equation.
   * @param normalEquationInfo Normal equation information.
   * @param eqn GNSS observation equation.
   * @param A Output design matrix.
   */
  void   designMatrix             (const GnssNormalEquationInfo &normalEquationInfo, const GnssObservationEquation &eqn, GnssDesignMatrix &A) const;

  /**
   * @brief Applies the constraints for a specific epoch.
   * @param normalEquationInfo Normal equation information.
   * @param idEpoch Epoch ID.
   * @param normals Output distributed normal equations.
   * @param n Output list of normal matrices.
   * @param lPl Output weighted sum of squared residuals.
   * @param obsCount Output number of observations.
  */
  void   constraintsEpoch         (const GnssNormalEquationInfo &normalEquationInfo, UInt idEpoch, MatrixDistributed &normals, std::vector<Matrix> &n, Double &lPl, UInt &obsCount) const;

  /**
   * @brief Applies the constraints.
   * @param normalEquationInfo Normal equation information.
   * @param normals Output distributed normal equations.
   * @param n Output list of normal matrices.
   * @param lPl Output weighted sum of squared residuals.
   * @param obsCount Output number of observations.
   */
  void   constraints              (const GnssNormalEquationInfo &normalEquationInfo, MatrixDistributed &normals, std::vector<Matrix> &n, Double &lPl, UInt &obsCount) const;

  /**
   * @brief Resolves the integer ambiguities.
   * @param normalEquationInfo Normal equation information.
   * @param normals Distributed normal equations.
   * @param n List of normal matrices.
   * @param lPl Weighted sum of squared residuals.
   * @param obsCount Number of observations.
   * @param selectedTransmitters List of selected transmitters.
   * @param selectedReceivers List of selected receivers.
   * @param searchInteger Function to search for integer ambiguities.
   * @return The maximum change in the parameters after ambiguity resolution.
   */
  Double ambiguityResolve         (const GnssNormalEquationInfo &normalEquationInfo, MatrixDistributed &normals, std::vector<Matrix> &n, Double &lPl, UInt &obsCount,
                                   const std::vector<Byte> &selectedTransmitters, const std::vector<Byte> &selectedReceivers,
                                   const std::function<Vector(const_MatrixSliceRef xFloat, MatrixSliceRef W, const_MatrixSliceRef d, Vector &xInt, Double &sigma)> &searchInteger);

  /**
   * @brief Updates the parameter vector.
   * @param normalEquationInfo Normal equation information.
   * @param x Float solution vector.
   * @param Wz Weight matrix of the observations.
   * @return The maximum change in the parameters after the update.
   */
  Double updateParameter          (const GnssNormalEquationInfo &normalEquationInfo, const_MatrixSliceRef x, const_MatrixSliceRef Wz);

  /**
   * @brief Updates the covariance matrix.
   * @param normalEquationInfo Normal equation information.
   * @param covariance Distributed covariance matrix.
   */
  void   updateCovariance         (const GnssNormalEquationInfo &normalEquationInfo, const MatrixDistributed &covariance);

  /**
   * @brief Writes the results to a file.
   * @param normalEquationInfo Normal equation information.
   * @param suffix Suffix for the output file.
   */
  void   writeResults             (const GnssNormalEquationInfo &normalEquationInfo, const std::string &suffix="");

  /**
   * @brief Sorts the list of used GNSS observation types with a given mask.
   * @param mask Mask for filtering observation types.
   * @return List of sorted used GNSS observation types.
   */
  std::vector<GnssType> types(const GnssType mask=GnssType::ALL) const;

  /**
   * @brief Selects the GNSS transmitters based on the platform selector.
   * @param selector Platform selector.
   * @return List of selected transmitters.
   */
  std::vector<Byte> selectTransmitters(PlatformSelectorPtr selector);

  /**
   * @brief Selects the GNSS receivers based on the platform selector.
   * @param selector Platform selector.
   * @return List of selected receivers.
   */
  std::vector<Byte> selectReceivers(PlatformSelectorPtr selector);

  /**
   * @brief Information about parameter changes.
   */
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

    /**
     * @brief Updates the parameter change information.
     * @param change Change in the parameter.
     * @return True if the update was successful.
     */
    Bool update(Double change);

    /**
     * @brief Synchronizes and prints the parameter change information.
     * @param comm Communication object.
     * @param convertToMeter Conversion factor to meters.
     * @param maxChangeTotal Reference to the total maximum change.
     */
    void synchronizeAndPrint(Parallel::CommunicatorPtr comm, Double convertToMeter, Double &maxChangeTotal);
  };
};

/// @}

/***********************************************/

#endif /* __GROOPS___ */
