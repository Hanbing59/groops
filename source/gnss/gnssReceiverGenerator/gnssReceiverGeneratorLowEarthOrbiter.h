/***********************************************/
/**
* @file gnssReceiverGeneratorLowEarthOrbiter.h
*
* @brief GNSS receiver generator for Low Earth Orbiter (LEO).
*
* @author Torsten Mayer-Guerr
* @author Sebastian Strasser
* @date 2021-02-25
*
*/
/***********************************************/

#ifndef __GROOPS_GNSSRECEIVERGENERATORLOWEARTHORBITER__
#define __GROOPS_GNSSRECEIVERGENERATORLOWEARTHORBITER__

// Latex documentation
#ifdef DOCSTRING_GnssReceiverGenerator
static const char *docstringGnssReceiverGeneratorLowEarthOrbiter = R"(
\subsection{LowEarthOrbiter}\label{gnssReceiverGeneratorType:lowEarthOrbiter}
A single low-Earth orbiting (LEO) satellite with an onboard GNSS receiver.
An apriori orbit is needed as \configFile{inputfileOrbit}{instrument}.
Attitude data must be provided via \configFile{inputfileStarCamera}{instrument}.
If no attitude data is available from the satellite operator,
the star camera data can be simulated by using \program{SimulateStarCamera}.
)";
#endif

/***********************************************/

#include "config/config.h"
#include "classes/noiseGenerator/noiseGenerator.h"
#include "gnss/gnssReceiverGenerator/gnssReceiverGenerator.h"

/***** CLASS ***********************************/

/** @brief GNSS receiver generator for Low Earth Orbiter (LEO).
* @ingroup gnssReceiverGeneratorGroup
* @see GnssReceiverGenerator */
class GnssReceiverGeneratorLowEarthOrbiter : public GnssReceiverGeneratorBase
{
  FileName              fileNameStationInfo;
  FileName              fileNameAntennaDef;
  FileName              fileNameReceiverDef;
  FileName              fileNameAccuracyDef;
  FileName              fileNameObs;
  FileName              fileNameOrbit;
  FileName              fileNameStarCamera;
  ExpressionVariablePtr exprSigmaPhase, exprSigmaCode;
  Bool                  integerAmbiguities;
  Double                wavelengthFactor;
  Angle                 elevationCutOff;
  std::vector<GnssType> useType, ignoreType;
  GnssAntennaDefinition::NoPatternFoundAction noPatternFoundAction;
  /// Whether to print preprocessing statistics
  Bool                  printInfo;
  Double                huber, huberPower;
  Double                codeMaxPosDiff;
  UInt                  minObsCountPerTrack;
  Double                denoisingLambda;
  UInt                  tecWindowSize;
  Double                tecSigmaFactor;
  FileName              fileNameTrackBefore, fileNameTrackAfter;
  GnssReceiverPtr       recv;

public:
  /** @brief Constructor. */
  GnssReceiverGeneratorLowEarthOrbiter(Config &config);

  /** @brief Initializes the GNSS receiver generator for Low Earth Orbiter (LEO). */
  void init(std::vector<GnssType> simulationTypes, const std::vector<Time> &times, const Time &timeMargin,
            const std::vector<GnssTransmitterPtr> &transmitters, EarthRotationPtr earthRotation,
            Parallel::CommunicatorPtr comm, std::vector<GnssReceiverPtr> &receivers) override;

  void preprocessing(Gnss *gnss, Parallel::CommunicatorPtr comm) override;

  void simulation(NoiseGeneratorPtr noiseClock, NoiseGeneratorPtr noiseObs,
                  Gnss *gnss, Parallel::CommunicatorPtr comm) override;
};

/***********************************************/

#endif
