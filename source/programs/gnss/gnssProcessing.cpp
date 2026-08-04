/***********************************************/
/**
* @file gnssProcessing.cpp
*
* @brief GNSS/LEO satellite orbit determination, station network analysis, PPP.
*
* @author Torsten Mayer-Guerr
* @author Sebastian Strasser
* @date 2010-08-04
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
This program processes GNSS observations. It calculates the linearized observation equations,
accumulates them into a system of normal equations and solves it.

The primary use cases of this program are:
\begin{itemize}
  \item \reference{GNSS satellite orbit determination and station network analysis}{cookbook.gnssNetwork}
  \item \reference{Kinematic orbit determination of LEO satellites}{cookbook.kinematicOrbit}
  \item \reference{GNSS precise point positioning (PPP)}{cookbook.gnssPpp}
\end{itemize}

The observation epochs are defined by \configClass{timeSeries}{timeSeriesType}
and only observations at these epochs (within a \config{timeMargin}) are considered.

To calculate observation equations from the tracks, the model parameters or unknown parameters need to be
defined beforehand. These unknown parameters can be chosen arbitrarily by the user with an adequate list of defined
\configClass{parametrization}{gnssParametrizationType}.
Some of the \configClass{parametrization}{gnssParametrizationType} also include a priori models.

Lastly it is required to define the process flow of the gnssProcessing. This is accomplished
with a list of \configClass{processingSteps}{gnssProcessingStepType}.
Each step is processed consecutively. Some steps allow the selection of parameters, epochs,
or the normal equation structure, which affects all subsequent steps.
A minimal example consists of following steps:
\begin{itemize}
  \item \configClass{estimate}{gnssProcessingStepType:estimate}: iterative float solution with outlier downeighting
  \item \configClass{resolveAmbiguities}{gnssProcessingStepType:resolveAmbiguities}:
        fix ambiguities to integer and remove them from the normals
  \item \configClass{estimate}{gnssProcessingStepType:estimate}: few iteration for final outlier downweighting
  \item \configClass{writeResults}{gnssProcessingStepType:writeResults}:
        write the output files defined in \configClass{parametrization}{gnssParametrizationType}
\end{itemize}

If the program is run on multiple processes the \configClass{receiver}{gnssReceiverGeneratorType}s
(stations or LEO satellites) are distributed over the processes.

See also \program{GnssSimulateReceiver}.
)";

/***********************************************/

#include "programs/program.h"
#include "config/configRegister.h"
#include "classes/timeSeries/timeSeries.h"
#include "classes/earthRotation/earthRotation.h"
#include "gnss/gnss.h"
#include "gnss/gnssReceiverGenerator/gnssReceiverGenerator.h"
#include "gnss/gnssTransmitterGenerator/gnssTransmitterGenerator.h"
#include "gnss/gnssParametrization/gnssParametrization.h"
#include "gnss/gnssProcessingStep/gnssProcessingStep.h"

/***** CLASS ***********************************/

/** @brief GNSS/LEO satellite orbit determination, station network analysis, PPP.
* @ingroup programsGroup */
class GnssProcessing
{
public:
  void run(Config &config, Parallel::CommunicatorPtr comm);
};

GROOPS_REGISTER_PROGRAM(GnssProcessing, PARALLEL, "GNSS/LEO satellite orbit determination, station network analysis, PPP", Gnss)

/***********************************************/

void GnssProcessing::run(Config &config, Parallel::CommunicatorPtr comm)
{
  try
  {
    TimeSeriesPtr               timeSeries;
    Double                      marginSeconds;
    GnssTransmitterGeneratorPtr transmitterGenerator;
    GnssReceiverGeneratorPtr    receiverGenerator;
    GnssParametrizationPtr      gnssParametrization;
    EarthRotationPtr            earthRotation;
    GnssProcessingStepPtr       processingSteps;

    readConfig(config, "timeSeries",      timeSeries,           Config::MUSTSET,  "",    "defines observation epochs");
    readConfig(config, "timeMargin",      marginSeconds,        Config::DEFAULT,  "0.1", "[seconds] margin to consider two times identical");
    readConfig(config, "transmitter",     transmitterGenerator, Config::MUSTSET,  "",    "constellation of GNSS satellites");
    readConfig(config, "receiver",        receiverGenerator,    Config::MUSTSET,  "",    "ground station network or LEO satellite");
    readConfig(config, "earthRotation",   earthRotation,        Config::MUSTSET,  "",    "apriori earth rotation");
    readConfig(config, "parametrization", gnssParametrization,  Config::MUSTSET,  "",    "models and parameters");
    readConfig(config, "processingStep",  processingSteps,      Config::MUSTSET,  "",    "steps are processed consecutively");
    if(isCreateSchema(config)) return;

    // ============================

    // init the GNSS system
    // --------------------
    logInfo<<"Init GNSS"<<Log::endl;
    std::vector<Time> times = timeSeries->times();
    GnssPtr gnss = std::make_shared<Gnss>();
    gnss->init({}, times, seconds2time(marginSeconds), transmitterGenerator, receiverGenerator, earthRotation, gnssParametrization, comm);
    receiverGenerator->preprocessing(gnss.get(), comm);
    gnss->synchronizeTransceivers(comm);
    transmitterGenerator = nullptr;
    receiverGenerator    = nullptr;
    gnssParametrization  = nullptr;
    earthRotation        = nullptr;

    // Get sorted list of useable transmitters
    std::set<GnssType> transmitterSet;
    for(auto trans : gnss->transmitters)
      if(trans->useable())
        transmitterSet.insert(trans->PRN());
    if(!transmitterSet.size())
    {
      logWarningOnce<<times.front().dateTimeStr()<<" - "<<times.back().dateTimeStr()<<": no useable transmitters"<<Log::endl;
      return;
    }
    std::stringstream ss;
    ss<<"transmitters: "<<transmitterSet.size();
    for(auto it = transmitterSet.begin(); it != transmitterSet.end(); ++it)
      ss<<" "<<it->prnStr();
    logInfo<<ss.str()<<Log::endl;

    // Get sorted list of useable receivers
    std::set<std::string> receiverSet;
    for(auto recv : gnss->receivers)
      if(recv->useable())
        receiverSet.insert(recv->name());
    if(!receiverSet.size())
    {
      logWarningOnce<<times.front().dateTimeStr()<<" - "<<times.back().dateTimeStr()<<": no useable receivers"<<Log::endl;
      return;
    }
    ss.str("");
    ss<<"receivers   : "<<receiverSet.size();
    for(auto it = receiverSet.begin(); it != receiverSet.end(); ++it)
      ss<<" "<<*it;
    logInfo<<ss.str()<<Log::endl;

    // count observation types for each transmitter
    logInfo<<"types and numbers of observations:"<<Log::endl;
    std::vector<GnssType> allTypes = gnss->types(~(GnssType::PRN + GnssType::FREQ_NO));
    // Count only RANGE and PHASE observations
    std::vector<GnssType> types;
    for(const auto &type : allTypes)
      if((type & GnssType::TYPE) == GnssType::RANGE || (type & GnssType::TYPE) == GnssType::PHASE)
        types.push_back(type);
    std::vector<GnssType> satellites(transmitterSet.begin(), transmitterSet.end());
    Vector countTypes(satellites.size() * types.size());
    for(auto recv : gnss->receivers)
      if(recv->isMyRank())
        for(UInt idEpoch=0; idEpoch<recv->idEpochSize(); idEpoch++)
          for(UInt idTrans=0; idTrans<recv->idTransmitterSize(idEpoch); idTrans++)
          {
            auto obs = recv->observation(idTrans, idEpoch);
            if(obs)
            {
              // Get satellite PRN of this transmitter. Is the idTrans can be indexed into the satellites vector?
              GnssType satPrn = gnss->transmitters.at(idTrans)->PRN();
              UInt idSat = std::distance(satellites.begin(), std::find(satellites.begin(), satellites.end(), satPrn));
              if(idSat >= satellites.size())
                continue;

              for(UInt idType=0; idType<obs->size(); idType++)
              {
                const UInt idx = GnssType::index(types, obs->at(idType).type);
                if(idx == NULLINDEX)
                  continue;
                countTypes(idSat * types.size() + idx)++;
              }
            }
          }
    Parallel::reduceSum(countTypes, 0, comm);

    // Table header
    ss.str("");
    ss<<"    ";
    for(const auto &type : types)
      ss<<"     "<<type.str().substr(0,4);
    logInfo<<ss.str()<<Log::endl;
    // Total number of observations for each satellite
    std::vector<UInt> colSum(types.size(), 0);
    // Print each row of the table
    for(UInt idSat=0; idSat<satellites.size(); idSat++)
    {
      ss.str("");
      ss<<" "<<satellites.at(idSat).prnStr();
      UInt rowSum = 0;
      for(UInt idType=0; idType<types.size(); idType++)
      {
        UInt count = countTypes(idSat * types.size() + idType);
        ss<<" "<<count%"%8i"s;
        rowSum += count;
        colSum[idType] += count;
      }
      ss<<" | "<<rowSum%"%10i"s;
      logInfo<<ss.str()<<Log::endl;
    }
    ss.str("");
    ss<<" SUM";
    UInt allSum = 0;
    for(UInt idType=0; idType<types.size(); idType++)
    {
      ss<<" "<<colSum[idType]%"%8i"s;
      allSum += colSum[idType];
    }
    ss<<" | "<<allSum%"%10i"s;
    logInfo<<ss.str()<<Log::endl;

    UInt countTracks = 0;
    for(auto recv : gnss->receivers)
      if(recv->isMyRank())
        countTracks += recv->tracks.size();
    Parallel::reduceSum(countTracks, 0, comm);
    logInfo<<"  number of tracks: "<<countTracks<<Log::endl;

    // Processing steps
    // ----------------
    GnssProcessingStep::State state(gnss, comm);
    processingSteps->process(state);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
