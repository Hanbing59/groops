/***********************************************/
/**
* @file slrProcessingStep.h
*
* @brief Processing steps for SLR normal equations.
*
* @author Torsten Mayer-Guerr
* @date 2022-04-28
*
*/
/***********************************************/

#ifndef __GROOPS_SLRPROCESSINGSTEP__
#define __GROOPS_SLRPROCESSINGSTEP__

// Latex documentation
#ifdef DOCSTRING_SlrProcessingStep
static const char *docstringSlrProcessingStep = R"(
\section{SlrProcessingStep}\label{slrProcessingStepType}
Processing step in \program{SlrProcessing}.

Processing steps enable a dynamic definition of the consecutive steps performed during any kind of SLR processing.
The most common steps are \configClass{estimate}{slrProcessingStepType:estimate}, which performs an iterative least
squares adjustment, and \configClass{writeResults}{slrProcessingStepType:writeResults}, which writes all output files
defined in \program{SlrProcessing} and is usually the last step.
Some steps such as \configClass{selectParametrizations}{slrProcessingStepType:selectParametrizations}
and \configClass{selectStations}{slrProcessingStepType:selectStations} affect all subsequent steps.
In case these steps are used within a \configClass{group}{slrProcessingStepType:group} step,
they only affect the steps within this level.
)";
#endif

/***********************************************/

#include "base/import.h"
#include "parallel/matrixDistributed.h"
#include "slr/slr.h"

/**
* @defgroup slrProcessingStepGroup SlrProcessingStep
* @brief Processing steps in @ref SlrProcessing.
* @ingroup slrGroup
* The interface is given by @ref SlrProcessingStep. */
/// @{

/***** TYPES ***********************************/

class Slr;
class SlrProcessingStep;
class SlrProcessingStepBase;
typedef std::shared_ptr<SlrProcessingStep> SlrProcessingStepPtr;

/***** CLASS ***********************************/

/** @brief Provides a list of SLR processing steps.
* An Instance of this class can be created by the function @ref readConfig(). */
class SlrProcessingStep
{
  /// The list of configured SLR processing steps
  std::vector<SlrProcessingStepBase*> bases;

public:
  /** @brief Stores and manages the state of the SLR processing. */
  class State
  {
  public:
    SlrPtr                 slr;
    /// Information about the normal equations, e.g. parameter names and block structure.
    SlrNormalEquationInfo  normalEquationInfo;
    /// Whether the normal equation info has changed and the normal equations need to be rebuilt.
    Bool                   changedNormalEquationInfo;
    /// Normal equations matrix.
    MatrixDistributed      normals;
    /// Right hand side of the normal equations, stored at master after solution.
    std::vector<Matrix>    n;
    /// lPl value, stored at master after solution.
    Vector                 lPl;
    /// Number of observations, stored at master after solution.
    UInt                   obsCount;
    /// Sigma factor for each station
    Vector                 sigmaFactor;

    /** @brief Constructor. */
    State(SlrPtr slr);

    /** @brief Regularize unused parameters by setting their corresponding diagonal elements in the normal matrix to 1. 
     * @param[in] blockStart Index of the first block to regularize.
     * @param[in] blockCount Number of blocks to regularize starting from blockStart.
    */
    void   regularizeNotUsedParameters(UInt blockStart, UInt blockCount);

    /** @brief Accumulate the normal equations
     * @param constraintsOnly Whether consider only constranits
     */
    void   buildNormals(Bool constraintsOnly);

    /** @brief Do the solution estimation.
     * @param computeResiduals
     * @param computeWeights
     * @param adjustSigma0
     * @param huber
     * @param huberPower
     */
    Double estimateSolution(Bool computeResiduals,  Bool computeWeights, Bool adjustSigma0, Double huber, Double huberPower);

    /** @brief Do the statistics for residuals of a station-satellite pair.
     * @param[in] idStat
     * @param[in] idSat
     * @param[out] ePe
     * @param[out] redundancy
     * @param[out] obsCount
     * @param[out] outlierCount
     */
    void   residualsStatistics(UInt idStat, UInt idSat, Double &ePe, Double &redundancy, UInt &obsCount, UInt &outlierCount);
  };

  /** @brief Constructor from configuration node @a config. */
  SlrProcessingStep(Config &config, const std::string &name);

  /** @brief Destructor. */
 ~SlrProcessingStep();

  /** @brief Performs those configured processing steps one by one. */
  void process(State &state);

  /** @brief Creates a derived instance of this class. */
  static SlrProcessingStepPtr create(Config &config, const std::string &name) {return SlrProcessingStepPtr(new SlrProcessingStep(config, name));}
};

/***** FUNCTIONS *******************************/

/** @brief Creates an instance of the class SlrProcessingStep.
* Search for a node with @a name in the Config node.
* if @a name is not found the function returns FALSE and @a var is untouched.
* @param config The config node which includes the node with the options for this class
* @param name Tag name in the config.
* @param[out] var Created class.
* @param mustSet If is MUSTSET and @a name is not found, this function throws an exception instead of returning with FALSE.
* @param defaultValue Ignored at the moment.
* @param annotation Description of the function of this class.
* @relates SlrProcessingStep */
template<> Bool readConfig(Config &config, const std::string &name, SlrProcessingStepPtr &var, Config::Appearance mustSet, const std::string &defaultValue, const std::string &annotation);

/// @}

/***** CLASS ***********************************/

/** @brief Internal class, the base class for different SLR processing steps. */
class SlrProcessingStepBase
{
public:
  virtual ~SlrProcessingStepBase() {}

  /** @brief Execute this processing step. */
  virtual void process(SlrProcessingStep::State &state) = 0;
  /** @brief Whether this processing step expects the parameters to be initialized. */
  virtual Bool expectInitializedParameters() const {return TRUE;}
};

/***********************************************/

#endif /* __GROOPS___ */
