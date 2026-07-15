/***********************************************/
/**
* @file gnssParametrizationTecBiases.h
*
* @brief TEC biases.
* @see GnssParametrization
*
* @author Torsten Mayer-Guerr
* @date 2021-01-23
*
*/
/***********************************************/

#ifndef __GROOPS_GNSSPARAMETRIZATIONTECBIASES__
#define __GROOPS_GNSSPARAMETRIZATIONTECBIASES__

// Latex documentation
#ifdef DOCSTRING_GnssParametrization
static const char *docstringGnssParametrizationTecBiases = R"(
\subsection{TecBiases}\label{gnssParametrizationType:tecBiases}
Each code observation (e.g \verb|C1C| or \verb|C2W|) contains a bias at transmitter/receiver level
\begin{equation}
  [C\nu a]_r^s(t) = \dots + \text{bias}[C\nu a]^s + \text{bias}[C\nu a]_r + \ldots
\end{equation}
This parametrization represents the linear combination of signal biases
which completely depend on the STEC parameters. Ignoring these bias combinations would result
in a biased STEC estimation (all other parameters are nearly unaffected).
To determine this part of the signal biases
the \configClass{parametrization:ionosphereSTEC}{gnssParametrizationType:ionosphereSTEC} should be constrained.
Furthermore, additional information about the ionosphere is required from
\configClass{parametrization:ionosphereVTEC}{gnssParametrizationType:ionosphereVTEC} or
\configClass{parametrization:ionosphereMap}{gnssParametrizationType:ionosphereMap}.

Rank deficiencies due to the signal bias parameters may occur if biases of
transmitters and receivers are estimated together.
The minimum norm nullspace is formulated as zero constraint equations and added with
a standard deviation of \config{sigmaZeroMeanConstraint}.

The accumulated estimated result can be written to files in
\configClass{parametrization:signalBiases}{gnssParametrizationType:signalBiases}.

The \file{parameter names}{parameterName} are \verb|<station or prn>:tecBias0<index><combi of gnssTypes>::|.
)";
#endif

/***********************************************/

#include "base/import.h"
#include "config/config.h"
#include "classes/platformSelector/platformSelector.h"
#include "gnss/gnss.h"
#include "gnss/gnssParametrization/gnssParametrization.h"

/***** CLASS ***********************************/

/** @brief TEC biases.
* @ingroup gnssParametrizationGroup
* @see GnssParametrization */
class GnssParametrizationTecBiases : public GnssParametrizationBase
{
  /** @brief Class for a TEC bias parameter. */
  class Parameter
  {
  public:
    /// Transmitter of this bias parameter
    GnssTransmitterPtr trans;
    /// Receiver of this bias parameter
    GnssReceiverPtr    recv;
    /// Index of this bias parameter in the normal equation
    GnssParameterIndex index;
    /// the linear combinations of signal biases that cannot be estimated from the observations
    Matrix             Bias;
  };

  Gnss                    *gnss;
  std::string              name, nameConstraint;
  PlatformSelectorPtr      selectTransmitters, selectReceivers;
  FileName                 fileNameTransmitter, fileNameReceiver;
  /// Whether assumes a linear dependence on the frequency channel for GLONASS
  Bool                     isLinearBias;
  /// Whether to apply the zero mean constraint on the bias parameters
  Bool                     applyConstraint;
  Double                   sigmaZeroMean;
  /// List of parameters for transmitter biases
  std::vector<Parameter*>  paraTrans;
  /// List of parameters for receiver biases
  std::vector<Parameter*>  paraRecv;
  /// indices in zeroMean matrix
  std::vector<UInt>        idxBiasTrans, idxBiasRecv;
  /// zero mean observation equations
  Matrix                   zeroMeanDesign;

public:
  GnssParametrizationTecBiases(Config &config);
 ~GnssParametrizationTecBiases();

  /** @brief Initializes the bias parameters lists. */
  void   init(Gnss *gnss, Parallel::CommunicatorPtr comm) override;

  /** @brief Initializes the TEC bias parameters in the normal equation. */
  void   initParameter(GnssNormalEquationInfo &normalEquationInfo) override;

  /** @brief Gets the a priori values of the TEC bias parameters. */
  void   aprioriParameter(const GnssNormalEquationInfo &normalEquationInfo, MatrixSliceRef x0) const override;

  /** @brief Forms the design matrix for the TEC bias parameters. */
  void   designMatrix(const GnssNormalEquationInfo &normalEquationInfo, const GnssObservationEquation &eqn, GnssDesignMatrix &A) const override;

  /** @brief Applies the constraints on the TEC bias parameters. */
  void   constraints(const GnssNormalEquationInfo &normalEquationInfo, MatrixDistributed &normals, std::vector<Matrix> &n, Double &lPl, UInt &obsCount) const override;

  /** @brief Updates the TEC bias parameters */
  Double updateParameter(const GnssNormalEquationInfo &normalEquationInfo, const_MatrixSliceRef x, const_MatrixSliceRef Wz) override;
};

/***********************************************/

#endif
