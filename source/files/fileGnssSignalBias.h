/***********************************************/
/**
* @file fileGnssSignalBias.h
*
* @brief File for GNSS code/phase biases.
*
* @author Torsten Mayer-Guerr
* @date 2013-08-11
*
*/
/***********************************************/

#ifndef __GROOPS_GNSSSIGNALBIAS__
#define __GROOPS_GNSSSIGNALBIAS__

// Latex documentation
#ifdef DOCSTRING_FILEFORMAT_GnssSignalBias
static const char *docstringGnssSignalBias = R"(
Signal biases of GNSS transmitters or receivers for different \configClass{gnssType}{gnssType}.

\begin{verbatim}
groops gnssSignalBias version=20200123
          5 # number of signals
# type   bias [m]
# ===============================
 C1CG06 -1.752461109688110974e-01
 C1WG06  4.005884595055994590e-02
 C2WG06  6.597469378913034532e-02
 L1*G06 -2.736169875580296909e-02
 L2*G06  3.422596762686257871e-02
 \end{verbatim}

See also \program{GnssProcessing}, \program{GnssSimulateReceiver}, \program{GnssSignalBias2Matrix}, \program{GnssSignalBias2SinexBias}.
)";
#endif

/***********************************************/

#include "base/gnssType.h"
#include "inputOutput/fileName.h"
#include "inputOutput/fileArchive.h"

/** @addtogroup filesGroup */
/// @{

/***** CONSTANTS ********************************/

const char *const FILE_GNSSSIGNALBIAS_TYPE    = "gnssSignalBias";
constexpr UInt    FILE_GNSSSIGNALBIAS_VERSION = std::max(UInt(20200123), FILE_BASE_VERSION);

/***** TYPES ***********************************/

class GnssSignalBias;
typedef std::shared_ptr<GnssSignalBias> GnssSignalBiasPtr;

/***** CLASS ***********************************/

/** @brief GNSS code/phase biases for a transmitter or receiver. */
class GnssSignalBias
{
  public:
  /// GNSS signal types
  std::vector<GnssType> types;
  /// GNSS signal biases in meters
  std::vector<Double>   biases;

  /** @brief Gets the signal biases for a list of specified GNSS signal types. */
  Vector compute(const std::vector<GnssType> &types) const;
};

/***** FUNCTIONS *******************************/

template<> void save(OutArchive &ar, const GnssSignalBias &x);
template<> void load(InArchive  &ar, GnssSignalBias &x);

/** @brief Writes the GNSS signal biases of a transmitter or receiver into a GnssSignalBias file. */
void writeFileGnssSignalBias(const FileName &fileName, const GnssSignalBias &x);

/** @brief Reads the GNSS signal biases of a transmitter or receiver from a GnssSignalBias file. */
void readFileGnssSignalBias(const FileName &fileName, GnssSignalBias &x);

/// @}

/***********************************************/

#endif /* __GROOPS___ */
