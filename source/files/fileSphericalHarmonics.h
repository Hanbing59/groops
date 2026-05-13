/***********************************************/
/**
* @file fileSphericalHarmonics.h
*
* @brief Read/Write a spherical harmonics coefficients file.
*
* @author Torsten Mayer-Guerr
* @date 2005-01-14
*
*/
/***********************************************/

#ifndef __GROOPS_FILESPHERICALHARMONICS__
#define __GROOPS_FILESPHERICALHARMONICS__

// Latex documentation
#ifdef DOCSTRING_FILEFORMAT_PotentialCoefficients
static const char *docstringPotentialCoefficients = R"(
The standard \verb|.gfc| format as defined by the ICGEM is used in the ASCII format.
Only the static part is used and temporal variations (e.g. trend) are ignored.
To write additional information and temporal variations, use \program{PotentialCoefficients2Icgem}.
)";
#endif

/***********************************************/

#include "base/exception.h"
#include "base/sphericalHarmonics.h"
#include "inputOutput/fileArchive.h"

/** @addtogroup filesGroup */
/// @{

/***** CONSTANTS ********************************/

const char *const FILE_POTENTIALCOEFFICIENTS_TYPE    = "potentialCoefficients";
constexpr UInt    FILE_POTENTIALCOEFFICIENTS_VERSION = std::max(UInt(20200123), FILE_BASE_VERSION);

/***** FUNCTIONS *******************************/

/** @brief Writes an SphericalHarmonics instance into a file. */
void writeFileSphericalHarmonics(const FileName &fileName, const SphericalHarmonics &x);

/** @brief Reads an SphericalHarmonics instance from a file. */
void readFileSphericalHarmonics(const FileName &fileName, SphericalHarmonics &x);

/// @}

/***********************************************/

#endif /* __GROOPS__ */
