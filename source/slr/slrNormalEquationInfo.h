/***********************************************/
/**
* @file slrNormalEquationInfo.h
*
* @brief SLR normal equations.
*
* @author Torsten Mayer-Guerr
* @date 2022-04-28
*
*/
/***********************************************/

#ifndef __GROOPS_SLRNORMALEQUATIONINFO__
#define __GROOPS_SLRNORMALEQUATIONINFO__

#include <regex>
#include "parallel/parallel.h"
#include "parallel/matrixDistributed.h"
#include "base/parameterName.h"

/** @addtogroup slrGroup */
/// @{

/***** CLASS ***********************************/

/** @brief SLR parameter index. */
class SlrParameterIndex
{
  UInt index;

public:
  SlrParameterIndex(UInt idx=NULLINDEX) : index(idx) {}
  explicit operator bool() const {return index != NULLINDEX;}

  friend class SlrNormalEquationInfo;
  friend class SlrDesignMatrix;
};

/***** CLASS ***********************************/

/** @brief Information about the SLR normal equation. */
class SlrNormalEquationInfo
{
public:
  /// wildcard matching, enable/disable
  std::vector<std::pair<std::regex, Bool>> enableParametrizations;
  std::vector<Byte>                        estimateStation;
  std::vector<Byte>                        estimateSatellite;
  UInt                                     defaultBlockSize;
  std::vector<SlrParameterIndex>           indexGravity;

  /** @brief Constructor. */
  SlrNormalEquationInfo(UInt countStations, UInt countSatellites);

  /** @brief Initializes the status of parameters. */
  void initNewParameterNames();

  /** @brief Adds station-specific parameters. */
  SlrParameterIndex parameterNamesStation  (UInt idStat, const std::vector<ParameterName> &parameterNames) {return addParameters(idStat,    NULLINDEX, parameterNames);}
  
  /** @brief Adds satellite-specific parameters. */
  SlrParameterIndex parameterNamesSatellite(UInt idSat,  const std::vector<ParameterName> &parameterNames) {return addParameters(NULLINDEX, idSat,     parameterNames);}
  
  /** @brief Adds other parameters. */
  SlrParameterIndex parameterNamesOther    (const std::vector<ParameterName> &parameterNames)              {return addParameters(NULLINDEX, NULLINDEX, parameterNames);}
  
  /** @brief Calculate indexes of normal matrix blocks. */
  void calculateIndex();

  UInt block(const SlrParameterIndex &index) const {return block_.at(index.index);}
  UInt index(const SlrParameterIndex &index) const {return index_.at(index.index);}
  UInt count(const SlrParameterIndex &index) const {return count_.at(index.index);}

  /** @brief Gets the list of parameters. */
  const std::vector<ParameterName> &parameterNames() const {return parameterNames_;}

  /** @brief Gets the list of start indexes of normal matrix blocks and the total parameter count as the last element. */
  const std::vector<UInt>          &blockIndices()   const {return blockIndices_;}

  /** @brief Gets the number of parameters. 
   * Number of rows/columns (dimension) of distributed matrix
  */
  UInt parameterCount()   const {return blockIndices_.back();}

  /** @brief Gets the start index of block @a i. 
   * Note that the last element of blockIndices() is the total parameter count.
  */
  UInt blockIndex(UInt i) const {return blockIndices_.at(i);}

  /** @brief Gets the size of block @a i. */
  UInt blockSize(UInt i)  const {return blockIndices_.at(i+1)-blockIndices_.at(i);}
  
  /** @brief Gets the number of blocks.
   * Note that the block count is blockIndices().size()-1, as the last element of blockIndices() is the total parameter count.
  */
  UInt blockCount()       const {return blockIndices_.size()-1;}

private:
  /// idStat, idSat, idx, name
  std::vector<std::tuple<UInt, UInt, UInt, std::vector<ParameterName>>> parameters;
  std::vector<UInt>          block_, index_, count_;
  /// Parameters list
  std::vector<ParameterName> parameterNames_;
  /// List of start indexes of normal matrix blocks and the total parameter count as the last element
  std::vector<UInt>          blockIndices_;
  std::vector<UInt>          blockCountEpoch_;

  SlrParameterIndex addParameters(UInt idStat, UInt idSat, const std::vector<ParameterName> &parameterNames);
};

/// @}

/***********************************************/

#endif
