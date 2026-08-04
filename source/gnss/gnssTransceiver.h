/***********************************************/
/**
* @file gnssTransceiver.h
*
* @brief GNSS receiver or transmitter.
*
* @author Torsten Mayer-Guerr
* @date 2021-01-30
*
*/
/***********************************************/

#ifndef __GROOPS_GNSSTRANSCEIVER__
#define __GROOPS_GNSSTRANSCEIVER__

#include "base/gnssType.h"
#include "files/fileGnssSignalBias.h"
#include "files/filePlatform.h"

/** @addtogroup gnssGroup */
/// @{

/***** TYPES ***********************************/

class GnssTransceiver;
typedef std::shared_ptr<GnssTransceiver> GnssTransceiverPtr;

/***** CLASS ***********************************/

/** @brief Abstract class for GNSS receiver or transmitter. */
class GnssTransceiver
{
  /// List of usable epochs of the platform
  Vector      useableEpochs;
  /// Number of usable epochs of the platform
  UInt        countUseableEpochs;
  /// Action to be taken if no antenna pattern is found for a given frequency
  GnssAntennaDefinition::NoPatternFoundAction noPatternFoundAction;

public:
  /// Index of the GNSS transceiver, set by Gnss::init()
  UInt           id_;
  /// Platform for the GNSS transceiver
  Platform       platform;
  /// Signal biases for the GNSS transceiver
  GnssSignalBias signalBias;

public:
  /** @brief Constructor. */
  GnssTransceiver(const Platform &platform, GnssAntennaDefinition::NoPatternFoundAction noPatternFoundAction, const Vector &useableEpochs);

  /** @brief Destructor. */
  virtual ~GnssTransceiver() {}

  /** @brief Returns the name of the transceiver (name of the platform). */
  std::string name() const {return platform.name;}

  /**
   * @brief Checks if the platform is usable at one given epoch or any epoch.
   * @param idEpoch Index of the epoch to check. If NULLINDEX, checks if there are any usable epochs.
   * @return True if the platform is usable at the specified epoch or any epoch if idEpoch is NULLINDEX.
   */
  Bool useable(UInt idEpoch=NULLINDEX) const {return countUseableEpochs && ((idEpoch == NULLINDEX) || useableEpochs(idEpoch));}

  /** @brief Disables a given epoch of the transceiver for a specific reason.
   * If no useable epochs remain after this disabling, the transceiver will be disabled.
  */
  virtual void disable(UInt idEpoch, const std::string &reason);

  /** @brief Disables the transceiver for a specific reason. */
  virtual void disable(const std::string &reason);

  /** @brief Returns the allowed signal types at a given time,
   * which will be empty if no GNSS receiver was found or
   * no signal / receiver definition was provided for this platform. */
  std::vector<GnssType> definedTypes(const Time &time) const;

  /**
   * @brief Returns the signal direction-dependent phase center corrections for GNSS range measurements at a given time.
   * Signal biases are included in the correction.
   * @param time Time at which to calculate the correction.
   * @param azimuth Azimuth angle in the antenna frame.
   * @param elevation Elevation angle in the antenna frame.
   * @param type Vector of GNSS signal types.
   * @return Vector of correction values.
   */
  Vector antennaVariations(const Time &time, Angle azimut, Angle elevation, const std::vector<GnssType> &type) const;

  /**
   * @brief Returns the signal direction-dependent standard deviation for GNSS range measurements at a given time.
   * @param time Time at which to calculate the accuracy.
   * @param azimut Azimuth angle in the (left-handed) antenna frame.
   * @param elevation Elevation angle in the (left-handed) antenna frame.
   * @param type Vector of GNSS signal types.
   * @return Vector of standard deviation values.
   */
  Vector accuracy(const Time &time, Angle azimut, Angle elevation, const std::vector<GnssType> &type) const;

  void save(OutArchive &oa) const;
  void load(InArchive  &ia);
};

/***********************************************/

inline GnssTransceiver::GnssTransceiver(const Platform &platform, GnssAntennaDefinition::NoPatternFoundAction noPatternFoundAction, const Vector &useableEpochs)
  : useableEpochs(useableEpochs), countUseableEpochs(sum(useableEpochs)), noPatternFoundAction(noPatternFoundAction), platform(platform) {}

/***********************************************/

inline void GnssTransceiver::disable(UInt idEpoch, const std::string &reason)
{
  try
  {
    if(useableEpochs(idEpoch))
      countUseableEpochs--;
    useableEpochs(idEpoch) = FALSE;
    if(countUseableEpochs == 0)
      disable(reason);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

inline void GnssTransceiver::disable(const std::string &/*reason*/)
{
  try
  {
    countUseableEpochs = 0;
    useableEpochs.setNull();
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

inline std::vector<GnssType> GnssTransceiver::definedTypes(const Time &time) const
{
  try
  {
    auto receiver = platform.findEquipment<PlatformGnssReceiver>(time);
    if(!receiver || !receiver->receiverDef)
      return std::vector<GnssType>();
    return receiver->receiverDef->types;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

inline Vector GnssTransceiver::antennaVariations(const Time &time, Angle azimut, Angle elevation, const std::vector<GnssType> &types) const
{
  try
  {
    Vector corr(types.size());
    // Have those signal biases been initialized at this moment?
    corr += signalBias.compute(types);

    auto antenna = platform.findEquipment<PlatformGnssAntenna>(time);
    if(!antenna)
      throw(Exception(platform.markerName+"."+platform.markerNumber+": no antenna definition found at "+time.dateTimeStr()));
    if(!antenna->antennaDef)
      throw(Exception("no antenna definition for "+antenna->str()));
    corr += antenna->antennaDef->antennaVariations(azimut, elevation, types, noPatternFoundAction);

    return corr;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

inline Vector GnssTransceiver::accuracy(const Time &time, Angle azimut, Angle elevation, const std::vector<GnssType> &types) const
{
  try
  {
    auto antenna = platform.findEquipment<PlatformGnssAntenna>(time);
    if(!antenna)
      throw(Exception(platform.markerName+"."+platform.markerNumber+": no antenna accuracy found at "+time.dateTimeStr()));
    if(!antenna->accuracyDef)
      throw(Exception("no accuracy definition for "+antenna->str()));
    return antenna->accuracyDef->antennaVariations(azimut, elevation, types, noPatternFoundAction);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

inline void GnssTransceiver::save(OutArchive &oa) const
{
  oa<<nameValue("useableEpochs",      useableEpochs);
  oa<<nameValue("countUseableEpochs", countUseableEpochs);
  oa<<nameValue("signalBias",         signalBias);
}

/***********************************************/

inline void GnssTransceiver::load(InArchive  &ia)
{
  ia>>nameValue("useableEpochs",      useableEpochs);
  ia>>nameValue("countUseableEpochs", countUseableEpochs);
  ia>>nameValue("signalBias",         signalBias);
}

/***********************************************/

/// @}

#endif /* __GROOPS___ */
