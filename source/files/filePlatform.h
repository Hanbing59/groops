/***********************************************/
/**
* @file filePlatform.h
*
* @brief Platform equipped with instruments.
*
* @author Torsten Mayer-Guerr
* @date 2022-11-07
*
*/
/***********************************************/

#ifndef __GROOPS_FILEPLATFORM__
#define __GROOPS_FILEPLATFORM__

// Latex documentation
#ifdef DOCSTRING_FILEFORMAT_Platform
static const char *docstringPlatform = R"(
Defines a platform with a local coordinate frame equipped with instruments.
The platform might be a reference station, a low Earth satellite,
or a transmitting GNSS satellite and is referenced by a marker name and number.
The reference point (marker or center of mass (CoM)) can change in time
relative to the local frame.

Each equipped instrument is described at least by the following information
\begin{itemize}
\item name
\item serial number
\item coordinates in the local frame
\item a time interval in which the instrument was active
\item the orientation for antennas and reflectors.
\end{itemize}

For GNSS satellites the platform defines the PRN. The different assigned SVNs
are defined by the transmitting antennas.

Platforms for GNSS stations can be created from station log files with
\program{GnssStationLog2Platform}. Platforms for GNSS satellites
can be created from an ANTEX file with \program{GnssAntex2AntennaDefinition}.

See also \program{PlatformCreate}.

\fig{!hb}{0.8}{fileFormatPlatform}{fig:fileFormatPlatform}{Platform for stations, LEOs, and GNSS satellites.}
)";
#endif

/***********************************************/

#include "base/import.h"
#include "inputOutput/fileName.h"
#include "inputOutput/fileArchive.h"
#include "files/fileGnssAntennaDefinition.h"
#include "files/fileGnssReceiverDefinition.h"

/** @addtogroup filesGroup */
/// @{

/***** CONSTANTS ********************************/

const char *const FILE_PLATFORM_TYPE    = "platform";
constexpr UInt    FILE_PLATFORM_VERSION = std::max(UInt(20200123), FILE_BASE_VERSION);

/***** TYPES ***********************************/

class Platform;
class PlatformEquipment;
typedef std::shared_ptr<PlatformEquipment> PlatformEquipmentPtr;

class GnssAntennaDefinition;
class GnssReceiverDefinition;
typedef std::shared_ptr<GnssAntennaDefinition>  GnssAntennaDefinitionPtr;
typedef std::shared_ptr<GnssReceiverDefinition> GnssReceiverDefinitionPtr;

/***** CLASS ***********************************/

/** @brief Class for a platform. */
class Platform
{
public:
  /** @brief Class for a reference point of a platform, e.g.
   * Center of mass (CoM) for satellites
  */
  class ReferencePoint
  {
  public:
    std::string comment;
    Vector3d    pointStart, pointEnd;
    Time        timeStart, timeEnd;
  };
  // Platform name: PRN for transmitters, marker name for stations, satellite name for LEOs
  std::string                       name;
  std::string                       markerName, markerNumber;
  std::string                       comment;
  Vector3d                          approxPosition;
  /// time sorted list of reference points
  std::vector<ReferencePoint>       referencePoints;
  /// list of equipments on this platform
  std::vector<PlatformEquipmentPtr> equipments;

  /** @brief Gets the reference point positions at a given epoch. */
  Vector3d referencePoint(const Time &time) const;

  /** @brief Gets a specified type of equipment at a given epoch. */
  template<typename T> std::shared_ptr<T> findEquipment(const Time &time) const;

  /** @brief Sets the antenna definition for the GNSS antenna of the platform. */
  void fillGnssAntennaDefinition (const std::vector<GnssAntennaDefinitionPtr> &antennaList);

  /** @brief Sets the antenna accuracy for the GNSS antenna of the platform. */
  void fillGnssAccuracyDefinition(const std::vector<GnssAntennaDefinitionPtr> &antennaList);

  /** @brief Sets the receiver definition for the GNSS receiver of the platform. */
  void fillGnssReceiverDefinition(const std::vector<GnssReceiverDefinitionPtr> &receiverList);
};

/***** CLASS ***********************************/

/** @brief Base class for an equipment on a platform. */
class PlatformEquipment
{
public:
  /// List of all equipment types
  enum Type : Int {UNDEFINED           = 0,
                   OTHER               = 1,
                   GNSSANTENNA         = 2,
                   GNSSRECEIVER        = 3,
                   SLRSTATION          = 4,
                   LASERRETROREFLECTOR = 5,
                   SATELLITEIDENTIFIER = 6};

  /// Type of this equipment
  static constexpr Type TYPE = OTHER;
  std::string comment;
  std::string name, serial;
  Time        timeStart, timeEnd;
  /// positions of this instrument in North, East, Up or vehicle system
  Vector3d    position;

  virtual ~PlatformEquipment() {}

  /** @brief Creates an equipment instance of a given type (as shared_ptr). */
  static PlatformEquipmentPtr create(Type type);

  /** @brief Gets the type of this equipment, e.g. GNSSANTENNA, SLRREFLECTOR. */
  virtual Type getType() const {return TYPE;}

  /** @brief Returns the string ID of the equipment consisting of its name and serial number. */
  virtual std::string str() const {return name+"|"+serial;};

  virtual void save(OutArchive &oa) const;
  virtual void load(InArchive  &ia);
};

/***** CLASS ***********************************/

/** @brief Class for a GNSS antenna equipment. */
class PlatformGnssAntenna : public PlatformEquipment
{
public:
  static constexpr Type TYPE = GNSSANTENNA;
  std::string              radome;
  /// rotation from the north, east, up or vehicle system to the antenna system
  Transform3d              local2antennaFrame;
  /// Antenna definition for phase center patterns
  GnssAntennaDefinitionPtr antennaDef;
  /// Antenna definition for observation accuracy patterns
  GnssAntennaDefinitionPtr accuracyDef;

  Type getType() const override {return TYPE;}
  std::string str() const override {return GnssAntennaDefinition::str(name, serial, radome);}
  void save(OutArchive &oa) const override;
  void load(InArchive  &ia) override;
};

/***** CLASS ***********************************/

/** @brief Class for a GNSS receiver equipment. */
class PlatformGnssReceiver : public PlatformEquipment
{
public:
  static constexpr Type TYPE = GNSSRECEIVER;
  // software version
  std::string               version;
  GnssReceiverDefinitionPtr receiverDef;

  Type getType() const override {return TYPE;}
  std::string str() const override {return GnssReceiverDefinition::str(name, serial, version);}
  void save(OutArchive &oa) const override;
  void load(InArchive  &ia) override;
};

/***** CLASS ***********************************/

/** @brief Class for a SLR station equipment. */
class PlatformSlrStation : public PlatformEquipment
{
public:
  static constexpr Type TYPE = SLRSTATION;
  Type getType() const override {return TYPE;}
  std::string str() const override {return name;}
  void save(OutArchive &oa) const override;
  void load(InArchive  &ia) override;
};

/***** CLASS ***********************************/

/** @brief Class for a SLR LRR equipment. */
class PlatformLaserRetroReflector : public PlatformEquipment
{
public:
  static constexpr Type TYPE = LASERRETROREFLECTOR;
  /// rotation from the satellite system to the reflector system
  Transform3d platform2reflectorFrame;
  Angle       dZenit;
  /// range variations (azimut(0..360) x zenit(0..dZenit*rows))
  Matrix      range;

  Type getType() const override {return TYPE;}
  std::string str() const override {return name;}
  void save(OutArchive &oa) const override;
  void load(InArchive  &ia) override;
};

/***** CLASS ***********************************/

/** @brief Class for a satellite ID equipment. */
class PlatformSatelliteIdentifier : public PlatformEquipment
{
public:
  static constexpr Type TYPE = SATELLITEIDENTIFIER;
  std::string cospar, norad, sic, sp3;

  Type getType() const override {return TYPE;}
  std::string str() const override {return name;}
  void save(OutArchive &oa) const override;
  void load(InArchive  &ia) override;
};

/***** FUNCTIONS *******************************/

template<> void save(OutArchive &ar, const Platform &x);
template<> void load(InArchive  &ar, Platform &x);

/** @brief Write a platform into a Platform file. */
void writeFilePlatform(const FileName &fileName, const Platform &x);

/** @brief Read a platform from a Platform file. */
void readFilePlatform(const FileName &fileName, Platform &x);

/***********************************************/
/***** INLINES   *******************************/
/***********************************************/

template<typename T> inline std::shared_ptr<T> Platform::findEquipment(const Time &time) const
{
  try
  {
    auto iter = std::find_if(equipments.begin(), equipments.end(), [&](const auto &x)
                            {return (x->getType() == T::TYPE) && (x->timeStart <= time) && (time < x->timeEnd);});
    if(iter == equipments.end())
      return std::shared_ptr<T>(nullptr);
    return std::dynamic_pointer_cast<T>(*iter);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

/// @}

#endif /* __GROOPS___ */
