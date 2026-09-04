/***********************************************/
/**
* @file fileSatelliteModel.h
*
* @brief Satellite model.
*
* @author Torsten Mayer-Guerr
* @author Sebastian Strasser
* @date 2015-05-25
*
*/
/***********************************************/

#ifndef __GROOPS_FILESATELLITEMODEL__
#define __GROOPS_FILESATELLITEMODEL__

// Latex documentation
#ifdef DOCSTRING_FILEFORMAT_SatelliteModel
static const char *docstringSatelliteModel = R"(
Properties of a satellite to model non-conservative forces (e.g. \configClass{miscAccelerations}{miscAccelerationsType}).
The file may contain surface properties, mass, drag coefficients, and antenna thrust values.

See \program{SatelliteModelCreate} and \program{SinexMetadata2SatelliteModel}.

\begin{verbatim}
<?xml version="1.0" encoding="UTF-8"?>
<groops type="satelliteModel" version="20190429">
   <satelliteCount>1</satelliteCount>
   <satellite>
       <satelliteName>GALILEO-2</satelliteName>
       <mass>7.00000000000000e+02</mass>
       <coefficientDrag>0.00000000000000e+00</coefficientDrag>
       <surfaceCount>15</surfaceCount>
       <surface>
           <type>0</type>
           <normal>
               <x>-1.00000000000000e+00</x>
               <y>0.00000000000000e+00</y>
               <z>0.00000000000000e+00</z>
           </normal>
           <area>4.40000000000000e-01</area>
           <reflexionVisible>0.00000000000000e+00</reflexionVisible>
           <diffusionVisible>7.00000000000000e-02</diffusionVisible>
           <absorptionVisible>9.30000000000000e-01</absorptionVisible>
           <reflexionInfrared>1.00000000000000e-01</reflexionInfrared>
           <diffusionInfrared>1.00000000000000e-01</diffusionInfrared>
           <absorptionInfrared>8.00000000000000e-01</absorptionInfrared>
           <hasThermalReemission>1</hasThermalReemission>
       </surface>
       ...
       <modulCount>2</modulCount>
       <modul>
           <type>1</type>
           <rotationAxis>
               <x>0.00000000000000e+00</x>
               <y>1.00000000000000e+00</y>
               <z>0.00000000000000e+00</z>
           </rotationAxis>
           <normal>
               <x>0.00000000000000e+00</x>
               <y>0.00000000000000e+00</y>
               <z>1.00000000000000e+00</z>
           </normal>
           <surface>
               <count>4</count>
               <cell>11</cell>
               <cell>12</cell>
               <cell>13</cell>
               <cell>14</cell>
           </surface>
       </modul>
       <modul>
           <type>2</type>
           <antennaThrust>
               <x>0.00000000000000e+00</x>
               <y>0.00000000000000e+00</y>
               <z>2.65000000000000e+02</z>
           </antennaThrust>
       </modul>
   </satellite>
</groops>
\end{verbatim}
)";
#endif

/***********************************************/

#include "inputOutput/fileName.h"
#include "inputOutput/fileArchive.h"
#include "inputOutput/logging.h"

/** @addtogroup filesGroup */
/// @{

/***** CONSTANTS ********************************/

const char *const FILE_SATELLITEMODEL_TYPE    = "satelliteModel";
constexpr UInt    FILE_SATELLITEMODEL_VERSION = std::max(UInt(20200123), FILE_BASE_VERSION);

/***** TYPES ***********************************/

class SatelliteModel;
class SatelliteModelModule;
typedef std::shared_ptr<SatelliteModel> SatelliteModelPtr;
typedef std::shared_ptr<SatelliteModelModule> SatelliteModelModulePtr;

/***** CLASS ***********************************/

/** @brief Satellite model. */
class SatelliteModel
{
public:
  /** @brief One surface of the satellite macro model. */
  class Surface
  {
  public:
    enum Type {PLATE=0, SPHERE=1, CYLINDER=2, GRACESHADOW=3};
    /// Type of the surface
    Type     type;
    /// Normal vector of the surface in the satellite reference frame (SRF)
    Vector3d normal;
    /// Area of the surface in [m^2]
    Double   area;
    Double   absorptionVisible, absorptionInfrared;
    Double   diffusionVisible,  diffusionInfrared;
    Double   reflexionVisible,  reflexionInfrared;
    /// Property about thermal reemission: 0, no thermal radiation; -1: spontaneous reemission of absorbed radiation [Ws/K/m^2]
    Double   specificHeatCapacity;

    Surface();
  };

  SatelliteModel();                                            //!< Constructor.
 ~SatelliteModel();                                            //!< Destructor.
  SatelliteModel &operator=(const SatelliteModel &x) = delete; //!< Disallow copying.
  SatelliteModel(const SatelliteModel &x) = delete;            //!< Disallow copy constructor.

  /// Name of the satellite.
  std::string                          satelliteName;
  /// Mass of satellite [kg].
  Double                               mass;
  /// Drag coefficient of satellite.
  Double                               coefficientDrag;
  /// Thrust power of satellite, e.g. of GNSS transmitting antenna in SRF [W].
  Vector3d                             thrustPower;
  /// Surfaces of satellite.
  std::vector<Surface>                 surfaces;
  /// Modules of the satellite
  std::vector<SatelliteModelModulePtr> modules;

  /**
   * @brief Computes the new state of the satellite.
   * E.g. moves parts of satellite surfaces, updates its mass.
   * @param time Time.
   * @param position in CRF [m].
   * @param velocity in CRF [m/s].
   * @param positionSun in CRF [m].
   * @param rotSat   Sat -> CRF
   * @param rotEarth CRF -> TRF */
  void changeState(const Time &time, const Vector3d &position, const Vector3d &velocity,
                   const Vector3d &positionSun, const Rotary3d &rotSat, const Rotary3d &rotEarth);
};

/***** FUNCTIONS *******************************/

template<> void save(OutArchive &ar, const SatelliteModelPtr &x);
template<> void load(InArchive  &ar, SatelliteModelPtr &x);

/** @brief Write one satellite model into a SatelliteModel file. */
void writeFileSatelliteModel(const FileName &fileName, const SatelliteModelPtr &x);

/** @brief Write multiple satellite models into a SatelliteModel file. */
void writeFileSatelliteModel(const FileName &fileName, const std::vector<SatelliteModelPtr> &x);

/** @brief Read one satellite model from a SatelliteModel file. */
void readFileSatelliteModel(const FileName &fileName, SatelliteModelPtr &x);

/** @brief Read multiple satellite models from a SatelliteModel file. */
void readFileSatelliteModel(const FileName &fileName, std::vector<SatelliteModelPtr> &x);

/// @}

/***** CLASS ***********************************/

/** @brief Base class for satellite model modules. */
class SatelliteModelModule
{
public:
  enum Type {SOLARPANEL = 1, ANTENNATHRUST = 2, MASSCHANGE = 3, SPECIFICHEATCAPACITY = 4};

  virtual ~SatelliteModelModule() {}
  /// @brief Returns the type of the module. */
  virtual Type type() const = 0;
  virtual void changeState(SatelliteModel &/*satellite*/, const Time &/*time*/, const Vector3d &/*position*/, const Vector3d &/*velocity*/,
                           const Vector3d &/*positionSun*/, const Rotary3d &/*rotSat*/, const Rotary3d &/*rotEarth*/) {}
  virtual void save(OutArchive &ar) const = 0;
  virtual void load(InArchive  &ar) = 0;

  /** @brief Create a satellite model module of given type (with new). */
  static SatelliteModelModulePtr create(Type type);
};

/***********************************************/

/** @brief Solar panel module of a satellite model. */
class SatelliteModelModuleSolarPanel : public SatelliteModelModule
{
public:
  /// Rotation axis of the solar panel in the satellite reference frame (SRF)
  Vector3d          rotationAxis;
  Vector3d          normal;
  /// Indexes of the solar panel surfaces in the satellite model surfaces vector, starting with 0.
  std::vector<UInt> indexSurface;

  Type type() const {return SOLARPANEL;}
  void changeState(SatelliteModel &satellite, const Time &time, const Vector3d &position, const Vector3d &velocity,
                   const Vector3d &positionSun, const Rotary3d &rotSat, const Rotary3d &rotEarth);
  void save(OutArchive &oa) const;
  void load(InArchive  &ia);
};

/***********************************************/

/** @brief Antenna thrust module of a satellite model. */
class SatelliteModelModuleAntennaThrust : public SatelliteModelModule
{
public:
  /// Flag if the thrust is applied or not.
  Bool     applied;
  /// Thrust vector in the satellite reference frame (SRF)
  Vector3d thrust;

  Type type() const {return ANTENNATHRUST;}
  void changeState(SatelliteModel &satellite, const Time &time, const Vector3d &position, const Vector3d &velocity,
                   const Vector3d &positionSun, const Rotary3d &rotSat, const Rotary3d &rotEarth);
  void save(OutArchive &oa) const;
  void load(InArchive  &ia);
};

/***********************************************/

/** @brief Mass change module of a satellite model. */
class SatelliteModelModuleMassChange : public SatelliteModelModule
{
public:
  std::vector<Time>   times;
  std::vector<Double> mass;

  Type type() const {return MASSCHANGE;}
  void changeState(SatelliteModel &satellite, const Time &time, const Vector3d &position, const Vector3d &velocity,
                   const Vector3d &positionSun, const Rotary3d &rotSat, const Rotary3d &rotEarth);
  void save(OutArchive &oa) const;
  void load(InArchive  &ia);
};

/***********************************************/

/** @brief Specific heat capacity module of a satellite model. */
class SatelliteModelModuleSetSpecificHeatCapacity : public SatelliteModelModule
{
public:
  Bool                applied;
  std::vector<Double> specificHeatCapacity;
  /// Indexes of the surfaces in the satellite model surfaces vector, starting with 0.
  std::vector<UInt>   indexSurface;

  Type type() const {return SPECIFICHEATCAPACITY;}
  void changeState(SatelliteModel &satellite, const Time &time, const Vector3d &position, const Vector3d &velocity,
                   const Vector3d &positionSun, const Rotary3d &rotSat, const Rotary3d &rotEarth);
  void save(OutArchive &oa) const;
  void load(InArchive  &ia);
};

/***********************************************/

#endif /* __GROOPS__ */
