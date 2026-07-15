/***********************************************/
/**
* @file orbitPropagator.h
*
* @brief Propagate orbit using different algorithms.
*
* @author Matthias Ellmer
* @date 2017-01-19
*
*/
/***********************************************/

#ifndef __GROOPS_ORBITPROPAGATOR__
#define __GROOPS_ORBITPROPAGATOR__

// Latex documentation
#ifdef DOCSTRING_OrbitPropagator
static const char *docstringOrbitPropagator = R"(
\section{OrbitPropagator}\label{orbitPropagatorType}
Implements the propagation of a satellite orbit under
the influence of \configClass{forces}{forcesType} as
used in \program{SimulateOrbit}
(dynamic orbits from numerical orbit integration).
)";
#endif

/***********************************************/

#include "base/import.h"
#include "config/config.h"
#include "files/fileInstrument.h"
#include "files/fileSatelliteModel.h"
#include "classes/forces/forces.h"

/**
* @defgroup orbitPropagatorGroup OrbitPropagator
* @brief Propagate a satellite orbit.
* @ingroup classesGroup
* The interface is given by @ref OrbitPropagator. */
/// @{

/***** TYPES ***********************************/

class OrbitPropagator;
typedef std::shared_ptr<OrbitPropagator> OrbitPropagatorPtr;

/***** CLASS ***********************************/

/** @brief Base class for satellite orbit propagators.
 * An instance of this class can be created by @ref readConfig. */
class OrbitPropagator
{
public:
  /** @brief Destructor. */
  virtual ~OrbitPropagator() {}

  /**
   * @brief Propagates an orbit arc through specified force fields.
   * @param startEpoch Initial epoch and state vector of the integration
   * @param sampling Time difference between epochs
   * @param posCount Total number of integration epochs (including @a startEpoch)
   * @param forces Force models to be considered
   * @param satellite Satellite macro model to consider
   * @param earthRotation Rotations between celestial and terrestrial reference frames
   * @param ephemerides Planetary ephemerides.
   * @param timing Whether to show a timer.
   * @returns OrbitArc Integrated orbit arc of size @a posCount, with @a startEpoch at index 0  */
  virtual OrbitArc integrateArc(const OrbitEpoch &startEpoch, const Time &sampling, UInt posCount, ForcesPtr forces, SatelliteModelPtr satellite,
                                EarthRotationPtr earthRotation, EphemeridesPtr ephemerides, Bool timing=TRUE) const = 0;

  /**
   * @brief Generates the rotation from the orbit local system to CRF.
   *
   * The orbit local system is defined as follows:
   * - x-axis, along the velocity vector (tangential)
   * - y-axis, along the orbit normal vector (normal)
   * - z-axis, orthogonal to x- and y-axis according to the right-hand rule,
   *           pointing towards the center of the orbit (radial)
   * @param time Time of evaluation
   * @param position Position of the satellite
   * @param velocity Velocity of the satellite
   * @param satellite Satellite macro model
   * @returns Rotary3d Rotation from the orbit local system to CRF. */
  virtual Rotary3d orientation(const Time &time, const Vector3d &position, const Vector3d &velocity, SatelliteModelPtr satellite) const;

  /** @brief Flips the arc, i.e. re-arranges the orbit in the order of reversed epochs. */
  static Arc flip(const Arc &arc);

  /** @brief Creates a derived instance of this class. */
  static OrbitPropagatorPtr create(Config &config, const std::string &name);

protected:
  /**
   * @brief Computes the accelerations of a satellite at a given epoch.
   * @param epoch Epoch with time, position and velocity
   * @param forces Force models to be considered
   * @param satellite Satellite macro model to consider
   * @param earthRotation Rotations between celestial and terrestrial reference frames
   * @param ephemerides Planetary ephemerides.
   * @returns Acceleration of the satellite in CRF [m/s^2] */
  Vector3d acceleration(const OrbitEpoch &epoch, ForcesPtr forces, SatelliteModelPtr satellite, EarthRotationPtr earthRotation, EphemeridesPtr ephemerides) const;
};

/***** FUNCTIONS *******************************/

/** @brief Creates an instance of the class OrbitPropagator.
* Search for a node with @a name in the Config node.
* if @a name is not found the function returns FALSE and @a orbitPropagator is untouched.
* @param config The config node which includes the node with the options for this class
* @param name Tag name in the config.
* @param[out] orbitPropagator Created class.
* @param mustSet If is MUSTSET and @a name is not found, this function throws an exception instead of returning with FALSE.
* @param defaultValue Ignored at the moment.
* @param annotation Description of the function of this class.
* @relates OrbitPropagator */
template<> Bool readConfig(Config &config, const std::string &name, OrbitPropagatorPtr &orbitPropagator, Config::Appearance mustSet, const std::string &defaultValue, const std::string &annotation);

/// @}

#endif
