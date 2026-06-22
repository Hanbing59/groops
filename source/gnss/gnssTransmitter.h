/***********************************************/
/**
* @file gnssTransmitter.h
*
* @brief GNSS transmitter.
*
* @author Torsten Mayer-Guerr
* @author Sebastian Strasser
* @date 2013-06-28
*
*/
/***********************************************/

#ifndef __GROOPS_GNSSTRANSMITTER__
#define __GROOPS_GNSSTRANSMITTER__

#include "base/polynomial.h"
#include "base/gnssType.h"
#include "gnss/gnssTransceiver.h"

/** @addtogroup gnssGroup */
/// @{

/***** TYPES ***********************************/

class GnssTransmitter;
typedef std::shared_ptr<GnssTransmitter> GnssTransmitterPtr;

/***** CLASS ***********************************/

/** @brief Abstract class for GNSS transmitters. */
class GnssTransmitter : public GnssTransceiver
{
  /// PRN
  GnssType                 type;
  /// The interpolation polynomial for positions and velocities
  Polynomial               polynomial;
  /// Clock errors
  std::vector<Double>      clk;
  std::vector<Double>      scale;
  /// Offsets between CoM and Antenna Reference Point (ARP) in SRF
  std::vector<Vector3d>    offset;
  /// Rotation from the Celestial Reference Frame (CRF) to the Satellite Reference Frame (SRF)
  std::vector<Transform3d> crf2srf;
  /// Rotation from the Satellite Reference Frame (SRF) to the left-handed Antenna Reference Frame (ARF)
  std::vector<Transform3d> srf2arf;

public:
  std::vector<Time> timesPosVel;
  /// CoM in CRF (epoch times (x,y,z))
  Matrix            pos;
  Matrix            vel;

  /** @brief Constructor. */
  GnssTransmitter(GnssType prn, const Platform &platform,
                  GnssAntennaDefinition::NoPatternFoundAction noPatternFoundAction,
                  const Vector &useableEpochs, const std::vector<Double> &clock, const std::vector<Double> &scale, const std::vector<Vector3d> &offset,
                  const std::vector<Transform3d> &crf2srf, const std::vector<Transform3d> &srf2arf,
                  const std::vector<Time> &timesPosVel, const_MatrixSliceRef position, const_MatrixSliceRef velocity, UInt interpolationDegree)
  : GnssTransceiver(platform, noPatternFoundAction, useableEpochs),
    type(prn), polynomial(timesPosVel, interpolationDegree, TRUE/*throwException*/, FALSE/*leastSquares*/, -(interpolationDegree+1.1), -1.1, 1e-7),
    clk(clock), scale(scale), offset(offset), crf2srf(crf2srf), srf2arf(srf2arf), timesPosVel(timesPosVel), pos(position), vel(velocity) {}

  /** @brief Destructor. */
  virtual ~GnssTransmitter() {}

  /** @brief Gets ID of the transmitter. */
  UInt idTrans() const {return id_;}

  /** @brief Gets PRN of the transmitter. */
  GnssType PRN() const {return type;}

  /** @brief Gets the clock error at a given epoch.
  * error = clock time - system time [s] */
  Double clockError(UInt idEpoch) const {return clk.at(idEpoch);}

  /** @brief Gets the scale due to frequency offset/clock drift at a given epoch.
  * observed = scaleFactor * true_range */
  Double scaleFactor(UInt idEpoch) const {return scale.at(idEpoch);}

  /** @brief Accumulates the clock error by a given amount.
  * error = observed clock time - system time [s] */
  void updateClockError(UInt idEpoch, Double deltaClock) {clk.at(idEpoch) += deltaClock;}

  /** @brief Returns the CoM positions in CRF at a given epoch. */
  Vector3d positionCoM(const Time &time) const;

  /** @brief Returns the ARP positions in CRF at a given epoch. */
  Vector3d position(UInt idEpoch, const Time &time) const {return positionCoM(time) + crf2srf.at(idEpoch).inverseTransform(offset.at(idEpoch));}

  /** @brief Returns the velocities in CRF [m/s]. */
  Vector3d velocity(const Time &time) const;

  /** @brief Gets the rotation from CRF to the left-handed ARF. */
  Transform3d celestial2antennaFrame(UInt idEpoch, const Time &/*time*/) const {return srf2arf.at(idEpoch) * crf2srf.at(idEpoch);}
};

/***********************************************/

inline Vector3d GnssTransmitter::positionCoM(const Time &time) const
{
  try
  {
    return Vector3d(polynomial.interpolate({time}, pos));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

inline Vector3d GnssTransmitter::velocity(const Time &time) const
{
  try
  {
    return Vector3d(polynomial.interpolate({time}, vel));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

/// @}

#endif /* __GROOPS___ */
