/***********************************************/
/**
* @file miscAccelerationsAntennaThrust.h
*
* @brief Antenna thrust (simple model).
* @see MiscAccelerations
*
* @author Sebastian Strasser
* @date 2015-03-04
*
*/
/***********************************************/

#ifndef __GROOPS_MISCACCELERATIONSANTENNATHRUST__
#define __GROOPS_MISCACCELERATIONSANTENNATHRUST__

// Latex documentation
#ifdef DOCSTRING_MiscAccelerations
static const char *docstringMiscAccelerationsAntennaThrust = R"(
\subsection{Antenna thrust}\label{miscAccelerationsType:antennaThrust}
The thrust (acceleration) in the opposite direction the antenna is facing
which is generated by satellite antenna broadcasts.
The thrust is defined in the satellite macro model.
)";
#endif

/***********************************************/

#include "classes/miscAccelerations/miscAccelerations.h"

/***** CLASS ***********************************/

/** @brief Antenna thrust (simple model).
* @ingroup miscAccelerationsGroup
* @see MiscAccelerations */
class MiscAccelerationsAntennaThrust : public MiscAccelerationsBase
{
  Double factor;

public:
  MiscAccelerationsAntennaThrust(Config &config);

  Vector3d acceleration(SatelliteModelPtr satellite, const Time &time, const Vector3d &position, const Vector3d &velocity,
                        const Rotary3d &rotSat, const Rotary3d &rotEarth, EphemeridesPtr ephemerides) override;
};

/***********************************************/

inline MiscAccelerationsAntennaThrust::MiscAccelerationsAntennaThrust(Config &config)
{
  try
  {

    readConfig(config, "factor", factor, Config::DEFAULT, "1.0", "the result is multplied by this factor");
    if(isCreateSchema(config)) return;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

inline Vector3d MiscAccelerationsAntennaThrust::acceleration(SatelliteModelPtr satellite, const Time &/*time*/,
                                                             const Vector3d &/*pos*/, const Vector3d &/*vel*/,
                                                             const Rotary3d &rotSat, const Rotary3d &rotEarth, EphemeridesPtr /*ephemerides*/)
{
  try
  {
    if(!satellite)
      throw(Exception("No satellite model given"));

    return factor * rotEarth.rotate(rotSat.rotate(satellite->accelerationThrust()));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

#endif