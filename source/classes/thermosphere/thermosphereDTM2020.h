/***********************************************/
/**
* @file thermosphereDTM2020.h
*
* @brief Density, temperature and velocity.
*
* @author Andreas Kvas
* @date 2022-04-20
*
*/
/***********************************************/

#ifndef __GROOPS_THERMOSPHEREDTM2020__
#define __GROOPS_THERMOSPHEREDTM2020__

// Latex documentation
#ifdef DOCSTRING_Thermosphere
static const char *docstringThermosphereDTM2020 = R"(
\subsection{DTM2020}
Thermosphere parameters from the DTM2020 operational model:

Documentation of the model can be found at https://swami-h2020-eu.github.io/mcm/
and in

Sean Bruinsma, Claude Boniface.
The operational and research DTM-2020 thermosphere models.
Journal of Space Weather and Space Climate, EDP sciences, 2021. 10.1051/swsc/2021032
)";
#endif

/***********************************************/

#include "base/planets.h"
#include "base/polynomial.h"
#include "external/dtm2020/dtm2020.h"
#include "classes/thermosphere/thermosphere.h"
#include "files/fileMatrix.h"

/***** CLASS ***********************************/

/** @brief Density, temperature and velocity.
* @ingroup thermosphereGroup
* @see Thermosphere */
class ThermosphereDTM2020 : public Thermosphere
{
private:
  Matrix f107Data, f107barData, apData;
  Polynomial polynomial;

  Matrix apKpLookupTable;

  Double ap2kp(const Double& ap) const
  {
    UInt idInterval;
    for(idInterval=0; idInterval+1<apKpLookupTable.rows(); idInterval++)
      if(ap < apKpLookupTable(idInterval+1, 0))
        break;

    return apKpLookupTable(idInterval, 1) + (ap - apKpLookupTable(idInterval, 0)) * (apKpLookupTable(idInterval+1, 1) - apKpLookupTable(idInterval, 1)) / (apKpLookupTable(idInterval+1, 0) - apKpLookupTable(idInterval, 0));
  }

public:
  inline ThermosphereDTM2020(Config &config);

  inline void state(const Time &time, const Vector3d &position, Double &density, Double &temperature, Vector3d &velocity) const override;
};

/***********************************************/

inline ThermosphereDTM2020::ThermosphereDTM2020(Config &config)
{
  try
  {
    FileName fileNameDtmConfig, fileNameMagnetic3hAp, fileNameFluxAndAp;

    readConfig(config, "inputfileDtm2020Data",  fileNameDtmConfig,    Config::MUSTSET,  "{groopsDataDir}/thermosphere/dtm2020/DTM_2020_F107_Kp.dat",   "");
    readConfig(config, "inputfileFluxAndAp",    fileNameFluxAndAp,    Config::OPTIONAL, "{groopsDataDir}/thermosphere/dtm2020/inputDTM2020.txt", "");
    readConfig(config, "inputfileMagnetic3hAp", fileNameMagnetic3hAp, Config::OPTIONAL, "{groopsDataDir}/thermosphere/hwm14/apActivity.txt", "indicies for wind model");
    readConfig(config, "hwm14DataDirectory",    fileNameHwm14Path,    Config::OPTIONAL, "{groopsDataDir}/thermosphere/hwm14",                "directory containing dwm07b104i.dat, gd2qd.dat, hwm123114.bin");
    if(isCreateSchema(config)) return;

#ifdef GROOPS_DISABLE_DTM2020
    throw(Exception("Compiled without DTM2020 library"));
#else
    dtm2020initWrapper(fileNameDtmConfig.c_str()); // Load DTM configuration
    Matrix tmp;
    readFileMatrix(fileNameFluxAndAp, tmp);
    std::vector<Time> fluxApTimes;
    fluxApTimes.reserve(tmp.rows());
    for(UInt k = 0; k < tmp.rows(); k++)
      fluxApTimes.push_back(mjd2time(tmp(k, 0)));
    f107Data = tmp.column(1);
    f107barData = tmp.column(2);
    apData = tmp.column(3);
    polynomial.init(fluxApTimes, 1, TRUE/*throw exception*/);

    apKpLookupTable = Matrix(28, 2);
    copy(Vector({0, 2, 3, 4, 5, 6, 7, 9, 12, 15, 18, 22, 27, 32, 39, 48, 56, 67, 80, 94, 111, 132, 154, 179, 207, 236, 300, 400 + 1e-5}), apKpLookupTable.column(0));
    copy(Vector({0, 0.33, 0.67, 1, 1.33, 1.67, 2, 2.33, 2.67, 3, 3.33, 3.67, 4, 4.33, 4.67, 5, 5.33, 5.67, 6, 6.33, 6.67, 7, 7.33, 7.67, 8, 8.33, 8.67, 9}), apKpLookupTable.column(1));
#endif

    magnetic3hAp = InstrumentFile::read(fileNameMagnetic3hAp);
#ifdef GROOPS_DISABLE_HWM14
    if(!fileNameHwm14Path.empty())
      logWarningOnce<<"Compiled without HWM14 wind model sources -> thermospheric wind is not calculated"<<Log::endl;
#endif
 }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

inline void ThermosphereDTM2020::state(const Time &time, const Vector3d &position, Double &density, Double &temperature, Vector3d &velocity) const
{
  try
  {
#ifndef GROOPS_DISABLE_DTM2020
    Ellipsoid ellipsoid;
    Angle     lon, lat;
    Double    height;
    ellipsoid(position, lon, lat, height);

    Time epoch_f107 = time - seconds2time(86400);
    Time epoch_f107bar = time;
    Time epoch_akp_3h = time - seconds2time(3*3600);
    std::vector<Time> timesApMean;
    timesApMean.reserve(24);
    for(UInt k = 0; k<24; k++)
      timesApMean.push_back(epoch_f107 + seconds2time(k * 3600.0));

    F77Float f107 = polynomial.interpolate({epoch_f107}, f107Data)(0, 0);
    F77Float f107bar = polynomial.interpolate({epoch_f107bar}, f107barData)(0, 0);
    F77Float kp3hDelay = ap2kp(polynomial.interpolate({epoch_akp_3h}, apData)(0, 0));
    F77Float meanKpLast24h = ap2kp(mean(polynomial.interpolate(timesApMean, apData)));

    UInt year, month, day, hour, minute;
    Double seconds;
    time2date(time, year, month, day, hour, minute, seconds);
    Time t0 = date2time(year, 1, 1);
    F77Float dayOfYear = (time - t0).mjdInt() + 1;
    F77Float hl = std::fmod(time.mjdMod()*2*PI + lon, 2*PI);

    F77Float ro, tinf, tz, wmm, d[6];
    dtm2020calcWrapper(dayOfYear, f107/*f*/, f107bar/*fbar*/, kp3hDelay/*kp*/, meanKpLast24h, height * 1e-3, hl/*hl*/, lat, lon, tz, tinf, ro, d, wmm);
    density     = 1000*ro; // g/cm^3 -> kg/m^3
    temperature = tz;
    velocity    = wind(time, position);
#endif
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

#endif
