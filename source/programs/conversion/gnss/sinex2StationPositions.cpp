/***********************************************/
/**
* @file sinex2StationPosition.cpp
*
* @brief Convert station positions from SINEX to InstrumentVector3d.
*
* @author Sebastian Strasser
* @author Torsten Mayer-Guerr
* @date 2016-12-07
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
Extracts station positions from \config{inputfileSinexSolution}
(\href{http://www.iers.org/IERS/EN/Organization/AnalysisCoordinator/SinexFormat/sinex.html}{SINEX format description})
and writes an \configFile{outputfileInstrument}{instrument} of type VECTOR3D
for each station.  Positions will be computed at \configClass{timeSeries}{timeSeriesType} based on position and velocity
of each provided interval in the SINEX file.
With \config{inputfileSinexDiscontinuities} the bounds of these time spans are adjusted to the exact epochs of discontinuities.
The \config{inputfileSinexPostSeismicDeformations} adds the ITRF post-seismic deformation model to the affected stations.
The \config{inputfileSinexFrequencies} adds annual and semi-annual frequencies.

If \config{extrapolateBackward} or \config{extrapolateForward} are provided, positions will also be computed for epochs
before the first interval/after the last interval, based on the position and velocity of the first/last interval.
Position extrapolation will stop at the first discontinuity before the first interval/after the last interval.

Stations can be limited via \config{stationName}, otherwise all stations in \config{inputfileSinexSolution} will be used.
)";

/***********************************************/

#include "programs/program.h"
#include "base/string.h"
#include "classes/timeSeries/timeSeries.h"
#include "inputOutput/fileSinex.h"
#include "files/fileInstrument.h"

/***** CLASS ***********************************/

/** @brief Convert GNSS station positions from SINEX to InstrumentVector3d.
* @ingroup programsConversionGroup */
class Sinex2StationPositions
{
public:
  void run(Config &config, Parallel::CommunicatorPtr comm);
};

GROOPS_REGISTER_PROGRAM(Sinex2StationPositions, SINGLEPROCESS, "Convert station positions from SINEX to InstrumentVector3d.", Conversion, Gnss, Instrument)

/***********************************************/

void Sinex2StationPositions::run(Config &config, Parallel::CommunicatorPtr /*comm*/)
{
  try
  {
    FileName                 fileNameInstrument;
    std::string              variableLoopStation;
    FileName                 fileNameSinex, fileNameSinexDisc, fileNameSinexPsd, fileNameSinexFreq;
    TimeSeriesPtr            timeSeries;
    Bool                     extrapolateForward, extrapolateBackward;
    std::vector<std::string> stationNames;

    readConfig(config, "outputfileInstrument",                  fileNameInstrument,  Config::MUSTSET,   "stationPosition.{station}.dat", "loop variable is replaced with station name (e.g. wtzz)");
    readConfig(config, "variableLoopStation",                   variableLoopStation, Config::DEFAULT,   "station", "variable name for station loop");
    readConfig(config, "inputfileSinexSolution",                fileNameSinex,       Config::MUSTSET,   "ITRF2020-IGS-TRF.SSC",            "SINEX file");
    readConfig(config, "inputfileSinexDiscontinuities",         fileNameSinexDisc,   Config::OPTIONAL,  "ITRF2020-soln-gnss.snx",          "SINEX file");
    readConfig(config, "inputfileSinexPostSeismicDeformations", fileNameSinexPsd,    Config::OPTIONAL,  "ITRF2020-psd-gnss.snx",           "SINEX file");
    readConfig(config, "inputfileSinexFrequencies",             fileNameSinexFreq,   Config::OPTIONAL,  "ITRF2020-Frequencies-XYZ-CF.snx", "SINEX file (XYZ or ENU)");
    readConfig(config, "timeSeries",                            timeSeries,          Config::MUSTSET,   "",  "compute positions for these epochs based on velocity");
    readConfig(config, "extrapolateForward",                    extrapolateForward,  Config::DEFAULT,   "1", "also compute positions for epochs after last interval defined in SINEX file");
    readConfig(config, "extrapolateBackward",                   extrapolateBackward, Config::DEFAULT,   "0", "also compute positions for epochs before first interval defined in SINEX file");
    readConfig(config, "stationName",                           stationNames,        Config::OPTIONAL,  "",  "convert only these stations");
    if(isCreateSchema(config)) return;

    // Change input station names to lower case for case-insensitive comparison with SINEX station names.
    for(std::string &name : stationNames)
      name = String::lowerCase(String::trim(name));

    struct Interval
    {
      std::string pointCode;
      std::string solutionId;
      Time        referenceTime, timeStart, timeEnd;
      Vector3d    position, velocity;
      Vector3dArc arc;
      std::map<std::string, Vector3d> ampl3d;
    };

    struct Station
    {
      std::vector<Interval> intervals;
      std::vector<Interval> discontinuities;
      std::string           pointCode, dome;
      // approximate position from the SITE/ID block of the SINEX solution file
      Vector3d              position;
      Transform3d           lnof2trf;
    };

    std::map<std::string, std::vector<Station>> stations;
    std::vector<Time> times = timeSeries->times();

    // =========================================================

    if(!fileNameSinex.empty())
    {
      logStatus<<"read SINEX solution file <"<<fileNameSinex<<">"<<Log::endl;
      Sinex sinex;
      readFileSinex(fileNameSinex, sinex);

      // SITE/ID
      // -------
      for(std::string &line : sinex.findBlock("SITE/ID")->lines)
      {
        // *         1         2         3         4         5         6         7         8
        // *12345678901234567890123456789012345678901234567890123456789012345678901234567890
        // *Code PT __DOMES__ T _STATION DESCRIPTION__ APPROX_LON_ APPROX_LAT_ _APP_H_
        const std::string name = String::lowerCase(String::trim(line.substr(1, 4)));
        if(stationNames.size() && std::find(stationNames.begin(), stationNames.end(), name) == stationNames.end())
          continue;

        const std::string pointCode = String::trim(line.substr(6, 2));
        // 4-char station name + 2-char point code should be unique
        const auto it = std::find_if(stations[name].begin(), stations[name].end(), [&](const Station &s) {return s.pointCode == pointCode;});
        if(it != stations[name].end())
        {
          logWarning<<"SITE/ID duplicate: "<<line<<Log::endl;
          continue;
        }
        const Double longitude   = String::toDouble(line.substr(44, 3)) + String::toDouble(line.substr(48, 2))/60 + String::toDouble(line.substr(51, 4))/3600;
        const Double latitude    = String::toDouble(line.substr(56, 3)) + (String::startsWith(String::trim(line.substr(56, 3)), "-") ? -1 : 1) * std::fabs(String::toDouble(line.substr(60, 2))/60 + String::toDouble(line.substr(62, 5))/3600);
        const Double height      = String::toDouble(line.substr(68, 7));
        Station station;
        station.pointCode = pointCode;
        station.dome      = String::trim(line.substr(9, 9));
        station.position  = Ellipsoid()(Angle(DEG2RAD*longitude), Angle(DEG2RAD*latitude), height);
        station.lnof2trf  = localNorthEastUp(polar(Angle(DEG2RAD*longitude), Angle(DEG2RAD*latitude), 1.));
        stations[name].push_back(station);
      }

      // SOLUTION/EPOCHS
      // ---------------
      for(std::string &line : sinex.findBlock("SOLUTION/EPOCHS")->lines)
      {
        const std::string name = String::lowerCase(String::trim(line.substr(1, 4)));
        if(stationNames.size() && std::find(stationNames.begin(), stationNames.end(), name) == stationNames.end())
          continue;
        if(stations.find(name) == stations.end())
        {
          logWarning<<"SOLUTION/EPOCHS without matching SITE/ID 4-char name: "<<line<<Log::endl;
          continue;
        }

        const std::string pointCode  = String::trim(line.substr(6, 2));
        const auto it = std::find_if(stations[name].begin(), stations[name].end(), [&](const Station &s) {return s.pointCode == pointCode;});
        if(it == stations[name].end())
        {
          logWarning<<"SOLUTION/EPOCHS without matching SITE/ID point code: "<<line<<Log::endl;
          continue;
        }
        Interval interval;
        interval.pointCode  = pointCode;
        interval.solutionId = String::trim(line.substr(9, 4));
        interval.timeStart  = Sinex::str2time(line, 16, FALSE);
        interval.timeEnd    = Sinex::str2time(line, 29, TRUE);
        it->intervals.push_back(interval);
      }

      // SOLUTION/ESTIMATE
      // -----------------
      UInt count = 0;
      for(std::string &line : sinex.findBlock("SOLUTION/ESTIMATE")->lines)
      {
        // *         1         2         3         4         5         6         7         8
        // *12345678901234567890123456789012345678901234567890123456789012345678901234567890
        // *INDEX _TYPE_ CODE PT SOLN _REF_EPOCH__ UNIT S ___ESTIMATED_VALUE___ __STD_DEV__
        //      1 STAX   7894  A    1 15:001:00000 m    2 -.219677819584269E+07 0.11095E+00
        const std::string name = String::lowerCase(String::trim(line.substr(14, 4)));
        if(stationNames.size() && std::find(stationNames.begin(), stationNames.end(), name) == stationNames.end())
          continue;

        const std::string parameterType = String::trim(line.substr(7, 6));
        const std::string pointCode     = String::trim(line.substr(19, 2));
        const std::string solutionId    = String::trim(line.substr(22, 4));
        const Double      value         = String::toDouble(line.substr(47, 21));

        if((parameterType != "STAX") && (parameterType != "STAY") && (parameterType != "STAZ") &&
           (parameterType != "VELX") && (parameterType != "VELY") && (parameterType != "VELZ"))
          continue;
        if(stations.find(name) == stations.end())
        {
          logWarning<<"SOLUTION/ESTIMATE without matching SITE/ID 4-char name: "<<line<<Log::endl;
          continue;
        }
        const auto station = std::find_if(stations[name].begin(), stations[name].end(), [&](const Station &s) {return s.pointCode == pointCode;});
        if(station == stations[name].end())
        {
          logWarning<<"SOLUTION/ESTIMATE without matching SITE/ID point code: "<<line<<Log::endl;
          continue;
        }
        const auto interval = std::find_if(station->intervals.begin(), station->intervals.end(), [&](const Interval &i) {return i.solutionId == solutionId;});
        if(interval == station->intervals.end())
        {
          logWarning<<"SOLUTION/ESTIMATE without matching SOLUTION/EPOCHS solution: "<<line<<Log::endl;
          continue;
        }

        interval->referenceTime = Sinex::str2time(line, 27, FALSE);
        Bool used = TRUE;
        if(     parameterType == "STAX") interval->position.x() = value;
        else if(parameterType == "STAY") interval->position.y() = value;
        else if(parameterType == "STAZ") interval->position.z() = value;
        else if(parameterType == "VELX") interval->velocity.x() = value;
        else if(parameterType == "VELY") interval->velocity.y() = value;
        else if(parameterType == "VELZ") interval->velocity.z() = value;
        else
          used = FALSE;
        count += used;
      }
      logInfo<<"  "<<count<<" of "<<sinex.findBlock("SOLUTION/ESTIMATE")->lines.size()<<" parameters used"<<Log::endl;
    }

    // =========================================================

    // test intervals
    // --------------
    for(auto &station : stations)
      for(auto &point : station.second)
      {
        auto &intervals = point.intervals;
        if(!intervals.size())
        {
          logWarning<<station.first<<" "<<point.pointCode<<" without solution intervals"<<Log::endl;
          continue;
        }

        std::stable_sort(point.intervals.begin(), point.intervals.end(), [](auto &a, auto &b) {return a.timeEnd < b.timeEnd;});

        // repair overlapping periods
        for(UInt i=0; i<intervals.size()-1; i++)
          if(intervals.at(i).timeEnd > intervals.at(i+1).timeStart)
          {
            if((intervals.at(i).timeEnd-intervals.at(i+1).timeStart).mjd() > 0.25)
              logWarning<<station.first<<" "<<point.pointCode<<" interval ends at "<<intervals.at(i).timeEnd.dateTimeStr()<<", next intervals starts at "<<intervals.at(i+1).timeStart.dateTimeStr()<<Log::endl;
            intervals.at(i).timeEnd = intervals.at(i+1).timeStart = 0.5*(intervals.at(i).timeEnd+intervals.at(i+1).timeStart);
          }

        // fill gaps
        for(UInt i=0; i<intervals.size()-1; i++)
          if(intervals.at(i).timeEnd < intervals.at(i+1).timeStart)
            intervals.at(i).timeEnd = intervals.at(i+1).timeStart;

        if(extrapolateBackward)
          intervals.front().timeStart = Time();
        if(extrapolateForward)
          intervals.back().timeEnd = date2time(2500, 1, 1);
      }

    // =========================================================

    if(!fileNameSinexDisc.empty())
    {
      logStatus<<"read SINEX discontinuities file <"<<fileNameSinexDisc<<">"<<Log::endl;
      Sinex sinex;
      readFileSinex(fileNameSinexDisc, sinex);

      // SOLUTION/DISCONTINUITY
      // ----------------------
      for(std::string &line : sinex.findBlock("SOLUTION/DISCONTINUITY")->lines)
      {
        // *         1         2         3         4         5         6         7         8
        // *12345678901234567890123456789012345678901234567890123456789012345678901234567890
        //  ADEA  A    1 D 00:000:00000 98:084:11545 P - 1998/03/25-Earthquake_M8.1_Balleny_Islands_(616km)
        const std::string name = String::lowerCase(String::trim(line.substr(1, 4)));
        if(stationNames.size() && std::find(stationNames.begin(), stationNames.end(), name) == stationNames.end())
          continue;
        if(stations.find(name) == stations.end())
          continue;

        const std::string pointCode = String::trim(line.substr(6, 2));
        const auto station = std::find_if(stations[name].begin(), stations[name].end(), [&](const Station &s) {return s.pointCode == pointCode;});
        if(station == stations[name].end())
        {
          logWarning<<"SOLUTION/DISCONTINUITY without matching SITE/ID point code: "<<line<<Log::endl;
          continue;
        }
        const std::string solutionId = String::trim(line.substr(9, 4));
        const auto interval = std::find_if(station->intervals.begin(), station->intervals.end(), [&](const Interval &i) {return i.solutionId == solutionId;});
        if(interval == station->intervals.end())
        {
          logWarning<<"SOLUTION/DISCONTINUITY without matching SOLUTION/EPOCHS solution: "<<line<<Log::endl;
          continue;
        }

        Interval discontinuity;
        discontinuity.pointCode  = pointCode;
        discontinuity.solutionId = solutionId;
        discontinuity.timeStart  = Sinex::str2time(line, 16, FALSE);
        discontinuity.timeEnd    = Sinex::str2time(line, 29, TRUE);

        if(line.substr(42,1) == "P") // position changes only
          station->discontinuities.push_back(discontinuity);
      }

      // adjust intervals
      // ----------------
      for(auto &station : stations)
        for(auto &point : station.second)
        {
          auto &intervals       = point.intervals;
          auto &discontinuities = point.discontinuities;
          if(!discontinuities.size())
            continue;

          std::stable_sort(discontinuities.begin(), discontinuities.end(), [](auto &a, auto &b) {return a.timeEnd < b.timeEnd;});

          // repair non-connected periods
          for(UInt i=0; i<discontinuities.size()-1; i++)
            if(discontinuities.at(i).timeEnd != discontinuities.at(i+1).timeStart)
            {
              logWarning<<station.first<<" "<<point.pointCode<<" Non-connected discontinuities: one ends at "
                        <<discontinuities.at(i).timeEnd.dateTimeStr()<<" while the next starts at "
                        <<discontinuities.at(i+1).timeStart.dateTimeStr()<<Log::endl;
              discontinuities.at(i).timeEnd = discontinuities.at(i+1).timeStart = 0.5*(discontinuities.at(i).timeEnd+discontinuities.at(i+1).timeStart);
            }

          // adjust solution intervals
          Time timeStart = intervals.front().timeStart;
          Time timeEnd   = intervals.back().timeEnd;
          for(UInt i=0; i<intervals.size(); i++)
          {
            auto discontinuity = std::find_if(discontinuities.begin(), discontinuities.end(), [&] (auto &x) {return x.solutionId == intervals.at(i).solutionId;});
            if(discontinuity != discontinuities.end())
            {
              // this should never happen
              if(intervals.at(i).pointCode != discontinuity->pointCode)
                logWarning<<station.first<<" point code differ "<<intervals.at(i).pointCode<<" != "<<discontinuity->pointCode<<Log::endl;
              intervals.at(i).pointCode = discontinuity->pointCode;
              intervals.at(i).timeStart = std::max(timeStart, discontinuity->timeStart);
              intervals.at(i).timeEnd   = std::min(timeEnd,   discontinuity->timeEnd);
            }
            else
              logWarning<<station.first<<" "<<point.pointCode<<" solution '"<<intervals.at(i).solutionId<<"' without discontinuity interval"<<Log::endl;
          }
        }
    }

    // =========================================================

    // create arcs
    // -----------
    for(auto &station : stations)
      for(auto &point : station.second)
        for(auto &interval : point.intervals)
        {
          logInfo<<station.first<<" "<<point.pointCode<<" solution "<<interval.solutionId
                <<" from "<<interval.timeStart.dateTimeStr()<<" to "<<interval.timeEnd.dateTimeStr()<<Log::endl;
          for(auto &time : times)
            if(time.isInInterval(interval.timeStart, interval.timeEnd))
            {
              Vector3dEpoch epoch;
              epoch.time     = time;
              epoch.vector3d = interval.position + interval.velocity * (time-interval.referenceTime).mjd()/365.25;
              interval.arc.push_back(epoch);
            }
        }


    // =========================================================

    if(!fileNameSinexPsd.empty())
    {
      logStatus<<"read SINEX post seismic deformation file <"<<fileNameSinexPsd<<">"<<Log::endl;
      Sinex sinex;
      readFileSinex(fileNameSinexPsd, sinex);

      // SITE/ID
      // -------
      for(std::string &line : sinex.findBlock("SITE/ID")->lines)
      {
        // *         1         2         3         4         5         6         7         8
        // *12345678901234567890123456789012345678901234567890123456789012345678901234567890
        // *Code PT __DOMES__ T _STATION DESCRIPTION__ APPROX_LON_ APPROX_LAT_ _APP_H_
        const std::string name = String::lowerCase(String::trim(line.substr(1, 4)));
        if(stationNames.size() && std::find(stationNames.begin(), stationNames.end(), name) == stationNames.end())
          continue;
        if(stations.find(name) == stations.end())
          continue;

        const std::string pointCode = String::trim(line.substr(6, 2));
        const auto station = std::find_if(stations[name].begin(), stations[name].end(), [&](const Station &s) {return s.pointCode == pointCode;});
        if(station == stations[name].end())
          continue;

        const Double longitude = String::toDouble(line.substr(44, 3)) + String::toDouble(line.substr(48, 2))/60 + String::toDouble(line.substr(51, 4))/3600;
        const Double latitude  = String::toDouble(line.substr(56, 3)) + (String::startsWith(String::trim(line.substr(56, 3)), "-") ? -1 : 1) * std::fabs(String::toDouble(line.substr(60, 2))/60 + String::toDouble(line.substr(62, 5))/3600);
        const Double height    = String::toDouble(line.substr(68, 7));
        const Double dist      = (station->position - Ellipsoid()(Angle(DEG2RAD*longitude), Angle(DEG2RAD*latitude), height)).r();
        if(dist > 10)
          logWarning<<name<<" "<<pointCode<<": approx position differs with that from the solution file by "<<dist<<" m"<<Log::endl;
        station->lnof2trf = localNorthEastUp(polar(Angle(DEG2RAD*longitude), Angle(DEG2RAD*latitude), 1.));
      }

      UInt count = 0;
      const std::vector<std::string> &lines = sinex.findBlock("SOLUTION/ESTIMATE")->lines;
      for(UInt i=0; i<lines.size(); i+=2)
      {
        // *         1         2         3         4         5         6         7         8
        // *12345678901234567890123456789012345678901234567890123456789012345678901234567890
        // *INDEX _TYPE_ CODE PT SOLN _REF_EPOCH__ UNIT S ___ESTIMATED_VALUE___ __STD_DEV__
        //      1 ALOG_E AREB  A ---- 01:174:73993 m    2 -1.63377151539946e-01 1.85315e-02
        //      2 TLOG_E AREB  A ---- 01:174:73993 m    2  8.15228736135023e+00 1.23120e+00
        const std::string name = String::lowerCase(String::trim(lines.at(i).substr(14, 4)));
        if(stationNames.size() && std::find(stationNames.begin(), stationNames.end(), name) == stationNames.end())
          continue;
        if(stations.find(name) == stations.end())
          continue;
        const std::string pointCode  = String::trim(lines.at(i).substr(19, 2));
        const auto station = std::find_if(stations[name].begin(), stations[name].end(), [&](const Station &s) {return s.pointCode == pointCode;});
        if(station == stations[name].end())
          continue;

        // parameter type
        const Bool isLog = (String::upperCase(lines.at(i  ).substr(7, 5)) == "ALOG_") &&
                           (String::upperCase(lines.at(i+1).substr(7, 5)) == "TLOG_");
        if(!isLog && ((String::upperCase(lines.at(i  ).substr(7, 5)) != "AEXP_") ||
                      (String::upperCase(lines.at(i+1).substr(7, 5)) != "TEXP_")))
          continue; // other parameter type

        const Double A     = String::toDouble(lines.at(i  ).substr(47, 21));
        const Double tau   = String::toDouble(lines.at(i+1).substr(47, 21)) * 365.25;
        const Time   time0 = Sinex::str2time(lines.at(i), 27, FALSE);

        Vector3d ampl3d;
        const Char c = String::upperCase(lines.at(i).substr(12, 1)).at(0);
        if(     c == 'N')                 ampl3d = station->lnof2trf.transform(Vector3d(A, 0, 0));
        else if(c == 'E')                 ampl3d = station->lnof2trf.transform(Vector3d(0, A, 0));
        else if((c == 'U') || (c == 'H')) ampl3d = station->lnof2trf.transform(Vector3d(0, 0, A));

        Bool used = FALSE;
        for(auto &interval : station->intervals)
          for(UInt i=0; i<interval.arc.size(); i++)
          {
            if(interval.arc.at(i).time >= time0)
            {
              used = TRUE;
              if(isLog)
                interval.arc.at(i).vector3d += ampl3d * std::log(1.+(interval.arc.at(i).time-time0).mjd()/tau);
              else
                interval.arc.at(i).vector3d += ampl3d * (1.-std::exp(-(interval.arc.at(i).time-time0).mjd()/tau));
            }
          }
        if(!used)
          logWarning<<"Unused: "<<lines.at(i)<<Log::endl;
        count += 2*used;
      }
      logInfo<<"  "<<count<<" of "<<sinex.findBlock("SOLUTION/ESTIMATE")->lines.size()<<" parameters used"<<Log::endl;
    }

    // =========================================================

    if(!fileNameSinexFreq.empty())
    {
      logStatus<<"read SINEX frequencies file <"<<fileNameSinexFreq<<">"<<Log::endl;
      Sinex sinex;
      readFileSinex(fileNameSinexFreq, sinex);

      for(std::string &line : sinex.findBlock("SITE/ID")->lines)
      {
        // *         1         2         3         4         5         6         7         8
        // *12345678901234567890123456789012345678901234567890123456789012345678901234567890
        // *Code PT __DOMES__ T _STATION DESCRIPTION__ APPROX_LON_ APPROX_LAT_ _APP_H_
        const std::string name = String::lowerCase(String::trim(line.substr(1, 4)));
        if(stationNames.size() && std::find(stationNames.begin(), stationNames.end(), name) == stationNames.end())
          continue;
        if(stations.find(name) == stations.end())
          continue;

        const std::string pointCode = String::trim(line.substr(6, 2));
        const auto station = std::find_if(stations[name].begin(), stations[name].end(), [&](const Station &s) {return s.pointCode == pointCode;});
        if(station == stations[name].end())
          continue;

        const Double longitude = String::toDouble(line.substr(44, 3)) + String::toDouble(line.substr(48, 2))/60 + String::toDouble(line.substr(51, 4))/3600;
        const Double latitude  = String::toDouble(line.substr(56, 3)) + (String::startsWith(String::trim(line.substr(56, 3)), "-") ? -1 : 1) * std::fabs(String::toDouble(line.substr(60, 2))/60 + String::toDouble(line.substr(62, 5))/3600);
        const Double height    = String::toDouble(line.substr(68, 7));
        const std::string dome = String::trim(line.substr(9, 9));
        const Double dist      = (station->position - Ellipsoid()(Angle(DEG2RAD*longitude), Angle(DEG2RAD*latitude), height)).r();
        if(dist > 10)
          logWarning<<name<<" "<<pointCode<<": approx position differs with that from the solution file by "<<dist<<" m"<<Log::endl;
        // stations[name].lnof2trf = localNorthEastUp(polar(Angle(DEG2RAD*longitude), Angle(DEG2RAD*latitude), 1.));
      }

      UInt count = 0;
      for(const std::string &line : sinex.findBlock("SOLUTION/ESTIMATE")->lines)
      {
        // *         1         2         3         4         5         6         7         8
        // *12345678901234567890123456789012345678901234567890123456789012345678901234567890
        // *INDEX _TYPE_ CODE PT SOLN _REF_EPOCH__ UNIT S ___ESTIMATED_VALUE___ __STD_DEV__
        //      1 A1COSX ALBH  A    1 15:001:00000 m    2  2.77340030328052e-05 1.64130e-04
        const std::string name = String::lowerCase(String::trim(line.substr(14, 4)));
        if(stationNames.size() && std::find(stationNames.begin(), stationNames.end(), name) == stationNames.end())
          continue;
        if(stations.find(name) == stations.end())
          continue;
        const std::string pointCode  = String::trim(line.substr(19, 2));
        const auto station = std::find_if(stations[name].begin(), stations[name].end(), [&](const Station &s) {return s.pointCode == pointCode;});
        if(station == stations[name].end())
          continue;

        const std::string solutionId = String::trim(line.substr(22, 4));
        const auto interval = std::find_if(station->intervals.begin(), station->intervals.end(), [&](const Interval &i) {return i.solutionId == solutionId;});
        if(interval == station->intervals.end())
          {
            logWarning<<"SOLUTION/ESTIMATE without matching SOLUTION/EPOCHS solution: "<<line<<Log::endl;
            continue;
          }

        // parameter type
        std::string para = String::upperCase(line.substr(7, 5));
        if(!String::startsWith(para, "A") || !(String::endsWith(para, "COS") || String::endsWith(para, "SIN")))
          continue;

        Vector3d ampl3d;
        const Double A = String::toDouble(line.substr(47, 21));
        const Char   c = String::upperCase(line.substr(12, 1)).at(0);
        if(     c == 'X') ampl3d = Vector3d(A, 0, 0);
        else if(c == 'Y') ampl3d = Vector3d(0, A, 0);
        else if(c == 'Z') ampl3d = Vector3d(0, 0, A);
        else if(c == 'N') ampl3d = station->lnof2trf.transform(Vector3d(A, 0, 0));
        else if(c == 'E') ampl3d = station->lnof2trf.transform(Vector3d(0, A, 0));
        else if(c == 'U') ampl3d = station->lnof2trf.transform(Vector3d(0, 0, A));
        interval->ampl3d[para] += ampl3d;
        count++;
      }
      logInfo<<"  "<<count<<" of "<<sinex.findBlock("SOLUTION/ESTIMATE")->lines.size()<<" parameters used"<<Log::endl;

      for(auto &station : stations)
        for(auto &point : station.second)
        {
          // find oscillations in one of the intervals
          std::map<std::string, Vector3d> ampl3d;
          for(auto &interval : point.intervals)
            if(interval.ampl3d.size())
              ampl3d = interval.ampl3d;
          // check if oscillations in all intervals the same
          Bool same = TRUE;
          for(auto &interval : point.intervals)
            if(interval.ampl3d.size())
              same = same && std::equal(ampl3d.begin(), ampl3d.end(), interval.ampl3d.begin(), [](auto &a, auto &b)
                                      {return (a.first == b.first) && ((a.second-b.second).r() < 0.0001);});
          if(same) // can extend to missing intervals
            for(auto &interval : point.intervals)
              if(!interval.ampl3d.size())
                interval.ampl3d = ampl3d;

          // check for empty intervals
          for(auto &interval : point.intervals)
            if(!interval.ampl3d.size())
              logWarning<<station.first<<" "<<point.pointCode<<" interval "<<interval.timeStart.dateTimeStr()<<" -- "<<interval.timeEnd.dateTimeStr()<<" without frequencies"<<Log::endl;

          // apply oscillations
          for(auto &interval : point.intervals)
            for(auto &A : interval.ampl3d)
            {
              const Bool isCos   = (A.first.substr(2, 3) == "COS");
              const Double omega = 2*PI*String::toDouble(A.first.substr(1, 1)); // annual or semiannual
              if(!isCos && (A.first.substr(2, 3) != "SIN"))
                throw(Exception("COS/SIN expected"));
              for(UInt i=0; i<interval.arc.size(); i++)
                interval.arc.at(i).vector3d += A.second * (isCos ? std::cos(omega*interval.arc.at(i).time.decimalYear()) : std::sin(omega*interval.arc.at(i).time.decimalYear()));
            }
        }
    } // if(!fileNameSinexFreq.empty())

    // =========================================================

    VariableList fileNameVariableList;
    fileNameVariableList.setVariable(variableLoopStation, "****");
    logStatus<<"write instrument files <"<<fileNameInstrument(fileNameVariableList)<<">"<<Log::endl;
    UInt count = 0;
    for(auto &station : stations)
      for(auto &point : station.second)
      {
        if(&point == &station.second.back())
          fileNameVariableList.setVariable(variableLoopStation, station.first);
        else
          fileNameVariableList.setVariable(variableLoopStation, station.first+"_"s+point.pointCode);

        std::vector<Vector3dArc> arcs;
        for(auto &interval : point.intervals)
          if(interval.arc.size())
            arcs.push_back(interval.arc);
        if(arcs.size())
        {
          Double dist = 0;
          for(auto &arc : arcs)
            for(auto &epoch : arc)
              dist = std::max(dist, (epoch.vector3d - point.position).r());
          if(dist > 10)
            logWarning<<station.first<<": max diff of approx position with its solution reaches "<<dist<<" m"<<Log::endl;

          InstrumentFile::write(fileNameInstrument(fileNameVariableList), arcs);
          count++;
          logStatus<<"  "<<count%"%3i"s<<" "<<station.first<<" "<<point.pointCode<<" "<<fileNameInstrument(fileNameVariableList)<<Log::endl;
        }
      }
    logInfo<<"  "<<count<<" station files written"<<Log::endl;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
