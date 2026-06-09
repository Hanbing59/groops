/***********************************************/
/**
* @file sentinel2StarCamera.cpp
*
* @brief read Sentinel-1/2/3/6 star camera data.
*
* @author Barbara Suesser-Rechberger
* @date 2020-10-21
*
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
This program reads in Sentinel-1/2/3/6 star camera data, i.e. AUX_PROQUA quaternions
(each line: YYYY/MM/DD HH:MM:SS.SSS  Q_COMPR QCOMP1 QCOMP2 QCOMP3 ...)
and then writes them to an GROOPS instrument file.
)";

/***********************************************/

#include "programs/program.h"
#include "inputOutput/system.h"
#include "inputOutput/file.h"
#include "files/fileInstrument.h"
#include "base/string.h"

/***** CLASS ***********************************/

/** @brief read Sentinel-1/2/3/6 star camera data.
* @ingroup programsConversionGroup */
class Sentinel2StarCamera
{
public:
  void run(Config &config, Parallel::CommunicatorPtr comm);
};

GROOPS_REGISTER_PROGRAM(Sentinel2StarCamera, SINGLEPROCESS, "read Sentinel-1/2/3/6 star camera data", Conversion, Instrument)

/***********************************************/

void Sentinel2StarCamera::run(Config &config, Parallel::CommunicatorPtr /*comm*/)
{
  try
  {
    FileName              outNameStar;
    std::vector<FileName> fileNamesIn;

    readConfig(config, "outputfileStarCamera", outNameStar, Config::MUSTSET,  "", "");
    readConfig(config, "inputfile",            fileNamesIn, Config::MUSTSET,  "", "Supports wildcards: * and ?");
    if(isCreateSchema(config)) return;

    logStatus<<"read input files"<<Log::endl;
    StarCameraArc arc;

    for(const auto &fileName : fileNamesIn)
    {
      try
      {
        logStatus<<"process file <"<<fileName<<">"<<Log::endl;
        // Expand wildcards
        const std::vector<FileName> fileList = System::fileList(fileName);
        for(const auto &fileOne : fileList)
        {
          logStatus<<"   read file <"<<fileOne<<">"<<Log::endl;
          InFile file(fileOne);
          std::string line;
          while(std::getline(file, line))
          {
            if(line.empty())
              continue;

            std::string lineID1 = line.substr(0,1);
            if(lineID1 == "#")
              continue;

            // Time format in AUX_PROQUA file: YYYY/MM/DD HH:MM:SS.SSS  Q_COMPR QCOMP1 QCOMP2 QCOMP3
            UInt   year   = String::toInt(line.substr(0, 4));
            UInt   month  = String::toInt(line.substr(5, 2));
            UInt   day    = String::toInt(line.substr(8, 2));
            UInt   hour   = String::toInt(line.substr(11, 2));
            UInt   minute = String::toInt(line.substr(14, 2));
            Double second = String::toDouble(line.substr(17, 6));

            std::stringstream ss(line.substr(23));
            ss.exceptions(std::ios::badbit | std::ios::failbit);

            Vector q(4);
            ss>>q(0)>>q(1)>>q(2)>>q(3);

            StarCameraEpoch epoch;
            // in GPS time
            epoch.time   = date2time(year, month, day, hour, minute, second);
            epoch.rotary = Rotary3d(q);
            arc.push_back(epoch);
          }
        }
      }
      catch(std::exception &e)
      {
        logError<<e.what()<<": continue..."<<Log::endl;
      }
    }

    logStatus<<"write star camera data to file <"<<outNameStar<<">"<<Log::endl;
    Arc::printStatistics(arc);
    InstrumentFile::write(outNameStar, arc);
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
