/***********************************************/
/**
* @file iersPotential2DoodsonHarmonics.cpp
*
* @brief Read ocean tide file in IERS format.
*
* @author Torsten Mayer-Guerr
* @date 2019-11-02
*
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
Read ocean tide file in IERS format.
)";

/***********************************************/

#include "programs/program.h"
#include "base/doodson.h"
#include "files/fileDoodsonHarmonic.h"

/***** CLASS ***********************************/

/** @brief Reads ocean tide file in IERS format and writes out a \b DoodsonHarmonic file.
 * 
 * This is the inverse conversion of the program \b DoodsonHarmonics2IersPotential.
 * The input file contains the geopotential harmonic amplitudes for different tide constituents 
 * of an ocean tide model (e.g. FES2004). Its data format should refer to that provided 
 * by the IERS Convention 2010 for the model FES2004 at 
 * ftp://tai.bipm.org/iers/conv2010/chapter6/tidemodels/fes2004_Cnm-Snm.dat.
 * Refer to the IERS Conventions 2010, Section 6.3, for details.
* @ingroup programsConversionGroup */
class IersPotential2DoodsonHarmonics
{
public:
  void run(Config &config, Parallel::CommunicatorPtr comm);
};

GROOPS_REGISTER_PROGRAM(IersPotential2DoodsonHarmonics, SINGLEPROCESS, "Read ocean tide file in IERS format", Conversion, DoodsonHarmonics)

/***********************************************/

void IersPotential2DoodsonHarmonics::run(Config &config, Parallel::CommunicatorPtr /*comm*/)
{
  try
  {
    FileName  fileNameOut;
    FileName  fileNameIn, fileNameTGP;
    UInt      countHeader;
    Double    GM, R;
    UInt      minDegree, maxDegree = INFINITYDEGREE;

    renameDeprecatedConfig(config, "outputfileDoodsonHarmoncis", "outputfileDoodsonHarmonics", date2time(2026, 7, 6));

    readConfig(config, "outputfileDoodsonHarmonics", fileNameOut, Config::MUSTSET, "",  "");
    readConfig(config, "inputfile",                  fileNameIn,  Config::MUSTSET, "",  "");
    readConfig(config, "headerLines",                countHeader, Config::MUSTSET, "4", "skip number of header lines");
    readConfig(config, "minDegree",                  minDegree,   Config::DEFAULT, "0", "");
    readConfig(config, "maxDegree",                  maxDegree,   Config::MUSTSET, "",  "");
    readConfig(config, "GM",                         GM,          Config::DEFAULT, STRING_DEFAULT_GM, "Geocentric gravitational constant");
    readConfig(config, "R",                          R,           Config::DEFAULT, STRING_DEFAULT_R,  "reference radius");
    if(isCreateSchema(config)) return;

    // ==============================

    logStatus<<"read file from <"<<fileNameIn<<">"<<Log::endl;
    InFile file(fileNameIn);

    // skip headerlines
    std::string line;
    for(UInt i=0; i<countHeader; i++)
      std::getline(file, line);

    std::vector<Doodson> doodson;
    std::vector<Matrix>  cnmCos, cnmSin, snmCos, snmSin;
    for(;;)
    {
      std::getline(file, line);
      if(file.eof())
        break;
      if(line.empty())
        continue;

      // Coefficients to compute variations in normalized Stokes coefficients (unit = 10^-11)
      // Ocean tide model: FES2004 normalized model (fev. 2004) up to (100,100)
      // (long period from FES2002 up to (50,50) + equilibrium Om1/Om2, atmospheric tide NOT included)
      // Doodson Darw  l   m    DelC+     DelS+       DelC-     DelS-
      // 255.555 M2    2   1 -12.07164  -4.38919    -3.09008   1.50139
      std::stringstream ss(line);
      std::string doodstring, name;
      UInt   n, m;
      Double cPlus, sPlus, cMinus, sMinus;
      // the geopotential harmonic amplitudes for each tide constituent
      ss>>doodstring>>name>>n>>m>>cPlus>>sPlus>>cMinus>>sMinus;
      if(doodstring.size() == 6)
        doodstring = '0'+doodstring;

      // Constituent already exist?
      UInt idx = std::distance(doodson.begin(), std::find(doodson.begin(), doodson.end(), Doodson(doodstring)));
      if(idx >= doodson.size())
      {
        doodson.push_back(Doodson(doodstring));
        logStatus<<doodstring<<" "<<name<<Log::endl;
        cnmCos.push_back(Matrix(maxDegree+1, Matrix::TRIANGULAR, Matrix::LOWER));
        cnmSin.push_back(Matrix(maxDegree+1, Matrix::TRIANGULAR, Matrix::LOWER));
        snmCos.push_back(Matrix(maxDegree+1, Matrix::TRIANGULAR, Matrix::LOWER));
        snmSin.push_back(Matrix(maxDegree+1, Matrix::TRIANGULAR, Matrix::LOWER));
      }

      // Refer to IERS conventions 2010, eq. (6.15)
      // Expand the right-hand side of eq. (6.15) by Euler's formula and group the terms
      // into real (Delta C) and imaginary (Delta S) parts.
      if((n>=minDegree) && (n<=maxDegree))
      {
        // coefficients for cosine terms of the real part (Delta C)
        cnmCos.at(idx)(n,m) =  (cPlus + cMinus) * 1e-11;
        // coefficients for cosine terms of the imaginary part (Delta S)
        snmCos.at(idx)(n,m) =  (sPlus - sMinus) * 1e-11;
        // coefficients for sine terms of the real part (Delta C)
        cnmSin.at(idx)(n,m) =  (sPlus + sMinus) * 1e-11;
        // coefficients for sine terms of the imaginary part (Delta S)
        snmSin.at(idx)(n,m) = -(cPlus - cMinus) * 1e-11;
      }
    }

    // write results
    // -------------
    logStatus<<"write tides to <"<<fileNameOut<<">"<<Log::endl;
    writeFileDoodsonHarmonic(fileNameOut, DoodsonHarmonic(GM, R, doodson, cnmCos, snmCos, cnmSin, snmSin));
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
