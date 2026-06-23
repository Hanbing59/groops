/***********************************************/
/**
* @file preprocessingVariationalEquationOrbitFit.cpp
*
* @brief Fit variational equations to orbit observations.
*
* @author Torsten Mayer-Guerr
* @date 2012-05-30
*
*/
/***********************************************/

// Latex documentation
#define DOCSTRING docstring
static const char *docstring = R"(
This program fits an \configFile{inputfileVariational}{variationalEquation} to an observed \configFile{inputfileOrbit}{instrument} by estimating parameters
in a least squares adjustment. Additional to the initial satellite state for each arc, these parameters can be
\configClass{parametrizationGravity}{parametrizationGravityType}, satellite \configClass{parametrizationAcceleration}{parametrizationAccelerationType}
and stochastic pulses (velocity jumps) at given times, \configClass{stochasticPulse}{timeSeriesType}. The estimated parameters can be stored with
\configFile{outputfileSolution}{matrix} and an extra file with the parameter names is created. The fitted orbit is written
as new reference in \configFile{outputfileVariational}{variationalEquation} and additionally in \configFile{outputfileOrbit}{instrument}.

The observed orbit positions (\configFile{inputfileOrbit}{instrument}) together with the epoch-wise covariance matrix
(\configFile{inputfileCovariancePodEpoch}{instrument}) must be split in the same arcs as the variational equations but not
necessarily uniformly distributed (use irregularData in \program{InstrumentSynchronize}). An iterative downweighting of
outliers is performed by M-Huber method.

The observation equations (parameter sensitivity matrix) are computed by integration of the variational equations
(\configFile{inputfileVariational}{variationalEquation}) using a polynomial with \config{integrationDegree} and interpolated to the
observation epochs using a polynomial with \config{interpolationDegree}.

All parameters used here must be re-estimated in the full least squares adjustment
for the gravity field determination to get a solution which is not biased towards the reference field.
The solutions of additional estimations are relative (deltas) as the parameters are already used as Taylor point
in the reference orbit.

See also \program{PreprocessingVariationalEquation}.
)";

/***********************************************/

#include "programs/program.h"
#include "base/polynomial.h"
#include "files/fileMatrix.h"
#include "files/fileInstrument.h"
#include "files/fileVariationalEquation.h"
#include "files/fileParameterName.h"
#include "classes/ephemerides/ephemerides.h"
#include "classes/parametrizationGravity/parametrizationGravity.h"
#include "classes/parametrizationAcceleration/parametrizationAcceleration.h"
#include "classes/timeSeries/timeSeries.h"
#include "misc/varianceComponentEstimation.h"
#include "misc/observation/variationalEquationFromFile.h"

/***** CLASS ***********************************/

/** @brief Fit variational equations to orbit observations
* @ingroup programsGroup */
class PreprocessingVariationalEquationOrbitFit
{
public:
  VariationalEquationFromFile    variationalEquationFromFile;
  InstrumentFile                 podFile;
  InstrumentFile                 covPodEpochFile;
  EphemeridesPtr                 ephemerides;
  ParametrizationAccelerationPtr parameterAcceleration;
  ParametrizationGravityPtr      parameterGravity;
  UInt                           interpolationDegree;
  /// Number of arcs
  UInt                           arcCount;

  /// indexes of all outliers detected in each arc
  std::vector<std::vector<UInt>> outlierIndexes;
  struct ArcResult
  {
    /// post-fit residuals
    OrbitArc          residualArc;
    /// indexes of outliers detected in this arc
    std::vector<UInt> outlierIndexesArc;

    void save(OutArchive &oa) const
    {
      oa << nameValue("residualArc", residualArc);
      oa << nameValue("outlierIndexesArc", outlierIndexesArc);
    }

    void load(InArchive &ia)
    {
      ia >> nameValue("residualArc", residualArc);
      ia >> nameValue("outlierIndexesArc", outlierIndexesArc);
    }
  };

  // normal equations
  // ----------------
  /// \f$ A^TPA \f$, Normal matrix
  Matrix N;
  /// \f$ A^TPl \f$, right hand side
  Vector n;
  /// \f$ l^TPl \f$, weighted norm of the observations
  Double lPl;
  /// number of observations
  UInt   obsCount;
  /// number of outliers detected in all iterations
  UInt   outlierCount;
  /// number of outliers detected in the current iteration
  UInt   outlierCountNew;
  /// the solution vector
  Vector x;
  /// the post-fit standard deviation of unit weight
  Double sigma0;
  /// the post-fit standard deviation of unit weight from the previous iteration
  Double sigma0Last;

  void run(Config &config, Parallel::CommunicatorPtr comm);

  /** @brief Build normal equations for a specific arc
   * @param arcNo arc number
   */
  ArcResult buildNormals(UInt arcNo);
};

GROOPS_REGISTER_PROGRAM(PreprocessingVariationalEquationOrbitFit, PARALLEL, "fit variational equations to orbit observations", Preprocessing, VariationalEquation)

/***********************************************/
/***********************************************/

void PreprocessingVariationalEquationOrbitFit::run(Config &config, Parallel::CommunicatorPtr comm)
{
  try
  {
    FileName fileNameOutVariational, fileNameOutOrbit, fileNameOutSolution, fileNameOutResiduals;
    FileName fileNameInVariational;
    FileName fileNameInOrbit, fileNameInOrbitCov;
    UInt              integrationDegree;
    UInt              iterCount;
    std::vector<Time> stochasticPulse;
    TimeSeriesPtr     stochasticPulsePtr;

    renameDeprecatedConfig(config, "representation", "parametrizationGravity",      date2time(2020, 6, 3));
    renameDeprecatedConfig(config, "parameter",      "parametrizationAcceleration", date2time(2020, 6, 3));

    readConfig(config, "outputfileVariational",       fileNameOutVariational, Config::MUSTSET,  "",    "approximate position and integrated state matrix");
    readConfig(config, "outputfileOrbit",             fileNameOutOrbit,       Config::OPTIONAL, "",    "integrated orbit");
    readConfig(config, "outputfileSolution",          fileNameOutSolution,    Config::OPTIONAL, "",    "estimated calibration and state parameters");
    readConfig(config, "outputfileResiduals",         fileNameOutResiduals,   Config::OPTIONAL, "",    "post-fit residuals of orbit fit");
    readConfig(config, "inputfileVariational",        fileNameInVariational,  Config::MUSTSET,  "",    "approximate position and integrated state matrix");
    readConfig(config, "inputfileOrbit",              fileNameInOrbit,        Config::MUSTSET,  "",    "kinematic positions of satellite as observations");
    readConfig(config, "inputfileCovariancePodEpoch", fileNameInOrbitCov,     Config::OPTIONAL, "",    "3x3 epoch wise covariances");
    readConfig(config, "ephemerides",                 ephemerides,            Config::OPTIONAL, "jpl", "may be needed by parametrizationAcceleration");
    readConfig(config, "parametrizationGravity",      parameterGravity,       Config::DEFAULT,  "",    "gravity field parametrization");
    readConfig(config, "parametrizationAcceleration", parameterAcceleration,  Config::DEFAULT,  "",    "orbit force parameters");
    readConfig(config, "stochasticPulse",             stochasticPulsePtr,     Config::DEFAULT,  "",    "");
    readConfig(config, "integrationDegree",           integrationDegree,      Config::DEFAULT,  "7",   "integration of forces by polynomial approximation of degree n");
    readConfig(config, "interpolationDegree",         interpolationDegree,    Config::DEFAULT,  "7",   "orbit interpolation by polynomial approximation of degree n");
    readConfig(config, "iterationCount",              iterCount,              Config::DEFAULT,  "10",  "maximum number of iterations for outlier downweighting");
    if(isCreateSchema(config)) return;

    if(integrationDegree%2 == 0)
      throw(Exception("polnomial degree for integration must be odd."));

    if(stochasticPulsePtr)
      stochasticPulse = stochasticPulsePtr->times();

    // init
    // ----
    podFile.open(fileNameInOrbit);
    covPodEpochFile.open(fileNameInOrbitCov);
    InstrumentFile::checkArcCount({podFile, covPodEpochFile});
    variationalEquationFromFile.open(fileNameInVariational, parameterGravity, parameterAcceleration, stochasticPulse, ephemerides, integrationDegree);

    // =============================================

    x = Vector(variationalEquationFromFile.parameterCount());
    sigma0 = 1;
    outlierIndexes.resize(podFile.arcCount());
    std::vector<ArcResult> arcResults(podFile.arcCount());
    for(UInt iter=0; iter<iterCount; iter++)
    {
      // build normals
      // -------------
      logStatus<<"  "<<iter+1<<"-th iteration of "<<iterCount<<" (maximum), accumulate normal equations"<<Log::endl;

      N            = Matrix(variationalEquationFromFile.parameterCount(), Matrix::SYMMETRIC);
      n            = Vector(variationalEquationFromFile.parameterCount());
      lPl          = 0;
      obsCount     = 0;
      outlierCount = 0;
      outlierCountNew = 0;

      Parallel::forEach(arcResults, [this](UInt arcNo) {return buildNormals(arcNo);}, comm);

      Parallel::reduceSum(N,            0, comm);
      Parallel::reduceSum(n,            0, comm);
      Parallel::reduceSum(lPl,          0, comm);
      Parallel::reduceSum(obsCount,     0, comm);
      Parallel::reduceSum(outlierCount, 0, comm);
      Parallel::reduceSum(outlierCountNew, 0, comm);
      Parallel::broadCast(outlierCountNew, 0, comm);

      if(Parallel::isMaster(comm))
      {
        for(UInt arcNo=0; arcNo<podFile.arcCount(); arcNo++)
        {
          /// gather the indexes of detected outliers from each arc
          outlierIndexes.at(arcNo) = arcResults[arcNo].outlierIndexesArc;
          /// @todo Add statistics, like residuals RMS for each arc
        }
      }
      Parallel::broadCast(outlierIndexes, 0, comm);

      // Estimate parameters
      // -------------------
      if(Parallel::isMaster(comm))
      {
        // Regularize unused parameters
        for(UInt i=0; i<N.rows(); i++)
          if(N(i,i) == 0)
            N(i,i) = 1.0;

        x = solve(N,n);
        sigma0Last = sigma0;
        sigma0 = Vce::standardDeviation(lPl-inner(n,x), obsCount-x.rows(), 2.5/*huber*/, 1./*huberPower*/);
        logInfo<<"  "<<variationalEquationFromFile.parameterCount()%"%8i"s<<" par, "
               <<"last sigma0="<<sigma0Last%"%12.8f"s<<", "<<"this sigma0="<<sigma0%"%12.8f"s<<", "
               <<"delta sigma0="<<(sigma0Last-sigma0)%"%12.8f"s<<", "
               <<outlierCount%"%8i"s<<" ("<<outlierCountNew%"%8i"s<<" new) outliers among "
               <<obsCount%"%10i"s<<" obs ("<<(100.*outlierCount/obsCount)%"%6.2f"s<<"%)"<<Log::endl;
      }
      Parallel::broadCast(x, 0, comm);
      Parallel::broadCast(sigma0, 0, comm);
      Parallel::broadCast(sigma0Last, 0, comm);
      if((iter>0) && (outlierCount==0))
        break;
    }

    if(Parallel::isMaster(comm) && !fileNameOutSolution.empty())
    {
      logStatus<<"  write solution to file <"<<fileNameOutSolution<<">"<<Log::endl;
      writeFileMatrix(fileNameOutSolution, x);

      std::vector<ParameterName> parameterName;
      variationalEquationFromFile.parameterName(parameterName);
      writeFileParameterName(fileNameOutSolution.replaceFullExtension(".parameterName.txt"), parameterName);
    }

    if(Parallel::isMaster(comm) && !fileNameOutResiduals.empty())
    {
      logStatus<<"  write residuals to file <"<<fileNameOutResiduals<<">"<<Log::endl;
      std::list<OrbitArc> arcList;
      for(const ArcResult &arcResult : arcResults)
        arcList.push_back(arcResult.residualArc);
      InstrumentFile::write(fileNameOutResiduals, arcList);
    }

    // =============================================================================

    // Improve orbits with estimated parameters
    // ----------------------------------------
    std::vector<VariationalEquationArc> arcs(variationalEquationFromFile.arcCount());
    Parallel::forEach(arcs, [this](UInt arcNo) {return variationalEquationFromFile.refineVariationalEquationArc(arcNo, x);}, comm);

    // =============================================================================

    if(Parallel::isMaster(comm) && !fileNameOutVariational.empty())
    {
      logStatus<<"  write variational equation to file <"<<fileNameOutVariational<<">"<<Log::endl;
      writeFileVariationalEquation(fileNameOutVariational, variationalEquationFromFile.satellite(), arcs);
    }

    // =============================================

    if(Parallel::isMaster(comm) && !fileNameOutOrbit.empty())
    {
      logStatus<<"  write orbit to file <"<<fileNameOutOrbit<<">"<<Log::endl;
      std::list<Arc> arcList;
      for(UInt arcNo=0; arcNo<arcs.size(); arcNo++)
        arcList.push_back( arcs.at(arcNo).orbitArc() );
      InstrumentFile::write(fileNameOutOrbit, arcList);
    }

    // =============================================
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/

PreprocessingVariationalEquationOrbitFit::ArcResult PreprocessingVariationalEquationOrbitFit::buildNormals(UInt arcNo)
{
  try
  {
    ArcResult arcResult;
    OrbitArc pod = podFile.readArc(arcNo);
    if(pod.size() == 0)
      return arcResult;

    Vector l(3*pod.size());
    for(UInt k=0; k<pod.size(); k++)
    {
      l(3*k+0) = pod.at(k).position.x();
      l(3*k+1) = pod.at(k).position.y();
      l(3*k+2) = pod.at(k).position.z();
    }

    std::vector<Time> timePod = pod.times();
    VariationalEquationFromFile::ObservationEquation eqn = variationalEquationFromFile.integrateArc(timePod.front(), timePod.back(), TRUE/*position*/, FALSE/*velocity*/);
    Polynomial polynomial(eqn.times, interpolationDegree);
    l -= polynomial.interpolate(timePod, eqn.pos0, 3); // reference orbit
    Matrix A = polynomial.interpolate(timePod, eqn.PosDesign, 3);

    // decorrelation
    Covariance3dArc covPod = covPodEpochFile.readArc(arcNo);
    Arc::checkSynchronized({pod, covPod});
    for(UInt i=0; i<covPod.size(); i++)
    {
      Matrix W = covPod.at(i).covariance.matrix();
      W.setType(Matrix::SYMMETRIC);
      cholesky(W);

      triangularSolve(1., W.trans(), l.row(3*i,3));
      triangularSolve(1., W.trans(), A.row(3*i,3));
    }

    // downweight outliers
    if(quadsum(x))
    {
      const Double huber = 2.5;
      Vector e = l;
      // compute residuals, e = l - A*x
      matMult(-1, A, x, e);
      for(UInt k=0; k<pod.size(); k++)
      {
        // store residuals for output
        OrbitEpoch epoch = pod.at(k);
        epoch.position = Vector3d(e(3*k), e(3*k+1), e(3*k+2));
        epoch.velocity = Vector3d();
        epoch.acceleration = Vector3d();
        arcResult.residualArc.push_back(epoch);

        Double s = sqrt(quadsum(e.row(3*k,3))/3);
        if(s>huber*sigma0)
        {
          l.row(3*k,3) *= huber*sigma0/s;
          A.row(3*k,3) *= huber*sigma0/s;
          outlierCount += 3;
          /// Check if this is a new outlier
          if(std::find(outlierIndexes.at(arcNo).begin(), outlierIndexes.at(arcNo).end(), k) == outlierIndexes.at(arcNo).end())
          {
            outlierIndexes.at(arcNo).push_back(k);
            outlierCountNew += 3;
          }
        }
      }
    }

    lPl += quadsum(l);
    obsCount += l.rows();
    rankKUpdate(1., A, N);
    matMult(1., A.trans(), l, n);

    arcResult.outlierIndexesArc = outlierIndexes.at(arcNo);
    return arcResult;
  }
  catch(std::exception &e)
  {
    GROOPS_RETHROW(e)
  }
}

/***********************************************/
