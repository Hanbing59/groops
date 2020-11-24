/***********************************************/
/**
* @file observationMiscPodVariational.h
*
* @brief Precise Orbit data (variational equations).
*
* @author Torsten Mayer-Guerr
* @date 2015-06-02
*
*/
/***********************************************/

#ifndef __GROOPS_OBSERVATIONMSICPODVARIATIONAL__
#define __GROOPS_OBSERVATIONMSICPODVARIATIONAL__

/***********************************************/

#include "base/polynomial.h"
#include "files/fileInstrument.h"
#include "classes/ephemerides/ephemerides.h"
#include "classes/parametrizationGravity/parametrizationGravity.h"
#include "classes/parametrizationAcceleration/parametrizationAcceleration.h"
#include "misc/observation/variationalEquationFromFile.h"
#include "misc/observation/observationMiscPod.h"

/***** TYPES ***********************************/

class ObservationMiscPodVariational;
typedef std::shared_ptr<ObservationMiscPodVariational> ObservationMiscPodVariationalPtr;

/***** CLASS ***********************************/

/** @brief Precise Orbit data (variational equations).
* @ingroup miscGroup
* @see Observation */
class ObservationMiscPodVariational : public ObservationMiscPod
{
  InstrumentFile                 podFile;
  VariationalEquationFromFile    variationalEquation;
  Polynomial                     polynomial;
  UInt                           countArc;
  EphemeridesPtr                 ephemerides;
  ParametrizationGravityPtr      parameterGravity;
  ParametrizationAccelerationPtr parameterAcceleration;
  Bool                           accelerateComputation;

  // Indicies for design matrix A
  UInt countAParameter;
  UInt idxGravity, gravityCount;
  UInt idxState;

public:
  ObservationMiscPodVariational(Config &config);
 ~ObservationMiscPodVariational() {}

  Bool setInterval(const Time &timeStart, const Time &timeEnd);
  UInt parameterCount()          const {return countAParameter;}
  UInt gravityParameterCount()   const {return gravityCount;}
  UInt rightSideCount()          const {return 1;}
  UInt arcCount()                const {return countArc;}
  void parameterName(std::vector<ParameterName> &name) const;

  Arc computeArc(UInt arcNo, CovariancePodPtr cov=CovariancePodPtr(nullptr));
};

/***********************************************/

#endif /* __GROOPS__ */
