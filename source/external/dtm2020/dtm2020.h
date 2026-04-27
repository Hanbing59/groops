/***********************************************/
/**
* @file dtm2020.h
*
* @brief Fortran Wrapper.
*
* @author Andreas Kvas
* @date   2022-04-19
*
*/
/***********************************************/

#ifndef __GROOPS_DTM2020__
#define __GROOPS_DTM2020__

#include "external/fortran.h"

extern "C"
{
  void dtm2020initWrapper(const char *fileName);
  void dtm2020calcWrapper(const F77Float &day, const F77Float &f, const F77Float &fbar, const F77Float &akp, const F77Float &akpMean, const F77Float &alti, const F77Float &hl, const F77Float &alat, const F77Float &xlon, F77Float &tz, F77Float &tinf, F77Float &ro, F77Float d[6], F77Float &wmm);
}

/***********************************************/

#endif
