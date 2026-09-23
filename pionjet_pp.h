#ifndef __APPROXIMATE_PIONJET_PP_H__
#define __APPROXIMATE_PIONJET_PP_H__

#include <common_sidis.h>

// cuba library
#include "cuba.h"

extern bool NLO_JET; // NLO CALCULATION OR NOT...

//typedef int (*integrand_t)(const int *ndim, const double x[],
//  const int *ncomp, double f[], const double userdata[]);

double calc_pionjet(double zc, double pt, double jt, double y_min, double y_max,  double jt_min, double jt_max);

double calc_pionjet1(void);


double calc_pionjet_NLO(double zc, double pt, double jt, double y_min, double y_max,  double jt_min, double jt_max);


extern int NDIM_JET;
extern int NCOMP_JET;


double Pqq( double x);

double Pgq( double x);

double deltaPqq( double x);



#endif //#ifndef  __APPROXIMATE_PIONJET_PP_H__ 
