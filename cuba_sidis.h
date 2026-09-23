/* /home/prokudin/cpp/collins08 my headers for fitaut */
/* Alexei Prokudin 31/10/2008        */

#ifndef __CUBA_SIVERS_SIDIS_H__
#define __CUBA_SIVERS_SIDIS_H__

#include <common_sidis.h>

#include "cuba.h"

//typedef void (*integrand_t)(const int *, const double [],
//  const int *, double [], void *userdata);

typedef int (*integrand_t)(const int *ndim, const double x[],
  const int *ncomp, double f[], void *userdata);


/* double Asymmetry_SIDIS(double x, integrand_t integrand); */
void Asymmetry_SIDIS(const int ncomputations, const int ndim, double x, integrand_t integrand, double asymmetry[]); // this depend on how many effects we want to sum for BM effect for example: BM + CAHN -> ncomputations = 2


// Integrand for Cuba Sivers, BM etc. with evolution!
int Integrand_EVOLV_SIDIS(const int *ndim, const double xxx[],
			     const int *ncomp, double ff[], void *userdata);


// Integrand for Cuba Sivers, BM etc.
int Integrand_SIDIS(const int *ndim, const double xxx[],
			     const int *ncomp, double ff[], void *userdata);

// Integrand for LP --> PI X
int Integrand_LP_SIDIS(const int *ndim, const double xxx[],
			     const int *ncomp, double ff[], void *userdata);

// Integrand for LP --> JET X
int Integrand_LPJET_SIDIS(const int *ndim, const double xxx[],
			     const int *ncomp, double ff[], void *userdata);

int Integrand_SIDIS_PTWEIGHTED_WW(const int *ndim, const double xxx[],
				const int *ncomp, double ff[], void *userdata); // Used for pt weighted integrations in WW approximation (one additional variable  is added!) ...

int Integrand_DIS(const int *ndim, const double xxx[],
		   const int *ncomp, double ff[], void *userdata); // Used for DIS, G2 etc...


int Integrand_SIDIS_WW(const int *ndim, const double xxx[],
			const int *ncomp, double ff[], void *userdata); // Used for pt UNweighted integrations in WW approximation (one additional variable  is added!) ...


//Integrand for PP --> H X twist-3 approach
int Integrand_PP(const int *ndim, const double xxx[],
                      const int *ncomp, double ff[], void *userdata);
//Integrand for PP --> JET X twist-3 approach
int Integrand_PP_JET(const int *ndim, const double xxx[],
                      const int *ncomp, double ff[], void *userdata);

// FOR UNPOLARISED CROSS SECTIONS... NO EVOLUTION
int Integrand_SIDIS_UNPOLARISED(const int *ndim, const double xxx[],
                      const int *ncomp, double ff[], void * userdata); 

// FOR UNPOLARISED CROSS SECTIONS...NORMALIZATION... NO EVOLUTION  
int Integrand_DIS_NORMALIZATION_UNPOLARISED(const int *ndim, const double xxx[],
                      const int *ncomp, double ff[], void * userdata); 



// FOR UNPOLARISED CROSS SECTIONS... WITH TMD EVOLUTION  
int Integrand_SIDIS_EVOLV_UNPOLARISED(const int *ndim, const double xxx[],
                      const int *ncomp, double ff[], void * userdata); 


#endif //#ifndef  __CUBA_SIVERS_SIDIS_H__
