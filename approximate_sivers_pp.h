#ifndef __APPROXIMATE_SIVERS_PP_H__
#define __APPROXIMATE_SIVERS_PP_H__


double jet_numerator_denominator(double xF); 

double jet_numerator_denominator_gauss(double xF);

double photon_numerator_denominator_gauss(double xF); // direct photon only

 

// C-------------------------------------------------------------------
// C     spin-dependent partonic cross section DIRECT PHOTON
// C-------------------------------------------------------------------
      double SivSIG_PHOTON(double s,double t,double u,
		double HU1,double HD1,double HUB1,double HDB1,double HS1,double HSB1,double HGL1,
		double U2,double D2,double UB2,double DB2,double S2,double SB2,double GL2);


    int Sivpp_PHOTON(double s,double t,double u,double* WQ1,double* WQ1b,double* WQ2,double* WQ3,
	  double* WQ4);
 

// C-------------------------------------------------------------------
// C     unpolarized partonic cross section PHOTON
// C-------------------------------------------------------------------

      double UUSIG_PHOTON(double s,double t,double u,
		double U1,double D1,double UB1,double DB1,double S1,double SB1,double GL1,
		double U2,double D2,double UB2,double DB2,double S2,double SB2,double GL2);

int Upp_PHOTON(double s,double t,double u, double* WQ1,double* WQ2,double* WQ3,
		double* WQ4);

#endif //#ifndef  __APPROXIMATE_SIVERS_PP_H__ 
