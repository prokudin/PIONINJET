/* /home/prokudin/cpp/collins08 my headers for fitaut */
/* Alexei Prokudin 31/10/2008        */

#ifndef __TWIST3_H__
#define __TWIST3_H__


 
using namespace std;

double alphasQCD( double q);

int Upp(double s,double t,double u, double* WQ1,double* WQ2,double* WQ3,
		double* WQ4,double* WQ5,double* WQ6,double* WQ7,double* WQ8);

      double UUSIG(double s,double t,double u,
		double U1,double D1,double UB1,double DB1,double S1,double SB1,double GL1,
		double U2,double D2,double UB2,double DB2,double S2,double SB2,double GL2,
		double fu,double fub,double fd,double fdb,double fs,double fsb,double fg
);


      int Sivpp(double s,double t,double u,double* WQ1,double* WQ2,double* WQ3,
	  double* WQ4,double* WQ5,double* WQ6,double* WQ7,double* WQ8,double* WQ9,double* WQ10,
	  double* WQ11,double* WQ12,double* WQ13,double* WQ14,double* WQ15,double* WQ16,double*  WQ17,
	  double* WQ18,double* WQ19,double* WQ20,double* WQ21,double* WQ22,double* WQ23,double* WQ24);


      double SivSIG(double s,double t,double u,
		double HU1,double HD1,double HUB1,double HDB1,double HS1,double HSB1,double HGL1,
		double U2,double D2,double UB2,double DB2,double S2,double SB2,double GL2,
		double fu,double fub,double fd,double fdb,double fs,double fsb,double fg
);


int Colpp(double s,double t,double u,double* WQ1,double* WQ2,double* WQ3,
	  double* WQ4,double* WQ5,double* WQ6,double* WQ7,double* WQ8);


      double ColSIG(double s,double t,double u,
		double HU1,double HD1,double HUB1,double HDB1,double HS1,double HSB1,double HGL1,
		double U2,double D2,double UB2,double DB2,double S2,double SB2,double GL2,
		double fu,double fub,double fd,double fdb,double fs,double fsb,double fg
);

#endif //#ifndef   __TWIST3_H__
