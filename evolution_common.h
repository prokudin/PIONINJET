#ifndef  __EVOLUTION_COMMON
#define __EVOLUTION_COMMON

// some useful functions for TMD evolution
 
double AlphaS(double mu);

double gammaK(double mu);

double gammafatone(double mu);

double Ktilde(double mu, double bT);

double PertEvolFeng(double mu, double bT);

 
double PertEvol(double mu, double bT, double xi,double Q0);

double NonPertEvol(double bt, double xi, double g2, double Q0);
 
double Sudakovkernel( double mu2 , double Q2);

void set_hoppet( void ); // this function will set hoppet correctly 


#endif // #ifndef __EVOLUTION_COMMOSN