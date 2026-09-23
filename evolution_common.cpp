//========================================================= includes
#include <iostream>
#include <iomanip>
#include <fstream>
#include <math.h>

#include <evolution_common.h>
#include <constants.h>


// Hoppet
#ifndef __HOPPET_V1__
#include "hoppet_v1.h"
#endif

using namespace std;



double AlphaS(double mu){

  double nf = 5; //2014
  double LamQCD = 0.225; //2014

  double beta_0 = 11.0 - 2.0*nf/3.0;
  double beta_1 = 102. - nf*38.0/3.0;



// Uncomment below for NLO expression

/*
return fourpi/beta_0/log(mu*mu/LamQCD/LamQCD) * ( 1 - beta_1/beta_0/beta_0*log(log(mu*mu/LamQCD/LamQCD))/log(mu*mu/LamQCD/LamQCD));
*/

//return fourpi/beta_0/log(mu*mu/LamQCD/LamQCD); // MAIN

return fourpi/beta_0/log(mu*mu/LamQCD/LamQCD) * ( 1. - beta_1/beta_0/beta_0*log(log(mu*mu/LamQCD/LamQCD))/log(mu*mu/LamQCD/LamQCD)); // PENG 2014

}


double gammaK(double mu){

return AlphaS(mu)/PI*2.0*CF;
}



double gammafatone(double mu){
return AlphaS(mu)/PI*CF*(3.0/2.0);
}


double Ktilde(double mu, double bT){
return -AlphaS(mu)/PI*CF*(2.0*log(mu*bT)-log(4.0)+2.0*EulerGamma);
}


//////////
//////////
// PROJECT WITH PENG AND FENG 2014 BEGIN
//////////
//////////

// from PENG AND FENG 2014: integral of Eq.(4)
double Sudakovkernel( double mu2 , double Q2) // It is integral of Spert from Eq.(4) Spert = int_mub^Q d mu/mu (A ln(Q/mu) + B) at NLL
{
double Bb0=(11.0-(2.0/3.)*nf);
double Bb1=102.0-38.0*nf/3;
double CA=3.0;
double A1=CF;
double A2= CF*((67.0/36.0-PI*PI/12.0)*CA-5.0/18.0*nf);
double B1=-(3.0/2.0)*CF;
 // double B2=CF*CF*(pi*pi/4.0-3.0/16.0-3*re3)+CF*Nc*(11.0/36.0*pi*pi-193.0/48.0+3.0/2.0*re3)+CF*nf*(-1.0/18.0*pi*pi+17.0/24.0);
double B2=0;
double LQCD2 = pow2(LamQCD);

return (2*(54*A1*pow(Bb0,5)*pow(log(mu2/LQCD2),4)*(-1. + log(log(mu2/LQCD2))) -
       8.*pow(Bb1,2)*(B2 + A2*log(Q2/mu2))*(2. + 6.*log(log(mu2/LQCD2)) + 9.*pow(log(log(mu2/LQCD2)),2)) +
       27.*pow(Bb0,3)*pow(log(mu2/LQCD2),3)*(-8.*A2*Bb0 + 2.*A1*Bb1 + 2.*(-4.*A2*Bb0 + B1*pow(Bb0,2) + A1*Bb1
 + A1*pow(Bb0,2)*log(Q2/mu2))*log(log(mu2/LQCD2)) +
          A1*Bb1*pow(log(log(mu2/LQCD2)),2)) + 2.*Bb1*log(mu2/LQCD2)*
        (54.*B2*pow(Bb0,2) + 19.*A2*Bb1 + 6.*(18.*B2*pow(Bb0,2) + 5.*A2*Bb1)*log(log(mu2/LQCD2)) + 18.*A2*Bb1*pow(log(log(mu2/LQCD2)),2) +
          54.*A2*pow(Bb0,2)*log(Q2/mu2)*(1. + 2.*log(log(mu2/LQCD2)))) -
       54.*pow(Bb0,2)*pow(log(mu2/LQCD2),2)*(4.*B2*pow(Bb0,2) + 6.*A2*Bb1 - B1*Bb0*Bb1 + (4.*A2 - B1*Bb0)*Bb1*log(log(mu2/LQCD2)) +
						Bb0*log(Q2/mu2)*(4.*A2*Bb0 - A1*Bb1
 - A1*Bb1*log(log(mu2/LQCD2))))))/(27.*pow(Bb0,6)*pow(log(mu2/LQCD2),3));

};


//FENG PENG 2014
double PertEvolFeng(double mu, double bt){

double bstar = bt/sqrt(1.0+ pow(bt,2)/pow(bmax,2));
double mub2 = pow(C1/bstar,2);
double Q2 =  pow2(mu);


double Ssud =  Sudakovkernel( Q2 , Q2) - Sudakovkernel( mub2 , Q2);

// Test no evolution
//Ssud = 0.;

if (   -Ssud/2.  <= -20.)
	return  exp(-20.) ;
else
	return  exp(-Ssud/2.) ;
};





//////////
//////////
// PROJECT WITH PENG AND FENG 2014 END
//////////
//////////



double PertEvol(double mu, double bT, double xi,double Q0){



double C1 = 1.123;
double bstar = bT/sqrt(1.0+ pow(bT,2)/pow(bmax,2));
//double mub2 = pow(C1/bstar,2);
double mub = C1/bstar;
double mu0 = Q0;





double A1 = (2*PI*log(log(pow(mu,2)/pow(LamQCD,2))/
       log(pow(mu0,2)/pow(LamQCD,2))))/beta_0;




double B1 = 0.0;

// For NLO result uncomment below
// double B1 = (2*beta_1*PI*(-((1 + log(log(pow(mu,2)/pow(LamQCD,2))))/
//           log(pow(mu,2)/pow(LamQCD,2))) +
//        (1 + log(log(pow(mu0,2)/pow(LamQCD,2))))/
//         log(pow(mu0,2)/pow(LamQCD,2))))/pow(beta_0,3);


double A2 = (-2*PI*log(mu/mu0))/beta_0 + (2*PI*log(sqrt(xi)/LamQCD)*
      log(log(pow(mu,2)/pow(LamQCD,2))/
        log(pow(mu0,2)/pow(LamQCD,2))))/beta_0 ;

double B2 = 0.0 ;


// For NLO result uncomment below

// double B2 = -(beta_1*PI*(pow(log(log(pow(mu,2)/pow(LamQCD,2))),2) -
//          pow(log(log(pow(mu0,2)/pow(LamQCD,2))),2)))/
//     (2.*pow(beta_0,3)) + (2*beta_1*PI*log(sqrt(xi)/LamQCD)*
//       (-((1 + log(log(pow(mu,2)/pow(LamQCD,2))))/
//            log(pow(mu,2)/pow(LamQCD,2))) +
//         (1 + log(log(pow(mu0,2)/pow(LamQCD,2))))/
//          log(pow(mu0,2)/pow(LamQCD,2))))/pow(beta_0,3);


double A3 = (2*PI*log(log(pow(mub,2)/pow(LamQCD,2))/
       log(pow(mu0,2)/pow(LamQCD,2))))/beta_0;


double B3 = 0.0 ;

// For NLO result uncomment below
/*
double B3 = (2*beta_1*PI*(-((1 + log(log(pow(mub,2)/pow(LamQCD,2))))/
          log(pow(mub,2)/pow(LamQCD,2))) +
       (1 + log(log(pow(mu0,2)/pow(LamQCD,2))))/
        log(pow(mu0,2)/pow(LamQCD,2))))/pow(beta_0,3);
*/


//cout << Ktilde(mub,bstar) << endl;

// NO SUDAKOV ENHANCEMENT:
double result = ( log(sqrt(xi)/mu0)*Ktilde(mub,bstar) + 3.0*CF/2.0/PI *(A1-B1) - 2.0*CF/PI*(A2-B2) +
log(sqrt(xi)/mu0) *2.0*CF/PI * (A3-B3)*0 );
//if (mu0 >= mu) return 1.; // test with Alessandro
//if (mub >= mu) return 1.; // test with Alessandro
// if ( result  <= -20. )
// 	result = -20.;
//  else if (result >= 0. ) // NO SUDAKOV ENHANCEMENT
//    result = 0.;
//
// KANG AND ZHONGBO CLAIM THAT SUDAKOV ENHANCEMENT IS APPROPRIATE, SO LET'S TRY TO INCLUDE THAT TOO:
//if ( result  <= -20. )
//	result = -20.;

//result = 0; // NO EVOLUTION TEST WITH STEFANO

//return exp(result); // MAIN


// //FENG AND PENG 2014:
return exp(-0.5*( Sudakovkernel( mu*mu , mu*mu) - Sudakovkernel( mub*mub , mu*mu)) ); //0.5 because I have S/2 for F and D

}




// THIS IS MAIN FUNCTION FOR OUR PRL PAPER!!!
double NonPertEvol(double bt, double xi, double g2, double Q0) {

// cout << " res = " << exp(-g2*bt*bt*log(sqrt(xi)/(Q0*QM))/2.0) << endl;
//NonPertEvNonPertEvol(bt,xi,g2,Q0)
// cout << " sqrt(xi) = " << sqrt(xi) << " Q0*QM = " << Q0*QM << endl;
// cout << " res = " << log(sqrt(xi)/(Q0*QM)) << endl;
 double result;
 if((result = -g2*log(sqrt(xi)/(Q0))*bt*bt/2.0) <= -20.)
  result = -20.;

// return exp(result);
return exp(result);

}
