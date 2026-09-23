//========================================================= includes
#include <iostream>
#include <iomanip>
#include <fstream>
#include <math.h>
#include <cfortran.h>
#include <stdio.h>
#include <sidis.h>
#include <cmath>

// Cuba library
#include <cuba_sidis.h>



using namespace std;
using namespace sidis_namespace;


int NDIM = 2;
int NCOMP = 2;
//#define NDIM 2
//#define NCOMP 2
#define USERDATA NULL
double EPSREL = 1.e-3; // MAIN
double EPSABS = 1.e-12;
/*double EPSREL = 1.e-3;
double EPSABS = 1.e-12;*/
//#define EPSREL 1e-2
//#define EPSABS 1e-16


int VERBOSE = 2;


//#define VERBOSE 2
 
 
 
 
#define LAST 4
#define SEED 0
#define NVEC 1

//#define MINEVAL 0
//#define MAXEVAL 50000
 
//#define NSTART 1000
//#define NINCREASE 500

 
#define NBATCH 100000
#define GRIDNO 0
#define STATEFILE NULL
#define SPIN NULL

#define NNEW 1000
#define FLATNESS 25.

#define KEY1 47
#define KEY2 1
#define KEY3 1
#define MAXPASS 5
#define BORDER 0.
#define MAXCHISQ 10.
#define MINDEVIATION .25
#define NGIVEN 0
#define LDXGIVEN NDIM
#define NEXTRA 0

#define KEY 0
 
// int MINEVAL = 50;  // TEST FITTING PP GLOBAL NOTE VERY ROUGHT
// int MAXEVAL = 500; // TEST FITTING PP GLOBAL
// int NSTART  = 50; // TEST FITTING PP GLOBAL
// int NINCREASE = 50; // TEST FITTING PP GLOBAL


// int MINEVAL = 100;  // TEST FITTING PP GLOBAL NOTE slightly better 
// int MAXEVAL = 2500; // TEST FITTING PP GLOBAL
// int NSTART  = 100; // TEST FITTING PP GLOBAL
// int NINCREASE = 500; // TEST FITTING PP GLOBAL


// int MINEVAL = 10000;  // MAIN FOR CUHRE
// int MAXEVAL = 15000; // TEST  
// int NSTART  = 10000; // TEST  
// int NINCREASE = 5000; // TEST  


//  int MINEVAL = 10000; // MAIN FOR VEGAS
//  int MAXEVAL = 1500000; // MAIN
//  int NSTART  = 20000; // MAIN
//  int NINCREASE = 20000; // MAIN

 
 int MINEVAL = 100000; // MAIN FOR VEGAS
 int MAXEVAL = 1500000; // MAIN
 int NSTART  = 100000; // MAIN
 int NINCREASE = 10000; // MAIN


// int MINEVAL = 500;  // MAIN FOR CUHRE
// int MAXEVAL = 15000; // TEST  
// int NSTART  = 500; // TEST  
// int NINCREASE = 500; // TEST  


void Asymmetry_SIDIS(const int ncomputations, const int ndim, double x, integrand_t integrand, double asymmetry[]){
  
  int comp, nregions, neval, fail;
  
  sidis& SIDIS = *sidis::Instance();

  NCOMP = ncomputations + 1; // Important!!! Check the subroutine  Integrand_SIDIS ncomputations defines how many effects we want to sum + 1 = also unpolarised
  NDIM  = ndim;

  if( cuba_verbose )
    { VERBOSE = 2;}
  else 
    { VERBOSE = 0;}



 

  switch(SIDIS.dependence){
  case xf_dependence:
    xf_used = x;
    if ( cuba_verbose ) cout << "*************************************" << endl;
    if ( cuba_verbose ) cout << "XF = " << xf_used << endl;
    break;
   case pt_dependence:
    pt_used = x;  
    if ( cuba_verbose ) cout << "*************************************" << endl;
    if ( cuba_verbose ) cout << "PT = " << pt_used << endl;
    break;
   default:
    printf("Unknown dependence in Asymmetry_sidis!\n");
    break;
  }

  double integral[NCOMP], error[NCOMP], prob[NCOMP];
  double numerator[NCOMP];
  double denominator = 0;

 
    cubacores(0,1000); // uses only ONE CORE FOR CALCULATIONS

    Vegas(NDIM, NCOMP, integrand, USERDATA, NVEC,
    EPSREL, EPSABS, VERBOSE, SEED,
    MINEVAL, MAXEVAL, NSTART, NINCREASE, NBATCH,
    GRIDNO, STATEFILE, SPIN,
    &neval, &fail, integral, error, prob);


//    Cuba(1, NDIM, NCOMP, integrand, integral, error, prob);
    

//   Suave(NDIM, NCOMP, integrand,
//     EPSREL, EPSABS, VERBOSE, MINEVAL, MAXEVAL,
// 	NNEW, FLATNESS, &nregions,
//     &neval, &fail, integral, error, prob);

//   Cuhre(NDIM, NCOMP, integrand, USERDATA, NVEC,
//     EPSREL, EPSABS, VERBOSE | LAST,
//     MINEVAL, MAXEVAL, KEY,
//     STATEFILE, SPIN,
//     &nregions, &neval, &fail, integral, error, prob);


//   Divonne(NDIM, NCOMP, integrand, USERDATA, NVEC,
//     EPSREL, EPSABS, VERBOSE, SEED,
//     MINEVAL, MAXEVAL, KEY1, KEY2, KEY3, MAXPASS,
//     BORDER, MAXCHISQ, MINDEVIATION,
//     NGIVEN, LDXGIVEN, NULL, NEXTRA, NULL,
//     STATEFILE, SPIN,
//     &nregions, &neval, &fail, integral, error, prob);



  if ( cuba_verbose ) printf("VEGAS RESULT:\tneval %d\tfail %d\n",
    neval, fail);

  if ( cuba_verbose ) for( int comp = 0; comp < NCOMP; ++comp ){
    cout << "VEGAS RESULT:" << integral[comp] << " +- " << error[comp] << " p = " << prob[comp] << endl;
  }  

  

  denominator = integral[0];
  if ( cuba_verbose ) cout << "UNP = " << denominator << " +- " << error[0] << endl;

  for( int comp = 0; comp < NCOMP - 1; ++comp ){ // polarised part
    numerator[comp] = integral[comp + 1];
    if ( cuba_verbose ) cout << "N: " << comp << " POL = " << numerator[comp] << " +- " << error[comp + 1] << endl;

  }


  for( int comp = 0; comp < NCOMP - 1; ++comp ){ // polarised part
    
    asymmetry[comp] = numerator[comp]/denominator; // MAIN
    
// UNPOLARISED// UNPOLARISED// UNPOLARISED// UNPOLARISED
    if(SIDIS.calculation == no_effect || SIDIS.calculation == normalization){ // UNPOLARISED
       for( int comp = 0; comp < NCOMP - 1; ++comp ){ // polarised part
         asymmetry[comp] = denominator; // UNPOLARISED
       if ( cuba_verbose ) cout << "UNPOLARISED: " << asymmetry[comp] << endl;
       };  
    }; // USE THIS TRICK FOR UNPOLARISED 


    if ( cuba_verbose ) cout << "N: " << comp << " ASYMMETRY = " << asymmetry[comp] << endl;
  }
  


}






  

//Integrand for PP --> H X twist-3 approach
int Integrand_PP(const int *ndim, const double xxx[],
                      const int *ncomp, double ff[], void *userdata)
{
  sidis& SIDIS = *sidis::Instance();

  double xa = 1., xb = 1., zc = 1., pt = 1., y = 1., xF = 1.; // y here stands for rapidity!
 
  double SS = 0, TT = 0, UU = 0, sh = 0, uh = 0, th = 0;

  double s = SIDIS.GetEnergy();
  double sqrts = sqrt(s); 

  double jacobian = 1.;

  double Q2 = 0;

  

 
  for( int i = 0; i < *ncomp; i++) ff[i] = 1.e-30; // Important!!!! otherwise the noise will count to the integral
 
  // xxx is n dim cube!!!
  switch(SIDIS.dependence){
  case pt_dependence:
    pt = pt_used;
    
    xF = SIDIS.Average.x; // this is average xF for STAR measurement
 

    xa = xxx[0];
    xb = xxx[1];
 

    Q2 = pt*pt;
 
 
    //calculate rapidity
    y = asinh( sqrts * xF/ (2. * pt) ); // pseudorapidity
     

    SS = s;
    TT = -sqrts * pt * exp(-y);
    UU = -sqrts * pt * exp(y);

    zc = -(TT/xb + UU/xa)/s;

    sh = xa * xb * s;
    uh = -pt/zc * xb * sqrts * exp(y);
    th = -pt/zc * xa * sqrts * exp(-y);


    break;
  case y_dependence: // rapidity TEST NOT WORKING !!!!!!!!!!!!!! TEST
    //pt = pt_used;
    //xF = xf_used;
    y = y_used;   

    xa = xxx[0];
    xb = xxx[1];
 

    xF = xa - xb;

    Q2 = pt*pt;

    pt = xF * sqrts/( exp(y) - exp(-y) ); 

    SS = s;
    TT = -sqrts * pt * exp(-y);
    UU = -sqrts * pt * exp(y);

    zc = -(TT/xb + UU/xa)/s;

    sh = xa * xb * s;
    uh = -pt/zc * xb * sqrts * exp(y);
    th = -pt/zc * xa * sqrts * exp(-y);


 
    break;
  case xf_dependence:
    xF = xf_used;


    xa = xxx[0];
    xb = xxx[1];
    
    y = SIDIS.Average.y; // this is rapidity!


    if(EXPERIMENT == RHIC_BRAHMS){
	y = - log( tan( (SIDIS.Average.y * PI/180.) / 2. ) ); // y is a \theta angle in grad this case...
    }

    pt = xF * sqrts/( exp(y) - exp(-y) ); 


    Q2 = pt*pt;
//    Q2 = 4.*pt*pt; // Q = 2 pt
//    Q2 = pt*pt/4.; // Q = pt/2.

    SS = s;
    TT = -sqrts * pt * exp(-y);
    UU = -sqrts * pt * exp(y);

    zc = -(TT/xb + UU/xa)/s;

    sh = xa * xb * s;
    uh = -pt/zc * xb * sqrts * exp(y);
    th = -pt/zc * xa * sqrts * exp(-y);

//cout << "here" << xa << " " << xb << " " << zc <<endl; 
    break;
  default:
    printf("Unknown dependence in Integrand_PP cuba_sidis.cpp !\n");
    break;
  }



  if(
//      SIDIS.x_cut(xa) &&
     (-UU/(SS+TT) < xa && xa  < 1.) &&
     (-TT/(SS+UU/xa) < xb && xb  < 1.) &&
     (0. < zc && zc  < 1.) &&
     SIDIS.Q2_cut(Q2) ){
     
     

    if ( cuba_verbose ) {
      switch(SIDIS.dependence){
      case kinematical_dependence:
      case pt_dependence:
      case pt2_dependence:
      case z_dependence:
      case x_dependence:
      case xf_dependence:
      case y_dependence:
// 	Q2_tree = Q2;
// 	xf_tree  = xF;
//  	x_tree  = xa;
// 	y_tree  = xb;
// 	z_tree  = zc;
// 	pt_tree = pt;
// 	tree->Fill();
	break;
      default:
	printf("Unknown dependence in Integrand_PP cuba_sivers_sidis.cpp !\n");
	break;
      }
    }




    if(SIDIS.calculation ==  collins_effect_pp){  //COLLINS
 
       //Collins
       SIDIS.effect = collins_effect; // NEED TO BE MORE FLEXIBLE
//       if( ( SIDIS.approximate_collins_pp_calculation( xa, xb, zc, y, pt, Q2) ) || 0){
// 	cerr << "cuba_sidis.cpp: ERROR wrong result!" << endl;     
//       }
       if( ( SIDIS.approximate_collins_pp_calculation_NLO( xa, xb, zc, y, pt, Q2) ) || 0){
	cerr << "cuba_sidis.cpp: ERROR wrong result!" << endl;     
      }

     
      
      ff[1] =  SIDIS.GetCrossSection() * jacobian;

      
      //Unpolarized
      SIDIS.effect = no_effect;
//       if( (  SIDIS.approximate_collins_pp_calculation( xa, xb, zc, y, pt, Q2 ) ) || 0){
// 	cerr << "cuba__sidis.cpp: ERROR wrong result!" << endl;     
//       }
      if( (  SIDIS.approximate_collins_pp_calculation_NLO( xa, xb, zc, y, pt, Q2 ) ) || 0){
	cerr << "cuba__sidis.cpp: ERROR wrong result!" << endl;     
      }

      ff[0] =  SIDIS.GetCrossSection() * jacobian;

    }




  }
return 0;
}













