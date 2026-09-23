#ifndef __CONSTANTS_H
#define __CONSTANTS_H

#include <cmath>

const double PI                = 3.1415926543; // Pi
const double mpr               = 0.93827203 ; // the proton mass
//const double mpr                  = 0.938 ; // TEST NOTE STEFANO the proton mass
//const double mdeutron          = mpr ; // the proton mass
const double mdeutron          = 1.875613/2.   ; // the deutron mass/2. !!!!!!!!!!!
const double mneutron          = 0.93956536;    // the neutron mass
const double mpion             = 0.13956995;	// Mass of pion
const double mkaon             = 0.493667;	// Mass of kaon+-
const double mmuon             = 0.105658;	// Mass of the muon
const double e                 = exp( 1. );
const double eps               = 1.e-3      ; // Epsilon
extern  double kt2_average;    // <kt^2>
const  double alpha_em         =  1./137.035999679; // aplha_em0
const  double eu               =  2./3.; // Up quark charge
const  double ed               = -1./3.; // Down quark charge
const  double es               = -1./3.; // Strange quark charge

const double epsilon = 1.e-6; // Our precision constant...

const double convert_gev_barn = 0.389379292e-27;

 
// Fragmentation
// ffset   1,2,3 means K, KKP, BFGW
// fforder 0,1 is LO, NLO(MSbar)
// ihadron 1,2,3,4,5 is pi,K,h,p,n
// icp     1,2,3 chooses between particle, anti-particle or sum of both 
// ipi     1,2,3 is a flag for BFGW; inactive for K, KKP  
// icharge 0,1,2,3 is 0,+,-,+&-
const int pion       = 1;
const int hadron     = 2;
const int kaon       = 3;
const int proton     = 4;
const int deutron    = 5;
const int neutron    = 6;
const int antiproton = 7;



const double positive    = +1.;
const double negative    = -1.;
const double neutral     =  0.;

const int helicity_plus = +1;
const int helicity_minus = -1;

inline double pow2(double a) {
  return a*a;
};

inline double mod3(double x, double y, double z = 0) {
  return sqrt( pow2(x) + pow2(y) + pow2(z) );
};


// CKM matrix elements
const double Vud = 0.97427;
const double Vus = 0.22534;
const double Vcd = 0.22520;
const double Vcs = 0.97344;


// masses W and Z
const double MZ = 91.187621;
const double MW = 80.39923; // PDG

// Branching ratios
const double Zepem = 3.363e-2; // Z -> e+ e-
 
//Weinberg angle at MZ
const double sin2theta = 0.2311613;

//Fermi constant GF
const double GF = 1.166371e-5;

// vector and axial couplings for up and down quarks...
 const double Vu =  1./2. - 4./3. * sin2theta;
 const double Vd = -1./2. + 2./3. * sin2theta;
 const double Au =  1./2.;
 const double Ad = -1./2.;
 


//All constants and common functions between DY and SIDIS will be here:

const double EulerGamma = 0.5772156649015329 ; 
const double fourpi = 4.0 * PI ;
const double zeta2 = pow(PI,2)/6.0 ; 
const double CF = 4.0/3.0;
const double TF = 1.0/2.0;
const double C1 = 2.*exp(-EulerGamma);
//const double nf = 5; // MAIN
//FENG:
// const double nf = 4; //FENG
// const double LamQCD = 0.192; //FENG


// ROGERS:
// const double nf = 3; //TED&LEO
// const double LamQCD = 0.2123; //TED&LEO
// 
// const double beta_0 = 11.0 - 2.0*nf/3.0;
// const double beta_1 = 102. - nf*38.0/3.0;

//FENG AND PENG
const double nf = 5; //2014
const double LamQCD = 0.225; //2014

const double beta_0 = 11.0 - 2.0*nf/3.0;
const double beta_1 = 102. - nf*38.0/3.0;






//const double Q0 = sqrt(2.4); // this is the inital scale!
//const double Q0 = sqrt(2.4); // MAIN this is the inital scale!
//const double QM = sqrt(2.4); // this is the matching scale!
//const double Q0 = 3.2; //  this is the inital scale BLNY FITS!
//const double Q0 = 5.; //  this is test initial scale!
//const double Q0 = 1.; //  this is the inital scale BLNY FITS!
//const double Q0 = sqrt(1.69); //  this is the inital scale lowest of CTEQ6L
const double Q0 = sqrt(2.4); //  this is the inital scale lowest of MSTW...
//const double Q0 = 1.6*2.; //  Stefano

 


const double g2 = 0.160; 

  
// BMAX parameter:
//const double bmax = 1.;
//const double bmax = 0.5;
const double bmax = 1.5;








#endif
