// /home/prokudin/cpp/dy/ 
// Functions for DY.
// Alexei Prokudin 26/04/2011

//========================================================= includes
#include <iostream>
#include <iomanip>
#include <fstream>
#include <math.h>
#include <cfortran.h>

#include <sidis.h>
#include <twist3.h>
  

using namespace std;


// c--------------------------------------------------------------
// c     alphas(q) -- strong coupling constant
// c--------------------------------------------------------------      

double alphasQCD( double q){
      int nf;
      double lambda, alphasQCD;
      double mb=4.5;

      if(q <= mb){
         nf=4;
         lambda = 0.326;}
      else{
         nf=5;
         lambda=0.226;
      };

      double b0=11.-2./3. * nf;
      double b1=51.-19./3. * nf;
      double q2=q*q;

      double lambda2=lambda*lambda;
      double tt=log(q2/lambda2);

      alphasQCD=4. * PI/(b0*tt)*(1.-2.*b1/(b0*b0)*log(tt)/tt);

      return alphasQCD;
};


// C-------------------------------------------------------------------

int Upp(double s,double t,double u, double* WQ1,double* WQ2,double* WQ3,
		double* WQ4,double* WQ5,double* WQ6,double* WQ7,double* WQ8){
       
      double Nc=3.;

// c     qq' ->qq'
      *WQ1=(pow2(Nc)-1.)/(2.*pow2(Nc))*(s*s+u*u)/(t*t);
// c     qq  ->qq
      *WQ2=(pow2(Nc)-1.)/(2.*pow2(Nc))*( (s*s+u*u)/(t*t)+(s*s+t*t)/(u*u) )
          -(pow2(Nc)-1.)/pow(Nc,3)*(s*s)/(t*u);
// c     qqb ->q'qb'
      *WQ3=(pow2(Nc)-1.)/(2.*pow2(Nc))*(t*t+u*u)/(s*s); 
// c     qqb ->qqb
      *WQ4=(pow2(Nc)-1.)/(2.*pow2(Nc))*( (s*s+u*u)/(t*t)+(t*t+u*u)/(s*s) )
          -(pow2(Nc)-1.)/pow(Nc,3)*(u*u)/(s*t);
// c     qqb ->gg
      *WQ5=pow2(pow2(Nc)-1.)/(2.*pow(Nc,3))*(u/t+t/u)
          -(pow2(Nc)-1.)/Nc*(t*t+u*u)/(s*s); // 
// c     gg  ->qqb
      *WQ6=1./(2.*Nc)*(t/u+u/t)-Nc/(pow2(Nc)-1.)*(t*t+u*u)/(s*s);
// c     qg  ->qg
      *WQ7=(pow2(Nc)-1.)/(2.*pow2(Nc))*(-s/u-u/s)+(s*s+u*u)/(t*t);
// c     gg  ->gg
      *WQ8=4.*pow2(Nc)/(pow2(Nc)-1.)*(3.-t*u/(s*s)-s*u/(t*t)-s*t/(u*u));


      return 0;
};




// C-------------------------------------------------------------------
// C     unpolarized partonic cross section
// C-------------------------------------------------------------------

      double UUSIG(double s,double t,double u,
		double U1,double D1,double UB1,double DB1,double S1,double SB1,double GL1,
		double U2,double D2,double UB2,double DB2,double S2,double SB2,double GL2,
		double fu,double fub,double fd,double fdb,double fs,double fsb,double fg
){
       
      //real*8 HU1,HUB1,HD1,HDB1,HS1,HSB1
      double SIG1,SIG2,SIG3,SIG4;
      double WQ1,WQ2,WQ3,WQ4,WQ5,WQ6,WQ7,WQ8;
      double WT1,WT2,WT3,WT4,WT5,WT6,WT7,WT8;
 

      Upp(s,t,u,&WQ1,&WQ2,&WQ3,&WQ4,&WQ5,&WQ6,&WQ7,&WQ8);
      Upp(s,u,t,&WT1,&WT2,&WT3,&WT4,&WT5,&WT6,&WT7,&WT8);

      SIG1=WQ2*(U1*U2*fu+UB1*UB2*fub+D1*D2*fd+DB1*DB2*fdb
              +S1*S2*fs+SB1*SB2*fsb)
         +WQ3*(U1*UB2*(fd+fs)+D1*DB2*(fu+fs)+S1*SB2*(fu+fd)
              +UB1*U2*(fdb+fsb)+DB1*D2*(fub+fsb)
              +SB1*S2*(fub+fdb))
         +WT3*(UB1*U2*(fd+fs)+DB1*D2*(fu+fs)+SB1*S2*(fu+fd)
              +U1*UB2*(fdb+fsb)+D1*DB2*(fub+fsb)
              +S1*SB2*(fub+fdb));

      SIG2=WQ4*(U1*UB2*fu+D1*DB2*fd+S1*SB2*fs
              +UB1*U2*fub+DB1*D2*fdb+SB1*S2*fsb)
         +WT4*(UB1*U2*fu+DB1*D2*fd+SB1*S2*fs
              +U1*UB2*fub+D1*DB2*fdb+S1*SB2*fsb)
         +WQ1*(U1*(D2+DB2+S2+SB2)*fu
              +D1*(U2+UB2+S2+SB2)*fd
              +UB1*(D2+DB2+S2+SB2)*fub
              +DB1*(U2+UB2+S2+SB2)*fdb
              +S1*(U2+UB2+D2+DB2)*fs
              +SB1*(U2+UB2+D2+DB2)*fsb);

      SIG3=WT1*((D1+DB1+S1+SB1)*U2*fu
              +(U1+UB1+S1+SB1)*D2*fd
              +(D1+DB1+S1+SB1)*UB2*fub
              +(U1+UB1+S1+SB1)*DB2*fdb
              +(U1+UB1+D1+DB1)*S2*fs
              +(U1+UB1+D1+DB1)*SB2*fsb)
         +WQ5*(U1*UB2+D1*DB2+S1*SB2)*fg
         +WT5*(UB1*U2+DB1*D2+SB1*S2)*fg;

      SIG4=WQ6*GL1*GL2*(fu+fd+fs)
         +WT6*GL1*GL2*(fub+fdb+fsb)
         +WQ7*(GL2*(U1*fu+D1*fd+S1*fs+UB1*fub+DB1*fdb+SB1*fsb)
              +GL1*fg*(U2+UB2+D2+DB2+S2+SB2))
         +WT7*(GL1*(U2*fu+D2*fd+S2*fs+UB2*fub+DB2*fdb+SB2*fsb)
              +GL2*fg*(U1+UB1+D1+DB1+S1+SB1))
         +WQ8*GL1*GL2*fg;

      double SIG = SIG1+SIG2+SIG3+SIG4;

      return SIG;
};    


 
 
