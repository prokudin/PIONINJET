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
#include "gauss_2d_sphere.h"


using namespace std;


 

//C-------------------------------------------------------------------

      int Colpp(double s,double t,double u,double* WQ1,double* WQ2,double* WQ3,
	  double* WQ4,double* WQ5,double* WQ6,double* WQ7,double* WQ8){
 
//c     qq'->qq'
      *WQ1=8./9.*(s*u)/(-t*t);
//c     qq->qq
      *WQ2=8./9.*( (s*u)/(-t*t)-1./3.*s/(-t));
//c     qqb->q'qb'
      *WQ3=0.;
//c     qqb->qqb
      *WQ4=8./9.*( (s*u)/(-t*t)+1./3.*u/t);
//c     qqb->gg
      *WQ5=0;
//c     gg->qqb
      *WQ6=0.;
//c     qg->qg
      *WQ7=8./9.+2.*(s*u)/(-t*t);
//c     gg->gg
      *WQ8=0.; 
 
      return 0;
};





// C-------------------------------------------------------------------
// C     spin-dependent partonic cross section
// C-------------------------------------------------------------------
      double ColSIG(double s,double t,double u,
		double U1,double D1,double UB1,double DB1,double S1,double SB1,double GL1,
		double U2,double D2,double UB2,double DB2,double S2,double SB2,double GL2,
		double fu,double fub,double fd,double fdb,double fs,double fsb,double fg
){
       
      //real*8 HU1,HUB1,HD1,HDB1,HS1,HSB1
      double SIG1,SIG2,SIG3,SIG4;
      double WQ1,WQ2,WQ3,WQ4,WQ5,WQ6,WQ7,WQ8;
      double WT1,WT4;
//      double WT1,WT2,WT3,WT4,WT5,WT6,WT7,WT8;
 
     Colpp(s,t,u,&WQ1,&WQ2,&WQ3,&WQ4,&WQ5,&WQ6,&WQ7,&WQ8);
 
//C----q'q->qq'
      WT1 = 0.;
//C----qbq->qqb
      WT4 = -8./27.;

      SIG1=WQ2*(U1*U2*fu+UB1*UB2*fub+D1*D2*fd+DB1*DB2*fdb+
              S1*S2*fs+S1*S2*fsb);

      SIG2=WQ1*((U1*D2*fu+U1*DB2*fu+2.*U1*S2*fu)+
              2.*D1*S2*fd+(UB1*D2*fub+UB1*DB2*fub+
              2.*UB1*S2*fub)+2.*DB1*S2*fdb+
              S1*fs*(U2+D2+UB2+DB2)+
              S1*fsb*(U2+D2+UB2+DB2)+
              D1*fd*(U2+UB2)+DB1*fdb*(U2+UB2));

      SIG3=WQ4*(U1*UB2*fu+D1*DB2*fd+UB1*U2*fub+
              DB1*D2*fdb+S1*S2*fsb)
          +WT4*(UB1*U2*fu+DB1*D2*fd+S1*S2*fs+U1*UB2*fub+
              D1*DB2*fdb+S1*S2*fsb);

      SIG4=WQ7*(GL2*(U1*fu+D1*fd+S1*fs+UB1*fub+DB1*fdb+S1*fsb)); 
 
      double SIG=SIG1+SIG2+SIG3+SIG4;
 
      return SIG;
       
}







//C------------------It is the difference of hard factors SH1perp and SHFU-------------------------------------------------

      int ColppPITONYAK(double s,double t,double u,double* WQ1,double* WQ2,double* WQ3,
	  double* WQ4,double* WQ5,double* WQ6,double* WQ7,double* WQ8){
 
//c     qq'->qq'
      *WQ1=1./9. * s * (2. * t - u)/(t * t * t) + (t + u)/(t * t);
//c     qq->qq
      *WQ2=1./27. * s * (u - t)/(t * t * u)  + 1./9. * s * (2. * t -u)/(t * t * t) - s /(t * t);
//c     qqb->q'qb'
      *WQ3=0.;
//c     qqb->qqb
      *WQ4=1./27. * (u + t)/(t * t)  - 1./9. * s * (u  + 2. * t)/(t * t * t) + 1./3. * 1./t;
      //what about  qbq->qbq
//c     qqb->gg
      *WQ5=0;
//c     gg->qqb
      *WQ6=0.;
//c     qg->qg
      *WQ7=1./9. * 1./t + 1./8. * s * (s-u)/(t * t * t) + s * s / (t * t * u);
//c     gg->gg
      *WQ8=0.; 
 
      return 0;
};


// C-------------------------------------------------------------------
// C     spin-dependent partonic cross section
// C-------------------------------------------------------------------
      double ColSIGPITONYAK(double s,double t,double u,
		double U1,double D1,double UB1,double DB1,double S1,double SB1,double GL1,
		double U2,double D2,double UB2,double DB2,double S2,double SB2,double GL2,
		double fu,double fub,double fd,double fdb,double fs,double fsb,double fg
){
       
      //real*8 HU1,HUB1,HD1,HDB1,HS1,HSB1
      double SIG1,SIG2,SIG3,SIG4,SIG5;
      double WQ1,WQ2,WQ3,WQ4,WQ5,WQ6,WQ7,WQ8,WQ55;
      double WT1,WT4;
//      double WT1,WT2,WT3,WT4,WT5,WT6,WT7,WT8;
 
     ColppPITONYAK(s,t,u,&WQ1,&WQ2,&WQ3,&WQ4,&WQ5,&WQ6,&WQ7,&WQ8);
 
//C----q'q->qq'
      WT1 = 0.;
//C----qbq->qqb
      WT4 = 1./27. * (u + t)/(t * u) + 1./3. * 1./t;
//C----qqb'->qqb'  
	  WQ55 = -1./9. * (u + 2. * t)/(t * t * t) * s;    

      SIG1=WQ2*(U1*U2*fu+UB1*UB2*fub+D1*D2*fd+DB1*DB2*fdb+
              S1*S2*fs+S1*S2*fsb);

      SIG2=WQ1*((U1*D2*fu+U1*S2*fu)+
              D1*S2*fd+(UB1*D2*fub+UB1*DB2*fub+
              2.*UB1*S2*fub)+2.*DB1*S2*fdb+
              S1*fs*(U2+D2)+
              S1*fsb*(U2+D2+UB2+DB2)+
              D1*fd*U2+DB1*fdb*(U2+UB2));
      // contribution from   qqb'->qqb'      
      SIG5 = WQ55 * (  U1*DB2*fu + U1*S2*fu +
      		  D1*S2*fd + S1*fs*(UB2+DB2) +
      		  D1*fd*UB2 );

      SIG3=WQ4*(U1*UB2*fu+D1*DB2*fd+UB1*U2*fub+
              DB1*D2*fdb+S1*S2*fsb)
          +WT4*(UB1*U2*fu+DB1*D2*fd+S1*S2*fs+U1*UB2*fub+
              D1*DB2*fdb+S1*S2*fsb);

      SIG4=WQ7*(GL2*(U1*fu+D1*fd+S1*fs+UB1*fub+DB1*fdb+S1*fsb)); 
 
      double SIG=SIG1+SIG2+SIG3+SIG4+SIG5;

      return SIG;
       
}

//C------------------It is the difference of hard factors SH and SHFU-------------------------------------------------

      int ColppPITONYAK_SH(double s,double t,double u,double* WQ1,double* WQ2,double* WQ3,
	  double* WQ4,double* WQ5,double* WQ6,double* WQ7,double* WQ8){
 
//c     qq'->qq'
      *WQ1=1./9. * s * ( t - 2. * u)/(t * t * t) + s*u/(t * t * t);
//c     qq->qq
      *WQ2=-1./27. * s * (t - 3. * u)/(2.* t * t * u)  + 1./9. * s * (t - 2. * u)/(2. * t * t * u) + s * u /(t * t * t);
//c     qqb->q'qb'
      *WQ3=0.;
//c     qqb->qqb
      *WQ4=1./27. * (4. * t + 3. * u)/(t * t)  - 1./9. * s * (3. * t  + 2. * u)/(t * t * t) - 1./3. * u/ (2.*t*t) + s*(t+u)/(t*t*t);
      //what about  qbq->qbq
//c     qqb->gg
      *WQ5=0;
//c     gg->qqb
      *WQ6=0.;
//c     qg->qg
      *WQ7=1./8. * s * (s-u)/(t * t * t) + 1./9. * (u-s) / (2.* t * u) + (u-s)*(t*t-2.*t*u-2.*u*u)/(2.*t*t*t*u);
//c     gg->gg
      *WQ8=0.; 
 
      return 0;
};

// C-------------------------------------------------------------------
// C     spin-dependent partonic cross section for H term
// C-------------------------------------------------------------------
      double ColSIGPITONYAK_SH(double s,double t,double u,
		double U1,double D1,double UB1,double DB1,double S1,double SB1,double GL1,
		double U2,double D2,double UB2,double DB2,double S2,double SB2,double GL2,
		double fu,double fub,double fd,double fdb,double fs,double fsb,double fg
){
       
      //real*8 HU1,HUB1,HD1,HDB1,HS1,HSB1
      double SIG1,SIG2,SIG3,SIG4,SIG5;
      double WQ1,WQ2,WQ3,WQ4,WQ5,WQ6,WQ7,WQ8,WQ55;
      double WT1,WT4;
//      double WT1,WT2,WT3,WT4,WT5,WT6,WT7,WT8;
 
     ColppPITONYAK_SH(s,t,u,&WQ1,&WQ2,&WQ3,&WQ4,&WQ5,&WQ6,&WQ7,&WQ8);
 
//C----q'q->qq'
      WT1 = 0.;
//C----qbq->qqb
      WT4 = ( 1./27. + 1./3. ) * (t + 2.*u)/(2.* t * u);
//C----qqb'->qqb'  
	  WQ55 = -1./9. * s * (3. * t + 2. * u)/(t * t * t) + s * (t + u)/(t*t*t);    

      SIG1=WQ2*(U1*U2*fu+UB1*UB2*fub+D1*D2*fd+DB1*DB2*fdb+
              S1*S2*fs+S1*S2*fsb);

      SIG2=WQ1*((U1*D2*fu+U1*S2*fu)+
              D1*S2*fd+(UB1*D2*fub+UB1*DB2*fub+
              2.*UB1*S2*fub)+2.*DB1*S2*fdb+
              S1*fs*(U2+D2)+
              S1*fsb*(U2+D2+UB2+DB2)+
              D1*fd*U2+DB1*fdb*(U2+UB2));
      // contribution from   qqb'->qqb'      
      SIG5 = WQ55 * (  U1*DB2*fu + U1*S2*fu +
      		  D1*S2*fd + S1*fs*(UB2+DB2) +
      		  D1*fd*UB2 );

      SIG3=WQ4*(U1*UB2*fu+D1*DB2*fd+UB1*U2*fub+
              DB1*D2*fdb+S1*S2*fsb)
          +WT4*(UB1*U2*fu+DB1*D2*fd+S1*S2*fs+U1*UB2*fub+
              D1*DB2*fdb+S1*S2*fsb);

      SIG4=WQ7*(GL2*(U1*fu+D1*fd+S1*fs+UB1*fub+DB1*fdb+S1*fsb)); 
 
      double SIG=SIG1+SIG2+SIG3+SIG4+SIG5;

      return SIG;
       
}



//========================================================= Approximate collins in PP...
int sidis::approximate_collins_pp_calculation( double xa, double xb, double zc, double y, double pt, double Q2 )
{
//  double integrand = 0.;

  double s = sidis::GetEnergy();
  double sqrts = sqrt(s);

//  double MASS = Target.GetMass();


  double Q  = sqrt( Q2 );

  double FU = 0,FD = 0,FUB = 0,FDB = 0,FS = 0,FSB = 0,FGL = 0; 
  double U1 = 0,D1 = 0,UB1 = 0,DB1 = 0,S1 = 0,SB1 = 0,GL1 = 0; 
  double U2 = 0,D2 = 0,UB2 = 0,DB2 = 0,S2 = 0,SB2 = 0,GL2 = 0.; 

  double PU1,PD1,PUB1,PDB1,PS1,PSB1,PGL1; 
  double MU1,MD1,MUB1,MDB1,MS1,MSB1,MGL1; 
  double DU1,DD1,DUB1,DDB1,DS1,DSB1,DGL1; 

  double cross_section;
 
  double  sh = xa * xb * s;
  double  uh = -pt/zc * xb * sqrts * exp(y);
  double  th = -pt/zc * xa * sqrts * exp(-y);

  double coeff, coeff_col, coeff_unp;

  double delta = 1.e-3;
  double zcp = zc;
  double zcm = zc;
  int count = 10;
//  bool condition;

  switch( effect ){
  case collins_effect: // COLLINS

    switch( Beam.GetPolarization() ){
    case transverse:

    TransversityDistribution( Beam , xa, Q2 );  // 1st hadron transversity

      U1  = Beam.partcontent.up;
      D1  = Beam.partcontent.down;     
      UB1 = Beam.partcontent.anti_up;   
      DB1 = Beam.partcontent.anti_down;  
      S1  = Beam.partcontent.strange;   
      SB1 = Beam.partcontent.anti_strange;  
      GL1 = Beam.partcontent.glu; 


    Unpolarised( Target , xb, Q2 );  // 2nd hadron

      U2  = Target.partcontent.up;
      D2  = Target.partcontent.down;     
      UB2 = Target.partcontent.anti_up;   
      DB2 = Target.partcontent.anti_down;  
      S2  = Target.partcontent.strange;   
      SB2 = Target.partcontent.anti_strange;  
      GL2 = Target.partcontent.glu; 

//Collins fragmentation
//C-------------------------------------------------------------------
//C     -z*d(H(z)/z^2)/dz
//C-------------------------------------------------------------------      
 
      zcp = zc*(1.+delta);
      zcm = zc*(1.-delta);

      count = 0;
      while( zcp > 1. ){
	delta = delta/2.;
        zcp = zc*(1.+delta);
	count = count +1;
	if (count >= 10) zcp = 1.; 
      };
 
      CollinsDistributionFirstMoment( Produced, zcp, Q2 ); //  FIRST MOMENT OF COLLINS
// see definition of Kang H(z) = -2 z mpion H1perp(1)
      PU1  = (-2.*zcp*mpion) * Produced.fragmentation.up;
      PD1  = (-2.*zcp*mpion) * Produced.fragmentation.down;     
      PUB1 = (-2.*zcp*mpion) * Produced.fragmentation.anti_up;   
      PDB1 = (-2.*zcp*mpion) * Produced.fragmentation.anti_down;  
      PS1  = (-2.*zcp*mpion) * Produced.fragmentation.strange;   
      PSB1 = (-2.*zcp*mpion) * Produced.fragmentation.anti_strange;  
      PGL1 = (-2.*zcp*mpion) * Produced.fragmentation.glu; 

      CollinsDistributionFirstMoment( Produced, zcm, Q2 ); //  FIRST MOMENT OF COLLINS

      MU1  = (-2.*zcm*mpion) * Produced.fragmentation.up;
      MD1  = (-2.*zcm*mpion) * Produced.fragmentation.down;     
      MUB1 = (-2.*zcm*mpion) * Produced.fragmentation.anti_up;   
      MDB1 = (-2.*zcm*mpion) * Produced.fragmentation.anti_down;  
      MS1  = (-2.*zcm*mpion) * Produced.fragmentation.strange;   
      MSB1 = (-2.*zcm*mpion) * Produced.fragmentation.anti_strange;  
      MGL1 = (-2.*zcm*mpion) * Produced.fragmentation.glu; 

      DU1  = -zc * (PU1/pow2(zcp) - MU1/pow2(zcm))  / (2.*zc*delta);
      DD1  = -zc * (PD1/pow2(zcp) - MD1/pow2(zcm))  / (2.*zc*delta);     
      DUB1 = -zc * (PUB1/pow2(zcp)- MUB1/pow2(zcm)) / (2.*zc*delta);   
      DDB1 = -zc * (PDB1/pow2(zcp)- MDB1/pow2(zcm)) / (2.*zc*delta) ;  
      DS1  = -zc * (PS1/pow2(zcp) - MS1/pow2(zcm))  / (2.*zc*delta);   
      DSB1 = -zc * (PSB1/pow2(zcp)- MSB1/pow2(zcm)) / (2.*zc*delta);  
      DGL1 = -zc * (PGL1/pow2(zcp)- MGL1/pow2(zcm)) / (2.*zc*delta); 

      FU  =   DU1;  
      FD  =   DD1;     
      FUB =   UB1;   
      FDB =   DDB1;  
      FS  =   DS1;   
      FSB =   DSB1;  
      FGL =   DGL1; 
  
      break;
    default:
	cerr << " No unp or long here!!  int sidis::approximate_collins_pp_calculation( double xa, double xb, double zc, double y, double pt, double Q2 ) " << endl;
    break;

    };
    
       cross_section = ColSIG( sh,  th,  uh,
		  U1,  D1,  UB1,  DB1,  S1,  SB1,  GL1,
		  U2,  D2,  UB2,  DB2,  S2,  SB2,  GL2,
		  FU,  FUB,  FD,  FDB,  FS,  FSB,  FGL);


    switch( dependence ){
    case pt_dependence:
//       coeff_col = (xa - xb)/( xa*(-uh) + xb*(-th) )/(xa*xb*zc*zc) * (-pt) ;
       coeff_col = 1./zc/(xa*xb) * (xa - xb)/( xa*(-uh) + xb*(-th) ) * (-pt) ; //

       break;
    default: // x, z dependence
//       coeff_col = (xa - xb)/( xa*(-uh) + xb*(-th) )/(xa*xb*zc*zc) * (-pt) ;
       coeff_col = 1./zc/(xa*xb) * (xa - xb)/( xa*(-uh) + xb*(-th) ) * (-pt) ; //
     break;
    }
            
    coeff = coeff_col;

    break;

  case no_effect: // UNPOLARISED


    Unpolarised( Beam , xa, Q2 );  // 1st hadron

      U1  = Beam.partcontent.up;
      D1  = Beam.partcontent.down;     
      UB1 = Beam.partcontent.anti_up;   
      DB1 = Beam.partcontent.anti_down;  
      S1  = Beam.partcontent.strange;   
      SB1 = Beam.partcontent.anti_strange;  
      GL1 = Beam.partcontent.glu; 

    Unpolarised( Target , xb, Q2 );  // 2nd hadron

      U2  = Target.partcontent.up;
      D2  = Target.partcontent.down;     
      UB2 = Target.partcontent.anti_up;   
      DB2 = Target.partcontent.anti_down;  
      S2  = Target.partcontent.strange;   
      SB2 = Target.partcontent.anti_strange;  
      GL2 = Target.partcontent.glu; 


    Fragmentation( Produced, zc, Q2);


      FU  = Produced.fragmentation.up;
      FD  = Produced.fragmentation.down;     
      FUB = Produced.fragmentation.anti_up;   
      FDB = Produced.fragmentation.anti_down;  
      FS  = Produced.fragmentation.strange;   
      FSB = Produced.fragmentation.anti_strange;  
      FGL = Produced.fragmentation.glu; 


  
      cross_section = UUSIG( sh,  th,  uh,
		  U1,  D1,  UB1,  DB1,  S1,  SB1,  GL1,
		  U2,  D2,  UB2,  DB2,  S2,  SB2,  GL2,
		  FU,  FUB,  FD,  FDB,  FS,  FSB,  FGL);

    switch( dependence ){
    case pt_dependence:
      coeff_unp = 1./zc/(xa*xb) ; 
       break;
    default: // x, z dependence
      coeff_unp = 1./zc/(xa*xb) ; 
     break;
    }

    coeff = coeff_unp;

    break;
  default: // 
    
    cout << "Error in sidis::approximate_collins_pp_calculation( double xa, double xb, double z, double pt, double Q2 ) no effect is defined" << endl;
    return 1;
    
    break;
  };



  sidis::cross_section = cross_section;  

  sidis::cross_section *= coeff * pow2( alphasQCD( Q ) / s ); 

  return 0;
}




//========================================================= Approximate collins in PP at NLO...
int sidis::approximate_collins_pp_calculation_NLO( double xa, double xb, double zc, double y, double pt, double Q2 )
{
//  double integrand = 0.;

  double s = sidis::GetEnergy();
  double sqrts = sqrt(s);

//  double MASS = Target.GetMass();


  double Q  = sqrt( Q2 );

  double FU = 0,FD = 0,FUB = 0,FDB = 0,FS = 0,FSB = 0,FGL = 0; 
  double FU1 = 0,FD1 = 0,FUB1 = 0,FDB1 = 0,FS1 = 0,FSB1 = 0,FGL1 = 0; 

  double U1 = 0,D1 = 0,UB1 = 0,DB1 = 0,S1 = 0,SB1 = 0,GL1 = 0; 
  double U2 = 0,D2 = 0,UB2 = 0,DB2 = 0,S2 = 0,SB2 = 0,GL2 = 0.; 

  double PU1,PD1,PUB1,PDB1,PS1,PSB1,PGL1; 
  double MU1,MD1,MUB1,MDB1,MS1,MSB1,MGL1; 
  double DU1,DD1,DUB1,DDB1,DS1,DSB1,DGL1; 

  double cross_section;
 
  double  sh = xa * xb * s;
  double  uh = -pt/zc * xb * sqrts * exp(y);
  double  th = -pt/zc * xa * sqrts * exp(-y);

  double coeff, coeff_col, coeff_unp;

  //double delta = 1.e-3;
  double delta = 1.e-3;
  double zcp = zc;
  double zcm = zc;
  int count = 10;
//  bool condition;

  switch( effect ){
  case collins_effect: // COLLINS

    switch( Beam.GetPolarization() ){
    case transverse:

    TransversityDistributionNLO( Beam , xa, Q2 );  // 1st hadron transversity

      U1  = Beam.partcontent.up;
      D1  = Beam.partcontent.down;     
      UB1 = Beam.partcontent.anti_up;   
      DB1 = Beam.partcontent.anti_down;  
      S1  = Beam.partcontent.strange;   
      SB1 = Beam.partcontent.anti_strange;  
      GL1 = Beam.partcontent.glu; 
      


    Unpolarised( Target , xb, Q2 );  // 2nd hadron

      U2  = Target.partcontent.up;
      D2  = Target.partcontent.down;     
      UB2 = Target.partcontent.anti_up;   
      DB2 = Target.partcontent.anti_down;  
      S2  = Target.partcontent.strange;   
      SB2 = Target.partcontent.anti_strange;  
      GL2 = Target.partcontent.glu; 

//Collins fragmentation
//C-------------------------------------------------------------------
//C     -z*d(H(z)/z^2)/dz
//C-------------------------------------------------------------------      
 
      zcp = zc*(1.+delta);
      zcm = zc*(1.-delta);

      count = 0;
      
      while( zcp > 1. ){
	    delta = delta/2.;
        zcp = zc*(1.+delta);
        zcm = zc*(1.-delta);
	    count = count + 1;
	    if (count >= 10) {
	        zcp = zcm;
	  		DU1  = 0;
      		DD1  = 0;     
      		DUB1 = 0;   
      		DDB1 = 0;  
      		DS1  = 0;   
      		DSB1 = 0;  
      		DGL1 = 0; 
	       }; 
      };
 
 
      CollinsDistributionNLO( Produced, zc, Q2 ); //  FIRST MOMENT OF COLLINS
 
// see definition of Kang or FengH(z) = -2 z mpion H1perp(1) Eq(58) from our paper with Feng
//This formula uses exactly H3: (-2 Mh) is already in the formula by Metz Pitonyak
      FU1  =   (-2.)/zc * Produced.fragmentation.up; // (-1.)/(zc) H(3) = (-2.) * (-2 H1perp)
      FD1  =   (-2.)/zc * Produced.fragmentation.down;     
      FUB1 =   (-2.)/zc * Produced.fragmentation.anti_up;   
      FDB1 =   (-2.)/zc * Produced.fragmentation.anti_down;  
      FS1  =   (-2.)/zc * Produced.fragmentation.strange;   
      FSB1 =   (-2.)/zc * Produced.fragmentation.anti_strange;  
      FGL1 =   (-2.)/zc * Produced.fragmentation.glu; 

 
 
      if(count < 10){
      CollinsDistributionNLO( Produced, zcp, Q2 ); //  FIRST MOMENT OF COLLINS
 
// see definition of Kang or FengH(z) = -2 z mpion H1perp(1) Eq(58) from our paper with Feng
//This formula uses exactly H3:
      PU1  =   Produced.fragmentation.up;
      PD1  =   Produced.fragmentation.down;     
      PUB1 =   Produced.fragmentation.anti_up;   
      PDB1 =   Produced.fragmentation.anti_down;  
      PS1  =   Produced.fragmentation.strange;   
      PSB1 =   Produced.fragmentation.anti_strange;  
      PGL1 =   Produced.fragmentation.glu; 

      CollinsDistributionNLO( Produced, zcm, Q2 ); //  FIRST MOMENT OF COLLINS

      MU1  =   Produced.fragmentation.up;
      MD1  =   Produced.fragmentation.down;     
      MUB1 =   Produced.fragmentation.anti_up;   
      MDB1 =   Produced.fragmentation.anti_down;  
      MS1  =   Produced.fragmentation.strange;   
      MSB1 =   Produced.fragmentation.anti_strange;  
      MGL1 =   Produced.fragmentation.glu; 

      DU1  = -zc * (PU1/pow2(zcp) - MU1/pow2(zcm))  / (2.*zc*delta);
      DD1  = -zc * (PD1/pow2(zcp) - MD1/pow2(zcm))  / (2.*zc*delta);     
      DUB1 = -zc * (PUB1/pow2(zcp)- MUB1/pow2(zcm)) / (2.*zc*delta);   
      DDB1 = -zc * (PDB1/pow2(zcp)- MDB1/pow2(zcm)) / (2.*zc*delta) ;  
      DS1  = -zc * (PS1/pow2(zcp) - MS1/pow2(zcm))  / (2.*zc*delta);   
      DSB1 = -zc * (PSB1/pow2(zcp)- MSB1/pow2(zcm)) / (2.*zc*delta);  
      DGL1 = -zc * (PGL1/pow2(zcp)- MGL1/pow2(zcm)) / (2.*zc*delta); 
      };

      FU  =   DU1;  
      FD  =   DD1;     
      FUB =   UB1;   
      FDB =   DDB1;  
      FS  =   DS1;   
      FSB =   DSB1;  
      FGL =   DGL1; 
  
      break;
    default:
	cerr << " No unp or long here!!  int sidis::approximate_collins_pp_calculation( double xa, double xb, double zc, double y, double pt, double Q2 ) " << endl;
    break;

    };



    if(!VALENCE_SIDIS){  
//       cross_section = ColSIG( sh,  th,  uh,
// 		  U1,  D1,  UB1,  DB1,  S1,  SB1,  GL1,
// 		  U2,  D2,  UB2,  DB2,  S2,  SB2,  GL2,
// 		  FU,  FUB,  FD,  FDB,  FS,  FSB,  FGL); // original calculation, only H1perp
	  cross_section = 1.* ColSIGPITONYAK( sh,  th,  uh,
		  U1,  D1,  UB1,  DB1,  S1,  SB1,  GL1,
		  U2,  D2,  UB2,  DB2,  S2,  SB2,  GL2,
		  FU,  FUB,  FD,  FDB,  FS,  FSB,  FGL)   +
		  1.* ColSIGPITONYAK_SH( sh,  th,  uh,  // If multiplied by "0" it is the original version, it corresponds to the 3rd term in Fengs paper, not suppressed by A
		  U1,  D1,  UB1,  DB1,  S1,  SB1,  GL1,
		  U2,  D2,  UB2,  DB2,  S2,  SB2,  GL2,
		  FU1,  FUB1,  FD1,  FDB1,  FS1,  FSB1,  FGL1) ;
    } else if(VALENCE_SIDIS){ // WE USE VALENCE STANDARD PARAMETRIZATION!!!!
//       cross_section = ColSIG( sh,  th,  uh,
// 		  U1+UB1,  D1+DB1,  UB1,  DB1,  S1,  SB1,  GL1,
// 		  U2+UB2,  D2+DB2,  UB2,  DB2,  S2,  SB2,  GL2,
// 		  FU,  FUB,  FD,  FDB,  FS,  FSB,  FGL);       
      cross_section = ColSIGPITONYAK( sh,  th,  uh,
		  U1+UB1,  D1+DB1,  UB1,  DB1,  S1,  SB1,  GL1,
		  U2+UB2,  D2+DB2,  UB2,  DB2,  S2,  SB2,  GL2,
		  FU,  FUB,  FD,  FDB,  FS,  FSB,  FGL)    +
		  ColSIGPITONYAK_SH( sh,  th,  uh,
		  U1+UB1,  D1+DB1,  UB1,  DB1,  S1,  SB1,  GL1,
		  U2+UB2,  D2+DB2,  UB2,  DB2,  S2,  SB2,  GL2,
		  FU1,  FUB1,  FD1,  FDB1,  FS1,  FSB1,  FGL1)  ;       

    }

    

    switch( dependence ){
    case pt_dependence:
       // coeff_col = 1./zc/(xa*xb) * (xa - xb)/( xa*(-uh) + xb*(-th) ) * (-pt) ; // IF CALCULATION WITH STANDARD
       coeff_col = 1./zc/(xa*xb) * pt ;
       break;
    default: // x, z dependence
       // coeff_col = 1./zc/(xa*xb) * (xa - xb)/( xa*(-uh) + xb*(-th) ) * (-pt) ; // IF CALCULATION WITH STANDARD
       coeff_col = 1./zc/(xa*xb) * pt ; // MAIN: CALCULATION WITH LIRS THEN (xa - xb)/( xa*(-uh) + xb*(-th) ) cancels out
     break;
    }
            
    coeff = coeff_col;

    break;

  case no_effect: // UNPOLARISED


    Unpolarised( Beam , xa, Q2 );  // 1st hadron

      U1  = Beam.partcontent.up;
      D1  = Beam.partcontent.down;     
      UB1 = Beam.partcontent.anti_up;   
      DB1 = Beam.partcontent.anti_down;  
      S1  = Beam.partcontent.strange;   
      SB1 = Beam.partcontent.anti_strange;  
      GL1 = Beam.partcontent.glu; 

    Unpolarised( Target , xb, Q2 );  // 2nd hadron

      U2  = Target.partcontent.up;
      D2  = Target.partcontent.down;     
      UB2 = Target.partcontent.anti_up;   
      DB2 = Target.partcontent.anti_down;  
      S2  = Target.partcontent.strange;   
      SB2 = Target.partcontent.anti_strange;  
      GL2 = Target.partcontent.glu; 


    Fragmentation( Produced, zc, Q2);


      FU  = Produced.fragmentation.up;
      FD  = Produced.fragmentation.down;     
      FUB = Produced.fragmentation.anti_up;   
      FDB = Produced.fragmentation.anti_down;  
      FS  = Produced.fragmentation.strange;   
      FSB = Produced.fragmentation.anti_strange;  
      FGL = Produced.fragmentation.glu; 


  
		  
	if(!VALENCE_SIDIS){  
      cross_section = UUSIG( sh,  th,  uh,
		  U1,  D1,  UB1,  DB1,  S1,  SB1,  GL1,
		  U2,  D2,  UB2,  DB2,  S2,  SB2,  GL2,
		  FU,  FUB,  FD,  FDB,  FS,  FSB,  FGL);
    } else if(VALENCE_SIDIS){ // WE USE VALENCE STANDARD PARAMETRIZATION!!!!
      cross_section = UUSIG( sh,  th,  uh,
		  U1+UB1,  D1+DB1,  UB1,  DB1,  S1,  SB1,  GL1,
		  U2+UB2,  D2+DB2,  UB2,  DB2,  S2,  SB2,  GL2,
		  FU,  FUB,  FD,  FDB,  FS,  FSB,  FGL);       
    }
	  
		  
		  

    switch( dependence ){
    case pt_dependence:
      coeff_unp = 1./zc/(xa*xb) ; // MAKE IT THE SAME AS IN BOER MULDERS
       break;
    default: // x, z dependence
//      coeff_unp = 1./zc/(pow2(xa)*pow2(xb)) ;  // MAIN ?
      coeff_unp = 1./zc/(xa*xb) ; // CORRECT 
     break;
    }

    coeff = coeff_unp;

    break;
  default: // 
    
    cout << "Error in sidis::approximate_collins_pp_calculation( double xa, double xb, double z, double pt, double Q2 ) no effect is defined" << endl;
    return 1;
    
    break;
  };



  sidis::cross_section = cross_section;  

  sidis::cross_section *= coeff * pow2( alphasQCD( Q ) / s ); 

  return 0;
}


// used in numerical integration by Holoborodko's method:
  struct PAR_FOR_AN { 
    double xF, zc, y, pt, Q2;
  };


// Let us also implement methods of numerical integration by Pavel Holoborodko
// Let us set up the integrand...
double an_gauss (double xa, double xb, void * params) {
  sidis& SIDIS = *sidis::Instance();


  struct PAR_FOR_AN * pars = (struct PAR_FOR_AN *) params;
   
  double SS = 0, TT = 0, UU = 0;


  double xF = pars->xF;
  double zc = pars->zc;
  double y  = pars->y;
  double pt = pars->pt;
  double Q2 = pars->Q2;
 
  int test;
  double result = 0.;
  
  double s = SIDIS.GetEnergy();
  double sqrts = sqrt(s); 


  SS = s;
  TT = -sqrts * pt * exp(-y);
  UU = -sqrts * pt * exp(y);

  zc = -(TT/xb + UU/xa)/SS;
 
  // We integrate in the disc (x-Xc)^2+(y-Yc)^2<=R^2, using R =1/sqrt(2) such that   
  // The sphere contains the region I am interested in:  xamin < xa < 1, xbmin < xb < 1
  // We should filter all variables, however...


  if(
     (-UU/(SS+TT) < xa && xa  < 1.) &&
     (-TT/(SS+UU/xa) < xb && xb  < 1.) &&
     (0. < zc && zc  < 1.) &&
     SIDIS.Q2_cut(Q2) ){
      
    
     if( ( SIDIS.approximate_collins_pp_calculation_NLO( xa, xb, zc, y, pt, Q2) ) || 0) {
		cerr << "approximate_collins_pp.cpp: ERROR wrong result in Holoborodko!" << endl;     
      }

    result =  SIDIS.GetCrossSection();
  //cout << " xa = " << xa << " xb = " << xb << " zc = " << zc << " y = " << y << " pt = " << pt << " Q2 = " << Q2<< " res = "  << result << endl;
   
   }
   
   return result;
   
   
};



// calculates numerator/denominator for AN depends on settings of sidis.effect
double an_numerator_denominator_gauss(double xF) 
{
  sidis& SIDIS = *sidis::Instance();
    
  double result;
  struct PAR_FOR_AN pars;
  
  double xa = 1., xb = 1., zc = 1., pt = 1., Q2 = 0.; // y here stands for rapidity!
 
  double SS = 0, TT = 0, UU = 0;

  
  double s = SIDIS.GetEnergy();
  double sqrts = sqrt(s); 

  double y = SIDIS.Average.y; // this is rapidity!

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


 
  pars.xF = xF; 
  pars.zc = zc; 
  pars.y  = y;
  pars.pt = pt;
  pars.Q2 = Q2;


  //int key = 50;
  //int key = 704;
  int key = 256;
 
 
/* 2D Numerical computation of int(f(x,y), (x-Xc)^2+(y-Yc)^2<=R^2) by Gauss-type n-th order high precision quadrature 
   Cubature uses n^2 points and exact on all x^i*y^j with i+j<=2*n-1 

	[in]n       - quadrature order
	[in]f       - integrand
	[in]data    - pointer on user-defined data which will 
				  be passed to f every time it called (as third parameter).
	[in]R       - radius of the sphere
	[in](Xc,Yc) - center of the sphere. Unit disk is R=1, Xc=Yc=0. 

return:
		-computed integral value or -1.0 if n order quadrature is not supported
*/

// The sphere of radius 2 does contain the region I am interested in:  xamin < xa < 1, xbmin < xb < 1...

  result = gauss_product_2D_sphere(key, an_gauss, &pars, 1./sqrt(2.), -0.5, -0.5);

  //result = gauss_product_2D_sphere(key, an_gauss, &pars, sqrt(2.), 0., 0.);
    

  return  result; // 
}



// calculates numerator/denominator for AN depends on settings of sidis.effect PT dependence
double an_numerator_denominator_gauss_pt(double pt) 
{
  sidis& SIDIS = *sidis::Instance();
    
  double result;
  struct PAR_FOR_AN pars;
  
  double xa = 1., xb = 1., zc = 1., Q2 = 0.; // y here stands for rapidity!
 
  double SS = 0, TT = 0, UU = 0;

  
  double s = SIDIS.GetEnergy();
  double sqrts = sqrt(s); 


  double xF = SIDIS.Average.x; // this is xF
  
  //calculate rapidity
  double y = asinh( sqrts * xF/ (2. * pt) ); // pseudorapidity


  Q2 = pt*pt;
//    Q2 = 4.*pt*pt; // Q = 2 pt
//    Q2 = pt*pt/4.; // Q = pt/2.

  SS = s;
  TT = -sqrts * pt * exp(-y);
  UU = -sqrts * pt * exp(y);

 
  pars.xF = xF; 
  pars.zc = zc; 
  pars.y  = y;
  pars.pt = pt;
  pars.Q2 = Q2;


  //int key = 50;
  //int key = 704;
  int key = 256;
 
/* 2D Numerical computation of int(f(x,y), (x-Xc)^2+(y-Yc)^2<=R^2) by Gauss-type n-th order high precision quadrature 
   Cubature uses n^2 points and exact on all x^i*y^j with i+j<=2*n-1 

	[in]n       - quadrature order
	[in]f       - integrand
	[in]data    - pointer on user-defined data which will 
				  be passed to f every time it called (as third parameter).
	[in]R       - radius of the sphere
	[in](Xc,Yc) - center of the sphere. Unit disk is R=1, Xc=Yc=0. 

return:
		-computed integral value or -1.0 if n order quadrature is not supported
*/

// The sphere of radius 2 does contain the region I am interested in:  xamin < xa < 1, xbmin < xb < 1...

  result = gauss_product_2D_sphere(key, an_gauss, &pars, 1./sqrt(2.), -0.5, -0.5);

  //result = gauss_product_2D_sphere(key, an_gauss, &pars, sqrt(2.), 0., 0.);
    

  return  result; // 
}

