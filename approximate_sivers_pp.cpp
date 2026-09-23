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
//#include <gsl/gsl_sf_erf.h>
#include <gauss_legendre.h>
#include <approximate_sivers_pp.h>


using namespace std;







//C-------------------------------------------------------------------

      int Sivpp(double s,double t,double u,double* WQ1,double* WQ2,double* WQ3,
	  double* WQ4,double* WQ5,double* WQ6,double* WQ7,double* WQ8,double* WQ9,double* WQ10,
	  double* WQ11,double* WQ12,double* WQ13,double* WQ14,double* WQ15,double* WQ16,double*  WQ17,
	  double* WQ18,double* WQ19,double* WQ20,double* WQ21,double* WQ22,double* WQ23,double* WQ24){
      double Nc,htu;

      Nc  = 3.;
      htu = 1. + u/t;

//c     qq' ->qq'
      *WQ1=-1./pow2(Nc)*(s*s+u*u)/(t*t)
          -1./2./pow2(Nc)*(s*s+u*u)/(t*t)*htu ;
 //c     qqb'->qqb'
      *WQ2 = -(pow2(Nc)-2.)/(2.*pow2(Nc))*(s*s+u*u)/(t*t)
          -1./2./pow2(Nc)*(s*s+u*u)/(t*t)*htu ;
//c     qbq'->qbq'
      *WQ3 = -*WQ2;
//     (pow2(Nc)-2.)/(2.*pow2(Nc))*(s*s+u*u)/(t*t)
//          +1./2./pow2(Nc)*(s*s+u*u)/(t*t)*htu ;
//c     qbqb'->qbqb'
      *WQ4 = -*WQ1;
//    1./pow2(Nc)*(s*s+u*u)/(t*t)
//          +1./2./pow2(Nc)*(s*s+u*u)/(t*t)*htu ;
//c     qq'->q'q
      *WQ5 = -1./pow2(Nc)*(s*s+t*t)/(u*u)
          +(pow2(Nc)-2.)/(2.*pow2(Nc))*(s*s+t*t)/(u*u)*htu ;
//c     qbq'->q'qb
      *WQ6 = (pow2(Nc)-2.)/(2.*pow2(Nc))*(s*s+t*t)/(u*u)
          -1./pow2(Nc)*(s*s+t*t)/(u*u)*htu ;
//c     qqb'->qb'q
      *WQ7 = -*WQ6;
//-(pow2(Nc)-2.)/(2.*pow2(Nc))*(s*s+t*t)/(u*u)
//          +1./pow2(Nc)*(s*s+t*t)/(u*u)*htu ;
//c     qbqb'->qb'qb
      *WQ8= -*WQ5;
//     1./pow2(Nc)*(s*s+t*t)/(u*u)
//          -(pow2(Nc)-2.)/(2.*pow2(Nc))*(s*s+t*t)/(u*u)*htu ;
//c     qq ->qq
      *WQ9=-1./pow2(Nc)*((s*s+u*u)/(t*t)+(s*s+t*t)/(u*u))
                 +(pow2(Nc)+1.)/pow(Nc,3)*(s*s)/(t*u)
          +(-1./2./pow2(Nc)*(s*s+u*u)/(t*t)
                 +(pow2(Nc)-2.)/(2.*pow2(Nc))*(s*s+t*t)/(u*u)
                 +1./pow(Nc,3)*(s*s)/(t*u))*htu ;
//c     qbqb->qbqb
      *WQ10 = -*WQ9;
//1./pow2(Nc)*((s*s+u*u)/(t*t)+(s*s+t*t)/(u*u))
//                 -(pow2(Nc)+1.)/pow(Nc,3)*(s*s)/(t*u)
//          +(1./2./pow2(Nc)*(s*s+u*u)/(t*t)
//                 -(pow2(Nc)-2.)/(2.*pow2(Nc))*(s*s+t*t)/(u*u)
//                 -1./pow(Nc,3)*(s*s)/(t*u))*htu ;
//c     qqb->q'qb'
      *WQ11=1./2./pow2(Nc)*(t*t+u*u)/(s*s)
          +(pow2(Nc)-2.)/(2.*pow2(Nc))*(t*t+u*u)/(s*s)*htu ;
//c     qbq->q'qb'
      *WQ12=-1./2./pow2(Nc)*(t*t+u*u)/(s*s)
          -1./pow2(Nc)*(t*t+u*u)/(s*s)*htu ;
//c     qqb->qb'q'
      *WQ13 = -*WQ12;
//1./2./pow2(Nc)*(t*t+u*u)/(s*s)
//          +1./pow2(Nc)*(t*t+u*u)/(s*s)*htu ;
//c     qbq->qb'q'
      *WQ14 = -*WQ11;
//-1./2./pow2(Nc)*(t*t+u*u)/(s*s)
//         -(pow2(Nc)-2.)/(2.*pow2(Nc))*(t*t+u*u)/(s*s)*htu ;
//c     qqb->qqb
      *WQ15 = -(pow2(Nc)-2.)/(2.*pow2(Nc))*(s*s+u*u)/(t*t)
                +1./2./pow2(Nc)*(t*t+u*u)/(s*s)-1./pow(Nc,3)*(u*u)/(s*t)
          +(-1./2./pow2(Nc)*(s*s+u*u)/(t*t)
                +(pow2(Nc)-2.)/(2.*pow2(Nc))*(t*t+u*u)/(s*s)
                +1./pow(Nc,3)*(u*u)/(s*t))*htu ;
//c     qqb->qbq
      *WQ16 = -(pow2(Nc)-2.)/(2.*pow2(Nc))*(s*s+u*u)/(t*t)
                +1./2./pow2(Nc)*(t*t+u*u)/(s*s)-1./pow(Nc,3)*(t*t)/(s*u)
          +(1./pow2(Nc)*(s*s+u*u)/(t*t)
                +1./pow2(Nc)*(t*t+u*u)/(s*s)
                -(pow2(Nc)+1.)/pow(Nc,3)*(t*t)/(s*u))*htu ;
//c     qbq->qqb
      *WQ17 = -*WQ16;
//(pow2(Nc)-2.)/(2.*pow2(Nc))*(s*s+u*u)/(t*t)
//               -1./2./pow2(Nc)*(t*t+u*u)/(s*s)+1./pow(Nc,3)*(t*t)/(s*u)
//          +(-1./pow2(Nc)*(s*s+u*u)/(t*t)
//                -1./pow2(Nc)*(t*t+u*u)/(s*s)
//                +(pow2(Nc)+1.)/pow(Nc,3)*(t*t)/(s*u))*htu ;
//c     qbq->qbq
      *WQ18 = -*WQ15;
//(pow2(Nc)-2.)/(2.*pow2(Nc))*(s*s+u*u)/(t*t)
//                -1./2./pow2(Nc)*(t*t+u*u)/(s*s)+1./pow(Nc,3)*(u*u)/(s*t)
//          +(1./2./pow2(Nc)*(s*s+u*u)/(t*t)
//                -(pow2(Nc)-2.)/(2.*pow2(Nc))*(t*t+u*u)/(s*s)
//                -1./pow(Nc,3)*(u*u)/(s*t))*htu ;
//c     qg ->qg
      *WQ19=1./2./(pow2(Nc)-1.)*(-s/u-u/s)
                +pow2(Nc)/2./(pow2(Nc)-1.)*(s*s+u*u)/(t*t)*u/s
          +(1./2./pow2(Nc)/(pow2(Nc)-1.)*(-s/u-u/s)
                -1./(pow2(Nc)-1.)*(s*s+u*u)/(t*t))*htu ;
//c     qbg->qbg
      *WQ20 = -*WQ19;
//-1./2./(pow2(Nc)-1.)*(-s/u-u/s)
//               -pow2(Nc)/2./(pow2(Nc)-1.)*(s*s+u*u)/(t*t)*u/s
//          +(-1./2./pow2(Nc)/(pow2(Nc)-1.)*(-s/u-u/s)
//                +1./(pow2(Nc)-1.)*(s*s+u*u)/(t*t))*htu ;
//c     qg ->gq
      *WQ21 = 1./2./(pow2(Nc)-1.)*(-s/t-t/s)
                +pow2(Nc)/2./(pow2(Nc)-1.)*(s*s+t*t)/(u*u)*t/s
          +(-1./2./(pow2(Nc)-1.)*(-s/t-t/s)
                -pow2(Nc)/2./(pow2(Nc)-1.)*(s*s+t*t)/(u*u)*s/t)*htu ;
//c     qbg->gqb
      *WQ22 = -*WQ21;
//-1./2./(pow2(Nc)-1.)*(-s/t-t/s)
//                -pow2(Nc)/2./(pow2(Nc)-1.)*(s*s+t*t)/(u*u)*t/s
//          +(1./2./(pow2(Nc)-1.)*(-s/t-t/s)
//                +pow2(Nc)/2./(pow2(Nc)-1.)*(s*s+t*t)/(u*u)*s/t)*htu ;
//c     qqb->gg
      *WQ23 = -1./2./pow(Nc,3)*(t/u+u/t)-1./Nc*(t*t+u*u)/(s*s)
          +(-1./2./Nc*(t/u+u/t)+Nc/2.*(t*t+u*u)/(s*s)*u/t)*htu ;
//c     qbq->qq
      *WQ24 = -*WQ23;
//1./2./pow(Nc,3)*(t/u+u/t)+1./Nc*(t*t+u*u)/(s*s)
//         +(1./2./Nc*(t/u+u/t)-Nc/2.*(t*t+u*u)/(s*s)*u/t)*htu ;


      return 0;
};





// C-------------------------------------------------------------------
// C     spin-dependent partonic cross section
// C-------------------------------------------------------------------
      double SivSIG(double s,double t,double u,
		double HU1,double HD1,double HUB1,double HDB1,double HS1,double HSB1,double HGL1,
		double U2,double D2,double UB2,double DB2,double S2,double SB2,double GL2,
		double fu,double fub,double fd,double fdb,double fs,double fsb,double fg
){

      //real*8 HU1,HUB1,HD1,HDB1,HS1,HSB1
      double SIG1,SIG2,SIG3,SIG4;
      double WQ1,WQ2,WQ3,WQ4,WQ5,WQ6,WQ7,WQ8;
//      double WT1,WT2,WT3,WT4,WT5,WT6,WT7,WT8,
      double WQ9,WQ10,
        WQ11,WQ12,WQ13,WQ14,WQ15,WQ16,WQ17,WQ18,
        WQ19,WQ20,WQ21,WQ22,WQ23,WQ24;

     Sivpp(s,t,u,&WQ1,&WQ2,&WQ3,&WQ4,&WQ5,&WQ6,&WQ7,&WQ8,&WQ9,&WQ10,
          &WQ11,&WQ12,&WQ13,&WQ14,&WQ15,&WQ16,&WQ17,&WQ18,&WQ19,&WQ20,&WQ21,&WQ22,
          &WQ23,&WQ24);



      SIG1=WQ9*(HU1*U2*fu+HD1*D2*fd+HS1*S2*fs)
          +WQ10*(HUB1*UB2*fub+HDB1*DB2*fdb+HSB1*SB2*fsb)
          +WQ11*(HU1*UB2*(fd+fs)+HD1*DB2*(fu+fs)+HS1*SB2*(fu+fd))
          +WQ14*(HUB1*U2*(fdb+fsb)+HDB1*D2*(fub+fsb)+HSB1*S2*(fub+fdb))
          +WQ12*(HUB1*U2*(fd+fs)+HDB1*D2*(fu+fs)+HSB1*S2*(fu+fd))
          +WQ13*(HU1*UB2*(fdb+fsb)+HD1*DB2*(fub+fsb)+HS1*SB2*(fub+fdb));


      SIG2=WQ15*(HU1*UB2*fu+HD1*DB2*fd+HS1*SB2*fs)
           +WQ18*(HUB1*U2*fub+HDB1*D2*fdb+HSB1*S2*fsb)
           +WQ16*(HU1*UB2*fub+HD1*DB2*fdb+HS1*SB2*fsb)
           +WQ17*(HUB1*U2*fu+HDB1*D2*fd+HSB1*S2*fs)
           +WQ1*(HU1*(D2+S2)*fu+HD1*(U2+S2)*fd+HS1*(U2+D2)*fs)
           +WQ4*(HUB1*(DB2+SB2)*fub+HDB1*(UB2+SB2)*fdb
                +HSB1*(UB2+DB2)*fsb)
           +WQ2*(HU1*(DB2+SB2)*fu+HD1*(UB2+SB2)*fd+HS1*(UB2+DB2)*fs)
           +WQ3*(HUB1*(D2+S2)*fub+HDB1*(U2+S2)*fdb
                +HSB1*(U2+D2)*fsb);

      SIG3=WQ5*((HD1+HS1)*U2*fu+(HU1+HS1)*D2*fd
                +(HU1+HD1)*S2*fs)
            +WQ8*((HDB1+HSB1)*UB2*fub+(HUB1+HSB1)*DB2*fdb
                +(HUB1+HDB1)*SB2*fsb)
            +WQ6*((HDB1+HSB1)*U2*fu+(HUB1+HSB1)*D2*fd
                +(HUB1+HDB1)*S2*fs)
            +WQ7*((HD1+HS1)*UB2*fub+(HU1+HS1)*DB2*fdb
                +(HU1+HD1)*SB2*fsb)
            +WQ23*(HU1*UB2+HD1*DB2+HS1*SB2)*fg
            +WQ24*(HUB1*U2+HDB1*D2+HSB1*S2)*fg;

      SIG4=WQ19*GL2*(HU1*fu+HD1*fd+HS1*fs)
            +WQ20*(HUB1*fub+HDB1*fdb+HSB1*fsb)
            +WQ21*GL2*fg*(HU1+HD1+HS1)
            +WQ22*GL2*fg*(HUB1+HDB1+HSB1);

      double SIG=SIG1+SIG2+SIG3+SIG4;

      return SIG;

}




//========================================================= Approximate sivers in PP...
int sidis::approximate_sivers_pp_calculation( double xa, double xb, double zc, double y, double pt, double Q2 )
{
//  double integrand = 0.;

  double s = sidis::GetEnergy();
  double sqrts = sqrt(s);

//  double MASS = Target.GetMass();


  double Q  = sqrt( Q2 );

  double FU =0.,FD =0.,FUB =0.,FDB =0.,FS =0.,FSB =0.,FGL =0.;
  double U1 =0.,D1 =0.,UB1 =0.,DB1 =0.,S1 =0.,SB1 =0.,GL1 =0.;
  double U2 =0.,D2 =0.,UB2 =0.,DB2 =0.,S2 =0.,SB2 =0.,GL2 =0.;

  double PU1,PD1,PUB1,PDB1,PS1,PSB1,PGL1;
  double MU1,MD1,MUB1,MDB1,MS1,MSB1,MGL1;
  double DU1,DD1,DUB1,DDB1,DS1,DSB1,DGL1;

  double cross_section;

  double  sh = xa * xb * s;
  double  uh = -pt/zc * xb * sqrts * exp(y);
  double  th = -pt/zc * xa * sqrts * exp(-y);

  double coeff, coeff_siv, coeff_unp;

  double delta = 1.e-3;
  double xap = xa;
  double xam = xa;
  int count = 10;
//  bool condition;

  switch( effect ){
  case sivers_effect: // SIVERS

    switch( Beam.GetPolarization() ){
    case transverse:
      F1TDistributionFirstMoment( Beam , xa, Q2 ); //  FIRST MOMENT OF SIVERS

      U1  = Beam.partcontent.up;
      D1  = Beam.partcontent.down;
      UB1 = Beam.partcontent.anti_up;
      DB1 = Beam.partcontent.anti_down;
      S1  = Beam.partcontent.strange;
      SB1 = Beam.partcontent.anti_strange;
      GL1 = Beam.partcontent.glu;


      xap = xa*(1.+delta);
      xam = xa*(1.-delta);

      count = 0;
      while( xap > 1. ){
	delta = delta/2.;
        xap = xa*(1.+delta);
	count = count +1;
	if (count >= 10) xap = 1.;
      };

      F1TDistributionFirstMoment( Beam , xap, Q2 ); //  FIRST MOMENT OF SIVERS

      PU1  = Beam.partcontent.up;
      PD1  = Beam.partcontent.down;
      PUB1 = Beam.partcontent.anti_up;
      PDB1 = Beam.partcontent.anti_down;
      PS1  = Beam.partcontent.strange;
      PSB1 = Beam.partcontent.anti_strange;
      PGL1 = Beam.partcontent.glu;

      F1TDistributionFirstMoment( Beam , xam, Q2 ); //  FIRST MOMENT OF SIVERS

      MU1  = Beam.partcontent.up;
      MD1  = Beam.partcontent.down;
      MUB1 = Beam.partcontent.anti_up;
      MDB1 = Beam.partcontent.anti_down;
      MS1  = Beam.partcontent.strange;
      MSB1 = Beam.partcontent.anti_strange;
      MGL1 = Beam.partcontent.glu;

      DU1  = (PU1-MU1) / (2.*xa*delta);
      DD1  = (PD1-MD1) / (2.*xa*delta);
      DUB1 = (PUB1-MUB1) / (2.*xa*delta);
      DDB1 = (PDB1-MDB1) / (2.*xa*delta) ;
      DS1  = (PS1-MS1) / (2.*xa*delta);
      DSB1 = (PSB1-MSB1) / (2.*xa*delta);
      DGL1 = (PGL1-MGL1) / (2.*xa*delta);


//       U1  = (PU1+MU1)   / 2.;
//       D1  = (PD1+MD1)   / 2.;
//       UB1 = (PUB1+MUB1) / 2.;
//       DB1 = (PDB1+MDB1) / 2.;
//       S1  = (PS1+MS1)   / 2.;
//       SB1 = (PSB1+MSB1) / 2.;
//       GL1 = (PGL1+MGL1) / 2.;


      U1 =   U1 -  xa * DU1   ;
      D1 =   D1 -  xa * DD1   ;
      UB1 = UB1 -  xa * DUB1  ;
      DB1 = DB1 -  xa * DDB1  ;
      S1  = S1  -  xa * DS1   ;
      SB1 = SB1 -  xa * DSB1  ;
      GL1 = GL1 -  xa * DGL1  ;



      U1  *=  -2. * mpr; // see definition of Kang
      D1  *=  -2. * mpr;
      UB1 *=  -2. * mpr;
      DB1 *=  -2. * mpr;
      S1  *=  -2. * mpr;
      SB1 *=  -2. * mpr;
      GL1 *=  -2. * mpr;

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

      break;
    default:
	cerr << " No unp or long here!! int sidis::approximate_sivers_pp_calculation( double xa, double xb, double zc, double y, double pt, double Q2 )  " << endl;
    break;
    };

    if(!VALENCE_SIDIS){
      cross_section = SivSIG( sh,  th,  uh,
		  U1,  D1,  UB1,  DB1,  S1,  SB1,  GL1,
		  U2,  D2,  UB2,  DB2,  S2,  SB2,  GL2,
		  FU,  FUB,  FD,  FDB,  FS,  FSB,  FGL);
    } else if(VALENCE_SIDIS){ // WE USE VALENCE STANDARD PARAMETRIZATION!!!!
      cross_section = SivSIG( sh,  th,  uh,
		  U1,  D1,  2.*UB1,  2.*DB1,  S1,  SB1,  GL1,
		  U2,  D2,  2.*UB2,  2.*DB2,  S2,  SB2,  GL2,
		  FU,  FUB,  FD,  FDB,  FS,  FSB,  FGL);
    }


    switch( dependence ){
    case pt_dependence:
       coeff_siv = 1. ;
       break;
    default: // x, z dependence
       coeff_siv = 1./zc/(pow2(xa)*pow2(xb))*(-pt/(zc*uh)) ;
     break;
    }

    coeff = coeff_siv;

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
		  U1,  D1,  2.*UB1,  2.*DB1,  S1,  SB1,  GL1,
		  U2,  D2,  2.*UB2,  2.*DB2,  S2,  SB2,  GL2,
		  FU,  FUB,  FD,  FDB,  FS,  FSB,  FGL);
    }

    switch( dependence ){
    case pt_dependence:
      coeff_unp = 1. ; // MAKE IT THE SAME AS IN BOER MULDERS
       break;
    default: // x, z dependence
      coeff_unp = 1./zc/(pow2(xa)*pow2(xb)) ;
     break;
    }

    coeff = coeff_unp;

    break;
  default: //

    cout << "Error in sidis::approximate_sivers_pp_calculation( double xa, double xb, double z, double pt, double Q2 ) no effect is defined" << endl;
    return 1;

    break;
  };



  sidis::cross_section = cross_section;

  sidis::cross_section *= coeff * pow2( alphasQCD( Q ) / s );

  return 0;
}



//========================================================= Approximate sivers in PP for jet production...
int sidis::approximate_sivers_pp_jet_calculation( double xb, double y, double pt, double Q2 )
{
//  double integrand = 0.;

  double s = sidis::GetEnergy();
  double sqrts = sqrt(s);

//  double MASS = Target.GetMass();


  double Q  = sqrt( Q2 );

  double FU,FD,FUB,FDB,FS,FSB,FGL;
  double U1 =0.,D1 =0.,UB1 =0.,DB1 =0.,S1 =0.,SB1 =0.,GL1 =0.;
  double U2 =0.,D2 =0.,UB2 =0.,DB2 =0.,S2 =0.,SB2 =0.,GL2 =0.;

  double PU1,PD1,PUB1,PDB1,PS1,PSB1,PGL1;
  double MU1,MD1,MUB1,MDB1,MS1,MSB1,MGL1;
  double DU1,DD1,DUB1,DDB1,DS1,DSB1,DGL1;



//   double Umax1,Dmax1,UBmax1,DBmax1,Smax1,SBmax1,GLmax1; //unpol in 1
//   double SUmax1,SDmax1,SUBmax1,SDBmax1,SSmax1,SSBmax1,SGLmax1; //siv in 1
//   double Umin1,Dmin1,UBmin1,DBmin1,Smin1,SBmin1,GLmin1; //unpol in xa_min
//   double SUmin1,SDmin1,SUBmin1,SDBmin1,SSmin1,SSBmin1,SGLmin1; //siv in xa_min


  double cross_section;

  double SS = s;
  double UU = - pt * sqrts * exp(y);
  double TT = - pt * sqrts * exp(-y);

  double xa = -xb*UU/(xb*SS +TT);
//  double xa_min = -UU/(SS +TT);

//   cout << "xb =" << xb << endl;
//   cout << "xa =" << xa << endl;


  if (xa >= 0.9 || xb >= 0.9) {sidis::cross_section = 0.; return 0;}; // otherwise we will hit x=1 point and sometimes (1-x)^-N --> infinity

  double  sh = xa * xb * SS;
  double  uh = xb * UU;
  double  th = xa * TT;

  double coeff, coeff_siv, coeff_unp;

  double delta = 1.e-3;
//  double delta = 1.e-6;
  double xap = xa;
  double xam = xa;
  int count = 10;
//  bool condition;

  switch( effect ){
  case sivers_effect: // SIVERS

    switch( Beam.GetPolarization() ){
    case transverse:
      F1TDistributionFirstMoment( Beam , xa, Q2 ); //  FIRST MOMENT OF SIVERS

      U1  = Beam.partcontent.up;
      D1  = Beam.partcontent.down;
      UB1 = Beam.partcontent.anti_up;
      DB1 = Beam.partcontent.anti_down;
      S1  = Beam.partcontent.strange;
      SB1 = Beam.partcontent.anti_strange;
      GL1 = Beam.partcontent.glu;

      xap = xa*(1.+delta);
      xam = xa*(1.-delta);

      count = 0;
      while( xap > 1. ){
	delta = delta/2.;
        xap = xa*(1.+delta);
	count = count +1;
	if (count >= 10) xap = 1.;
      };



      F1TDistributionFirstMoment( Beam , xap, Q2 ); //  FIRST MOMENT OF SIVERS

      PU1  = Beam.partcontent.up;
      PD1  = Beam.partcontent.down;
      PUB1 = Beam.partcontent.anti_up;
      PDB1 = Beam.partcontent.anti_down;
      PS1  = Beam.partcontent.strange;
      PSB1 = Beam.partcontent.anti_strange;
      PGL1 = Beam.partcontent.glu;

      F1TDistributionFirstMoment( Beam , xam, Q2 ); //  FIRST MOMENT OF SIVERS

      MU1  = Beam.partcontent.up;
      MD1  = Beam.partcontent.down;
      MUB1 = Beam.partcontent.anti_up;
      MDB1 = Beam.partcontent.anti_down;
      MS1  = Beam.partcontent.strange;
      MSB1 = Beam.partcontent.anti_strange;
      MGL1 = Beam.partcontent.glu;

      DU1  = (PU1-MU1) / (2.*xa*delta);
      DD1  = (PD1-MD1) / (2.*xa*delta);
      DUB1 = (PUB1-MUB1) / (2.*xa*delta);
      DDB1 = (PDB1-MDB1) / (2.*xa*delta) ;
      DS1  = (PS1-MS1) / (2.*xa*delta);
      DSB1 = (PSB1-MSB1) / (2.*xa*delta);
      DGL1 = (PGL1-MGL1) / (2.*xa*delta);


// cout << DU1 << endl;
// cout << DD1 << endl;
// cout << DUB1 << endl;
// cout << DDB1 << endl;
// cout << DS1 << endl;
// cout << DSB1 << endl;
// cout << DGL1 << endl;


//Kang Gamberg derivative contribution
       U1 =  U1 -  xa * DU1   ;
       D1 =  D1 -  xa * DD1   ;
       UB1 = UB1 -  xa * DUB1  ;
       DB1 = DB1 -  xa * DDB1  ;
       S1  = S1  -  xa * DS1   ;
       SB1 = SB1 -  xa * DSB1  ;
       GL1 = GL1 -  xa * DGL1  ;


// Anikin Teryaev no derivatives:
//        U1 =  U1 -  0*xa * DU1   ;
//        D1 =  D1 -  0*xa * DD1   ;
//        UB1 = UB1 -  0*xa * DUB1  ;
//        DB1 = DB1 -  0*xa * DDB1  ;
//        S1  = S1  -  0*xa * DS1   ;
//        SB1 = SB1 -  0*xa * DSB1  ;
//        GL1 = GL1 -  0*xa * DGL1  ;


//  cout << "---------" << endl;
//  cout << "xa =" << xa << endl;
//  cout << U1 << endl;
//  cout << D1 << endl;
//  cout << UB1 << endl;
//  cout << DB1 << endl;
//  cout << S1 << endl;
//  cout << SB1 << endl;
//  cout << GL1 << endl;


// let's do partial integration

//       F1TDistributionFirstMoment( Target , 0.99, Q2 ); //  FIRST MOMENT OF SIVERS at 1...
//
//       SUmax1  = Target.partcontent.up;
//       SDmax1  = Target.partcontent.down;
//       SUBmax1 = Target.partcontent.anti_up;
//       SDBmax1 = Target.partcontent.anti_down;
//       SSmax1  = Target.partcontent.strange;
//       SSBmax1 = Target.partcontent.anti_strange;
//       SGLmax1 = Target.partcontent.glu;
//
//
//       cout << SUmax1 << endl;
// //
//       Unpolarised( Target , 1., Q2 ); //  Unpolarised at 1...
//
//       Umax1  = Target.partcontent.up;
//       Dmax1  = Target.partcontent.down;
//       UBmax1 = Target.partcontent.anti_up;
//       DBmax1 = Target.partcontent.anti_down;
//       Smax1  = Target.partcontent.strange;
//       SBmax1 = Target.partcontent.anti_strange;
//       GLmax1 = Target.partcontent.glu;
//
//
//       cout << Umax1 << endl;
//
// //
//       F1TDistributionFirstMoment( Target , xa_min, Q2 ); //  FIRST MOMENT OF SIVERS at xmin...
//
//       SUmin1  = Target.partcontent.up;
//       SDmin1  = Target.partcontent.down;
//       SUBmin1 = Target.partcontent.anti_up;
//       SDBmin1 = Target.partcontent.anti_down;
//       SSmin1  = Target.partcontent.strange;
//       SSBmin1 = Target.partcontent.anti_strange;
//       SGLmin1 = Target.partcontent.glu;
//
//       Unpolarised( Target , xa_min, Q2 ); //  Unpolarised at xmin...
//
//       Umin1  = Target.partcontent.up;
//       Dmin1  = Target.partcontent.down;
//       UBmin1 = Target.partcontent.anti_up;
//       DBmin1 = Target.partcontent.anti_down;
//       Smin1  = Target.partcontent.strange;
//       SBmin1 = Target.partcontent.anti_strange;
//       GLmin1 = Target.partcontent.glu;
//
//
// //derivative of unpolarised function!
//       Unpolarised( Target , xap, Q2 ); //  Unpolarised
//
//       PU1  = Target.partcontent.up;
//       PD1  = Target.partcontent.down;
//       PUB1 = Target.partcontent.anti_up;
//       PDB1 = Target.partcontent.anti_down;
//       PS1  = Target.partcontent.strange;
//       PSB1 = Target.partcontent.anti_strange;
//       PGL1 = Target.partcontent.glu;
//
//       Unpolarised( Target , xam, Q2 ); //  Unpolarised
//
//       MU1  = Target.partcontent.up;
//       MD1  = Target.partcontent.down;
//       MUB1 = Target.partcontent.anti_up;
//       MDB1 = Target.partcontent.anti_down;
//       MS1  = Target.partcontent.strange;
//       MSB1 = Target.partcontent.anti_strange;
//       MGL1 = Target.partcontent.glu;
//
//       DU1  = (PU1-MU1) / (2.*xa*delta);
//       DD1  = (PD1-MD1) / (2.*xa*delta);
//       DUB1 = (PUB1-MUB1) / (2.*xa*delta);
//       DDB1 = (PDB1-MDB1) / (2.*xa*delta) ;
//       DS1  = (PS1-MS1) / (2.*xa*delta);
//       DSB1 = (PSB1-MSB1) / (2.*xa*delta);
//       DGL1 = (PGL1-MGL1) / (2.*xa*delta);
//
// //   double Umax1,Dmax1,UBmax1,DBmax1,Smax1,SBmax1,GLmax1; //unpol in 1
// //   double SUmax1,SDmax1,SUBmax1,SDBmax1,SSmax1,SSBmax1,SGLmax1; //siv in 1
// //   double Umin1,Dmin1,UBmin1,DBmin1,Smin1,SBmin1,GLmin1; //unpol in xa_min
// //   double SUmin1,SDmin1,SUBmin1,SDBmin1,SSmin1,SSBmax1,SGLmin1; //siv in xa_min
//
//       double interval = 1 - xa_min;
//
//        U1 =   U1 - 0*Umax1*SUmax1/interval + 0*Umin1*SUmin1/interval + 0*xa * DU1   ;
//        D1 =   D1 -  xa * DD1   ;
//        UB1 = UB1 -  xa * DUB1  ;
//        DB1 = DB1 -  xa * DDB1  ;
//        S1  = S1  -  xa * DS1   ;
//        SB1 = SB1 -  xa * DSB1  ;
//        GL1 = GL1 -  xa * DGL1  ;

/* cout << U1 << endl;
 cout << D1 << endl;
 cout << UB1 << endl;
 cout << DB1 << endl;
 cout << S1 << endl;
 cout << SB1 << endl;
 cout << GL1 << endl; */

      U1  *=  -2. * mpr; // see definition of Kang
      D1  *=  -2. * mpr;
      UB1 *=  -2. * mpr;
      DB1 *=  -2. * mpr;
      S1  *=  -2. * mpr;
      SB1 *=  -2. * mpr;
      GL1 *=  -2. * mpr;

    Unpolarised( Target , xb, Q2 );  // 2nd hadron

      U2  = Target.partcontent.up;
      D2  = Target.partcontent.down;
      UB2 = Target.partcontent.anti_up;
      DB2 = Target.partcontent.anti_down;
      S2  = Target.partcontent.strange;
      SB2 = Target.partcontent.anti_strange;
      GL2 = Target.partcontent.glu;




      // NO fragmentation!
      FU  = 1.;
      FD  = 1.;
      FUB = 1.;
      FDB = 1.;
      FS  = 1.;
      FSB = 1.;
      FGL = 1.;

      break;
    default:
	cerr << " No unp or long here!! int sidis::approximate_sivers_pp_jet_calculation( double xb, double y, double pt, double Q2 )   " << endl;
    break;

    };

    if(!VALENCE_SIDIS){
      cross_section = SivSIG( sh,  th,  uh,
		  U1,  D1,  UB1,  DB1,  S1,  SB1,  GL1,
		  U2,  D2,  UB2,  DB2,  S2,  SB2,  GL2,
		  FU,  FUB,  FD,  FDB,  FS,  FSB,  FGL);
    } else if(VALENCE_SIDIS){ // WE USE VALENCE STANDARD PARAMETRIZATION!!!!
       cross_section = SivSIG( sh,  th,  uh,
		  U1,  D1,  2.*UB1,  2.*DB1,  S1,  SB1,  GL1,
		  U2,  D2,  2.*UB2,  2.*DB2,  S2,  SB2,  GL2,
		  FU,  FUB,  FD,  FDB,  FS,  FSB,  FGL);
    }


    switch( dependence ){
    case pt_dependence:
       coeff_siv = 1. ;
       break;
    default: // x, z dependence
       coeff_siv = 1./(xb*SS+TT)/(xb*xa)*(-pt/(uh)) ;
     break;
    }

    coeff = coeff_siv;

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


    // NO Fragmentation


      FU  = 1.;
      FD  = 1.;
      FUB = 1.;
      FDB = 1.;
      FS  = 1.;
      FSB = 1.;
      FGL = 1.;


    if(!VALENCE_SIDIS){
      cross_section = UUSIG( sh,  th,  uh,
		  U1,  D1,  UB1,  DB1,  S1,  SB1,  GL1,
		  U2,  D2,  UB2,  DB2,  S2,  SB2,  GL2,
		  FU,  FUB,  FD,  FDB,  FS,  FSB,  FGL);
    } else if(VALENCE_SIDIS){ // WE USE VALENCE STANDARD PARAMETRIZATION!!!!
      cross_section = UUSIG( sh,  th,  uh,
		  U1,  D1,  2.*UB1,  2.*DB1,  S1,  SB1,  GL1,
		  U2,  D2,  2.*UB2,  2.*DB2,  S2,  SB2,  GL2,
		  FU,  FUB,  FD,  FDB,  FS,  FSB,  FGL);
    }

    switch( dependence ){
    case pt_dependence:
      coeff_unp = 1. ; // MAKE IT THE SAME AS IN BOER MULDERS
       break;
    default: // x, z dependence
      coeff_unp = 1./(xb*SS+TT)/(xb*xa) ;
     break;
    }

    coeff = coeff_unp;

    break;
  default: //

    cout << "Error in sidis::approximate_sivers_pp_calculation( double xa, double xb, double z, double pt, double Q2 ) no effect is defined" << endl;
    return 1;

    break;
  };



  sidis::cross_section = cross_section;

  sidis::cross_section *= coeff * pow2( alphasQCD( Q ) ) / s ;

  return 0;
}


// Let us integrate jet cross section using another method

struct PAR_FOR_JET {
  double y, pt, Q2;
};


// we will use GSL integration
double pp_jet (double xb, void * params) {
  sidis& SIDIS = *sidis::Instance();

  struct PAR_FOR_JET * pars = (struct PAR_FOR_JET *) params;

  double y = pars->y;
  double pt = pars->pt;
  double Q2 = pars->Q2;

  SIDIS.approximate_sivers_pp_jet_calculation( xb, y, pt, Q2 );

  return SIDIS.GetCrossSection();
};




// THE SAME USING GAUS LEGENDRE TEST TEST TEST
// we will use  GAUS LEGENDRE gaus_legendre.c TEST TEST TEST
double pp_jet_gauss (double xb, void * params) {
  sidis& SIDIS = *sidis::Instance();

//  double s = SIDIS.GetEnergy();
//  double sqrts = sqrt(s);

  struct PAR_FOR_JET * pars = (struct PAR_FOR_JET *) params;

  double y = pars->y;
  double pt = pars->pt;
  double Q2 = pars->Q2;


  SIDIS.approximate_sivers_pp_jet_calculation( xb, y, pt, Q2 );

  return SIDIS.GetCrossSection();
};

#ifndef FABS
	#define FABS(a) ((a)>=0?(a):-(a))
#endif

/* n = 12 */
static double W[12] = {0.1012285,0.2223810,0.3137067,0.3623838,0.0271525,
            0.0622535,0.0951585,0.1246290,0.1495960,0.1691565,
            0.1826034,0.1894506 };
static double X[12] = {0.9602899,0.7966665,0.5255324,0.1834346,0.9894009,
            0.9445750,0.8656312,0.7554044,0.6178762,0.4580168,
            0.2816036,0.0950125};



// let me rewrite the same subroutine used by Zhongbo...
double gaussp(double (*f)(double,void*), void* data, double A, double B, double ACC){

      double GAUSSP,CONST,DELTA,AA,BB,Y;
      double C1,C2,S8,S16,U,SPLIT;
//      int I;

      CONST = 1.E-5;
      SPLIT = 2.;


      DELTA=CONST*FABS(A-B);

      GAUSSP=0.0;

      AA=A;
      Y=B-AA;
      BB=AA+Y;

      do{

      //if(FABS(Y) <= DELTA) return GAUSSP;

      C1=0.5*(AA+BB);
      C2=C1-AA;
      S8=0.0;
      S16=0.0;

      for(int I=0;I<=3;I++){
      //cout<< I << " " << X[I] << " " << W[I] << endl;
      U=X[I]*C2;
      S8=S8+W[I]*((*f)(C1+U,data)+(*f)(C1-U,data));
      }

      for(int I=4;I<=11;I++){
      //cout<< I << " " << X[I] << " " << W[I] << endl;
      U=X[I]*C2;
      S16=S16+W[I]*((*f)(C1+U,data)+(*f)(C1-U,data));
      }


      S8=S8*C2;
      S16=S16*C2;

      if( FABS(S16-S8) > ACC*FABS(S16) ){
        Y=Y/SPLIT;
        if(FABS(Y) > DELTA){
          BB = AA +Y;
        } else {
          GAUSSP=GAUSSP+S16;
          AA = BB;
          BB = AA +Y;
        };
      } else {
          GAUSSP=GAUSSP+S16;
          AA = BB;
          Y  = B - AA;
          BB = AA + Y;
      };

      } while (FABS(Y) >= DELTA);

      return GAUSSP;
}


// calculates numerator/denominator for jet depends on settings of sidis.effect
double jet_numerator_denominator_gauss(double xF)
{
  sidis& SIDIS = *sidis::Instance();

  double result;
  struct PAR_FOR_JET pars;
  double s = SIDIS.GetEnergy();
  double sqrts = sqrt(s);

  double y = SIDIS.Average.y; // this is rapidity!

  double pt = xF * sqrts/( exp(y) - exp(-y) );

  double Q2 = pt*pt;

  double SS = s;
  double TT = -sqrts * pt * exp(-y);
  double UU = -sqrts * pt * exp(y);


  pars.y  = y;
  pars.pt = pt;
  pars.Q2 = Q2;


  int key = 1024;
//  int key = 256;

  result = gauss_legendre (key, pp_jet_gauss, &pars, -TT/(SS+UU), 0.95);
//  result =   gaussp (pp_jet_gauss, &pars, -TT/(SS+UU), 0.95, 1.e-2);



  return  result; //
}




// THE SAME USING GAUS LEGENDRE CALCULATE DIRECT PHOTONS!!!! TEST TEST TEST
// we will use  GAUS LEGENDRE gaus_legendre.c TEST TEST TEST
double pp_photon_gauss (double xb, void * params) {
  sidis& SIDIS = *sidis::Instance();

//  double s = SIDIS.GetEnergy();
//  double sqrts = sqrt(s);

  struct PAR_FOR_JET * pars = (struct PAR_FOR_JET *) params;

  double y = pars->y;
  double pt = pars->pt;
  double Q2 = pars->Q2;


  SIDIS.approximate_sivers_pp_photon_calculation( xb, y, pt, Q2 ); // direct photon

  return SIDIS.GetCrossSection();
};


// calculates numerator/denominator for jet depends on settings of sidis.effect DIRECT PHOTON
double photon_numerator_denominator_gauss(double xF)
{
  sidis& SIDIS = *sidis::Instance();

  double result;
  struct PAR_FOR_JET pars;
  double s = SIDIS.GetEnergy();
  double sqrts = sqrt(s);

  double y = SIDIS.Average.y; // this is rapidity!

  double pt = xF * sqrts/( exp(y) - exp(-y) );

  double Q2 = pt*pt;

  double SS = s;
  double TT = -sqrts * pt * exp(-y);
  double UU = -sqrts * pt * exp(y);


  pars.y  = y;
  pars.pt = pt;
  pars.Q2 = Q2;


  int key = 1024;
//  int key = 256;

  result = gauss_legendre (key, pp_photon_gauss, &pars, -TT/(SS+UU), 0.99);




  return  result; //
}



//========================================================= Approximate sivers in PP for direct photon production...
int sidis::approximate_sivers_pp_photon_calculation( double xb, double y, double pt, double Q2 )
{
//  double integrand = 0.;

  double s = sidis::GetEnergy();
  double sqrts = sqrt(s);

//  double MASS = Target.GetMass();


  double Q  = sqrt( Q2 );

//  double FU,FD,FUB,FDB,FS,FSB,FGL;

  double U1 =0.,D1 =0.,UB1 =0.,DB1 =0.,S1 =0.,SB1 =0.,GL1 =0.;
  double U2 =0.,D2 =0.,UB2 =0.,DB2 =0.,S2 =0.,SB2 =0.,GL2 =0.;

  double PU1,PD1,PUB1,PDB1,PS1,PSB1,PGL1;
  double MU1,MD1,MUB1,MDB1,MS1,MSB1,MGL1;
  double DU1,DD1,DUB1,DDB1,DS1,DSB1,DGL1;



//   double Umax1,Dmax1,UBmax1,DBmax1,Smax1,SBmax1,GLmax1; //unpol in 1
//   double SUmax1,SDmax1,SUBmax1,SDBmax1,SSmax1,SSBmax1,SGLmax1; //siv in 1
//   double Umin1,Dmin1,UBmin1,DBmin1,Smin1,SBmin1,GLmin1; //unpol in xa_min
//   double SUmin1,SDmin1,SUBmin1,SDBmin1,SSmin1,SSBmin1,SGLmin1; //siv in xa_min


  double cross_section;

  double SS = s;
  double UU = - pt * sqrts * exp(y);
  double TT = - pt * sqrts * exp(-y);

  double xa = -xb*UU/(xb*SS +TT);
//  double xa_min = -UU/(SS +TT);

//   cout << "xb =" << xb << endl;
//   cout << "xa =" << xa << endl;


  if (xa >= 0.9 || xb >= 0.99) {sidis::cross_section = 0.; return 0;}; // otherwise we will hit x=1 point and sometimes (1-x)^-N --> infinity

  double  sh = xa * xb * SS;
  double  uh = xb * UU;
  double  th = xa * TT;

  double coeff, coeff_siv, coeff_unp;

  double delta = 1.e-3;
//  double delta = 1.e-6;
  double xap = xa;
  double xam = xa;
  int count = 10;
//  bool condition;

  switch( effect ){
  case sivers_effect: // SIVERS

    switch( Target.GetPolarization() ){
    case transverse:
      F1TDistributionFirstMoment( Target , xa, Q2 ); //  FIRST MOMENT OF SIVERS

      U1  = Target.partcontent.up;
      D1  = Target.partcontent.down;
      UB1 = Target.partcontent.anti_up;
      DB1 = Target.partcontent.anti_down;
      S1  = Target.partcontent.strange;
      SB1 = Target.partcontent.anti_strange;
      GL1 = Target.partcontent.glu;

      xap = xa*(1.+delta);
      xam = xa*(1.-delta);

      count = 0;
      while( xap > 1. ){
	delta = delta/2.;
        xap = xa*(1.+delta);
	count = count +1;
	if (count >= 10) xap = 1.;
      };



      F1TDistributionFirstMoment( Target , xap, Q2 ); //  FIRST MOMENT OF SIVERS

      PU1  = Target.partcontent.up;
      PD1  = Target.partcontent.down;
      PUB1 = Target.partcontent.anti_up;
      PDB1 = Target.partcontent.anti_down;
      PS1  = Target.partcontent.strange;
      PSB1 = Target.partcontent.anti_strange;
      PGL1 = Target.partcontent.glu;

      F1TDistributionFirstMoment( Target , xam, Q2 ); //  FIRST MOMENT OF SIVERS

      MU1  = Target.partcontent.up;
      MD1  = Target.partcontent.down;
      MUB1 = Target.partcontent.anti_up;
      MDB1 = Target.partcontent.anti_down;
      MS1  = Target.partcontent.strange;
      MSB1 = Target.partcontent.anti_strange;
      MGL1 = Target.partcontent.glu;

      DU1  = (PU1-MU1) / (2.*xa*delta);
      DD1  = (PD1-MD1) / (2.*xa*delta);
      DUB1 = (PUB1-MUB1) / (2.*xa*delta);
      DDB1 = (PDB1-MDB1) / (2.*xa*delta) ;
      DS1  = (PS1-MS1) / (2.*xa*delta);
      DSB1 = (PSB1-MSB1) / (2.*xa*delta);
      DGL1 = (PGL1-MGL1) / (2.*xa*delta);


// Kang Gamberg -> derivative terms
       U1 =  U1 -  xa * DU1    ;
       D1 =  D1 -  xa * DD1    ;
       UB1 = UB1 -  xa * DUB1   ;
       DB1 = DB1 -  xa * DDB1   ;
       S1  = S1  -  xa * DS1    ;
       SB1 = SB1 -  xa * DSB1   ;
       GL1 = GL1 -  xa * DGL1 *0 ;

// Anikin Teryaev -> No derivative terms
//        U1 =  U1      ;
//        D1 =  D1     ;
//        UB1 = UB1     ;
//        DB1 = DB1     ;
//        S1  = S1       ;
//        SB1 = SB1     ;
//        GL1 = GL1   ;

      U1  *=  -2. * mpr; // see definition of Kang
      D1  *=  -2. * mpr;
      UB1 *=  -2. * mpr;
      DB1 *=  -2. * mpr;
      S1  *=  -2. * mpr;
      SB1 *=  -2. * mpr;
      GL1 *=  -2. * mpr;

    Unpolarised( Target , xb, Q2 );  // 2nd hadron

      U2  = Target.partcontent.up;
      D2  = Target.partcontent.down;
      UB2 = Target.partcontent.anti_up;
      DB2 = Target.partcontent.anti_down;
      S2  = Target.partcontent.strange;
      SB2 = Target.partcontent.anti_strange;
      GL2 = Target.partcontent.glu;





      break;
    default:
	cerr << " No unp or long here!! int sidis::approximate_sivers_pp_photon_calculation( double xb, double y, double pt, double Q2 ) " << endl;
    break;

    };

    if(!VALENCE_SIDIS){
      cross_section = SivSIG_PHOTON( sh,  th,  uh,
		  U1,  D1,  UB1,  DB1,  S1,  SB1,  GL1,
		  U2,  D2,  UB2,  DB2,  S2,  SB2,  GL2 );
    } else if(VALENCE_SIDIS){ // WE USE VALENCE STANDARD PARAMETRIZATION!!!!
       cross_section = SivSIG_PHOTON( sh,  th,  uh,
		  U1,  D1,  2.*UB1,  2.*DB1,  S1,  SB1,  GL1,
		  U2,  D2,  2.*UB2,  2.*DB2,  S2,  SB2,  GL2);
    }


    switch( dependence ){
    case pt_dependence:
       coeff_siv = 1. ;
       break;
    default: // x, z dependence
       coeff_siv = 1./(xb*SS+TT)/(xb*xa)*(-pt/(uh)) ;
     break;
    }

    coeff = coeff_siv;

    break;

  case no_effect: // UNPOLARISED


    Unpolarised( Target , xa, Q2 );  // 1st hadron

      U1  = Target.partcontent.up;
      D1  = Target.partcontent.down;
      UB1 = Target.partcontent.anti_up;
      DB1 = Target.partcontent.anti_down;
      S1  = Target.partcontent.strange;
      SB1 = Target.partcontent.anti_strange;
      GL1 = Target.partcontent.glu;

    Unpolarised( Target , xb, Q2 );  // 2nd hadron

      U2  = Target.partcontent.up;
      D2  = Target.partcontent.down;
      UB2 = Target.partcontent.anti_up;
      DB2 = Target.partcontent.anti_down;
      S2  = Target.partcontent.strange;
      SB2 = Target.partcontent.anti_strange;
      GL2 = Target.partcontent.glu;




    if(!VALENCE_SIDIS){
      cross_section = UUSIG_PHOTON( sh,  th,  uh,
		  U1,  D1,  UB1,  DB1,  S1,  SB1,  GL1,
		  U2,  D2,  UB2,  DB2,  S2,  SB2,  GL2);
    } else if(VALENCE_SIDIS){ // WE USE VALENCE STANDARD PARAMETRIZATION!!!!
      cross_section = UUSIG_PHOTON( sh,  th,  uh,
		  U1,  D1,  2.*UB1,  2.*DB1,  S1,  SB1,  GL1,
		  U2,  D2,  2.*UB2,  2.*DB2,  S2,  SB2,  GL2);
    }

    switch( dependence ){
    case pt_dependence:
      coeff_unp = 1. ; // MAKE IT THE SAME AS IN BOER MULDERS
       break;
    default: // x, z dependence
      coeff_unp = 1./(xb*SS+TT)/(xb*xa) ;
     break;
    }

    coeff = coeff_unp;

    break;
  default: //

    cout << "Error in sidis::approximate_sivers_pp_photon_calculation( double xa, double xb, double z, double pt, double Q2 ) no effect is defined" << endl;
    return 1;

    break;
  };



  sidis::cross_section = cross_section;

  sidis::cross_section *= coeff * pow2( alphasQCD( Q ) ) * alpha_em/ s ;

  return 0;
}


//C-------------------------------------------------------------------

      int Sivpp_PHOTON(double s,double t,double u,double* WQ1,double* WQ1b,double* WQ2,double* WQ3,
	  double* WQ4){
      double Nc;

      Nc  = 3.;



//c     qg ->photon+q
      *WQ1=-Nc/(pow2(Nc)-1.)*(-t/s-s/t);
//c     qbg->photon+qb
      *WQ1b=Nc/(pow2(Nc)-1.)*(-t/s-s/t);
//c     gq ->photon+q
      *WQ2=0.;
//c     qqb->photon+g
      *WQ3=1./(pow2(Nc))*(u/t+t/u);
//c     qbq->photon+g
      *WQ4=-1./(pow2(Nc))*(u/t+t/u);




      return 0;
};


// C-------------------------------------------------------------------
// C     spin-dependent partonic cross section DIRECT PHOTON
// C-------------------------------------------------------------------
      double SivSIG_PHOTON(double s,double t,double u,
		double HU1,double HD1,double HUB1,double HDB1,double HS1,double HSB1,double HGL1,
		double U2,double D2,double UB2,double DB2,double S2,double SB2,double GL2
){


      double SIG;
      double WQ1,WQ1b,WQ2,WQ3,WQ4;

      double eu2 = 4./9.;
      double ed2 = 1./9.;
      double es2 = 1./9.;

     Sivpp_PHOTON(s,t,u,&WQ1,&WQ1b,&WQ2,&WQ3,&WQ4);



      SIG=WQ1*(eu2*HU1+ed2*HD1+es2*HS1)*GL2
          +WQ1b*(eu2*HUB1+ed2*HDB1+es2*HSB1)*GL2
          +WQ3*(eu2*HU1*UB2+ed2*HD1*DB2+es2*HS1*SB2)
          +WQ4*(eu2*HUB1*U2+ed2*HDB1*D2+es2*HSB1*S2);

      return SIG;

}



// C-------------------------------------------------------------------

int Upp_PHOTON(double s,double t,double u, double* WQ1,double* WQ2,double* WQ3,
		double* WQ4){

      double Nc=3.;


//c     qg ->photon+q
      *WQ1=1./Nc*(-t/s-s/t);
//c     gq ->photon+q
      *WQ2=1./Nc*(-u/s-s/u);
//c     qqb->photon+g
      *WQ3=(pow2(Nc)-1.)/(pow2(Nc))*(u/t+t/u);
//c     qbq->photon+g
      *WQ4=(pow2(Nc)-1.)/(pow2(Nc))*(u/t+t/u);



      return 0;
};

// C-------------------------------------------------------------------
// C     unpolarized partonic cross section PHOTON
// C-------------------------------------------------------------------

      double UUSIG_PHOTON(double s,double t,double u,
		double U1,double D1,double UB1,double DB1,double S1,double SB1,double GL1,
		double U2,double D2,double UB2,double DB2,double S2,double SB2,double GL2
){


      double SIG ;
      double WQ1,WQ2,WQ3,WQ4;

      double eu2 = 4./9.;
      double ed2 = 1./9.;
      double es2 = 1./9.;

      Upp_PHOTON(s,t,u,&WQ1,&WQ2,&WQ3,&WQ4);

      SIG = WQ1*(eu2*U1+ed2*D1+eu2*UB1+ed2*DB1+es2*S1+es2*SB1)*GL2
          +WQ2*GL1*(eu2*U2+ed2*D2+eu2*UB2+ed2*DB2+es2*S2+es2*SB2)
          +WQ3*(eu2*U1*UB2+ed2*D1*DB2+es2*S1*SB2)
          +WQ4*(eu2*UB1*U2+ed2*DB1*D2+es2*SB1*S2);


      return SIG;
};
