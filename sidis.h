/* /home/prokudin/cpp/dy/sidis.h my headers for SIDIS functions */
/* Alexei Prokudin 31/01/2009        */

#ifndef __SIDIS_H__
#define __SIDIS_H__

#include <cfortran.h>
#include <LHAPDF.h>  // THIS IS CHANGED WITH RESPECT TO MAC TO:
#include "TPaveStats.h"
#include "TText.h"



#ifndef  __HADRON_H_
#include "hadron.h"
#endif

#ifndef  __LEPTON_H_
#include "lepton.h"
#endif

#ifndef  __PARTON_H_
#include "parton.h"
#endif

#ifndef  __PARAMETERS_H_
#include "parameters.h"
#endif

#ifndef  __CONSTANTS_H_
#include "constants.h"
#endif

#ifndef  __COMMON_SIDIS_H_
#include "common_sidis.h"
#endif

#ifndef  __AVERAGE_SIDIS_H_
#include "average_sidis.h"
#endif


// TWIST-3 X SECTION...
#ifndef  __TWIST3_H_
#include "twist3.h"
#endif

// Hoppet
#ifndef __HOPPET_V1__
#include "hoppet_v1.h"
#endif


#ifndef  __EVOLUTION_COMMON_
#include "evolution_common.h"
#endif


#include <TComplex.h>
#include <TCanvas.h>
#include <TLorentzVector.h>

using namespace std;
using namespace sidis_namespace;


/* *    POLARIZED RADIATIVELY GENERATED LO AND NLO PARTON DENSITIES    * */
/* *                                                                   * */
/* *         M. GLUCK, E. REYA, M. STRATMANN AND W. VOGELSANG,         * */
/* *              Phys.Rev.D63:094005,2001,    hep-ph/0011215          * */
/* PROTOCCALLSFSUB11(PARPOL,parpol,INT,DOUBLE,DOUBLE,PDOUBLE,PDOUBLE,PDOUBLE,PDOUBLE,PDOUBLE,PDOUBLE,PDOUBLE,PDOUBLE) */
/* #define  PARPOL(ISET, X, Q2, U, D, UB, DB, ST, GL, G1P, G1N) \ */
/* CCALLSFSUB11(PARPOL,parpol,INT,DOUBLE,DOUBLE,PDOUBLE,PDOUBLE,PDOUBLE,PDOUBLE,PDOUBLE,PDOUBLE,PDOUBLE,PDOUBLE,ISET, X, Q2, U, D, UB, DB, ST, GL, G1P, G1N) */

/* ********************************************************************* */
/* *                                                                   * */
/* *                    TRANSVERSITY DENSITIES                         * */
/* *                                                                   * */
/* *   INPUT:   ISET = number of the parton set :                      * */
/* *              ISET = 1  NEXT-TO-LEADING ORDER  (MS-bar)            *  */
/* *                        (DATA FILE 'transmaxnlo.grid' UNIT=11, TO  * */
/* *                         BE DEFINED BY THE USER )                  *  */
/* *              ISET = 2  LEADING ORDER                              *  */
/* *                        (DATA FILE 'transmaxlo_new.grid' UNIT=22   * */
/* *                                                                   * */
/* *            X  = Bjorken-x       (between  1.E-4  and  1)          * */
/* *            Q2 = scale in GeV**2 (between  0.8  and   1.E6)        * */
/* *             (for values outside the allowed range the program     * */
/* *              writes a warning and extrapolates to the x and       * */
/* *              Q2 values requested)                                 * */
/* *                                                                   * */
/* *   OUTPUT:  U = x * DELTA u                                        *  */
/* *            D = x * DELTA d                                        * */
/* *            UB = x * UBAR                                          * */
/* *            DB = x * DBAR                                          *    */
/* *            ST = x * DELTA STRANGE = x * DELTA STRANGE(BAR)        *      */
/* *                                                                   * */
/* *          (  For the parton distributions always x times           * */
/* *                   the distribution is returned   )                * */
/* *                                                                   * */
/* *                                                                   * */
/* *   COMMON:  The main program or the calling routine has to have    * */
/* *            a common block  COMMON / INTINI / IINI , and  IINI     * */
/* *            has always to be zero when PARPOL is called for the    * */
/* *            first time or when 'ISET' has been changed.            * */
/* *                                                                   * */
/* ********************************************************************* */
/* PROTOCCALLSFSUB8(PARPOLT,parpolt,INT,DOUBLE,DOUBLE,PDOUBLE,PDOUBLE,PDOUBLE,PDOUBLE,PDOUBLE) */
/* #define  PARPOLT(ISET, X, Q2, U, D, UB, DB, ST) \ */
/* CCALLSFSUB8(PARPOLT,parpolt,INT,DOUBLE,DOUBLE,PDOUBLE,PDOUBLE,PDOUBLE,PDOUBLE,PDOUBLE,ISET, X, Q2, U, D, UB, DB, ST) */

/* C***************************************************************C */
/* C                                                               C */
/* C    PION DISTRIBUTION FUNCTIONS EXTRACTED FROM DRELL-YAN AND   C */
/* C    PROMPT PHOTON PION BEAM DATA, USING (NUCLEAR TARGET        C */
/* C    CORRECTED) HMRS(B) DISTRIBUTIONS FOR THE NUCLEON.          C */
/* C    MAIN UNCERTAINTY IS IN SEA, HENCE FOLLOWING THREE SETS...  C */
/* C    MODE=2 IS THE "BEST FIT" SET.                              C */
/* C                                                               C */
/* C    MODE 1 CORRESPONDS TO A 10% MOMENTUM SEA    (0.11213)      C */
/* C    MODE 2 CORRESPONDS TO A 15% MOMENTUM SEA    (0.16119)      C */
/* C    MODE 3 CORRESPONDS TO A 20% MOMENTUM SEA    (0.23785)      C */
/* C                                                               C */
/* C    (NUMBER IN BRACKETS IS LAST ENTRY IN FIRST ROW OF GRID)    C */
/* C                                                               C */
/* C                         -*-                                   C */
/* C                                                               C */
/* C    REFERENCE: A.D. MARTIN. R.G. ROBERTS. W.J. STIRLING        C */
/* C    AND P.J. SUTTON, PHYS. REV. D45 (1992) 2349.               C */
/* C                         -*-                                   C */
/* C                                                               C */
/* C    FOR THESE SETS....                                         C */
/* C    * LAMBDA(MSBAR,NF=4) = 190 MEV                             C */
/* C    * UV=DV, STR=USEA=DSEA, BTM=0                              C */
/* C                                                               C */
/* C    NOTE THAT X TIMES THE PARTON DISTRIBUTION FUNCTION         C */
/* C    IS RETURNED I.E. G(X) = GLU/X ETC, AND THAT "SEA"          C */
/* C    IS THE QUARK SEA / 6   I.E. UBAR(X)=DBAR(X)                C */
/* C    = SEA/X FOR A PION. IF IN DOUBT, CHECK THE                 C */
/* C    MOMENTUM SUM RULE! NOTE ALSO THAT SCALE=Q IN GEV           C */
/* C                                                               C */
/* C                         -*-                                   C */
/* C                                                               C */
/* C     THE RANGE OF APPLICABILITY IS CURRENTLY:                  C */
/* C     10**-5 < X < 1  AND  5 < Q**2 < 1.31 * 10**6              C */
/* C     HIGHER Q**2 VALUES CAN BE SUPPLIED ON REQUEST             C */
/* C     - PROBLEMS, COMMENTS ETC TO w.j.stirling@durham.ac.uk     C */
/* C                                                               C */
/* C                                                               C */
/* C***************************************************************C */
/* ********************************************************************* */
PROTOCCALLSFSUB10(PION,pion,DOUBLE,DOUBLE,INT,PDOUBLE,PDOUBLE,PDOUBLE,PDOUBLE,PDOUBLE,PDOUBLE,PDOUBLE)
#define  PION(X,SCALE,MODE,UPV,DNV,SEA,STR,CHM,BOT,GLU) \
CCALLSFSUB10(PION,pion,DOUBLE,DOUBLE,INT,PDOUBLE,PDOUBLE,PDOUBLE,PDOUBLE,PDOUBLE,PDOUBLE,PDOUBLE,X,SCALE,MODE,UPV,DNV,SEA,STR,CHM,BOT,GLU)
PROTOCCALLSFSUB9(GRVPILO,grvpilo,DOUBLE,DOUBLE,PDOUBLE,PDOUBLE,PDOUBLE,PDOUBLE,PDOUBLE,PDOUBLE,PDOUBLE)
#define  GRVPILO(X,Q,ZUV,ZDV,ZUDB,ZSB,ZCB,ZBB,ZGL) \
CCALLSFSUB9(GRVPILO,grvpilo,PDOUBLE,DOUBLE,DOUBLE,PDOUBLE,PDOUBLE,PDOUBLE,PDOUBLE,PDOUBLE,PDOUBLE,X,Q,ZUV,ZDV,ZUDB,ZSB,ZCB,ZBB,ZGL)



/* * */
/*       SUBROUTINE GRV98PA (ISET, X, Q2, UV, DV, US, DS, SS, GL) */
/* ********************************************************************* */
/* *                                                                   * */
/* *   THE PARTON ROUTINE.                                             * */
/* *                                     __                            * */
/* *   INPUT:   ISET =  1 (LO),  2 (NLO, MS), or  3 (NLO, DIS)         * */
/* *            X  =  Bjorken-x        (between  1.E-9 and 1.)         * */
/* *            Q2 =  scale in GeV**2  (between  0.8 and 1.E6)         * */
/* *                                                                   * */
/* *   OUTPUT:  UV = u - u(bar),  DV = d - d(bar),  US = u(bar),       * */
/* *            DS = d(bar),  SS = s = s(bar),  GL = gluon.            * */
/* *            Always x times the distribution is returned.           * */
/* *                                                                   * */
/* *   COMMON:  The main program or the calling routine has to have    * */
/* *            a common block  COMMON / INTINIP / IINIP , and the     * */
/* *            integer variable  IINIP  has always to be zero when    * */
/* *            GRV98PA is called for the first time or when  ISET     * */
/* *            has been changed.                                      * */
/* *                                                                   * */
/* *   GRIDS:   1. grv98lo.grid, 2. grv98nlm.grid, 3. grv98nld.grid,   * */
/* *            (1+1809 lines with 6 columns, 4 significant figures)   * */
/* *                                                                   * */
/* *******************************************************i************* */
/* * */
PROTOCCALLSFSUB9(GRV98PA,grv98pa,INT,DOUBLE,DOUBLE,PDOUBLE,PDOUBLE,PDOUBLE,PDOUBLE,PDOUBLE,PDOUBLE)
#define GRV98PA(ISET, x, Q2, UV, DV, US, DS, SS, GL) \
CCALLSFSUB9(GRV98PA,grv98pa,INT,DOUBLE,DOUBLE,PDOUBLE,PDOUBLE,PDOUBLE,PDOUBLE,PDOUBLE,PDOUBLE,ISET, x, Q2, UV, DV, US, DS, SS, GL)


// FRAGMENTATION LIB
PROTOCCALLSFSUB9(DLIB,dlib,DOUBLE,DOUBLE,PVOID,INT,INT,INT,INT,INT,INT)
#define DLIB(Z,Q2,DFF,FFSET,FFORDER,IHADRON,ICHARGE,ICP,IPI) \
  CCALLSFSUB9(MYDLIB,dlib,DOUBLE,DOUBLE,PVOID,INT,INT,INT,INT,INT,INT,\
	      Z,Q2,DFF,FFSET,FFORDER,IHADRON,ICHARGE,ICP,IPI)


/* *    POLARIZED RADIATIVELY GENERATED LO AND NLO PARTON DENSITIES    * */
/* *                                                                   * */
/* *         M. GLUCK, E. REYA, M. STRATMANN AND W. VOGELSANG,         * */
/* *              Phys.Rev.D63:094005,2001,    hep-ph/0011215          * */
PROTOCCALLSFSUB11(PARPOL,parpol,INT,DOUBLE,DOUBLE,PDOUBLE,PDOUBLE,PDOUBLE,PDOUBLE,PDOUBLE,PDOUBLE,PDOUBLE,PDOUBLE)
#define  PARPOL(ISET, X, Q2, U, D, UB, DB, ST, GL, G1P, G1N) \
CCALLSFSUB11(PARPOL,parpol,INT,DOUBLE,DOUBLE,PDOUBLE,PDOUBLE,PDOUBLE,PDOUBLE,PDOUBLE,PDOUBLE,PDOUBLE,PDOUBLE,ISET, X, Q2, U, D, UB, DB, ST, GL, G1P, G1N)


/* CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC */
/* C                                                                            C */
/* C       ------ DSSV POLARIZED PARTON DISTRIBUTIONS ------                    C */
/* C                                                                            C */
/* C       REFERENCES:                                                          C */
/* C           "GLOBAL ANALYSIS OF HELICITY PARTON DENSITIES                    C */
/* C            AND THEIR UNCERTAINTIES",                                       C */
/* C            arXiv:0804:0422 [hep-ph]                                        C */
/* C           "EXTRACTION OF SPIN-DEPENDENT PARTON DENSITIES                   C */
/* C            AND THEIR UNCERTAINTIES",                                       C */
/* C            arXiv:0904.3821 [hep-ph]                                        C  */
/* C                                                                            C */
/* C       INPUT:                                                               C */
/* C           X     :  BJORKEN-X BETWEEN 10**(-5)< X < 1.0                     C */
/* C           Q2    :  SCALE**2 IN GeV**2 BETWEEN  1 < Q**2 < 10**5 GeV**2     C */
/* C                                                                            C */
/* C       OUTPUT:                                                              C */
/* C           DUV   :  X * U VALENCE DISTRIBUTION                              C */
/* C           DDV   :  X * D VALENCE DISTRIBUTION                              C */
/* C           DUBAR :  X * UBAR DISTRIBUTION                                   C */
/* C           DDBAR :  X * DBAR DISTRIBUTION                                   C */
/* C           DSTR  :  X * STRANGE DISTRIBUTION                                C */
/* C           DGLU  :  X * GLUON DISTRIBUTION                                  C */
/* C                                                                            C */
/* C       IMPORTANT:                                                           C */
/* C           ALWAYS X*DISTRIBUTION IS RETURNED !!!                            C */
/* C           ALL PDFs ARE IN THE MSbar SCHEME                                 C */
/* C                                                                            C */
/* C           BEFORE CALLING THE SUBROUTINE 'DSSVFIT' FOR THE FIRST TIME,      C */
/* C           THE SUBROUTINE 'DSSVINI' MUST BE CALLED (ONLY ONCE) TO           C */
/* C           INITIALIZE THE GRIDS !!                                          C */
/* C                                                                            C */
/* C       HESSIAN EIGENVECTOR SETS FOR DELTA CHI**2 = 1 UNCERTAINTIES:         C */
/* C           38 EIGENVECTOR SETS AVAILABLE TO COMPUTE DELTA CHI**2=1 UNCERT.  C */
/* C           SETS ARE INITIALIZED BY CALLING THE SUBROUTINE DSSVINI           C */
/* C           CALL DSSVINI(0) : OPTIMUM/BEST FIT OF DSSV                       C */
/* C           CALL DSSVINI(-19 .... 19) : EIGENVECTOR SETS WITH STEP IN        C */
/* C                                       MINUS OR PLUS "i" DIRECTION          C */
/* C           TO COMPUTE UNCERTAINTIES OF AN OBSERVABLE O USE                  C */
/* C           EQUATION (11) IN arXiv:0904.3821 [hep-ph]                        C */
/* C                                                                            C */
/* C       IN CASE OF PROBLEMS, DOUBTS, ETC, PLEASE E-MAIL US:                  C */
/* C           D. de Florian  deflo@df.uba.ar                                   C */
/* C           R. Sassot      sassot@df.uba.ar                                  C */
/* C           M. Stratmann   marco@ribf.riken.jp                               C */
/* C           W. Vogelsang   vogelsan@quark.phy.bnl.gov                        C */
/* C                                                                            C */
/* CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC */
PROTOCCALLSFSUB8(DSSVFIT,dssvfit,DOUBLE,DOUBLE,PDOUBLE,PDOUBLE,PDOUBLE,PDOUBLE,PDOUBLE,PDOUBLE)
#define  DSSVFIT(X,Q2,DUV,DDV,DUBAR,DDBAR,DSTR,DGLU) \
  CCALLSFSUB8(DSSVFIT,dssvfit,DOUBLE,DOUBLE,PDOUBLE,PDOUBLE,PDOUBLE,PDOUBLE,PDOUBLE,PDOUBLE,X,Q2,DUV,DDV,DUBAR,DDBAR,DSTR,DGLU)
//INITIALIZATION OF THE DENSITY...DSSVINI(ISET)
PROTOCCALLSFSUB1(DSSVINI,dssvini,INT)
#define DSSVINI(ISET) \
  CCALLSFSUB1(DSSVINI,dssvini,INT,ISET)


/* ********************************************************************* */
/* *                                                                   * */
/* *                    TRANSVERSITY DENSITIES                         * */
/* *                                                                   * */
/* *   INPUT:   ISET = number of the parton set :                      * */
/* *              ISET = 1  NEXT-TO-LEADING ORDER  (MS-bar)            *  */
/* *                        (DATA FILE 'transmaxnlo.grid' UNIT=11, TO  * */
/* *                         BE DEFINED BY THE USER )                  *  */
/* *              ISET = 2  LEADING ORDER                              *  */
/* *                        (DATA FILE 'transmaxlo_new.grid' UNIT=22   * */
/* *                                                                   * */
/* *            X  = Bjorken-x       (between  1.E-4  and  1)          * */
/* *            Q2 = scale in GeV**2 (between  0.8  and   1.E6)        * */
/* *             (for values outside the allowed range the program     * */
/* *              writes a warning and extrapolates to the x and       * */
/* *              Q2 values requested)                                 * */
/* *                                                                   * */
/* *   OUTPUT:  U = x * DELTA u                                        *  */
/* *            D = x * DELTA d                                        * */
/* *            UB = x * UBAR                                          * */
/* *            DB = x * DBAR                                          *    */
/* *            ST = x * DELTA STRANGE = x * DELTA STRANGE(BAR)        *      */
/* *                                                                   * */
/* *          (  For the parton distributions always x times           * */
/* *                   the distribution is returned   )                * */
/* *                                                                   * */
/* *                                                                   * */
/* *   COMMON:  The main program or the calling routine has to have    * */
/* *            a common block  COMMON / INTINI / IINI , and  IINI     * */
/* *            has always to be zero when PARPOL is called for the    * */
/* *            first time or when 'ISET' has been changed.            * */
/* *                                                                   * */
/* ********************************************************************* */
PROTOCCALLSFSUB8(PARPOLT,parpolt,INT,DOUBLE,DOUBLE,PDOUBLE,PDOUBLE,PDOUBLE,PDOUBLE,PDOUBLE)
#define  PARPOLT(ISET, X, Q2, U, D, UB, DB, ST) \
CCALLSFSUB8(PARPOLT,parpolt,INT,DOUBLE,DOUBLE,PDOUBLE,PDOUBLE,PDOUBLE,PDOUBLE,PDOUBLE,ISET, X, Q2, U, D, UB, DB, ST)


/*ADAPTIVE GAUSSIAN INTEGRATION IN ../CERNLIB/
C     RES = Estimated Integral of F from A to B,
C     ERR = Estimated absolute error on RES.
C     NSEG  specifies how the adaptation is to be done:
C        =0   means use previous binning,
C        =1   means fully automatic, adapt until tolerance attained.
C        =n>1 means first split interval into n equal segments,
C             then adapt as necessary to attain tolerance.
C     The specified tolerances are:
C            relative: RELTOL ;  absolute: ABSTOL.
C        It stops when one OR the other is satisfied, or number of
C        segments exceeds NDIM.  Either TOLA or TOLR (but not both!)
C        can be set to zero, in which case only the other is used.
C     RES = Estimated Integral of F from A to B,
C     ERR = Estimated absolute error on RES.
C     NSEG  specifies how the adaptation is to be done:
C        =0   means use previous binning,
C        =1   means fully automatic, adapt until tolerance attained.
C        =n>1 means first split interval into n equal segments,
C             then adapt as necessary to attain tolerance.
C     The specified tolerances are:
C            relative: RELTOL ;  absolute: ABSTOL.
C        It stops when one OR the other is satisfied, or number of
C        segments exceeds NDIM.  Either TOLA or TOLR (but not both!)
C        can be set to zero, in which case only the other is used.
      SUBROUTINE DADAPT(F,A,B,NSEG,RELTOL,ABSTOL,RES,ERR)
*/
PROTOCCALLSFSUB8(DADAPT,dadapt,DOUBLEV,DOUBLE,DOUBLE,INT,DOUBLE,DOUBLE,PDOUBLE,PDOUBLE)
#define  DADAPT(F,A,B,NSEG,RELTOL,ABSTOL,RES,ERR)  \
CCALLSFSUB8(DADAPT,dadapt,DOUBLEV,DOUBLE,DOUBLE,INT,DOUBLE,DOUBLE,PDOUBLE,PDOUBLE,F,A,B,NSEG,RELTOL,ABSTOL,RES,ERR)



const int LIN = 0;
const int LOG = 1;


const int standart_nlo = 1; // for polarised PDFs
const int  valence_nlo = 2;
const int standart_lo  = 3;
const int  valence_lo  = 4;


const int transversity_nlo = 1; // for transversity
const int transversity_lo  = 2;



// Class  description
class sidis {
  static sidis* instance;
  double energy;
  double cross_section;

  string experiment_name;


  private:
  PARTONCONTENT unpolarised(double x, double Q2);
  //  PARTONCONTENT unpolarised_pion(double x, double Q2);


  double gauss( double kt, double width );

  void fragmentation(HADRON& h, int charge, double z, double Q2); // unp fragmentation
  void fragmentationCollins(HADRON& h, int charge, double z, double Q2); // collins frag ala transversity

  sidis();

  public:
  static sidis* Instance();

  // Beam + Target --> Produced + X
  HADRON Beam;
  HADRON Target;
  HADRON Produced;

  PARTON Parton_a;
  PARTON Parton_c;


  PARAMETERS Params;

  enum EFFECT effect;
  enum EFFECT calculation;
  enum DEPENDENCE dependence;

  struct LIMITS Limits;
  AVERAGE Average;

  void SetName(string Name) { experiment_name = Name; };
  string GetName(void) {return experiment_name;};



  // Experimental cuts
  bool xF_cut( double x){
    Limits.xF_cut_used = true;
    return x >= Limits.xF_min && x <= Limits.xF_max;
  };

  bool x_cut( double x){
    Limits.x_cut_used = true;
    return x >= Limits.x_min && x <= Limits.x_max;
  };

  bool z_cut( double x){
    Limits.z_cut_used = true;
    return x >= Limits.z_min && x <= Limits.z_max;
  };

  bool y_cut( double x){
    Limits.y_cut_used = true;
    return x >= Limits.y_min && x <= Limits.y_max;
  };

  bool pt_cut( double x){
    Limits.pt_cut_used = true;
    return x >= Limits.pt_min && x <= Limits.pt_max;
  };

  bool kt_cut( double x){
    Limits.kt_cut_used = true;
    return x >= Limits.kt_min && x <= Limits.kt_max;
  };

  bool Q2_cut( double x){
    Limits.Q2_cut_used = true;
    return x >= Limits.Q2_min && x <= Limits.Q2_max;
  };

  bool W2_cut( double x){
    Limits.W2_cut_used = true;
    return x >= Limits.W2_min && x <= Limits.W2_max;
  };

  bool Eh_cut( double x){
    Limits.Eh_cut_used = true;
    return x >= Limits.Eh_min && x <= Limits.Eh_max;
  };

  //description
  string process(void);
  string print_energy(void);
  string limit_Y(void);
  string limit_x(void);
  string limit_z(void);
  string limit_xF(void);
  string limit_PT(void);
  string limit_kt(void);
  string limit_Q2(void);
  string limit_W2(void);
  string limit_Eh(void);
  TPaveText* description( void );

  string HadronDescription( HADRON h );

  // Settings
  double SetEnergy(double Energy) { return energy = Energy; }
  double GetEnergy() const { return energy; }
  double GetCrossSection() const {return cross_section;};

  // Additional
  void Deutron( struct PARTONCONTENT& x );
  void Neutron( struct PARTONCONTENT& x );
  void Antiproton( struct PARTONCONTENT& x );
  void PionNegative( struct PARTONCONTENT& partcontent );




  //Helicity distribution
  PARTONCONTENT helicity(double x, double Q2); // helicity distribution...
  void HelicityDistribution(HADRON& h,double x, double Q2);


  //Unpolarised distributions
  void Unpolarised(HADRON& h, double x, double Q2);

  // Pion distribution
  PARTONCONTENT unpolarised_pion(double x, double Q2);

  // Unpolarised Fragmentation
   void Fragmentation(HADRON& h, double z, double Q2);



  //Collins
  void FragmentationCollins(HADRON& h, double z, double Q2); // ala transversity
  double collins_z_dependence(double z, double a, double b, double n);
  PARTONCONTENT collins( double z ); // BAD WORK ON REALIZATION! should know if pion or kaon
  PARTONCONTENT collinsFirstMoment( double z ); // BAD WORK ON REALIZATION!
  void FragmentationCollinsNLO(HADRON& h, double z, double Q2); // using hoppet u,d--> pi+
  void CollinsDistributionNLO( HADRON& h, double z, double Q2); // at NLO

  // Sivers
  double sivers_x_dependence(double x, double a, double b, double n, double gamma, double eta);
  PARTONCONTENT sivers( double x );
  void SiversDistributionFirstMoment(HADRON& h, double x, double Q2);
  void F1TDistributionFirstMoment( HADRON& h,  double x, double Q2);



  double GetFavoured();          // first favoured
  double GetAntiFavoured();      // seconf favoured
  double GetUnFavoured();        // first unfavoured
  double GetUnFavouredStrange(); // second unfavoured
  // some auxiliary function
  double sum_collins( PARTONCONTENT az1, PARTONCONTENT bz2, PARTONCONTENT bz1, PARTONCONTENT az2 ); // some auxiliary function


  void collins_distribution( HADRON& h, PARTONCONTENT a); // USEFUL FUNCTION TO RE_SHUFFLE THINGS FOR COLLINS DISTRIBUTION...


  void CollinsDistribution( HADRON& h, double z, double Q2);   // BAD WORK ON REALIZATION!
  void CollinsDistributionFirstMoment( HADRON& h, double z, double Q2);


  // Transversity
  double transversity_x_dependence(double x, double a, double b, double n, double gamma, double eta);
  PARTONCONTENT transversity( double x );
  PARTONCONTENT SofferBound( double x, double Q2);

  void TransversityDistributionNLO( HADRON& h, double x, double Q2); //using hoppet
  void TransversityDistribution(HADRON& h,double x, double Q2);
  void TransversityDistributionFirstMoment(HADRON& h,double x, double Q2);




//========================================================= Approximate sivers&collins in PP...
  int approximate_collins_pp_calculation( double xa, double xb, double z, double y, double pt, double Q2 );
  int approximate_collins_pp_calculation_NLO( double xa, double xb, double zc, double y, double pt, double Q2 ); // NLO




// Everything needed for Puion in Jet calculation
  int approximate_pionjet_pp_calculation_NLO( double xb, double y, double zc, double zhat, double bt, double pt, double jt, double Q2 );
  int approximate_pionjet_pp_calculation( double xb, double y, double zc, double pt, double jt, double Q2 );


//========================================================= Approximate sivers&collins in PP...
  int approximate_sivers_pp_calculation( double xa, double xb, double z, double y, double pt, double Q2 );
  int approximate_sivers_pp_jet_calculation( double xb, double y, double pt, double Q2 );

//========================================================= Approximate sivers in PP for direct photon production...
  int approximate_sivers_pp_photon_calculation( double xb, double y, double pt, double Q2 );


  PARTONCONTENT  collinsff_btpart_nonpert_Feng( HADRON& h, double bt, double x);
  PARTONCONTENT  unpolarisedff_btpart_nonpert_Feng( HADRON& h, double bt, double x);
  PARTONCONTENT  unpolarised_btpart_nonpert_Feng(double bt, double x);

  double unpolarisedwidth_x_dependence_Feng(double x, double a, double b, double n, double gamma, double eta);


};




void hold(bool);

// the initial condition WILL BE USED IN PDF FITTING
void  pdf_init(const double & x,
                   const double & Q,
                   double * pdf);

void set_hoppet( void );

void reset_hoppet( void );

#endif //#ifndef  __SIDIS_H__
