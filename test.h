#ifndef __TEST_H__
#define __TEST_H__

#include <cfortran.h>


// FRAGMENTATION LIB
PROTOCCALLSFSUB9(DLIB,dlib,DOUBLE,DOUBLE,PVOID,INT,INT,INT,INT,INT,INT)
#define DLIB(Z,Q2,DFF,FFSET,FFORDER,IHADRON,ICHARGE,ICP,IPI) \
  CCALLSFSUB9(MYDLIB,dlib,DOUBLE,DOUBLE,PVOID,INT,INT,INT,INT,INT,INT,\
	      Z,Q2,DFF,FFSET,FFORDER,IHADRON,ICHARGE,ICP,IPI)


  
  
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



#endif // #ifndef __TEST_H__