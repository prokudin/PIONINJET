	DOUBLE PRECISION FUNCTION F(X)
	DOUBLE PRECISION X
           F =  DEXP(-X*X/0.25D0)
        RETURN
        END


	PROGRAM TOY

        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
 

        EXTERNAL F
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
 

  	A = 0.
	B = 4.
	NSEG = 1
	RELTOL = 1.E-4
	ABSTOL = 0
        CALL DADAPT(F,A,B,NSEG,RELTOL,ABSTOL,RES,ERR)

        PRINT*,RES,ERR
        END
