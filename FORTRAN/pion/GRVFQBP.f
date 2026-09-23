*
* $Id: grvfqbp.F,v 1.1 2010/02/19 16:42:50 decowski Exp $
*
* $Log: grvfqbp.F,v $
* Revision 1.1  2010/02/19 16:42:50  decowski
*
* [decowski]
*
* CERN lib libpdf804 with CTEQ6 structure functions instead the default
* CTEQ5 (necessary for simulating high E interactions).
*
* Modified by Aart Heijboer and Teresa Montaruli
*
* Revision 1.1.1.2  1996/10/30 08:28:37  cernlib
* Version 7.04
*
* Revision 1.1.1.1  1996/04/12 15:29:23  plothow
* Version 7.01
*
*
* ALEXEI PROKUDIN 2/9/2017 IMPLICIT REAL*8 (A - Z)

       FUNCTION GRVFQBP (X, S, ST, AL, BE, AK, AG, B, D, E, ES)
c     IMPLICIT REAL (A - Z)
       IMPLICIT REAL*8 (A - Z)
       DX = DSQRT (X)
       LX = DLOG (1./X)
       IF (S .LE. ST) THEN
          GRVFQBP = 0.0
       ELSE
          GRVFQBP = (S-ST)**AL / LX**AK * (1.+ AG*DX + B*X) * (1.- X)**D
     1           * EXP (-E + DSQRT (ES * S**BE * LX))
       END IF
       RETURN
       END
