*
* $Id: grvfgp.F,v 1.1 2010/02/19 16:42:50 decowski Exp $
*
* $Log: grvfgp.F,v $
* Revision 1.1  2010/02/19 16:42:50  decowski
*
* [decowski]
*
* CERN lib libpdf804 with CTEQ6 structure functions instead the default
* CTEQ5 (necessary for simulating high E interactions).
*
* Modified by Aart Heijboer and Teresa Montaruli
*
* Revision 1.1.1.2  1996/10/30 08:28:36  cernlib
* Version 7.04
*
* Revision 1.1.1.1  1996/04/12 15:29:23  plothow
* Version 7.01
*
*
* ALEXEI PROKUDIN 2/9/2017 IMPLICIT REAL*8 (A - Z)
      
       FUNCTION GRVFGP (X, S, AL, BE, AK, BK, AG, BG, C, D, E, ES)
c       IMPLICIT REAL (A - Z)
       IMPLICIT REAL*8 (A - Z)
       DX = DSQRT (X)
       LX = DLOG (1./X)
       GRVFGP = (X**AK * (AG + BG*DX + C*X) * LX**BK + S**AL
     1       * EXP (-E + SQRT (ES * S**BE * LX))) * (1.- X)**D
       RETURN
       END
