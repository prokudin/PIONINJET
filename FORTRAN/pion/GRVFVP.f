*
* $Id: grvfvp.F,v 1.1 2010/02/19 16:42:50 decowski Exp $
*
* $Log: grvfvp.F,v $
* Revision 1.1  2010/02/19 16:42:50  decowski
*
* [decowski]
*
* CERN lib libpdf804 with CTEQ6 structure functions instead the default
* CTEQ5 (necessary for simulating high E interactions).
*
* Modified by Aart Heijboer and Teresa Montaruli
*
* Revision 1.1.1.2  1996/10/30 08:28:38  cernlib
* Version 7.04
*
* Revision 1.1.1.1  1996/04/12 15:29:23  plothow
* Version 7.01
*
* ALEXEI PROKUDIN 2/9/2017 IMPLICIT REAL*8 (A - Z)
       FUNCTION GRVFVP (X, N, AK, AG, D)
c       IMPLICIT REAL (A - Z)
       IMPLICIT REAL*8 (A - Z)
       DX = SQRT (X)
       GRVFVP = N * X**AK * (1.+ AG*DX) * (1.- X)**D
       RETURN
       END
