C *********************************************************************
C  aac03.f  Version 2.2                                  Jan. 21, 2005
C
C  [Package for the AAC polarized parton distribution functions]
C *********************************************************************
C  AAC Polarized Parton Distribution Functions.
C      (Asymmetry Analysis Collaboration)
C
C  This package contains AAC03-NLO PDFs and gradient terms of 
C  these distributions. 
C
C  REFERENCE:
C    "Determination of polarized parton distribution functions 
C     and their uncertainties",
C     M. Hirai, S. Kumano, N. Saito, 
C     (Asymmetry Analysis Collaboration),
C     Phys. Rev. D69, 054021 (2004)
C
C  QUESTION OR COMMENT TO:
C    Masanori Hirai     [mhirai@rarfaxp.riken.jp]
C ---------------------------------------------------------------------
C  SUBROUTINE AACPDFE(Q2,X,ISET,POLPDF,GRADPDF):
C
C   Subroutine AACPDF returns the values of AAC polarized parton
C   distribution functions at specified Q^2 and x point 
C   by interpolating the grid data.
C   [ Log(Q^2): LINEAR INTERPOLATION, x: CUBIC SPLINE INTERPOLATION ]
C
C   INPUT:
C
C     Q2, X   ... Q^2 and x values at which the distributions are
C                 calculated. Available range: 10^-9 <= X <= 1.0,
C                 1.0 GeV^2 <= Q2 <= 10^8 GeV^2.
C     ISET=1: AAC03 polarized PDFs 
C          2: Gradient terms of AAC PDFs: d x Delta f(x)/d a_i
C          3: Polarized PDFs and their gradient terms 
C
C   OUTPUT: Arrays POLPDF(-3:3)
C
C     POLPDF(I) --> AAC03 polarized PDFs.
C      I = -3 ... s-bar quark (x Delta s-bar = x Delta s)
C          -2 ... d-bar quark (x Delta d-bar = x Delta d_sea)
C          -1 ... u-bar quark (x Delta u-bar = x Delta u_sea)
C           0 ... gluon (x Delta g)
C           1 ... u quark (x Delta u = x Delta u_v + x Delta u_sea)
C           2 ... d quark (x Delta d = x Delta d_v + x Delta d_sea)
C           3 ... s quark (x Delta s = x Delta s-bar)
C
C     GRADPDF(I) --> Gradient terms of AAC03 polarized PDFs
C      I = 1..4: gluon (delta, nu, kappa, mu)
C          5, 6: u_v (delta, mu)
C          7, 8: d_v (delta, mu)
C         9..11: antiquark (delta, nu, kappa)
C
C   NOTE: The returned values are the distributions multiplied by x.
C
C      *  Error matrix can be used by declaring a common block:
C         COMMON/ERRM/EM(11,11). This matrix is defined as multiplying
C         Delta chi^2 by inverse matrix of Hessian: 
C         EM(i,j)=Delta chi^2*H_ij^-1 (Delta ch^2=12.647)
C
C *********************************************************************
      SUBROUTINE AACPDFE(Q2,X,ISET,POLPDF,GRADPDF)
C ---------------------------------------------------------------------
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER (NQ=33, NX=117, ND=18, NSET=3, IFILE=10)
      DIMENSION IREAD(2),QG(NQ), XG(NX),PDFJ1(ND), PDFJ2(ND)
     +         ,POLPDF(-3:3), GRADPDF(11)
     +         ,BXG(NX,NQ,ND), CXG(NX,NQ,ND), DXG(NX,NQ,ND)
     +         ,PDFG(NX,NQ,ND)

      COMMON/ERRM/EM(11,11)
      SAVE IREAD, BXG, CXG, DXG
      DATA IREAD / 1, 1 /

C Q2 AND X GRID.
      DATA QG /
     +  1.000000D+00, 1.467799D+00, 2.154435D+00,
     +  3.162278D+00, 4.641589D+00, 6.812921D+00,
     +  1.000000D+01, 1.467799D+01, 2.154435D+01,
     +  3.162278D+01, 4.641589D+01, 6.812921D+01,
     +  1.000000D+02, 1.778279D+02, 3.162278D+02, 5.623413D+02,
     +  1.000000D+03, 1.778279D+03, 3.162278D+03, 5.623413D+03,
     +  1.000000D+04, 1.778279D+04, 3.162278D+04, 5.623413D+04,
     +  1.000000D+05, 1.778279D+05, 3.162278D+05, 5.623413D+05,
     +  1.000000D+06, 4.641589D+06, 
     +  1.000000D+07, 4.641589D+07,  
     +  1.000000D+08  /

      DATA XG / 
     +  1.000000D-09, 1.333521D-09, 1.778279D-09, 2.371374D-09,
     +  3.162278D-09, 4.216965D-09, 5.623413D-09, 7.498942D-09,
     +  1.000000D-08, 1.333521D-08, 1.778279D-08, 2.371374D-08,
     +  3.162278D-08, 4.216965D-08, 5.623413D-08, 7.498942D-08,
     +  1.000000D-07, 1.333521D-07, 1.778279D-07, 2.371374D-07,
     +  3.162278D-07, 4.216965D-07, 5.623413D-07, 7.498942D-07,
     +  1.000000D-06, 1.333521D-06, 1.778279D-06, 2.371374D-06,
     +  3.162278D-06, 4.216965D-06, 5.623413D-06, 7.498942D-06,
     +  1.000000D-05, 1.333521D-05, 1.778279D-05, 2.371374D-05,
     +  3.162278D-05, 4.216965D-05, 5.623413D-05, 7.498942D-05,
     +  1.000000D-04, 1.333521D-04, 1.778279D-04, 2.371374D-04,
     +  3.162278D-04, 4.216965D-04, 5.623413D-04, 7.498942D-04,
     +  1.000000D-03, 1.154782D-03, 1.333521D-03, 1.539927D-03,
     +  1.778279D-03, 2.053525D-03, 2.371374D-03, 2.738420D-03,
     +  3.162278D-03, 3.651741D-03, 4.216965D-03, 4.869675D-03,
     +  5.623413D-03, 6.493816D-03, 7.498942D-03, 8.659643D-03,
     +  1.000000D-02, 1.154782D-02, 1.333521D-02, 1.539927D-02,
     +  1.778279D-02, 2.053525D-02, 2.371374D-02, 2.738420D-02,
     +  3.162278D-02, 3.651741D-02, 4.216965D-02, 4.869675D-02,
     +  5.623413D-02, 6.493816D-02, 7.498942D-02, 8.659643D-02,
     +  1.000000D-1, 1.250000D-1, 1.500000D-1, 1.750000D-1,
     +  2.000000D-1, 2.250000D-1, 2.500000D-1, 2.750000D-1,
     +  3.000000D-1, 3.250000D-1, 3.500000D-1, 3.750000D-1,
     +  4.000000D-1, 4.250000D-1, 4.500000D-1, 4.750000D-1, 
     +  5.000000D-1, 5.250000D-1, 5.500000D-1, 5.750000D-1,
     +  6.000000D-1, 6.250000D-1, 6.500000D-1, 6.750000D-1,
     +  7.000000D-1, 7.250000D-1, 7.500000D-1, 7.750000D-1,
     +  8.000000D-1, 8.250000D-1, 8.500000D-1, 8.750000D-1,
     +  9.000000D-1, 9.250000D-1, 9.500000D-1, 9.750000D-1,
     +  1.000000D+0 /

C ERROR MATRIX EM(i,j)=Dchi^2*H_ij^-1 (Dch^2=12.647)
      DATA ((EM(I,J),J=1,11), I=1,11)/
     + 1.985579D+02,-2.061461D+01,-2.200578D+01,-6.057913D+00,
     + 6.930556D-01, 2.428224D+00, 2.289107D+00,-8.878194D-01,
     +-5.223211D+00, 2.592635D+00,-1.277347D+00,
     +-2.061461D+01, 1.504993D+01, 2.529400D+00,-8.979370D+00,
     +-8.081433D-02,-3.161750D-01,-2.567341D-01, 3.073221D-01,
     + 1.644110D-01,-4.692037D-01, 1.024407D-01,
     +-2.200578D+01, 2.529400D+00, 6.399382D+00, 1.745286D+00,
     +-1.454405D-01,-4.944977D-01,-5.463504D-01, 1.656757D-01,
     + 1.033260D+00,-5.248505D-01, 2.832928D-01,
     +-6.057913D+00,-8.979370D+00, 1.745286D+00, 8.916135D+00,
     +-4.375862D-02,-1.492346D-01,-1.479699D-01, 5.450857D-02,
     + 2.301754D-01,-1.429111D-01, 7.638788D-02,
     + 6.930556D-01,-8.081433D-02,-1.454405D-01,-4.375862D-02,
     + 1.236877D-01, 2.934104D-01,-8.903488D-02,-3.781453D-02,
     + 4.944977D-01, 2.782340D-02, 3.288220D-03,
     + 2.428224D+00,-3.161750D-01,-4.944977D-01,-1.492346D-01,
     + 2.934104D-01, 8.941429D-01, 2.643223D-02,-2.883516D-01,
     + 1.884403D+00, 5.855561D-02, 2.757046D-02,
     + 2.289107D+00,-2.567341D-01,-5.463504D-01,-1.479699D-01,
     +-8.903488D-02, 2.643223D-02, 1.795874D+00,-1.568228D+00,
     +-3.326161D-01, 1.656757D-01,-5.147329D-02,
     +-8.878194D-01, 3.073221D-01, 1.656757D-01, 5.450857D-02,
     +-3.781453D-02,-2.883516D-01,-1.568228D+00, 2.744399D+00,
     +-1.732639D+00,-6.854674D-02,-3.515866D-02,
     +-5.223211D+00, 1.644110D-01, 1.033260D+00, 2.301754D-01,
     + 4.944977D-01, 1.884403D+00,-3.326161D-01,-1.732639D+00,
     + 1.248259D+01,-1.378523D+00, 5.400269D-01,
     + 2.592635D+00,-4.692037D-01,-5.248505D-01,-1.429111D-01,
     + 2.782340D-02, 5.855561D-02, 1.656757D-01,-6.854674D-02,
     +-1.378523D+00, 3.958511D-01,-9.207016D-02,
     +-1.277347D+00, 1.024407D-01, 2.832928D-01, 7.638788D-02,
     + 3.288220D-03, 2.757046D-02,-5.147329D-02,-3.515866D-02,
     + 5.400269D-01,-9.207016D-02, 3.604395D-02/

C CALCULATE SPLINE COEFFICIENTS.
      IF(IREAD(1).NE.1) GO TO 20

C READ GRID DATA AND CALCULATE SPLINE COEFFICIENTS. 
      IF((ISET.GT.0).AND.(ISET.LT.4))THEN
        OPEN(UNIT=IFILE,FILE='aacnlo.grd',STATUS='OLD')
        OPEN(UNIT=IFILE+1,FILE='gradn.grd',STATUS='OLD')
      ELSE
        WRITE(*,1010) ISET
 1010   FORMAT(' ','AACPDF ERROR: ISET =', I3)
        STOP
      END IF

      DO J=1,NQ
        DO K=1,NX-1
          READ(IFILE,1020) (PDFG(K,J,I), I=1,6)
          READ(IFILE+1,1020) (PDFG(K,J,I), I=7,17)
 1020     FORMAT(1X,11(1PE14.6))
        ENDDO
      ENDDO
      CLOSE(IFILE)
      CLOSE(IFILE+1)

      DO I=1,17
        DO J=1,NQ
          PDFG(NX,J,I)=0.D0 ! x=1 PolPDF=0.D0
          CALL SPLINE(NX,XG,PDFG,BXG,CXG,DXG,ISET,I,J)
        ENDDO
      ENDDO

      IREAD(1)=2
   20 CONTINUE

      DO I=1,7
        POLPDF(I-4)=0.D0
      END DO
      DO I=1,11
        GRADPDF(I)=0.D0
      ENDDO

      IF(ISET.EQ.1) THEN
        NDS=1
        NDE=6
      ELSE IF (ISET.EQ.2) THEN
        NDS=7
        NDE=17
      ELSE 
        NDS=1
        NDE=17
      ENDIF

C CHECK X AND Q2 VALUES.
      IF((X.LT.1.D-9).OR.(X.GT.1.D0)) THEN
        WRITE(*,1030) X
 1030   FORMAT (' ','AACPDF WARNING: OUT OF RANGE --> X =', 1PE12.3)
        STOP
      ENDIF
      IF((Q2.LT.1.D0).OR.(Q2.GT.1.D8)) THEN
        WRITE(*,1040) Q2
 1040   FORMAT (' ','AACPDF WARNING: OUT OF RANGE --> Q2 =', 1PE12.3)
        STOP
      ENDIF

C INTERPOLATION.
C X: CUBIC SPLINE INTERPOLATION, LOG(Q2): LINEAR INTERPOLATION.
      J=ISERCH(NQ,QG,Q2)
      IF(J.EQ.NQ) J=NQ-1
      K=ISERCH(NX,XG,X)
      DO I=NDS,NDE
        DX=X-XG(K)
        PDFJ1(I)=PDFG(K,J,I)
     1       +DX*(BXG(K,J,I)+DX*(CXG(K,J,I)+DX*DXG(K,J,I)))
        PDFJ2(I)=PDFG(K,J+1,I)
     1       +DX*(BXG(K,J+1,I)+DX*(CXG(K,J+1,I)+DX*DXG(K,J+1,I)))
      ENDDO

      T=(DLOG(Q2)-DLOG(QG(J)))/(DLOG(QG(J+1))-DLOG(QG(J)))
      IF(ISET.NE.1) THEN
        DO I=1,11
          GRADPDF(I)=(1.D0-T)*PDFJ1(I+6)+T*PDFJ2(I+6)
        ENDDO

        IF(ISET.EQ.3) THEN
          DO I=1,3
            POLPDF(I-1)=(1.D0-T)*PDFJ1(I)+T*PDFJ2(I)
            POLPDF(-I)=(1.D0-T)*PDFJ1(I+3)+T*PDFJ2(I+3)
          ENDDO
          POLPDF(3)=POLPDF(-3)
        END IF

      ELSE
        DO I=1,3
          POLPDF(I-1)=(1.D0-T)*PDFJ1(I)+T*PDFJ2(I)
          POLPDF(-I)=(1.D0-T)*PDFJ1(I+3)+T*PDFJ2(I+3)
        ENDDO
        POLPDF(3)=POLPDF(-3)
      ENDIF


      RETURN
      END
C ---------------------------------------------------------------------
      SUBROUTINE SPLINE(N,X,Y,B,C,D,ISET,I,J)
C ---------------------------------------------------------------------
C CALCULATE THE COEFFICIENTS B,C,D IN A CUBIC SPLINE INTERPOLATION.
C INTERPOLATION SUBROUTINES ARE TAKEN FROM
C G.E. FORSYTHE, M.A. MALCOLM AND C.B. MOLER,
C COMPUTER METHODS FOR MATHEMATICAL COMPUTATIONS (PRENTICE-HALL, 1977).
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER (NQ=33, NX=117, ND=18)
      DIMENSION Y(NX,NQ,ND),B(NX,NQ,ND),C(NX,NQ,ND),D(NX,NQ,ND)
     1         ,X(NX) 
      NM1=N-1
      IF(N.LT.2) RETURN
      IF(N.LT.3) GO TO 250
      D(1,J,I)=X(2)-X(1)
      C(2,J,I)=(Y(2,J,I)-Y(1,J,I))/D(1,J,I)
      DO 210 K=2,NM1
        D(K,J,I)=X(K+1)-X(K)
        B(K,J,I)=2.0D0*(D(K-1,J,I)+D(K,J,I))
        C(K+1,J,I)=(Y(K+1,J,I)-Y(K,J,I))/D(K,J,I)
        C(K,J,I)=C(K+1,J,I)-C(K,J,I)
  210 CONTINUE
      B(1,J,I)=-D(1,J,I)
      B(N,J,I)=-D(N-1,J,I)
      C(1,J,I)=0.0D0
      C(N,J,I)=0.0D0
      IF(N.EQ.3) GO TO 215
      C(1,J,I)=C(3,J,I)/(X(4)-X(2))-C(2,J,I)/(X(3)-X(1))
      C(N,J,I)=C(N-1,J,I)/(X(N)-X(N-2))-C(N-2,J,I)/(X(N-1)-X(N-3))
      C(1,J,I)=C(1,J,I)*D(1,J,I)**2.0D0/(X(4)-X(1))
      C(N,J,I)=-C(N,J,I)*D(N-1,J,I)**2.0D0/(X(N)-X(N-3))
  215 CONTINUE
      DO 220 K=2,N
        T=D(K-1,J,I)/B(K-1,J,I)
        B(K,J,I)=B(K,J,I)-T*D(K-1,J,I)
        C(K,J,I)=C(K,J,I)-T*C(K-1,J,I)
  220 CONTINUE
      C(N,J,I)=C(N,J,I)/B(N,J,I)
      DO 230 IB=1,NM1
        K=N-IB
        C(K,J,I)=(C(K,J,I)-D(K,J,I)*C(K+1,J,I))/B(K,J,I)
  230 CONTINUE
      B(N,J,I)=(Y(N,J,I)-Y(NM1,J,I))/D(NM1,J,I)
     1        +D(NM1,J,I)*(C(NM1,J,I)+2.0D0*C(N,J,I))
      DO 240 K=1,NM1
        B(K,J,I)=(Y(K+1,J,I)-Y(K,J,I))/D(K,J,I)
     1          -D(K,J,I)*(C(K+1,J,I)+2.0D0*C(K,J,I))
        D(K,J,I)=(C(K+1,J,I)-C(K,J,I))/D(K,J,I)
        C(K,J,I)=3.0D0*C(K,J,I)
  240 CONTINUE
      C(N,J,I)=3.0D0*C(N,J,I)
      D(N,J,I)=D(N-1,J,I)
      RETURN
  250 CONTINUE
      B(1,J,I)=(Y(2,J,I)-Y(1,J,I))/(X(2)-X(1))
      C(1,J,I)=0.0D0
      D(1,J,I)=0.0D0
      B(2,J,I)=B(1,J,I)
      C(2,J,I)=0.0D0
      D(2,J,I)=0.0D0
      RETURN
      END
C ---------------------------------------------------------------------
      INTEGER FUNCTION ISERCH(N,X,Y)
C ---------------------------------------------------------------------
C THIS FUNCTION SEARCHES "I" WHICH SATISFIES THE RELATION
C X(I) <= Y < X(I+1) BY USING A BINARY SEARCH.
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION X(118)

      MIN=1
      MAX=N+1

   10 CONTINUE
      MID=(MIN+MAX)/2
      IF(Y.LT.X(MID)) THEN
        MAX=MID
      ELSE
        MIN=MID
      END IF
      IF((MAX-MIN).GT.1) GO TO 10

      ISERCH=MIN

      RETURN
      END
C *********************************************************************
C THE END OF THE PROGRAM.
C *********************************************************************
