C *********************************************************************
      PROGRAM SAMPLE
C ---------------------------------------------------------------------
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER (NX=200)
      DIMENSION XPDF(-3:3),E(5),AA(7),GRADPDF(11)

      double precision U, D, UB, DB, ST, GL, G1P, G1N
      integer iset

      DATA XMIN,XMAX/1.D-4,1.0D0/
      DATA Q2/1.D0/           ! Q^2[GeV^2]

      COMMON/ERRM/EM(11,11)   ! Error matrix: Dchi^2*H_ij^-1
                              ! EM(11,11) is defined in aac03e.f

      XLSTEP=(DLOG10(XMAX)-DLOG10(XMIN))/DFLOAT(NX)
      DO L=1,NX+1
        DLOGX=DFLOAT(L-1)*XLSTEP+DLOG10(XMIN)
        X=10.D0**(DLOGX) 

                   ! ISET=3: AAC Pol-PDFs & Gradient terms
        Call AACPDFE(Q2,X,3,XPDF,GRADPDF)  

        iset = 3  ! 'STANDARD' SCENARIO LO


        call PARPOL(ISET, X, Q2, U, D, UB, DB, ST, GL, G1P, G1N)
c        write(14,1020) X, Q2, U, D, UB, DB, ST, GL, G1P, G1N
        write(14,1020) X, GL, U-UB, D-DB, UB, DB, ST, GL, G1P, G1N
  
  
C Plarized PDFs by AAC03 analysis
C XPDF(I) (0: gluon, 1:u=uv+usea, 2:d=dv+dsea, 3:strange, -I:antiquarks)
        write(12,1010) X,XPDF(0)                           ! xDg
     +                  ,XPDF(1)-XPDF(-1)                  ! xDu_v
     +                  ,XPDF(2)-XPDF(-2)                  ! xDd_v
     +                  ,(XPDF(-1)+XPDF(-2)+XPDF(-3))/3.D0 ! xDqb
     +     ,XPDF(1)+XPDF(-1)+XPDF(2)+XPDF(-2)+2.D0*XPDF(3) ! xDsigma

C Uncertainties of AAC polarized PDFs
C GRADPDF(I) (1..4:gluon, 5&6:uv, 7&8:dv, 9..11:antiquark)
        DO N=1,5
          IF(N.EQ.1)THEN      ! xDg
            NP=4;NS=0
          ELSE IF(N.EQ.2)THEN ! xDuv
            NP=2;NS=4
          ELSE IF(N.EQ.3)THEN ! xDdv
            NP=2;NS=6
          ELSE IF(N.EQ.4)THEN ! xDqb
            NP=3;NS=8
          ELSE                ! singlet(= xDuv + xDdv + 6*xDqb)
            NP=7;NS=4
            DO M=9,11
              GRADPDF(M)=GRADPDF(M)*6.D0
            END DO
          END IF

          DO I=1,7
            AA(I)=0.D0
          ENDDO
          EE=0.D0

          DO I=1,NP     ! PDF uncertainty by the Hessian method 
            DO J=1,NP
              AA(I)=AA(I)+EM(NS+I,NS+J)*GRADPDF(NS+J)
            ENDDO
            EE=EE+GRADPDF(NS+I)*AA(I)
          ENDDO
          E(N)=DSQRT(EE)
        ENDDO
                       !  eDg,  eDuv, eDdv, eDb,  eDsigma
        WRITE(13,1010) X, E(1), E(2), E(3), E(4), E(5)
      ENDDO

 1010 FORMAT(' ', 6(1PE16.7))
 1020 FORMAT(' ', 10(1PE16.7))
      STOP
      END
C *********************************************************************
C THE END OF THE PROGRAM.
C *********************************************************************
