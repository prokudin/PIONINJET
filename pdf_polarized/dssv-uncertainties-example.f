C --- EXAMPLE HOW TO USE HESSIAN EIGENVECTOR SETS
C --- COMPUTES (TRUNCATED) 1ST MOMENT OF (DELTA U + DELTA UBAR)
C --- AND ITS HESSIAN UNCERTAINTY
C --- 
      PROGRAM EXAMPLE
      IMPLICIT DOUBLE PRECISION (A-Z)
      integer iset
      external func

      sum = 0.D0

      call dssvini(0)           !COMPUTE BEST FIT VALUE
      y0 = dinteg(func,0.001d0,1.d0,1.d-3)

      do 1 iset=1, 19           !LOOP THROUGH EIGENVECTOR SETS
      call dssvini(iset)
      y1 = dinteg(func,0.001d0,1.d0,1.d-3) !COMPUTE TRUNCATED 1ST MOMENT
      call dssvini(-iset)
      y2 = dinteg(func,0.001d0,1.d0,1.d-3)

      sum= sum + (y1-y2)**2     !COMPUTE UNCERTAINTY

      write(6,3) iset,y0,y1,y2
 3    format(I3, 3(1pe12.4,1x))
 1    continue 

      write(6,4) y0,0.5*sqrt(sum)
 4    format(3(1pe12.4,1x))

      end

      double precision function func(x)
      implicit double precision (a-z)

      q2=10.d0                  !CHOOSE Q2 VALUE

      CALL DSSVFIT(X,Q2,DUV,DDV,DUBAR,DDBAR,DSTR,DGLU)

      func=(duv+2.d0*dubar)/x   !CHOOSE INTEGRAND, HERE U+UBAR
      return
      end

C
C...DOUBLE PRECISE GAUSS INTEGRATION :
       FUNCTION DINTEG (F, ALFA, BETA, EPS)
       IMPLICIT DOUBLE PRECISION (A-H, O-Z)
       DIMENSION W(12), X(12)
       DATA CONST / 1.0 D-12 /
       DATA W
     1  /0.10122 85362 90376, 0.22238 10344 53374, 0.31370 66458 77887,
     2   0.36268 37833 78362, 0.02715 24594 11754, 0.06225 35239 38647,
     3   0.09515 85116 82492, 0.12462 89712 55533, 0.14959 59888 16576,
     4   0.16915 65193 95002, 0.18260 34150 44923, 0.18945 06104 55068/
       DATA X
     1  /0.96028 98564 97536, 0.79666 64774 13627, 0.52553 24099 16329,
     2   0.18343 46424 95650, 0.98940 09349 91649, 0.94457 50230 73232,
     3   0.86563 12023 87831, 0.75540 44083 55003, 0.61787 62444 02643,
     4   0.45801 67776 57227, 0.28160 35507 79258, 0.09501 25098 37637/
       DINTEG = 0.0 D0
       IF ( ALFA . EQ. BETA ) RETURN
       A = ALFA
       B = BETA
       DELTA = CONST * (DABS(A-B))
       AA = A
    1  Y = B - AA
       IF( DABS(Y) .LE. DELTA ) RETURN
    2  BB = AA + Y
       C1 = 0.5 D0 * (AA + BB)
       C2 = C1 - AA
       S8 = 0.0 D0
       S16 = 0.0 D0
       DO 15 I = 1, 4
          C3 = X(I) * C2
          S8 = S8 + W(I) * (F(C1+C3) + F(C1-C3))
   15  CONTINUE
       DO 16 I = 5, 12
          C3 = X(I) * C2
          S16 = S16 + W(I) * (F(C1+C3) + F(C1-C3))
  16   CONTINUE
       S8 = S8 * C2
       S16= S16 * C2
       IF( DABS(S16-S8) .GT. EPS * DABS(S8)) THEN
          Y = 0.5 * Y
          IF ( DABS(Y) .LE. DELTA ) THEN
             DINTEG = 0.0
             WRITE (*,10)
  10         FORMAT (1X,' DINTEG : TOO HIGH ACCURACY ')
          ELSE
             GOTO 2
          END IF
       ELSE
          DINTEG = DINTEG + S16
          AA = BB
          GOTO 1
       END IF
       RETURN
       END
