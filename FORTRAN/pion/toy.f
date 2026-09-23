      program toy
      IMPLICIT REAL*8 (A - Y)
      Double Precision
     +     ZX,ZQ,ZUV,ZDV,ZUDB,ZSB,ZCB,ZBB,ZGL
      real*8  CHM,BOT,GLU
      integer n1
      real*8 x, Q2
      integer nd, ni, nj, MODE
      MODE = 2 ! best fit

 
      x  = 0.1
      Q = 9.

      xmin = 0.0001
      xmax = 0.9999
      Q2min = 1.0000
      Q2max = 90.

      nd = 100
      stepx = (xmax - xmin)/nd
      stepQ2 = (Q2max - Q2min)/nd
      x = xmin
      Q2 = Q2min
      

      open (unit = 6, file = "pion_mrss.dat")
      do ni = 0, nd, 1
      x = xmin + stepx*ni
      do nj = 0, nd, 1
      Q2 = Q2min + stepQ2*nj
c      call GRVPILO (x, dsqrt(Q),ZUV,ZDV,ZUDB,ZSB,ZCB,ZBB,ZGL) ! GRV
c     write(6,*) x, Q, ZUV,ZDV,ZSB,ZUDB,ZUDB,ZSB
      call PION(X,dsqrt(Q),MODE,ZUV,ZDV,ZUDB,ZSB,CHM,BOT,GLU) ! MRSS
      write(6,*) x, Q, ZUV/x,ZDV/x,ZSB/x,ZUDB/x,ZUDB/x,ZSB/x
      enddo
      enddo

      end
