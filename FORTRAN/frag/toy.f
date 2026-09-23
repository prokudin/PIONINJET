      program toy
      implicit double precision ( a-z )
      dimension dff(-5:5)
      integer ffset, fforder, iparton, ipi, ifini, old1, 
     >        iset, ihadron, icharge, old2, old3, ih, icp
      integer n, n1
      data ipi / 1 /  


      ffset    = 6 ! HKNS
      fforder  = 1
      ihadron  = 1
      icharge  = 3
      icp      = 1

      z  = 0.2d0
      Q2 = 1.d4

      do n1 = 0, 10
      z = z + 0.01 * n1
      call dlib(z,Q2,dff,ffset,fforder,ihadron,icharge,icp,ipi) 
      write(6,*) dff

      su3sing = 0.d0
      do 10 n = -3, 3
      if (n .eq. 0) goto 10
      su3sing = su3sing + dff(n)
 10   continue
      write(6,*) su3sing
      enddo

      end
