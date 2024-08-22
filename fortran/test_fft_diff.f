      program main
      implicit real*8 (a-h,o-z)
      parameter(n=32)
      dimension x(n), a(n), b(n), ax(n)

      pi = 4.0d0*atan(1.0d0)

      xmin = 0.0d0
      xmax = 2.0d0*pi
      dx = (xmax - xmin)/n

      do i=1,n
         x(i) = xmin + (i-1)*dx
         a(i) = exp(sin(x(i)))    ! function
         ax(i) = a(i) * cos(x(i)) ! derivative
         b(i) = 0.0d0
      enddo

      call tofour(a,b,n)
      call difffour(a,b,n)
      call fromfour(a,b,n)
      scale = 2.0d0/(xmax-xmin)
      do i=1,n
         a(i) = scale*a(i)
      enddo

      open(10,file='diff.txt')
      do i=1,n
         write(10,*) x(i),a(i),a(i)-ax(i)
      enddo
      close(10)

      stop
      end program main
