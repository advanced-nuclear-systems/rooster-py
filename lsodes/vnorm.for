      double precision function vnorm (n, v, w)
      integer n, i
      double precision v, w,   sum
      dimension v(n), w(n)
      sum = 0.0d0
      do i = 1,n
         sum = sum + (v(i)*w(i))**2
      end do
      vnorm = dsqrt(sum/dfloat(n))
      return
      end
