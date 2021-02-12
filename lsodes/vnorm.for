
! this function routine computes the weighted root-mean-square norm of the vector of length n contained in the array v, 
! with weights contained in the array w of length n..
! vnorm = sqrt( (1/n) * sum( v(i)*w(i) )**2 )

      double precision function vnorm(n, v, w)
      integer n, i
      double precision v, w, sum
      dimension v(n), w(n)
      sum = 0.0d0
      do i = 1,n
         sum = sum + (v(i)*w(i))**2
      end do
      vnorm = dsqrt(sum/dfloat(n))
      return
      end
