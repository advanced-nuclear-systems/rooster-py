
!     numerical solution of sparse nonsymmetric system of linear equations given ldu-factorization (compressed pointer storage)

      subroutine nnsc(n, r, c, il, jl, ijl, l, d, iu, ju, iju, u, z, b, tmp)
      integer r(1), c(1), il(n), jl(1), ijl(1), iu(1), ju(1), iju(1)
      double precision  l(1), d(1), u(1), b(1), z(1), tmp(1), tmpk, sum

!     set tmp to reordered b
      do k=1,n
         tmp(k) = b(r(k))
      end do
!     solve  ly = b  by forward substitution
      do k=1,n
         jmin = il(k)
         jmax = il(k+1) - 1
         tmpk = -d(k) * tmp(k)
         tmp(k) = -tmpk
         if(jmin .le. jmax)then
            ml = ijl(k) - jmin
            do j=jmin,jmax
               tmp(jl(ml+j)) = tmp(jl(ml+j)) + tmpk * l(j)
            end do
         end if
      end do
!     solve  ux = y  by back substitution
      k = n
      do i=1,n
         sum = -tmp(k)
         jmin = iu(k)
         jmax = iu(k+1) - 1
         if(jmin .le. jmax)then
            mu = iju(k) - jmin
            do j=jmin,jmax
               sum = sum + u(j) * tmp(ju(mu+j))
            end do
         end if
         tmp(k) = -sum
         z(c(k)) = -sum
         k = k - 1
      end do
      return
      end
