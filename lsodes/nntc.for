
!     numeric solution of the transpose of a sparse nonsymmetric system of linear equations given lu-factorization (compressed pointer storage)

      subroutine nntc(n, r, c, il, jl, ijl, l, d, iu, ju, iju, u, z, b, tmp)
      integer r(1), c(1), il(n), jl(1), ijl(1), iu(n), ju(1), iju(1)
      double precision l(1), d(1), u(1), b(1), z(1), tmp(1), tmpk,sum

!     set tmp to reordered b
      do k=1,n
         tmp(k) = b(c(k))
      end do
!     solve  ut y = b  by forward substitution
      do k=1,n
         jmin = iu(k)
         jmax = iu(k+1) - 1
         tmpk = -tmp(k)
         if(jmin .le. jmax)then
            mu = iju(k) - jmin
            do j=jmin,jmax
               tmp(ju(mu+j)) = tmp(ju(mu+j)) + tmpk * u(j)
            end do
         end if
      end do
!     solve  lt x = y  by back substitution
      k = n
      do i=1,n
         sum = -tmp(k)
         jmin = il(k)
         jmax = il(k+1) - 1
         if(jmin .le. jmax)then
            ml = ijl(k) - jmin
            do j=jmin,jmax
               sum = sum + l(j) * tmp(jl(ml+j))
            end do
         end if
         tmp(k) = -sum * d(k)
         z(r(k)) = tmp(k)
         k = k - 1
      end do
      return
      end
