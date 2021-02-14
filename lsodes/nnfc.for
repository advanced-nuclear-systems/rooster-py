
! numerical ldu-factorization of sparse nonsymmetric matrix and solution of system of linear equations (compressed pointer storage)

      subroutine nnfc(n, r, c, ic, ia, ja, a, z, b, lmax, il, jl, ijl, l, d, umax, iu, ju, iju, u, row, tmp, irl, jrl, flag)
      integer rk, umax
      integer r(1), c(1), ic(1), ia(1), ja(1), il(1), jl(1), ijl(1)
      integer iu(1), ju(1), iju(1), irl(1), jrl(1), flag
      double precision a(1), l(1), d(1), u(1), z(1), b(1), row(1)
      double precision tmp(1), lki, sum, dk

!     initialize pointers and test storage
      if(il(n+1)-1 .gt. lmax)then
!        error.. insufficient storage for l
         flag = 4*n + 1
         return
      end if
      if(iu(n+1)-1 .gt. umax)then
!        error.. insufficient storage for u
         flag = 7*n + 1
         return
      end if
      do k=1,n
         irl(k) = il(k)
         jrl(k) = 0
      end do

!     for each row
      do k=1,n
!        reverse jrl and zero row where kth row of l will fill in
         row(k) = 0
         i1 = 0
         if(jrl(k) .ne. 0)then
            i = jrl(k)
            do while(.true.)
               i2 = jrl(i)
               jrl(i) = i1
               i1 = i
               row(i) = 0
               i = i2
               if(i .eq. 0) exit
            end do
         end if
!        set row to zero where u will fill in
         jmin = iju(k)
         jmax = jmin + iu(k+1) - iu(k) - 1
         if(jmin .le. jmax)then
            do j=jmin,jmax
               row(ju(j)) = 0
            end do
         end if
!        place kth row of a in row
         rk = r(k)
         jmin = ia(rk)
         jmax = ia(rk+1) - 1
         do j=jmin,jmax
            row(ic(ja(j))) = a(j)
         end do
!        initialize sum, and link through jrl
         sum = b(rk)
         i = i1
         if(i .ne. 0)then
            do while(.true.)
!              assign the kth row of l and adjust row, sum
               lki = -row(i)
!              if l is not required, then comment out the following line
               l(irl(i)) = -lki
               sum = sum + lki * tmp(i)
               jmin = iu(i)
               jmax = iu(i+1) - 1
               if(jmin .le. jmax)then
                  mu = iju(i) - jmin
                  do j=jmin,jmax
                     row(ju(mu+j)) = row(ju(mu+j)) + lki * u(j)
                  end do
               end if
               i = jrl(i)
               if(i .eq. 0) exit
            end do
         end if

!        assign kth row of u and diagonal d, set tmp(k)
         if(row(k) .eq. 0.0d0)then
!           error.. zero pivot
            flag = 8*n + k
            return
         end if
         dk = 1.0d0 / row(k)
         d(k) = dk
         tmp(k) = sum * dk
         if(k .ne. n)then
            jmin = iu(k)
            jmax = iu(k+1) - 1
            if(jmin .le. jmax)then
               mu = iju(k) - jmin
               do j=jmin,jmax
                  u(j) = row(ju(mu+j)) * dk
               end do
            end if
            
!           update irl and jrl, keeping jrl in decreasing order
            i = i1
            if(i .ne. 0)then
               do while(.true.)
                  irl(i) = irl(i) + 1
                  i1 = jrl(i)
                  if(irl(i) .lt. il(i+1))then
                     ijlb = irl(i) - il(i) + ijl(i)
                     j = jl(ijlb)
                     do while(i .le. jrl(j))
                        j = jrl(j)
                     end do
                     jrl(i) = jrl(j)
                     jrl(j) = i
                  end if
                  i = i1
                  if(i .eq. 0) exit
               end do
            end if
            if(irl(k) .lt. il(k+1))then
               j = jl(ijl(k))
               jrl(k) = jrl(j)
               jrl(j) = k
            end if
         end if
      end do

!     solve  ux = tmp  by back substitution
      k = n
      do i=1,n
         sum =  tmp(k)
         jmin = iu(k)
         jmax = iu(k+1) - 1
         if(jmin .le. jmax)then
            mu = iju(k) - jmin
            do j=jmin,jmax
               sum = sum - u(j) * tmp(ju(mu+j))
            end do
         end if
         tmp(k) =  sum
         z(c(k)) =  sum
      end do
      k = k-1
      flag = 0
      return

      end
