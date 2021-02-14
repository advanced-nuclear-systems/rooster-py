
! reorders rows of a, leaving row order unchanged

      subroutine nroc (n, ic, ia, ja, a, jar, ar, p, flag)
      integer  ic(1), ia(n), ja(1), jar(1), p(1), flag
      double precision  a(1), ar(1)

      do k=1,n
         jmin = ia(k)
         jmax = ia(k+1) - 1
         if(jmin .le. jmax)then
            p(n+1) = n + 1
            do j=jmin,jmax
               newj = ic(ja(j))
               i = n + 1
               do while(.true.)
                  if(p(i) .ge. newj) exit
                  i = p(i)
               end do
               if(p(i) .eq. newj)then
                  flag = n + k
                  return
               end if
               p(newj) = p(i)
               p(i) = newj
               jar(newj) = ja(j)
               ar(newj) = a(j)
            end do
            i = n + 1
            do j=jmin,jmax
               i = p(i)
               ja(j) = jar(i)
               a(j) = ar(i)
            end do
         end if
      end do
      flag = 0
      return
      end
