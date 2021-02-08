      subroutine cntnzu (n, ia, ja, nzsut)
      integer n, ia, ja, nzsut
      dimension ia(n), ja(n)
      integer ii, jj, j, jmin, jmax, k, kmin, kmax, num

      num = 0
      do ii = 1,n
         jmin = ia(ii)
         jmax = ia(ii+1) - 1
         if(jmin .le. jmax)then
            do 40 j = jmin,jmax
              if (ja(j) - ii) 10, 40, 30
 10           jj =ja(j)
              kmin = ia(jj)
              kmax = ia(jj+1) - 1
              if (kmin .gt. kmax) go to 30
              do k = kmin,kmax
                 if (ja(k) .eq. ii) go to 40
              end do
 30           num = num + 1
 40           continue
         end if
      end do
      nzsut = num
      return
      end
