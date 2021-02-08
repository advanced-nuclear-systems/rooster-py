      subroutine jgroup (n,ia,ja,maxg,ngrp,igp,jgp,incl,jdone,ier)
      integer n, ia, ja, maxg, ngrp, igp, jgp, incl, jdone, ier
      dimension ia(1), ja(1), igp(1), jgp(n), incl(n), jdone(n)
      integer i, j, k, kmin, kmax, ncol, ng

      ier = 0
      do j = 1,n
         jdone(j) = 0
      end do
      ncol = 1
      do ng = 1,maxg
          igp(ng) = ncol
          do i = 1,n
             incl(i) = 0
          end do
          do j = 1,n
             if(jdone(j) .ne. 1)then
                kmin = ia(j)
                kmax = ia(j+1) - 1
                do k = kmin,kmax
                   i = ja(k)
                   if(incl(i) .eq. 1) exit
                end do
                jgp(ncol) = j
                ncol = ncol + 1
                jdone(j) = 1
                do k = kmin,kmax
                   i = ja(k)
                   incl(i) = 1
                end do
             end if
          end do
          if(ncol .eq. igp(ng))then
             ngrp = ng - 1
             return
          end if
      end do
      if(ncol .gt. n)then
         ng = maxg
         ngrp = ng - 1
         return
      end if
      ier = 1
      return
      end
