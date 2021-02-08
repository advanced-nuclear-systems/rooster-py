      subroutine sro(n, ip, ia, ja, a, q, r, dflag)

      integer ip(1), ia(n), ja(1),  q(1), r(1)
      double precision a(1),  ak
      logical dflag

      do i=1,n
         q(i) = 0
      end do

      do i=1,n
        jmin = ia(i)
        jmax = ia(i+1) - 1
        if(jmin .le. jmax)then
           do j=jmin,jmax
             k = ja(j)
             if (ip(k).lt.ip(i))  ja(j) = i
             if (ip(k).ge.ip(i))  k = i
             r(j) = k
             q(k) = q(k) + 1
           end do
        end if
      end do
      do i=1,n
         ia(i+1) = ia(i) + q(i)
         q(i) = ia(i+1)
      end do

      ilast = 0
      jmin = ia(1)
      jmax = ia(n+1) - 1
      j = jmax
      do jdummy=jmin,jmax
        i = r(j)
        if(.not.dflag .or. ja(j).ne.i .or. i.eq.ilast)then
           q(i) = q(i) - 1
           r(j) = q(i)
        else
           r(j) = ia(i)
           ilast = i
        end if
        j = j - 1
      end do


      do j = jmin, jmax
         do while(r(j) .ne. j)
            k = r(j)
            r(j) = r(k)
            r(k) = k
            jak = ja(k)
            ja(k) = ja(j)
            ja(j) = jak
            ak = a(k)
            a(k) = a(j)
            a(j) = ak
        end do
      end do

      return
      end
