
! purge inactive elements and do mass elimination

      subroutine mdp(k, ek, tail, v, l, head, last, next, mark)
      integer ek, tail, v(1), l(1), head(1), last(1), next(1), mark(1), tag, free, li, vi, lvi, evi, s, ls, es, ilp, ilpmax

!     initialize tag
      tag = mark(ek)

!     for each vertex vi in ek
      li = ek
      ilpmax = last(ek)
      if(ilpmax .gt. 0)then
         do ilp=1,ilpmax
            i = li
            li = l(i)
            vi = v(li)
            
!           remove vi from degree list
            if(last(vi) .gt. 0)then
               next(last(vi)) = next(vi)
            else if(last(vi) .lt. 0)then
               head(-last(vi)) = next(vi)
            end if
            if(next(vi) .gt. 0) last(next(vi)) = last(vi)

!           remove inactive items from element list of vi
            ls = vi
            do while(.true.)
               s = ls
               ls = l(s)
               if(ls .eq. 0)exit
               es = v(ls)
               if(mark(es) .ge. tag)then
                  free = ls
                  l(s) = l(ls)
                  ls = s
               end if
            end do
            
!           if vi is interior vertex, then remove from list and eliminate
            lvi = l(vi)
            if(lvi .eq. 0)then
               l(i) = l(li)
               li = i
               
               k = k+1
               next(vi) = -k
               last(ek) = last(ek) - 1
            else
!              else classify vertex vi
               if(l(lvi) .eq. 0)then
                  evi = v(lvi)
                  if(next(evi) .lt. 0)then
                     if(mark(evi).lt.0)then
!                       if vi is duplicate vertex, then mark as such and adjust overlap count for corresponding element
                        last(vi) = 0
                        mark(evi) = mark(evi) - 1
                     else
!                       else if vi is prototype vertex, then mark as such, initialize overlap count for corresponding element, and move vi to end of boundary list
                        last(vi) = evi
                        mark(evi) = -1
                        l(tail) = li
                        tail = li
                        l(i) = l(li)
                        li = i
                     end if
                  end if
               else
!                 else mark vi to compute degree
                  last(vi) = -ek
               end if
            
!             insert ek in element list of vi
              v(free) = ek
              l(free) = l(vi)
              l(vi) = free
            end if
         end do
      end if
!     terminate boundary list
      l(tail) = 0
      return
      end
