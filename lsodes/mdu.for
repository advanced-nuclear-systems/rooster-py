
!     update degrees of uneliminated vertices in ek

      subroutine mdu(ek, dmin, v, l, head, last, next, mark)
      integer ek, dmin, v(1), l(1), head(1), last(1), next(1), mark(1), tag, vi, evi, dvi, s, ves, b, vb, ilp, ilpmax, blp, blpmax

!     initialize tag
      tag = mark(ek) - last(ek)

!     for each vertex vi in ek
      i = ek
      ilpmax = last(ek)
      if(ilpmax.gt.0)then
         do ilp=1,ilpmax
            i = l(i)
            vi = v(i)
            if(last(vi) .lt. 0)then
!              if vi neither prototype nor duplicate vertex, then merge elements to compute degree
               tag = tag + 1
               dvi = last(ek)
               
!              for each vertex/element ves in element list of vi
               s = l(vi)
               do while(.true.)
                  s = l(s)
                  if(s .eq. 0)then
!                    insert vi in appropriate degree list
                     next(vi) = head(dvi)
                     head(dvi) = vi
                     last(vi) = -dvi
                     if(next(vi).gt.0)  last(next(vi)) = vi
                     if(dvi.lt.dmin)  dmin = dvi
                     exit
                  end if
                  ves = v(s)
                  if(next(ves).ge.0)then
!                    if ves is uneliminated vertex, then tag and adjust degree
                     mark(ves) = tag
                     dvi = dvi + 1
                  else
!                    if ves is active element, then expand check for outmatched vertex
                     if(mark(ves).lt.0)then
!                       else if vi is outmatched vertex, then adjust overlaps but do not compute degree
                        last(vi) = 0
                        mark(ves) = mark(ves) - 1
                        do while(.true.)
                           s = l(s)
                           if(s .eq. 0) exit
                           ves = v(s)
                           if(mark(ves).lt.0) mark(ves) = mark(ves) - 1
                        end do
                        exit
                     else
!                       for each vertex vb in es
                        b = ves
                        blpmax = last(ves)
                        do blp=1,blpmax
                           b = l(b)
                           vb = v(b)
!                          if vb is untagged, then tag and adjust degree
                           if(mark(vb) .lt. tag)then
                              mark(vb) = tag
                              dvi = dvi + 1
                           end if
                        end do
                     end if
                  end if
               end do
            else if(last(vi) .gt. 0)then
!              else if vi is prototype vertex, then calculate degree by inclusion/exclusion and reset overlap count
               evi = last(vi)
               dvi = last(ek) + last(evi) + mark(evi)
               mark(evi) = 0

!              insert vi in appropriate degree list
               next(vi) = head(dvi)
               head(dvi) = vi
               last(vi) = -dvi
               if(next(vi).gt.0)  last(next(vi)) = vi
               if(dvi.lt.dmin)  dmin = dvi
            end if
         end do
      end if
      return
      end
