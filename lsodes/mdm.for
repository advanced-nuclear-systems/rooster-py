
!  mdm -- form element from uneliminated neighbors of vk

      subroutine mdm(vk,tail, v,l, last,next, mark)
      integer  vk, tail, v(1), l(1), last(1), next(1), mark(1), tag, s, ls, vs, es, b, lb, vb, blp, blpmax
      equivalence  (vs, es)

!     initialize tag and list of uneliminated neighbors
      tag = mark(vk)
      tail = vk
!     for each vertex/element vs/es in element list of vk
      ls = l(vk)
      s = ls
      do while(s.ne.0)
         ls = l(s)
         vs = v(s)
         if(next(vs).lt.0)then
!        if es is active element, then for each vertex vb in boundary list of element es
            lb = l(es)
            blpmax = last(es)
            do blp=1,blpmax
               b = lb
               lb = l(b)
               vb = v(b)
!              if vb is untagged vertex, then tag and append to list of uneliminated neighbors
               if(mark(vb).lt.tag)then
                  mark(vb) = tag
                  l(tail) = b
                  tail = b
               end if
            end do
!           mark es inactive
            mark(es) = tag
         else
!           if vs is uneliminated vertex, then tag and append to list of uneliminated neighbors
            mark(vs) = tag
            l(tail) = s
            tail = s
         end if
         s = ls
      end do

!     terminate list of uneliminated neighbors
      l(tail) = 0

      return
      end
