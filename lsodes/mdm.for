      subroutine mdm(vk,tail, v,l, last,next, mark)
      integer  vk, tail, v(1), l(1), last(1), next(1), mark(1), tag, s, ls, vs, es, b, lb, vb, blp, blpmax
      equivalence  (vs, es)

      tag = mark(vk)
      tail = vk
      ls = l(vk)
      s = ls
      do while(s.ne.0)
         ls = l(s)
         vs = v(s)
         if(next(vs).lt.0)then
            lb = l(es)
            blpmax = last(es)
            do blp=1,blpmax
               b = lb
               lb = l(b)
               vb = v(b)
               if(mark(vb).lt.tag)then
                  mark(vb) = tag
                  l(tail) = b
                  tail = b
               end if
            end do
            mark(es) = tag
         else
            mark(vs) = tag
            l(tail) = s
            tail = s
         end if
         s = ls
      end do

      l(tail) = 0

      return
      end
