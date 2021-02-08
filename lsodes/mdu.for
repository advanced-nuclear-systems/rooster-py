
!     update degrees of uneliminated vertices in ek

      subroutine mdu(ek,dmin, v,l, head, last, next, mark)
      integer  ek, dmin,  v(1), l(1),  head(1), last(1), next(1), mark(1),  tag, vi,evi,dvi, s, vs, es, b, vb, ilp, ilpmax, blp, blpmax
      equivalence  (vs, es)

      tag = mark(ek) - last(ek)

      i = ek
      ilpmax = last(ek)
      if (ilpmax.le.0) return
      do ilp=1,ilpmax
         i = l(i)
         vi = v(i)
         if (last(vi))  1, 10, 8
         
   1       tag = tag + 1
           dvi = last(ek)
         
           s = l(vi)
   2       s = l(s)
           if (s.eq.0)  go to 9
             vs = v(s)
             if (next(vs).lt.0)  go to 3
         
               mark(vs) = tag
               dvi = dvi + 1
               go to 5
         
   3           if (mark(es).lt.0)  go to 6
         
               b = es
               blpmax = last(es)
               do 4 blp=1,blpmax
                 b = l(b)
                 vb = v(b)
         
                 if (mark(vb).ge.tag)  go to 4
                   mark(vb) = tag
                   dvi = dvi + 1
   4             continue
         
   5         go to 2
         
   6       last(vi) = 0
           mark(es) = mark(es) - 1
   7       s = l(s)
           if (s.eq.0)  go to 10
             es = v(s)
             if (mark(es).lt.0)  mark(es) = mark(es) - 1
             go to 7
         
   8       evi = last(vi)
           dvi = last(ek) + last(evi) + mark(evi)
           mark(evi) = 0
         
   9     next(vi) = head(dvi)
         head(dvi) = vi
         last(vi) = -dvi
         if (next(vi).gt.0)  last(next(vi)) = vi
         if (dvi.lt.dmin)  dmin = dvi
         
  10     continue
      end do
      end
