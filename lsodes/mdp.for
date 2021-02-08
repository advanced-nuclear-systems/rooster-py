      subroutine mdp(k,ek,tail, v,l, head,last,next, mark)
      integer  ek, tail,  v(1), l(1),  head(1), last(1), next(1), mark(1),  tag, free, li,vi,lvi,evi, s,ls,es, ilp,ilpmax

      tag = mark(ek)

      li = ek
      ilpmax = last(ek)
      if (ilpmax.le.0)  go to 12
      do 11 ilp=1,ilpmax
        i = li
        li = l(i)
        vi = v(li)

        if (last(vi).eq.0)  go to 3
          if (last(vi).gt.0)  go to 1
            head(-last(vi)) = next(vi)
            go to 2
   1        next(last(vi)) = next(vi)
   2      if (next(vi).gt.0)  last(next(vi)) = last(vi)

   3    ls = vi
   4    s = ls
        ls = l(s)
        if (ls.eq.0)  go to 6
          es = v(ls)
          if (mark(es).lt.tag)  go to 5
            free = ls
            l(s) = l(ls)
            ls = s
   5      go to 4

   6    lvi = l(vi)
        if (lvi.ne.0)  go to 7
          l(i) = l(li)
          li = i

          k = k+1
          next(vi) = -k
          last(ek) = last(ek) - 1
          go to 11

   7      if (l(lvi).ne.0)  go to 9
            evi = v(lvi)
            if (next(evi).ge.0)  go to 9
              if (mark(evi).lt.0)  go to 8
                last(vi) = evi
                mark(evi) = -1
                l(tail) = li
                tail = li
                l(i) = l(li)
                li = i
                go to 10

   8            last(vi) = 0
                mark(evi) = mark(evi) - 1
                go to 10

   9            last(vi) = -ek

  10      v(free) = ek
          l(free) = l(vi)
          l(vi) = free
  11    continue

  12  l(tail) = 0

      return
      end
