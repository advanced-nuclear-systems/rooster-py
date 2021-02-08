      subroutine mdi(n, ia, ja, max, v, l, head, last, next, mark, tag, flag)

      integer ia(n), ja(1), v(1), l(1), head(1), last(1), next(1), mark(1), tag, flag, sfs, vi, dvi, vj

      do vi=1,n
         mark(vi) = 1
         l(vi) = 0
         head(vi) = 0
      end do
      sfs = n+1
      do vi=1,n
        jmin = ia(vi)
        jmax = ia(vi+1) - 1
        if(jmin.le.jmax)then
           do j=jmin,jmax
              vj = ja(j)
              if(vj < vi)then
                 lvk = vi
                 kmax = mark(vi) - 1
                 if(kmax .ne. 0)then
                    do k=1,kmax
                       lvk = l(lvk)
                       if(v(lvk).eq.vj) exit
                    end do
                 end if
              else if(vj > vi)then
                 if(sfs.ge.max)then
                    flag = 9*n + vi
                    return
                 end if
                 
                 mark(vi) = mark(vi) + 1
                 v(sfs) = vj
                 l(sfs) = l(vi)
                 l(vi) = sfs
                 sfs = sfs+1
                 
                 mark(vj) = mark(vj) + 1
                 v(sfs) = vi
                 l(sfs) = l(vj)
                 l(vj) = sfs
                 sfs = sfs+1
              end if
           end do
        end if
      end do

      do vi=1,n
         dvi = mark(vi)
         next(vi) = head(dvi)
         head(dvi) = vi
         last(vi) = -dvi
         nextvi = next(vi)
         if (nextvi.gt.0) last(nextvi) = vi
         mark(vi) = tag
      end do
      return
      end
