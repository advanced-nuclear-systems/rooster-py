!
!     md finds a minimum degree ordering of the rows and columns of a general sparse matrix m stored in (ia,ja,a) format.
!     when the structure of m is nonsymmetric, the ordering is that obtained for the symmetric matrix  m + m-transpose.
!
!     additional parameters
!
!     max  - declared dimension of the one-dimensional arrays v and l.
!            max must be at least  n+2k,  where k is the number of nonzeroes in the strict upper triangle of m + m-transpose
!     v    - integer one-dimensional work array.  dimension = max
!     l    - integer one-dimensional work array.  dimension = max
!     head - integer one-dimensional work array.  dimension = n
!     last - integer one-dimensional array used to return the permutation of the rows and columns of m corresponding to the minimum degree ordering.  
!            dimension = n
!     next - integer one-dimensional array used to return the inverse of the permutation returned in last.  dimension = n
!     mark - integer one-dimensional work array (may be the same as v). dimension = n
!     flag - integer error flag.  values and their meanings are -
!            0 no errors detected
!            9n+k  insufficient storage in md

      subroutine md(n, ia, ja, max, v, l, head, last, next, mark, flag)
      integer ia(1), ja(1), v(1), l(1), head(1), last(1), next(1), mark(1), flag, tag, dmin, vk, ek, tail
      equivalence (vk,ek)

!     initialization
      tag = 0
      call mdi(n, ia,ja, max,v,l, head,last,next, mark,tag, flag)
      if (flag.ne.0)  return

      k = 0
      dmin = 1

      do while(k.lt.n)
!        search for vertex of minimum degree
         do while(head(dmin).le.0)
           dmin = dmin + 1
         end do
!        remove vertex vk of minimum degree from degree list
         vk = head(dmin)
         head(dmin) = next(vk)
         if (head(dmin).gt.0)  last(head(dmin)) = -dmin

!        number vertex vk, adjust tag, and tag vk
         k = k + 1
         next(vk) = -k
         last(ek) = dmin - 1
         tag = tag + last(ek)
         mark(vk) = tag

!        form element ek from uneliminated neighbors of vk
         call mdm(vk, tail, v, l, last, next, mark)
!        purge inactive elements and do mass elimination
         call mdp(k, ek, tail, v, l, head, last, next, mark)
!        update degrees of uneliminated vertices in ek
         call mdu(ek, dmin, v, l, head, last, next, mark)
      end do

!     generate inverse permutation from permutation
      do k=1,n
         next(k) = -next(k)
         last(next(k)) = k
      end do
      return
      end
