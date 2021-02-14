
!  symmetric reordering of sparse symmetric matrix.
!
!  the nonzero entries of the matrix m are assumed to be stored symmetrically in (ia,ja,a) format (i.e., not both m(i,j) and m(j,i) are stored if i ne j).
!
!  sro does not rearrange the order of the rows, but does move nonzeroes from one row to another to ensure that if m(i,j) will be
!  in the upper triangle of m with respect to the new ordering, then m(i,j) is stored in row i (and thus m(j,i) is not stored),  whereas
!  if m(i,j) will be in the strict lower triangle of m, then m(j,i) is stored in row j (and thus m(i,j) is not stored).
!
!  additional parameters
!  q - integer one-dimensional work array.  dimension = n
!  r - integer one-dimensional work array.  dimension = number of nonzero entries in the upper triangle of m
!  dflag - logical variable.  if dflag = .true., then store nonzero diagonal elements at the beginning of the row

      subroutine sro(n, ip, ia, ja, a, q, r, dflag)

      integer ip(1), ia(n), ja(1),  q(1), r(1)
      double precision a(1),  ak
      logical dflag

!     phase 1 -- find row in which to store each nonzero initialize count of nonzeroes to be stored in each row
      do i=1,n
         q(i) = 0
      end do

!     for each nonzero element a(j)
      do i=1,n
        jmin = ia(i)
        jmax = ia(i+1) - 1
        if(jmin .le. jmax)then
           do j=jmin,jmax
!            find row (=r(j)) and column (=ja(j)) in which to store a(j) ...
             k = ja(j)
             if (ip(k).lt.ip(i))  ja(j) = i
             if (ip(k).ge.ip(i))  k = i
             r(j) = k
!            ... and increment count of nonzeroes (=q(r(j)) in that row
             q(k) = q(k) + 1
           end do
        end if
      end do

!     phase 2 -- find new ia and permutation to apply to (ja,a) determine pointers to delimit rows in permuted (ja,a)
      do i=1,n
         ia(i+1) = ia(i) + q(i)
         q(i) = ia(i+1)
      end do

!     determine where each (ja(j),a(j)) is stored in permuted (ja,a) for each nonzero element (in reverse order)
      ilast = 0
      jmin = ia(1)
      jmax = ia(n+1) - 1
      j = jmax
      do jdummy=jmin,jmax
        i = r(j)
        if(.not.dflag .or. ja(j).ne.i .or. i.eq.ilast)then
!          put (off-diagonal) nonzero in last unused location in row
           q(i) = q(i) - 1
           r(j) = q(i)
        else
!          if dflag, then put diagonal nonzero at beginning of row
           r(j) = ia(i)
           ilast = i
        end if
        j = j - 1
      end do

!     phase 3: permute (ja,a) to upper triangular form (wrt new ordering)
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
