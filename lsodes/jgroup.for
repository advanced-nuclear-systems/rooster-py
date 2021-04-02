
! this subroutine constructs groupings of the column indices of the jacobian matrix, 
! used in the numerical evaluation of the jacobian by finite differences.
! input..
! n      = the order of the matrix.
! ia,ja  = sparse structure descriptors of the matrix by rows.
! maxg   = length of available storate in the igp array.
!
! output..
! ngrp   = number of groups.
! jgp    = array of length n containing the column indices by groups.
! igp    = pointer array of length ngrp + 1 to the locations in jgp of the beginning of each group.
! ier    = error indicator.  ier = 0 if no error occurred, or 1 if maxg was insufficient.
!
! incl and jdone are working arrays of length n.

      subroutine jgroup (n,ia,ja,maxg,ngrp,igp,jgp,incl,jdone,ier)
      integer n, ia, ja, maxg, ngrp, igp, jgp, incl, jdone, ier
      dimension ia(1), ja(1), igp(1), jgp(n), incl(n), jdone(n)
      integer i, j, k, kmin, kmax, ncol, ng
      logical exit1

      exit1 = .false.
      ier = 0
      do j = 1,n
         jdone(j) = 0
      end do
      ncol = 1
      do ng = 1,maxg
          igp(ng) = ncol
          do i = 1,n
             incl(i) = 0
          end do
          do j = 1,n
!            reject column j if it is already in a group.
             if(jdone(j) .ne. 1)then
                kmin = ia(j)
                kmax = ia(j+1) - 1
                do k = kmin,kmax
!                  reject column j if it overlaps any column already in this group.
                   i = ja(k)
                   if(incl(i) .eq. 1)then
                      exit1 = .true.
                      exit
                   end if
                end do
                if(exit1)then
                   exit1 = .false.
                else
!                  accept column j into group ng.
                   jgp(ncol) = j
                   ncol = ncol + 1
                   jdone(j) = 1
                   do k = kmin,kmax
                      i = ja(k)
                      incl(i) = 1
                   end do
                end if
             end if
          end do
!         stop if this group is empty (grouping is complete).
          if(ncol .eq. igp(ng))then
             ngrp = ng - 1
!             write(*,*)'ngrp ',ngrp
!             write(*,*)'igp ',(igp(jj),jj=1,ng)
!             write(*,*)'jgp ',jgp
             RETURN
          end if
      end do
!     error return if not all columns were chosen (maxg too small).
      if(ncol .gt. n)then
         ng = maxg
         ngrp = ng - 1
         RETURN
      end if
      ier = 1
      RETURN
      end
