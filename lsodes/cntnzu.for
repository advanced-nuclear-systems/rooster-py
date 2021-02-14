
! this routine counts the number of nonzero elements in the strict upper triangle of the matrix m + m(transpose), 
! where the sparsity structure of m is given by pointer arrays ia and ja.
! this is needed to compute the storage requirements for the sparse matrix reordering operation in odrv.

! if matrix elements are held in a(k), k=1,...,nz,
! ja(k) holds the column number of the element held in a(k).
! ia(i) contains the address of the first element of row i and ia(n+1) contains the address of the first unused element in a.

      subroutine cntnzu(n, ia, ja, nzsut)
      integer n, nzsut
      integer ia(n), ja(n)
      integer ii, jj, j, jmin, jmax, k, kmin, kmax, num
      logical flag

!     counter of nonzeros
      num = 0
!     loop over columns
      do ii = 1,n
!        addresses of the first and last elements of column ii
         jmin = ia(ii) 
         jmax = ia(ii+1) - 1
         if(jmin .le. jmax)then
            do j = jmin,jmax
               if(ja(j) .lt. ii)then
                  jj =ja(j)
                  kmin = ia(jj)
                  kmax = ia(jj+1) - 1
                  if(kmin .le. kmax)then
                     flag = .false.
                     do k = kmin,kmax
                        if(ja(k) .eq. ii)then
                           flag = .true.
                           exit
                        end if
                     end do
                  end if
                  if(.not. flag) num = num + 1
               else if(ja(j) .gt. ii)then
                  num = num + 1
               end if
            end do
         end if
      end do
      nzsut = num
      write(*,*)'nzsut ', nzsut
      return
      end
