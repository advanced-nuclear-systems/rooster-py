      subroutine solve_eigenvalue_problem(n, m, a, b)
       
      implicit none
      integer n, m, i, j
      real*8 a(:,:)
      real*8 b(:)

      do i = 1, n
         write(*,*)(a(i,j),j=1,m)
      end do
      do i = 1, m
         write(*,*)b(i)
      end do
      a(1,1) = 200.
      b(1) = 200.
      end
