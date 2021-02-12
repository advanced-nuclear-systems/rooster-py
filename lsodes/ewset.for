
! this subroutine sets the error weight vector ewt according to ewt(i) = rtol(i)*abs(ycur(i)) + atol(i),  i = 1,...,n,
! with the subscript on rtol and/or atol possibly replaced by 1 above, depending on the value of itol.

      subroutine ewset(n, itol, rtol, atol, ycur, ewt)
      integer n, itol
      integer i
      double precision rtol, atol, ycur, ewt
      dimension rtol(1), atol(1), ycur(n), ewt(n)

      if(itol .eq. 1)then
         do i = 1,n
            ewt(i) = rtol(1)*dabs(ycur(i)) + atol(1)
         end do
      else if(itol .eq. 2)then
         do i = 1,n
            ewt(i) = rtol(1)*dabs(ycur(i)) + atol(i)
         end do
      else if(itol .eq. 3)then
         do i = 1,n
            ewt(i) = rtol(i)*dabs(ycur(i)) + atol(1)
         end do
      else if(itol .eq. 4)then
          do i = 1,n
             ewt(i) = rtol(i)*dabs(ycur(i)) + atol(i)
          end do
      end if
      RETURN
      end
