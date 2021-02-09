      subroutine ewset (n, itol, rtol, atol, ycur, ewt)
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
