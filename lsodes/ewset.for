      subroutine ewset (n, itol, rtol, atol, ycur, ewt)
      integer n, itol
      integer i
      double precision rtol, atol, ycur, ewt
      dimension rtol(1), atol(1), ycur(n), ewt(n)

      if(itol == 10)then
         do i = 1,n
            ewt(i) = rtol(1)*dabs(ycur(i)) + atol(1)
         end do
         return
      else if(itol == 20)then
         do i = 1,n
            ewt(i) = rtol(1)*dabs(ycur(i)) + atol(i)
         end do
         return
      else if(itol == 30)then
         do i = 1,n
            ewt(i) = rtol(i)*dabs(ycur(i)) + atol(1)
         end do
         return
      else if(itol == 40)then
          do i = 1,n
             ewt(i) = rtol(i)*dabs(ycur(i)) + atol(i)
          end do
          return
      end if
      end
