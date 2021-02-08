!
! cfode is called by the integrator routine to set coefficients needed there.
! the coefficients for the current method, as given by the value of meth, are set for all orders and saved.
! the maximum order assumed here is 12 if meth = 1 and 5 if meth = 2.
! (a smaller value of the maximum order is also allowed.)
! cfode is called once at the beginning of the problem, and is not called again unless and until meth is changed.
!
! the elco array contains the basic method coefficients.
! the coefficients el(i), 1 .le. i .le. nq+1, for the method of order nq are stored in elco(i,nq).
! they are given by a genetrating polynomial, i.e., l(x) = el(1) + el(2)*x + ... + el(nq+1)*x**nq.
! for the implicit adams methods, l(x) is given by dl/dx = (x+1)*(x+2)*...*(x+nq-1)/factorial(nq-1),    l(-1) = 0.
! for the bdf methods, l(x) is given by l(x) = (x+1)*(x+2)* ... *(x+nq)/k, where k = factorial(nq)*(1 + 1/2 + ... + 1/nq).
!
! the tesco array contains test constants used for the local error test and the selection of step size and/or order.
! at order nq, tesco(k,nq) is used for the selection of step size at order nq - 1 if k = 1, at order nq if k = 2, and at order nq + 1 if k = 3.
!
      subroutine cfode(meth, elco, tesco)

      integer meth
      integer i, ib, nq, nqm1, nqp1
      double precision elco, tesco
      double precision agamq, fnq, fnqm1, pc, pint, ragq, rqfac, rq1fac, tsign, xpin
      dimension elco(13,12), tesco(3,12)
      dimension pc(12)

      if(meth .eq. 1)then
         elco(1,1) = 1.0d0
         elco(2,1) = 1.0d0
         tesco(1,1) = 0.0d0
         tesco(2,1) = 2.0d0
         tesco(1,2) = 1.0d0
         tesco(3,12) = 0.0d0
         pc(1) = 1.0d0
         rqfac = 1.0d0
         do nq = 2,12
!           the pc array will contain the coefficients of the polynomial p(x) = (x+1)*(x+2)*...*(x+nq-1). initially, p(x) = 1.
            rq1fac = rqfac
            rqfac = rqfac/dfloat(nq)
            nqm1 = nq - 1
            fnqm1 = dfloat(nqm1)
            nqp1 = nq + 1
!           form coefficients of p(x)*(x+nq-1).
            pc(nq) = 0.0d0
            do ib = 1,nqm1
               i = nqp1 - ib
               pc(i) = pc(i-1) + fnqm1*pc(i)
            end do
            pc(1) = fnqm1*pc(1)
!           compute integral, -1 to 0, of p(x) and x*p(x).
            pint = pc(1)
            xpin = pc(1)/2.0d0
            tsign = 1.0d0
            do i = 2,nq
               tsign = -tsign
               pint = pint + tsign*pc(i)/dfloat(i)
               xpin = xpin + tsign*pc(i)/dfloat(i+1)
            end do
!           store coefficients in elco and tesco.
            elco(1,nq) = pint*rq1fac
            elco(2,nq) = 1.0d0
            do i = 2,nq
               elco(i+1,nq) = rq1fac*pc(i)/dfloat(i)
            end do
            agamq = rqfac*xpin
            ragq = 1.0d0/agamq
            tesco(2,nq) = ragq
            if (nq .lt. 12) tesco(1,nqp1) = ragq*rqfac/dfloat(nqp1)
            tesco(3,nqm1) = ragq
         end do

      else ! meth == 2
         pc(1) = 1.0d0
         rq1fac = 1.0d0
         do nq = 1,5
!           the pc array will contain the coefficients of the polynomial p(x) = (x+1)*(x+2)*...*(x+nq). initially, p(x) = 1.
            fnq = dfloat(nq)
            nqp1 = nq + 1
!           form coefficients of p(x)*(x+nq).
            pc(nqp1) = 0.0d0
            do ib = 1,nq
               i = nq + 2 - ib
               pc(i) = pc(i-1) + fnq*pc(i)
            end do
            pc(1) = fnq*pc(1)
!           store coefficients in elco and tesco.
            do i = 1,nqp1
               elco(i,nq) = pc(i)/pc(2)
            end do
            elco(2,nq) = 1.0d0
            tesco(1,nq) = rq1fac
            tesco(2,nq) = dfloat(nqp1)/elco(1,nq)
            tesco(3,nq) = dfloat(nq+2)/elco(1,nq)
            rq1fac = rq1fac/fnq
         end do
      end if
      return
      end
