!
! stode performs one step of the integration of an initial value problem for a system of ordinary differential equations.
! note.. stode is independent of the value of the iteration method indicator miter, when this is .ne. 0, and hence is independent
! of the type of chord method used, or the jacobian structure.
!
! communication with stode is done with the following variables..
!
! n      = the number of first-order differential equations.
! neq    = integer array containing problem size in neq(1), and passed as the neq argument in all calls to f and jac.
! y      = an array of length .ge. n used as the y argument in all calls to f and jac.
! yh     = an nnyh by maxord+1 array containing the dependent variables and their approximate scaled derivatives, where.  
!          yh(i,j+1) contains the approximate j-th derivative of y(i), scaled by h**j/factorial(j) (j = 0,1,...,nq).  
!          on entry for the first step, the first two columns of yh must be set from the initial values.
! nnyh   = a constant integer .ge. n, the first dimension of yh.
! yh1    = a one-dimensional array occupying the same space as yh.
! ewt    = an array of length n containing multiplicative weights for local error measurements.  
!          local errors in y(i) are compared to 1.0/ewt(i) in various error tests.
! savf   = an array of working storage, of length n. also used for input of yh(*,maxord+2) when jstart = -1
!          and maxord .lt. the current order nq.
! acor   = a work array of length n, used for the accumulated corrections.  
!          on a successful return, acor(i) contains the estimated one-step local error in y(i).
! wm,iwm = real and integer work arrays associated with matrix operations in chord iteration (miter .ne. 0).
! fjac   = name of routine to evaluate and preprocess jacobian matrix and p = i - h*el0*jac, if a chord method is being used.
! slss   = name of routine to solve linear system in chord iteration.
! h      = the step size to be attempted on the next step.
!          h is altered by the error control algorithm during the problem.  
!          h can be either positive or negative, but its sign must remain constant throughout the problem.
! hmin   = the minimum absolute value of the step size h to be used.
! hmxi   = inverse of the maximum absolute value of h to be used. hmxi = 0.0 is allowed and corresponds to an infinite hmax.
!          hmin and hmxi may be changed at any time, but will not take effect until the next change of h is considered.
! tn     = the independent variable. tn is updated on each step taken.
! jstart = an integer used for input only, with the following values and meanings..
!               0  perform the first step.
!           .gt.0  take a new step continuing from the last.
!              -1  take the next step with a new value of h, maxord, n, meth, miter, and/or matrix parameters.
!              -2  take the next step with a new value of h, but with other inputs unchanged.
!          on return, jstart is set to 1 to facilitate continuation.
! kflag  = a completion code with the following meanings..
!               0  the step was succesful.
!              -1  the requested error could not be achieved.
!              -2  corrector convergence could not be achieved.
!              -3  fatal error in fjac or slss.
!          a return with kflag = -1 or -2 means either abs(h) = hmin or 10 consecutive failures occurred.
!          on a return with kflag negative, the values of tn and the yh array are as of the beginning of the last step, 
!          and h is the last step size attempted.
! maxord = the maximum order of integration method to be allowed.
! maxcor = the maximum number of corrector iterations allowed.
! meth   = the basic linear multistep method..
!               1  means the implicit adams method.
!               2  means the method based on backward differentiation formulas (bdf-s).
! miter  = the corrector iteration method..
!               0  means functional iteration (no jacobian matrix is involved).
!               1  means chord iteration with a user-supplied sparse jacobian, given by subroutine jac.
!               2  means chord iteration with an internally generated (difference quotient) sparse jacobian
!                  (using ngp extra calls to f per df/dy value, where ngp is an optional output described below.)
!               3  means chord iteration with an internally generated diagonal jacobian approximation. (using 1 extra call to f per df/dy evaluation).
!          if miter = 1 the user must supply a subroutine jac (the name is arbitrary) as described in lsodes.
!          otherwise, a dummy argument can be used.
!
      subroutine stode(neq, y, yh, nnyh, yh1, ewt, savf, acor, wm, iwm)

      integer neq, nnyh, iwm
      integer i, i1, j, jb, m, ncf, newq
      double precision y, yh, yh1, ewt, savf, acor, wm
      double precision dcon, ddn, del, delp, dsm, dup, exdn, exsm, exup, r, rh, rhdn, rhsm, rhup, told, vnorm
      dimension neq(1), y(1), yh(nnyh,1), yh1(1), ewt(1), savf(1), acor(1), wm(1), iwm(1)

      double precision conit, crate, el, elco, hold, rmax, tesco, 
     +   ccmax, el0, h, hmin, hmxi, hu, rc, tn, uround
      integer illin, init, lyh, lewt, lacor, lsavf, lwm, liwm, mxstep, mxhnil, nhnil, ntrep, nslast, nyh,
     +   ialth, ipup, lmax, meo, nqnyh, nslp, 
     +   icf, ierpj, iersl, jcur, jstart, kflag, l, meth, miter, 
     +   maxord, maxcor, msbp, n, nq, nst, nfe, nje, nqu
      common /ls0001/ conit, crate, el(13), elco(13,12), hold, rmax, tesco(3,12),
     +   ccmax, el0, h, hmin, hmxi, hu, rc, tn, uround, 
     +   illin, init, lyh, lewt, lacor, lsavf, lwm, liwm, mxstep, mxhnil, nhnil, ntrep, nslast, nyh, 
     +   ialth, ipup, lmax, meo, nqnyh, nslp, 
     +   icf, ierpj, iersl, jcur, jstart, kflag, l, meth, miter, 
     +   maxord, maxcor, msbp, mxncf, n, nq, nst, nfe, nje, nqu

      logical exit1, exit2, exit3

      kflag = 0
      told = tn
      ncf = 0
      ierpj = 0
      iersl = 0
      jcur = 0
      icf = 0
      delp = 0.0d0

      if(jstart .eq. 0)then
!        on the first call, the order is set to 1, and other variables are initialized.  
         nq = 1
!        rmax is the maximum ratio by which h can be increased in a single step.  
!        it is initially 1.e4 to compensate for the small initial h, but then is normally equal to 10.
!        if a failure occurs (in corrector convergence or error test), rmax is set at 2 for the next increase.
         rmax = 10000.0d0
         rc = 0.0d0
         el0 = 1.0d0
         crate = 0.7d0
         hold = h
         meo = meth
         nslp = 0
         ipup = miter
!        cfode is called to reset the coefficients of the method.
         call cfode(meth, elco, tesco)
!        el vector and related constants are reset at the start of the problem.
         call order_changed()

      else if(jstart .eq. -1)then
!        preliminaries needed when jstart = -1. 
!        if h is to be changed, yh must be rescaled.
!        if h or meth is being changed, ialth is reset to nq + 1 to prevent further changes in h for that many steps.
!        ipup is set to miter to force a matrix update.
         ipup = miter
!        if an order increase is about to be considered (ialth = 1), ialth is reset to 2 to postpone consideration one more step.
         if(ialth .eq. 1) ialth = 2
!        if the caller has changed meth, cfode is called to reset the coefficients of the method.
         if(meth .ne. meo)then
            call cfode(meth, elco, tesco)
            meo = meth
            if(nq .le. maxord)then
!              el vector and related constants are reset when the order nq is changed.
               call order_changed()
            end if
         end if
         if(nq .le. maxord)then
!           if h is being changed, the h ratio rh is reset.  
            if(h .ne. hold)then
               rh = h/hold
               h = hold
            end if
         else
!           if the caller has changed maxord to a value less than the current order nq, nq is reduced to maxord, and a new h chosen accordingly.
            nq = maxord
            call order_changed()
            ddn = vnorm(n, savf, ewt)/tesco(1,nq+1)
            exdn = 1.0d0/dfloat(nq + 1)
            rhdn = 1.0d0/(1.3d0*ddn**exdn + 0.0000013d0)
            rh = dmin1(rhdn,1.0d0)
            if(h .eq. hold)then
               rh = dmax1(rh,hmin/dabs(h))
            else
               rh = dmin1(rh,dabs(h/hold))
               h = hold
            end if
         end if
!        if h is being changed, the h ratio rh is checked against rmax, hmin, and hmxi, and the yh array rescaled.  
         if(h .ne. hold)then
            call step_changed(rh, nyh, yh)
         end if

      else if(jstart .eq. -2 .and. h .ne. hold)then
!        if only h is changed, the h ratio rh is reset, checked against rmax, hmin, and hmxi, and the yh array rescaled.
         rh = h/hold
         h = hold
         call step_changed(rh, nyh, yh)
      end if

      do while(.true.)
200      continue
!        rc is the ratio of new to old values of the coefficient h*el(1).
!        when rc differs from 1 by more than 0.3, ipup is set to miter to force fjac to be called, if a jacobian is involved.
!        in any case, fjac is called at least every 20 steps. 0.3 is maximum relative change in h*el0 before fjac is called.
         if(dabs(rc - 1.0d0) .gt. 0.3 .or. nst .ge. nslp + 20) ipup = miter

         tn = tn + h

!        the predicted values is computed by effectively multiplying the yh array by the pascal triangle matrix.
         i1 = nqnyh + 1
         do jb = 1,nq
            i1 = i1 - nnyh
            do i = i1,nqnyh
               yh1(i) = yh1(i) + yh1(i+nnyh)
            end do
         end do

         exit1 = .false.
         do while(.not. exit1)
!           up to maxcor corrector iterations are taken.  
!           a convergence test is made on the r.m.s. norm of each correction, weighted by the error weight vector ewt.
!           the sum of the corrections is accumulated in the vector acor(i). the yh array is not altered in the corrector loop.
            m = 0
            do i = 1,n
               y(i) = yh(i,1)
            end do
            call rhs(neq, tn, y, savf)
            nfe = nfe + 1
            if(ipup .gt. 0)then
!              if indicated, the matrix p = i - h*el(1)*j is reevaluated and preprocessed before starting the corrector iteration.  
               call prjs(neq, y, yh, nnyh, ewt, acor, savf, wm, iwm)
!              ipup is set to 0 as an indicator that the matrix was evaluated.
               ipup = 0
               rc = 1.0d0
               nslp = nst
               crate = 0.7d0
               if(ierpj .ne. 0)then
                  icf = 2
                  ncf = ncf + 1
!                 the values of tn and the yh array are as of the beginning of the last step.
                  call set_back(told, yh1)
!                 kflag  = a completion code -3: fatal error in prjs.
                  kflag = -3
                  hold = h
                  jstart = 1
                  RETURN
!                 ------
               end if
            end if
            
            do i = 1,n
               acor(i) = 0.0d0
            end do
            
            exit2 = .false.
            exit3 = .false.
            do while(.not. exit2)
270            continue
               if(miter .eq. 0)then
!                 in case of functional iteration, update y directly from the result of the last function evaluation.
                  do i = 1,n
                     savf(i) = h*savf(i) - yh(i,2)
                     y(i) = savf(i) - acor(i)
                  end do
                  del = vnorm (n, y, ewt)
                  do i = 1,n
                     y(i) = yh(i,1) + el(1)*savf(i)
                     acor(i) = savf(i)
                  end do
               else
!                 in case of chord method, compute corrector error, and solve linear system with that as right-hand side and p as coefficient matrix.
                  do i = 1,n
                     y(i) = h*savf(i) - (yh(i,2) + acor(i))
                  end do
                  call slss(wm, iwm, y, savf)
!                 slss failed
                  if(iersl .lt. 0)then
                     icf = 2
                     ncf = ncf + 1
!                    the values of tn and the yh array are set as of the beginning of the last step.
                     call set_back(told, yh1)
!                    kflag  = a completion code -3: fatal error in slss.
                     kflag = -3
                     hold = h
                     jstart = 1
                     RETURN
!                    ------
                  end if

!                 slss was successful
                  if(iersl .eq. 0)then
                     del = vnorm(n, y, ewt)
                     do i = 1,n
                        acor(i) = acor(i) + y(i)
                        y(i) = yh(i,1) + el(1)*acor(i)
                     end do
                  end if
               end if

               if(iersl .eq. 0)then
!                 slss succesful or miter = 0. test for convergence.
!                 if m.gt.0, an estimate of the convergence rate constant is stored in crate, and this is used in the test.
                  if(m .gt. 0) crate = dmax1(0.2d0*crate,del/delp)
                  dcon = del*dmin1(1.0d0,1.5d0*crate)/(tesco(2,nq)*conit)
                  
                  if(dcon .le. 1.0d0)then
!                    the corrector has converged.  jcur is set to 0 to signal that the jacobian involved may need updating later.
!                    the local error test is made
                     jcur = 0
                     if(m .eq. 0)then
                        dsm = del/tesco(2,nq)
                     else
                        dsm = vnorm(n, acor, ewt)/tesco(2,nq)
                     end if

                     if(dsm .le. 1.0d0)then

!                       a successful step, update the yh array. consider changing h if ialth = 1.  otherwise decrease ialth by 1.
!                       if ialth is then 1 and nq .lt. maxord, then acor is saved for use in a possible order increase on the next step.
!                       if a change in h is considered, an increase or decrease in order by one is considered also.
!                       a change in h is made only if it is by afactor of at least 1.1. if not, ialth is set to 3 to prevent testing for that many steps.
                        kflag = 0
                        nst = nst + 1
                        hu = h
                        nqu = nq
                        do j = 1,nq+1
                           do i = 1,n
                              yh(i,j) = yh(i,j) + el(j)*acor(i)
                           end do
                        end do
                        ialth = ialth - 1
                        if(ialth .eq. 0)then
!                          the success of the step, factors rhdn, rhsm, and rhup are computed, 
!                          by which h could be multiplied at order nq - 1, order nq, or order nq + 1, respectively.
!                          the largest of these is determined and the new order chosen accordingly.
!                          if the order is to be increased, we compute one additional scaled derivative.
                           rhup = 0.0d0
                           if(nq .lt. maxord)then
                              do i = 1,n
                                 savf(i) = acor(i) - yh(i,maxord+1)
                              end do
                              dup = vnorm(n, savf, ewt)/tesco(3,nq)
                              exup = 1.0d0/dfloat(nq+2)
                              rhup = 1.0d0/(1.4d0*dup**exup + 0.0000014d0)
                           end if
                           
                           exsm = 1.0d0/dfloat(nq + 1)
                           rhsm = 1.0d0/(1.2d0*dsm**exsm + 0.0000012d0)
                           rhdn = 0.0d0
                           if(nq .ne. 1)then
                              ddn = vnorm (n, yh(1,nq+1), ewt)/tesco(1,nq)
                              exdn = 1.0d0/dfloat(nq)
                              rhdn = 1.0d0/(1.3d0*ddn**exdn + 0.0000013d0)
                           end if

                           rmax = max(rhdn, rhsm, rhup)  
                           if(rmax == rhdn)then
                              newq = nq - 1
                              rh = rhdn
                           else if(rmax == rhsm)then
                              newq = nq
                              rh = rhsm
                           else if(rmax == rhup)then
                              newq = nq + 1
                              rh = rhup
                              r = el(newq)/dfloat(newq)
                              do i = 1,n
                                 yh(i,newq+1) = acor(i)*r
                              end do
                           end if
!                          if h is being changed, the h ratio rh is checked against rmax, hmin, and hmxi, and the yh array rescaled.  
                           if(rh .ne. 1.0d0)then
                              rh = dmax1(rh,hmin/dabs(h))
!                             the h ratio rh is checked against rmax, hmin, and hmxi, and the yh array rescaled.  
                              call step_changed(rh, nyh, yh)
                              if(rh .lt. 1.1d0) ialth = 3
                              rmax = 10.0d0
                           end if

                        else if(ialth .eq. 1 .and. nq .lt. maxord)then
                           do i = 1,n
                              yh(i,maxord+1) = acor(i)
                           end do
                        end if
                        
                        !successful step
                        r = 1.0d0/tesco(2,nqu)
                        do i = 1,n
                           acor(i) = acor(i)*r
                        end do
                        hold = h
                        jstart = 1
                        RETURN
!                       ------ 
   
                     else !dsm .gt. 1

!                       the local error test failed. kflag keeps track of multiple failures.
                        kflag = kflag - 1

!                       restore tn and the yh array to their previous values, and prepare to try the step again.
                        call set_back(told, yh1)
                        if(dabs(h) .le. hmin*1.00001d0)then
!                          kflag  = a completion code -1: the requested error could not be achieved. abs(h) = hmin.
                           kflag = -1
!                          h is saved in hold to allow the caller to change h on the next step.
                           hold = h
                           jstart = 1
                           RETURN
!                          ------
                        end if
!                       if 3 or more failures have occured...
                        if(kflag .le. -3)then
!                          if 10 failures have occurred, exit with kflag = -1.
                           if(kflag .eq. -10)then
!                             kflag = a completion code -1: the requested error could not be achieved. 10 consecutive failures occurred.
!                             the values of tn and the yh array are as of the beginning of the last step, and h is the last step size attempted.
!                             h is saved in hold to allow the caller to change h on the next step.
                              kflag = -1
                              hold = h
                              jstart = 1
                              RETURN
!                             ------
                           end if
!                          it is assumed that the derivatives that have accumulated in the yh array have errors of the wrong order.  
!                          hence the first derivative is recomputed, and the order is set to 1. 
                           do i = 1,n
                              y(i) = yh(i,1)
                           end do
                           call rhs(neq, tn, y, savf)
                           nfe = nfe + 1
                           do i = 1,n
                              yh(i,2) = h*savf(i)
                           end do
!                          then h is reduced by a factor of 10, and the step is retried, until it succeeds or h reaches hmin.
                           rh = 0.1d0
                           rh = dmax1(hmin/dabs(h),rh)
                           h = h*rh
                           ipup = miter
                           if(nq .ne. 1)then
                              nq = 1
!                             el vector and related constants are reset because the order nq is changed.
                              call order_changed()
                           end if
                           ialth = 5
                           GO TO 200
                           
                        else ! less than 3 failures...

!                          failure of the step, factors rhdn, rhsm, and rhup are computed, 
!                          by which h could be multiplied at order nq - 1, order nq, or order nq + 1, respectively.
!                          in the case of failure, rhup = 0.0 to avoid an order increase. the largest of these is determined and the new order chosen accordingly.
!                          if the order is to be increased, we compute one additional scaled derivative.
                           rhup = 0.0d0
                           exsm = 1.0d0/dfloat(nq + 1)
                           rhsm = 1.0d0/(1.2d0*dsm**exsm + 0.0000012d0)
                           rhdn = 0.0d0 
                           if(nq .ne. 1)then
                              ddn = vnorm (n, yh(1,nq+1), ewt)/tesco(1,nq)
                              exdn = 1.0d0/dfloat(nq)
                              rhdn = 1.0d0/(1.3d0*ddn**exdn + 0.0000013d0)
                           end if

                           rmax = max(rhdn, rhsm)  
                           if(rmax == rhdn)then
                              newq = nq - 1
                              rh = rhdn
                           else if(rmax == rhsm)then
                              newq = nq
                              rh = rhsm
                           end if
                           if(rh .gt. 1.0d0) rh = 1.0d0
                           if(kflag .eq. -2) rh = dmin1(rh,0.2d0)

!                          if h is being changed, the h ratio rh is checked against rmax, hmin, and hmxi, and the yh array rescaled.
                           if(rh .lt. 1.0d0)then
                              rh = dmax1(rh,hmin/dabs(h))
                              call step_changed(rh, nyh, yh)
                           end if

!                          if there is a change of order, reset nq and the coefficients.
!                          h is reset according to rh and the yh array is rescaled. then redo the step otherwise.
                           if(newq .ne. nq)then
                              nq = newq
!                             el vector and related constants are reset when the order nq is changed.
                              call order_changed()
                           end if
                        end if ! number of failures
                     end if ! local error test (dsm)
                     
                  else ! dcon > 1
   
!                    the corrector iteration failed to converge.
                     m = m + 1
                     if(m .lt. maxcor .and. (m .eq. 1 .or. del .le. 2.0d0*delp))then
                        delp = del
                        call rhs(neq, tn, y, savf)
                        nfe = nfe + 1
                        exit3 = .true.
                     end if

                  end if ! corrector convergence test (dcon <= 1.0d0)

               end if ! iersl .eq. 0

               if(exit3)then
                  exit3 = .false.
               else
                  if(miter .ne. 0 .and. jcur .eq. 0)then
!                    if miter .ne. 0 and the jacobian is out of date, pjac is called for the next try.
                     icf = 1
                     ipup = miter
                     exit2 = .true.
                  else
                     exit2 = .true.
                     exit1 = .true.
                  end if
               end if
            end do !exit2
         end do !exit1

!        the yh array is retracted to its values before prediction, and h is reduced, if possible.
         icf = 2
         ncf = ncf + 1

         call set_back(told, yh1)

         if(dabs(h) .le. hmin*1.00001d0 .or. ncf .eq. 10)then
!           kflag  = a completion code -2: corrector convergence could not be achieved. either abs(h) = hmin or 10 consecutive failures occurred.
            kflag = -2
!           h is saved in hold to allow the caller to change h on the next step.
            hold = h
            jstart = 1
            RETURN
!           ------
         end if

         rh = 0.25d0
         ipup = miter
         rh = dmax1(rh,hmin/dabs(h))
!        if h is being changed, the h ratio rh is checked against rmax, hmin, and hmxi, and the yh array rescaled.
         call step_changed(rh, nyh, yh)
      end do
      end


      subroutine order_changed()
      integer i

      double precision conit, crate, el, elco, hold, rmax, tesco, 
     +   ccmax, el0, h, hmin, hmxi, hu, rc, tn, uround
      integer illin, init, lyh, lewt, lacor, lsavf, lwm, liwm, mxstep, mxhnil, nhnil, ntrep, nslast, nyh,
     +   ialth, ipup, lmax, meo, nqnyh, nslp, 
     +   icf, ierpj, iersl, jcur, jstart, kflag, l, meth, miter, 
     +   maxord, maxcor, msbp, n, nq, nst, nfe, nje, nqu
      common /ls0001/ conit, crate, el(13), elco(13,12), hold, rmax, tesco(3,12),
     +   ccmax, el0, h, hmin, hmxi, hu, rc, tn, uround, 
     +   illin, init, lyh, lewt, lacor, lsavf, lwm, liwm, mxstep, mxhnil, nhnil, ntrep, nslast, nyh, 
     +   ialth, ipup, lmax, meo, nqnyh, nslp, 
     +   icf, ierpj, iersl, jcur, jstart, kflag, l, meth, miter, 
     +   maxord, maxcor, msbp, mxncf, n, nq, nst, nfe, nje, nqu

!     el vector and related constants are reset because the order nq is changed.
      do i = 1,nq+1
         el(i) = elco(i,nq)
      end do
      nqnyh = nq*nnyh
      rc = rc*el(1)/el0
      el0 = el(1)
      conit = 0.5d0/dfloat(nq+2)
!     ialth is set to nq + 1 to prevent a change of h for that many steps.
      ialth = nq + 1
      return
      end


      subroutine step_changed(rh, nnyh, yh)

      integer i, j, nnyh
      double precision rh, r, yh(nnyh,1)

      double precision conit, crate, el, elco, hold, rmax, tesco, 
     +   ccmax, el0, h, hmin, hmxi, hu, rc, tn, uround
      integer illin, init, lyh, lewt, lacor, lsavf, lwm, liwm, mxstep, mxhnil, nhnil, ntrep, nslast, nyh,
     +   ialth, ipup, lmax, meo, nqnyh, nslp, 
     +   icf, ierpj, iersl, jcur, jstart, kflag, l, meth, miter, 
     +   maxord, maxcor, msbp, n, nq, nst, nfe, nje, nqu
      common /ls0001/ conit, crate, el(13), elco(13,12), hold, rmax, tesco(3,12),
     +   ccmax, el0, h, hmin, hmxi, hu, rc, tn, uround, 
     +   illin, init, lyh, lewt, lacor, lsavf, lwm, liwm, mxstep, mxhnil, nhnil, ntrep, nslast, nyh, 
     +   ialth, ipup, lmax, meo, nqnyh, nslp, 
     +   icf, ierpj, iersl, jcur, jstart, kflag, l, meth, miter, 
     +   maxord, maxcor, msbp, mxncf, n, nq, nst, nfe, nje, nqu

!     if h is being changed, the h ratio rh is checked against rmax, hmin, and hmxi, and the yh array rescaled.
      rh = dmin1(rh,rmax)
      rh = rh/dmax1(1.0d0,dabs(h)*hmxi*rh)
      r = 1.0d0
      do j = 2,nq+1
         r = r*rh
         do i = 1,n
            yh(i,j) = yh(i,j)*r
         end do
      end do
      h = h*rh
      rc = rc*rh
!     ialth is set to nq + 1 to prevent a change of h for that many steps.
      ialth = nq + 1
      return
      end




      subroutine set_back(told, yh1)

      integer i1, jb, i
      double precision told, yh1(1)

      double precision conit, crate, el, elco, hold, rmax, tesco, 
     +   ccmax, el0, h, hmin, hmxi, hu, rc, tn, uround
      integer illin, init, lyh, lewt, lacor, lsavf, lwm, liwm, mxstep, mxhnil, nhnil, ntrep, nslast, nyh,
     +   ialth, ipup, lmax, meo, nqnyh, nslp, 
     +   icf, ierpj, iersl, jcur, jstart, kflag, l, meth, miter, 
     +   maxord, maxcor, msbp, n, nq, nst, nfe, nje, nqu
      common /ls0001/ conit, crate, el(13), elco(13,12), hold, rmax, tesco(3,12),
     +   ccmax, el0, h, hmin, hmxi, hu, rc, tn, uround, 
     +   illin, init, lyh, lewt, lacor, lsavf, lwm, liwm, mxstep, mxhnil, nhnil, ntrep, nslast, nyh, 
     +   ialth, ipup, lmax, meo, nqnyh, nslp, 
     +   icf, ierpj, iersl, jcur, jstart, kflag, l, meth, miter, 
     +   maxord, maxcor, msbp, mxncf, n, nq, nst, nfe, nje, nqu

!     the maximum ratio by which h can be increased in a single step.  
      rmax = 2.0d0
      
!     the values of tn and the yh array are set back as of the beginning of the last step.
      tn = told
      
      i1 = nqnyh + 1
      do jb = 1,nq
         i1 = i1 - nyh
         do i = i1,nqnyh
            yh1(i) = yh1(i) - yh1(i+nyh)
         end do
      end do
      return
      end
