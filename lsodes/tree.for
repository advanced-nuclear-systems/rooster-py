lsodes:

!     block a.
!     this code block is executed on every call. it tests istate and itask for legality and branches appropriately.
!     if istate .gt. 1 but the flag init shows that initialization has not yet been done, an error return occurs.
!     if istate = 1 and tout = t, jump to block g and return immediately.

      if(istate .eq. 1 .or. istate .eq. 3)then
!        block b.
!        the next code block is executed for the initial call (istate = 1), or for a continuation call with parameter changes (istate = 3).
!        it contains checking of all inputs and various initializations. 
!        if istate = 1, the final setting of work space pointers, the matrix preprocessing, and other initializations are done in block c.
!        first check legality of the non-optional inputs neq, itol, iopt, mf, ml, and mu.

         if(istate .eq. 1)then
!           block c.
!           the next block is for the initial call only (istate = 1). it contains all remaining initializations, the initial call to f,
!           the sparse matrix preprocessing (miter = 1 or 2), and the calculation of the initial step size. the error weights in ewt are inverted after being loaded.

!           initial call to f.
            call rhs(neq, t, y, rwork(lf0))

!           load and invert the ewt array.  (h is temporarily set to 1.0.)
            call ewset(n, itol, rtol, atol, rwork(lyh), rwork(lewt))

            if(miter .eq. 1 .or. miter .eq. 2)then
!              iprep and prep do sparse matrix preprocessing if miter = 1 or 2.
               call iprep(neq, y, rwork, iwork(lia), iwork(lja), ipflag)
            end if

         else ! istate .eq. 3

            if(miter .eq. 1 .or. miter .eq. 2)then
               if(moss .eq. 2)then
!                 temporarily load ewt if miter = 1 or 2 and moss = 2.
                  call ewset (n, itol, rtol, atol, rwork(lyh), rwork(lewt))
                  ...
                  end do
               end if
!              iprep and prep do sparse matrix preprocessing if miter = 1 or 2.
               ...
               call iprep(neq, y, rwork, iwork(lia), iwork(lja), ipflag)
                    └─ ! do matrix preprocessing operations.
                       call prep(neq, y, rwork(lyh), rwork(lsavf), rwork(lewt), rwork(lacor), ia, ja, rwork(lwm), rwork(lwm), ipflag)
                            └─ 

               ...
            end if
            ...
         end if
      end if

      if(istate .eq. 2 .or. istate .eq. 3)then
!        block d.
!        the next code block is for continuation calls only (istate = 2 or 3) and is to check stop conditions before taking a step.
         nslast = nst
         if(itask .eq. 1)then
            ...
               ! computes interpolated values of the vector yh and stores it in y
               call intdy (tout, 0, rwork(lyh), nyh, y, iflag)
            ...
            end if
         else if(itask .eq. 3)then
            ...
         else if(itask .eq. 4)then
            ...
               ! computes interpolated values of the vector yh and stores it in y
               call intdy (tout, 0, rwork(lyh), nyh, y, iflag)
            ...
         else if(itask .eq. 5)then
            ...
         end if

         if(itask .eq. 4 .or. itask .eq. 5)then
            ...
         end if

!     block e.
!     the next block is normally executed for all calls and contains the call to the one-step core integrator stode.
!     this is a looping point for the integration steps.
!     first check for too many steps being taken, update ewt (if not at start of problem), check for too much accuracy being requested, and
!     check for h below the roundoff level in t.
      do while(.true.)
         ...
         call ewset (n, itol, rtol, atol, rwork(lyh), rwork(lewt))
         ...
         call stode (neq, y, rwork(lyh), nyh, rwork(lyh), rwork(lewt), rwork(lsavf), rwork(lacor), rwork(lwm), rwork(lwm))
         
!        block f.
!        the following block handles the case of a successful return from the core integrator (kflag = 0).  test for stop conditions.
         ...
         if(itask .eq. 1) then
         ...
         else if(itask .eq. 2) then
         ...
         else if(itask .eq. 3) then
         ...
         else if(itask .eq. 4) then
!           see if tout or tcrit was reached. adjust h if necessary.
         ...
               call intdy (tout, 0, rwork(lyh), nyh, y, iflag)
         ...
         else if(itask .eq. 5) then
      end do
