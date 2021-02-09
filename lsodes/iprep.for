
! this routine serves as an interface between the driver and subroutine prep. it is called only if miter is 1 or 2.
! tasks performed here are..
!  * call prep,
!  * reset the required wm segment length lenwk,
!  * move yh back to its final location (following wm in rwork),
!  * reset pointers for yh, savf, ewt, and acor, and
!  * move ewt to its new position if istate = 1.
! ipflag is an output error indication flag.  ipflag = 0 if there was no trouble, and ipflag is the value of the prep error flag ipper
! if there was trouble in subroutine prep.

      subroutine iprep (neq, y, rwork, ia, ja, ipflag)

      integer neq, ia, ja, ipflag
      integer i, imax, lewtn, lyhd, lyhn
      double precision y, rwork
      dimension neq(1), y(1), rwork(1), ia(1), ja(1)

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

      double precision con0, conmin, ccmxj, psmall, rbig, seth
      integer iplost, iesp, istatc, iys, iba, ibian, ibjan, ibjgp, 
     +   ipian, ipjan, ipjgp, ipigp, ipr, ipc, ipic, ipisp, iprsp, ipa, 
     +   lenyh, lenyhm, lenwk, lreq, lrat, lrest, lwmin, moss, msbj,
     +   nslj, ngp, nlu, nnz, nsp, nzl, nzu
      common /lss001/ con0, conmin, ccmxj, psmall, rbig, seth,
     +   iplost, iesp, istatc, iys, iba, ibian, ibjan, ibjgp,
     +   ipian, ipjan, ipjgp, ipigp, ipr, ipc, ipic, ipisp, iprsp, ipa,
     +   lenyh, lenyhm, lenwk, lreq, lrat, lrest, lwmin, moss, msbj,
     +   nslj, ngp, nlu, nnz, nsp, nzl, nzu

      ipflag = 0
!     call prep to do matrix preprocessing operations.
      call prep(neq, y, rwork(lyh), rwork(lsavf), rwork(lewt), rwork(lacor), ia, ja, rwork(lwm), rwork(lwm), ipflag)
      lenwk = max0(lreq,lwmin)
      if(ipflag .lt. 0) RETURN

!     if prep was successful, move yh to end of required space for wm.
      lyhn = lwm + lenwk
      if(lyhn .gt. lyh) RETURN

      lyhd = lyh - lyhn
      if(lyhd .ne. 0)then
         imax = lyhn - 1 + lenyhm
         do i = lyhn,imax
            rwork(i) = rwork(i+lyhd)
         end do
         lyh = lyhn
      end if

!     reset pointers for savf, ewt, and acor.
      lsavf = lyh + lenyh
      lewtn = lsavf + n
      lacor = lewtn + n
!     if istate = 1, move ewt (left) to its new position.
      if(lewtn .gt. lewt) RETURN
      do i = 1,n
         rwork(i+lewtn-1) = rwork(i+lewt-1)
      end do
      lewt = lewtn
      RETURN
      end