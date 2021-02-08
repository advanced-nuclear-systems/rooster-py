      subroutine lsodes(neq, y, t, tout, itol, rtol, atol, itask, istate, iopt, rwork, lrw, iwork, liw, jac, mf)
      external jac
      integer neq, itol, itask, iopt, lrw, iwork, liw, mf
      double precision y, t, tout, rtol, atol, rwork
      dimension y(1), rtol(1), atol(1), rwork(lrw), iwork(liw)

      integer illin, init, lyh, lewt, lacor, lsavf, lwm, liwm, mxstep, mxhnil, nhnil, nyh, iowns
      integer icf, ierpj, iersl, jcur, jstart, kflag, l, meth, miter, maxord, maxcor, n, nq, nst, nfe, nje, nqu
      integer iplost, iesp, iys, iba, ibian, ibjan, ibjgp, ipian, ipjan, ipjgp, ipigp, ipr, ipc, ipic, ipisp, iprsp, ipa, lenyh, lenyhm, lenwk, lreq, lrest, lwmin, moss, nslj, ngp, nlu, nnz, nsp, nzl, nzu
      integer i, i1, i2, iflag, imax, imul, imxer, ipflag, ipgo, irem, j, kgo, lenyht, leniw, lenrw, lf0, lia, lja, lrtem, lwtem, lyhd, lyhn, mf1, mord, mxhnl0, mxstp0, ncolm
      double precision rowns, el0, h, hmin, hmxi, hu, rc, tn, uround
      double precision con0, conmin, ccmxj, psmall, rbig, seth
      double precision atoli, ayi, h0, hmax, hmx, rh, rtoli, tdist, tnext, tol, tp, size, sum, w0, d1mach, vnorm
      dimension mord(2)
      common /ls0001/ rowns(209), el0, h, hmin, hmxi, hu, rc, tn, uround, illin, init, lyh, lewt, lacor, lsavf, lwm, liwm, mxstep, mxhnil, nhnil, nyh, iowns(6), icf, ierpj, iersl, jcur, jstart, kflag, l, meth, miter, maxord, maxcor, n, nq, nst, nfe, nje, nqu
      common /lss001/ con0, conmin, ccmxj, psmall, rbig, seth, iplost, iesp, iys, iba, ibian, ibjan, ibjgp, ipian, ipjan, ipjgp, ipigp, ipr, ipc, ipic, ipisp, iprsp, ipa, lenyh, lenyhm, lenwk, lreq, lrest, lwmin, moss, nslj, ngp, nlu, nnz, nsp, nzl, nzu

      data mord(1),mord(2)/12,5/, mxstp0/500/, mxhnl0/10/

      if(istate .eq. 1)then
         init = 0
         n = neq
         moss = mf/100
         mf1 = mf - 100*moss
         meth = mf1/10
         miter = mf1 - 10*meth
         if(miter .eq. 0 .or. miter .eq. 3) moss = 0
         if(iopt .eq. 1)then
            maxord = iwork(5)
            if(maxord .eq. 0) maxord = 100
            maxord = min0(maxord,mord(meth))
            mxstep = iwork(6)
            if(mxstep .eq. 0) mxstep = mxstp0
            mxhnil = iwork(7)
            if(mxhnil .eq. 0) mxhnil = mxhnl0
            h0 = rwork(5)
            hmax = rwork(6)
            hmxi = 0.0d0
            if(hmax .gt. 0.0d0) hmxi = 1.0d0/hmax
            hmin = rwork(7)
            seth = rwork(8)
         else
            maxord = mord(meth)
            mxstep = mxstp0
            mxhnil = mxhnl0
            h0 = 0.0d0
            hmxi = 0.0d0
            hmin = 0.0d0
            seth = 0.0d0
         end if

         rtoli = rtol(1)
         atoli = atol(1)
         do i = 1,n
           if(itol .ge. 3) rtoli = rtol(i)
           if(itol .eq. 2 .or. itol .eq. 4) atoli = atol(i)
         end do

         nyh = n
         lwmin = 0
         if(miter .eq. 1) lwmin = 4*n + 10*n/2
         if(miter .eq. 2) lwmin = 4*n + 11*n/2
         if(miter .eq. 3) lwmin = n + 2
         lenyh = (maxord+1)*nyh
         lrest = lenyh + 3*n
         lenrw = 20 + lwmin + lrest
         iwork(17) = lenrw
         leniw = 30
         if(moss .eq. 0 .and. miter .ne. 0 .and. miter .ne. 3) leniw = leniw + n + 1
         iwork(18) = leniw
         lia = 31
         if(moss .eq. 0 .and. miter .ne. 0 .and. miter .ne. 3) leniw = leniw + iwork(lia+n) - 1
         iwork(18) = leniw
         lja = lia + n + 1
         lia = min0(lia,liw)
         lja = min0(lja,liw)
         lwm = 21
         nq = 1
         ncolm = min0(nq+1,maxord+2)
         lenyhm = ncolm*nyh
         lenyht = lenyh
         if(miter .eq. 1 .or. miter .eq. 2) lenyht = lenyhm
         imul = 2
         if(moss .eq. 2) imul = 3
         lrtem = lenyht + imul*n
         lwtem = lwmin
         if(miter .eq. 1 .or. miter .eq. 2) lwtem = lrw - 20 - lrtem
         lenwk = lwtem
         lyhn = lwm + lwtem
         lsavf = lyhn + lenyht
         lewt = lsavf + n
         lacor = lewt + n
         lyh = lyhn
         iwork(22) = lyh
         tn = t
         nst = 0
         h = 1.0d0
         nnz = 0
         ngp = 0
         nzl = 0
         nzu = 0
         do i = 1,n
            rwork(i+lyh-1) = y(i)
         end do
         lf0 = lyh + nyh
         call rhs(neq, t, y, rwork(lf0))
         nfe = 1
!        ewset sets the error weight vector ewt before each step.
         call ewset(n, itol, rtol, atol, rwork(lyh), rwork(lewt))
         do i = 1,n
            rwork(i+lewt-1) = 1.0d0/rwork(i+lewt-1)
         end do
         if(miter .eq. 1 .or. miter .eq. 2)then
            lacor = min0(lacor,lrw)
!           interface between lsodes and prep, and also does adjusting of work space pointers and work arrays
            call iprep(neq, y, rwork, iwork(lia), iwork(lja), ipflag)
            lenrw = lwm - 1 + lenwk + lrest
            iwork(17) = lenrw
            if(ipflag .ne. -1) iwork(23) = ipian
            if(ipflag .ne. -1) iwork(24) = ipjan
            ipgo = -ipflag + 1
            iwork(22) = lyh
         end if

         uround = 1.0e-14
         jstart = 0
         if(miter .ne. 0) rwork(lwm) = dsqrt(uround)
         nslj = 0
         ccmxj = 0.2d0
         psmall = 1000.0d0*uround
         rbig = 0.01d0/psmall
         nhnil = 0
         nje = 0
         nlu = 0
         hu = 0.0d0
         nqu = 0
         maxcor = 3
         lf0 = lyh + nyh
         if(h0 .eq. 0.0d0)then
            tdist = dabs(tout - t)
            w0 = dmax1(dabs(t),dabs(tout))
            tol = rtol(1)
            if(itol .gt. 2)then
               do i = 1,n
                  tol = dmax1(tol,rtol(i))
               end do
            end if
            if(tol .le. 0.0d0)then
               atoli = atol(1)
               do i = 1,n
                  if(itol .eq. 2 .or. itol .eq. 4) atoli = atol(i)
                  ayi = dabs(y(i))
                  if(ayi .ne. 0.0d0) tol = dmax1(tol,atoli/ayi)
               end do
            end if
            tol = dmax1(tol,100.0d0*uround)
            tol = dmin1(tol,0.001d0)
            sum = vnorm (n, rwork(lf0), rwork(lewt))
            sum = 1.0d0/(tol*w0*w0) + tol*sum**2
            h0 = 1.0d0/dsqrt(sum)
            h0 = dmin1(h0,tdist)
            h0 = dsign(h0,tout-t)
         end if
         rh = dabs(h0)*hmxi
         if(rh .gt. 1.0d0) h0 = h0/rh
         h = h0
         do i = 1,n
            rwork(i+lf0-1) = h0*rwork(i+lf0-1)
         end do
      end if
         
!     ewset sets the error weight vector ewt before each step.
      call ewset (n, itol, rtol, atol, rwork(lyh), rwork(lewt))
      do i = 1,n
         rwork(i+lewt-1) = 1.0d0/rwork(i+lewt-1)
      end do

!     core integrator, which does one step of the integration and the associated error control
      call stode(neq, y, rwork(lyh), nyh, rwork(lyh), rwork(lewt), rwork(lsavf), rwork(lacor), rwork(lwm), rwork(lwm))
      kgo = 1 - kflag

      init = 1
 
      t = tout

      first_step = .false.
      illin = 0
      rwork(11) = hu
      rwork(12) = h
      rwork(13) = tn
      iwork(11) = nst
      iwork(12) = nfe
      iwork(13) = nje
      iwork(14) = nqu
      iwork(15) = nq
      iwork(19) = nnz
      iwork(20) = ngp
      iwork(21) = nlu
      iwork(25) = nzl
      iwork(26) = nzu
      return
      end