!
! prjs is called to compute and process the matrix p = i - h*el(1)*j , where j is an approximation to the jacobian.
! j is computed by columns, either by the user-supplied routine jac if miter = 1, or by finite differencing if miter = 2.
! if miter = 3, a diagonal approximation to j is used.
! if miter = 1 or 2, and if the existing value of the jacobian (as contained in p) is considered acceptable, then a new value of p is reconstructed from the old value.
! in any case, when miter is 1 or 2, the p matrix is subjected to lu decomposition in cdrv. p and its lu decomposition are stored (separately) in wk.
!
! in addition to variables described previously, communication with prjs uses the following..
! y     = array containing predicted values on entry.
! ftem  = work array of length n (acor in stode).
! savf  = array containing f evaluated at predicted y.
! wk    = real work space for matrices.  on output it contains the inverse diagonal matrix if miter = 3, and p and its sparse lu decomposition if miter is 1 or 2.
!         storage of matrix elements starts at wk(3). wk also contains the following matrix-related data..
!         wk(1) = sqrt(uround), used in numerical jacobian increments.
!         wk(2) = h*el0, saved for later use if miter = 3.
! iwk   = integer work space for matrix-related data, assumed to be equivalenced to wk.  in addition, wk(iprsp) and iwk(ipisp) are assumed to have identical locations.
! el0   = el(1) (input).
! ierpj = output error flag (in common).
!       = 0 if no error.
!       = 1  if zero pivot found in cdrv.
!       = 2  if a singular matrix arose with miter = 3.
!       = -1 if insufficient storage for cdrv (should not occur here).
!       = -2 if other error found in cdrv (should not occur here).
! jcur  = output flag = 1 to indicate that the jacobian matrix (or approximation) is now current.
!
      subroutine prjs(neq, y, yh, nnyh, ewt, ftem, savf, wk, iwk)

      integer neq, nnyh, iwk
      integer i, imul, j, jj, jok, jmax, jmin, k, kmax, kmin, ng
      double precision y, yh, ewt, ftem, savf, wk
      double precision con, di, fac, hl0, pij, r, r0, rcon, rcont, srur, vnorm
      dimension y(1), yh(nnyh,1), ewt(1), ftem(1), savf(1), wk(1), iwk(1)                                                        

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

      hl0 = h*el0
      con = -hl0

      if(miter .eq. 3)then
!        if miter = 3, construct a diagonal approximation to j and p.
         jcur = 1
         nje = nje + 1
         wk(2) = hl0
         ierpj = 0
         r = el0*0.1d0
         do i = 1,n
            y(i) = y(i) + r*(h*savf(i) - yh(i,2))
         end do
         call rhs(neq, tn, y, wk(3))
         nfe = nfe + 1
         do i = 1,n
            r0 = h*savf(i) - yh(i,2)
            di = 0.1d0*r0 - h*(wk(i+2) - savf(i))
            wk(i+2) = 1.0d0
            if(dabs(r0) .ge. uround/ewt(i))then
               if(dabs(di) .eq. 0.0d0)then
                  ierpj = 2
                  return
               end if
               wk(i+2) = 0.1d0*r0/di
            end if
         end do
         return
      end if

!     see whether jacobian should be reevaluated (jok = 0) or not (jok = 1).
      jok = 1
      if(nst .eq. 0 .or. nst .ge. nslj+50) jok = 0
      if(icf .eq. 1 .and. dabs(rc - 1.0d0) .lt. ccmxj) jok = 0
      if(icf .eq. 2) jok = 0
      if(jok .eq. 1 .and. dabs(con)/conmin .gt. rbig .and. iplost .eq. 1) jok = 0

      if(jok .eq. 0)then
!        miter = 1 or 2, and the jacobian is to be reevaluated.
         jcur = 1
         nje = nje + 1
         nslj = nst
         iplost = 0
         conmin = dabs(con)
         if(miter .eq. 1)then
!           if miter = 1, call jac, multiply by scalar, and add identity.
            kmin = iwk(ipian)
            do j = 1, n
               kmax = iwk(ipian+j) - 1
               do i = 1,n
                  ftem(i) = 0.0d0
               end do
               call fjac(neq, tn, y, j, iwk(ipian), iwk(ipjan), ftem)
               do k = kmin, kmax
                  i = iwk(ibjan+k)
                  wk(iba+k) = ftem(i)*con
                  if (i .eq. j) wk(iba+k) = wk(iba+k) + 1.0d0
               end do
               kmin = kmax + 1
            end do
         else ! miter == 2
!           if miter = 2, make ngp calls to f to approximate j and p.
            fac = vnorm(n, savf, ewt)
            r0 = 1000.0d0 * dabs(h) * uround * dfloat(n) * fac
            if (r0 .eq. 0.0d0) r0 = 1.0d0
            srur = wk(1)
            jmin = iwk(ipigp)
            do ng = 1,ngp
               jmax = iwk(ipigp+ng) - 1
               do j = jmin,jmax
                  jj = iwk(ibjgp+j)
                  r = dmax1(srur*dabs(y(jj)),r0/ewt(jj))
                  y(jj) = y(jj) + r
               end do
               call rhs(neq, tn, y, ftem)
               do j = jmin,jmax
                  jj = iwk(ibjgp+j)
                  y(jj) = yh(jj,1)
                  r = dmax1(srur*dabs(y(jj)),r0/ewt(jj))
                  fac = -hl0/r
                  kmin =iwk(ibian+jj)
                  kmax =iwk(ibian+jj+1) - 1
                  do k = kmin,kmax
                    i = iwk(ibjan+k)
                    wk(iba+k) = (ftem(i) - savf(i))*fac
                    if (i .eq. jj) wk(iba+k) = wk(iba+k) + 1.0d0
                  end do
               end do
               jmin = jmax + 1
            end do
            nfe = nfe + ngp
         end if

      else ! jok .eq. 1

!        if jok = 1, reconstruct new p from old p.
         jcur = 0
         kmin = iwk(ipian)
         do j = 1,n
           kmax = iwk(ipian+j) - 1
           do k = kmin,kmax
             i = iwk(ibjan+k)
             pij = wk(iba+k)
             if(i .eq. j)then
                pij = pij - 1.0d0
                if(dabs(pij) .lt. psmall)then
                   iplost = 1
                   conmin = dmin1(dabs(con0),conmin)
                end if
             end if
             rcon = con/con0
             pij = pij*rcon
             if (i .eq. j) pij = pij + 1.0d0
             wk(iba+k) = pij
           end do
           kmin = kmax + 1
         end do
      end if

!     do numerical factorization of p matrix.
      nlu = nlu + 1
      con0 = con
      ierpj = 0
      do i = 1,n
         ftem(i) = 0.0d0
      end do
      call cdrv(n,iwk(ipr),iwk(ipc),iwk(ipic),iwk(ipian),iwk(ipjan), wk(ipa),ftem,ftem,nsp,iwk(ipisp),wk(iprsp),iesp,2,iys)
      if(iys .ne. 0) then
         imul = (iys - 1)/n
         ierpj = -2
         if(imul .eq. 8) ierpj = 1
         if(imul .eq. 10) ierpj = -1
      end if
      return
      end
