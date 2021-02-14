!
! this routine performs preprocessing related to the sparse linear systems that must be solved if miter = 1 or 2.
! the operations that are performed here are..
!  * compute sparseness structure of jacobian according to moss,
!  * compute grouping of column indices (miter = 2),
!  * compute a new ordering of rows and columns of the matrix,
!  * reorder ja corresponding to the new ordering,
!  * perform a symbolic lu factorization of the matrix, and
!  * set pointers for segments of the iwk/wk array.
! in addition to variables described previously, prep uses the following for communication..
! yh     = the history array.  only the first column, containing the current y vector, is used. used only if moss .ne. 0.
! savf   = a work array of length neq, used only if moss .ne. 0.
! ewt    = array of length neq containing (inverted) error weights. used only if moss = 2 or if istate = moss = 1.
! ftem   = a work array of length neq, identical to acor in the driver, used only if moss = 2.
! wk     = a real work array of length lenwk, identical to wm in the driver.
! iwk    = integer work array, assumed to occupy the same space as wk.
! lenwk  = the length of the work arrays wk and iwk.
! istatc = a copy of the driver input argument istate (= 1 on the first call, = 3 on a continuation call).
! iys    = flag value from odrv or cdrv.
! ipper  = output error flag with the following values and meanings..
!          0  no error.
!         -1  insufficient storage for internal structure pointers.
!         -2  insufficient storage for jgroup.
!         -3  insufficient storage for odrv.
!         -4  other error flag from odrv (should never occur).
!         -5  insufficient storage for cdrv.
!         -6  other error flag from cdrv.
! moss   = the method to be used to obtain the sparsity structure of the jacobian matrix if miter = 1 or 2..
!          0  means the user has supplied ia and ja (see descriptions under iwork in the lsodes).
!          1  means the user has supplied jac and the structure will be obtained from neq initial calls to jac.
!          2  means the structure will be obtained from neq+1 initial calls to f.
! miter  = the corrector iteration method..
!          0  means functional iteration (no jacobian matrix is involved).
!          1  means chord iteration with a user-supplied sparse jacobian, given by subroutine jac.
!          2  means chord iteration with an internally generated (difference quotient) sparse jacobian
!             (using ngp extra calls to f per df/dy value, where ngp is an optional output described below.)
!          3  means chord iteration with an internally generated diagonal jacobian approximation. (using 1 extra call to f per df/dy evaluation).
!
      subroutine prep (neq, y, yh, savf, ewt, ftem, ia, ja, wk, iwk, ipper)

      integer neq, ia, ja, iwk, ipper
      integer i, ibr, ier, ipil, ipiu, iptt1, iptt2, j, jfound, k, knew, kmax, kmin, ldif, lenigp, liwk, maxg, np1, nzsut
      double precision y, yh, savf, ewt, ftem, wk
      double precision dq, dyj, erwt, fac, yj
      dimension neq(1), y(1), yh(1), savf(1), ewt(1), ftem(1), ia(1), ja(1), wk(1), iwk(1)

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

      ibian = 4
      ipian = ibian + 1
      np1 = n + 1
      ipjan = ipian + np1
      ibjan = ipjan - 1
      liwk = lenwk*2
      if(ipjan+n-1 .gt. liwk)then
         ipper = -1
         lreq = 2 + (2*n + 1)/2
         lreq = max0(lenwk+1,lreq)
         RETURN
      end if

      if(moss .eq. 0)then
!        moss = 0. process user-s ia,ja. add diagonal entries if necessary.
         knew = ipjan
         kmin = ia(1)
         iwk(ipian) = 1
         do j = 1,n
            jfound = 0
            kmax = ia(j+1) - 1
            if(kmin .le. kmax)then
               do k = kmin,kmax
                  i = ja(k)
                  if(i .eq. j) jfound = 1
                  if(knew .gt. liwk)then
                     ipper = -1
                     lreq = 2 + (2*n + 1)/2
                     lreq = max0(lenwk+1,lreq)
                     RETURN
                  end if
                  iwk(knew) = i
                  knew = knew + 1
               end do
            end if
            if(jfound .eq. 0)then
               if(knew .gt. liwk)then
                  ipper = -1
                  lreq = 2 + (2*n + 1)/2
                  lreq = max0(lenwk+1,lreq)
                  RETURN
               end if
               iwk(knew) = j
               knew = knew + 1
            end if
            iwk(ipian+j) = knew + 1 - ipjan
            kmin = kmax + 1
         end do

      else ! moth .ne. 0

         if(istatc .eq. 3)then
!           istate = 3 and moss .ne. 0. load y from yh(*,1).
            do i = 1,n
               y(i) = yh(i)
            end do
         else
!           istate = 1 and moss .ne. 0. perturb y for structure determination.
            do i = 1,n
              erwt = 1.0d0/ewt(i)
              fac = 1.0d0 + 1.0d0/(dfloat(i)+1.0d0)
              y(i) = y(i) + fac*dsign(erwt,y(i))
            end do
         end if
      end if
      
      if(moss .eq. 1)then
         continue
!        a dummy call to rhs allows user to create temporaries for use in jac.
         call rhs(neq, tn, y, savf)
         k = ipjan
         iwk(ipian) = 1
         do j = 1,n
            if(k .gt. liwk)then
               ipper = -1
               lreq = 2 + (2*n + 1)/2
               lreq = max0(lenwk+1,lreq)
               RETURN
            end if
            iwk(k) = j
            k = k + 1
            do i = 1,n
               savf(i) = 0.0d0
            end do
            call fjac(neq, tn, y, j, iwk(ipian), iwk(ipjan), savf)
            do i = 1,n
               if(dabs(savf(i)) .gt. seth .and. i .ne. j)then
                  if(k .gt. liwk)then
                     ipper = -1
                     lreq = 2 + (2*n + 1)/2
                     lreq = max0(lenwk+1,lreq)
                     RETURN
                  end if
                  iwk(k) = i
                  k = k + 1
               end if
            end do
            iwk(ipian+j) = k + 1 - ipjan
         end do
      end if

      if(moss .eq. 2)then
!         moss = 2. compute structure from results of n + 1 calls to f.
         k = ipjan
         iwk(ipian) = 1
         call rhs(neq, tn, y, savf)
         do j = 1,n
            if(k .gt. liwk)then
               ipper = -1
               lreq = 2 + (2*n + 1)/2
               lreq = max0(lenwk+1,lreq)
               RETURN
            end if
            iwk(k) = j
            k = k + 1
            yj = y(j)
            erwt = 1.0d0/ewt(j)
            dyj = dsign(erwt,yj)
            y(j) = yj + dyj
            call rhs(neq, tn, y, ftem)
            y(j) = yj
            do i = 1,n
               dq = (ftem(i) - savf(i))/dyj
               if(dabs(dq) .gt. seth .and. i .ne. j)then
                  if(k .gt. liwk)then
                     ipper = -1
                     lreq = 2 + (2*n + 1)/2
                     lreq = max0(lenwk+1,lreq)
                     RETURN
                  end if
                  iwk(k) = i
                  k = k + 1
               end if
            end do
            iwk(ipian+j) = k + 1 - ipjan
         end do
      end if

      if(moss .ne. 0 .and. istatc .eq. 1)then ! CHECK that istatc .eq. 1 is correct here
!        if istate = 1 and moss .ne. 0, restore y from yh.
         do i = 1,n
            y(i) = yh(i)
         end do
      end if
      nnz = iwk(ipian+n) - 1
      lenigp = 0
      ipigp = ipjan + nnz
      if(miter .eq. 2)then
!        compute grouping of column indices (miter = 2).
         maxg = np1
         ipjgp = ipjan + nnz
         ibjgp = ipjgp - 1
         ipigp = ipjgp + n
         iptt1 = ipigp + np1
         iptt2 = iptt1 + n
         lreq = iptt2 + n - 1
         if(lreq .gt. liwk)then
            ipper = -2
            lreq = (lreq - 1)/2 + 1
            RETURN
         end if
         call jgroup (n, iwk(ipian), iwk(ipjan), maxg, ngp, iwk(ipigp), iwk(ipjgp), iwk(iptt1), iwk(iptt2), ier)
         if(ier .ne. 0)then
            ipper = -2
            lreq = (lreq - 1)/2 + 1
            RETURN
         end if
         lenigp = ngp + 1
      end if
!     compute new ordering of rows/columns of jacobian.
      ipr = ipigp + lenigp
      ipc = ipr
      ipic = ipc + n
      ipisp = ipic + n
      iprsp = (ipisp - 2)/2 + 2
      iesp = lenwk + 1 - iprsp
      if(iesp .lt. 0)then
         ipper = -3
         call cntnzu(n, iwk(ipian), iwk(ipjan), nzsut)
         lreq = lenwk - iesp + (3*n + 4*nzsut - 1)/2 + 1
         RETURN
      end if
      ibr = ipr - 1
      do i = 1,n
         iwk(ibr+i) = i
      end do
      nsp = liwk + 1 - ipisp
      call odrv(n, iwk(ipian), iwk(ipjan), wk, iwk(ipr), iwk(ipic), nsp, iwk(ipisp), 1, iys)
      if(iys .eq. 11*n+1)then
         ipper = -4
         RETURN
      end if
      if(iys .ne. 0)then
         ipper = -3
         call cntnzu(n, iwk(ipian), iwk(ipjan), nzsut)
         lreq = lenwk - iesp + (3*n + 4*nzsut - 1)/2 + 1
         RETURN
      end if

!     reorder jan and do symbolic lu factorization of matrix.
      ipa = lenwk + 1 - nnz
      nsp = ipa - iprsp
      lreq = max0(12*n/2, 6*n/2+2*n+nnz) + 3
      lreq = lreq + iprsp - 1 + nnz
      if(lreq .gt. lenwk)then
         ipper = -5
         RETURN
      end if
      iba = ipa - 1
      do i = 1,nnz
         wk(iba+i) = 0.0d0
      end do
      ipisp = 2*(iprsp - 1) + 1
      call cdrv(n,iwk(ipr),iwk(ipc),iwk(ipic),iwk(ipian),iwk(ipjan),wk(ipa),wk(ipa),wk(ipa),nsp,iwk(ipisp),wk(iprsp),iesp,5,iys)
      lreq = lenwk - iesp
      if(iys .eq. 10*n+1)then
         ipper = -5
         RETURN
      end if
      if(iys .ne. 0)then
         ipper = -6
         lreq = lenwk
         RETURN
      end if
      ipil = ipisp
      ipiu = ipil + 2*n + 1
      nzu = iwk(ipil+n) - iwk(ipil)
      nzl = iwk(ipiu+n) - iwk(ipiu)
      if(nnz .eq. n) lreq = lreq + 1
      nsp = nsp + lreq - lenwk
      ipa = lreq + 1 - nnz
      iba = ipa - 1
      ipper = 0
      RETURN
      end
