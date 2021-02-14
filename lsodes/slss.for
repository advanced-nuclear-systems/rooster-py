
! this routine manages the solution of the linear system arising from a chord iteration.
! it is called if miter .ne. 0. if miter is 1 or 2, it calls cdrv to accomplish this.
! if miter = 3 it updates the coefficient h*el0 in the diagonal matrix, and then computes the solution.
! communication with slss uses the following variables..
! wk    = real work space containing the inverse diagonal matrix if miter = 3 and the lu decomposition of the matrix otherwise.
!         storage of matrix elements starts at wk(3). wk also contains the following matrix-related data..
!         wk(1) = sqrt(uround) (not used here),
!         wk(2) = hl0, the previous value of h*el0, used if miter = 3.
! iwk   = integer work space for matrix-related data, assumed to be equivalenced to wk.  
!         in addition, wk(iprsp) and iwk(ipisp) are assumed to have identical locations.
! x     = the right-hand side vector on input, and the solution vector on output, of length n.
! tem   = vector of work space of length n, not used in this version.
! iersl = output flag (in common).
!         iersl = 0  if no trouble occurred.
!         iersl = -1 if cdrv returned an error flag (miter = 1 or 2). this should never occur and is considered fatal.
!         iersl = 1  if a singular matrix arose with miter = 3.
! this routine also uses other variables in common.

      subroutine slss(wk, iwk, x, tem)

      integer iwk
      integer i
      double precision wk, x, tem
      double precision di, hl0, phl0, r
      dimension wk(1), iwk(1), x(1), tem(1)

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

      iersl = 0
      if(miter .eq. 3)then
         phl0 = wk(2)
         hl0 = h*el0
         wk(2) = hl0
         if(hl0 .ne. phl0)then
            r = hl0/phl0
            do i = 1,n
              di = 1.0d0 - r*(1.0d0 - 1.0d0/wk(i+2))
              if(dabs(di) .eq. 0.0d0)then
                 iersl = 1
                 return
              end if
              wk(i+2) = 1.0d0/di
            end do
         end if
         do i = 1,n
            x(i) = wk(i+2)*x(i)
         end do
      else
         call cdrv(n,iwk(ipr),iwk(ipc),iwk(ipic),iwk(ipian),iwk(ipjan),wk(ipa),x,x,nsp,iwk(ipisp),wk(iprsp),iesp,4,iersl)
         if (iersl .ne. 0) iersl = -1
      end if
      return
      end
