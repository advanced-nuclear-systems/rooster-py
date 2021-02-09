      subroutine slss(wk, iwk, x, tem)

      integer iwk
      integer i
      double precision wk, x, tem
      double precision di, hl0, phl0, r
      dimension wk(n), iwk(1), x(1), tem(1)                                   

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
