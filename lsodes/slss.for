      subroutine slss (wk, iwk, x, tem)
      integer iwk
      integer iownd, iowns, icf, ierpj, iersl, jcur, jstart, kflag, l, meth, miter, maxord, maxcor, msbp, n, nq, nst, nfe, nje, nqu
      integer iplost, iesp, iys, iba, ibian, ibjan, ibjgp, ipian, ipjan, ipjgp, ipigp, ipr, ipc, ipic, ipisp, iprsp, ipa, lenyh, lenyhm, lenwk, lreq, lrest, lwmin, moss, nslj, ngp, nlu, nnz, nsp, nzl, nzu
      integer i
      double precision wk, x, tem
      double precision rowns, el0, h, hmin, hmxi, hu, rc, tn, uround
      double precision rlss
      double precision di, hl0, phl0, r
      dimension wk(n), iwk(1), x(1), tem(1)                                   
      common /ls0001/ rowns(209), el0, h, hmin, hmxi, hu, rc, tn, uround, iownd(14), iowns(6), icf, ierpj, iersl, jcur, jstart, kflag, l, meth, miter, maxord, maxcor, msbp, n, nq, nst, nfe, nje, nqu
      common /lss001/ rlss(6), iplost, iesp, iys, iba, ibian, ibjan, ibjgp, ipian, ipjan, ipjgp, ipigp, ipr, ipc, ipic, ipisp, iprsp, ipa, lenyh, lenyhm, lenwk, lreq, lrest, lwmin, moss, nslj, ngp, nlu, nnz, nsp, nzl, nzu

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
