      subroutine iprep (neq, y, rwork, ia, ja, ipflag)

      integer neq, ia, ja, ipflag
      integer illin, init, lyh, lewt, lacor, lsavf, lwm, liwm,
     1   mxstep, mxhnil, nhnil, ntrep, nslast, nyh, iowns
      integer icf, ierpj, iersl, jcur, jstart, kflag, l, meth, miter,
     1   maxord, maxcor, n, nq, nst, nfe, nje, nqu
      integer iplost, iesp, iys, iba, ibian, ibjan, ibjgp,
     1   ipian, ipjan, ipjgp, ipigp, ipr, ipc, ipic, ipisp, iprsp, ipa,
     2   lenyh, lenyhm, lenwk, lreq, lrest, lwmin, moss, nslj, ngp, nlu, nnz, nsp, nzl, nzu
      integer i, imax, lewtn, lyhd, lyhn
      double precision y, rwork
      double precision rowns,
     1   el0, h, hmin, hmxi, hu, rc, tn, uround
      double precision rlss
      dimension neq(1), y(1), rwork(1), ia(1), ja(1)
      common /ls0001/ rowns(209), el0, h, hmin, hmxi, hu, rc, tn, uround, illin, init, lyh, lewt, lacor, lsavf, lwm, liwm, mxstep, mxhnil, nhnil, ntrep, nslast, nyh, iowns(6), icf, ierpj, iersl, jcur, jstart, kflag, l, meth, miter, maxord, maxcor, n, nq, nst, nfe, nje, nqu
      common /lss001/ rlss(6),
     1   iplost, iesp, iys, iba, ibian, ibjan, ibjgp,
     2   ipian, ipjan, ipjgp, ipigp, ipr, ipc, ipic, ipisp, iprsp, ipa,
     3   lenyh, lenyhm, lenwk, lreq, lrest, lwmin, moss, nslj, ngp, nlu, nnz, nsp, nzl, nzu
      ipflag = 0
      call prep (neq, y, rwork(lyh), rwork(lsavf), rwork(lewt), rwork(lacor), ia, ja, rwork(lwm), rwork(lwm), ipflag)
      lenwk = max0(lreq,lwmin)
      if (ipflag .lt. 0) return
      lyhn = lwm + lenwk
      if (lyhn .gt. lyh) return
      lyhd = lyh - lyhn
      if (lyhd .eq. 0) go to 20
      imax = lyhn - 1 + lenyhm
      do 10 i = lyhn,imax
 10     rwork(i) = rwork(i+lyhd)
      lyh = lyhn
 20   lsavf = lyh + lenyh
      lewtn = lsavf + n
      lacor = lewtn + n
      if (lewtn .gt. lewt) return
      do 30 i = 1,n
 30     rwork(i+lewtn-1) = rwork(i+lewt-1)
      lewt = lewtn
      return
      end
