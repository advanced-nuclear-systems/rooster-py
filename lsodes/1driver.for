!
      double precision atol, rtol, rwork, t, tout, y
      dimension y(1), rwork(500), iwork(300)
      data lrw/500/, liw/300/

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

      neq = 1
      y(1) = 1.0d0
      t = 0.0d0
      tout = 0.001d0
      itol = 1
      rtol = 1.0d-4
      atol = 1.0d-6
      itask = 1
      istate = 1
      iopt = 0
      mf = 121
      do iout = 1,1
        call lsodes(neq, y, t, tout, itol, rtol, atol, itask, istate, iopt, rwork, lrw, iwork, liw, mf)
        write(*,*)t,(y(i),i=1,neq)
        if(istate .lt. 0)then
           write(6,90)istate
 90        format(///22h error halt.. istate =,i3)
           stop
        end if
        tout = tout*10.0d0
      end do
      lenrw = iwork(17)
      leniw = iwork(18)
      nst = iwork(11)
      nfe = iwork(12)
      nje = iwork(13)
      nlu = iwork(21)
      nnz = iwork(19)
      nnzlu = iwork(25) + iwork(26) + neq
      write (6,70) lenrw,leniw,nst,nfe,nje,nlu,nnz,nnzlu
 70   format(//22h required rwork size =,i4,15h   iwork size =,i4/
     +   12h no. steps =,i4,12h   no. f-s =,i4,12h   no. j-s =,i4,
     +   13h   no. lu-s =,i4/23h no. of nonzeros in j =,i5,
     +   26h   no. of nonzeros in lu =,i5)
      stop
      end
     
      subroutine rhs(neq, t, y, ydot)
      double precision t, y(neq), ydot(neq)

      ydot(1) = y(1)
      return
      end
     
      subroutine fjac(neq, t, y, j, ia, ja, pdj)
      double precision t, y, pdj
      dimension y(100), ia(100), ja(100), pdj(100)
      return
      end
