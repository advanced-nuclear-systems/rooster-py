!
! intdy computes interpolated values of the k-th derivative of the dependent variable vector y, and stores it in dky.
! this routine is called within the package with k = 0 and t = tout, but may also be called by the user for any k up to the current order.
! (see detailed instructions in the usage documentation.)
!
! the computed values in dky are gotten by interpolation using the nordsieck history array yh.
! this array corresponds uniquely to a vector-valued polynomial of degree nqcur or less, and dky is set to the k-th derivative of this polynomial at t.
! the formula for dky is..
!              q
!  dky(i)  =  sum  c(j,k) * (t - tn)**(j-k) * h**(-j) * yh(i,j+1)
!             j=k
! where  c(j,k) = j*(j-1)*...*(j-k+1), q = nqcur, tn = tcur, h = hcur.
! the quantities  nq = nqcur, l = nq+1, n = neq, tn, and h are communicated by common.  the above sum is done in reverse order.
! iflag is returned negative if either k or t is out of bounds.

      subroutine intdy (t, k, yh, nnyh, dky, iflag)

      integer k, nnyh, iflag
      integer i, ic, j, jb, jb2, jj, jj1, jp1
      double precision t, yh, dky
      double precision c, r, s, tp
      dimension yh(nnyh,1), dky(1)

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

      iflag = 0
      if(k .lt. 0 .or. k .gt. nq)then
         call xerrwv("intdy--  k (=i1) illegal      ", 30, 51, 0, 1, k, 0, 0, 0.0d0, 0.0d0)
         iflag = -1
         RETURN
      end if
      tp = tn - hu -  100.0d0*uround*(tn + hu)
      if((t-tp)*(t-tn) .gt. 0.0d0)then
         call xerrwv("intdy--  t (=r1) illegal      ", 30, 52, 0, 0, 0, 0, 1, t, 0.0d0)
         call xerrwv("      t not in interval tcur - hu (= r1) to tcur (=r2)      ", 60, 52, 0, 0, 0, 0, 2, tp, tn)
         iflag = -2
         RETURN
      end if

      s = (t - tn)/h
      ic = 1
      if(k .ne. 0)then
         jj1 = l - k
         do jj = jj1,nq
            ic = ic*jj
         end do
      end if
      c = dfloat(ic)
      do i = 1,n
         dky(i) = c*yh(i,l)
      end do
      if(k .ne. nq)then
         jb2 = nq - k
         do jb = 1,jb2
            j = nq - jb
            jp1 = j + 1
            ic = 1
            if(k .ne. 0)then
               jj1 = jp1 - k
               do jj = jj1,j
                  ic = ic*jj
               end do
            end if
            c = dfloat(ic)
            do i = 1,n
               dky(i) = c*yh(i,jp1) + s*dky(i)
            end do
         end do
         if (k .eq. 0) RETURN
      end if
      r = h**(-k)
      do i = 1,n
         dky(i) = r*dky(i)
      end do
      RETURN
      end