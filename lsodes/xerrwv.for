      subroutine xerrwv (msg, nmes, nerr, level, ni, i1, i2, nr, r1, r2)
      character msg
      integer nmes, nerr, level, ni, i1, i2, nr, i, lun, lunit, mesflg, ncpw, nch, nwds
      double precision r1, r2
      dimension msg(nmes)
      common /eh0001/ mesflg, lunit

      mesflg=1
      lunit=6
      data ncpw/4/
      if(mesflg .ne. 0)then
         lun = lunit
         nch = min0(nmes,60)
         nwds = nch/ncpw
         if (nch .ne. nwds*ncpw) nwds = nwds + 1
         write (lun, *) msg
         if (ni .eq. 1) write (lun, 20) i1
 20      format(6x,23hin above message,  i1 =,i10)
         if (ni .eq. 2) write (lun, 30) i1,i2
 30      format(6x,23hin above message,  i1 =,i10,3x,4hi2 =,i10)
         if (nr .eq. 1) write (lun, 40) r1
 40      format(6x,23hin above message,  r1 =,d21.13)
         if (nr .eq. 2) write (lun, 50) r1,r2
 50      format(6x,15hin above,  r1 =,d21.13,3x,4hr2 =,d21.13)
      end if
      if (level .ne. 2) return
      stop
      end
