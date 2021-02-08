
      integer istate, iopt, itask, itol, liw, lrw, mf, neq
      parameter(liw = 100)
      integer iwork(liw)
      parameter(lrw = 1000, maxeq = 1000)
      real*8 atol, rtol, t, tout 
      real*8 rwork(lrw), y(maxeq)
      logical first_step
      external fjac

      neq = 1
      y(1) = 1.
      t = 0.
      tout = 100.
      itol = 1
      rtol = 1.e-6
      atol = 1.e-6
      itask = 1
      first_step = .true.
      iopt = 0
      mf = 222
      call lsodes(neq, y, t, tout, itol, rtol, atol, itask, first_step, iopt, rwork, lrw, iwork, liw, fjac, mf)
      return
      end

      subroutine rhs(neq, t, y, ydot)
      double precision y, t, ydot
      dimension y(neq), ydot(neq)
      ydot(1) = y(1)
      return
      end

      subroutine fjac(neq, t, y, j, ian, jan, pdj)
      implicit real*8 (a-h,o-z)
      dimension y(1), ian(1), jan(1), pdj(1)
      return
      end
