!-----------------------------------------------------------------------
! this is the march 30, 1987 version of
! lsodes.. livermore solver for ordinary differential equations
!          with general sparse jacobian matrices.
! this version is in double precision.
!
! lsodes solves the initial value problem for stiff or nonstiff
! systems of first order ode-s,
!     dy/dt = f(t,y) ,  or, in component form,
!     dy(i)/dt = f(i) = f(i,t,y(1),y(2),...,y(neq)) (i = 1,...,neq).
! lsodes is a variant of the lsode package, and is intended for
! problems in which the jacobian matrix df/dy has an arbitrary
! sparse structure (when the problem is stiff).
!
! authors..      alan c. hindmarsh,
!                computing and mathematics research division, l-316
!                lawrence livermore national laboratory
!                livermore, ca 94550.
!
! and            andrew h. sherman
!                j. s. nolen and associates
!                houston, tx 77084
!-----------------------------------------------------------------------
! references..
! 1.  alan c. hindmarsh,  odepack, a systematized collection of ode
!     solvers, in scientific computing, r. s. stepleman et al. (eds.),
!     north-holland, amsterdam, 1983, pp. 55-64.
!
! 2.  s. c. eisenstat, m. c. gursky, m. h. schultz, and a. h. sherman,
!     yale sparse matrix package.. i. the symmetric codes,
!     int. j. num. meth. eng., 18 (1982), pp. 1145-1151.
!
! 3.  s. c. eisenstat, m. c. gursky, m. h. schultz, and a. h. sherman,
!     yale sparse matrix package.. ii. the nonsymmetric codes,
!     research report no. 114, dept. of computer sciences, yale
!     university, 1977.
!-----------------------------------------------------------------------
! summary of usage.
!
! communication between the user and the lsodes package, for normal
! situations, is summarized here.  this summary describes only a subset
! of the full set of options available.  see the full description for
! details, including optional communication, nonstandard options,
! and instructions for special situations.  see also the example
! problem (with program and output) following this summary.
!
! a. first provide a subroutine of the form..
!               subroutine f (neq, t, y, ydot)
!               dimension y(neq), ydot(neq)
! which supplies the vector function f by loading ydot(i) with f(i).
!
! b. next determine (or guess) whether or not the problem is stiff.
! stiffness occurs when the jacobian matrix df/dy has an eigenvalue
! whose real part is negative and large in magnitude, compared to the
! reciprocal of the t span of interest.  if the problem is nonstiff,
! use a method flag mf = 10.  if it is stiff, there are two standard
! for the method flag, mf = 121 and mf = 222.  in both cases, lsodes
! requires the jacobian matrix in some form, and it treats this matrix
! in general sparse form, with sparsity structure determined internally.
! (for options where the user supplies the sparsity structure, see
! the full description of mf below.)
!
! c. if the problem is stiff, you are encouraged to supply the jacobian
! directly (mf = 121), but if this is not feasible, lsodes will
! compute it internally by difference quotients (mf = 222).
! if you are supplying the jacobian, provide a subroutine of the form..
!               subroutine jac (neq, t, y, j, ian, jan, pdj)
!               dimension y(1), ian(1), jan(1), pdj(1)
! here neq, t, y, and j are input arguments, and the jac routine is to
! load the array pdj (of length neq) with the j-th column of df/dy.
! i.e., load pdj(i) with df(i)/dy(j) for all relevant values of i.
! the arguments ian and jan should be ignored for normal situations.
! lsodes will call the jac routine with j = 1,2,...,neq.
! only nonzero elements need be loaded.  usually, a crude approximation
! to df/dy, possibly with fewer nonzero elements, will suffice.
!
! d. write a main program which calls subroutine lsodes once for
! each point at which answers are desired.  this should also provide
! for possible use of logical unit 6 for output of error messages
! by lsodes.  on the first call to lsodes, supply arguments as follows..
! f      = name of subroutine for right-hand side vector f.
!          this name must be declared external in calling program.
! neq    = number of first order ode-s.
! y      = array of initial values, of length neq.
! t      = the initial value of the independent variable.
! tout   = first point where output is desired (.ne. t).
! itol   = 1 or 2 according as atol (below) is a scalar or array.
! rtol   = relative tolerance parameter (scalar).
! atol   = absolute tolerance parameter (scalar or array).
!          the estimated local error in y(i) will be controlled so as
!          to be roughly less (in magnitude) than
!             ewt(i) = rtol*abs(y(i)) + atol     if itol = 1, or
!             ewt(i) = rtol*abs(y(i)) + atol(i)  if itol = 2.
!          thus the local error test passes if, in each component,
!          either the absolute error is less than atol (or atol(i)),
!          or the relative error is less than rtol.
!          use rtol = 0.0 for pure absolute error control, and
!          use atol = 0.0 (or atol(i) = 0.0) for pure relative error
!          control.  caution.. actual (global) errors may exceed these
!          local tolerances, so choose them conservatively.
! itask  = 1 for normal computation of output values of y at t = tout.
! istate = integer flag (input and output).  set istate = 1.
! iopt   = 0 to indicate no optional inputs used.
! rwork  = real work array of length at least..
!             20 + 16*neq            for mf = 10,
!             20 + (2 + 1./2)*nnz + (11 + 9./2)*neq
!                                    for mf = 121 or 222,
!          where..
!          nnz    = the number of nonzero elements in the sparse
!                   jacobian (if this is unknown, use an estimate), and
!          in any case, the required size of rwork cannot generally
!          be predicted in advance if mf = 121 or 222, and the value
!          above is a rough estimate of a crude lower bound.  some
!          experimentation with this size may be necessary.
!          (when known, the correct required length is an optional
!          output, available in iwork(17).)
! lrw    = declared length of rwork (in user-s dimension).
! iwork  = integer work array of length at least 30.
! liw    = declared length of iwork (in user-s dimension).
! jac    = name of subroutine for jacobian matrix (mf = 121).
!          if used, this name must be declared external in calling
!          program.  if not used, pass a dummy name.
! mf     = method flag.  standard values are..
!          10  for nonstiff (adams) method, no jacobian used.
!          121 for stiff (bdf) method, user-supplied sparse jacobian.
!          222 for stiff method, internally generated sparse jacobian.
! note that the main program must declare arrays y, rwork, iwork,
! and possibly atol.
!
! e. the output from the first call (or any call) is..
!      y = array of computed values of y(t) vector.
!      t = corresponding value of independent variable (normally tout).
! istate = 2  if lsodes was successful, negative otherwise.
!          -1 means excess work done on this call (perhaps wrong mf).
!          -2 means excess accuracy requested (tolerances too small).
!          -3 means illegal input detected (see printed message).
!          -4 means repeated error test failures (check all inputs).
!          -5 means repeated convergence failures (perhaps bad jacobian
!             supplied or wrong choice of mf or tolerances).
!          -6 means error weight became zero during problem. (solution
!             component i vanished, and atol or atol(i) = 0.)
!          -7 means a fatal error return flag came from the sparse
!             solver cdrv by way of prjs or slss.  should never happen.
!          a return with istate = -1, -4, or -5 may result from using
!          an inappropriate sparsity structure, one that is quite
!          different from the initial structure.  consider calling
!          lsodes again with istate = 3 to force the structure to be
!          reevaluated.  see the full description of istate below.
!
! f. to continue the integration after a successful return, simply
! reset tout and call lsodes again.  no other parameters need be reset.
!
!-----------------------------------------------------------------------
! example problem.
!
! the following is a simple example problem, with the coding
! needed for its solution by lsodes.  the problem is from chemical
! kinetics, and consists of the following 12 rate equations..
!    dy1/dt  = -rk1*y1
!    dy2/dt  = rk1*y1 + rk11*rk14*y4 + rk19*rk14*y5
!                - rk3*y2*y3 - rk15*y2*y12 - rk2*y2
!    dy3/dt  = rk2*y2 - rk5*y3 - rk3*y2*y3 - rk7*y10*y3
!                + rk11*rk14*y4 + rk12*rk14*y6
!    dy4/dt  = rk3*y2*y3 - rk11*rk14*y4 - rk4*y4
!    dy5/dt  = rk15*y2*y12 - rk19*rk14*y5 - rk16*y5
!    dy6/dt  = rk7*y10*y3 - rk12*rk14*y6 - rk8*y6
!    dy7/dt  = rk17*y10*y12 - rk20*rk14*y7 - rk18*y7
!    dy8/dt  = rk9*y10 - rk13*rk14*y8 - rk10*y8
!    dy9/dt  = rk4*y4 + rk16*y5 + rk8*y6 + rk18*y7
!    dy10/dt = rk5*y3 + rk12*rk14*y6 + rk20*rk14*y7
!                + rk13*rk14*y8 - rk7*y10*y3 - rk17*y10*y12
!                - rk6*y10 - rk9*y10
!    dy11/dt = rk10*y8
!    dy12/dt = rk6*y10 + rk19*rk14*y5 + rk20*rk14*y7
!                - rk15*y2*y12 - rk17*y10*y12
!
! with rk1 = rk5 = 0.1,  rk4 = rk8 = rk16 = rk18 = 2.5,
!      rk10 = 5.0,  rk2 = rk6 = 10.0,  rk14 = 30.0,
!      rk3 = rk7 = rk9 = rk11 = rk12 = rk13 = rk19 = rk20 = 50.0,
!      rk15 = rk17 = 100.0.
!
! the t interval is from 0 to 1000, and the initial conditions
! are y1 = 1, y2 = y3 = ... = y12 = 0.  the problem is stiff.
!
! the following coding solves this problem with lsodes, using mf = 121
! and printing results at t = .1, 1., 10., 100., 1000.  it uses
! itol = 1 and mixed relative/absolute tolerance controls.
! during the run and at the end, statistical quantities of interest
! are printed (see optional outputs in the full description below).
!
!     external fex, jex
!     double precision atol, rtol, rwork, t, tout, y
!     dimension y(12), rwork(500), iwork(30)
!     data lrw/500/, liw/30/
!     neq = 12
!     do 10 i = 1,neq
! 10    y(i) = 0.0d0
!     y(1) = 1.0d0
!     t = 0.0d0
!     tout = 0.1d0
!     itol = 1
!     rtol = 1.0d-4
!     atol = 1.0d-6
!     itask = 1
!     istate = 1
!     iopt = 0
!     mf = 121
!     do 40 iout = 1,5
!       call lsodes (fex, neq, y, t, tout, itol, rtol, atol,
!    1     itask, istate, iopt, rwork, lrw, iwork, liw, jex, mf)
!       write(6,30)t,iwork(11),rwork(11),(y(i),i=1,neq)
! 30    format(//7h at t =,e11.3,4x,
!    1    12h no. steps =,i5,4x,12h last step =,e11.3/
!    2    13h  y array =  ,4e14.5/13x,4e14.5/13x,4e14.5)
!       if(istate .lt. 0) go to 80
!       tout = tout*10.0d0
! 40    continue
!     lenrw = iwork(17)
!     leniw = iwork(18)
!     nst = iwork(11)
!     nfe = iwork(12)
!     nje = iwork(13)
!     nlu = iwork(21)
!     nnz = iwork(19)
!     nnzlu = iwork(25) + iwork(26) + neq
!     write (6,70) lenrw,leniw,nst,nfe,nje,nlu,nnz,nnzlu
! 70  format(//22h required rwork size =,i4,15h   iwork size =,i4/
!    1   12h no. steps =,i4,12h   no. f-s =,i4,12h   no. j-s =,i4,
!    2   13h   no. lu-s =,i4/23h no. of nonzeros in j =,i5,
!    3   26h   no. of nonzeros in lu =,i5)
!     stop
! 80  write(6,90)istate
! 90  format(///22h error halt.. istate =,i3)
!     stop
!     end
!
!     subroutine fex (neq, t, y, ydot)
!     double precision t, y, ydot
!     double precision rk1, rk2, rk3, rk4, rk5, rk6, rk7, rk8, rk9,
!    1   rk10, rk11, rk12, rk13, rk14, rk15, rk16, rk17
!     dimension y(12), ydot(12)
!     data rk1/0.1d0/, rk2/10.0d0/, rk3/50.0d0/, rk4/2.5d0/, rk5/0.1d0/,
!    1   rk6/10.0d0/, rk7/50.0d0/, rk8/2.5d0/, rk9/50.0d0/, rk10/5.0d0/,
!    2   rk11/50.0d0/, rk12/50.0d0/, rk13/50.0d0/, rk14/30.0d0/,
!    3   rk15/100.0d0/, rk16/2.5d0/, rk17/100.0d0/, rk18/2.5d0/,
!    4   rk19/50.0d0/, rk20/50.0d0/
!     ydot(1)  = -rk1*y(1)
!     ydot(2)  = rk1*y(1) + rk11*rk14*y(4) + rk19*rk14*y(5)
!    1           - rk3*y(2)*y(3) - rk15*y(2)*y(12) - rk2*y(2)
!     ydot(3)  = rk2*y(2) - rk5*y(3) - rk3*y(2)*y(3) - rk7*y(10)*y(3)
!    1           + rk11*rk14*y(4) + rk12*rk14*y(6)
!     ydot(4)  = rk3*y(2)*y(3) - rk11*rk14*y(4) - rk4*y(4)
!     ydot(5)  = rk15*y(2)*y(12) - rk19*rk14*y(5) - rk16*y(5)
!     ydot(6)  = rk7*y(10)*y(3) - rk12*rk14*y(6) - rk8*y(6)
!     ydot(7)  = rk17*y(10)*y(12) - rk20*rk14*y(7) - rk18*y(7)
!     ydot(8)  = rk9*y(10) - rk13*rk14*y(8) - rk10*y(8)
!     ydot(9)  = rk4*y(4) + rk16*y(5) + rk8*y(6) + rk18*y(7)
!     ydot(10) = rk5*y(3) + rk12*rk14*y(6) + rk20*rk14*y(7)
!    1           + rk13*rk14*y(8) - rk7*y(10)*y(3) - rk17*y(10)*y(12)
!    2           - rk6*y(10) - rk9*y(10)
!     ydot(11) = rk10*y(8)
!     ydot(12) = rk6*y(10) + rk19*rk14*y(5) + rk20*rk14*y(7)
!    1           - rk15*y(2)*y(12) - rk17*y(10)*y(12)
!     return
!     end
!
!     subroutine jex (neq, t, y, j, ia, ja, pdj)
!     double precision t, y, pdj
!     double precision rk1, rk2, rk3, rk4, rk5, rk6, rk7, rk8, rk9,
!    1   rk10, rk11, rk12, rk13, rk14, rk15, rk16, rk17
!     dimension y(1), ia(1), ja(1), pdj(1)
!     data rk1/0.1d0/, rk2/10.0d0/, rk3/50.0d0/, rk4/2.5d0/, rk5/0.1d0/,
!    1   rk6/10.0d0/, rk7/50.0d0/, rk8/2.5d0/, rk9/50.0d0/, rk10/5.0d0/,
!    2   rk11/50.0d0/, rk12/50.0d0/, rk13/50.0d0/, rk14/30.0d0/,
!    3   rk15/100.0d0/, rk16/2.5d0/, rk17/100.0d0/, rk18/2.5d0/,
!    4   rk19/50.0d0/, rk20/50.0d0/
!     go to (1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12), j
! 1   pdj(1) = -rk1
!     pdj(2) = rk1
!     return
! 2   pdj(2) = -rk3*y(3) - rk15*y(12) - rk2
!     pdj(3) = rk2 - rk3*y(3)
!     pdj(4) = rk3*y(3)
!     pdj(5) = rk15*y(12)
!     pdj(12) = -rk15*y(12)
!     return
! 3   pdj(2) = -rk3*y(2)
!     pdj(3) = -rk5 - rk3*y(2) - rk7*y(10)
!     pdj(4) = rk3*y(2)
!     pdj(6) = rk7*y(10)
!     pdj(10) = rk5 - rk7*y(10)
!     return
! 4   pdj(2) = rk11*rk14
!     pdj(3) = rk11*rk14
!     pdj(4) = -rk11*rk14 - rk4
!     pdj(9) = rk4
!     return
! 5   pdj(2) = rk19*rk14
!     pdj(5) = -rk19*rk14 - rk16
!     pdj(9) = rk16
!     pdj(12) = rk19*rk14
!     return
! 6   pdj(3) = rk12*rk14
!     pdj(6) = -rk12*rk14 - rk8
!     pdj(9) = rk8
!     pdj(10) = rk12*rk14
!     return
! 7   pdj(7) = -rk20*rk14 - rk18
!     pdj(9) = rk18
!     pdj(10) = rk20*rk14
!     pdj(12) = rk20*rk14
!     return
! 8   pdj(8) = -rk13*rk14 - rk10
!     pdj(10) = rk13*rk14
!     pdj(11) = rk10
! 9   return
! 10  pdj(3) = -rk7*y(3)
!     pdj(6) = rk7*y(3)
!     pdj(7) = rk17*y(12)
!     pdj(8) = rk9
!     pdj(10) = -rk7*y(3) - rk17*y(12) - rk6 - rk9
!     pdj(12) = rk6 - rk17*y(12)
! 11  return
! 12  pdj(2) = -rk15*y(2)
!     pdj(5) = rk15*y(2)
!     pdj(7) = rk17*y(10)
!     pdj(10) = -rk17*y(10)
!     pdj(12) = -rk15*y(2) - rk17*y(10)
!     return
!     end
!
! the output of this program (on a cray-1 in single precision)
! is as follows..
!
!
! at t =  1.000e-01     no. steps =   12     last step =  1.515e-02
!  y array =     9.90050e-01   6.28228e-03   3.65313e-03   7.51934e-07
!                1.12167e-09   1.18458e-09   1.77291e-12   3.26476e-07
!                5.46720e-08   9.99500e-06   4.48483e-08   2.76398e-06
!
!
! at t =  1.000e+00     no. steps =   33     last step =  7.880e-02
!  y array =     9.04837e-01   9.13105e-03   8.20622e-02   2.49177e-05
!                1.85055e-06   1.96797e-06   1.46157e-07   2.39557e-05
!                3.26306e-05   7.21621e-04   5.06433e-05   3.05010e-03
!
!
! at t =  1.000e+01     no. steps =   48     last step =  1.239e+00
!  y array =     3.67876e-01   3.68958e-03   3.65133e-01   4.48325e-05
!                6.10798e-05   4.33148e-05   5.90211e-05   1.18449e-04
!                3.15235e-03   3.56531e-03   4.15520e-03   2.48741e-01
!
!
! at t =  1.000e+02     no. steps =   91     last step =  3.764e+00
!  y array =     4.44981e-05   4.42666e-07   4.47273e-04  -3.53257e-11
!                2.81577e-08  -9.67741e-11   2.77615e-07   1.45322e-07
!                1.56230e-02   4.37394e-06   1.60104e-02   9.52246e-01
!
!
! at t =  1.000e+03     no. steps =  111     last step =  4.156e+02
!  y array =    -2.65492e-13   2.60539e-14  -8.59563e-12   6.29355e-14
!               -1.78066e-13   5.71471e-13  -1.47561e-12   4.58078e-15
!                1.56314e-02   1.37878e-13   1.60184e-02   9.52719e-01
!
!
! required rwork size = 442   iwork size =  30
! no. steps = 111   no. f-s = 142   no. j-s =   2   no. lu-s =  20
! no. of nonzeros in j =   44   no. of nonzeros in lu =   50
!-----------------------------------------------------------------------
! full description of user interface to lsodes.
!
! the user interface to lsodes consists of the following parts.
!
! i.   the call sequence to subroutine lsodes, which is a driver
!      routine for the solver.  this includes descriptions of both
!      the call sequence arguments and of user-supplied routines.
!      following these descriptions is a description of
!      optional inputs available through the call sequence, and then
!      a description of optional outputs (in the work arrays).
!
! ii.  descriptions of other routines in the lsodes package that may be
!      (optionally) called by the user.  these provide the ability to
!      alter error message handling, save and restore the internal
!      common, and obtain specified derivatives of the solution y(t).
!
! iii. descriptions of common blocks to be declared in overlay
!      or similar environments, or to be saved when doing an interrupt
!      of the problem and continued solution later.
!
! iv.  description of two routines in the lsodes package, either of
!      which the user may replace with his own version, if desired.
!      these relate to the measurement of errors.
!
!-----------------------------------------------------------------------
! part i.  call sequence.
!
! the call sequence parameters used for input only are
!     f, neq, tout, itol, rtol, atol, itask, iopt, lrw, liw, jac, mf,
! and those used for both input and output are
!     y, t, istate.
! the work arrays rwork and iwork are also used for conditional and
! optional inputs and optional outputs.  (the term output here refers
! to the return from subroutine lsodes to the user-s calling program.)
!
! the legality of input parameters will be thoroughly checked on the
! initial call for the problem, but not checked thereafter unless a
! change in input parameters is flagged by istate = 3 on input.
!
! the descriptions of the call arguments are as follows.
!
! f      = the name of the user-supplied subroutine defining the
!          ode system.  the system must be put in the first-order
!          form dy/dt = f(t,y), where f is a vector-valued function
!          of the scalar t and the vector y.  subroutine f is to
!          compute the function f.  it is to have the form
!               subroutine f (neq, t, y, ydot)
!               dimension y(1), ydot(1)
!          where neq, t, and y are input, and the array ydot = f(t,y)
!          is output.  y and ydot are arrays of length neq.
!          (in the dimension statement above, 1 is a dummy
!          dimension.. it can be replaced by any value.)
!          subroutine f should not alter y(1),...,y(neq).
!          f must be declared external in the calling program.
!
!          subroutine f may access user-defined quantities in
!          neq(2),... and/or in y(neq(1)+1),... if neq is an array
!          (dimensioned in f) and/or y has length exceeding neq(1).
!          see the descriptions of neq and y below.
!
!          if quantities computed in the f routine are needed
!          externally to lsodes, an extra call to f should be made
!          for this purpose, for consistent and accurate results.
!          if only the derivative dy/dt is needed, use intdy instead.
!
! neq    = the size of the ode system (number of first order
!          ordinary differential equations).  used only for input.
!          neq may be decreased, but not increased, during the problem.
!          if neq is decreased (with istate = 3 on input), the
!          remaining components of y should be left undisturbed, if
!          these are to be accessed in f and/or jac.
!
!          normally, neq is a scalar, and it is generally referred to
!          as a scalar in this user interface description.  however,
!          neq may be an array, with neq(1) set to the system size.
!          (the lsodes package accesses only neq(1).)  in either case,
!          this parameter is passed as the neq argument in all calls
!          to f and jac.  hence, if it is an array, locations
!          neq(2),... may be used to store other integer data and pass
!          it to f and/or jac.  subroutines f and/or jac must include
!          neq in a dimension statement in that case.
!
! y      = a real array for the vector of dependent variables, of
!          length neq or more.  used for both input and output on the
!          first call (istate = 1), and only for output on other calls.
!          on the first call, y must contain the vector of initial
!          values.  on output, y contains the computed solution vector,
!          evaluated at t.  if desired, the y array may be used
!          for other purposes between calls to the solver.
!
!          this array is passed as the y argument in all calls to
!          f and jac.  hence its length may exceed neq, and locations
!          y(neq+1),... may be used to store other real data and
!          pass it to f and/or jac.  (the lsodes package accesses only
!          y(1),...,y(neq).)
!
! t      = the independent variable.  on input, t is used only on the
!          first call, as the initial point of the integration.
!          on output, after each call, t is the value at which a
!          computed solution y is evaluated (usually the same as tout).
!          on an error return, t is the farthest point reached.
!
! tout   = the next value of t at which a computed solution is desired.
!          used only for input.
!
!          when starting the problem (istate = 1), tout may be equal
!          to t for one call, then should .ne. t for the next call.
!          for the initial t, an input value of tout .ne. t is used
!          in order to determine the direction of the integration
!          (i.e. the algebraic sign of the step sizes) and the rough
!          scale of the problem.  integration in either direction
!          (forward or backward in t) is permitted.
!
!          if itask = 2 or 5 (one-step modes), tout is ignored after
!          the first call (i.e. the first call with tout .ne. t).
!          otherwise, tout is required on every call.
!
!          if itask = 1, 3, or 4, the values of tout need not be
!          monotone, but a value of tout which backs up is limited
!          to the current internal t interval, whose endpoints are
!          tcur - hu and tcur (see optional outputs, below, for
!          tcur and hu).
!
! itol   = an indicator for the type of error control.  see
!          description below under atol.  used only for input.
!
! rtol   = a relative error tolerance parameter, either a scalar or
!          an array of length neq.  see description below under atol.
!          input only.
!
! atol   = an absolute error tolerance parameter, either a scalar or
!          an array of length neq.  input only.
!
!             the input parameters itol, rtol, and atol determine
!          the error control performed by the solver.  the solver will
!          control the vector e = (e(i)) of estimated local errors
!          in y, according to an inequality of the form
!                      rms-norm of ( e(i)/ewt(i) )   .le.   1,
!          where       ewt(i) = rtol(i)*abs(y(i)) + atol(i),
!          and the rms-norm (root-mean-square norm) here is
!          rms-norm(v) = sqrt(sum v(i)**2 / neq).  here ewt = (ewt(i))
!          is a vector of weights which must always be positive, and
!          the values of rtol and atol should all be non-negative.
!          the following table gives the types (scalar/array) of
!          rtol and atol, and the corresponding form of ewt(i).
!
!             itol    rtol       atol          ewt(i)
!              1     scalar     scalar     rtol*abs(y(i)) + atol
!              2     scalar     array      rtol*abs(y(i)) + atol(i)
!              3     array      scalar     rtol(i)*abs(y(i)) + atol
!              4     array      array      rtol(i)*abs(y(i)) + atol(i)
!
!          when either of these parameters is a scalar, it need not
!          be dimensioned in the user-s calling program.
!
!          if none of the above choices (with itol, rtol, and atol
!          fixed throughout the problem) is suitable, more general
!          error controls can be obtained by substituting
!          user-supplied routines for the setting of ewt and/or for
!          the norm calculation.  see part iv below.
!
!          if global errors are to be estimated by making a repeated
!          run on the same problem with smaller tolerances, then all
!          components of rtol and atol (i.e. of ewt) should be scaled
!          down uniformly.
!
! itask  = an index specifying the task to be performed.
!          input only.  itask has the following values and meanings.
!          1  means normal computation of output values of y(t) at
!             t = tout (by overshooting and interpolating).
!          2  means take one step only and return.
!          3  means stop at the first internal mesh point at or
!             beyond t = tout and return.
!          4  means normal computation of output values of y(t) at
!             t = tout but without overshooting t = tcrit.
!             tcrit must be input as rwork(1).  tcrit may be equal to
!             or beyond tout, but not behind it in the direction of
!             integration.  this option is useful if the problem
!             has a singularity at or beyond t = tcrit.
!          5  means take one step, without passing tcrit, and return.
!             tcrit must be input as rwork(1).
!
!          note..  if itask = 4 or 5 and the solver reaches tcrit
!          (within roundoff), it will return t = tcrit (exactly) to
!          indicate this (unless itask = 4 and tout comes before tcrit,
!          in which case answers at t = tout are returned first).
!
! istate = an index used for input and output to specify the
!          the state of the calculation.
!
!          on input, the values of istate are as follows.
!          1  means this is the first call for the problem
!             (initializations will be done).  see note below.
!          2  means this is not the first call, and the calculation
!             is to continue normally, with no change in any input
!             parameters except possibly tout and itask.
!             (if itol, rtol, and/or atol are changed between calls
!             with istate = 2, the new values will be used but not
!             tested for legality.)
!          3  means this is not the first call, and the
!             calculation is to continue normally, but with
!             a change in input parameters other than
!             tout and itask.  changes are allowed in
!             neq, itol, rtol, atol, iopt, lrw, liw, mf,
!             the conditional inputs ia and ja,
!             and any of the optional inputs except h0.
!             in particular, if miter = 1 or 2, a call with istate = 3
!             will cause the sparsity structure of the problem to be
!             recomputed (or reread from ia and ja if moss = 0).
!          note..  a preliminary call with tout = t is not counted
!          as a first call here, as no initialization or checking of
!          input is done.  (such a call is sometimes useful for the
!          purpose of outputting the initial conditions.)
!          thus the first call for which tout .ne. t requires
!          istate = 1 on input.
!
!          on output, istate has the following values and meanings.
!           1  means nothing was done, as tout was equal to t with
!              istate = 1 on input.  (however, an internal counter was
!              set to detect and prevent repeated calls of this type.)
!           2  means the integration was performed successfully.
!          -1  means an excessive amount of work (more than mxstep
!              steps) was done on this call, before completing the
!              requested task, but the integration was otherwise
!              successful as far as t.  (mxstep is an optional input
!              and is normally 500.)  to continue, the user may
!              simply reset istate to a value .gt. 1 and call again
!              (the excess work step counter will be reset to 0).
!              in addition, the user may increase mxstep to avoid
!              this error return (see below on optional inputs).
!          -2  means too much accuracy was requested for the precision
!              of the machine being used.  this was detected before
!              completing the requested task, but the integration
!              was successful as far as t.  to continue, the tolerance
!              parameters must be reset, and istate must be set
!              to 3.  the optional output tolsf may be used for this
!              purpose.  (note.. if this condition is detected before
!              taking any steps, then an illegal input return
!              (istate = -3) occurs instead.)
!          -3  means illegal input was detected, before taking any
!              integration steps.  see written message for details.
!              note..  if the solver detects an infinite loop of calls
!              to the solver with illegal input, it will cause
!              the run to stop.
!          -4  means there were repeated error test failures on
!              one attempted step, before completing the requested
!              task, but the integration was successful as far as t.
!              the problem may have a singularity, or the input
!              may be inappropriate.
!          -5  means there were repeated convergence test failures on
!              one attempted step, before completing the requested
!              task, but the integration was successful as far as t.
!              this may be caused by an inaccurate jacobian matrix,
!              if one is being used.
!          -6  means ewt(i) became zero for some i during the
!              integration.  pure relative error control (atol(i)=0.0)
!              was requested on a variable which has now vanished.
!              the integration was successful as far as t.
!          -7  means a fatal error return flag came from the sparse
!              solver cdrv by way of prjs or slss (numerical
!              factorization or backsolve).  this should never happen.
!              the integration was successful as far as t.
!
!          note.. an error return with istate = -1, -4, or -5 and with
!          miter = 1 or 2 may mean that the sparsity structure of the
!          problem has changed significantly since it was last
!          determined (or input).  in that case, one can attempt to
!          complete the integration by setting istate = 3 on the next
!          call, so that a new structure determination is done.
!
!          note..  since the normal output value of istate is 2,
!          it does not need to be reset for normal continuation.
!          also, since a negative input value of istate will be
!          regarded as illegal, a negative output value requires the
!          user to change it, and possibly other inputs, before
!          calling the solver again.
!
! iopt   = an integer flag to specify whether or not any optional
!          inputs are being used on this call.  input only.
!          the optional inputs are listed separately below.
!          iopt = 0 means no optional inputs are being used.
!                   default values will be used in all cases.
!          iopt = 1 means one or more optional inputs are being used.
!
! rwork  = a work array used for a mixture of real (double precision)
!          and integer work space.
!          the length of rwork (in real words) must be at least
!             20 + nyh*(maxord + 1) + 3*neq + lwm    where
!          nyh    = the initial value of neq,
!          maxord = 12 (if meth = 1) or 5 (if meth = 2) (unless a
!                   smaller value is given as an optional input),
!          lwm = 0                                    if miter = 0,
!          lwm = 2*nnz + 2*neq + (nnz+9*neq)/2   if miter = 1,
!          lwm = 2*nnz + 2*neq + (nnz+10*neq)/2  if miter = 2,
!          lwm = neq + 2                              if miter = 3.
!          in the above formulas,
!          nnz    = number of nonzero elements in the jacobian matrix.
!          (see the mf description for meth and miter.)
!          thus if maxord has its default value and neq is constant,
!          the minimum length of rwork is..
!             20 + 16*neq        for mf = 10,
!             20 + 16*neq + lwm  for mf = 11, 111, 211, 12, 112, 212,
!             22 + 17*neq        for mf = 13,
!             20 +  9*neq        for mf = 20,
!             20 +  9*neq + lwm  for mf = 21, 121, 221, 22, 122, 222,
!             22 + 10*neq        for mf = 23.
!          if miter = 1 or 2, the above formula for lwm is only a
!          crude lower bound.  the required length of rwork cannot
!          be readily predicted in general, as it depends on the
!          sparsity structure of the problem.  some experimentation
!          may be necessary.
!
!          the first 20 words of rwork are reserved for conditional
!          and optional inputs and optional outputs.
!
!          the following word in rwork is a conditional input..
!            rwork(1) = tcrit = critical value of t which the solver
!                       is not to overshoot.  required if itask is
!                       4 or 5, and ignored otherwise.  (see itask.)
!
! lrw    = the length of the array rwork, as declared by the user.
!          (this will be checked by the solver.)
!
! iwork  = an integer work array.  the length of iwork must be at least
!             31 + neq + nnz   if moss = 0 and miter = 1 or 2, or
!             30               otherwise.
!          (nnz is the number of nonzero elements in df/dy.)
!
!          in lsodes, iwork is used only for conditional and
!          optional inputs and optional outputs.
!
!          the following two blocks of words in iwork are conditional
!          inputs, required if moss = 0 and miter = 1 or 2, but not
!          otherwise (see the description of mf for moss).
!            iwork(30+j) = ia(j)     (j=1,...,neq+1)
!            iwork(31+neq+k) = ja(k) (k=1,...,nnz)
!          the two arrays ia and ja describe the sparsity structure
!          to be assumed for the jacobian matrix.  ja contains the row
!          indices where nonzero elements occur, reading in columnwise
!          order, and ia contains the starting locations in ja of the
!          descriptions of columns 1,...,neq, in that order, with
!          ia(1) = 1.  thus, for each column index j = 1,...,neq, the
!          values of the row index i in column j where a nonzero
!          element may occur are given by
!            i = ja(k),  where   ia(j) .le. k .lt. ia(j+1).
!          if nnz is the total number of nonzero locations assumed,
!          then the length of the ja array is nnz, and ia(neq+1) must
!          be nnz + 1.  duplicate entries are not allowed.
!
! liw    = the length of the array iwork, as declared by the user.
!          (this will be checked by the solver.)
!
! note..  the work arrays must not be altered between calls to lsodes
! for the same problem, except possibly for the conditional and
! optional inputs, and except for the last 3*neq words of rwork.
! the latter space is used for internal scratch space, and so is
! available for use by the user outside lsodes between calls, if
! desired (but not for use by f or jac).
!
! jac    = name of user-supplied routine (miter = 1 or moss = 1) to
!          compute the jacobian matrix, df/dy, as a function of
!          the scalar t and the vector y.  it is to have the form
!               subroutine jac (neq, t, y, j, ian, jan, pdj)
!               dimension y(1), ian(1), jan(1), pdj(1)
!          where neq, t, y, j, ian, and jan are input, and the array
!          pdj, of length neq, is to be loaded with column j
!          of the jacobian on output.  thus df(i)/dy(j) is to be
!          loaded into pdj(i) for all relevant values of i.
!          here t and y have the same meaning as in subroutine f,
!          and j is a column index (1 to neq).  ian and jan are
!          undefined in calls to jac for structure determination
!          (moss = 1).  otherwise, ian and jan are structure
!          descriptors, as defined under optional outputs below, and
!          so can be used to determine the relevant row indices i, if
!          desired.  (in the dimension statement above, 1 is a
!          dummy dimension.. it can be replaced by any value.)
!               jac need not provide df/dy exactly.  a crude
!          approximation (possibly with greater sparsity) will do.
!               in any case, pdj is preset to zero by the solver,
!          so that only the nonzero elements need be loaded by jac.
!          calls to jac are made with j = 1,...,neq, in that order, and
!          each such set of calls is preceded by a call to f with the
!          same arguments neq, t, and y.  thus to gain some efficiency,
!          intermediate quantities shared by both calculations may be
!          saved in a user common block by f and not recomputed by jac,
!          if desired.  jac must not alter its input arguments.
!          jac must be declared external in the calling program.
!               subroutine jac may access user-defined quantities in
!          neq(2),... and y(neq(1)+1),... if neq is an array
!          (dimensioned in jac) and y has length exceeding neq(1).
!          see the descriptions of neq and y above.
!
! mf     = the method flag.  used only for input.
!          mf has three decimal digits-- moss, meth, miter--
!             mf = 100*moss + 10*meth + miter.
!          moss indicates the method to be used to obtain the sparsity
!          structure of the jacobian matrix if miter = 1 or 2..
!            moss = 0 means the user has supplied ia and ja
!                     (see descriptions under iwork above).
!            moss = 1 means the user has supplied jac (see below)
!                     and the structure will be obtained from neq
!                     initial calls to jac.
!            moss = 2 means the structure will be obtained from neq+1
!                     initial calls to f.
!          meth indicates the basic linear multistep method..
!            meth = 1 means the implicit adams method.
!            meth = 2 means the method based on backward
!                     differentiation formulas (bdf-s).
!          miter indicates the corrector iteration method..
!            miter = 0 means functional iteration (no jacobian matrix is involved).
!            miter = 1 means chord iteration with a user-supplied sparse jacobian, given by subroutine jac.
!            miter = 2 means chord iteration with an internally generated (difference quotient) sparse jacobian
!                      (using ngp extra calls to f per df/dy value, where ngp is an optional output described below.)
!            miter = 3 means chord iteration with an internally generated diagonal jacobian approximation.
!                      (using 1 extra call to f per df/dy evaluation).
!          if miter = 1 or moss = 1, the user must supply a subroutine
!          jac (the name is arbitrary) as described above under jac.
!          otherwise, a dummy argument can be used.
!
!          the standard choices for mf are..
!            mf = 10  for a nonstiff problem,
!            mf = 21 or 22 for a stiff problem with ia/ja supplied
!                     (21 if jac is supplied, 22 if not),
!            mf = 121 for a stiff problem with jac supplied,
!                     but not ia/ja,
!            mf = 222 for a stiff problem with neither ia/ja nor
!                     jac supplied.
!          the sparseness structure can be changed during the
!          problem by making a call to lsodes with istate = 3.
!-----------------------------------------------------------------------
! optional inputs.
!
! the following is a list of the optional inputs provided for in the
! call sequence.  (see also part ii.)  for each such input variable,
! this table lists its name as used in this documentation, its
! location in the call sequence, its meaning, and the default value.
! the use of any of these inputs requires iopt = 1, and in that
! case all of these inputs are examined.  a value of zero for any
! of these optional inputs will cause the default value to be used.
! thus to use a subset of the optional inputs, simply preload
! locations 5 to 10 in rwork and iwork to 0.0 and 0 respectively, and
! then set those of interest to nonzero values.
!
! name    location      meaning and default value
!
! h0      rwork(5)  the step size to be attempted on the first step.
!                   the default value is determined by the solver.
!
! hmax    rwork(6)  the maximum absolute step size allowed.
!                   the default value is infinite.
!
! hmin    rwork(7)  the minimum absolute step size allowed.
!                   the default value is 0.  (this lower bound is not
!                   enforced on the final step before reaching tcrit
!                   when itask = 4 or 5.)
!
! seth    rwork(8)  the element threshhold for sparsity determination
!                   when moss = 1 or 2.  if the absolute value of
!                   an estimated jacobian element is .le. seth, it
!                   will be assumed to be absent in the structure.
!                   the default value of seth is 0.
!
! maxord  iwork(5)  the maximum order to be allowed.  the default
!                   value is 12 if meth = 1, and 5 if meth = 2.
!                   if maxord exceeds the default value, it will
!                   be reduced to the default value.
!                   if maxord is changed during the problem, it may
!                   cause the current order to be reduced.
!
! mxstep  iwork(6)  maximum number of (internally defined) steps
!                   allowed during one call to the solver.
!                   the default value is 500.
!
! mxhnil  iwork(7)  maximum number of messages printed (per problem)
!                   warning that t + h = t on a step (h = step size).
!                   this must be positive to result in a non-default
!                   value.  the default value is 10.
!-----------------------------------------------------------------------
! optional outputs.
!
! as optional additional output from lsodes, the variables listed
! below are quantities related to the performance of lsodes
! which are available to the user.  these are communicated by way of
! the work arrays, but also have internal mnemonic names as shown.
! except where stated otherwise, all of these outputs are defined
! on any successful return from lsodes, and on any return with
! istate = -1, -2, -4, -5, or -6.  on an illegal input return
! (istate = -3), they will be unchanged from their existing values
! (if any), except possibly for tolsf, lenrw, and leniw.
! on any error return, outputs relevant to the error will be defined,
! as noted below.
!
! name    location      meaning
!
! hu      rwork(11) the step size in t last used (successfully).
!
! hcur    rwork(12) the step size to be attempted on the next step.
!
! tcur    rwork(13) the current value of the independent variable
!                   which the solver has actually reached, i.e. the
!                   current internal mesh point in t.  on output, tcur
!                   will always be at least as far as the argument
!                   t, but may be farther (if interpolation was done).
!
! tolsf   rwork(14) a tolerance scale factor, greater than 1.0,
!                   computed when a request for too much accuracy was
!                   detected (istate = -3 if detected at the start of
!                   the problem, istate = -2 otherwise).  if itol is
!                   left unaltered but rtol and atol are uniformly
!                   scaled up by a factor of tolsf for the next call,
!                   then the solver is deemed likely to succeed.
!                   (the user may also ignore tolsf and alter the
!                   tolerance parameters in any other way appropriate.)
!
! nst     iwork(11) the number of steps taken for the problem so far.
!
! nfe     iwork(12) the number of f evaluations for the problem so far,
!                   excluding those for structure determination
!                   (moss = 2).
!
! nje     iwork(13) the number of jacobian evaluations for the problem
!                   so far, excluding those for structure determination
!                   (moss = 1).
!
! nqu     iwork(14) the method order last used (successfully).
!
! nqcur   iwork(15) the order to be attempted on the next step.
!
! imxer   iwork(16) the index of the component of largest magnitude in
!                   the weighted local error vector ( e(i)/ewt(i) ),
!                   on an error return with istate = -4 or -5.
!
! lenrw   iwork(17) the length of rwork actually required.
!                   this is defined on normal returns and on an illegal
!                   input return for insufficient storage.
!
! leniw   iwork(18) the length of iwork actually required.
!                   this is defined on normal returns and on an illegal
!                   input return for insufficient storage.
!
! nnz     iwork(19) the number of nonzero elements in the jacobian
!                   matrix, including the diagonal (miter = 1 or 2).
!                   (this may differ from that given by ia(neq+1)-1
!                   if moss = 0, because of added diagonal entries.)
!
! ngp     iwork(20) the number of groups of column indices, used in
!                   difference quotient jacobian aproximations if
!                   miter = 2.  this is also the number of extra f
!                   evaluations needed for each jacobian evaluation.
!
! nlu     iwork(21) the number of sparse lu decompositions for the
!                   problem so far.
!
! lyh     iwork(22) the base address in rwork of the history array yh,
!                   described below in this list.
!
! ipian   iwork(23) the base address of the structure descriptor array
!                   ian, described below in this list.
!
! ipjan   iwork(24) the base address of the structure descriptor array
!                   jan, described below in this list.
!
! nzl     iwork(25) the number of nonzero elements in the strict lower
!                   triangle of the lu factorization used in the chord
!                   iteration (miter = 1 or 2).
!
! nzu     iwork(26) the number of nonzero elements in the strict upper
!                   triangle of the lu factorization used in the chord
!                   iteration (miter = 1 or 2).
!                   the total number of nonzeros in the factorization
!                   is therefore nzl + nzu + neq.
!
! the following four arrays are segments of the rwork array which
! may also be of interest to the user as optional outputs.
! for each array, the table below gives its internal name,
! its base address, and its description.
! for yh and acor, the base addresses are in rwork (a real array).
! the integer arrays ian and jan are to be obtained by declaring an
! integer array iwk and identifying iwk(1) with rwork(21), using either
! an equivalence statement or a subroutine call.  then the base
! addresses ipian (of ian) and ipjan (of jan) in iwk are to be obtained
! as optional outputs iwork(23) and iwork(24), respectively.
! thus ian(1) is iwk(ipian), etc.
!
! name    base address      description
!
! ian    ipian (in iwk)  structure descriptor array of size neq + 1.
! jan    ipjan (in iwk)  structure descriptor array of size nnz.
!         (see above)    ian and jan together describe the sparsity
!                        structure of the jacobian matrix, as used by
!                        lsodes when miter = 1 or 2.
!                        jan contains the row indices of the nonzero
!                        locations, reading in columnwise order, and
!                        ian contains the starting locations in jan of
!                        the descriptions of columns 1,...,neq, in
!                        that order, with ian(1) = 1.  thus for each
!                        j = 1,...,neq, the row indices i of the
!                        nonzero locations in column j are
!                        i = jan(k),  ian(j) .le. k .lt. ian(j+1).
!                        note that ian(neq+1) = nnz + 1.
!                        (if moss = 0, ian/jan may differ from the
!                        input ia/ja because of a different ordering
!                        in each column, and added diagonal entries.)
!
! yh      lyh            the nordsieck history array, of size nyh by
!          (optional     (nqcur + 1), where nyh is the initial value
!          output)       of neq.  for j = 0,1,...,nqcur, column j+1
!                        of yh contains hcur**j/factorial(j) times
!                        the j-th derivative of the interpolating
!                        polynomial currently representing the solution,
!                        evaluated at t = tcur.  the base address lyh
!                        is another optional output, listed above.
!
! acor     lenrw-neq+1   array of size neq used for the accumulated
!                        corrections on each step, scaled on output
!                        to represent the estimated local error in y
!                        on the last step.  this is the vector e in
!                        the description of the error control.  it is
!                        defined only on a successful return from
!                        lsodes.
!
!-----------------------------------------------------------------------
! part ii.  other routines callable.
!
! the following are optional calls which the user may make to
! gain additional capabilities in conjunction with lsodes.
! (the routines xsetun and xsetf are designed to conform to the
! slatec error handling package.)
!
!     form of call                  function
!   call xsetun(lun)          set the logical unit number, lun, for
!                             output of messages from lsodes, if
!                             the default is not desired.
!                             the default value of lun is 6.
!
!   call xsetf(mflag)         set a flag to control the printing of
!                             messages by lsodes.
!                             mflag = 0 means do not print. (danger..
!                             this risks losing valuable information.)
!                             mflag = 1 means print (the default).
!
!                             either of the above calls may be made at
!                             any time and will take effect immediately.
!
!   call srcms(rsav,isav,job) saves and restores the contents of
!                             the internal common blocks used by
!                             lsodes (see part iii below).
!                             rsav must be a real array of length 224
!                             or more, and isav must be an integer
!                             array of length 75 or more.
!                             job=1 means save common into rsav/isav.
!                             job=2 means restore common from rsav/isav.
!                                srcms is useful if one is
!                             interrupting a run and restarting
!                             later, or alternating between two or
!                             more problems solved with lsodes.
!
!   call intdy(,,,,,)         provide derivatives of y, of various
!        (see below)          orders, at a specified point t, if
!                             desired.  it may be called only after
!                             a successful return from lsodes.
!
! the detailed instructions for using intdy are as follows.
! the form of the call is..
!
!   lyh = iwork(22)
!   call intdy (t, k, rwork(lyh), nyh, dky, iflag)
!
! the input parameters are..
!
! t         = value of independent variable where answers are desired
!             (normally the same as the t last returned by lsodes).
!             for valid results, t must lie between tcur - hu and tcur.
!             (see optional outputs for tcur and hu.)
! k         = integer order of the derivative desired.  k must satisfy
!             0 .le. k .le. nqcur, where nqcur is the current order
!             (see optional outputs).  the capability corresponding
!             to k = 0, i.e. computing y(t), is already provided
!             by lsodes directly.  since nqcur .ge. 1, the first
!             derivative dy/dt is always available with intdy.
! lyh       = the base address of the history array yh, obtained
!             as an optional output as shown above.
! nyh       = column length of yh, equal to the initial value of neq.
!
! the output parameters are..
!
! dky       = a real array of length neq containing the computed value
!             of the k-th derivative of y(t).
! iflag     = integer flag, returned as 0 if k and t were legal,
!             -1 if k was illegal, and -2 if t was illegal.
!             on an error return, a message is also written.
!-----------------------------------------------------------------------
! part iii.  common blocks.
!
! if lsodes is to be used in an overlay situation, the user
! must declare, in the primary overlay, the variables in..
!   (1) the call sequence to lsodes,
!   (2) the three internal common blocks
!         /ls0001/  of length  257  (218 double precision words
!                         followed by 39 integer words),
!         /lss001/  of length  40    ( 6 double precision words
!                         followed by 34 integer words),
!         /eh0001/  of length  2 (integer words).
!
! if lsodes is used on a system in which the contents of internal
! common blocks are not preserved between calls, the user should
! declare the above three common blocks in his main program to insure
! that their contents are preserved.
!
! if the solution of a given problem by lsodes is to be interrupted
! and then later continued, such as when restarting an interrupted run
! or alternating between two or more problems, the user should save,
! following the return from the last lsodes call prior to the
! interruption, the contents of the call sequence variables and the
! internal common blocks, and later restore these values before the
! next lsodes call for that problem.  to save and restore the common
! blocks, use subroutine srcms (see part ii above).
!
!-----------------------------------------------------------------------
! part iv.  optionally replaceable solver routines.
!
! below are descriptions of two routines in the lsodes package which
! relate to the measurement of errors.  either routine can be
! replaced by a user-supplied version, if desired.  however, since such
! a replacement may have a major impact on performance, it should be
! done only when absolutely necessary, and only with great caution.
! (note.. the means by which the package version of a routine is
! superseded by the user-s version may be system-dependent.)
!
! (a) ewset.
! the following subroutine is called just before each internal
! integration step, and sets the array of error weights, ewt, as
! described under itol/rtol/atol above..
!     subroutine ewset (neq, itol, rtol, atol, ycur, ewt)
! where neq, itol, rtol, and atol are as in the lsodes call sequence,
! ycur contains the current dependent variable vector, and
! ewt is the array of weights set by ewset.
!
! if the user supplies this subroutine, it must return in ewt(i)
! (i = 1,...,neq) a positive quantity suitable for comparing errors
! in y(i) to.  the ewt array returned by ewset is passed to the
! vnorm routine (see below), and also used by lsodes in the computation
! of the optional output imxer, the diagonal jacobian approximation,
! and the increments for difference quotient jacobians.
!
! in the user-supplied version of ewset, it may be desirable to use
! the current values of derivatives of y.  derivatives up to order nq
! are available from the history array yh, described above under
! optional outputs.  in ewset, yh is identical to the ycur array,
! extended to nq + 1 columns with a column length of nyh and scale
! factors of h**j/factorial(j).  on the first call for the problem,
! given by nst = 0, nq is 1 and h is temporarily set to 1.0.
! the quantities nq, nyh, h, and nst can be obtained by including
! in ewset the statements..
!     double precision h, rls
!     common /ls0001/ rls(218),ils(39)
!     nq = ils(35)
!     nyh = ils(14)
!     nst = ils(36)
!     h = rls(212)
! thus, for example, the current value of dy/dt can be obtained as
! ycur(nyh+i)/h  (i=1,...,neq)  (and the division by h is
! unnecessary when nst = 0).
!
! (b) vnorm.
! the following is a real function routine which computes the weighted
! root-mean-square norm of a vector v..
!     d = vnorm (n, v, w)
! where..
!   n = the length of the vector,
!   v = real array of length n containing the vector,
!   w = real array of length n containing weights,
!   d = sqrt( (1/n) * sum(v(i)*w(i))**2 ).
! vnorm is called with n = neq and with w(i) = 1.0/ewt(i), where
! ewt is as set by subroutine ewset.
!
! if the user supplies this function, it should return a non-negative
! value of vnorm suitable for use in the error control in lsodes.
! none of the arguments should be altered by vnorm.
! for example, a user-supplied vnorm routine might..
!   -substitute a max-norm of (v(i)*w(i)) for the rms-norm, or
!   -ignore some components of v in the norm, with the effect of
!    suppressing the error control on those components of y.
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
! other routines in the lsodes package.
!
! in addition to subroutine lsodes, the lsodes package includes the
! following subroutines and function routines..
!  iprep    acts as an iterface between lsodes and prep, and also does
!           adjusting of work space pointers and work arrays.
!  prep     is called by iprep to compute sparsity and do sparse matrix
!           preprocessing if miter = 1 or 2.
!  jgroup   is called by prep to compute groups of jacobian column
!           indices for use when miter = 2.
!  adjlr    adjusts the length of required sparse matrix work space.
!           it is called by prep.
!  cntnzu   is called by prep and counts the nonzero elements in the
!           strict upper triangle of j + j-transpose, where j = df/dy.
!  intdy    computes an interpolated value of the y vector at t = tout.
!  stode    is the core integrator, which does one step of the
!           integration and the associated error control.
!  cfode    sets all method coefficients and test constants.
!  prjs     computes and preprocesses the jacobian matrix j = df/dy
!           and the newton iteration matrix p = i - h*l0*j.
!  slss     manages solution of linear system in chord iteration.
!  ewset    sets the error weight vector ewt before each step.
!  vnorm    computes the weighted r.m.s. norm of a vector.
!  srcms    is a user-callable routine to save and restore
!           the contents of the internal common blocks.
!  odrv     constructs a reordering of the rows and columns of
!           a matrix by the minimum degree algorithm.  odrv is a
!           driver routine which calls subroutines md, mdi, mdm,
!           mdp, mdu, and sro.  see ref. 2 for details.  (the odrv
!           module has been modified since ref. 2, however.)
!  cdrv     performs reordering, symbolic factorization, numerical
!           factorization, or linear system solution operations,
!           depending on a path argument ipath.  cdrv is a
!           driver routine which calls subroutines nroc, nsfc,
!           nnfc, nnsc, and nntc.  see ref. 3 for details.
!           lsodes uses cdrv to solve linear systems in which the
!           coefficient matrix is  p = i - con*j, where i is the
!           identity, con is a scalar, and j is an approximation to
!           the jacobian df/dy.  because cdrv deals with rowwise
!           sparsity descriptions, cdrv works with p-transpose, not p.
!  d1mach   computes the unit roundoff in a machine-independent manner.
!  xerrwv, xsetun, and xsetf   handle the printing of all error
!           messages and warnings.  xerrwv is machine-dependent.
! note..  vnorm and d1mach are function routines.
! all the others are subroutines.
!
! the intrinsic and external routines used by lsodes are..
! dabs, dmax1, dmin1, dfloat, max0, min0, mod, dsign, dsqrt, and write.
!
! a block data subprogram is also included with the package,
! for loading some of the variables in internal common.
!
!-----------------------------------------------------------------------
! the following card is for optimized compilation on lll compilers.
!lll. optimize
!-----------------------------------------------------------------------
      subroutine lsodes(neq, y, t, tout, itol, rtol, atol, itask, istate, iopt, rwork, lrw, iwork, liw, mf)

      integer neq, itol, itask, istate, iopt, lrw, iwork, liw, mf
      double precision y, t, tout, rtol, atol, rwork
      dimension neq(1), y(1), rtol(1), atol(1), rwork(lrw), iwork(liw)
      integer i, i1, i2, iflag, imax, imul, imxer, ipflag, ipgo, irem, j, kgo, 
     +        lenyht, leniw, lenrw, lf0, lia, lja, lrtem, lwtem, lyhd, lyhn, mf1, mord, mxhnl0, ncolm
      double precision atoli, ayi, big, ewti, h0, hmax, hmx, rh, rtoli, tcrit, tdist, tnext, tol, tolsf, tp, size, sum, w0, vnorm
      dimension mord(2)
      logical ihit

!     the following two internal common blocks contain
!     (a) variables which are local to any subroutine but whose values must be preserved between calls to the routine (own variables), and
!     (b) variables which are communicated between subroutines. the structure of each block is as follows..  
!         all real variables are listed first, followed by all integers.  
!         within each type, the variables are grouped with those local to subroutine lsodes first,
!         then those local to subroutine stode or subroutine prjs (no other routines have own variables), and finally those used for communication.
!     the block ls0001 is declared in subroutines lsodes, iprep, prep, intdy, stode, prjs, and slss.  
!     the block lss001 is declared in subroutines lsodes, iprep, prep, prjs, and slss.
!     groups of variables are replaced by dummy arrays in the common declarations in routines where those variables are not used.

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

      data mord(1),mord(2)/12,5/, mxhnl0/10/

!     block a.
!     this code block is executed on every call. it tests istate and itask for legality and branches appropriately.
!     if istate .gt. 1 but the flag init shows that initialization has not yet been done, an error return occurs.
!     if istate = 1 and tout = t, jump to block g and return immediately.
      if(istate .lt. 1 .or. istate .gt. 3)then
         call xerrwv("lsodes-- istate (=i1) illegal ", 30, 1, 0, 1, istate, 0, 0, 0.0d0, 0.0d0)
         STOP
      end if
      if(itask .lt. 1 .or. itask .gt. 5)then
         call xerrwv("lsodes-- itask (=i1) illegal  ", 30, 2, 0, 1, itask, 0, 0, 0.0d0, 0.0d0)
         STOP
      end if
      if(istate .ne. 1 .and. init .eq. 0)then
         call xerrwv("lsodes-- istate .gt. 1 but lsodes not initialized ", 50, 3, 0, 0, 0, 0, 0, 0.0d0, 0.0d0)
         STOP
      end if
      if(istate .eq. 1)then
         init = 0
         if(tout .eq. t)then
            ntrep = ntrep + 1
            if(ntrep .lt. 5) return
            call xerrwv("lsodes-- repeated calls with istate = 1 and tout = t (=r1)  ", 60, 301, 0, 0, 0, 0, 1, t, 0.0d0)
            call xerrwv("lsodes-- run aborted.. apparent infinite loop     ", 50, 303, 2, 0, 0, 0, 0, 0.0d0, 0.0d0)
!           stop
         end if
      end if

      if(istate .eq. 1 .or. istate .eq. 3)then
         ntrep = 0
         
!        block b.
!        the next code block is executed for the initial call (istate = 1), or for a continuation call with parameter changes (istate = 3).
!        it contains checking of all inputs and various initializations. 
!        if istate = 1, the final setting of work space pointers, the matrix preprocessing, and other initializations are done in block c.
!        first check legality of the non-optional inputs neq, itol, iopt, mf, ml, and mu.
         if(neq(1) .le. 0)then
            call xerrwv("lsodes-- neq (=i1) .lt. 1     ", 30, 4, 0, 1, neq(1), 0, 0, 0.0d0, 0.0d0)
            STOP
         end if
         if(istate .ne. 1 .and. neq(1) .gt. n)then
            call xerrwv("lsodes-- istate = 3 and neq increased (i1 to i2)  ", 50, 5, 0, 2, n, neq(1), 0, 0.0d0, 0.0d0)
            STOP
         end if
         n = neq(1)
         if(itol .lt. 1 .or. itol .gt. 4)then
            call xerrwv("lsodes-- itol (=i1) illegal   ", 30, 6, 0, 1, itol, 0, 0, 0.0d0, 0.0d0)
            STOP
         end if
         if(iopt .lt. 0 .or. iopt .gt. 1)then
            call xerrwv("lsodes-- iopt (=i1) illegal   ", 30, 7, 0, 1, iopt, 0, 0, 0.0d0, 0.0d0)
            STOP
         end if
         moss = mf/100
         mf1 = mf - 100*moss
         meth = mf1/10
         miter = mf1 - 10*meth
         if(moss .lt. 0 .or. moss .gt. 2)then
            call xerrwv("lsodes-- mf (=i1) illegal     ", 30, 8, 0, 1, mf, 0, 0, 0.0d0, 0.0d0)
            STOP
         end if
         if(meth .lt. 1 .or. meth .gt. 2)then
            call xerrwv("lsodes-- mf (=i1) illegal     ", 30, 8, 0, 1, mf, 0, 0, 0.0d0, 0.0d0)
            STOP
         end if
         if(miter .lt. 0 .or. miter .gt. 3)then
            call xerrwv("lsodes-- mf (=i1) illegal     ", 30, 8, 0, 1, mf, 0, 0, 0.0d0, 0.0d0)
            STOP
         end if
         if(miter .eq. 0 .or. miter .eq. 3) moss = 0
!        next process and check the optional inputs.
         if(iopt .eq. 0)then
            maxord = mord(meth)
            mxstep = 1000000
            mxhnil = mxhnl0
            if(istate .eq. 1) h0 = 0.0d0
            hmxi = 0.0d0
            hmin = 0.0d0
            seth = 0.0d0
         else ! iopt == 1
            maxord = iwork(5)
            if(maxord .lt. 0)then
               call xerrwv("lsodes-- maxord (=i1) .lt. 0  ", 30, 11, 0, 1, maxord, 0, 0, 0.0d0, 0.0d0)
               STOP
            end if
            if(maxord .eq. 0) maxord = 100
            maxord = min0(maxord,mord(meth))
            mxstep = iwork(6)
            if(mxstep .lt. 0)then
               call xerrwv("lsodes-- mxstep (=i1) .lt. 0  ", 30, 12, 0, 1, mxstep, 0, 0, 0.0d0, 0.0d0)
               STOP
            end if
            if(mxstep .eq. 0) mxstep = 500
            mxhnil = iwork(7)
            if(mxhnil .lt. 0)then
               call xerrwv("lsodes-- mxhnil (=i1) .lt. 0  ", 30, 13, 0, 1, mxhnil, 0, 0, 0.0d0, 0.0d0)
               STOP
            end if
            if(mxhnil .eq. 0) mxhnil = mxhnl0
            if(istate .eq. 1)then
               h0 = rwork(5)
               if((tout - t)*h0 .lt. 0.0d0)then
                  call xerrwv("lsodes-- tout (=r1) behind t (=r2)      ", 40, 14, 0, 0, 0, 0, 2, tout, t)
                  call xerrwv("      integration direction is given by h0 (=r1)  ", 50, 14, 0, 0, 0, 0, 1, h0, 0.0d0)
                  STOP
               end if
            end if
            hmax = rwork(6)
            if(hmax .lt. 0.0d0)then
               call xerrwv("lsodes-- hmax (=r1) .lt. 0.0  ", 30, 15, 0, 0, 0, 0, 1, hmax, 0.0d0)
               STOP
            end if
            hmxi = 0.0d0
            if(hmax .gt. 0.0d0) hmxi = 1.0d0/hmax
            hmin = rwork(7)
            if(hmin .lt. 0.0d0)then
               call xerrwv("lsodes-- hmin (=r1) .lt. 0.0  ", 30, 16, 0, 0, 0, 0, 1, hmin, 0.0d0)
               STOP
            end if
            seth = rwork(8)
            if(seth .lt. 0.0d0)then
               call xerrwv("lsodes-- seth (=r1) .lt. 0.0  ", 30, 9, 0, 0, 0, 0, 1, seth, 0.0d0)
               STOP
            end if
         end if
         
!        check rtol and atol for legality.
         rtoli = rtol(1)
         atoli = atol(1)
         do i = 1,n
            if(itol .ge. 3) rtoli = rtol(i)
            if(itol .eq. 2 .or. itol .eq. 4) atoli = atol(i)
            if(rtoli .lt. 0.0d0)then
               call xerrwv("lsodes-- rtol(i1) is r1 .lt. 0.0        ", 40, 19, 0, 1, i, 0, 1, rtoli, 0.0d0)
               STOP
            end if
            if(atoli .lt. 0.0d0)then
               call xerrwv("lsodes-- atol(i1) is r1 .lt. 0.0        ", 40, 20, 0, 1, i, 0, 1, atoli, 0.0d0)
               STOP
            end if
         end do
         
!        compute required work array lengths, as far as possible, and test these against lrw and liw.  then set tentative pointers for work arrays.
!        pointers to rwork/iwork segments are named by prefixing l to the name of the segment.  e.g., the segment yh starts at rwork(lyh).
!        segments of rwork (in order) are denoted  wm, yh, savf, ewt, acor.
!        if miter = 1 or 2, the required length of the matrix work space wm is not yet known, 
!        and so a crude minimum value is used for the initial tests of lrw and liw, and yh is temporarily stored as far
!        to the right in rwork as possible, to leave the maximum amount of space for wm for matrix preprocessing.
!        thus if miter = 1 or 2 and moss .ne. 2, some of the segments of rwork are temporarily omitted, as they are not needed in the preprocessing.
!        these omitted segments are.. acor if istate = 1, ewt and acor if istate = 3 and moss = 1, and savf, ewt, and acor if istate = 3 and moss = 0.
         if(istate .eq. 1) nyh = n
         lwmin = 0
         if(miter .eq. 1) lwmin = 4*n + 10*n/2
         if(miter .eq. 2) lwmin = 4*n + 11*n/2
         if(miter .eq. 3) lwmin = n + 2
         lenyh = (maxord+1)*nyh
         lrest = lenyh + 3*n
         lenrw = 20 + lwmin + lrest
         iwork(17) = lenrw
         leniw = 30
         if(moss .eq. 0 .and. miter .ne. 0 .and. miter .ne. 3) leniw = leniw + n + 1
         iwork(18) = leniw
         if(lenrw .gt. lrw)then
            call xerrwv("lsodes-- rwork length is insufficient to proceed. ", 50, 17, 0, 0, 0, 0, 0, 0.0d0, 0.0d0)
            call xerrwv("        length needed is .ge. lenrw (=i1), exceeds lrw (=i2)", 60, 17, 0, 2, lenrw, lrw, 0, 0.0d0, 0.0d0)
            STOP
         end if
         if(leniw .gt. liw)then
            call xerrwv("lsodes-- iwork length is insufficient to proceed. ", 50, 18, 0, 0, 0, 0, 0, 0.0d0, 0.0d0)
            call xerrwv("        length needed is .ge. leniw (=i1), exceeds liw (=i2)", 60, 18, 0, 2, leniw, liw, 0, 0.0d0, 0.0d0)
            STOP
         end if
         lia = 31
         if(moss .eq. 0 .and. miter .ne. 0 .and. miter .ne. 3) leniw = leniw + iwork(lia+n) - 1
         iwork(18) = leniw
         if(leniw .gt. liw)then
            call xerrwv("lsodes-- iwork length is insufficient to proceed. ", 50, 18, 0, 0, 0, 0, 0, 0.0d0, 0.0d0)
            call xerrwv("        length needed is .ge. leniw (=i1), exceeds liw (=i2)", 60, 18, 0, 2, leniw, liw, 0, 0.0d0, 0.0d0)
            STOP
         end if
         lja = lia + n + 1
         lia = min0(lia,liw)
         lja = min0(lja,liw)
         lwm = 21
         if(istate .eq. 1) nq = 1
         ncolm = min0(nq+1,maxord+2)
         lenyhm = ncolm*nyh
         lenyht = lenyh
         if(miter .eq. 1 .or. miter .eq. 2) lenyht = lenyhm
         imul = 2
         if(istate .eq. 3) imul = moss
         if(moss .eq. 2) imul = 3
         lrtem = lenyht + imul*n
         lwtem = lwmin
         if(miter .eq. 1 .or. miter .eq. 2) lwtem = lrw - 20 - lrtem
         lenwk = lwtem
         lyhn = lwm + lwtem
         lsavf = lyhn + lenyht
         lewt = lsavf + n
         lacor = lewt + n
         istatc = istate
         if(istate .eq. 1)then
!           block c.
!           the next block is for the initial call only (istate = 1).
!           it contains all remaining initializations, the initial call to f,
!           the sparse matrix preprocessing (miter = 1 or 2), and the
!           calculation of the initial step size.
!           the error weights in ewt are inverted after being loaded.
            lyh = lyhn
            iwork(22) = lyh
            tn = t
            nst = 0
            h = 1.0d0
            nnz = 0
            ngp = 0
            nzl = 0
            nzu = 0
!           load the initial value vector in yh.
            do i = 1,n
               rwork(i+lyh-1) = y(i)
            end do
!           initial call to f.  (lf0 points to yh(*,2).)
            lf0 = lyh + nyh
            call rhs(neq, t, y, rwork(lf0))
            nfe = 1
!           load and invert the ewt array.  (h is temporarily set to 1.0.)
            call ewset(n, itol, rtol, atol, rwork(lyh), rwork(lewt))
            do i = 1,n
               if(rwork(i+lewt-1) .le. 0.0d0)then
                  ewti = rwork(lewt+i-1)
                  call xerrwv("lsodes-- ewt(i1) is r1 .le. 0.0         ", 40, 21, 0, 1, i, 0, 1, ewti, 0.0d0)
                  STOP
               end if
               rwork(i+lewt-1) = 1.0d0/rwork(i+lewt-1)
            end do
            if(miter .eq. 1 .or. miter .eq. 2)then
!              iprep and prep do sparse matrix preprocessing if miter = 1 or 2.
               lacor = min0(lacor,lrw)
               call iprep(neq, y, rwork, iwork(lia), iwork(lja), ipflag)
               lenrw = lwm - 1 + lenwk + lrest
               iwork(17) = lenrw
               if(ipflag .ne. -1) iwork(23) = ipian
               if(ipflag .ne. -1) iwork(24) = ipjan
               ipgo = -ipflag + 1
               if(ipgo .eq. 2)then
                  call xerrwv("lsodes-- rwork length insufficient (for subroutine prep).   ", 60, 28, 0, 0, 0, 0, 0, 0.0d0, 0.0d0)
                  call xerrwv("        length needed is .ge. lenrw (=i1), exceeds lrw (=i2)", 60, 28, 0, 2, lenrw, lrw, 0, 0.0d0, 0.0d0)
                  STOP
               else if(ipgo .eq. 3)then
                  call xerrwv("lsodes-- rwork length insufficient (for subroutine jgroup). ", 60, 29, 0, 0, 0, 0, 0, 0.0d0, 0.0d0)
                  call xerrwv("        length needed is .ge. lenrw (=i1), exceeds lrw (=i2)", 60, 29, 0, 2, lenrw, lrw, 0, 0.0d0, 0.0d0)
                  STOP
               else if(ipgo .eq. 4)then
                  call xerrwv("lsodes-- rwork length insufficient (for subroutine odrv).   ", 60, 30, 0, 0, 0, 0, 0, 0.0d0, 0.0d0)
                  call xerrwv("        length needed is .ge. lenrw (=i1), exceeds lrw (=i2)", 60, 30, 0, 2, lenrw, lrw, 0, 0.0d0, 0.0d0)
                  STOP
               else if(ipgo .eq. 5)then
                  call xerrwv("lsodes-- error from odrv in yale sparse matrix package      ", 60, 31, 0, 0, 0, 0, 0, 0.0d0, 0.0d0)
                  imul = (iys - 1)/n
                  irem = iys - imul*n
                  call xerrwv("      at t (=r1), odrv returned error flag = i1*neq + i2.   ", 60, 31, 0, 2, imul, irem, 1, tn, 0.0d0)
                  STOP
               else if(ipgo .eq. 6)then
                  call xerrwv("lsodes-- rwork length insufficient (for subroutine cdrv).   ", 60, 32, 0, 0, 0, 0, 0, 0.0d0, 0.0d0)
                  call xerrwv("        length needed is .ge. lenrw (=i1), exceeds lrw (=i2)", 60, 32, 0, 2, lenrw, lrw, 0, 0.0d0, 0.0d0)
                  STOP
               else if(ipgo .eq. 7)then
                  call xerrwv("lsodes-- error from cdrv in yale sparse matrix package      ", 60, 33, 0, 0, 0, 0, 0, 0.0d0, 0.0d0)
                  imul = (iys - 1)/n
                  irem = iys - imul*n
                  call xerrwv("      at t (=r1), cdrv returned error flag = i1*neq + i2.   ", 60, 33, 0, 2, imul, irem, 1, tn, 0.0d0)
                  if(imul .eq. 2) call xerrwv("        duplicate entry in sparsity structure descriptors   ", 60, 33, 0, 0, 0, 0, 0, 0.0d0, 0.0d0)
                  if(imul .eq. 3 .or. imul .eq. 6) call xerrwv("        insufficient storage for nsfc (called by cdrv)      ", 60, 33, 0, 0, 0, 0, 0, 0.0d0, 0.0d0)
                  STOP
               end if

               iwork(22) = lyh
               if(lenrw .gt. lrw)then
                  call xerrwv("lsodes-- rwork length is insufficient to proceed. ", 50, 17, 0, 0, 0, 0, 0, 0.0d0, 0.0d0)
                  call xerrwv("        length needed is .ge. lenrw (=i1), exceeds lrw (=i2)", 60, 17, 0, 2, lenrw, lrw, 0, 0.0d0, 0.0d0)
                  STOP
               end if
!              check tcrit for legality (itask = 4 or 5).
            end if
            if(itask .eq. 4 .or. itask .eq. 5)then
               tcrit = rwork(1)
               if((tcrit - tout)*(tout - t) .lt. 0.0d0) then
                  call xerrwv("lsodes-- itask = 4 or 5 and tcrit (=r1) behind tout (=r2)   ", 60, 25, 0, 0, 0, 0, 2, tcrit, tout)
                  STOP
               end if
               if(h0 .ne. 0.0d0 .and. (t + h0 - tcrit)*h0 .gt. 0.0d0) h0 = tcrit - t
            end if
!           initialize all remaining parameters.
            uround = 2.2204460492503131E-016
            jstart = 0
            if(miter .ne. 0) rwork(lwm) = dsqrt(uround)
            msbj = 50
            nslj = 0
            ccmxj = 0.2d0
            psmall = 1000.0d0*uround
            rbig = 0.01d0/psmall
            nhnil = 0
            nje = 0
            nlu = 0
            nslast = 0
            hu = 0.0d0
            nqu = 0
            ccmax = 0.3d0
            maxcor = 3
            msbp = 20
            mxncf = 10
!           the coding below computes the step size, h0, to be attempted on the first step, unless the user has supplied a value for this.
!           first check that tout - t differs significantly from zero.
!           a scalar tolerance quantity tol is computed, as max(rtol(i)) if this is positive, or max(atol(i)/abs(y(i))) otherwise, 
!           adjusted so as to be between 100*uround and 1.0e-3.
!           then the computed value h0 is given by..
!                                                neq
!             h0**2 = tol / ( w0**-2 + (1/neq) * sum ( f(i)/ywt(i) )**2  )
!                                                 1
!           where   w0     = max ( abs(t), abs(tout) ),
!                   f(i)   = i-th component of initial value of f,
!                   ywt(i) = ewt(i)/tol  (a weight for y(i)).
!           the sign of h0 is inferred from the initial values of tout and t.
            lf0 = lyh + nyh
            if(h0 .eq. 0.0d0)then
               tdist = dabs(tout - t)
               w0 = dmax1(dabs(t),dabs(tout))
               if(tdist .lt. 2.0d0*uround*w0)then
                  call xerrwv("lsodes-- tout (=r1) too close to t(=r2) to start integration", 60, 22, 0, 0, 0, 0, 2, tout, t)
                  STOP
               end if
               tol = rtol(1)
               if(itol .gt. 2)then
                  do i = 1,n
                     tol = dmax1(tol,rtol(i))
                  end do
               end if
               if(tol .le. 0.0d0)then
                  atoli = atol(1)
                  do i = 1,n
                     if(itol .eq. 2 .or. itol .eq. 4) atoli = atol(i)
                     ayi = dabs(y(i))
                     if(ayi .ne. 0.0d0) tol = dmax1(tol,atoli/ayi)
                  end do
               end if
               tol = dmax1(tol,100.0d0*uround)
               tol = dmin1(tol,0.001d0)
               sum = vnorm (n, rwork(lf0), rwork(lewt))
               sum = 1.0d0/(tol*w0*w0) + tol*sum**2
               h0 = 1.0d0/dsqrt(sum)
               h0 = dmin1(h0,tdist)
               h0 = dsign(h0,tout-t)
            end if
!           adjust h0 if necessary to meet hmax bound.
            rh = dabs(h0)*hmxi
            if(rh .gt. 1.0d0) h0 = h0/rh
!           load h with h0 and scale yh(*,2) by h0.
            h = h0
            do i = 1,n
               rwork(i+lf0-1) = h0*rwork(i+lf0-1)
            end do

         else ! istate .eq. 3
         
!           istate = 3. move yh to its new location.
!           note that only the part of yh needed for the next step, namely min(nq+1,maxord+2) columns, is actually moved.
!           a temporary error weight array ewt is loaded if moss = 2. sparse matrix processing is done in iprep/prep if miter = 1 or 2.
!           if maxord was reduced below nq, then the pointers are finally set so that savf is identical to yh(*,maxord+2).
            lyhd = lyh - lyhn
            imax = lyhn - 1 + lenyhm
!           move yh. branch for move right, no move, or move left.
            if(lyhd .eq. 1)then
!              move right
               do i = lyhn,imax
                  j = imax + lyhn - i
                  rwork(j) = rwork(j+lyhd)
               end do
            else if(lyhd .eq. 3)then
!              move left
               do i = lyhn,imax
                  rwork(i) = rwork(i+lyhd)
               end do
            end if
            lyh = lyhn
            iwork(22) = lyh
            if(miter .eq. 1 .or. miter .eq. 2)then
               if(moss .eq. 2)then
!                 temporarily load ewt if miter = 1 or 2 and moss = 2.
                  call ewset (n, itol, rtol, atol, rwork(lyh), rwork(lewt))
                  do i = 1,n
                     if(rwork(i+lewt-1) .le. 0.0d0)then
                        ewti = rwork(lewt+i-1)
                        call xerrwv("lsodes-- ewt(i1) is r1 .le. 0.0         ", 40, 21, 0, 1, i, 0, 1, ewti, 0.0d0)
                        STOP
                     end if
                     rwork(i+lewt-1) = 1.0d0/rwork(i+lewt-1)
                  end do
               end if
!              iprep and prep do sparse matrix preprocessing if miter = 1 or 2.
               lsavf = min0(lsavf,lrw)
               lewt = min0(lewt,lrw)
               lacor = min0(lacor,lrw)
               call iprep(neq, y, rwork, iwork(lia), iwork(lja), ipflag)
               lenrw = lwm - 1 + lenwk + lrest
               iwork(17) = lenrw
               if(ipflag .ne. -1) iwork(23) = ipian
               if(ipflag .ne. -1) iwork(24) = ipjan
               ipgo = -ipflag + 1

               if(ipgo .eq. 2)then
                  call xerrwv("lsodes-- rwork length insufficient (for subroutine prep).   ", 60, 28, 0, 0, 0, 0, 0, 0.0d0, 0.0d0)
                  call xerrwv("        length needed is .ge. lenrw (=i1), exceeds lrw (=i2)", 60, 28, 0, 2, lenrw, lrw, 0, 0.0d0, 0.0d0)
                  STOP
               else if(ipgo .eq. 3)then
                  call xerrwv("lsodes-- rwork length insufficient (for subroutine jgroup). ", 60, 29, 0, 0, 0, 0, 0, 0.0d0, 0.0d0)
                  call xerrwv("        length needed is .ge. lenrw (=i1), exceeds lrw (=i2)", 60, 29, 0, 2, lenrw, lrw, 0, 0.0d0, 0.0d0)
                  STOP
               else if(ipgo .eq. 4)then
                  call xerrwv("lsodes-- rwork length insufficient (for subroutine odrv).   ", 60, 30, 0, 0, 0, 0, 0, 0.0d0, 0.0d0)
                  call xerrwv("        length needed is .ge. lenrw (=i1), exceeds lrw (=i2)", 60, 30, 0, 2, lenrw, lrw, 0, 0.0d0, 0.0d0)
                  STOP
               else if(ipgo .eq. 5)then
                  call xerrwv("lsodes-- error from odrv in yale sparse matrix package      ", 60, 31, 0, 0, 0, 0, 0, 0.0d0, 0.0d0)
                  imul = (iys - 1)/n
                  irem = iys - imul*n
                  call xerrwv("      at t (=r1), odrv returned error flag = i1*neq + i2.   ", 60, 31, 0, 2, imul, irem, 1, tn, 0.0d0)
                  STOP
               else if(ipgo .eq. 6)then
                  call xerrwv("lsodes-- rwork length insufficient (for subroutine cdrv).   ", 60, 32, 0, 0, 0, 0, 0, 0.0d0, 0.0d0)
                  call xerrwv("        length needed is .ge. lenrw (=i1), exceeds lrw (=i2)", 60, 32, 0, 2, lenrw, lrw, 0, 0.0d0, 0.0d0)
                  STOP
               else if(ipgo .eq. 7)then
                  call xerrwv("lsodes-- error from cdrv in yale sparse matrix package      ", 60, 33, 0, 0, 0, 0, 0, 0.0d0, 0.0d0)
                  imul = (iys - 1)/n
                  irem = iys - imul*n
                  call xerrwv("      at t (=r1), cdrv returned error flag = i1*neq + i2.   ", 60, 33, 0, 2, imul, irem, 1, tn, 0.0d0)
                  if(imul .eq. 2) call xerrwv("        duplicate entry in sparsity structure descriptors   ", 60, 33, 0, 0, 0, 0, 0, 0.0d0, 0.0d0)
                  if(imul .eq. 3 .or. imul .eq. 6) call xerrwv("        insufficient storage for nsfc (called by cdrv)      ", 60, 33, 0, 0, 0, 0, 0, 0.0d0, 0.0d0)
                  STOP
               end if

               iwork(22) = lyh
               if(lenrw .gt. lrw)then
                  call xerrwv("lsodes-- rwork length is insufficient to proceed. ", 50, 17, 0, 0, 0, 0, 0, 0.0d0, 0.0d0)
                  call xerrwv("        length needed is .ge. lenrw (=i1), exceeds lrw (=i2)", 60, 17, 0, 2, lenrw, lrw, 0, 0.0d0, 0.0d0)
                  STOP
               end if
            end if

!           set flag to signal parameter changes to stode.
            jstart = -1
            if(n .ne. nyh)then
!              neq was reduced.  zero part of yh to avoid undefined references.
               i1 = lyh + (nq+1)*nyh
               i2 = lyh + (maxord + 1)*nyh - 1
               if(i1 .le. i2)then
                  do i = i1,i2
                     rwork(i) = 0.0d0
                  end do
               end if
            end if
         end if
      end if

      if(istate .eq. 2 .or. istate .eq. 3)then
!        block d.
!        the next code block is for continuation calls only (istate = 2 or 3)
!        and is to check stop conditions before taking a step.
         nslast = nst
         if(itask .eq. 1)then
            if((tn - tout)*h .ge. 0.0d0)then
               call intdy (tout, 0, rwork(lyh), nyh, y, iflag)
               if(iflag .ne. 0)then
                  call xerrwv("lsodes-- trouble from intdy. itask = i1, tout = r1", 50, 27, 0, 1, itask, 0, 1, tout, 0.0d0)
                  STOP
               end if
               t = tout
               istate = 2
               rwork(11:13) = (/ hu, h, tn /)
               iwork(11:15) = (/ nst, nfe, nje, nqu, nq /)
               iwork(19:21) = (/ nnz, ngp, nlu /)
               iwork(25:26) = (/ nzl, nzu /)
               RETURN
            end if
         
         else if(itask .eq. 3)then
            tp = tn - hu*(1.0d0 + 100.0d0*uround)
            if((tp - tout)*h .gt. 0.0d0)then
               call xerrwv("lsodes-- itask = i1 and tout (=r1) behind tcur - hu (= r2)  ", 60, 23, 0, 1, itask, 0, 2, tout, tp)
               STOP
            end if
            if((tn - tout)*h .ge. 0.0d0)then
               do i = 1,n
                  y(i) = rwork(i+lyh-1)
               end do
               t = tn
               istate = 2
               rwork(11:13) = (/ hu, h, tn /)
               iwork(11:15) = (/ nst, nfe, nje, nqu, nq /)
               iwork(19:21) = (/ nnz, ngp, nlu /)
               iwork(25:26) = (/ nzl, nzu /)
               RETURN
            end if
         
         else if(itask .eq. 4)then
            tcrit = rwork(1)
            if((tn - tcrit)*h .gt. 0.0d0)then
               call xerrwv("lsodes-- itask = 4 or 5 and tcrit (=r1) behind tcur (=r2)   ", 60, 24, 0, 0, 0, 0, 2, tcrit, tn)
               STOP
            end if
            if((tcrit - tout)*h .lt. 0.0d0)then
               call xerrwv("lsodes-- itask = 4 or 5 and tcrit (=r1) behind tout (=r2)   ", 60, 25, 0, 0, 0, 0, 2, tcrit, tout)
               STOP
            end if
            if((tn - tout)*h .ge. 0.0d0)then
               call intdy (tout, 0, rwork(lyh), nyh, y, iflag)
               if(iflag .ne. 0)then
                  call xerrwv("lsodes-- trouble from intdy. itask = i1, tout = r1", 50, 27, 0, 1, itask, 0, 1, tout, 0.0d0)
                  STOP
               end if
               t = tout
               istate = 2
               rwork(11:13) = (/ hu, h, tn /)
               iwork(11:15) = (/ nst, nfe, nje, nqu, nq /)
               iwork(19:21) = (/ nnz, ngp, nlu /)
               iwork(25:26) = (/ nzl, nzu /)
               RETURN
            end if
         
         else if(itask .eq. 5)then
            tcrit = rwork(1)
            if((tn - tcrit)*h .gt. 0.0d0)then
               call xerrwv("lsodes-- itask = 4 or 5 and tcrit (=r1) behind tcur (=r2)   ", 60, 24, 0, 0, 0, 0, 2, tcrit, tn)
               STOP
            end if
         end if

         if(itask .eq. 4 .or. itask .eq. 5)then
            hmx = dabs(tn) + dabs(h)
            ihit = dabs(tn - tcrit) .le. 100.0d0*uround*hmx
            if(ihit)then
               do i = 1,n
                  y(i) = rwork(i+lyh-1)
               end do
               t = tn
               if(ihit) t = tcrit
               istate = 2
               rwork(11:13) = (/ hu, h, tn /)
               iwork(11:15) = (/ nst, nfe, nje, nqu, nq /)
               iwork(19:21) = (/ nnz, ngp, nlu /)
               iwork(25:26) = (/ nzl, nzu /)
               RETURN
            end if
            tnext = tn + h*(1.0d0 + 4.0d0*uround)
            if((tnext - tcrit)*h .gt. 0.0d0)then
               h = (tcrit - tn)*(1.0d0 - 4.0d0*uround)
               if(istate .eq. 2) jstart = -2
            end if
         end if
      end if

!     block e.
!     the next block is normally executed for all calls and contains the call to the one-step core integrator stode.
!     this is a looping point for the integration steps.
!     first check for too many steps being taken, update ewt (if not at start of problem), check for too much accuracy being requested, and
!     check for h below the roundoff level in t.
      do while(.true.)
!         if(istate .ne. 1)then
            if((nst-nslast) .ge. mxstep)then
!              the maximum number of steps was taken before reaching tout.
               call xerrwv("lsodes-- at current t (=r1), mxstep (=i1) steps   ", 50, 201, 0, 0, 0, 0, 0, 0.0d0, 0.0d0)
               call xerrwv("      taken on this call before reaching tout     ", 50, 201, 0, 1, mxstep, 0, 1, tn, 0.0d0)
               istate = -1
!              set y vector, t and optional outputs.
               do i = 1,n
                  y(i) = rwork(i+lyh-1)
               end do
               t = tn
               rwork(11:13) = (/ hu, h, tn /)
               iwork(11:15) = (/ nst, nfe, nje, nqu, nq /)
               iwork(19:21) = (/ nnz, ngp, nlu /)
               iwork(25:26) = (/ nzl, nzu /)
               RETURN
            end if
            call ewset (n, itol, rtol, atol, rwork(lyh), rwork(lewt))
            do i = 1,n
               if(rwork(i+lewt-1) .le. 0.0d0)then
!                 ewt(i) .le. 0.0 for some i (not at start of problem).
                  ewti = rwork(lewt+i-1)
                  call xerrwv("lsodes-- at t (=r1), ewt(i1) has become r2 .le. 0.", 50, 202, 0, 1, i, 0, 2, tn, ewti)
                  istate = -6
!                 set y vector, t and optional outputs.
                  do j = 1,n
                     y(j) = rwork(j+lyh-1)
                  end do
                  t = tn
                  rwork(11:13) = (/ hu, h, tn /)
                  iwork(11:15) = (/ nst, nfe, nje, nqu, nq /)
                  iwork(19:21) = (/ nnz, ngp, nlu /)
                  iwork(25:26) = (/ nzl, nzu /)
                  RETURN
               end if
               rwork(i+lewt-1) = 1.0d0/rwork(i+lewt-1)
            end do
!         end if
         tolsf = uround*vnorm (n, rwork(lyh), rwork(lewt))
         if(tolsf .gt. 1.0d0)then
            tolsf = tolsf*2.0d0
            if(nst .eq. 0)then
               call xerrwv("lsodes-- at start of problem, too much accuracy   ", 50, 26, 0, 0, 0, 0, 0, 0.0d0, 0.0d0)
               call xerrwv("      requested for precision of machine..  see tolsf (=r1) ", 60, 26, 0, 0, 0, 0, 1, tolsf, 0.0d0)
               STOP
            end if
!           too much accuracy requested for machine precision.
            call xerrwv("lsodes-- at t (=r1), too much accuracy requested  ", 50, 203, 0, 0, 0, 0, 0, 0.0d0, 0.0d0)
            call xerrwv("      for precision of machine..  see tolsf (=r2) ", 50, 203, 0, 0, 0, 0, 2, tn, tolsf)
            rwork(14) = tolsf
            istate = -2
!           set y vector, t and optional outputs.
            do i = 1,n
               y(i) = rwork(i+lyh-1)
            end do
            t = tn
            rwork(11:13) = (/ hu, h, tn /)
            iwork(11:15) = (/ nst, nfe, nje, nqu, nq /)
            iwork(19:21) = (/ nnz, ngp, nlu /)
            iwork(25:26) = (/ nzl, nzu /)
            RETURN
         end if
         
         if((tn + h) .eq. tn)then
            nhnil = nhnil + 1
            if(nhnil .le. mxhnil)then
               call xerrwv("lsodes-- warning..internal t (=r1) and h (=r2) are", 50, 101, 0, 0, 0, 0, 0, 0.0d0, 0.0d0)
               call xerrwv("      such that in the machine, t + h = t on the next step  ", 60, 101, 0, 0, 0, 0, 0, 0.0d0, 0.0d0)
               call xerrwv("      (h = step size). solver will continue anyway", 50, 101, 0, 0, 0, 0, 2, tn, h)
               if(nhnil .eq. mxhnil)then
                  call xerrwv("lsodes-- above warning has been issued i1 times.  2", 50, 102, 0, 0, 0, 0, 0, 0.0d0, 0.0d0)
                  call xerrwv("      it will not be issued again for this problem", 50, 102, 0, 1, mxhnil, 0, 0, 0.0d0, 0.0d0)
               end if
            end if
         end if

         call stode (neq, y, rwork(lyh), nyh, rwork(lyh), rwork(lewt), rwork(lsavf), rwork(lacor), rwork(lwm), rwork(lwm))

         kgo = 1 - kflag
         
!        block f.
!        the following block handles the case of a successful return from the
!        core integrator (kflag = 0).  test for stop conditions.
         if(kgo .eq. 1)then
            init = 1

         else if(kgo .eq. 2)then
!           kflag = -1.  error test failed repeatedly or with abs(h) = hmin.
            call xerrwv("lsodes-- at t(=r1) and step size h(=r2), the error", 50, 204, 0, 0, 0, 0, 0, 0.0d0, 0.0d0)
            call xerrwv("      test failed repeatedly or with abs(h) = hmin", 50, 204, 0, 0, 0, 0, 2, tn, h)
            istate = -4
!           compute imxer.
            big = 0.0d0
            imxer = 1
            do i = 1,n
               size = dabs(rwork(i+lacor-1)*rwork(i+lewt-1))
               if(big .lt. size)then
                  big = size
                  imxer = i
               end if
            end do
            iwork(16) = imxer
!           set y vector, t and optional outputs.
            do i = 1,n
               y(i) = rwork(i+lyh-1)
            end do
            t = tn
            rwork(11:13) = (/ hu, h, tn /)
            iwork(11:15) = (/ nst, nfe, nje, nqu, nq /)
            iwork(19:21) = (/ nnz, ngp, nlu /)
            iwork(25:26) = (/ nzl, nzu /)
            RETURN

         else if(kgo .eq. 3)then
!           kflag = -2.  convergence failed repeatedly or with abs(h) = hmin.
            call xerrwv("lsodes-- at t (=r1) and step size h (=r2), the    ", 50, 205, 0, 0, 0, 0, 0, 0.0d0, 0.0d0)
            call xerrwv("      corrector convergence failed repeatedly     ", 50, 205, 0, 0, 0, 0, 0, 0.0d0, 0.0d0)
            call xerrwv("      or with abs(h) = hmin   ", 30, 205, 0, 0, 0, 0, 2, tn, h)
            istate = -5
!           compute imxer.
            big = 0.0d0
            imxer = 1
            do i = 1,n
               size = dabs(rwork(i+lacor-1)*rwork(i+lewt-1))
               if(big .lt. size)then
                  big = size
                  imxer = i
               end if
            end do
            iwork(16) = imxer
!           set y vector, t and optional outputs.
            do i = 1,n
               y(i) = rwork(i+lyh-1)
            end do
            t = tn
            rwork(11:13) = (/ hu, h, tn /)
            iwork(11:15) = (/ nst, nfe, nje, nqu, nq /)
            iwork(19:21) = (/ nnz, ngp, nlu /)
            iwork(25:26) = (/ nzl, nzu /)
            RETURN

         else if(kgo .eq. 4)then
!           kflag = -3.  fatal error flag returned by prjs or slss (cdrv).
            call xerrwv("lsodes-- at t (=r1) and step size h (=r2), a fatal", 50, 207, 0, 0, 0, 0, 0, 0.0d0, 0.0d0)
            call xerrwv("      error flag was returned by cdrv (by way of  ", 50, 207, 0, 0, 0, 0, 0, 0.0d0, 0.0d0)
            call xerrwv("      subroutine prjs or slss)", 30, 207, 0, 0, 0, 0, 2, tn, h)
            istate = -7
!           set y vector, t and optional outputs.
            do i = 1,n
               y(i) = rwork(i+lyh-1)
            end do
            t = tn
            rwork(11:13) = (/ hu, h, tn /)
            iwork(11:15) = (/ nst, nfe, nje, nqu, nq /)
            iwork(19:21) = (/ nnz, ngp, nlu /)
            iwork(25:26) = (/ nzl, nzu /)
            RETURN
         end if
         
         if(itask .eq. 1) then
!           if tout has been reached, interpolate.
            if((tn - tout)*h .ge. 0.0d0)then
               call intdy(tout, 0, rwork(lyh), nyh, y, iflag)
               t = tout
               istate = 2
               rwork(11:13) = (/ hu, h, tn /)
               iwork(11:15) = (/ nst, nfe, nje, nqu, nq /)
               iwork(19:21) = (/ nnz, ngp, nlu /)
               iwork(25:26) = (/ nzl, nzu /)
               RETURN
            end if
         
         else if(itask .eq. 2) then
!           jump to exit.
            do i = 1,n
               y(i) = rwork(i+lyh-1)
            end do
            t = tn
            istate = 2
            rwork(11:13) = (/ hu, h, tn /)
            iwork(11:15) = (/ nst, nfe, nje, nqu, nq /)
            iwork(19:21) = (/ nnz, ngp, nlu /)
            iwork(25:26) = (/ nzl, nzu /)
            RETURN
         
         else if(itask .eq. 3) then
!           jump to exit if tout was reached.
            if((tn - tout)*h .ge. 0.0d0)then
               do i = 1,n
                  y(i) = rwork(i+lyh-1)
               end do
               t = tn
               istate = 2
               rwork(11:13) = (/ hu, h, tn /)
               iwork(11:15) = (/ nst, nfe, nje, nqu, nq /)
               iwork(19:21) = (/ nnz, ngp, nlu /)
               iwork(25:26) = (/ nzl, nzu /)
               RETURN
            end if
         
         else if(itask .eq. 4) then
!           see if tout or tcrit was reached. adjust h if necessary.
            if((tn - tout)*h .ge. 0.0d0)then
               call intdy (tout, 0, rwork(lyh), nyh, y, iflag)
               t = tout
               istate = 2
               rwork(11:13) = (/ hu, h, tn /)
               iwork(11:15) = (/ nst, nfe, nje, nqu, nq /)
               iwork(19:21) = (/ nnz, ngp, nlu /)
               iwork(25:26) = (/ nzl, nzu /)
               RETURN
            else
               hmx = dabs(tn) + dabs(h)
               ihit = dabs(tn - tcrit) .le. 100.0d0*uround*hmx
               if(ihit)then
                  do i = 1,n
                     y(i) = rwork(i+lyh-1)
                  end do
                  t = tn
                  if(ihit) t = tcrit
                  istate = 2
                  rwork(11:13) = (/ hu, h, tn /)
                  iwork(11:15) = (/ nst, nfe, nje, nqu, nq /)
                  iwork(19:21) = (/ nnz, ngp, nlu /)
                  iwork(25:26) = (/ nzl, nzu /)
                  RETURN
               else
                  tnext = tn + h*(1.0d0 + 4.0d0*uround)
                  if((tnext - tcrit)*h .gt. 0.0d0)then
                     h = (tcrit - tn)*(1.0d0 - 4.0d0*uround)
                     jstart = -2
                  end if
               end if
            end if
         
         else if(itask .eq. 5) then
!           see if tcrit was reached and jump to exit.
            hmx = dabs(tn) + dabs(h)
            ihit = dabs(tn - tcrit) .le. 100.0d0*uround*hmx
            do i = 1,n
               y(i) = rwork(i+lyh-1)
            end do
            t = tn
            if(ihit) t = tcrit
            istate = 2
            rwork(11:13) = (/ hu, h, tn /)
            iwork(11:15) = (/ nst, nfe, nje, nqu, nq /)
            iwork(19:21) = (/ nnz, ngp, nlu /)
            iwork(25:26) = (/ nzl, nzu /)
            RETURN
         end if
      end do
      end
