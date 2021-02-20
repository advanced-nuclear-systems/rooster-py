#--------------------------------------------------------------------------------------------------
# this is the march 30, 1987 version of
# lsodes.. livermore solver for ordinary differential equations
#          with general sparse jacobian matrices.
# this version is in double precision.
#
# lsodes solves the initial value problem for stiff or nonstiff
# systems of first order ode-s,
#     dy/dt = f(t,y) , or, in component form,
#     dy(i)/dt = f(i) = f(i,t,y(1),y(2),...,y(neq)) (i = 1,...,neq).
# lsodes is a variant of the lsode package, and is intended for
# problems in which the jacobian matrix df/dy has an arbitrary
# sparse structure (when the problem is stiff).
#
# authors..      alan c. hindmarsh,
#                computing and mathematics research division, l-316
#                lawrence livermore national laboratory
#                livermore, ca 94550.
#
# and            andrew h. sherman
#                j. s. nolen and associates
#                houston, tx 77084
#-----------------------------------------------------------------------
# references..
# 1.  alan c. hindmarsh, odepack, a systematized collection of ode
#     solvers, in scientific computing, r. s. stepleman et al. (eds.),
#     north-holland, amsterdam, 1983, pp. 55-64.
#
# 2.  s. c. eisenstat, m. c. gursky, m. h. schultz, and a. h. sherman,
#     yale sparse matrix package.. i. the symmetric codes,
#     int. j. num. meth. eng., 18 (1982), pp. 1145-1151.
#
# 3.  s. c. eisenstat, m. c. gursky, m. h. schultz, and a. h. sherman,
#     yale sparse matrix package.. ii. the nonsymmetric codes,
#     research report no. 114, dept. of computer sciences, yale
#     university, 1977.
#-----------------------------------------------------------------------
# summary of usage.
#
# communication between the user and the lsodes package, for normal
# situations, is summarized here.  this summary describes only a subset
# of the full set of options available.  see the full description for
# details, including optional communication, nonstandard options,
# and instructions for special situations.  see also the example
# problem (with program and output) following this summary.
#
# a. first provide a subroutine of the form..
#               subroutine f (neq, t, y, ydot)
#               dimension y(neq), ydot(neq)
# which supplies the vector function f by loading ydot(i) with f(i).
#
# b. next determine (or guess) whether or not the problem is stiff.
# stiffness occurs when the jacobian matrix df/dy has an eigenvalue
# whose real part is negative and large in magnitude, compared to the
# reciprocal of the t span of interest.  if the problem is nonstiff,
# use a method flag mf = 10.  if it is stiff, there are two standard
# for the method flag, mf = 121 and mf = 222.  in both cases, lsodes
# requires the jacobian matrix in some form, and it treats this matrix
# in general sparse form, with sparsity structure determined internally.
# (for options where the user supplies the sparsity structure, see
# the full description of mf below.)
#
# c. if the problem is stiff, you are encouraged to supply the jacobian
# directly (mf = 121), but if this is not feasible, lsodes will
# compute it internally by difference quotients (mf = 222).
# if you are supplying the jacobian, provide a subroutine of the form..
#               subroutine jac (neq, t, y, j, ian, jan, pdj)
#               dimension y(1), ian(1), jan(1), pdj(1)
# here neq, t, y, and j are input arguments, and the jac routine is to
# load the array pdj (of length neq) with the j-th column of df/dy.
# i.e., load pdj(i) with df(i)/dy(j) for all relevant values of i.
# the arguments ian and jan should be ignored for normal situations.
# lsodes will call the jac routine with j = 1,2,...,neq.
# only nonzero elements need be loaded.  usually, a crude approximation
# to df/dy, possibly with fewer nonzero elements, will suffice.
#
# d. write a main program which calls subroutine lsodes once for
# each point at which answers are desired.  this should also provide
# for possible use of logical unit 6 for output of error messages
# by lsodes.  on the first call to lsodes, supply arguments as follows..
# f      = name of subroutine for right-hand side vector f.
#          this name must be declared external in calling program.
# neq    = number of first order ode-s.
# y      = array of initial values, of length neq.
# t      = the initial value of the independent variable.
# tout   = first point where output is desired (!= t).
# itol   = 1 or 2 according as atol (below) is a scalar or array.
# rtol   = relative tolerance parameter (scalar).
# atol   = absolute tolerance parameter (scalar or array).
#          the estimated local error in y(i) will be controlled so as
#          to be roughly less (in magnitude) than
#             ewt(i) = rtol*abs(y(i)) + atol     if itol = 1, or
#             ewt(i) = rtol*abs(y(i)) + atol(i)  if itol = 2.
#          thus the local error test passes if, in each component,
#          either the absolute error is less than atol (or atol(i)),
#          or the relative error is less than rtol.
#          use rtol = 0.0 for pure absolute error control, and
#          use atol = 0.0 (or atol(i) = 0.0) for pure relative error
#          control.  caution.. actual (global) errors may exceed these
#          local tolerances, so choose them conservatively.
# itask  = 1 for normal computation of output values of y at t = tout.
# istate = integer flag (input and output).  set istate = 1.
# iopt   = 0 to indicate no optional inputs used.
# rwork  = real work array of length at least..
#             20 + 16*neq            for mf = 10,
#             20 + (2 + 1./2)*nnz + (11 + 9./2)*neq
#                                    for mf = 121 or 222,
#          where..
#          nnz    = the number of nonzero elements in the sparse
#                   jacobian (if this is unknown, use an estimate), and
#          in any case, the required size of rwork cannot generally
#          be predicted in advance if mf = 121 or 222, and the value
#          above is a rough estimate of a crude lower bound.  some
#          experimentation with this size may be necessary.
#          (when known, the correct required length is an optional
#          output, available in iwork(17).)
# lrw    = declared length of rwork (in user-s dimension).
# iwork  = integer work array of length at least 30.
# liw    = declared length of iwork (in user-s dimension).
# jac    = name of subroutine for jacobian matrix (mf = 121).
#          if used, this name must be declared external in calling
#          program.  if not used, pass a dummy name.
# mf     = method flag.  standard values are..
#          10  for nonstiff (adams) method, no jacobian used.
#          121 for stiff (bdf) method, user-supplied sparse jacobian.
#          222 for stiff method, internally generated sparse jacobian.
# note that the main program must declare arrays y, rwork, iwork,
# and possibly atol.
#
# e. the output from the first call (or any call) is..
#      y = array of computed values of y(t) vector.
#      t = corresponding value of independent variable (normally tout).
# istate = 2  if lsodes was successful, negative otherwise.
#          -1 means excess work done on this call (perhaps wrong mf).
#          -2 means excess accuracy requested (tolerances too small).
#          -3 means illegal input detected (see printed message).
#          -4 means repeated error test failures (check all inputs).
#          -5 means repeated convergence failures (perhaps bad jacobian
#             supplied or wrong choice of mf or tolerances).
#          -6 means error weight became zero during problem. (solution
#             component i vanished, and atol or atol(i) = 0.)
#          -7 means a fatal error return flag came from the sparse
#             solver cdrv by way of prjs or slss.  should never happen.
#          a return with istate = -1, -4, or -5 may result from using
#          an inappropriate sparsity structure, one that is quite
#          different from the initial structure.  consider calling
#          lsodes again with istate = 3 to force the structure to be
#          reevaluated.  see the full description of istate below.
#
# f. to continue the integration after a successful return, simply
# reset tout and call lsodes again.  no other parameters need be reset.
#
#-----------------------------------------------------------------------
# full description of user interface to lsodes.
#
# the user interface to lsodes consists of the following parts.
#
# i.   the call sequence to subroutine lsodes, which is a driver
#      routine for the solver.  this includes descriptions of both
#      the call sequence arguments and of user-supplied routines.
#      following these descriptions is a description of
#      optional inputs available through the call sequence, and then
#      a description of optional outputs (in the work arrays).
#
# ii.  descriptions of other routines in the lsodes package that may be
#      (optionally) called by the user.  these provide the ability to
#      alter error message handling, save and restore the internal
#      common, and obtain specified derivatives of the solution y(t).
#
# iii. descriptions of common blocks to be declared in overlay
#      or similar environments, or to be saved when doing an interrupt
#      of the problem and continued solution later.
#
# iv.  description of two routines in the lsodes package, either of
#      which the user may replace with his own version, if desired.
#      these relate to the measurement of errors.
#
#-----------------------------------------------------------------------
# part i.  call sequence.
#
# the call sequence parameters used for input only are
#     f, neq, tout, itol, rtol, atol, itask, iopt, lrw, liw, jac, mf,
# and those used for both input and output are
#     y, t, istate.
# the work arrays rwork and iwork are also used for conditional and
# optional inputs and optional outputs.  (the term output here refers
# to the return from subroutine lsodes to the user-s calling program.)
#
# the legality of input parameters will be thoroughly checked on the
# initial call for the problem, but not checked thereafter unless a
# change in input parameters is flagged by istate = 3 on input.
#
# the descriptions of the call arguments are as follows.
#
# f      = the name of the user-supplied subroutine defining the
#          ode system.  the system must be put in the first-order
#          form dy/dt = f(t,y), where f is a vector-valued function
#          of the scalar t and the vector y.  subroutine f is to
#          compute the function f.  it is to have the form
#               subroutine f (neq, t, y, ydot)
#               dimension y(1), ydot(1)
#          where neq, t, and y are input, and the array ydot = f(t,y)
#          is output.  y and ydot are arrays of length neq.
#          (in the dimension statement above, 1 is a dummy
#          dimension.. it can be replaced by any value.)
#          subroutine f should not alter y(1),...,y(neq).
#          f must be declared external in the calling program.
#
#
#          if quantities computed in the f routine are needed
#          externally to lsodes, an extra call to f should be made
#          for this purpose, for consistent and accurate results.
#          if only the derivative dy/dt is needed, use intdy instead.
#
# neq    = the size of the ode system (number of first order
#          ordinary differential equations).  used only for input.
#          neq may be decreased, but not increased, during the problem.
#          if neq is decreased (with istate = 3 on input), the
#          remaining components of y should be left undisturbed, if
#          these are to be accessed in f and/or jac.
#
#          normally, neq is a scalar, and it is generally referred to
#          as a scalar in this user interface description.  however,
#          neq may be an array, with neq set to the system size.
#          (the lsodes package accesses only neq.)
#
# y      = a real array for the vector of dependent variables, of
#          length neq or more.  used for both input and output on the
#          first call (istate = 1), and only for output on other calls.
#          on the first call, y must contain the vector of initial
#          values.  on output, y contains the computed solution vector,
#          evaluated at t.  if desired, the y array may be used
#          for other purposes between calls to the solver.
#
#          this array is passed as the y argument in all calls to
#          f and jac.  hence its length may exceed neq, and locations
#          y(neq+1),... may be used to store other real data and
#          pass it to f and/or jac.  (the lsodes package accesses only
#          y(1),...,y(neq).)
#
# t      = the independent variable.  on input, t is used only on the
#          first call, as the initial point of the integration.
#          on output, after each call, t is the value at which a
#          computed solution y is evaluated (usually the same as tout).
#          on an error return, t is the farthest point reached.
#
# tout   = the next value of t at which a computed solution is desired.
#          used only for input.
#
#          when starting the problem (istate = 1), tout may be equal
#          to t for one call, then should != t for the next call.
#          for the initial t, an input value of tout != t is used
#          in order to determine the direction of the integration
#          (i.e. the algebraic sign of the step sizes) and the rough
#          scale of the problem.  integration in either direction
#          (forward or backward in t) is permitted.
#
#          if itask = 2 or 5 (one-step modes), tout is ignored after
#          the first call (i.e. the first call with tout != t).
#          otherwise, tout is required on every call.
#
#          if itask = 1, 3, or 4, the values of tout need not be
#          monotone, but a value of tout which backs up is limited
#          to the current internal t interval, whose endpoints are
#          tcur - hu and tcur (see optional outputs, below, for
#          tcur and hu).
#
# itol   = an indicator for the type of error control.  see
#          description below under atol.  used only for input.
#
# rtol   = a relative error tolerance parameter, either a scalar or
#          an array of length neq.  see description below under atol.
#          input only.
#
# atol   = an absolute error tolerance parameter, either a scalar or
#          an array of length neq.  input only.
#
#             the input parameters itol, rtol, and atol determine
#          the error control performed by the solver.  the solver will
#          control the vector e = (e(i)) of estimated local errors
#          in y, according to an inequality of the form
#                      rms-norm of ( e(i)/ewt(i) )   <=   1,
#          where       ewt(i) = rtol(i)*abs(y(i)) + atol(i),
#          and the rms-norm (root-mean-square norm) here is
#          rms-norm(v) = sqrt(sum v(i)**2 / neq).  here ewt = (ewt(i))
#          is a vector of weights which must always be positive, and
#          the values of rtol and atol should all be non-negative.
#          the following table gives the types (scalar/array) of
#          rtol and atol, and the corresponding form of ewt(i).
#
#             itol    rtol       atol          ewt(i)
#              1     scalar     scalar     rtol*abs(y(i)) + atol
#              2     scalar     array      rtol*abs(y(i)) + atol(i)
#              3     array      scalar     rtol(i)*abs(y(i)) + atol
#              4     array      array      rtol(i)*abs(y(i)) + atol(i)
#
#          when either of these parameters is a scalar, it need not
#          be dimensioned in the user-s calling program.
#
#          if none of the above choices (with itol, rtol, and atol
#          fixed throughout the problem) is suitable, more general
#          error controls can be obtained by substituting
#          user-supplied routines for the setting of ewt and/or for
#          the norm calculation.  see part iv below.
#
#          if global errors are to be estimated by making a repeated
#          run on the same problem with smaller tolerances, then all
#          components of rtol and atol (i.e. of ewt) should be scaled
#          down uniformly.
#
# itask  = an index specifying the task to be performed.
#          input only.  itask has the following values and meanings.
#          1  means normal computation of output values of y(t) at
#             t = tout (by overshooting and interpolating).
#          2  means take one step only and return.
#          3  means stop at the first internal mesh point at or
#             beyond t = tout and return.
#          4  means normal computation of output values of y(t) at
#             t = tout but without overshooting t = tcrit.
#             tcrit must be input as rwork(1).  tcrit may be equal to
#             or beyond tout, but not behind it in the direction of
#             integration.  this option is useful if the problem
#             has a singularity at or beyond t = tcrit.
#          5  means take one step, without passing tcrit, and return.
#             tcrit must be input as rwork(1).
#
#          note..  if itask = 4 or 5 and the solver reaches tcrit
#          (within roundoff), it will return t = tcrit (exactly) to
#          indicate this (unless itask = 4 and tout comes before tcrit,
#          in which case answers at t = tout are returned first).
#
# istate = an index used for input and output to specify the
#          the state of the calculation.
#
#          on input, the values of istate are as follows.
#          1  means this is the first call for the problem
#             (initializations will be done).  see note below.
#          2  means this is not the first call, and the calculation
#             is to continue normally, with no change in any input
#             parameters except possibly tout and itask.
#             (if itol, rtol, and/or atol are changed between calls
#             with istate = 2, the new values will be used but not
#             tested for legality.)
#          3  means this is not the first call, and the
#             calculation is to continue normally, but with
#             a change in input parameters other than
#             tout and itask.  changes are allowed in
#             neq, itol, rtol, atol, iopt, lrw, liw, mf,
#             the conditional inputs ia and ja,
#             and any of the optional inputs except h0.
#             in particular, if miter = 1 or 2, a call with istate = 3
#             will cause the sparsity structure of the problem to be
#             recomputed (or reread from ia and ja if moss = 0).
#          note..  a preliminary call with tout = t is not counted
#          as a first call here, as no initialization or checking of
#          input is done.  (such a call is sometimes useful for the
#          purpose of outputting the initial conditions.)
#          thus the first call for which tout != t requires
#          istate = 1 on input.
#
#          on output, istate has the following values and meanings.
#           1  means nothing was done, as tout was equal to t with
#              istate = 1 on input.  (however, an internal counter was
#              set to detect and prevent repeated calls of this type.)
#           2  means the integration was performed successfully.
#          -1  means an excessive amount of work (more than mxstep
#              steps) was done on this call, before completing the
#              requested task, but the integration was otherwise
#              successful as far as t.  (mxstep is an optional input
#              and is normally 500.)  to continue, the user may
#              simply reset istate to a value > 1 and call again
#              (the excess work step counter will be reset to 0).
#              in addition, the user may increase mxstep to avoid
#              this error return (see below on optional inputs).
#          -2  means too much accuracy was requested for the precision
#              of the machine being used.  this was detected before
#              completing the requested task, but the integration
#              was successful as far as t.  to continue, the tolerance
#              parameters must be reset, and istate must be set
#              to 3.  the optional output tolsf may be used for this
#              purpose.  (note.. if this condition is detected before
#              taking any steps, then an illegal input return
#              (istate = -3) occurs instead.)
#          -3  means illegal input was detected, before taking any
#              integration steps.  see written message for details.
#              note..  if the solver detects an infinite loop of calls
#              to the solver with illegal input, it will cause
#              the run to stop.
#          -4  means there were repeated error test failures on
#              one attempted step, before completing the requested
#              task, but the integration was successful as far as t.
#              the problem may have a singularity, or the input
#              may be inappropriate.
#          -5  means there were repeated convergence test failures on
#              one attempted step, before completing the requested
#              task, but the integration was successful as far as t.
#              this may be caused by an inaccurate jacobian matrix,
#              if one is being used.
#          -6  means ewt(i) became zero for some i during the
#              integration.  pure relative error control (atol(i)=0.0)
#              was requested on a variable which has now vanished.
#              the integration was successful as far as t.
#          -7  means a fatal error return flag came from the sparse
#              solver cdrv by way of prjs or slss (numerical
#              factorization or backsolve).  this should never happen.
#              the integration was successful as far as t.
#
#          note.. an error return with istate = -1, -4, or -5 and with
#          miter = 1 or 2 may mean that the sparsity structure of the
#          problem has changed significantly since it was last
#          determined (or input).  in that case, one can attempt to
#          complete the integration by setting istate = 3 on the next
#          call, so that a new structure determination is done.
#
#          note..  since the normal output value of istate is 2,
#          it does not need to be reset for normal continuation.
#          also, since a negative input value of istate will be
#          regarded as illegal, a negative output value requires the
#          user to change it, and possibly other inputs, before
#          calling the solver again.
#
# iopt   = an integer flag to specify whether or not any optional
#          inputs are being used on this call.  input only.
#          the optional inputs are listed separately below.
#          iopt = 0 means no optional inputs are being used.
#                   default values will be used in all cases.
#          iopt = 1 means one or more optional inputs are being used.
#
# rwork  = a work array used for a mixture of real (double precision)
#          and integer work space.
#          the length of rwork (in real words) must be at least
#             20 + nyh*(maxord + 1) + 3*neq + lwm    where
#          nyh    = the initial value of neq,
#          maxord = 12 (if meth = 1) or 5 (if meth = 2) (unless a
#                   smaller value is given as an optional input),
#          lwm = 0                                    if miter = 0,
#          lwm = 2*nnz + 2*neq + (nnz+9*neq)/2   if miter = 1,
#          lwm = 2*nnz + 2*neq + (nnz+10*neq)/2  if miter = 2,
#          lwm = neq + 2                              if miter = 3.
#          in the above formulas,
#          nnz    = number of nonzero elements in the jacobian matrix.
#          (see the mf description for meth and miter.)
#          thus if maxord has its default value and neq is constant,
#          the minimum length of rwork is..
#             20 + 16*neq        for mf = 10,
#             20 + 16*neq + lwm  for mf = 11, 111, 211, 12, 112, 212,
#             22 + 17*neq        for mf = 13,
#             20 +  9*neq        for mf = 20,
#             20 +  9*neq + lwm  for mf = 21, 121, 221, 22, 122, 222,
#             22 + 10*neq        for mf = 23.
#          if miter = 1 or 2, the above formula for lwm is only a
#          crude lower bound.  the required length of rwork cannot
#          be readily predicted in general, as it depends on the
#          sparsity structure of the problem.  some experimentation
#          may be necessary.
#
#          the first 20 words of rwork are reserved for conditional
#          and optional inputs and optional outputs.
#
#          the following word in rwork is a conditional input..
#            rwork(1) = tcrit = critical value of t which the solver
#                       is not to overshoot.  required if itask is
#                       4 or 5, and ignored otherwise.  (see itask.)
#
# lrw    = the length of the array rwork, as declared by the user.
#          (this will be checked by the solver.)
#
# iwork  = an integer work array.  the length of iwork must be at least
#             31 + neq + nnz   if moss = 0 and miter = 1 or 2, or
#             30               otherwise.
#          (nnz is the number of nonzero elements in df/dy.)
#
#          in lsodes, iwork is used only for conditional and
#          optional inputs and optional outputs.
#
#          the following two blocks of words in iwork are conditional
#          inputs, required if moss = 0 and miter = 1 or 2, but not
#          otherwise (see the description of mf for moss).
#            iwork(30+j) = ia(j)     (j=1,...,neq+1)
#            iwork(31+neq+k) = ja(k) (k=1,...,nnz)
#          the two arrays ia and ja describe the sparsity structure
#          to be assumed for the jacobian matrix.  ja contains the row
#          indices where nonzero elements occur, reading in columnwise
#          order, and ia contains the starting locations in ja of the
#          descriptions of columns 1,...,neq, in that order, with
#          ia(1) = 1.  thus, for each column index j = 1,...,neq, the
#          values of the row index i in column j where a nonzero
#          element may occur are given by
#            i = ja(k), where   ia(j) <= k < ia(j+1).
#          if nnz is the total number of nonzero locations assumed,
#          then the length of the ja array is nnz, and ia(neq+1) must
#          be nnz + 1.  duplicate entries are not allowed.
#
# liw    = the length of the array iwork, as declared by the user.
#          (this will be checked by the solver.)
#
# note..  the work arrays must not be altered between calls to lsodes
# for the same problem, except possibly for the conditional and
# optional inputs, and except for the last 3*neq words of rwork.
# the latter space is used for internal scratch space, and so is
# available for use by the user outside lsodes between calls, if
# desired (but not for use by f or jac).
#
# jac    = name of user-supplied routine (miter = 1 or moss = 1) to
#          compute the jacobian matrix, df/dy, as a function of
#          the scalar t and the vector y.  it is to have the form
#               subroutine jac (neq, t, y, j, ian, jan, pdj)
#               dimension y(1), ian(1), jan(1), pdj(1)
#          where neq, t, y, j, ian, and jan are input, and the array
#          pdj, of length neq, is to be loaded with column j
#          of the jacobian on output.  thus df(i)/dy(j) is to be
#          loaded into pdj(i) for all relevant values of i.
#          here t and y have the same meaning as in subroutine f,
#          and j is a column index (1 to neq).  ian and jan are
#          undefined in calls to jac for structure determination
#          (moss = 1).  otherwise, ian and jan are structure
#          descriptors, as defined under optional outputs below, and
#          so can be used to determine the relevant row indices i, if
#          desired.  (in the dimension statement above, 1 is a
#          dummy dimension.. it can be replaced by any value.)
#               jac need not provide df/dy exactly.  a crude
#          approximation (possibly with greater sparsity) will do.
#               in any case, pdj is preset to zero by the solver,
#          so that only the nonzero elements need be loaded by jac.
#          calls to jac are made with j = 1,...,neq, in that order, and
#          each such set of calls is preceded by a call to f with the
#          same arguments neq, t, and y.  thus to gain some efficiency,
#          intermediate quantities shared by both calculations may be
#          saved in a user common block by f and not recomputed by jac,
#          if desired.  jac must not alter its input arguments.
#          jac must be declared external in the calling program.
#
# mf     = the method flag.  used only for input.
#          mf has three decimal digits-- moss, meth, miter--
#             mf = 100*moss + 10*meth + miter.
#          moss indicates the method to be used to obtain the sparsity
#          structure of the jacobian matrix if miter = 1 or 2..
#            moss = 0 means the user has supplied ia and ja
#                     (see descriptions under iwork above).
#            moss = 1 means the user has supplied jac (see below)
#                     and the structure will be obtained from neq
#                     initial calls to jac.
#            moss = 2 means the structure will be obtained from neq+1
#                     initial calls to f.
#          meth indicates the basic linear multistep method..
#            meth = 1 means the implicit adams method.
#            meth = 2 means the method based on backward
#                     differentiation formulas (bdf-s).
#          miter indicates the corrector iteration method..
#            miter = 0 means functional iteration (no jacobian matrix is involved).
#            miter = 1 means chord iteration with a user-supplied sparse jacobian, given by subroutine jac.
#            miter = 2 means chord iteration with an internally generated (difference quotient) sparse jacobian
#                      (using ngp extra calls to f per df/dy value, where ngp is an optional output described below.)
#            miter = 3 means chord iteration with an internally generated diagonal jacobian approximation.
#                      (using 1 extra call to f per df/dy evaluation).
#          if miter = 1 or moss = 1, the user must supply a subroutine
#          jac (the name is arbitrary) as described above under jac.
#          otherwise, a dummy argument can be used.
#
#          the standard choices for mf are..
#            mf = 10  for a nonstiff problem,
#            mf = 21 or 22 for a stiff problem with ia/ja supplied
#                     (21 if jac is supplied, 22 if not),
#            mf = 121 for a stiff problem with jac supplied,
#                     but not ia/ja,
#            mf = 222 for a stiff problem with neither ia/ja nor
#                     jac supplied.
#          the sparseness structure can be changed during the
#          problem by making a call to lsodes with istate = 3.
#-----------------------------------------------------------------------
# optional inputs.
#
# the following is a list of the optional inputs provided for in the
# call sequence.  (see also part ii.)  for each such input variable,
# this table lists its name as used in this documentation, its
# location in the call sequence, its meaning, and the default value.
# the use of any of these inputs requires iopt = 1, and in that
# case all of these inputs are examined.  a value of zero for any
# of these optional inputs will cause the default value to be used.
# thus to use a subset of the optional inputs, simply preload
# locations 5 to 10 in rwork and iwork to 0.0 and 0 respectively, and
# then set those of interest to nonzero values.
#
# name    location      meaning and default value
#
# h0      rwork(5)  the step size to be attempted on the first step.
#                   the default value is determined by the solver.
#
# hmax    rwork(6)  the maximum absolute step size allowed.
#                   the default value is infinite.
#
# hmin    rwork(7)  the minimum absolute step size allowed.
#                   the default value is 0.  (this lower bound is not
#                   enforced on the final step before reaching tcrit
#                   when itask = 4 or 5.)
#
# seth    rwork(8)  the element threshhold for sparsity determination
#                   when moss = 1 or 2.  if the absolute value of
#                   an estimated jacobian element is <= seth, it
#                   will be assumed to be absent in the structure.
#                   the default value of seth is 0.
#
# maxord  iwork(5)  the maximum order to be allowed.  the default
#                   value is 12 if meth = 1, and 5 if meth = 2.
#                   if maxord exceeds the default value, it will
#                   be reduced to the default value.
#                   if maxord is changed during the problem, it may
#                   cause the current order to be reduced.
#
# mxstep  iwork(6)  maximum number of (internally defined) steps
#                   allowed during one call to the solver.
#                   the default value is 500.
#
# mxhnil  iwork(7)  maximum number of messages printed (per problem)
#                   warning that t + h = t on a step (h = step size).
#                   this must be positive to result in a non-default
#                   value.  the default value is 10.
#-----------------------------------------------------------------------
# optional outputs.
#
# as optional additional output from lsodes, the variables listed
# below are quantities related to the performance of lsodes
# which are available to the user.  these are communicated by way of
# the work arrays, but also have internal mnemonic names as shown.
# except where stated otherwise, all of these outputs are defined
# on any successful return from lsodes, and on any return with
# istate = -1, -2, -4, -5, or -6.  on an illegal input return
# (istate = -3), they will be unchanged from their existing values
# (if any), except possibly for tolsf, lenrw, and leniw.
# on any error return, outputs relevant to the error will be defined,
# as noted below.
#
# name    location      meaning
#
# hu      rwork(11) the step size in t last used (successfully).
#
# hcur    rwork(12) the step size to be attempted on the next step.
#
# tcur    rwork(13) the current value of the independent variable
#                   which the solver has actually reached, i.e. the
#                   current internal mesh point in t.  on output, tcur
#                   will always be at least as far as the argument
#                   t, but may be farther (if interpolation was done).
#
# tolsf   rwork(14) a tolerance scale factor, greater than 1.0,
#                   computed when a request for too much accuracy was
#                   detected (istate = -3 if detected at the start of
#                   the problem, istate = -2 otherwise).  if itol is
#                   left unaltered but rtol and atol are uniformly
#                   scaled up by a factor of tolsf for the next call,
#                   then the solver is deemed likely to succeed.
#                   (the user may also ignore tolsf and alter the
#                   tolerance parameters in any other way appropriate.)
#
# nst     iwork(11) the number of steps taken for the problem so far.
#
# nfe     iwork(12) the number of f evaluations for the problem so far,
#                   excluding those for structure determination
#                   (moss = 2).
#
# nje     iwork(13) the number of jacobian evaluations for the problem
#                   so far, excluding those for structure determination
#                   (moss = 1).
#
# nqu     iwork(14) the method order last used (successfully).
#
# nqcur   iwork(15) the order to be attempted on the next step.
#
# imxer   iwork(16) the index of the component of largest magnitude in
#                   the weighted local error vector ( e(i)/ewt(i) ),
#                   on an error return with istate = -4 or -5.
#
# lenrw   iwork(17) the length of rwork actually required.
#                   this is defined on normal returns and on an illegal
#                   input return for insufficient storage.
#
# leniw   iwork(18) the length of iwork actually required.
#                   this is defined on normal returns and on an illegal
#                   input return for insufficient storage.
#
# nnz     iwork(19) the number of nonzero elements in the jacobian
#                   matrix, including the diagonal (miter = 1 or 2).
#                   (this may differ from that given by ia(neq+1)-1
#                   if moss = 0, because of added diagonal entries.)
#
# ngp     iwork(20) the number of groups of column indices, used in
#                   difference quotient jacobian aproximations if
#                   miter = 2.  this is also the number of extra f
#                   evaluations needed for each jacobian evaluation.
#
# nlu     iwork(21) the number of sparse lu decompositions for the
#                   problem so far.
#
# lyh     iwork(22) the base address in rwork of the history array yh,
#                   described below in this list.
#
# ipian   iwork(23) the base address of the structure descriptor array
#                   ian, described below in this list.
#
# ipjan   iwork(24) the base address of the structure descriptor array
#                   jan, described below in this list.
#
# nzl     iwork(25) the number of nonzero elements in the strict lower
#                   triangle of the lu factorization used in the chord
#                   iteration (miter = 1 or 2).
#
# nzu     iwork(26) the number of nonzero elements in the strict upper
#                   triangle of the lu factorization used in the chord
#                   iteration (miter = 1 or 2).
#                   the total number of nonzeros in the factorization
#                   is therefore nzl + nzu + neq.
#
# the following four arrays are segments of the rwork array which
# may also be of interest to the user as optional outputs.
# for each array, the table below gives its internal name,
# its base address, and its description.
# for yh and acor, the base addresses are in rwork (a real array).
# the integer arrays ian and jan are to be obtained by declaring an
# integer array iwk and identifying iwk(1) with rwork(21), using either
# an equivalence statement or a subroutine call.  then the base
# addresses ipian (of ian) and ipjan (of jan) in iwk are to be obtained
# as optional outputs iwork(23) and iwork(24), respectively.
# thus ian(1) is iwk(ipian), etc.
#
# name    base address      description
#
# ian    ipian (in iwk)  structure descriptor array of size neq + 1.
# jan    ipjan (in iwk)  structure descriptor array of size nnz.
#         (see above)    ian and jan together describe the sparsity
#                        structure of the jacobian matrix, as used by
#                        lsodes when miter = 1 or 2.
#                        jan contains the row indices of the nonzero
#                        locations, reading in columnwise order, and
#                        ian contains the starting locations in jan of
#                        the descriptions of columns 1,...,neq, in
#                        that order, with ian(1) = 1.  thus for each
#                        j = 1,...,neq, the row indices i of the
#                        nonzero locations in column j are
#                        i = jan(k), ian(j) <= k < ian(j+1).
#                        note that ian(neq+1) = nnz + 1.
#                        (if moss = 0, ian/jan may differ from the
#                        input ia/ja because of a different ordering
#                        in each column, and added diagonal entries.)
#
# yh      lyh            the nordsieck history array, of size nyh by
#          (optional     (nqcur + 1), where nyh is the initial value
#          output)       of neq.  for j = 0,1,...,nqcur, column j+1
#                        of yh contains hcur**j/factorial(j) times
#                        the j-th derivative of the interpolating
#                        polynomial currently representing the solution,
#                        evaluated at t = tcur.  the base address lyh
#                        is another optional output, listed above.
#
# acor     lenrw-neq+1   array of size neq used for the accumulated
#                        corrections on each step, scaled on output
#                        to represent the estimated local error in y
#                        on the last step.  this is the vector e in
#                        the description of the error control.  it is
#                        defined only on a successful return from
#                        lsodes.
#
#-----------------------------------------------------------------------
# part ii.  other routines callable.
#
# the following are optional calls which the user may make to
# gain additional capabilities in conjunction with lsodes.
# (the routines xsetun and xsetf are designed to conform to the
# slatec error handling package.)
#
#     form of call                  function
#   call xsetun(lun)          set the logical unit number, lun, for
#                             output of messages from lsodes, if
#                             the default is not desired.
#                             the default value of lun is 6.
#
#   call xsetf(mflag)         set a flag to control the printing of
#                             messages by lsodes.
#                             mflag = 0 means make not print. (danger..
#                             this risks losing valuable information.)
#                             mflag = 1 means print (the default).
#
#                             either of the above calls may be made at
#                             any time and will take effect immediately.
#
#   call srcms(rsav,isav,job) saves and restores the contents of
#                             the internal common blocks used by
#                             lsodes (see part iii below).
#                             rsav must be a real array of length 224
#                             or more, and isav must be an integer
#                             array of length 75 or more.
#                             job=1 means save common into rsav/isav.
#                             job=2 means restore common from rsav/isav.
#                                srcms is useful if one is
#                             interrupting a run and restarting
#                             later, or alternating between two or
#                             more problems solved with lsodes.
#
#   call intdy(,,,,,)         provide derivatives of y, of various
#        (see below)          orders, at a specified point t, if
#                             desired.  it may be called only after
#                             a successful return from lsodes.
#
# the detailed instructions for using intdy are as follows.
# the form of the call is..
#
#   lyh = iwork(22)
#   call intdy (t, k, rwork(lyh), nyh, dky, iflag)
#
# the input parameters are..
#
# t         = value of independent variable where answers are desired
#             (normally the same as the t last returned by lsodes).
#             for valid results, t must lie between tcur - hu and tcur.
#             (see optional outputs for tcur and hu.)
# k         = integer order of the derivative desired.  k must satisfy
#             0 <= k <= nqcur, where nqcur is the current order
#             (see optional outputs).  the capability corresponding
#             to k = 0, i.e. computing y(t), is already provided
#             by lsodes directly.  since nqcur >= 1, the first
#             derivative dy/dt is always available with intdy.
# lyh       = the base address of the history array yh, obtained
#             as an optional output as shown above.
# nyh       = column length of yh, equal to the initial value of neq.
#
# the output parameters are..
#
# dky       = a real array of length neq containing the computed value
#             of the k-th derivative of y(t).
# iflag     = integer flag, returned as 0 if k and t were legal,
#             -1 if k was illegal, and -2 if t was illegal.
#             on an error return, a message is also written.
#-----------------------------------------------------------------------
# part iii.  common blocks.
#
# if lsodes is to be used in an overlay situation, the user
# must declare, in the primary overlay, the variables in..
#   (1) the call sequence to lsodes,
#   (2) the three internal common blocks
#         /ls0001/  of length  257  (218 double precision words
#                         followed by 39 integer words),
#         /lss001/  of length  40    ( 6 double precision words
#                         followed by 34 integer words),
#         /eh0001/  of length  2 (integer words).
#
# if lsodes is used on a system in which the contents of internal
# common blocks are not preserved between calls, the user should
# declare the above three common blocks in his main program to insure
# that their contents are preserved.
#
# if the solution of a given problem by lsodes is to be interrupted
# and then later continued, such as when restarting an interrupted run
# or alternating between two or more problems, the user should save,
# following the return from the last lsodes call prior to the
# interruption, the contents of the call sequence variables and the
# internal common blocks, and later restore these values before the
# next lsodes call for that problem.  to save and restore the common
# blocks, use subroutine srcms (see part ii above).
#
#-----------------------------------------------------------------------
# part iv.  optionally replaceable solver routines.
#
# below are descriptions of two routines in the lsodes package which
# relate to the measurement of errors.  either routine can be
# replaced by a user-supplied version, if desired.  however, since such
# a replacement may have a major impact on performance, it should be
# done only when absolutely necessary, and only with great caution.
# (note.. the means by which the package version of a routine is
# superseded by the user-s version may be system-dependent.)
#
# (a) ewset.
# the following subroutine is called just before each internal
# integration step, and sets the array of error weights, ewt, as
# described under itol/rtol/atol above..
#     subroutine ewset (neq, itol, rtol, atol, ycur, ewt)
# where neq, itol, rtol, and atol are as in the lsodes call sequence,
# ycur contains the current dependent variable vector, and
# ewt is the array of weights set by ewset.
#
# if the user supplies this subroutine, it must return in ewt(i)
# (i = 1,...,neq) a positive quantity suitable for comparing errors
# in y(i) to.  the ewt array returned by ewset is passed to the
# vnorm routine (see below), and also used by lsodes in the computation
# of the optional output imxer, the diagonal jacobian approximation,
# and the increments for difference quotient jacobians.
#
# in the user-supplied version of ewset, it may be desirable to use
# the current values of derivatives of y.  derivatives up to order nq
# are available from the history array yh, described above under
# optional outputs.  in ewset, yh is identical to the ycur array,
# extended to nq + 1 columns with a column length of nyh and scale
# factors of h**j/factorial(j).  on the first call for the problem,
# given by nst = 0, nq is 1 and h is temporarily set to 1.0.
# the quantities nq, nyh, h, and nst can be obtained by including
# in ewset the statements..
#     double precision h, rls
#     common /ls0001/ rls(218),ils(39)
#     nq = ils(35)
#     nyh = ils(14)
#     nst = ils(36)
#     h = rls(212)
# thus, for example, the current value of dy/dt can be obtained as
# ycur(nyh+i)/h  (i=1,...,neq)  (and the division by h is
# unnecessary when nst = 0).
#
# (b) vnorm.
# the following is a real function routine which computes the weighted
# root-mean-square norm of a vector v..
#     d = vnorm (n, v, w)
# where..
#   n = the length of the vector,
#   v = real array of length n containing the vector,
#   w = real array of length n containing weights,
#   d = sqrt( (1/n) * sum(v(i)*w(i))**2 ).
# vnorm is called with n = neq and with w(i) = 1.0/ewt(i), where
# ewt is as set by subroutine ewset.
#
# if the user supplies this function, it should return a non-negative
# value of vnorm suitable for use in the error control in lsodes.
# none of the arguments should be altered by vnorm.
# for example, a user-supplied vnorm routine might..
#   -substitute a max-norm of (v(i)*w(i)) for the rms-norm, or
#   -ignore some components of v in the norm, with the effect of
#    suppressing the error control on those components of y.
#-----------------------------------------------------------------------
#-----------------------------------------------------------------------
# other routines in the lsodes package.
#
# in addition to subroutine lsodes, the lsodes package includes the
# following subroutines and function routines..
#  iprep    acts as an iterface between lsodes and prep, and also does
#           adjusting of work space pointers and work arrays.
#  prep     is called by iprep to compute sparsity and make sparse matrix
#           preprocessing if miter = 1 or 2.
#  jgroup   is called by prep to compute groups of jacobian column
#           indices for use when miter = 2.
#  adjlr    adjusts the length of required sparse matrix work space.
#           it is called by prep.
#  cntnzu   is called by prep and counts the nonzero elements in the
#           strict upper triangle of j + j-transpose, where j = df/dy.
#  intdy    computes an interpolated value of the y vector at t = tout.
#  stode    is the core integrator, which does one step of the
#           integration and the associated error control.
#  cfode    sets all method coefficients and test constants.
#  prjs     computes and preprocesses the jacobian matrix j = df/dy
#           and the newton iteration matrix p = i - h*l0*j.
#  slss     manages solution of linear system in chord iteration.
#  ewset    sets the error weight vector ewt before each step.
#  vnorm    computes the weighted r.m.s. norm of a vector.
#  srcms    is a user-callable routine to save and restore
#           the contents of the internal common blocks.
#  odrv     constructs a reordering of the rows and columns of
#           a matrix by the minimum degree algorithm.  odrv is a
#           driver routine which calls subroutines md, mdi, mdm,
#           mdp, mdu, and sro.  see ref. 2 for details.  (the odrv
#           module has been modified since ref. 2, however.)
#  cdrv     performs reordering, symbolic factorization, numerical
#           factorization, or linear system solution operations,
#           depending on a path argument ipath.  cdrv is a
#           driver routine which calls subroutines nroc, nsfc,
#           nnfc, nnsc, and nntc.  see ref. 3 for details.
#           lsodes uses cdrv to solve linear systems in which the
#           coefficient matrix is  p = i - con*j, where i is the
#           identity, con is a scalar, and j is an approximation to
#           the jacobian df/dy.  because cdrv deals with rowwise
#           sparsity descriptions, cdrv works with p-transpose, not p.
#
#--------------------------------------------------------------------------------------------------
import math
import sys

def lsodes(neq, y, t, tout, itol, rtol, atol, itask, istate, iopt, rwork, lrw, iwork, liw, mf)

   integer neq, itol, itask, istate, iopt, lrw, iwork, liw, mf
   double precision y, t, tout, rtol, atol, rwork
   dimension y(1), rtol(1), atol(1), rwork(lrw), iwork(liw)
   integer i, i1, i2, iflag, imax, imul, imxer, ipflag, ipgo, irem, j, kgo, 
  +        lenyht, leniw, lenrw, lf0, lia, lja, lrtem, lwtem, lyhd, lyhn, mf1, mord, mxhnl0, ncolm
   double precision atoli, ayi, big, ewti, h0, hmax, hmx, rh, rtoli, tcrit, tdist, tnext, tol, tolsf, tp, size, sum, w0, vnorm
   dimension mord(2)
   logical ihit
 
#  the following two internal common blocks contain
#  (a) variables which are local to any subroutine but whose values must be preserved between calls to the routine (own variables), and
#  (b) variables which are communicated between subroutines. the structure of each block is as follows..  
#      all real variables are listed first, followed by all integers.  
#      within each type, the variables are grouped with those local to subroutine lsodes first,
#      then those local to subroutine stode or subroutine prjs (no other routines have own variables), and finally those used for communication.
#  the block ls0001 is declared in subroutines lsodes, iprep, prep, intdy, stode, prjs, and slss.  
#  the block lss001 is declared in subroutines lsodes, iprep, prep, prjs, and slss.
#  groups of variables are replaced by dummy arrays in the common declarations in routines where those variables are not used.
 
   double precision crate, el, elco, hold, rmax, tesco, 
  +   ccmax, el0, h, hmin, hmxi, hu, rc, tn, uround

   integer init, lyh, lewt, lacor, lsavf, lwm, liwm, mxstep, mxhnil, nhnil, ntrep, nslast, nyh,
  +   ialth, ipup, lmax, meo, nqnyh, nslp, 
  +   icf, ierpj, iersl, jcur, jstart, kflag, meth, miter, 
  +   maxord, maxcor, n, nq, nst, nfe, nje, nqu
   common /ls0001/ crate, el(13), elco(13,12), hold, rmax, tesco(3,12),
  +   ccmax, el0, h, hmin, hmxi, hu, rc, tn, uround,
  +   init, lyh, lewt, lacor, lsavf, lwm, liwm, mxstep, mxhnil, nhnil, ntrep, nslast, nyh, 
  +   ialth, ipup, lmax, meo, nqnyh, nslp, 
  +   icf, ierpj, iersl, jcur, jstart, kflag, meth, miter,
  +   maxord, maxcor, n, nq, nst, nfe, nje, nqu
 
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
 
#  this code block is executed on every call. it tests istate and itask for legality and branches appropriately.
#  if istate > 1 but the flag init shows that initialization has not yet been done, an error return occurs.
#  if istate = 1 and tout = t, return immediately.
   if istate < 1 or istate > 3 :
      write(*,*)'***ERROR lsodes: istate (=', istate, ') is illegal'
      sys.exit()
   if itask < 1 or itask > 5 :
      write(*,*)'***ERROR lsodes: itask (=', itask, ') is illegal'
      sys.exit()
   if istate != 1 and init == 0 :
      write(*,*)'***ERROR lsodes: istate > 1 but lsodes not initialized (init = 0)'
      sys.exit()
   if istate == 1 :
      init = 0
      if tout == t :
         write(*,*)'***ERROR lsodes: istate = 1 and tout = t (=', t, ')'
         sys.exit()
 
   if istate == 1 or istate == 3 :
      ntrep = 0
      
#     the next code block is executed for the initial call (istate = 1), or for a continuation call with parameter changes (istate = 3).
#     it contains checking of all inputs and various initializations. 
#     if istate = 1, the final setting of work space pointers, the matrix preprocessing, and other initializations are done.
#     first check legality of the non-optional inputs neq, itol, iopt, mf, ml, and mu.
      if neq <= 0 :
         write(*,*)'***ERROR lsodes: neq (=', neq, ') < 1'
         sys.exit()
      if istate != 1 and neq > n :
         write(*,*)'***ERROR lsodes: istate = 3 and neq increased (', n, ' to ', neq, ')'
         sys.exit()
      n = neq
      if itol < 1 or itol > 4 :
         write(*,*)'***ERROR lsodes: itol (=', itol, ') is illegal'
         sys.exit()
      if iopt < 0 or iopt > 1 :
         write(*,*)'***ERROR lsodes: iopt (=', iopt, ') is illegal'
         sys.exit()
      moss = mf/100
      mf1 = mf - 100*moss
      meth = mf1/10
      miter = mf1 - 10*meth
      if moss < 0 or moss > 2 :
         write(*,*)'***ERROR lsodes: mf (=', mf, ') is illegal'
         sys.exit()
      if meth < 1 or meth > 2 :
         write(*,*)'***ERROR lsodes: mf (=', mf, ') is illegal'
         sys.exit()
      if miter < 0 or miter > 3 :
         write(*,*)'***ERROR lsodes: mf (=', mf, ') is illegal'
         sys.exit()
      if miter == 0 or miter == 3 : moss = 0
#     next process and check the optional inputs.
      if iopt == 0 :
         maxord = mord(meth)
         mxstep = 1000000
         mxhnil = mxhnl0
         if istate == 1 : h0 = 0.0
         hmxi = 0.0
         hmin = 0.0
         seth = 0.0
      else : # iopt == 1
         maxord = iwork(5)
         if maxord < 0 :
            write(*,*)'***ERROR lsodes: maxord (=', maxord, ') < 0'
            sys.exit()
         if maxord == 0 : maxord = 100
         maxord = min(maxord,mord(meth))
         mxstep = iwork(6)
         if mxstep < 0 :
            write(*,*)'***ERROR lsodes: maxord (=', mxstep, ') < 0'
            sys.exit()
         if mxstep == 0 : mxstep = 500
         mxhnil = iwork(7)
         if mxhnil < 0 :
            write(*,*)'***ERROR lsodes: mxhnil (=', mxhnil, ') < 0'
            sys.exit()
         if mxhnil == 0 : mxhnil = mxhnl0
         if istate == 1 :
            h0 = rwork(5)
            if (tout - t)*h0 < 0.0 :
               write(*,*)'***ERROR lsodes: tout (=', tout, ') is behind t (', t, '). integration direction is given by h0 (=', h0, ')'
               sys.exit()
         hmax = rwork(6)
         if hmax < 0.0 :
            write(*,*)'***ERROR lsodes: hmax (=', hmax, ') < 0'
            sys.exit()
         hmxi = 0.0
         if hmax > 0.0 : hmxi = 1.0/hmax
         hmin = rwork(7)
         if hmin < 0.0 :
            write(*,*)'***ERROR lsodes: hmin (=', hmin, ') < 0'
            sys.exit()
         seth = rwork(8)
         if seth < 0.0 :
            write(*,*)'***ERROR lsodes: seth (=', seth, ') < 0'
            sys.exit()
      
#     check rtol and atol for legality.
      rtoli = rtol(1)
      atoli = atol(1)
      do i = 1,n
         if itol >= 3 : rtoli = rtol(i)
         if itol == 2 or itol == 4 : atoli = atol(i)
         if rtoli < 0.0 :
            write(*,*)'***ERROR lsodes: rtol(', i, ') = ', rtoli, '< 0'
            sys.exit()
         if atoli < 0.0 :
            write(*,*)'***ERROR lsodes: atol(', i, ') = ', atoli, '< 0'
            sys.exit()
      
#     compute required work array lengths, as far as possible, and test these against lrw and liw.  then set tentative pointers for work arrays.
#     pointers to rwork/iwork segments are named by prefixing l to the name of the segment.  e.g., the segment yh starts at rwork(lyh).
#     segments of rwork (in order) are denoted  wm, yh, savf, ewt, acor.
#     if miter = 1 or 2, the required length of the matrix work space wm is not yet known, 
#     and so a crude minimum value is used for the initial tests of lrw and liw, and yh is temporarily stored as far
#     to the right in rwork as possible, to leave the maximum amount of space for wm for matrix preprocessing.
#     thus if miter = 1 or 2 and moss != 2, some of the segments of rwork are temporarily omitted, as they are not needed in the preprocessing.
#     these omitted segments are.. acor if istate = 1, ewt and acor if istate = 3 and moss = 1, and savf, ewt, and acor if istate = 3 and moss = 0.
      if istate == 1 : nyh = n
      lwmin = 0
      if miter == 1 : lwmin = 4*n + 10*n/2
      if miter == 2 : lwmin = 4*n + 11*n/2
      if miter == 3 : lwmin = n + 2
      lenyh = (maxord+1)*nyh
      lrest = lenyh + 3*n
      lenrw = 20 + lwmin + lrest
      iwork(17) = lenrw
      leniw = 30
      if moss == 0 and miter != 0 and miter != 3 : leniw = leniw + n + 1
      iwork(18) = leniw
      if lenrw > lrw :
         write(*,*)'***ERROR lsodes: rwork length is insufficient to proceed. length needed is at least lenrw (=', lenrw, ') exceeds lrw (=', lrw, ')'
         sys.exit()
      if leniw > liw :
         write(*,*)'***ERROR lsodes: iwork length is insufficient to proceed. length needed is at least lenrw (=', leniw, ') exceeds lrw (=', liw, ')'
         sys.exit()
      lia = 31
      if moss == 0 and miter != 0 and miter != 3 : leniw = leniw + iwork(lia+n) - 1
      iwork(18) = leniw
      if leniw > liw :
         write(*,*)'***ERROR lsodes: iwork length is insufficient to proceed. length needed is at least lenrw (=', leniw, ') exceeds lrw (=', liw, ')'
         sys.exit()
      lja = lia + n + 1
      lia = min(lia,liw)
      lja = min(lja,liw)
      lwm = 21
      if istate == 1 : nq = 1
      ncolm = min(nq+1,maxord+2)
      lenyhm = ncolm*nyh
      lenyht = lenyh
      if miter == 1 or miter == 2 : lenyht = lenyhm
      imul = 2
      if istate == 3 : imul = moss
      if moss == 2 : imul = 3
      lrtem = lenyht + imul*n
      lwtem = lwmin
      if miter == 1 or miter == 2 : lwtem = lrw - 20 - lrtem
      lenwk = lwtem
      lyhn = lwm + lwtem
      lsavf = lyhn + lenyht
      lewt = lsavf + n
      lacor = lewt + n
      istatc = istate
      if istate == 1 :
#        the next block is for the initial call only (istate = 1).
#        it contains all remaining initializations, the initial call to f,
#        the sparse matrix preprocessing (miter = 1 or 2), and the
#        calculation of the initial step size.
#        the error weights in ewt are inverted after being loaded.
         lyh = lyhn
         iwork(22) = lyh
         tn = t
         nst = 0
         h = 1.0
         nnz = 0
         ngp = 0
         nzl = 0
         nzu = 0
#        load the initial value vector in yh.
         do i = 1,n
            rwork(i+lyh-1) = y(i)
#        initial call to f.  (lf0 points to yh(*,2).)
         lf0 = lyh + nyh
         call rhs(neq, t, y, rwork(lf0))
         nfe = 1
#        load and invert the ewt array. (h is temporarily set to 1.0.)
         call ewset(n, itol, rtol, atol, rwork(lyh), rwork(lewt))
         do i = 1,n
            if rwork(i+lewt-1) <= 0.0 :
               ewti = rwork(lewt+i-1)
               write(*,*)'***ERROR lsodes: ewt(', i, ') is (', ewti, ') non-positive.'
               sys.exit()
            rwork(i+lewt-1) = 1.0/rwork(i+lewt-1)
         if miter == 1 or miter == 2 :
#           iprep and prep make sparse matrix preprocessing if miter = 1 or 2.
            lacor = min(lacor,lrw)
            call iprep(neq, y, rwork, iwork(lia), iwork(lja), ipflag)
            lenrw = lwm - 1 + lenwk + lrest
            iwork(17) = lenrw
            if ipflag != -1 : iwork(23) = ipian
            if ipflag != -1 : iwork(24) = ipjan
            ipgo = -ipflag + 1
            if ipgo == 2 :
               write(*,*)'***ERROR lsodes: rwork length is insufficient for subroutine prep. length needed is at least lenrw (=', lenrw, ') exceeds lrw (=', lrw, ')'
               sys.exit()
            elif ipgo == 3 :
               write(*,*)'***ERROR lsodes: rwork length is insufficient for subroutine jgroup. length needed is at least lenrw (=', lenrw, ') exceeds lrw (=', lrw, ')'
               sys.exit()
            elif ipgo == 4 :
               write(*,*)'***ERROR lsodes: rwork length is insufficient for subroutine odrv. length needed is at least lenrw (=', lenrw, ') exceeds lrw (=', lrw, ')'
               sys.exit()
            elif ipgo == 5 :
               imul = (iys - 1)/n
               irem = iys - imul*n
               write(*,*)'***ERROR lsodes: error from odrv in yale sparse matrix package at t (=', tn, ') odrv returned error flag = ', imul, '*neq + ', irem, '.'
               sys.exit()
            elif ipgo == 6 :
               write(*,*)'***ERROR lsodes: rwork length is insufficient for subroutine cdrv. length needed is at least lenrw (=', lenrw, ') exceeds lrw (=', lrw, ')'
               sys.exit()
            elif ipgo == 7 :
               imul = (iys - 1)/n
               irem = iys - imul*n
               write(*,*)'***ERROR lsodes: error from cdrv in yale sparse matrix package at t (=', tn, ') cdrv returned error flag = ', imul, '*neq + ', irem, '.'
               if imul == 2 : write(*,*)'***ERROR lsodes: duplicate entry in sparsity structure descriptors'
               if imul == 3 or imul == 6 : write(*,*)'***ERROR lsodes: insufficient storage for nsfc (called by cdrv)'
               sys.exit()
 
            iwork(22) = lyh
            if lenrw > lrw :
               write(*,*)'***ERROR lsodes: rwork length is insufficient to proceed. length needed is at least lenrw (=', lenrw, ') exceeds lrw (=', lrw, ')'
               sys.exit()
#        check tcrit for legality (itask = 4 or 5).
         if itask == 4 or itask == 5 :
            tcrit = rwork(1)
            if (tcrit - tout)*(tout - t) < 0.0 :
               write(*,*)'***ERROR lsodes: itask = 4 or 5 and tcrit (=', tcrit, ') is behind tout (', tout, ').'
               sys.exit()
            if h0 != 0.0 and (t + h0 - tcrit)*h0 > 0.0 : h0 = tcrit - t
#        initialize all remaining parameters.
         uround = 2.2204460492503131E-016
         jstart = 0
         if miter != 0 : rwork(lwm) = math.sqrt(uround)
         msbj = 50
         nslj = 0
         ccmxj = 0.2
         psmall = 1000.0*uround
         rbig = 0.01/psmall
         nje = 0
         nlu = 0
         nslast = 0
         hu = 0.0
         nqu = 0
         ccmax = 0.3
         maxcor = 3
#        the coding below computes the step size, h0, to be attempted on the first step, unless the user has supplied a value for this.
#        first check that tout - t differs significantly from zero.
#        a scalar tolerance quantity tol is computed, as max(rtol(i)) if this is positive, or max(atol(i)/abs(y(i))) otherwise, 
#        adjusted so as to be between 100*uround and 1.0e-3.
#        then the computed value h0 is given by..
#                                             neq
#          h0**2 = tol / ( w0**-2 + (1/neq) * sum ( f(i)/ywt(i) )**2  )
#                                              1
#        where   w0     = max ( abs(t), abs(tout) ),
#                f(i)   = i-th component of initial value of f,
#                ywt(i) = ewt(i)/tol  (a weight for y(i)).
#        the sign of h0 is inferred from the initial values of tout and t.
         lf0 = lyh + nyh
         if h0 == 0.0 :
            tdist = abs(tout - t)
            w0 = max(abs(t), abs(tout))
            if tdist < 2.0*uround*w0 :
               write(*,*)'***ERROR lsodes: tout (=', tout, ') too close to t(=', t, ') to start integration.'
               sys.exit()
            tol = rtol(1)
            if itol > 2 :
               do i = 1,n
                  tol = max(tol,rtol(i))
            if tol <= 0.0 :
               atoli = atol(1)
               do i = 1,n
                  if itol == 2 or itol == 4 : atoli = atol(i)
                  ayi = abs(y(i))
                  if ayi != 0.0 : tol = max(tol,atoli/ayi)
            tol = max(tol,100.0*uround)
            tol = min(tol,0.001)
            sum = vnorm(n, rwork(lf0), rwork(lewt))
            sum = 1.0/(tol*w0*w0) + tol*sum**2
            h0 = 1.0/math.sqrt(sum)
            h0 = min(h0,tdist)
            h0 = math.copysign(h0,tout-t)
#        adjust h0 if necessary to meet hmax bound.
         rh = abs(h0)*hmxi
         if rh > 1.0 : h0 = h0/rh
#        load h with h0 and scale yh(*,2) by h0.
         h = h0
         do i = 1,n
            rwork(i+lf0-1) = h0*rwork(i+lf0-1)
 
      else : # istate == 3
      
#        istate = 3. move yh to its new location.
#        note that only the part of yh needed for the next step, namely min(nq+1,maxord+2) columns, is actually moved.
#        a temporary error weight array ewt is loaded if moss = 2. sparse matrix processing is done in iprep/prep if miter = 1 or 2.
#        if maxord was reduced below nq, then the pointers are finally set so that savf is identical to yh(*,maxord+2).
         lyhd = lyh - lyhn
         imax = lyhn - 1 + lenyhm
#        move yh. branch for move right, no move, or move left.
         if lyhd == 1 :
#           move right
            do i = lyhn,imax
               j = imax + lyhn - i
               rwork(j) = rwork(j+lyhd)
         elif lyhd == 3 :
#           move left
            do i = lyhn,imax
               rwork(i) = rwork(i+lyhd)
         lyh = lyhn
         iwork(22) = lyh
         if miter == 1 or miter == 2 :
            if moss == 2 :
#              temporarily load ewt if miter = 1 or 2 and moss = 2.
               call ewset (n, itol, rtol, atol, rwork(lyh), rwork(lewt))
               do i = 1,n
                  if rwork(i+lewt-1) <= 0.0 :
                     ewti = rwork(lewt+i-1)
                     write(*,*)'***ERROR lsodes: ewt (', i, ') is ', ewti, 'non-positive.'
                     sys.exit()
                  rwork(i+lewt-1) = 1.0/rwork(i+lewt-1)
#           iprep and prep make sparse matrix preprocessing if miter = 1 or 2.
            lsavf = min(lsavf,lrw)
            lewt = min(lewt,lrw)
            lacor = min(lacor,lrw)
            call iprep(neq, y, rwork, iwork(lia), iwork(lja), ipflag)
            lenrw = lwm - 1 + lenwk + lrest
            iwork(17) = lenrw
            if ipflag != -1 : iwork(23) = ipian
            if ipflag != -1 : iwork(24) = ipjan
            ipgo = -ipflag + 1
 
            if ipgo == 2 :
               write(*,*)'***ERROR lsodes: rwork length is insufficient for subroutine prep. length needed is at least lenrw (=', lenrw, ') exceeds lrw (=', lrw, ')'
               sys.exit()
            elif ipgo == 3 :
               write(*,*)'***ERROR lsodes: rwork length is insufficient for subroutine jgroup. length needed is at least lenrw (=', lenrw, ') exceeds lrw (=', lrw, ')'
               sys.exit()
            elif ipgo == 4 :
               write(*,*)'***ERROR lsodes: rwork length is insufficient for subroutine odrv. length needed is at least lenrw (=', lenrw, ') exceeds lrw (=', lrw, ')'
               sys.exit()
            elif ipgo == 5 :
               imul = (iys - 1)/n
               irem = iys - imul*n
               write(*,*)'***ERROR lsodes: error from odrv in yale sparse matrix package at t (=', tn, ') odrv returned error flag = ', imul, '*neq + ', irem, '.'
               sys.exit()
            elif ipgo == 6 :
               write(*,*)'***ERROR lsodes: rwork length is insufficient for subroutine cdrv. length needed is at least lenrw (=', lenrw, ') exceeds lrw (=', lrw, ')'
               sys.exit()
            elif ipgo == 7 :
               imul = (iys - 1)/n
               irem = iys - imul*n
               write(*,*)'***ERROR lsodes: error from cdrv in yale sparse matrix package at t (=', tn, ') cdrv returned error flag = ', imul, '*neq + ', irem, '.'
               if imul == 2 : write(*,*)'***ERROR lsodes: duplicate entry in sparsity structure descriptors'
               if imul == 3 or imul == 6 : write(*,*)'***ERROR lsodes: insufficient storage for nsfc (called by cdrv)'
               sys.exit()
 
            iwork(22) = lyh
            if lenrw > lrw :
               write(*,*)'***ERROR lsodes: rwork length is insufficient to proceed. length needed is at least lenrw (=', lenrw, ') exceeds lrw (=', lrw, ')'
               sys.exit()
 
#        set flag to signal parameter changes to stode.
         jstart = -1
         if n != nyh :
#           neq was reduced.  zero part of yh to avoid undefined references.
            i1 = lyh + (nq+1)*nyh
            i2 = lyh + (maxord + 1)*nyh - 1
            if i1 <= i2 :
               do i = i1,i2
                  rwork(i) = 0.0
 
   if istate == 2 or istate == 3 :
#     the next code block is for continuation calls only (istate = 2 or 3)
#     and is to check stop conditions before taking a step.
      nslast = nst
      if itask == 1 :
         if (tn - tout)*h >= 0.0 :
            call intdy (tout, 0, rwork(lyh), nyh, y, iflag)
            if iflag != 0 :
               write(*,*)'***ERROR lsodes: trouble from intdy. itask = ', itask, ', tout = ', tout, '.'
               sys.exit()
            t = tout
            istate = 2
            rwork(11:13) = (/ hu, h, tn /)
            iwork(11:15) = (/ nst, nfe, nje, nqu, nq /)
            iwork(19:21) = (/ nnz, ngp, nlu /)
            iwork(25:26) = (/ nzl, nzu /)
            return
      
      elif itask == 3 :
         tp = tn - hu*(1.0 + 100.0*uround)
         if (tp - tout)*h > 0.0 :
            write(*,*)'***ERROR lsodes: itask = ', itask, ' and tout (=', tout, ') behind tcur - hu (= ', tp, ').'
            sys.exit()
         if (tn - tout)*h >= 0.0 :
            do i = 1,n
               y(i) = rwork(i+lyh-1)
            t = tn
            istate = 2
            rwork(11:13) = (/ hu, h, tn /)
            iwork(11:15) = (/ nst, nfe, nje, nqu, nq /)
            iwork(19:21) = (/ nnz, ngp, nlu /)
            iwork(25:26) = (/ nzl, nzu /)
            return
      
      elif itask == 4 :
         tcrit = rwork(1)
         if (tn - tcrit)*h > 0.0 :
            write(*,*)'***ERROR lsodes: itask = 4 or 5 and tcrit (=', tcrit, ') is behind tcur (', tn, ').'
            sys.exit()
         if (tcrit - tout)*h < 0.0 :
            write(*,*)'***ERROR lsodes: itask = 4 or 5 and tcrit (=', tcrit, ') is behind tout (', tout, ').'
            sys.exit()
         if (tn - tout)*h >= 0.0 :
            call intdy (tout, 0, rwork(lyh), nyh, y, iflag)
            if iflag != 0 :
               write(*,*)'***ERROR lsodes: trouble from intdy. itask = ', itask, ', tout = ', tout, '.'
               sys.exit()
            t = tout
            istate = 2
            rwork(11:13) = (/ hu, h, tn /)
            iwork(11:15) = (/ nst, nfe, nje, nqu, nq /)
            iwork(19:21) = (/ nnz, ngp, nlu /)
            iwork(25:26) = (/ nzl, nzu /)
            return
      
      elif itask == 5 :
         tcrit = rwork(1)
         if (tn - tcrit)*h > 0.0 :
            write(*,*)'***ERROR lsodes: itask = 4 or 5 and tcrit (=', tcrit, ') is behind tcur (', tn, ').'
            sys.exit()
 
      if itask == 4 or itask == 5 :
         hmx = abs(tn) + abs(h)
         ihit = abs(tn - tcrit) <= 100.0*uround*hmx
         if ihit :
            do i = 1,n
               y(i) = rwork(i+lyh-1)
            t = tn
            if ihit : t = tcrit
            istate = 2
            rwork(11:13) = (/ hu, h, tn /)
            iwork(11:15) = (/ nst, nfe, nje, nqu, nq /)
            iwork(19:21) = (/ nnz, ngp, nlu /)
            iwork(25:26) = (/ nzl, nzu /)
            return
         tnext = tn + h*(1.0 + 4.0*uround)
         if (tnext - tcrit)*h > 0.0 :
            h = (tcrit - tn)*(1.0 - 4.0*uround)
            if istate == 2 : jstart = -2
 
#  the next block is normally executed for all calls and contains the call to the one-step core integrator stode.
#  this is a looping point for the integration steps.
#  first check for too many steps being taken, update ewt (if not at start of problem), check for too much accuracy being requested, and
#  check for h below the roundoff level in t.
   while True :
      if (nst-nslast) >= mxstep :
#        the maximum number of steps was taken before reaching tout.
         write(*,*)'***ERROR lsodes: at current t (=', mxstep, '), mxstep (=', tn, ') steps taken on this call before reaching tout.'
         sys.exit()
      call ewset (n, itol, rtol, atol, rwork(lyh), rwork(lewt))
      do i = 1,n
         if rwork(i+lewt-1) <= 0.0 :
#           ewt(i) <= 0.0 for some i (not at start of problem).
            ewti = rwork(lewt+i-1)
            write(*,*)'***ERROR lsodes: at t (=', tn, '), ewt(', i, ') has become non-positive (=', ewti, ').'
            sys.exit()
         rwork(i+lewt-1) = 1.0/rwork(i+lewt-1)
      tolsf = uround*vnorm(n, rwork(lyh), rwork(lewt))
      if tolsf > 1.0 :
         write(*,*)'***ERROR lsodes: at t (=', tn, '), too much accuracy requested for precision of machine. see tolsf (=', tolsf, ').'
         sys.exit()
      
      if (tn + h) == tn :
         write(*,*)'***ERROR lsodes: internal t (', tn, ') and h (=', h, ') are such that in the machine, t + h = t on the next step.'
         sys.exit()
 
      call stode (neq, y, rwork(lyh), nyh, rwork(lyh), rwork(lewt), rwork(lsavf), rwork(lacor), rwork(lwm), rwork(lwm))
 
      kgo = 1 - kflag
      
#     the following block handles the case of a successful return from the
#     core integrator (kflag = 0).  test for stop conditions.
      if kgo == 1 :
         init = 1
 
      elif kgo == 2 :
#        kflag = -1.  error test failed repeatedly or with abs(h) = hmin.
         write(*,*)'***ERROR lsodes: at t(=', tn, ') and step size h(=', h, '), the error test failed repeatedly or with abs(h) = hmin.'
         sys.exit()
 
      elif kgo == 3 :
#        kflag = -2.  convergence failed repeatedly or with abs(h) = hmin.
         write(*,*)'***ERROR lsodes: at t(=', tn, ') and step size h(=', h, '), the corrector convergence failed repeatedly or with abs(h) = hmin.'
         sys.exit()
 
      elif kgo == 4 :
#        kflag = -3.  fatal error flag returned by prjs or slss (cdrv).
         write(*,*)'***ERROR lsodes: at t(=', tn, ') and step size h(=', h, '), a fatal error flag was returned by cdrv (by way of subroutine prjs or slss).'
         sys.exit()
      
      if itask == 1 :
#        if tout has been reached, interpolate.
         if (tn - tout)*h >= 0.0 :
            call intdy(tout, 0, rwork(lyh), nyh, y, iflag)
            t = tout
            istate = 2
            rwork(11:13) = (/ hu, h, tn /)
            iwork(11:15) = (/ nst, nfe, nje, nqu, nq /)
            iwork(19:21) = (/ nnz, ngp, nlu /)
            iwork(25:26) = (/ nzl, nzu /)
            return
      
      elif itask == 2 :
#        jump to exit.
         do i = 1,n
            y(i) = rwork(i+lyh-1)
         t = tn
         istate = 2
         rwork(11:13) = (/ hu, h, tn /)
         iwork(11:15) = (/ nst, nfe, nje, nqu, nq /)
         iwork(19:21) = (/ nnz, ngp, nlu /)
         iwork(25:26) = (/ nzl, nzu /)
         return
      
      elif itask == 3 :
#        jump to exit if tout was reached.
         if (tn - tout)*h >= 0.0 :
            do i = 1,n
               y(i) = rwork(i+lyh-1)
            t = tn
            istate = 2
            rwork(11:13) = (/ hu, h, tn /)
            iwork(11:15) = (/ nst, nfe, nje, nqu, nq /)
            iwork(19:21) = (/ nnz, ngp, nlu /)
            iwork(25:26) = (/ nzl, nzu /)
            return
      
      elif itask == 4 :
#        see if tout or tcrit was reached. adjust h if necessary.
         if (tn - tout)*h >= 0.0 :
            call intdy (tout, 0, rwork(lyh), nyh, y, iflag)
            t = tout
            istate = 2
            rwork(11:13) = (/ hu, h, tn /)
            iwork(11:15) = (/ nst, nfe, nje, nqu, nq /)
            iwork(19:21) = (/ nnz, ngp, nlu /)
            iwork(25:26) = (/ nzl, nzu /)
            return
         else :
            hmx = abs(tn) + abs(h)
            ihit = abs(tn - tcrit) <= 100.0*uround*hmx
            if ihit :
               do i = 1,n
                  y(i) = rwork(i+lyh-1)
               t = tn
               if ihit : t = tcrit
               istate = 2
               rwork(11:13) = (/ hu, h, tn /)
               iwork(11:15) = (/ nst, nfe, nje, nqu, nq /)
               iwork(19:21) = (/ nnz, ngp, nlu /)
               iwork(25:26) = (/ nzl, nzu /)
               return
            else :
               tnext = tn + h*(1.0 + 4.0*uround)
               if (tnext - tcrit)*h > 0.0 :
                  h = (tcrit - tn)*(1.0 - 4.0*uround)
                  jstart = -2
      
      elif itask == 5 :
#        see if tcrit was reached and jump to exit.
         hmx = abs(tn) + abs(h)
         ihit = abs(tn - tcrit) <= 100.0*uround*hmx
         do i = 1,n
            y(i) = rwork(i+lyh-1)
         t = tn
         if ihit : t = tcrit
         istate = 2
         rwork(11:13) = (/ hu, h, tn /)
         iwork(11:15) = (/ nst, nfe, nje, nqu, nq /)
         iwork(19:21) = (/ nnz, ngp, nlu /)
         iwork(25:26) = (/ nzl, nzu /)
         return

#--------------------------------------------------------------------------------------------------
# driver for subroutines for solving sparse nonsymmetric systems of linear equations (compressed pointer storage)
#
# parameters
# class abbreviations are
#    n - integer variable
#    f - real variable
#    v - supplies a value to the driver
#    r - returns a result from the driver
#    i - used internally by the driver
#    a - array
#
# the nonzero entries of the coefficient matrix m are stored row-by-row in the array a.  to identify the individual nonzero entries in each row, 
# we need to know in which column each entry lies.
# the column indices which correspond to the nonzero entries of m are stored in the array ja.  i.e., if  a(k) = m(i,j), then ja(k) = j.  
# in addition, we need to know where each row starts and how long it is.
# the index positions in ja and a where the rows of m begin are stored in the array ia.  
# i.e., if m(i,j) is the first nonzero entry (stored) in the i-th row and a(k) = m(i,j), then# ia(i) = k.
# moreover, the index in ja and a of the first location following the last element in the last row is stored in ia(n+1).
# thus, the number of entries in the i-th row is given by ia(i+1) - ia(i), the nonzero entries of the i-th row are stored consecutively in
# a(ia(i)), a(ia(i)+1), ..., a(ia(i+1)-1), and the corresponding column indices are stored consecutively in ja(ia(i)), ja(ia(i)+1), ..., ja(ia(i+1)-1).
# for example, the 5 by 5 matrix
#             ( 1. 0. 2. 0. 0.)
#             ( 0. 3. 0. 0. 0.)
#         m = ( 0. 4. 5. 6. 0.)
#             ( 0. 0. 0. 7. 0.)
#             ( 0. 0. 0. 8. 9.)
# would be stored as
#         ia - 1  3  4  7  8 10
#         ja - 1  3  2  2  3  4  4  4  5
#          a - 1. 2. 3. 4. 5. 6. 7. 8. 9.         .
#
# nv    - n     - number of variables/equations.
# fva   - a     - nonzero entries of the coefficient matrix m, stored by rows. size = number of nonzero entries in m.
# nva   - ia    - pointers to delimit the rows in a. size = n+1.
# nva   - ja    - column numbers corresponding to the elements of a. size = size of a.
# fva   - b     - right-hand side b. b and z can the same array. size = n.
# fra   - z     - solution x.  b and z can be the same array. size = n.
#
# the rows and columns of the original matrix m can be reordered (e.g., to reduce fillin or ensure numerical stability) before calling the driver.
# if no reordering is done, then set r(i) = c(i) = ic(i) = i for i=1,...,n.  the solution z is returned in the original order.
# if the columns have been reordered (i.e., c(i)!=i  for some i), then the driver will call a subroutine (nroc) which rearranges
# each row of ja and a, leaving the rows in the original order, but placing the elements of each row in increasing order with respect to the new ordering.
# if  path!=1, then nroc is assumed to have been called already.
#
# nva   - r     - ordering of the rows of m. size = n.
# nva   - c     - ordering of the columns of m. size = n.
# nva   - ic    - inverse of the ordering of the columns of m.  i.e., ic(c(i)) = i  for i=1,...,n. size = n.
#
# the solution of the system of linear equations is divided into three stages:
# nsfc -- the matrix m is processed symbolically to determine where filling will occur during the numeric factorization.
# nnfc -- the matrix m is factored numerically into the product ldu of a unit lower triangular matrix l, a diagonal matrix d, and 
#         a unit upper triangular matrix u, and the system mx = b  is solved.
# nnsc -- the linear system  mx = b  is solved using the ldu or factorization from nnfc.
# nntc -- the transposed linear system  mt x = b  is solved using the ldu factorization from nnf.
#
# for several systems whose coefficient matrices have the same nonzero structure, nsfc need be done only once (for the first system).
# then nnfc is done once for each additional system. for several systems with the same coefficient matrix, nsfc and nnfc need be done only once (for the first system).
# then nnsc or nntc is done once for each additional right-hand side.
#
# nv    - path  - path specification.  values and their meanings are:
#    1  perform nroc, nsfc, and nnfc.
#    2  perform nnfc only  (nsfc is assumed to have been done in a manner compatible with the storage allocation used in the driver).
#    3  perform nnsc only  (nsfc and nnfc are assumed to have been done in a manner compatible with the storage allocation used in the driver).
#    4  perform nntc only  (nsfc and nnfc are assumed to have been done in a manner compatible with the storage allocation used in the driver).
#    5  perform nroc and nsfc.
#
# various errors are detected by the driver and the individual subroutines.
#
# nr  error flag.  values and their meanings are:
#         0     no errors detected
#         n+k   null row in a  --  row = k
#        2n+k   duplicate entry in a  --  row = k
#        3n+k   insufficient storage in nsfc  --  row = k
#        4n+1   insufficient storage in nnfc
#        5n+k   null pivot  --  row = k
#        6n+k   insufficient storage in nsfc  --  row = k
#        7n+1   insufficient storage in nnfc
#        8n+k   zero pivot  --  row = k
#       10n+1   insufficient storage in cdrv
#       11n+1   illegal path specification
#
# working storage is needed for the factored form of the matrix m plus various temporary vectors.
# the arrays isp and rsp should be equivalenced.  integer storage is allocated from the beginning of isp and real storage from the end of rsp.
#
# nv    - nsp   - declared dimension of rsp.  nsp generally must be larger than  8n+2 + 2k  (where  k = (number of nonzero entries in m)).
# nvira - isp   - integer working storage divided up into various arrays needed by the subroutines.  isp and rsp should be equivalenced. size = lratio*nsp.
# fvira - rsp   - real working storage divided up into various arrays needed by the subroutines.  isp and rsp should be equivalenced. size = nsp.
# nr    - esp   - if sufficient storage was available to perform the symbolic factorization (nsfc), then esp is set to the amount of excess storage provided 
#                 (negative if insufficient storage was available to perform the numeric factorization (nnfc)).
#
#
#--------------------------------------------------------------------------------------------------
def cdrv(n, r, c, ic, ia, ja, a, b, z, nsp, isp, rsp, esp, path, flag)

   integer r(1), c(1), ic(1), ia(1), ja(1), isp(1), esp, path, flag, d, u, q, row, tmp, ar, umax
   double precision a(1), b(1), z(1), rsp(1)
 
   if path < 1 or path > 5 :
      flag = 10*n + 1
      return
#  initialize and divide up temporary storage
   il = 1
   ijl = il + (n+1)
   iu = ijl + n
   iju = iu + (n+1)
   irl = iju + n
   jrl = irl + n
   jl = jrl + n
   if (path-1) * (path-5) == 0 :
      max = (2*nsp + 1 - jl) - (n+1) - 5*n
      jlmax = max/2
      q = jl + jlmax
      ira = q + (n+1)
      jra = ira + n
      irac = jra + n
      iru = irac + n
      jru = iru + n
      jutmp = jru + n
      jumax = 2*nsp + 1 - jutmp
      esp = max/2
      if jlmax<=0 or jumax<=0 :
#        error.. insufficient storage
         flag = 10*n + 1
         return
      do i = 1,n
         if c(i)!=i :
            ar = nsp + 1 - n
            call nroc(n, ic, ia, ja, a, isp(il), rsp(ar), isp(iu), flag)
            if flag != 0 :
#              error.. in nroc, nsfc, nnfc, or nnsc
               return 
            break
      call  nsfc(n, r, ic, ia, ja, jlmax, isp(il), isp(jl), isp(ijl), jumax, isp(iu), isp(jutmp), isp(iju), isp(q), isp(ira), isp(jra), isp(irac), isp(irl), isp(jrl), isp(iru), isp(jru), flag)
      if flag != 0 :
#        error.. in nroc, nsfc, nnfc, or nnsc
         return 
      jlmax = isp(ijl+n-1)
      ju = jl + jlmax
      jumax = isp(iju+n-1)
      if jumax > 0 :
         do j = 1,jumax
            isp(ju+j-1) = isp(jutmp+j-1)
   jlmax = isp(ijl+n-1)
   ju = jl  + jlmax
   jumax = isp(iju+n-1)
   l = (ju + jumax - 2 + 2)/2 + 1
   lmax = isp(il+n) - 1
   d = l + lmax
   u = d + n
   row = nsp + 1 - n
   tmp = row - n
   umax = tmp - u
   esp = umax - (isp(iu+n) - 1)
 
   if (path-1) * (path-2) == 0 :
      if umax<0 :
#        error.. insufficient storage
         flag = 10*n + 1
         return
      call nnfc(n, r, c, ic, ia, ja, a, z, b, lmax, isp(il), isp(jl), isp(ijl), rsp(l), rsp(d), umax, isp(iu), isp(ju), isp(iju), rsp(u), rsp(row), rsp(tmp), isp(irl), isp(jrl), flag)
      if flag != 0 :
#        error.. in nroc, nsfc, nnfc, or nnsc
         return 
   if (path-3) == 0 : call nnsc(n, r, c, isp(il), isp(jl), isp(ijl), rsp(l), rsp(d), isp(iu), isp(ju), isp(iju), rsp(u), z, b, rsp(tmp))
   if (path-4) == 0 : call nntc(n, r, c, isp(il), isp(jl), isp(ijl), rsp(l), rsp(d), isp(iu), isp(ju), isp(iju), rsp(u), z, b, rsp(tmp))
   return

#--------------------------------------------------------------------------------------------------
# cfode is called by the integrator routine to set coefficients needed there.
# the coefficients for the current method, as given by the value of meth, are set for all orders and saved.
# the maximum order assumed here is 12 if meth = 1 and 5 if meth = 2.
# (a smaller value of the maximum order is also allowed.)
# cfode is called once at the beginning of the problem, and is not called again unless and until meth is changed.
#
# the elco array contains the basic method coefficients.
# the coefficients el(i), 1 <= i <= nq+1, for the method of order nq are stored in elco(i,nq).
# they are given by a genetrating polynomial, i.e., l(x) = el(1) + el(2)*x + ... + el(nq+1)*x**nq.
# for the implicit adams methods, l(x) is given by dl/dx = (x+1)*(x+2)*...*(x+nq-1)/factorial(nq-1), l(-1) = 0.
# for the bdf methods, l(x) is given by l(x) = (x+1)*(x+2)* ... *(x+nq)/k, where k = factorial(nq)*(1 + 1/2 + ... + 1/nq).
#
# the tesco array contains test constants used for the local error test and the selection of step size and/or order.
# at order nq, tesco(k,nq) is used for the selection of step size at order nq - 1 if k = 1, at order nq if k = 2, and at order nq + 1 if k = 3.
#--------------------------------------------------------------------------------------------------
def cfode(meth, elco, tesco)

   integer meth
   integer i, ib, nq, nqm1, nqp1
   double precision elco, tesco
   double precision agamq, fnq, fnqm1, pc, pint, ragq, rqfac, rq1fac, tsign, xpin
   dimension elco(13,12), tesco(3,12)
   dimension pc(12)
 
   if meth == 1 :
      elco(1,1) = 1.0
      elco(2,1) = 1.0
      tesco(1,1) = 0.0
      tesco(2,1) = 2.0
      tesco(1,2) = 1.0
      tesco(3,12) = 0.0
      pc(1) = 1.0
      rqfac = 1.0
      do nq = 2,12
#        the pc array will contain the coefficients of the polynomial p(x) = (x+1)*(x+2)*...*(x+nq-1). initially, p(x) = 1.
         rq1fac = rqfac
         rqfac = rqfac/float(nq)
         nqm1 = nq - 1
         fnqm1 = float(nqm1)
         nqp1 = nq + 1
#        form coefficients of p(x)*(x+nq-1).
         pc(nq) = 0.0
         do ib = 1,nqm1
            i = nqp1 - ib
            pc(i) = pc(i-1) + fnqm1*pc(i)
         pc(1) = fnqm1*pc(1)
#        compute integral, -1 to 0, of p(x) and x*p(x).
         pint = pc(1)
         xpin = pc(1)/2.0
         tsign = 1.0
         do i = 2,nq
            tsign = -tsign
            pint = pint + tsign*pc(i)/float(i)
            xpin = xpin + tsign*pc(i)/float(i+1)
#        store coefficients in elco and tesco.
         elco(1,nq) = pint*rq1fac
         elco(2,nq) = 1.0
         do i = 2,nq
            elco(i+1,nq) = rq1fac*pc(i)/float(i)
         agamq = rqfac*xpin
         ragq = 1.0/agamq
         tesco(2,nq) = ragq
         if nq < 12 : tesco(1,nqp1) = ragq*rqfac/float(nqp1)
         tesco(3,nqm1) = ragq
 
   else : # meth == 2
      pc(1) = 1.0
      rq1fac = 1.0
      do nq = 1,5
#        the pc array will contain the coefficients of the polynomial p(x) = (x+1)*(x+2)*...*(x+nq). initially, p(x) = 1.
         fnq = float(nq)
         nqp1 = nq + 1
#        form coefficients of p(x)*(x+nq).
         pc(nqp1) = 0.0
         do ib = 1,nq
            i = nq + 2 - ib
            pc(i) = pc(i-1) + fnq*pc(i)
         pc(1) = fnq*pc(1)
#        store coefficients in elco and tesco.
         do i = 1,nqp1
            elco(i,nq) = pc(i)/pc(2)
         elco(2,nq) = 1.0
         tesco(1,nq) = rq1fac
         tesco(2,nq) = float(nqp1)/elco(1,nq)
         tesco(3,nq) = float(nq+2)/elco(1,nq)
         rq1fac = rq1fac/fnq
   return

#--------------------------------------------------------------------------------------------------
# this routine counts the number of nonzero elements in the strict upper triangle of the matrix m + m(transpose), 
# where the sparsity structure of m is given by pointer arrays ia and ja.
# this is needed to compute the storage requirements for the sparse matrix reordering operation in odrv.

# if matrix elements are held in a(k), k=1,...,nz,
# ja(k) holds the column number of the element held in a(k).
# ia(i) contains the address of the first element of row i and ia(n+1) contains the address of the first unused element in a.
#--------------------------------------------------------------------------------------------------
def cntnzu(n, ia, ja, nzsut)

   integer n, nzsut
   integer ia(n), ja(n)
   integer ii, jj, j, jmin, jmax, k, kmin, kmax, num
   logical flag
 
#  counter of nonzeros
   num = 0
#  loop over columns
   do ii = 1,n
#     addresses of the first and last elements of column ii
      jmin = ia(ii) 
      jmax = ia(ii+1) - 1
      if jmin <= jmax :
         do j = jmin,jmax
            if ja(j) < ii :
               jj =ja(j)
               kmin = ia(jj)
               kmax = ia(jj+1) - 1
               if kmin <= kmax :
                  flag = False
                  do k = kmin,kmax
                     if ja(k) == ii :
                        flag = True
                        break
               if not flag : num = num + 1
            elif ja(j) > ii :
               num = num + 1
   nzsut = num
   write(*,*)'nzsut ', nzsut
   return

#--------------------------------------------------------------------------------------------------
# this subroutine sets the error weight vector ewt according to ewt(i) = rtol(i)*abs(ycur(i)) + atol(i), i = 1,...,n,
# with the subscript on rtol and/or atol possibly replaced by 1 above, depending on the value of itol.
#--------------------------------------------------------------------------------------------------
def ewset(n, itol, rtol, atol, ycur, ewt)

   integer n, itol
   integer i
   double precision rtol, atol, ycur, ewt
   dimension rtol(1), atol(1), ycur(n), ewt(n)
 
   if itol == 1 :
      do i = 1,n
         ewt(i) = rtol(1)*abs(ycur(i)) + atol(1)
   elif itol == 2 :
      do i = 1,n
         ewt(i) = rtol(1)*abs(ycur(i)) + atol(i)
   elif itol == 3 :
      do i = 1,n
         ewt(i) = rtol(i)*abs(ycur(i)) + atol(1)
   elif itol == 4 :
       do i = 1,n
          ewt(i) = rtol(i)*abs(ycur(i)) + atol(i)
   return

#--------------------------------------------------------------------------------------------------
# intdy computes interpolated values of the k-th derivative of the dependent variable vector y, and stores it in dky.
# this routine is called within the package with k = 0 and t = tout, but may also be called by the user for any k up to the current order.
# (see detailed instructions in the usage documentation.)
#
# the computed values in dky are gotten by interpolation using the nordsieck history array yh.
# this array corresponds uniquely to a vector-valued polynomial of degree nqcur or less, 
# and dky is set to the k-th derivative of this polynomial at t.
# the formula for dky is..
#              q
#  dky(i)  =  sum  c(j,k) * (t - tn)**(j-k) * h**(-j) * yh(i,j+1)
#             j=k
# where  c(j,k) = j*(j-1)*...*(j-k+1), q = nqcur, tn = tcur, h = hcur.
# the quantities  nq = nqcur, l = nq+1, n = neq, tn, and h are communicated by common.  the above sum is done in reverse order.
# iflag is returned negative if either k or t is out of bounds.
#--------------------------------------------------------------------------------------------------
def intdy(t, k, yh, nnyh, dky, iflag)

   integer k, nnyh, iflag
   integer i, ic, j, jb, jb2, jj, jj1, jp1
   double precision t, yh, dky
   double precision c, r, s, tp
   dimension yh(nnyh,1), dky(1)
 
   double precision crate, el, elco, hold, rmax, tesco, 
  +   ccmax, el0, h, hmin, hmxi, hu, rc, tn, uround
   integer init, lyh, lewt, lacor, lsavf, lwm, liwm, mxstep, mxhnil, nhnil, ntrep, nslast, nyh,
  +   ialth, ipup, lmax, meo, nqnyh, nslp, 
  +   icf, ierpj, iersl, jcur, jstart, kflag, meth, miter, 
  +   maxord, maxcor, n, nq, nst, nfe, nje, nqu
   common /ls0001/ crate, el(13), elco(13,12), hold, rmax, tesco(3,12),
  +   ccmax, el0, h, hmin, hmxi, hu, rc, tn, uround,
  +   init, lyh, lewt, lacor, lsavf, lwm, liwm, mxstep, mxhnil, nhnil, ntrep, nslast, nyh, 
  +   ialth, ipup, lmax, meo, nqnyh, nslp, 
  +   icf, ierpj, iersl, jcur, jstart, kflag, meth, miter,
  +   maxord, maxcor, n, nq, nst, nfe, nje, nqu
 
   iflag = 0
   if k < 0 or k > nq :
      write(*,*)'***ERROR lsodes/intdy: k (=', k, ') is illegal'
      iflag = -1
      return
   tp = tn - hu -  100.0*uround*(tn + hu)
   if (t-tp)*(t-tn) > 0.0 :
      write(*,*)'***ERROR lsodes/intdy: t (=', t, ') is illegal, i.e. not in interval tcur - hu (=', tp, ') to tcur (=', tn, ')'
      iflag = -2
      return
 
   s = (t - tn)/h      
   ic = 1
   if k != 0 :
      jj1 = nq + 1 - k
      do jj = jj1,nq
         ic = ic*jj
   c = float(ic)
   do i = 1,n
      dky(i) = c*yh(i,nq+1)
   if k != nq :
      jb2 = nq - k
      do jb = 1,jb2
         j = nq - jb
         jp1 = j + 1
         ic = 1
         if k != 0 :
            jj1 = jp1 - k
            do jj = jj1,j
               ic = ic*jj
         c = float(ic)
         do i = 1,n
            dky(i) = c*yh(i,jp1) + s*dky(i)
      if k == 0 : return
   r = h**(-k)
   do i = 1,n
      dky(i) = r*dky(i)
   return

#--------------------------------------------------------------------------------------------------
# this routine serves as an interface between the driver and subroutine prep. it is called only if miter is 1 or 2.
# tasks performed here are..
#  * call prep,
#  * reset the required wm segment length lenwk,
#  * move yh back to its final location (following wm in rwork),
#  * reset pointers for yh, savf, ewt, and acor, and
#  * move ewt to its new position if istate = 1.
# ipflag is an output error indication flag.  ipflag = 0 if there was no trouble, and ipflag is the value of the prep error flag ipper
# if there was trouble in subroutine prep.
#--------------------------------------------------------------------------------------------------
def iprep(neq, y, rwork, ia, ja, ipflag)

   integer neq, ia, ja, ipflag
   integer i, imax, lewtn, lyhd, lyhn
   double precision y, rwork
   dimension y(1), rwork(1), ia(1), ja(1)
 
   double precision crate, el, elco, hold, rmax, tesco, 
  +   ccmax, el0, h, hmin, hmxi, hu, rc, tn, uround
   integer init, lyh, lewt, lacor, lsavf, lwm, liwm, mxstep, mxhnil, nhnil, ntrep, nslast, nyh,
  +   ialth, ipup, lmax, meo, nqnyh, nslp, 
  +   icf, ierpj, iersl, jcur, jstart, kflag, meth, miter, 
  +   maxord, maxcor, n, nq, nst, nfe, nje, nqu
   common /ls0001/ crate, el(13), elco(13,12), hold, rmax, tesco(3,12),
  +   ccmax, el0, h, hmin, hmxi, hu, rc, tn, uround,
  +   init, lyh, lewt, lacor, lsavf, lwm, liwm, mxstep, mxhnil, nhnil, ntrep, nslast, nyh, 
  +   ialth, ipup, lmax, meo, nqnyh, nslp, 
  +   icf, ierpj, iersl, jcur, jstart, kflag, meth, miter,
  +   maxord, maxcor, n, nq, nst, nfe, nje, nqu
 
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
 
   ipflag = 0
#  call prep to make matrix preprocessing operations.
   call prep(neq, y, rwork(lyh), rwork(lsavf), rwork(lewt), rwork(lacor), ia, ja, rwork(lwm), rwork(lwm), ipflag)
   lenwk = max0(lreq,lwmin)
   if ipflag < 0 : return
 
#  if prep was successful, move yh to end of required space for wm.
   lyhn = lwm + lenwk
   if lyhn > lyh : return
 
   lyhd = lyh - lyhn
   if lyhd != 0 :
      imax = lyhn - 1 + lenyhm
      do i = lyhn,imax
         rwork(i) = rwork(i+lyhd)
      lyh = lyhn
 
#  reset pointers for savf, ewt, and acor.
   lsavf = lyh + lenyh
   lewtn = lsavf + n
   lacor = lewtn + n
#  if istate = 1, move ewt (left) to its new position.
   if lewtn > lewt : return
   do i = 1,n
      rwork(i+lewtn-1) = rwork(i+lewt-1)
   lewt = lewtn
   return

#--------------------------------------------------------------------------------------------------
# this subroutine constructs groupings of the column indices of the jacobian matrix, 
# used in the numerical evaluation of the jacobian by finite differences.
# input..
# n      = the order of the matrix.
# ia,ja  = sparse structure descriptors of the matrix by rows.
# maxg   = length of available storate in the igp array.
#
# output..
# ngrp   = number of groups.
# jgp    = array of length n containing the column indices by groups.
# igp    = pointer array of length ngrp + 1 to the locations in jgp of the beginning of each group.
# ier    = error indicator.  ier = 0 if no error occurred, or 1 if maxg was insufficient.
#
# incl and jdone are working arrays of length n.
#--------------------------------------------------------------------------------------------------
def jgroup(n,ia,ja,maxg,ngrp,igp,jgp,incl,jdone,ier)

   integer n, ia, ja, maxg, ngrp, igp, jgp, incl, jdone, ier
   dimension ia(1), ja(1), igp(1), jgp(n), incl(n), jdone(n)
   integer i, j, k, kmin, kmax, ncol, ng
 
   ier = 0
   do j = 1,n
      jdone(j) = 0
   ncol = 1
   do ng = 1,maxg
       igp(ng) = ncol
       do i = 1,n
          incl(i) = 0
       do j = 1,n
#         reject column j if it is already in a group.
          if jdone(j) != 1 :
             kmin = ia(j)
             kmax = ia(j+1) - 1
             do k = kmin,kmax
#               reject column j if it overlaps any column already in this group.
                i = ja(k)
                if incl(i) == 1 : break
#            accept column j into group ng.
             jgp(ncol) = j
             ncol = ncol + 1
             jdone(j) = 1
             do k = kmin,kmax
                i = ja(k)
                incl(i) = 1
#      stop if this group is empty (grouping is complete).
       if ncol == igp(ng) :
          ngrp = ng - 1
          return
#  error return if not all columns were chosen (maxg too small).
   if ncol > n :
      ng = maxg
      ngrp = ng - 1
      return
   ier = 1
   return

#--------------------------------------------------------------------------------------------------
#     initialization
#--------------------------------------------------------------------------------------------------
def mdi(n, ia, ja, max, v, l, head, last, next, mark, tag, flag)

   integer ia(n), ja(1), v(1), l(1), head(1), last(1), next(1), mark(1), tag, flag, sfs, vi, dvi, vj
 
#  initialize degrees, element lists, and degree lists
   do vi = 1,n
      mark(vi) = 1
      l(vi) = 0
      head(vi) = 0
 
   sfs = n + 1
#  create nonzero structure for each nonzero entry a(vi,vj)
   do vi = 1,n
     jmin = ia(vi)
     jmax = ia(vi+1) - 1
     if jmin<=jmax :
        do j = jmin,jmax
           vj = ja(j)
           if vj < vi :
#             if a(vi,vj) is in strict lower triangle check for previous occurrence of a(vj,vi)
              lvk = vi
              kmax = mark(vi) - 1
              if kmax != 0 :
                 do k = 1,kmax
                    lvk = l(lvk)
                    if v(lvk)==vj : break
 
           elif vj > vi :
#             for unentered entries a(vi,vj)
              if sfs>=max :
#                error - insufficient storage
                 flag = 9*n + vi
                 return
 
#             enter vj in element list for vi
              mark(vi) = mark(vi) + 1
              v(sfs) = vj
              l(sfs) = l(vi)
              l(vi) = sfs
              sfs = sfs+1
 
#             enter vi in element list for vj
              mark(vj) = mark(vj) + 1
              v(sfs) = vi
              l(sfs) = l(vj)
              l(vj) = sfs
              sfs = sfs+1
 
#  сreate degree lists and initialize mark vector
   do vi = 1,n
      dvi = mark(vi)
      next(vi) = head(dvi)
      head(dvi) = vi
      last(vi) = -dvi
      nextvi = next(vi)
      if nextvi > 0 : last(nextvi) = vi
      mark(vi) = tag
   return

#--------------------------------------------------------------------------------------------------
# form element from uneliminated neighbors of vk
#--------------------------------------------------------------------------------------------------
def mdm(vk, tail, v, l, last, next, mark)

   integer  vk, tail, v(1), l(1), last(1), next(1), mark(1), tag, s, ls, vs, es, b, lb, vb, blp, blpmax
   equivalence  (vs, es)
 
#  initialize tag and list of uneliminated neighbors
   tag = mark(vk)
   tail = vk
#  for each vertex/element vs/es in element list of vk
   ls = l(vk)
   s = ls
   while s != 0 :
      ls = l(s)
      vs = v(s)
      if next(vs)<0 :
#        if es is active element, then for each vertex vb in boundary list of element es
         lb = l(es)
         blpmax = last(es)
         do blp = 1,blpmax
            b = lb
            lb = l(b)
            vb = v(b)
#           if vb is untagged vertex, then tag and append to list of uneliminated neighbors
            if mark(vb)<tag :
               mark(vb) = tag
               l(tail) = b
               tail = b
#        mark es inactive
         mark(es) = tag
      else :
#        if vs is uneliminated vertex, then tag and append to list of uneliminated neighbors
         mark(vs) = tag
         l(tail) = s
         tail = s
      s = ls
 
#  terminate list of uneliminated neighbors
   l(tail) = 0
 
   return

#--------------------------------------------------------------------------------------------------
# purge inactive elements and make mass elimination
#--------------------------------------------------------------------------------------------------
def mdp(k, ek, tail, v, l, head, last, next, mark)

   integer ek, tail, v(1), l(1), head(1), last(1), next(1), mark(1), tag, free, li, vi, lvi, evi, s, ls, es, ilp, ilpmax
 
#  initialize tag
   tag = mark(ek)
 
#  for each vertex vi in ek
   li = ek
   ilpmax = last(ek)
   if ilpmax > 0 :
      do ilp=1,ilpmax
         i = li
         li = l(i)
         vi = v(li)
         
#        remove vi from degree list
         if last(vi) > 0 :
            next(last(vi)) = next(vi)
         elif last(vi) < 0 :
            head(-last(vi)) = next(vi)
         if next(vi) > 0 : last(next(vi)) = last(vi)
 
#        remove inactive items from element list of vi
         ls = vi
         while True :
            s = ls
            ls = l(s)
            if ls == 0 : break
            es = v(ls)
            if mark(es) >= tag :
               free = ls
               l(s) = l(ls)
               ls = s
         
#        if vi is interior vertex, then remove from list and eliminate
         lvi = l(vi)
         if lvi == 0 :
            l(i) = l(li)
            li = i
            
            k = k+1
            next(vi) = -k
            last(ek) = last(ek) - 1
         else :
#           else classify vertex vi
            if l(lvi) == 0 :
               evi = v(lvi)
               if next(evi) < 0 :
                  if mark(evi)<0 :
#                    if vi is duplicate vertex, then mark as such and adjust overlap count for corresponding element
                     last(vi) = 0
                     mark(evi) = mark(evi) - 1
                  else :
#                    else if vi is prototype vertex, then mark as such, initialize overlap count for corresponding element, 
#                    and move vi to end of boundary list
                     last(vi) = evi
                     mark(evi) = -1
                     l(tail) = li
                     tail = li
                     l(i) = l(li)
                     li = i
            else :
#              else mark vi to compute degree
               last(vi) = -ek
         
#          insert ek in element list of vi
           v(free) = ek
           l(free) = l(vi)
           l(vi) = free
#  terminate boundary list
   l(tail) = 0
   return

#--------------------------------------------------------------------------------------------------
#     md finds a minimum degree ordering of the rows and columns of a general sparse matrix m stored in (ia,ja,a) format.
#     when the structure of m is nonsymmetric, the ordering is that obtained for the symmetric matrix  m + m-transpose.
#
#     additional parameters
#
#     max  - declared dimension of the one-dimensional arrays v and l.
#            max must be at least  n+2k, where k is the number of nonzeroes in the strict upper triangle of m + m-transpose
#     v    - integer one-dimensional work array.  dimension = max
#     l    - integer one-dimensional work array.  dimension = max
#     head - integer one-dimensional work array.  dimension = n
#     last - integer one-dimensional array used to return the permutation of the rows and columns of m corresponding to the minimum degree ordering.
#            dimension = n
#     next - integer one-dimensional array used to return the inverse of the permutation returned in last.  dimension = n
#     mark - integer one-dimensional work array (may be the same as v). dimension = n
#     flag - integer error flag.  values and their meanings are -
#            0 no errors detected
#            9n+k  insufficient storage in md
#--------------------------------------------------------------------------------------------------
def md(n, ia, ja, max, v, l, head, last, next, mark, flag)

   integer ia(1), ja(1), v(1), l(1), head(1), last(1), next(1), mark(1), flag, tag, dmin, vk, ek, tail
   equivalence (vk,ek)
 
#  initialization
   tag = 0
   call mdi(n, ia,ja, max,v,l, head,last,next, mark,tag, flag)
   if flag!=0 : return
 
   k = 0
   dmin = 1
 
   while k < n :
#     search for vertex of minimum degree
      while head(dmin) <= 0 :
        dmin = dmin + 1
#     remove vertex vk of minimum degree from degree list
      vk = head(dmin)
      head(dmin) = next(vk)
      if head(dmin) > 0 : last(head(dmin)) = -dmin
 
#     number vertex vk, adjust tag, and tag vk
      k = k + 1
      next(vk) = -k
      last(ek) = dmin - 1
      tag = tag + last(ek)
      mark(vk) = tag
 
#     form element ek from uneliminated neighbors of vk
      call mdm(vk, tail, v, l, last, next, mark)
#     purge inactive elements and make mass elimination
      call mdp(k, ek, tail, v, l, head, last, next, mark)
#     update degrees of uneliminated vertices in ek
      call mdu(ek, dmin, v, l, head, last, next, mark)
 
#  generate inverse permutation from permutation
   do k=1,n
      next(k) = -next(k)
      last(next(k)) = k
   return

#--------------------------------------------------------------------------------------------------
# update degrees of uneliminated vertices in ek
#--------------------------------------------------------------------------------------------------
def mdu(ek, dmin, v, l, head, last, next, mark)

   integer ek, dmin, v(1), l(1), head(1), last(1), next(1), mark(1), tag, vi, evi, dvi, s, ves, b, vb, ilp, ilpmax, blp, blpmax
 
#  initialize tag
   tag = mark(ek) - last(ek)
 
#  for each vertex vi in ek
   i = ek
   ilpmax = last(ek)
   if ilpmax > 0 :
      do ilp=1,ilpmax
         i = l(i)
         vi = v(i)
         if last(vi) < 0 :
#           if vi neither prototype nor duplicate vertex, then merge elements to compute degree
            tag = tag + 1
            dvi = last(ek)
            
#           for each vertex/element ves in element list of vi
            s = l(vi)
            while True :
               s = l(s)
               if s == 0 :
#                 insert vi in appropriate degree list
                  next(vi) = head(dvi)
                  head(dvi) = vi
                  last(vi) = -dvi
                  if next(vi) > 0 : last(next(vi)) = vi
                  if dvi < dmin : dmin = dvi
                  break
               ves = v(s)
               if next(ves) >= 0 :
#                 if ves is uneliminated vertex, then tag and adjust degree
                  mark(ves) = tag
                  dvi = dvi + 1
               else :
#                 if ves is active element, then expand check for outmatched vertex
                  if mark(ves) < 0 :
#                    else if vi is outmatched vertex, then adjust overlaps but do not compute degree
                     last(vi) = 0
                     mark(ves) = mark(ves) - 1
                     while True :
                        s = l(s)
                        if s == 0 : break
                        ves = v(s)
                        if mark(ves) < 0 : mark(ves) = mark(ves) - 1
                     break
                  else :
#                    for each vertex vb in es
                     b = ves
                     blpmax = last(ves)
                     do blp = 1,blpmax
                        b = l(b)
                        vb = v(b)
#                       if vb is untagged, then tag and adjust degree
                        if mark(vb) < tag :
                           mark(vb) = tag
                           dvi = dvi + 1
         elif last(vi) > 0 :
#           else if vi is prototype vertex, then calculate degree by inclusion/exclusion and reset overlap count
            evi = last(vi)
            dvi = last(ek) + last(evi) + mark(evi)
            mark(evi) = 0
 
#           insert vi in appropriate degree list
            next(vi) = head(dvi)
            head(dvi) = vi
            last(vi) = -dvi
            if next(vi) > 0 : last(next(vi)) = vi
            if dvi < dmin : dmin = dvi
   return

#--------------------------------------------------------------------------------------------------
# numerical ldu-factorization of sparse nonsymmetric matrix and solution of system of linear equations (compressed pointer storage)
#--------------------------------------------------------------------------------------------------
def nnfc(n, r, c, ic, ia, ja, a, z, b, lmax, il, jl, ijl, l, d, umax, iu, ju, iju, u, row, tmp, irl, jrl, flag)

   integer rk, umax
   integer r(1), c(1), ic(1), ia(1), ja(1), il(1), jl(1), ijl(1)
   integer iu(1), ju(1), iju(1), irl(1), jrl(1), flag
   double precision a(1), l(1), d(1), u(1), z(1), b(1), row(1)
   double precision tmp(1), lki, sum, dk
 
#  initialize pointers and test storage
   if il(n+1)-1 > lmax :
#     error.. insufficient storage for l
      flag = 4*n + 1
      return
   if iu(n+1)-1 > umax :
#     error.. insufficient storage for u
      flag = 7*n + 1
      return
   do k = 1,n
      irl(k) = il(k)
      jrl(k) = 0
 
#  for each row
   do k = 1,n
#     reverse jrl and zero row where kth row of l will fill in
      row(k) = 0
      i1 = 0
      if jrl(k) != 0 :
         i = jrl(k)
         while True :
            i2 = jrl(i)
            jrl(i) = i1
            i1 = i
            row(i) = 0
            i = i2
            if i == 0 : break
#     set row to zero where u will fill in
      jmin = iju(k)
      jmax = jmin + iu(k+1) - iu(k) - 1
      if jmin <= jmax :
         do j = jmin,jmax
            row(ju(j)) = 0
#     place kth row of a in row
      rk = r(k)
      jmin = ia(rk)
      jmax = ia(rk+1) - 1
      do j = jmin,jmax
         row(ic(ja(j))) = a(j)
#     initialize sum, and link through jrl
      sum = b(rk)
      i = i1
      if i != 0 :
         while True :
#           assign the kth row of l and adjust row, sum
            lki = -row(i)
#           if l is not required, then comment out the following line
            l(irl(i)) = -lki
            sum = sum + lki * tmp(i)
            jmin = iu(i)
            jmax = iu(i+1) - 1
            if jmin <= jmax :
               mu = iju(i) - jmin
               do j = jmin,jmax
                  row(ju(mu+j)) = row(ju(mu+j)) + lki * u(j)
            i = jrl(i)
            if i == 0 : break
 
#     assign kth row of u and diagonal d, set tmp(k)
      if row(k) == 0.0 :
#        error.. zero pivot
         flag = 8*n + k
         return
      dk = 1.0 / row(k)
      d(k) = dk
      tmp(k) = sum * dk
      if k != n :
         jmin = iu(k)
         jmax = iu(k+1) - 1
         if jmin <= jmax :
            mu = iju(k) - jmin
            do j = jmin,jmax
               u(j) = row(ju(mu+j)) * dk
         
#        update irl and jrl, keeping jrl in decreasing order
         i = i1
         if i != 0 :
            while True :
               irl(i) = irl(i) + 1
               i1 = jrl(i)
               if irl(i) < il(i+1) :
                  ijlb = irl(i) - il(i) + ijl(i)
                  j = jl(ijlb)
                  while i <= jrl(j) :
                     j = jrl(j)
                  jrl(i) = jrl(j)
                  jrl(j) = i
               i = i1
               if i == 0 : break
         if irl(k) < il(k+1) :
            j = jl(ijl(k))
            jrl(k) = jrl(j)
            jrl(j) = k
 
#  solve  ux = tmp  by back substitution
   k = n
   do i = 1,n
      sum =  tmp(k)
      jmin = iu(k)
      jmax = iu(k+1) - 1
      if jmin <= jmax :
         mu = iju(k) - jmin
         do j = jmin,jmax
            sum = sum - u(j) * tmp(ju(mu+j))
      tmp(k) =  sum
      z(c(k)) =  sum
   k = k-1
   flag = 0
   return
 
#--------------------------------------------------------------------------------------------------
#     numerical solution of sparse nonsymmetric system of linear equations given ldu-factorization (compressed pointer storage)
#--------------------------------------------------------------------------------------------------
def nnsc(n, r, c, il, jl, ijl, l, d, iu, ju, iju, u, z, b, tmp)

   integer r(1), c(1), il(n), jl(1), ijl(1), iu(1), ju(1), iju(1)
   double precision l(1), d(1), u(1), b(1), z(1), tmp(1), tmpk, sum
 
#  set tmp to reordered b
   do k = 1,n
      tmp(k) = b(r(k))
#  solve  ly = b  by forward substitution
   do k = 1,n
      jmin = il(k)
      jmax = il(k+1) - 1
      tmpk = -d(k) * tmp(k)
      tmp(k) = -tmpk
      if jmin <= jmax :
         ml = ijl(k) - jmin
         do j = jmin,jmax
            tmp(jl(ml+j)) = tmp(jl(ml+j)) + tmpk * l(j)
#  solve  ux = y  by back substitution
   k = n
   do i = 1,n
      sum = -tmp(k)
      jmin = iu(k)
      jmax = iu(k+1) - 1
      if jmin <= jmax :
         mu = iju(k) - jmin
         do j = jmin,jmax
            sum = sum + u(j) * tmp(ju(mu+j))
      tmp(k) = -sum
      z(c(k)) = -sum
      k = k - 1
   return

#--------------------------------------------------------------------------------------------------
#     numeric solution of the transpose of a sparse nonsymmetric system of linear equations given lu-factorization (compressed pointer storage)
#--------------------------------------------------------------------------------------------------
def nntc(n, r, c, il, jl, ijl, l, d, iu, ju, iju, u, z, b, tmp)

   integer r(1), c(1), il(n), jl(1), ijl(1), iu(n), ju(1), iju(1)
   double precision l(1), d(1), u(1), b(1), z(1), tmp(1), tmpk,sum
 
#  set tmp to reordered b
   do k = 1,n
      tmp(k) = b(c(k))
#  solve  ut y = b  by forward substitution
   do k = 1,n
      jmin = iu(k)
      jmax = iu(k+1) - 1
      tmpk = -tmp(k)
      if jmin <= jmax :
         mu = iju(k) - jmin
         do j = jmin,jmax
            tmp(ju(mu+j)) = tmp(ju(mu+j)) + tmpk * u(j)
#  solve  lt x = y  by back substitution
   k = n
   do i = 1,n
      sum = -tmp(k)
      jmin = il(k)
      jmax = il(k+1) - 1
      if jmin <= jmax :
         ml = ijl(k) - jmin
         do j = jmin,jmax
            sum = sum + l(j) * tmp(jl(ml+j))
      tmp(k) = -sum * d(k)
      z(r(k)) = tmp(k)
      k = k - 1
   return

#--------------------------------------------------------------------------------------------------
# reorders rows of a, leaving row order unchanged
#--------------------------------------------------------------------------------------------------
def nroc (n, ic, ia, ja, a, jar, ar, p, flag)

   integer  ic(1), ia(n), ja(1), jar(1), p(1), flag
   double precision  a(1), ar(1)
 
   do k = 1,n
      jmin = ia(k)
      jmax = ia(k+1) - 1
      if jmin <= jmax :
         p(n+1) = n + 1
         do j = jmin,jmax
            newj = ic(ja(j))
            i = n + 1
            while True :
               if p(i) >= newj : break
               i = p(i)
            if p(i) == newj :
               flag = n + k
               return
            p(newj) = p(i)
            p(i) = newj
            jar(newj) = ja(j)
            ar(newj) = a(j)
         i = n + 1
         do j = jmin,jmax
            i = p(i)
            ja(j) = jar(i)
            a(j) = ar(i)
   flag = 0
   return

#--------------------------------------------------------------------------------------------------
# symbolic ldu-factorization of nonsymmetric sparse matrix (compressed pointer storage)
#--------------------------------------------------------------------------------------------------
def nsfc(n, r, ic, ia, ja, jlmax, il, jl, ijl, jumax, iu, ju, iju, q, ira, jra, irac, irl, jrl, iru, jru, flag)

   integer cend, qm, rend, rk, vj
   integer ia(1), ja(1), ira(1), jra(1), il(1), jl(1), ijl(1)
   integer iu(1), ju(1), iju(1), irl(1), jrl(1), iru(1), jru(1)
   integer r(1), ic(1), q(1), irac(1), flag
   logical exit17, exit 34
 
   np1 = n + 1
   jlmin = 1
   jlptr = 0
   il(1) = 1
   jumin = 1
   juptr = 0
   iu(1) = 1
   do k = 1,n
      irac(k) = 0
      jra(k) = 0
      jrl(k) = 0
      jru(k) = 0
 
   do k = 1,n
      rk = r(k)
      iak = ia(rk)
      if iak >= ia(rk+1) :
         flag = n + rk
         return
      jaiak = ic(ja(iak))
      if jaiak > k :
         flag = 5*n + k
         return
      jra(k) = irac(jaiak)
      irac(jaiak) = k
      ira(k) = iak
 
   do k = 1,n
      q(np1) = np1
      luk = -1
      vj = irac(k)
      if vj != 0 :
         while True :
            qm = np1
            while True :
               m = qm
               qm =  q(m)
               if qm >= vj : break
            if qm == vj :
               flag = 2*n + rk
               return
            luk = luk + 1
            q(m) = vj
            q(vj) = qm
            vj = jra(vj)
            if vj == 0 : break
      lastid = 0
      lasti = 0
      ijl(k) = jlptr
      i = k
      while True :
         i = jru(i)
         if i == 0 : break
         qm = np1
         jmin = irl(i)
         jmax = ijl(i) + il(i+1) - il(i) - 1
         long = jmax - jmin
         if long >= 0 :
            jtmp = jl(jmin)
            if jtmp != k : long = long + 1
            if jtmp == k : r(i) = -r(i)
            if lastid < long :
               lasti = i
               lastid = long
            do j = jmin,jmax
               vj = jl(j)
               while True :
                  m = qm
                  qm = q(m)
                  if qm >= vj : break
               if qm != vj :
                  luk = luk + 1
                  q(m) = vj
                  q(vj) = qm
                  qm = vj
 
      qm = q(np1)
      if qm != k :
         flag = 5*n + k
         return
 
      exit17 = False
      if luk != 0 :
         if lastid == luk :
            irll = irl(lasti)
            ijl(k) = irll + 1
            if jl(irll) != k : ijl(k) = ijl(k) - 1
         else :
            if jlmin <= jlptr :
               qm = q(qm)
               do j = jlmin,jlptr
                  if jl(j) >= qm :
                     if jl(j) == qm :
                        ijl(k) = j
                        do i = j,jlptr
                           if jl(i) == qm :
                              qm = q(qm)
                              if qm > n :
                                 exit17 = True
                                 break
                           else :
                              jlmin = jlptr + 1
                              ijl(k) = jlmin
                              if luk != 0 :
                                 jlptr = jlptr + luk
                                 if jlptr > jlmax :
                                    flag = 3*n + k
                                    return
                                 qm = q(np1)
                                 do jj = jlmin,jlptr
                                    qm = q(qm)
                                    jl(jj) = qm
                              exit17 = True
                              break
                        if not exit17 : jlptr = j - 1
                     if not exit17 :
                        jlmin = jlptr + 1
                        ijl(k) = jlmin
                        if luk != 0 :
                           jlptr = jlptr + luk
                           if jlptr > jlmax :
                              flag = 3*n + k
                              return
                           qm = q(np1)
                           do jj = jlmin,jlptr
                              qm = q(qm)
                              jl(jj) = qm
                     exit17 = True
                     break
 
            if not exit17 :
               jlmin = jlptr + 1
               ijl(k) = jlmin
               if luk != 0 :
                  jlptr = jlptr + luk
                  if jlptr > jlmax :
                     flag = 3*n + k
                     return
                  qm = q(np1)
                  do jj = jlmin,jlptr
                     qm = q(qm)
                     jl(jj) = qm
 
      irl(k) = ijl(k)
      il(k+1) = il(k) + luk     
      q(np1) = np1
      luk = -1
      rk = r(k)
      jmin = ira(k)
      jmax = ia(rk+1) - 1
      if jmin <= jmax :
         do j = jmin,jmax
            vj = ic(ja(j))
            qm = np1
            while True :
               m = qm
               qm = q(m)
               if qm >= vj : break
            if qm == vj :
               flag = 2*n + rk
               return
            luk = luk + 1
            q(m) = vj
            q(vj) = qm
      lastid = 0
      lastid = 0
      lasti = 0
      iju(k) = juptr
      i = k
      i1 = jrl(k)
 
      while True :
         i = i1
         if i == 0 : break
         i1 = jrl(i)
         qm = np1
         jmin = iru(i)
         jmax = iju(i) + iu(i+1) - iu(i) - 1
         long = jmax - jmin
         if long >= 0 :
            jtmp = ju(jmin)
            if jtmp != k :
               long = long + 1
               cend = ijl(i) + il(i+1) - il(i)
               irl(i) = irl(i) + 1
               if irl(i) < cend :
                  j = jl(irl(i))
                  jrl(i) = jrl(j)
                  jrl(j) = i
            if lastid < long :
               lasti = i
               lastid = long
            do j = jmin,jmax
               vj = ju(j)
               while True :
                  m = qm
                  qm = q(m)
                  if qm >= vj : break
               if qm != vj :
                  luk = luk + 1
                  q(m) = vj
                  q(vj) = qm
                  qm = vj
 
      if il(k+1) > il(k) :
         j = jl(irl(k))
         jrl(k) = jrl(j)
         jrl(j) = k
      qm = q(np1)
      if qm != k :
         flag = 5*n + k
         return
 
      exit34 = False
      if luk != 0 :
         if lastid == luk :
            irul = iru(lasti)
            iju(k) = irul + 1
            if ju(irul) != k : iju(k) = iju(k) - 1
         else :
            if jumin <= juptr :
               qm = q(qm)
               do j = jumin,juptr
                  if ju(j) >= qm :
                     if ju(j) == qm :
                        iju(k) = j
                        do i = j,juptr
                           if ju(i) == qm :
                              qm = q(qm)
                              if qm > n :
                                 exit34 = True
                                 break
                           else :
                              jumin = juptr + 1
                              iju(k) = jumin
                              if luk != 0 :
                                 juptr = juptr + luk
                                 if juptr > jumax :
                                    flag = 3*n + k
                                    return
                                 qm = q(np1)
                                 do jj = jumin,juptr
                                    qm = q(qm)
                                    ju(jj) = qm
                              exit34 = True
                              break
                        if not exit34 : juptr = j - 1
                     if not exit34 :
                        jumin = juptr + 1
                        iju(k) = jumin
                        if luk != 0 :
                           juptr = juptr + luk
                           if juptr > jumax :
                              flag = 3*n + k
                              return
                           qm = q(np1)
                           do jj = jumin,juptr
                              qm = q(qm)
                              ju(jj) = qm
                     exit34 = True
                     break
 
            if not exit34 :
               jumin = juptr + 1
               iju(k) = jumin
               if luk != 0 :
                  juptr = juptr + luk
                  if juptr > jumax :
                     flag = 3*n + k
                     return
                  qm = q(np1)
                  do jj = jumin,juptr
                     qm = q(qm)
                     ju(jj) = qm
 
      iru(k) = iju(k)
      iu(k+1) = iu(k) + luk
      i = k
      while True :
         i1 = jru(i)
         if r(i) >= 0 :
            rend = iju(i) + iu(i+1) - iu(i)
            if iru(i) < rend :
               j = ju(iru(i))
               jru(i) = jru(j)
               jru(j) = i
         else :
            r(i) = -r(i)
         i = i1
         if i == 0 : break
         iru(i) = iru(i) + 1
      i = irac(k)
      if i != 0 :
         while True :
            i1 = jra(i)
            ira(i) = ira(i) + 1
            if ira(i) < ia(r(i)+1) :
               irai = ira(i)
               jairai = ic(ja(irai))
               if jairai <= i :
                  jra(i) = irac(jairai)
                  irac(jairai) = i
            i = i1
            if i == 0 : break
   ijl(n) = jlptr
   iju(n) = juptr
   flag = 0
   return

#--------------------------------------------------------------------------------------------------
# driver for sparse matrix reordering routines
#
# odrv finds a minimum degree ordering of the rows and columns of a matrix m stored in (ia,ja,a) format (see below).
# for the reordered matrix, the work and storage required to perform gaussian elimination is (usually) significantly less.
#
# note.. odrv and its subordinate routines have been modified to compute orderings for general matrices, 
# not necessarily having any symmetry.
# the miminum degree ordering is computed for the structure of the symmetric matrix  m + m-transpose.
# modifications to the original odrv module have been made in the coding in subroutine mdi, and in the initial comments in
# subroutines odrv and md.
#
# if only the nonzero entries in the upper triangle of m are being stored, then odrv symmetrically reorders (ia,ja,a),
# (optionally) with the diagonal entries placed first in each row.
# this is to ensure that if m(i,j) will be in the upper triangle of m with respect to the new ordering,
# then m(i,j) is stored in row i (and thus m(j,i) is not stored), whereas if m(i,j) will be in the strict lower triangle of m, 
# then m(j,i) is stored in row j (and thus m(i,j) is not stored).
#
# storage of sparse matrices.  the nonzero entries of the matrix m are stored row-by-row in the array a.
# to identify the individual nonzero entries in each row, we need to know in which column each entry lies.
# these column indices are stored in the array ja. i.e., if  a(k) = m(i,j), then ja(k) = j.
# to identify the individual rows, we need to know where each row starts.
# these row pointers are stored in the array ia. i.e.,
# if m(i,j) is the first nonzero entry (stored) in the i-th row and a(k) = m(i,j), then ia(i) = k.
# moreover, ia(n+1) points to the first location following the last element in the last row.
# thus, the number of entries in the i-th row is  ia(i+1) - ia(i), the nonzero entries in the i-th row are stored consecutively in
# a(ia(i)), a(ia(i)+1), ..., a(ia(i+1)-1), and the corresponding column indices are stored consecutively in
# ja(ia(i)), ja(ia(i)+1), ..., ja(ia(i+1)-1).
#
# since the coefficient matrix is symmetric, only the nonzero entries in the upper triangle need be stored.  for example, the matrix
#
#     ( 1  0  2  3  0 )
#     ( 0  4  0  0  0 )
# m = ( 2  0  5  6  0 )
#     ( 3  0  6  7  8 )
#     ( 0  0  0  8  9 )
#
# could be stored as
#
#    - 1  2  3  4  5  6  7  8  9 10 11 12 13
# ---+--------------------------------------
# ia - 1  4  5  8 12 14
# ja - 1  3  4  2  1  3  4  1  3  4  5  4  5
#  a - 1  2  3  4  2  5  6  3  6  7  8  8  9
#
# or (symmetrically) as
#
#    - 1  2  3  4  5  6  7  8  9
# ---+--------------------------
# ia - 1  4  5  7  9 10
# ja - 1  3  4  2  3  4  4  5  5
#  a - 1  2  3  4  5  6  7  8  9          .
#
#  parameters:
# n    - order of the matrix
#
# ia   - integer one-dimensional array containing pointers to delimit rows in ja and a.  dimension = n+1
#
# ja   - integer one-dimensional array containing the column indices corresponding to the elements of a.
#        dimension = number of nonzero entries in (the upper triangle of) m
#
# a    - real one-dimensional array containing the nonzero entries in (the upper triangle of) m, stored by rows.
#        dimension = number of nonzero entries in (the upper triangle of) m
#
# p    - integer one-dimensional array used to return the permutation of the rows and columns of m corresponding to the minimum degree ordering.  
#        dimension = n
#
# ip   - integer one-dimensional array used to return the inverse of the permutation returned in p.  dimension = n
#
# nsp  - declared dimension of the one-dimensional array isp.  nsp must be at least  3n+4k,
#        where k is the number of nonzeroes in the strict upper triangle of m
#
# isp  - integer one-dimensional array used for working storage.
#        dimension = nsp
#
# path - integer path specification.  values and their meanings are -
# 1  find minimum degree ordering only
# 2  find minimum degree ordering and reorder symmetrically stored matrix (used when only the nonzero entries in
#    the upper triangle of m are being stored)
# 3  reorder symmetrically stored matrix as specified by input permutation (used when an ordering has already
#    been determined and only the nonzero entries in the upper triangle of m are being stored)
# 4  same as 2 but put diagonal entries at start of each row
# 5  same as 3 but put diagonal entries at start of each row
#
# flag - integer error flag.  values and their meanings are:
# 0 no errors detected
# 9n+k insufficient storage in md
# 10n+1 insufficient storage in odrv
# 11n+1 illegal path specification
#--------------------------------------------------------------------------------------------------
def odrv(n, ia, ja, a, p, ip, nsp, isp, path, flag)

   integer ia(1), ja(1), p(1), ip(1), isp(1), path, flag, v, l, head, tmp, q
   double precision a(1)
   logical dflag
 
#  initialize error flag and validate path specification
   flag = 0
   if path < 1 or path > 5 :
#     error -- illegal path specified
      flag = 11*n + 1
      return
 
#  allocate storage and find minimum degree ordering
   if (path-1) * (path-2) * (path-4) == 0 :
      max = (nsp-n)/2
      v = 1
      l = v +  max
      head = l +  max
      next = head +  n
      if max < n :
         flag = 10*n + 1
         return
      call md(n, ia, ja, max, isp(v), isp(l), isp(head), p, ip, isp(v), flag)
#     error detected in md
      if flag != 0 : return
   if (path-2) * (path-3) * (path-4) * (path-5) == 0 :
      tmp = (nsp+1) - n
      q  = tmp - (ia(n+1)-1)
      if q < 1 :
#        error -- insufficient storage
         flag = 10*n + 1
         return
      dflag = path==4 or path==5
      call sro(n, ip, ia, ja, a, isp(tmp), isp(q), dflag)
   return

#--------------------------------------------------------------------------------------------------
# this routine performs preprocessing related to the sparse linear systems that must be solved if miter = 1 or 2.
# the operations that are performed here are..
#  * compute sparseness structure of jacobian according to moss,
#  * compute grouping of column indices (miter = 2),
#  * compute a new ordering of rows and columns of the matrix,
#  * reorder ja corresponding to the new ordering,
#  * perform a symbolic lu factorization of the matrix, and
#  * set pointers for segments of the iwk/wk array.
# in addition to variables described previously, prep uses the following for communication..
# yh     = the history array.  only the first column, containing the current y vector, is used. used only if moss != 0.
# savf   = a work array of length neq, used only if moss != 0.
# ewt    = array of length neq containing (inverted) error weights. used only if moss = 2 or if istate = moss = 1.
# ftem   = a work array of length neq, identical to acor in the driver, used only if moss = 2.
# wk     = a real work array of length lenwk, identical to wm in the driver.
# iwk    = integer work array, assumed to occupy the same space as wk.
# lenwk  = the length of the work arrays wk and iwk.
# istatc = a copy of the driver input argument istate (= 1 on the first call, = 3 on a continuation call).
# iys    = flag value from odrv or cdrv.
# ipper  = output error flag with the following values and meanings..
#          0  no error.
#         -1  insufficient storage for internal structure pointers.
#         -2  insufficient storage for jgroup.
#         -3  insufficient storage for odrv.
#         -4  other error flag from odrv (should never occur).
#         -5  insufficient storage for cdrv.
#         -6  other error flag from cdrv.
# moss   = the method to be used to obtain the sparsity structure of the jacobian matrix if miter = 1 or 2..
#          0  means the user has supplied ia and ja (see descriptions under iwork in the lsodes).
#          1  means the user has supplied jac and the structure will be obtained from neq initial calls to jac.
#          2  means the structure will be obtained from neq+1 initial calls to f.
# miter  = the corrector iteration method..
#          0  means functional iteration (no jacobian matrix is involved).
#          1  means chord iteration with a user-supplied sparse jacobian, given by subroutine jac.
#          2  means chord iteration with an internally generated (difference quotient) sparse jacobian
#             (using ngp extra calls to f per df/dy value, where ngp is an optional output described below.)
#          3  means chord iteration with an internally generated diagonal jacobian approximation. (using 1 extra call to f per df/dy evaluation).
#--------------------------------------------------------------------------------------------------
def prep (neq, y, yh, savf, ewt, ftem, ia, ja, wk, iwk, ipper)

   integer neq, ia, ja, iwk, ipper
   integer i, ibr, ier, ipil, ipiu, iptt1, iptt2, j, jfound, k, knew, kmax, kmin, ldif, lenigp, liwk, maxg, np1, nzsut
   double precision y, yh, savf, ewt, ftem, wk
   double precision dq, dyj, erwt, fac, yj
   dimension y(1), yh(1), savf(1), ewt(1), ftem(1), ia(1), ja(1), wk(1), iwk(1)
 
   double precision crate, el, elco, hold, rmax, tesco, 
  +   ccmax, el0, h, hmin, hmxi, hu, rc, tn, uround
   integer init, lyh, lewt, lacor, lsavf, lwm, liwm, mxstep, mxhnil, nhnil, ntrep, nslast, nyh,
  +   ialth, ipup, lmax, meo, nqnyh, nslp, 
  +   icf, ierpj, iersl, jcur, jstart, kflag, meth, miter, 
  +   maxord, maxcor, n, nq, nst, nfe, nje, nqu
   common /ls0001/ crate, el(13), elco(13,12), hold, rmax, tesco(3,12),
  +   ccmax, el0, h, hmin, hmxi, hu, rc, tn, uround,
  +   init, lyh, lewt, lacor, lsavf, lwm, liwm, mxstep, mxhnil, nhnil, ntrep, nslast, nyh, 
  +   ialth, ipup, lmax, meo, nqnyh, nslp, 
  +   icf, ierpj, iersl, jcur, jstart, kflag, meth, miter,
  +   maxord, maxcor, n, nq, nst, nfe, nje, nqu
 
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
   if ipjan+n-1 > liwk :
      ipper = -1
      lreq = 2 + (2*n + 1)/2
      lreq = max0(lenwk+1,lreq)
      return
 
   if moss == 0 :
#     moss = 0. process user-s ia,ja. add diagonal entries if necessary.
      knew = ipjan
      kmin = ia(1)
      iwk(ipian) = 1
      do j = 1,n
         jfound = 0
         kmax = ia(j+1) - 1
         if kmin <= kmax :
            do k = kmin,kmax
               i = ja(k)
               if i == j : jfound = 1
               if knew > liwk :
                  ipper = -1
                  lreq = 2 + (2*n + 1)/2
                  lreq = max0(lenwk+1,lreq)
                  return
               iwk(knew) = i
               knew = knew + 1
         if jfound == 0 :
            if knew > liwk :
               ipper = -1
               lreq = 2 + (2*n + 1)/2
               lreq = max0(lenwk+1,lreq)
               return
            iwk(knew) = j
            knew = knew + 1
         iwk(ipian+j) = knew + 1 - ipjan
         kmin = kmax + 1
 
   else : # moss != 0
 
      if istatc == 3 :
#        istate = 3 and moss != 0. load y from yh(*,1).
         do i = 1,n
            y(i) = yh(i)
      else :
#        istate = 1 and moss != 0. perturb y for structure determination.
         do i = 1,n
           erwt = 1.0/ewt(i)
           fac = 1.0 + 1.0/(float(i)+1.0)
           y(i) = y(i) + fac*math.copysign(erwt,y(i))
   
   if moss == 1 :
      continue
#     a dummy call to rhs allows user to create temporaries for use in jac.
      call rhs(neq, tn, y, savf)
      k = ipjan
      iwk(ipian) = 1
      do j = 1,n
         if k > liwk :
            ipper = -1
            lreq = 2 + (2*n + 1)/2
            lreq = max0(lenwk+1,lreq)
            return
         iwk(k) = j
         k = k + 1
         do i = 1,n
            savf(i) = 0.0
         call fjac(neq, tn, y, j, iwk(ipian), iwk(ipjan), savf)
         do i = 1,n
            if abs(savf(i)) > seth and i != j :
               if k > liwk :
                  ipper = -1
                  lreq = 2 + (2*n + 1)/2
                  lreq = max0(lenwk+1,lreq)
                  return
               iwk(k) = i
               k = k + 1
         iwk(ipian+j) = k + 1 - ipjan
 
   if moss == 2 :
#     moss = 2. compute structure from results of n + 1 calls to f.
      k = ipjan
      iwk(ipian) = 1
      call rhs(neq, tn, y, savf)
      do j = 1,n
         if k > liwk :
            ipper = -1
            lreq = 2 + (2*n + 1)/2
            lreq = max0(lenwk+1,lreq)
            return
         iwk(k) = j
         k = k + 1
         yj = y(j)
         erwt = 1.0/ewt(j)
         dyj = math.copysign(erwt,yj)
         y(j) = yj + dyj
         call rhs(neq, tn, y, ftem)
         y(j) = yj
         do i = 1,n
            dq = (ftem(i) - savf(i))/dyj
            if abs(dq) > seth and i != j :
               if k > liwk :
                  ipper = -1
                  lreq = 2 + (2*n + 1)/2
                  lreq = max0(lenwk+1,lreq)
                  return
               iwk(k) = i
               k = k + 1
         iwk(ipian+j) = k + 1 - ipjan
 
   if moss != 0 and istatc == 1 : # CHECK that istatc == 1 is correct here
#     if istate = 1 and moss != 0, restore y from yh.
      do i = 1,n
         y(i) = yh(i)
   nnz = iwk(ipian+n) - 1
   lenigp = 0
   ipigp = ipjan + nnz
   if miter == 2 :
#     compute grouping of column indices (miter = 2).
      maxg = np1
      ipjgp = ipjan + nnz
      ibjgp = ipjgp - 1
      ipigp = ipjgp + n
      iptt1 = ipigp + np1
      iptt2 = iptt1 + n
      lreq = iptt2 + n - 1
      if lreq > liwk :
         ipper = -2
         lreq = (lreq - 1)/2 + 1
         return
      call jgroup (n, iwk(ipian), iwk(ipjan), maxg, ngp, iwk(ipigp), iwk(ipjgp), iwk(iptt1), iwk(iptt2), ier)
      if ier != 0 :
         ipper = -2
         lreq = (lreq - 1)/2 + 1
         return
      lenigp = ngp + 1
#  compute new ordering of rows/columns of jacobian.
   ipr = ipigp + lenigp
   ipc = ipr
   ipic = ipc + n
   ipisp = ipic + n
   iprsp = (ipisp - 2)/2 + 2
   iesp = lenwk + 1 - iprsp
   if iesp < 0 :
      ipper = -3
      call cntnzu(n, iwk(ipian), iwk(ipjan), nzsut)
      lreq = lenwk - iesp + (3*n + 4*nzsut - 1)/2 + 1
      return
   ibr = ipr - 1
   do i = 1,n
      iwk(ibr+i) = i
   nsp = liwk + 1 - ipisp
   call odrv(n, iwk(ipian), iwk(ipjan), wk, iwk(ipr), iwk(ipic), nsp, iwk(ipisp), 1, iys)
   if iys == 11*n+1 :
      ipper = -4
      return
   if iys != 0 :
      ipper = -3
      call cntnzu(n, iwk(ipian), iwk(ipjan), nzsut)
      lreq = lenwk - iesp + (3*n + 4*nzsut - 1)/2 + 1
      return
 
#  reorder jan and make symbolic lu factorization of matrix.
   ipa = lenwk + 1 - nnz
   nsp = ipa - iprsp
   lreq = max0(12*n/2, 6*n/2+2*n+nnz) + 3
   lreq = lreq + iprsp - 1 + nnz
   if lreq > lenwk :
      ipper = -5
      return
   iba = ipa - 1
   do i = 1,nnz
      wk(iba+i) = 0.0
   ipisp = 2*(iprsp - 1) + 1
   call cdrv(n,iwk(ipr),iwk(ipc),iwk(ipic),iwk(ipian),iwk(ipjan),wk(ipa),wk(ipa),wk(ipa),nsp,iwk(ipisp),wk(iprsp),iesp,5,iys)
   lreq = lenwk - iesp
   if iys == 10*n+1 :
      ipper = -5
      return
   if iys != 0 :
      ipper = -6
      lreq = lenwk
      return
   ipil = ipisp
   ipiu = ipil + 2*n + 1
   nzu = iwk(ipil+n) - iwk(ipil)
   nzl = iwk(ipiu+n) - iwk(ipiu)
   if nnz == n : lreq = lreq + 1
   nsp = nsp + lreq - lenwk
   ipa = lreq + 1 - nnz
   iba = ipa - 1
   ipper = 0
   return
#
# prjs is called to compute and process the matrix p = i - h*el(1)*j , where j is an approximation to the jacobian.
# j is computed by columns, either by the user-supplied routine jac if miter = 1, or by finite differencing if miter = 2.
# if miter = 3, a diagonal approximation to j is used.
# if miter = 1 or 2, and if the existing value of the jacobian (as contained in p) is considered acceptable, then a new value of p is reconstructed 
# from the old value. in any case, when miter is 1 or 2, the p matrix is subjected to lu decomposition in cdrv. 
# p and its lu decomposition are stored (separately) in wk.
#
# in addition to variables described previously, communication with prjs uses the following..
# y     = array containing predicted values on entry.
# ftem  = work array of length n (acor in stode).
# savf  = array containing f evaluated at predicted y.
# wk    = real work space for matrices.  on output it contains the inverse diagonal matrix if miter = 3, and p and its sparse lu decomposition 
#         if miter is 1 or 2.
#         storage of matrix elements starts at wk(3). wk also contains the following matrix-related data..
#         wk(1) = sqrt(uround), used in numerical jacobian increments.
#         wk(2) = h*el0, saved for later use if miter = 3.
# iwk   = integer work space for matrix-related data, assumed to be equivalenced to wk.
#         in addition, wk(iprsp) and iwk(ipisp) are assumed to have identical locations.
# el0   = el(1) (input).
# ierpj = output error flag (in common).
#       = 0 if no error.
#       = 1  if zero pivot found in cdrv.
#       = 2  if a singular matrix arose with miter = 3.
#       = -1 if insufficient storage for cdrv (should not occur here).
#       = -2 if other error found in cdrv (should not occur here).
# jcur  = output flag = 1 to indicate that the jacobian matrix (or approximation) is now current.
#
#--------------------------------------------------------------------------------------------------
def prjs(neq, y, yh, nnyh, ewt, ftem, savf, wk, iwk)

   integer neq, nnyh, iwk
   integer i, imul, j, jj, jok, jmax, jmin, k, kmax, kmin, ng
   double precision y, yh, ewt, ftem, savf, wk
   double precision con, di, fac, hl0, pij, r, r0, rcon, rcont, srur, vnorm
   dimension y(1), yh(nnyh,1), ewt(1), ftem(1), savf(1), wk(1), iwk(1)                                                        
 
   double precision crate, el, elco, hold, rmax, tesco, 
  +   ccmax, el0, h, hmin, hmxi, hu, rc, tn, uround
   integer init, lyh, lewt, lacor, lsavf, lwm, liwm, mxstep, mxhnil, nhnil, ntrep, nslast, nyh,
  +   ialth, ipup, lmax, meo, nqnyh, nslp, 
  +   icf, ierpj, iersl, jcur, jstart, kflag, meth, miter, 
  +   maxord, maxcor, n, nq, nst, nfe, nje, nqu
   common /ls0001/ crate, el(13), elco(13,12), hold, rmax, tesco(3,12),
  +   ccmax, el0, h, hmin, hmxi, hu, rc, tn, uround,
  +   init, lyh, lewt, lacor, lsavf, lwm, liwm, mxstep, mxhnil, nhnil, ntrep, nslast, nyh, 
  +   ialth, ipup, lmax, meo, nqnyh, nslp, 
  +   icf, ierpj, iersl, jcur, jstart, kflag, meth, miter,
  +   maxord, maxcor, n, nq, nst, nfe, nje, nqu
 
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
 
   hl0 = h*el0
   con = -hl0
 
   if miter == 3 :
#     if miter = 3, construct a diagonal approximation to j and p.
      jcur = 1
      nje = nje + 1
      wk(2) = hl0
      ierpj = 0
      r = el0*0.1
      do i = 1,n
         y(i) = y(i) + r*(h*savf(i) - yh(i,2))
      call rhs(neq, tn, y, wk(3))
      nfe = nfe + 1
      do i = 1,n
         r0 = h*savf(i) - yh(i,2)
         di = 0.1*r0 - h*(wk(i+2) - savf(i))
         wk(i+2) = 1.0
         if abs(r0) >= uround/ewt(i) :
            if abs(di) == 0.0 :
               ierpj = 2
               return
            wk(i+2) = 0.1*r0/di
      return
 
#  see whether jacobian should be reevaluated (jok = 0) or not (jok = 1).
   jok = 1
   if nst == 0 or nst >= nslj+50 : jok = 0
   if icf == 1 and abs(rc - 1.0) < ccmxj : jok = 0
   if icf == 2 : jok = 0
   if jok == 1 and abs(con)/conmin > rbig and iplost == 1 : jok = 0
 
   if jok == 0 :
#     miter = 1 or 2, and the jacobian is to be reevaluated.
      jcur = 1
      nje = nje + 1
      nslj = nst
      iplost = 0
      conmin = abs(con)
      if miter == 1 :
#        if miter = 1, call jac, multiply by scalar, and add identity.
         kmin = iwk(ipian)
         do j = 1, n
            kmax = iwk(ipian+j) - 1
            do i = 1,n
               ftem(i) = 0.0
            call fjac(neq, tn, y, j, iwk(ipian), iwk(ipjan), ftem)
            do k = kmin, kmax
               i = iwk(ibjan+k)
               wk(iba+k) = ftem(i)*con
               if i == j : wk(iba+k) = wk(iba+k) + 1.0
            kmin = kmax + 1
      else : # miter == 2
#        if miter = 2, make ngp calls to f to approximate j and p.
         fac = vnorm(n, savf, ewt)
         r0 = 1000.0 * abs(h) * uround * float(n) * fac
         if r0 == 0.0 : r0 = 1.0
         srur = wk(1)
         jmin = iwk(ipigp)
         do ng = 1,ngp
            jmax = iwk(ipigp+ng) - 1
            do j = jmin,jmax
               jj = iwk(ibjgp+j)
               r = max(srur*abs(y(jj)),r0/ewt(jj))
               y(jj) = y(jj) + r
            call rhs(neq, tn, y, ftem)
            do j = jmin,jmax
               jj = iwk(ibjgp+j)
               y(jj) = yh(jj,1)
               r = max(srur*abs(y(jj)),r0/ewt(jj))
               fac = -hl0/r
               kmin =iwk(ibian+jj)
               kmax =iwk(ibian+jj+1) - 1
               do k = kmin,kmax
                 i = iwk(ibjan+k)
                 wk(iba+k) = (ftem(i) - savf(i))*fac
                 if i == jj : wk(iba+k) = wk(iba+k) + 1.0
            jmin = jmax + 1
         nfe = nfe + ngp
 
   else : # jok == 1
 
#     if jok = 1, reconstruct new p from old p.
      jcur = 0
      kmin = iwk(ipian)
      do j = 1,n
        kmax = iwk(ipian+j) - 1
        do k = kmin,kmax
          i = iwk(ibjan+k)
          pij = wk(iba+k)
          if i == j :
             pij = pij - 1.0
             if abs(pij) < psmall :
                iplost = 1
                conmin = min(abs(con0),conmin)
          rcon = con/con0
          pij = pij*rcon
          if i == j : pij = pij + 1.0
          wk(iba+k) = pij
        kmin = kmax + 1
 
#  make numerical factorization of p matrix.
   nlu = nlu + 1
   con0 = con
   ierpj = 0
   do i = 1,n
      ftem(i) = 0.0
   call cdrv(n,iwk(ipr),iwk(ipc),iwk(ipic),iwk(ipian),iwk(ipjan), wk(ipa),ftem,ftem,nsp,iwk(ipisp),wk(iprsp),iesp,2,iys)
   if iys != 0 :
      imul = (iys - 1)/n
      ierpj = -2
      if imul == 8 : ierpj = 1
      if imul == 10 : ierpj = -1
   return

#--------------------------------------------------------------------------------------------------
# this routine manages the solution of the linear system arising from a chord iteration.
# it is called if miter != 0. if miter is 1 or 2, it calls cdrv to accomplish this.
# if miter = 3 it updates the coefficient h*el0 in the diagonal matrix, and then computes the solution.
# communication with slss uses the following variables..
# wk    = real work space containing the inverse diagonal matrix if miter = 3 and the lu decomposition of the matrix otherwise.
#         storage of matrix elements starts at wk(3). wk also contains the following matrix-related data..
#         wk(1) = sqrt(uround) (not used here),
#         wk(2) = hl0, the previous value of h*el0, used if miter = 3.
# iwk   = integer work space for matrix-related data, assumed to be equivalenced to wk.  
#         in addition, wk(iprsp) and iwk(ipisp) are assumed to have identical locations.
# x     = the right-hand side vector on input, and the solution vector on output, of length n.
# tem   = vector of work space of length n, not used in this version.
# iersl = output flag (in common).
#         iersl = 0  if no trouble occurred.
#         iersl = -1 if cdrv returned an error flag (miter = 1 or 2). this should never occur and is considered fatal.
#         iersl = 1  if a singular matrix arose with miter = 3.
# this routine also uses other variables in common.
#--------------------------------------------------------------------------------------------------
def slss(wk, iwk, x, tem)

   integer iwk
   integer i
   double precision wk, x, tem
   double precision di, hl0, phl0, r
   dimension wk(1), iwk(1), x(1), tem(1)
 
   double precision crate, el, elco, hold, rmax, tesco, 
  +   ccmax, el0, h, hmin, hmxi, hu, rc, tn, uround
   integer init, lyh, lewt, lacor, lsavf, lwm, liwm, mxstep, mxhnil, nhnil, ntrep, nslast, nyh,
  +   ialth, ipup, lmax, meo, nqnyh, nslp, 
  +   icf, ierpj, iersl, jcur, jstart, kflag, meth, miter, 
  +   maxord, maxcor, n, nq, nst, nfe, nje, nqu
   common /ls0001/ crate, el(13), elco(13,12), hold, rmax, tesco(3,12),
  +   ccmax, el0, h, hmin, hmxi, hu, rc, tn, uround,
  +   init, lyh, lewt, lacor, lsavf, lwm, liwm, mxstep, mxhnil, nhnil, ntrep, nslast, nyh, 
  +   ialth, ipup, lmax, meo, nqnyh, nslp, 
  +   icf, ierpj, iersl, jcur, jstart, kflag, meth, miter,
  +   maxord, maxcor, n, nq, nst, nfe, nje, nqu
 
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
   if miter == 3 :
      phl0 = wk(2)
      hl0 = h*el0
      wk(2) = hl0
      if hl0 != phl0 :
         r = hl0/phl0
         do i = 1,n
           di = 1.0 - r*(1.0 - 1.0/wk(i+2))
           if abs(di) == 0.0 :
              iersl = 1
              return
           wk(i+2) = 1.0/di
      do i = 1,n
         x(i) = wk(i+2)*x(i)
   else :
      call cdrv(n,iwk(ipr),iwk(ipc),iwk(ipic),iwk(ipian),iwk(ipjan),wk(ipa),x,x,nsp,iwk(ipisp),wk(iprsp),iesp,4,iersl)
      if iersl != 0 : iersl = -1
   return

#--------------------------------------------------------------------------------------------------
#  symmetric reordering of sparse symmetric matrix.
#
#  the nonzero entries of the matrix m are assumed to be stored symmetrically in (ia,ja,a) format
#  (i.e., not both m(i,j) and m(j,i) are stored if i ne j).
#
#  sro does not rearrange the order of the rows, but does move nonzeroes from one row to another to ensure that if m(i,j) will be
#  in the upper triangle of m with respect to the new ordering, then m(i,j) is stored in row i (and thus m(j,i) is not stored), whereas
#  if m(i,j) will be in the strict lower triangle of m, then m(j,i) is stored in row j (and thus m(i,j) is not stored).
#
#  additional parameters
#  q - integer one-dimensional work array.  dimension = n
#  r - integer one-dimensional work array.  dimension = number of nonzero entries in the upper triangle of m
#  dflag - logical variable.  if dflag = True, then store nonzero diagonal elements at the beginning of the row
#--------------------------------------------------------------------------------------------------
def sro(n, ip, ia, ja, a, q, r, dflag)

   integer ip(1), ia(n), ja(1), q(1), r(1)
   double precision a(1), ak
   logical dflag
 
#  phase 1 -- find row in which to store each nonzero initialize count of nonzeroes to be stored in each row
   do i = 1,n
      q(i) = 0
 
#  for each nonzero element a(j)
   do i = 1,n
     jmin = ia(i)
     jmax = ia(i+1) - 1
     if jmin <= jmax :
        do j = jmin,jmax
#         find row (=r(j)) and column (=ja(j)) in which to store a(j) ...
          k = ja(j)
          if ip(k)<ip(i) : ja(j) = i
          if ip(k)>=ip(i) : k = i
          r(j) = k
#         ... and increment count of nonzeroes (=q(r(j)) in that row
          q(k) = q(k) + 1
 
#  phase 2 -- find new ia and permutation to apply to (ja,a) determine pointers to delimit rows in permuted (ja,a)
   do i = 1,n
      ia(i+1) = ia(i) + q(i)
      q(i) = ia(i+1)
 
#  determine where each (ja(j),a(j)) is stored in permuted (ja,a) for each nonzero element (in reverse order)
   ilast = 0
   jmin = ia(1)
   jmax = ia(n+1) - 1
   j = jmax
   do jdummy = jmin,jmax
     i = r(j)
     if not dflag or ja(j) != i or i == ilast :
#       put (off-diagonal) nonzero in last unused location in row
        q(i) = q(i) - 1
        r(j) = q(i)
     else :
#       if dflag, then put diagonal nonzero at beginning of row
        r(j) = ia(i)
        ilast = i
     j = j - 1
 
#  phase 3: permute (ja,a) to upper triangular form (wrt new ordering)
   do j = jmin, jmax
      while r(j) != j :
         k = r(j)
         r(j) = r(k)
         r(k) = k
         jak = ja(k)
         ja(k) = ja(j)
         ja(j) = jak
         ak = a(k)
         a(k) = a(j)
         a(j) = ak
   return

#--------------------------------------------------------------------------------------------------
# stode performs one step of the integration of an initial value problem for a system of ordinary differential equations.
# note.. stode is independent of the value of the iteration method indicator miter, when this is != 0, and hence is independent
# of the type of chord method used, or the jacobian structure.
#
# communication with stode is done with the following variables..
#
# n      = the number of first-order differential equations.
# neq    = integer array containing problem size in neq, and passed as the neq argument in all calls to f and jac.
# y      = an array of length >= n used as the y argument in all calls to f and jac.
# yh     = an nnyh by maxord+1 array containing the dependent variables and their approximate scaled derivatives, where.  
#          yh(i,j+1) contains the approximate j-th derivative of y(i), scaled by h**j/factorial(j) (j = 0,1,...,nq).  
#          on entry for the first step, the first two columns of yh must be set from the initial values.
# nnyh   = a constant integer >= n, the first dimension of yh.
# yh1    = a one-dimensional array occupying the same space as yh.
# ewt    = an array of length n containing multiplicative weights for local error measurements.  
#          local errors in y(i) are compared to 1.0/ewt(i) in various error tests.
# savf   = an array of working storage, of length n. also used for input of yh(*,maxord+2) when jstart = -1
#          and maxord < the current order nq.
# acor   = a work array of length n, used for the accumulated corrections.  
#          on a successful return, acor(i) contains the estimated one-step local error in y(i).
# wm,iwm = real and integer work arrays associated with matrix operations in chord iteration (miter != 0).
# fjac   = name of routine to evaluate and preprocess jacobian matrix and p = i - h*el0*jac, if a chord method is being used.
# slss   = name of routine to solve linear system in chord iteration.
# h      = the step size to be attempted on the next step.
#          h is altered by the error control algorithm during the problem.  
#          h can be either positive or negative, but its sign must remain constant throughout the problem.
# hmin   = the minimum absolute value of the step size h to be used.
# hmxi   = inverse of the maximum absolute value of h to be used. hmxi = 0.0 is allowed and corresponds to an infinite hmax.
#          hmin and hmxi may be changed at any time, but will not take effect until the next change of h is considered.
# tn     = the independent variable. tn is updated on each step taken.
# jstart = an integer used for input only, with the following values and meanings..
#               0  perform the first step.
#           >0  take a new step continuing from the last.
#              -1  take the next step with a new value of h, maxord, n, meth, miter, and/or matrix parameters.
#              -2  take the next step with a new value of h, but with other inputs unchanged.
#          on return, jstart is set to 1 to facilitate continuation.
# kflag  = a completion code with the following meanings..
#               0  the step was succesful.
#              -1  the requested error could not be achieved.
#              -2  corrector convergence could not be achieved.
#              -3  fatal error in fjac or slss.
#          a return with kflag = -1 or -2 means either abs(h) = hmin or 10 consecutive failures occurred.
#          on a return with kflag negative, the values of tn and the yh array are as of the beginning of the last step, 
#          and h is the last step size attempted.
# maxord = the maximum order of integration method to be allowed.
# maxcor = the maximum number of corrector iterations allowed.
# meth   = the basic linear multistep method..
#               1  means the implicit adams method.
#               2  means the method based on backward differentiation formulas (bdf-s).
# miter  = the corrector iteration method..
#               0  means functional iteration (no jacobian matrix is involved).
#               1  means chord iteration with a user-supplied sparse jacobian, given by subroutine jac.
#               2  means chord iteration with an internally generated (difference quotient) sparse jacobian
#                  (using ngp extra calls to f per df/dy value, where ngp is an optional output described below.)
#               3  means chord iteration with an internally generated diagonal jacobian approximation.
#                  (using 1 extra call to f per df/dy evaluation).
#          if miter = 1 the user must supply a subroutine jac (the name is arbitrary) as described in lsodes.
#          otherwise, a dummy argument can be used.
#--------------------------------------------------------------------------------------------------
def stode(neq, y, yh, nnyh, yh1, ewt, savf, acor, wm, iwm)

   integer neq, nnyh, iwm
   integer i, i1, j, jb, m, ncf, newq
   double precision y, yh, yh1, ewt, savf, acor, wm, rmx
   double precision dcon, ddn, del, delp, dsm, dup, exdn, exsm, exup, r, rh, rhdn, rhsm, rhup, told, vnorm
   dimension y(1), yh(nnyh,1), yh1(1), ewt(1), savf(1), acor(1), wm(1), iwm(1)
 
   double precision crate, el, elco, hold, rmax, tesco, 
  +   ccmax, el0, h, hmin, hmxi, hu, rc, tn, uround
   integer init, lyh, lewt, lacor, lsavf, lwm, liwm, mxstep, mxhnil, nhnil, ntrep, nslast, nyh,
  +   ialth, ipup, lmax, meo, nqnyh, nslp, 
  +   icf, ierpj, iersl, jcur, jstart, kflag, meth, miter, 
  +   maxord, maxcor, n, nq, nst, nfe, nje, nqu
   common /ls0001/ crate, el(13), elco(13,12), hold, rmax, tesco(3,12),
  +   ccmax, el0, h, hmin, hmxi, hu, rc, tn, uround, 
  +   init, lyh, lewt, lacor, lsavf, lwm, liwm, mxstep, mxhnil, nhnil, ntrep, nslast, nyh, 
  +   ialth, ipup, lmax, meo, nqnyh, nslp, 
  +   icf, ierpj, iersl, jcur, jstart, kflag, meth, miter, 
  +   maxord, maxcor, n, nq, nst, nfe, nje, nqu
 
   logical exit0, exit1, exit2, exit3
 
   kflag = 0
   told = tn
   ncf = 0
   ierpj = 0
   iersl = 0
   jcur = 0
   icf = 0
   delp = 0.0
 
   if jstart == 0 :
#     on the first call, the order is set to 1, and other variables are initialized.  
      nq = 1
      ialth = 2
#     rmax is the maximum ratio by which h can be increased in a single step.  
#     it is initially 1.e4 to compensate for the small initial h, but then is normally equal to 10.
#     if a failure occurs (in corrector convergence or error test), rmax is set at 2 for the next increase.
      rmax = 10000.0
      rc = 0.0
      el0 = 1.0
      crate = 0.7
      hold = h
      meo = meth
      nslp = 0
      ipup = miter
#     cfode is called to reset the coefficients of the method.
      call cfode(meth, elco, tesco)
#     el vector and related constants are reset at the start of the problem.
      call order_changed()
 
   elif jstart == -1 :
#     preliminaries needed when jstart = -1. 
#     if h is to be changed, yh must be rescaled.
#     if h or meth is being changed, ialth is reset to nq + 1 to prevent further changes in h for that many steps.
#     ipup is set to miter to force a matrix update.
      ipup = miter
#     if an order increase is about to be considered (ialth = 1), ialth is reset to 2 to postpone consideration one more step.
      if ialth == 1 : ialth = 2
#     if the caller has changed meth, cfode is called to reset the coefficients of the method.
      if meth != meo :
         call cfode(meth, elco, tesco)
         meo = meth
#        el vector and related constants are reset when the order nq is changed.
         if nq <= maxord :
            ialth = l
            call order_changed()
      if nq <= maxord :
#        if h is being changed, the h ratio rh is reset.  
         if h != hold :
            rh = h/hold
            h = hold
      else :
#        if the caller has changed maxord to a value less than the current order nq, nq is reduced to maxord, and a new h chosen accordingly.
         nq = maxord
         call order_changed()
         ddn = vnorm(n, savf, ewt)/tesco(1,nq+1)
         exdn = 1.0/float(nq + 1)
         rhdn = 1.0/(1.3*ddn**exdn + 0.0000013)
         rh = min(rhdn,1.0)
         if h == hold :
            rh = max(rh,hmin/abs(h))
         else :
            rh = min(rh,abs(h/hold))
            h = hold
#     if h is being changed, the h ratio rh is checked against rmax, hmin, and hmxi, and the yh array rescaled.  
      if h != hold :
         call step_changed(rh, nyh, yh)
 
   elif jstart == -2 and h != hold :
#     if only h is changed, the h ratio rh is reset, checked against rmax, hmin, and hmxi, and the yh array rescaled.
      rh = h/hold
      h = hold
      call step_changed(rh, nyh, yh)
 
   exit0 = False
   while True :
#     rc is the ratio of new to old values of the coefficient h*el(1).
#     when rc differs from 1 by more than 0.3, ipup is set to miter to force fjac to be called, if a jacobian is involved.
#     in any case, fjac is called at least every 20 steps. 0.3 is maximum relative change in h*el0 before fjac is called.
      if abs(rc - 1.0) > 0.3 or nst >= nslp + 20 : ipup = miter
 
      tn = tn + h
 
#     the predicted values is computed by effectively multiplying the yh array by the pascal triangle matrix.
      i1 = nqnyh + 1
      do jb = 1,nq
         i1 = i1 - nnyh
         do i = i1,nqnyh
            yh1(i) = yh1(i) + yh1(i+nnyh)
 
      exit1 = False
      while not exit1 :
#        up to maxcor corrector iterations are taken.  
#        a convergence test is made on the r.m.s. norm of each correction, weighted by the error weight vector ewt.
#        the sum of the corrections is accumulated in the vector acor(i). the yh array is not altered in the corrector loop.
         m = 0
         do i = 1,n
            y(i) = yh(i,1)
         call rhs(neq, tn, y, savf)
         nfe = nfe + 1
         if ipup > 0 :
#           if indicated, the matrix p = i - h*el(1)*j is reevaluated and preprocessed before starting the corrector iteration.  
            call prjs(neq, y, yh, nnyh, ewt, acor, savf, wm, iwm)
#           ipup is set to 0 as an indicator that the matrix was evaluated.
            ipup = 0
            rc = 1.0
            nslp = nst
            crate = 0.7
            if ierpj != 0 :
               icf = 2
               ncf = ncf + 1
#              the values of tn and the yh array are as of the beginning of the last step.
               call set_back(told, yh1)
#              kflag  = a completion code -3: fatal error in prjs.
               kflag = -3
               hold = h
               jstart = 1
               return
#              ------
         
         do i = 1,n
            acor(i) = 0.0
         
         exit2 = False
         exit3 = False
         while not exit2 :
            if miter == 0 :
#              in case of functional iteration, update y directly from the result of the last function evaluation.
               do i = 1,n
                  savf(i) = h*savf(i) - yh(i,2)
                  y(i) = savf(i) - acor(i)
               del = vnorm (n, y, ewt)
               do i = 1,n
                  y(i) = yh(i,1) + el(1)*savf(i)
                  acor(i) = savf(i)
            else :
#              in case of chord method, compute corrector error, and solve linear system with that as right-hand side and p as coefficient matrix.
               do i = 1,n
                  y(i) = h*savf(i) - (yh(i,2) + acor(i))
               call slss(wm, iwm, y, savf)
#              slss failed
               if iersl < 0 :
                  icf = 2
                  ncf = ncf + 1
#                 the values of tn and the yh array are set as of the beginning of the last step.
                  call set_back(told, yh1)
#                 kflag  = a completion code -3: fatal error in slss.
                  kflag = -3
                  hold = h
                  jstart = 1
                  return
#                 ------
 
#              slss was successful
               if iersl == 0 :
                  del = vnorm(n, y, ewt)
                  do i = 1,n
                     acor(i) = acor(i) + y(i)
                     y(i) = yh(i,1) + el(1)*acor(i)
 
            if iersl == 0 :
#              slss succesful or miter = 0. test for convergence.
#              if m>0, an estimate of the convergence rate constant is stored in crate, and this is used in the test.
               if m > 0 : crate = max(0.2*crate,del/delp)
               dcon = del*min(1.0,1.5*crate)/(tesco(2,nq)*0.5/float(nq+2))
 
               if dcon <= 1.0 :
#                 the corrector has converged.  jcur is set to 0 to signal that the jacobian involved may need updating later.
                  jcur = 0
#                 the local error test is made
                  if m == 0 :
                     dsm = del/tesco(2,nq)
                  else :
                     dsm = vnorm(n, acor, ewt)/tesco(2,nq)
 
                  if dsm <= 1.0 :
#                    a successful step, update the yh array. consider changing h if ialth = 1.  otherwise decrease ialth by 1.
#                    if ialth is then 1 and nq < maxord, then acor is saved for use in a possible order increase on the next step.
#                    if a change in h is considered, an increase or decrease in order by one is considered also.
#                    a change in h is made only if it is by a factor of at least 1.1. 
#                    if not, ialth is set to 3 to prevent testing for that many steps.
                     kflag = 0
                     nst = nst + 1
                     hu = h
                     nqu = nq
                     do j = 1,nq+1
                        do i = 1,n
                           yh(i,j) = yh(i,j) + el(j)*acor(i)
                     ialth = ialth - 1
                     if ialth == 0 :
#                       rhdn, rhsm, and rhup are computed, by which h could be multiplied at orders nq - 1, nq, or nq + 1, respectively.
                        rhup = 0.0
                        if nq < maxord :
                           do i = 1,n
                              savf(i) = acor(i) - yh(i,maxord+1)
                           dup = vnorm(n, savf, ewt)/tesco(3,nq)
                           exup = 1.0/float(nq+2)
                           rhup = 1.0/(1.4*dup**exup + 0.0000014)
                        
                        exsm = 1.0/float(nq + 1)
                        rhsm = 1.0/(1.2*dsm**exsm + 0.0000012)
                        rhdn = 0.0
                        if nq != 1 :
                           ddn = vnorm (n, yh(1,nq+1), ewt)/tesco(1,nq)
                           exdn = 1.0/float(nq)
                           rhdn = 1.0/(1.3*ddn**exdn + 0.0000013)
 
#                       the largest of rhdn, rhsm, and rhup is determined and the new order chosen accordingly.
                        rmx = max(rhdn, rhsm, rhup)
                        if rmx == rhdn :
                           newq = nq - 1
                           rh = rhdn
                        elif rmx == rhup :
                           newq = nq + 1
                           rh = rhup
                           if rh >= 1.1 :
                              r = el(newq)/float(newq)
                              do i = 1,n
                                 yh(i,newq+1) = acor(i)*r
                        else : 
                           newq = nq
                           rh = rhsm
 
                        if rh < 1.1 :
                           ialth = 3
                        else : 
                           if newq != nq :
#                             if there is a change of order, reset nq and the coefficients.
                              nq = newq
#                             el vector and related constants are reset when the order nq is changed.
                              call order_changed()
#                          if h is being changed, the h ratio rh is checked against rmax, hmin, and hmxi, and the yh array rescaled.  
                           if rh != 1.0 :
                              rh = max(rh,hmin/abs(h))
#                             the h ratio rh is checked against rmax, hmin, and hmxi, and the yh array rescaled.  
                              call step_changed(rh, nyh, yh)
                              rmax = 10.0
 
                     elif ialth == 1 and nq < maxord :
#                       if the order is to be increased, we compute one additional scaled derivative.
                        do i = 1,n
                           yh(i,maxord+1) = acor(i)
                     
                     #successful step
                     r = 1.0/tesco(2,nqu)
                     do i = 1,n
                        acor(i) = acor(i)*r
                     hold = h
                     jstart = 1
                     return
#                    ------ 
 
                  else : #dsm > 1
 
#                    the local error test failed. kflag is a count of failures.
                     kflag = kflag - 1
#                    restore tn and the yh array to their previous values, and prepare to try the step again.
                     call set_back(told, yh1)
                     if abs(h) <= hmin*1.00001 :
#                       kflag  = a completion code -1: the requested error could not be achieved. abs(h) = hmin.
                        kflag = -1
#                       h is saved in hold to allow the caller to change h on the next step.
                        hold = h
                        jstart = 1
                        return
#                       ------
#                    if 3 or more failures have occured...
                     if kflag <= -3 :
#                       if 10 failures have occurred, exit with kflag = -1.
                        if kflag == -10 :
#                          kflag = a completion code -1: the requested error could not be achieved. 10 consecutive failures occurred.
#                          the values of tn and the yh array are as of the beginning of the last step, and h is the last step size attempted.
#                          h is saved in hold to allow the caller to change h on the next step.
                           kflag = -1
                           hold = h
                           jstart = 1
                           return
#                          ------
#                       it is assumed that the derivatives that have accumulated in the yh array have errors of the wrong order.  
#                       hence the first derivative is recomputed, and the order is set to 1. 
                        do i = 1,n
                           y(i) = yh(i,1)
                        call rhs(neq, tn, y, savf)
                        nfe = nfe + 1
                        do i = 1,n
                           yh(i,2) = h*savf(i)
#                       then h is reduced by a factor of 10, and the step is retried, until it succeeds or h reaches hmin.
                        rh = 0.1
                        rh = max(hmin/abs(h),rh)
                        h = h*rh
                        ipup = miter
                        if nq != 1 :
                           nq = 1
#                          el vector and related constants are reset because the order nq is changed.
                           call order_changed()
                        ialth = 5
                        
                     else : # one or two failures...
 
#                       1 or 2 failures of the step
#                       factor rhsm by which h could be multiplied at order nq
                        exsm = 1.0/float(nq + 1)
                        rhsm = 1.0/(1.2*dsm**exsm + 0.0000012)
 
#                       factor rhdn by which h could be multiplied at order nq-1
                        rhdn = 0.0 
                        if nq != 1 :
                           ddn = vnorm(n, yh(1,nq+1), ewt)/tesco(1,nq)
                           exdn = 1.0/float(nq)
                           rhdn = 1.0/(1.3*ddn**exdn + 0.0000013)
 
#                       the largest of these is determined and the new order chosen accordingly.
                        if rhdn > rhsm :
                           newq = nq - 1
                           rh = rhdn
                           if rh > 1.0 : rh = 1.0
                        else :
                           newq = nq
                           rh = rhsm
                        if kflag == -2 : rh = min(rh,0.2)
 
#                       if there is a change of order, reset nq and the coefficients.
#                       h is reset according to rh and the yh array is rescaled. then remake the step otherwise.
                        if newq != nq :
                           nq = newq
#                          el vector and related constants are reset when the order nq is changed.
                           call order_changed()
 
#                       if h is being changed, the h ratio rh is checked against rmax, hmin, and hmxi, and the yh array rescaled.
                        if rh < 1.0 :
                           rh = max(rh,hmin/abs(h))
                           call step_changed(rh, nyh, yh)
                     exit0 = True
 
               else : # dcon > 1
 
#                 the corrector iteration failed to converge.
                  m = m + 1
                  if m < maxcor and (m == 1 or del <= 2.0*delp) :
                     delp = del
                     call rhs(neq, tn, y, savf)
                     nfe = nfe + 1
                     exit3 = True
 
            if exit0 :
               exit2 = True
               exit1 = True
            elif exit3 :
               exit3 = False
            else :
               if miter != 0 and jcur == 0 :
#                 if miter != 0 and the jacobian is out of date, pjac is called for the next try.
                  icf = 1
                  ipup = miter
                  exit2 = True
               else :
                  exit2 = True
                  exit1 = True
 
      if not exit0 :
#        the yh array is retracted to its values before prediction, and h is reduced, if possible.
         icf = 2
         ncf = ncf + 1
         
         call set_back(told, yh1)
         
         if abs(h) <= hmin*1.00001 or ncf == 10 :
#           kflag  = a completion code -2: corrector convergence could not be achieved. either abs(h) = hmin or 10 consecutive failures occurred.
            kflag = -2
#           h is saved in hold to allow the caller to change h on the next step.
            hold = h
            jstart = 1
            return
#           ------
         
         rh = 0.25
         ipup = miter
         rh = max(rh,hmin/abs(h))
#        if h is being changed, the h ratio rh is checked against rmax, hmin, and hmxi, and the yh array rescaled.
         call step_changed(rh, nyh, yh)
      else :
         exit0 = False


#--------------------------------------------------------------------------------------------------
def order_changed()

   integer i
 
   double precision crate, el, elco, hold, rmax, tesco, 
  +   ccmax, el0, h, hmin, hmxi, hu, rc, tn, uround
   integer init, lyh, lewt, lacor, lsavf, lwm, liwm, mxstep, mxhnil, nhnil, ntrep, nslast, nyh,
  +   ialth, ipup, lmax, meo, nqnyh, nslp, 
  +   icf, ierpj, iersl, jcur, jstart, kflag, meth, miter, 
  +   maxord, maxcor, n, nq, nst, nfe, nje, nqu
   common /ls0001/ crate, el(13), elco(13,12), hold, rmax, tesco(3,12),
  +   ccmax, el0, h, hmin, hmxi, hu, rc, tn, uround, 
  +   init, lyh, lewt, lacor, lsavf, lwm, liwm, mxstep, mxhnil, nhnil, ntrep, nslast, nyh, 
  +   ialth, ipup, lmax, meo, nqnyh, nslp, 
  +   icf, ierpj, iersl, jcur, jstart, kflag, meth, miter, 
  +   maxord, maxcor, n, nq, nst, nfe, nje, nqu
 
#  el vector and related constants are reset because the order nq is changed.
   do i = 1,nq+1
      el(i) = elco(i,nq)
   nqnyh = nq*nyh
   rc = rc*el(1)/el0
   el0 = el(1)
#  ialth is set to nq + 1 to prevent a change of h for that many steps.
   ialth = nq + 1
   return


#--------------------------------------------------------------------------------------------------
def step_changed(rh, nnyh, yh)

   integer i, j, nnyh
   double precision rh, r, yh(nnyh,1)
 
   double precision crate, el, elco, hold, rmax, tesco, 
  +   ccmax, el0, h, hmin, hmxi, hu, rc, tn, uround
   integer init, lyh, lewt, lacor, lsavf, lwm, liwm, mxstep, mxhnil, nhnil, ntrep, nslast, nyh,
  +   ialth, ipup, lmax, meo, nqnyh, nslp, 
  +   icf, ierpj, iersl, jcur, jstart, kflag, meth, miter, 
  +   maxord, maxcor, n, nq, nst, nfe, nje, nqu
   common /ls0001/ crate, el(13), elco(13,12), hold, rmax, tesco(3,12),
  +   ccmax, el0, h, hmin, hmxi, hu, rc, tn, uround, 
  +   init, lyh, lewt, lacor, lsavf, lwm, liwm, mxstep, mxhnil, nhnil, ntrep, nslast, nyh, 
  +   ialth, ipup, lmax, meo, nqnyh, nslp, 
  +   icf, ierpj, iersl, jcur, jstart, kflag, meth, miter, 
  +   maxord, maxcor, n, nq, nst, nfe, nje, nqu
 
#  if h is being changed, the h ratio rh is checked against rmax, hmin, and hmxi, and the yh array rescaled.
   rh = min(rh,rmax)
   rh = rh/max(1.0,abs(h)*hmxi*rh)
   r = 1.0
   do j = 2,nq+1
      r = r*rh
      do i = 1,n
         yh(i,j) = yh(i,j)*r
   h = h*rh
   rc = rc*rh
#  ialth is set to nq + 1 to prevent a change of h for that many steps.
   ialth = nq + 1
   return

#--------------------------------------------------------------------------------------------------
def set_back(told, yh1)

   integer i1, jb, i
   double precision told, yh1(1)

   double precision crate, el, elco, hold, rmax, tesco, 
  +   ccmax, el0, h, hmin, hmxi, hu, rc, tn, uround
   integer init, lyh, lewt, lacor, lsavf, lwm, liwm, mxstep, mxhnil, nhnil, ntrep, nslast, nyh,
  +   ialth, ipup, lmax, meo, nqnyh, nslp, 
  +   icf, ierpj, iersl, jcur, jstart, kflag, meth, miter, 
  +   maxord, maxcor, n, nq, nst, nfe, nje, nqu
   common /ls0001/ crate, el(13), elco(13,12), hold, rmax, tesco(3,12),
  +   ccmax, el0, h, hmin, hmxi, hu, rc, tn, uround, 
  +   init, lyh, lewt, lacor, lsavf, lwm, liwm, mxstep, mxhnil, nhnil, ntrep, nslast, nyh, 
  +   ialth, ipup, lmax, meo, nqnyh, nslp, 
  +   icf, ierpj, iersl, jcur, jstart, kflag, meth, miter, 
  +   maxord, maxcor, n, nq, nst, nfe, nje, nqu

#  the maximum ratio by which h can be increased in a single step.  
   rmax = 2.0
   
#  the values of tn and the yh array are set back as of the beginning of the last step.
   tn = told
   
   i1 = nqnyh + 1
   do jb = 1,nq
      i1 = i1 - nyh
      do i = i1,nqnyh
         yh1(i) = yh1(i) - yh1(i+nyh)
   return

#--------------------------------------------------------------------------------------------------
# this function routine computes the weighted root-mean-square norm of the vector of length n contained in the array v, 
# with weights contained in the array w of length n..
# vnorm = sqrt( (1/n) * sum( v(i)*w(i) )**2 )
#--------------------------------------------------------------------------------------------------
def vnorm(n, v, w)

   sum = 0.0
   for i in range(0,n)
      sum = sum + (v[i]*w[i])**2

   return math.sqrt(sum/float(n))
