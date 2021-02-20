import math
import sys

class Lsodes:

   MAXORD = 12

   def solve(self):
      solve_(self)

   def __init__(self, input):
      initialize(self, input)

#--------------------------------------------------------------------------------------------------
def initialize(lsodes, input):
#  initial values of unknowns
   lsodes.y = input['y']
   if len(lsodes.y) == 0 : sys.exit('***ERROR: lsodes : length of y is zero')

#  number of equations
   lsodes.neq = len(lsodes.y)

#  time
   lsodes.t = input['t']
   if lsodes.t < 0 : sys.exit('***ERROR: lsodes : initial time t is negative')

#  output time
   lsodes.tout = input['tout']
   if lsodes.tout < 0 : sys.exit('***ERROR: lsodes : output time tout is negative')
   if lsodes.tout < lsodes.t : sys.exit('***ERROR: lsodes : output time tout is less than time t')

#  relative tolerance (scalar or vector)
   if isinstance(input['rtol'], list): #vector
      if len(input['rtol']) != lsodes.neq : 
         sys.exit('***ERROR: lsodes : length of rtol not equal to number of equations neq')
      else:
         lsodes.rtol = input['rtol']
   else: #scalar
      lsodes.rtol = [input['rtol']] * lsodes.neq

#  absolute tolerance (scalar or vector)
   if isinstance(input['atol'], list): #vector
      if len(input['atol']) != lsodes.neq : 
         sys.exit('***ERROR: lsodes : length of atol not equal to number of equations neq')
      else:
         lsodes.rtol = input['atol']
   else: #scalar
      lsodes.atol = [input['atol']] * lsodes.neq

#  output time flag
   lsodes.itask = input['itask']
   if lsodes.itask < 1 or lsodes.itask > 5:
      sys.exit('***ERROR: lsodes : itask (1 to 5) is wrong: ' + str(lsodes.itask))

#  istate flag (1: first call, 2: not first call, 3: not first call with changed input parameters)
   lsodes.istate = input['istate']
   if lsodes.istate < 1 or lsodes.istate > 3:
      sys.exit('***ERROR: lsodes : istate (1 to 3) is wrong: ' + str(lsodes.istate))

#  the element threshhold for sparsity determination.
#  if the absolute value of an estimated jacobian element is <= seth, 
#  it will be assumed to be absent in the structure. the default value of seth is 0.
   if 'seth' in input:
      lsodes.seth = input['seth']
      if lsodes.seth < 0:
         sys.exit('***ERROR: lsodes : seth is negative: ' + str(lsodes.seth))
   else:
      lsodes.seth = 0

#     current time
      lsodes.tn = lsodes.t

#     time step
      lsodes.h = 1.0
      
#     counter of steps
      lsodes.nst = 0

#     counter of non-zeros in jacobian
      lsodes.nnz = 0

#     number of rhs evaluations
      lsodes.nfe = 0

#     number of extra rhs evaluations needed for each jacobian evaluation
      lsodes.ngp = 0

#     number of nonzero elements in the strict lower triangle of the lu factorization used in the chord iteration
      lsodes.nzl = 0

#     number of nonzero elements in the strict lower triangle of the lu factorization used in the chord iteration
#     the total number of nonzeros in the factorization is therefore nzl + nzu + neq.
      lsodes.nzu = 0

#     the nordsieck history array, of size nyh by (nqcur+1), where nyh is the initial value of neq.  
#     for j = 0,1,...,nqcur, column j+1 of yh contains hcur**j/factorial(j) times the j-th derivative of 
#     the interpolating polynomial currently representing the solution, evaluated at t = tcur.
      lsodes.yh = [list(lsodes.y)]

#     initial call to rhs
      lsodes.f0 = rhs(lsodes.t, lsodes.y)

#     error weight vector
      lsodes.ewt = []
      for i in range(0,lsodes.neq):
         e = lsodes.rtol[i]*abs(lsodes.y[i]) + lsodes.atol[i]
         if e <= 0.0:
            sys.exit('***ERROR: lsodes : ewt is non positive: ' + str(e))
         else:
            lsodes.ewt.append(1.0/e)
            fac = 1.0 + 1.0/(float(i+1)+1.0)
#           perturb y for structure determination.
            lsodes.y[i] += fac*math.copysign(e,lsodes.y[i])

#     the block computes jacobian structure from results of n + 1 calls to rhs
      lsodes.savf = rhs(lsodes.t, lsodes.y)
#     calculate ia and ja
      lsodes.ia = [0]
      lsodes.ja = []
      k = 0
#     j is a column index
      for j in range(0,lsodes.neq):
#        ja is a vector of length nnz: ja(k) is an index of row (0 ... neq-1) the non-zero a[k] belongs to
         lsodes.ja.append(j) # diagonal element
         k += 1
         e = 1.0/lsodes.ewt[j]
         dyj = math.copysign(e,lsodes.y[j])
         yj = lsodes.y[j]
         lsodes.y[j] += dyj
         lsodes.ftem = rhs(lsodes.t, lsodes.y)
         lsodes.y[j] = yj
#        i is a row index
         for i in range(0,lsodes.neq):
            dq = (lsodes.ftem[i] - lsodes.savf[i])/dyj
            if abs(dq) > lsodes.seth and i != j:
               lsodes.ja.append(i) # non-diagonal element
               k += 1
#        ia is a vector of length neq+1: ia(j) is index k (0 ... nnz-1) of first non-zero a[k] of column j. ia(neq+1) = nnz.
         lsodes.ia.append(k)
#     restore y from yh
      lsodes.y = lsodes.yh[0]

#     the block constructs groupings of the column indices of the jacobian matrix
      exit1 = False 
      ncol = 0
      igp = []
      jgp = []
      ng = 0
#     working array
      jdone = [0]*lsodes.neq
      while True :
         igp.append(ncol)
#        working array
         incl = [0]*lsodes.neq
#        j is column
         for j in range(0,lsodes.neq):
#           reject column j if it is already in a group
            if jdone[j] != 1 :
               kmin = lsodes.ia[j]
               kmax = lsodes.ia[j+1]
               for k in range(kmin,kmax):
#                 reject column j if it overlaps any column already in this group
                  i = lsodes.ja[k]
                  if incl[i] == 1 :
                     exit1 = True                  
                     break
               if exit1:
                  exit1 = False
               else:
#                 accept column j into group ng
                  jgp.append(j)
                  ncol += 1
                  jdone[j] = 1
                  for k in range(kmin,kmax):
                     i = lsodes.ja[k]
                     incl[i] = 1
#        stop if this group is empty (grouping is complete).
         if ncol == igp[ng] : break
         ng += 1
      ngrp = ng
      print('igp: ', igp)
      print('jgp: ', jgp)

#--------------------------------------------------------------------------------------------------
#   ---00--- ---01--- ---02--- ---03--- ---04--- ---05--- ---06--- ---07--- ---08--- ---09--- ---10--- ---11--- 
#00:   X        .        .        .        .        .        .        .        .        .        .        .
#01:   X        X        X        X        X        .        .        .        .        .        .        X
#02:   .        X        X        X        .        X        .        .        .        X        .        .
#03:   .        X        X        X        .        .        .        .        .        .        .        .
#04:   .        X        .        .        X        .        .        .        .        .        .        X
#05:   .        .        X        .        .        X        .        .        .        X        .        .
#06:   .        .        .        .        .        .        X        .        .        X        .        X
#07:   .        .        .        .        .        .        .        X        .        X        .        .
#08:   .        .        .        X        X        X        X        .        .        .        .        .
#09:   .        .        X        .        .        X        X        X        .        X        .        X
#10:   .        .        .        .        .        .        .        X        .        .        .        .
#11:   .        X        .        .        X        .        X        .        .        X        .        X
#--------------------------------------------------------------------------------------------------
def rhs(t, y):
   rk = [0.1, 10.0, 50.0, 2.5, 0.1, 10.0, 50.0, 2.5, 50.0, 5.0, 50.0, 50.0, 50.0, 30.0, 100.0, 2.5, 100.0, 2.5, 50.0, 50.0]
   ydot = [0.0]*len(y)
   ydot[0] = -rk[0]*y[0]
   ydot[1] = rk[0]*y[0] + rk[10]*rk[13]*y[3] + rk[18]*rk[13]*y[4] - rk[2]*y[1]*y[2] - rk[14]*y[1]*y[11] - rk[1]*y[1]
   ydot[2] = rk[1]*y[1] - rk[4]*y[2] - rk[2]*y[1]*y[2] - rk[6]*y[9]*y[2] + rk[10]*rk[13]*y[3] + rk[11]*rk[13]*y[5]
   ydot[3] = rk[2]*y[1]*y[2] - rk[10]*rk[13]*y[3] - rk[3]*y[3]
   ydot[4] = rk[14]*y[1]*y[11] - rk[18]*rk[13]*y[4] - rk[15]*y[4]
   ydot[5] = rk[6]*y[9]*y[2] - rk[11]*rk[13]*y[5] - rk[7]*y[5]
   ydot[6] = rk[16]*y[9]*y[11] - rk[19]*rk[13]*y[6] - rk[17]*y[6]
   ydot[7] = rk[8]*y[9] - rk[12]*rk[13]*y[7] - rk[9]*y[7]
   ydot[8] = rk[3]*y[3] + rk[15]*y[4] + rk[7]*y[5] + rk[17]*y[6]
   ydot[9] = rk[4]*y[2] + rk[11]*rk[13]*y[5] + rk[19]*rk[13]*y[6] + rk[12]*rk[13]*y[7] - rk[6]*y[9]*y[2] - rk[16]*y[9]*y[11] - rk[5]*y[9] - rk[8]*y[9]
   ydot[10] = rk[9]*y[7]
   ydot[11] = rk[5]*y[9] + rk[18]*rk[13]*y[4] + rk[19]*rk[13]*y[6] - rk[14]*y[1]*y[11] - rk[16]*y[9]*y[11]
#

   return ydot

#--------------------------------------------------------------------------------------------------
def solve_(lsodes):
#   print('lsodes.tout', lsodes.tout, 'lsodes.y:', lsodes.y)
   neq = len(lsodes.y)
