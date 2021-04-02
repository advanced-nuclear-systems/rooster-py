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
   lsodes.neq = len(lsodes.y)-1

#  time
   lsodes.t = input['t']
   if lsodes.t < 0 : sys.exit('***ERROR: lsodes : initial time t is negative')

#  output time
   lsodes.tout = input['tout']
   if lsodes.tout < 0 : sys.exit('***ERROR: lsodes : output time tout is negative')
   if lsodes.tout < lsodes.t : sys.exit('***ERROR: lsodes : output time tout is less than time t')

#  relative tolerance (scalar or vector)
   if isinstance(input['rtol'], list): #vector
      if len(input['rtol']) != lsodes.neq+1 : 
         sys.exit('***ERROR: lsodes : length of rtol not equal to number of equations neq')
      else:
         lsodes.rtol = input['rtol']
   else: #scalar
      lsodes.rtol = [input['rtol']] * (lsodes.neq+1)

#  absolute tolerance (scalar or vector)
   if isinstance(input['atol'], list): #vector
      if len(input['atol']) != lsodes.neq+1 : 
         sys.exit('***ERROR: lsodes : length of atol not equal to number of equations neq')
      else:
         lsodes.rtol = input['atol']
   else: #scalar
      lsodes.atol = [input['atol']] * (lsodes.neq+1)

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

#  current time
   lsodes.tn = lsodes.t

#  time step
   lsodes.h = 1.0
   
#  counter of steps
   lsodes.nst = 0

#  counter of non-zeros in jacobian
   lsodes.nnz = 0

#  number of rhs evaluations
   lsodes.nfe = 0

#  number of extra rhs evaluations needed for each jacobian evaluation
   lsodes.ngp = 0

#  number of nonzero elements in the strict lower triangle of the lu factorization used in the chord iteration
   lsodes.nzl = 0

#  number of nonzero elements in the strict lower triangle of the lu factorization used in the chord iteration
#  the total number of nonzeros in the factorization is therefore nzl + nzu + neq.
   lsodes.nzu = 0

#  the nordsieck history array, of size nyh by (nqcur+1), where nyh is the initial value of neq.  
#  for j = 0,1,...,nqcur, column j+1 of yh contains hcur**j/factorial(j) times the j-th derivative of 
#  the interpolating polynomial currently representing the solution, evaluated at t = tcur.
   lsodes.yh = [list(lsodes.y)]

#  initial call to rhs
   lsodes.f0 = rhs(lsodes.t, lsodes.y)

#  error weight vector
   lsodes.ewt = [0]
   for i in range(1,lsodes.neq+1):
      e = lsodes.rtol[i]*abs(lsodes.y[i]) + lsodes.atol[i]
      if e <= 0.0:
         sys.exit('***ERROR: lsodes : ewt is non-positive: ' + str(e))
      else:
         lsodes.ewt.append(1.0/e)
         fac = 1.0 + 1.0/(float(i)+1.0)
#        perturb y for structure determination.
         lsodes.y[i] += fac*math.copysign(e,lsodes.y[i])

#  the block computes jacobian structure from results of n + 1 calls to rhs
   lsodes.savf = rhs(lsodes.t, lsodes.y)
#  calculate ia and ja
   lsodes.ia = [0, 1]
   lsodes.ja = [0]
   k = 1
#  j is a column index
   for j in range(1,lsodes.neq+1):
#     ja is a vector of length nnz: ja(k) is an index of row (0 ... neq-1) the non-zero a[k] belongs to
      lsodes.ja.append(j) # diagonal element
      k += 1
      e = 1.0/lsodes.ewt[j]
      dyj = math.copysign(e,lsodes.y[j])
      yj = lsodes.y[j]
      lsodes.y[j] += dyj
      lsodes.ftem = rhs(lsodes.t, lsodes.y)
      lsodes.y[j] = yj
#     i is a row index
      for i in range(1,lsodes.neq+1):
         dq = (lsodes.ftem[i] - lsodes.savf[i])/dyj
         if abs(dq) > lsodes.seth and i != j:
            lsodes.ja.append(i) # non-diagonal element
            k += 1
#     ia is a vector of length neq+: ia(j) is index k (1 ... nnz) of first non-zero a[k] of column j. ia(neq+1) = nnz.
      lsodes.ia.append(k)
   lsodes.nnz = lsodes.ia[lsodes.neq]
#  restore y from yh
   lsodes.y = lsodes.yh[0]

#  the block constructs groupings of the column indices of the jacobian matrix (igp[neq] and jgp[ngrp+1])
   exit1 = False 
   ncol = 1
   igp = [0]
   jgp = [0]
   ng = 1
#  working array
   jdone = [0]*(lsodes.neq+1)
   while True :
      igp.append(ncol)
#     working array
      incl = [0]*(lsodes.neq+1)
#     j is column
      for j in range(1,lsodes.neq+1):
#        reject column j if it is already in a group
         if jdone[j] != 1 :
            kmin = lsodes.ia[j]
            kmax = lsodes.ia[j+1]
            for k in range(kmin,kmax):
#              reject column j if it overlaps any column already in this group
               i = lsodes.ja[k]
               if incl[i] == 1 :
                  exit1 = True                  
                  break
            if exit1:
               exit1 = False
            else:
#              accept column j into group ng
               jgp.append(j)
               ncol += 1
               jdone[j] = 1
               for k in range(kmin,kmax):
                  i = lsodes.ja[k]
                  incl[i] = 1
#     stop if this group is empty (grouping is complete).
      if ncol == igp[ng] : break
      ng += 1
   ngrp = ng - 1

#  the block finds a minimum degree ordering of the rows and columns of a matrix m stored in (ia,ja,a) format.
   tag = 0
#  initialize degrees, element lists, and degree lists
   mark = [0] + [1]*(lsodes.neq)
   l = [0]*(lsodes.neq + 1 + 2*lsodes.nnz)
   v = [0]*(lsodes.neq + 1 + 2*lsodes.nnz)
   head = [0]*(lsodes.neq+1)
   next = [0]*(lsodes.neq+1)
   last = list(range(lsodes.neq+1))
   sfs = lsodes.neq + 1
#  create nonzero structure for each nonzero entry a(vi,vj)
   for vi in range(1,lsodes.neq+1):
      jmin = lsodes.ia[vi]
      jmax = lsodes.ia[vi+1] - 1
      if jmin <= jmax :
         for j in range(jmin,jmax+1):
            vj = lsodes.ja[j]
#           if a(vi,vj) is in strict lower triangle check for previous occurrence of a(vj,vi)
            if vj < vi :
               lvk = vi
               kmax = mark[vi] - 1
               if kmax != 0 :
                  for k in range(1,kmax+1):
                     lvk = l[lvk]
                     if v[lvk] == vj : break
      
#           for unentered entries a(vi,vj) 
            elif vj > vi :
#              enter vj in element list for vi
               mark[vi] += 1
               v[sfs] = vj
               l[sfs] = l[vi]
               l[vi] = sfs
               sfs += 1
      
#              enter vi in element list for vj
               mark[vj] += 1
               v[sfs] = vi
               l[sfs] = l[vj]
               l[vj] = sfs
               sfs += 1

#  сreate degree lists and initialize mark vector
   for vi in range(1,lsodes.neq+1):
      dvi = mark[vi]
      next[vi] = head[dvi]
      head[dvi] = vi
      last[vi] = -dvi
      if next[vi] > 0 : 
         last[next[vi]] = vi
      mark[vi] = tag

   k = 0
   dmin = 1

   while k < lsodes.neq :

#     search for vertex of minimum degree
      while head[dmin] <= 0 :
         dmin += 1
#     remove vertex vk of minimum degree from degree list
      vk = head[dmin]
      head[dmin] = next[vk]
      if head[dmin] > 0 : last[head[dmin]] = -dmin
 
#     number vertex vk, adjust tag, and tag vk
      k += 1
      next[vk] = -k
      last[vk] = dmin - 1
      tag += last[vk]
      mark[vk] = tag

#     the block (mdm) forms element vk from uneliminated neighbors of vk
#     initialize tag2 and list of uneliminated neighbors
      tag2 = mark[vk]
      tail = vk
#     for each vertex/element vs/es in element list of vk
      ls = l[vk]
      s = ls
      while s != 0 :
         ls = l[s]
         vs = v[s]
         if next[vs] < 0 :
#           if vs is active element, then for each vertex vb in boundary list of element vs
            lb = l[vs]
            blpmax = last[vs]
            for blp in range(1,blpmax+1):
               b = lb
               lb = l[b]
               vb = v[b]
#              if vb is untagged vertex, then tag2 and append to list of uneliminated neighbors
               if mark[vb] < tag2 :
                  mark[vb] = tag2
                  l[tail] = b
                  tail = b
#           mark vs inactive
            mark[vs] = tag2
         else :
#           if vs is uneliminated vertex, then tag2 and append to list of uneliminated neighbors
            mark[vs] = tag2
            l[tail] = s
            tail = s
         s = ls
      
#     terminate list of uneliminated neighbors
      l[tail] = 0

#     the block (mdp) purges inactive elements and make mass elimination
#     initialize tag2
      tag2 = mark[vk]
      
#     for each vertex vi in vk
      li = vk
      ilpmax = last[vk]
      if ilpmax > 0 :
         for ilp in range(1,ilpmax+1):
            i = li
            li = l[i]
            vi = v[li]
            
#           remove vi from degree list
            if last[vi] > 0 :
               next[last[vi]] = next[vi]
            elif last[vi] < 0 :
               head[-last[vi]] = next[vi]
            if next[vi] > 0 : last[next[vi]] = last[vi]
      
#           remove inactive items from element list of vi
            ls = vi
            while True :
               s = ls
               ls = l[s]
               if ls == 0 : break
               vs = v[ls]
               if mark[vs] >= tag2 :
                  free = ls
                  l[s] = l[ls]
                  ls = s
            
#           if vi is interior vertex, then remove from list and eliminate
            lvi = l[vi]
            if lvi == 0 :
               l[i] = l[li]
               li = i
               
               k += 1
               next[vi] = -k
               last[vk] -= 1
            else :
#              else classify vertex vi
               if l[lvi] == 0 :
                  evi = v[lvi]
                  if next[evi] < 0 :
                     if mark[evi] < 0 :
#                       if vi is duplicate vertex, then mark as such and adjust overlap count for corresponding element
                        last[vi] = 0
                        mark[evi] -= 1
                     else :
#                       else if vi is prototype vertex, then mark as such, initialize overlap count for corresponding element, 
#                       and move vi to end of boundary list
                        last[vi] = evi
                        mark[evi] = -1
                        l[tail] = li
                        tail = li
                        l[i] = l[li]
                        li = i
               else :
#                 else mark vi to compute degree
                  last[vi] = -vk
            
#              insert vk in element list of vi
               v[free] = vk
               l[free] = l[vi]
               l[vi] = free
#     terminate boundary list
      l[tail] = 0

#     the block (mdu) updates degrees of uneliminated vertices in vk
#     initialize tag2
      tag2 = mark[vk] - last[vk]

#     for each vertex vi in vk
      i = vk
      ilpmax = last[vk]
      if ilpmax > 0 :
         for ilp in range(1,ilpmax+1):
            i = l[i]
            vi = v[i]
            if last[vi] < 0 :
#              if vi neither prototype nor duplicate vertex, then merge elements to compute degree
               tag2 += 1
               dvi = last[vk]
               
#              for each vertex/element ves in element list of vi
               s = l[vi]
               while True :
                  s = l[s]
                  if s == 0 :
#                    insert vi in appropriate degree list
                     next[vi] = head[dvi]
                     head[dvi] = vi
                     last[vi] = -dvi
                     if next[vi] > 0 : last[next[vi]] = vi
                     if dvi < dmin : dmin = dvi
                     break
                  ves = v[s]
                  if next[ves] >= 0 :
#                    if ves is uneliminated vertex, then tag2 and adjust degree
                     mark[ves] = tag2
                     dvi += 1
                  else :
#                    if ves is active element, then expand check for outmatched vertex
                     if mark[ves] < 0 :
#                       else if vi is outmatched vertex, then adjust overlaps but do not compute degree
                        last[vi] = 0
                        mark[ves] -= 1
                        while True :
                           s = l[s]
                           if s == 0 : break
                           ves = v[s]
                           if mark[ves] < 0 : mark[ves] -= 1
                        break
                     else :
#                       for each vertex vb in es
                        b = ves
                        blpmax = last[ves]
                        for blp in range(1,blpmax+1):
                           b = l[b]
                           vb = v[b]
#                          if vb is untagged, then tag2 and adjust degree
                           if mark[vb] < tag2 :
                              mark[vb] = tag2
                              dvi += 1
            elif last[vi] > 0 :
#              else if vi is prototype vertex, then calculate degree by inclusion/exclusion and reset overlap count
               evi = last[vi]
               dvi = last[vk] + last[evi] + mark[evi]
               mark[evi] = 0
     
#              insert vi in appropriate degree list
               next[vi] = head[dvi]
               head[dvi] = vi
               last[vi] = -dvi
               if next[vi] > 0 : last[next[vi]] = vi
               if dvi < dmin : dmin = dvi

#  generate inverse permutation from permutation
   for k in range(1,lsodes.neq+1):
      next[k] = -next[k]
      last[next[k]] = k
   print('mark ', mark)
   print('head ', head)
   print('last ', last)
   print('next ', next)

   p = [0]*(lsodes.neq+2)
   jar = [0]*len(lsodes.ja)
   for i in range(1,lsodes.neq+1):
      if last[i] != i :
#        reorders rows of a, leaving row order unchanged
         for k in range(1,lsodes.neq+1):
            jmin = lsodes.ia[k]
            jmax = lsodes.ia[k+1] - 1
            if jmin <= jmax :
               p[lsodes.neq+1] = lsodes.neq + 1
               for j in range(jmin,jmax+1):
                  newj = next[lsodes.ja[j]]
                  i = lsodes.neq + 1
                  while True :
                     if p[i] >= newj : break
                     i = p[i]
                  p[newj] = p[i]
                  p[i] = newj
                  jar[newj] = lsodes.ja[j]
               i = lsodes.neq + 1
               for j in range(jmin,jmax+1):
                  i = p[i]
                  lsodes.ja[j] = jar[i]
         break
# everything seems ok up to here (in nroc there is operation with a but a is always zero so I did not include this coding here)
# next is nsfc
#--------------------------------------------------------------------------------------------------
#   ---01--- ---02--- ---03--- ---04--- ---05--- ---06--- ---07--- ---08--- ---09--- ---10--- ---11--- ---12--- 
#01:   X        .        .        .        .        .        .        .        .        .        .        .
#02:   X        X        X        X        X        .        .        .        .        .        .        X
#03:   .        X        X        X        .        X        .        .        .        X        .        .
#04:   .        X        X        X        .        .        .        .        .        .        .        .
#05:   .        X        .        .        X        .        .        .        .        .        .        X
#06:   .        .        X        .        .        X        .        .        .        X        .        .
#07:   .        .        .        .        .        .        X        .        .        X        .        X
#08:   .        .        .        .        .        .        .        X        .        X        .        .
#09:   .        .        .        X        X        X        X        .        o        .        .        .
#10:   .        .        X        .        .        X        X        X        .        X        .        X
#11:   .        .        .        .        .        .        .        X        .        .        o        .
#12:   .        X        .        .        X        .        X        .        .        X        .        X
#--------------------------------------------------------------------------------------------------
def rhs(t, y):
   rk = [0.1, 10.0, 50.0, 2.5, 0.1, 10.0, 50.0, 2.5, 50.0, 5.0, 50.0, 50.0, 50.0, 30.0, 100.0, 2.5, 100.0, 2.5, 50.0, 50.0]
   ydot = [0.0]*len(y)
   ydot[1] = -rk[0]*y[1]
   ydot[2] = rk[0]*y[1] + rk[10]*rk[13]*y[4] + rk[18]*rk[13]*y[5] - rk[2]*y[2]*y[3] - rk[14]*y[2]*y[12] - rk[1]*y[2]
   ydot[3] = rk[1]*y[2] - rk[4]*y[3] - rk[2]*y[2]*y[3] - rk[6]*y[10]*y[3] + rk[10]*rk[13]*y[4] + rk[11]*rk[13]*y[6]
   ydot[4] = rk[2]*y[2]*y[3] - rk[10]*rk[13]*y[4] - rk[3]*y[4]
   ydot[5] = rk[14]*y[2]*y[12] - rk[18]*rk[13]*y[5] - rk[15]*y[5]
   ydot[6] = rk[6]*y[10]*y[3] - rk[11]*rk[13]*y[6] - rk[7]*y[6]
   ydot[7] = rk[16]*y[10]*y[12] - rk[19]*rk[13]*y[7] - rk[17]*y[7]
   ydot[8] = rk[8]*y[10] - rk[12]*rk[13]*y[8] - rk[9]*y[8]
   ydot[9] = rk[3]*y[4] + rk[15]*y[5] + rk[7]*y[6] + rk[17]*y[7]
   ydot[10] = rk[4]*y[3] + rk[11]*rk[13]*y[6] + rk[19]*rk[13]*y[7] + rk[12]*rk[13]*y[8] - rk[6]*y[10]*y[3] - rk[16]*y[10]*y[12] - rk[5]*y[10] - rk[8]*y[10]
   ydot[11] = rk[9]*y[8]
   ydot[12] = rk[5]*y[10] + rk[18]*rk[13]*y[5] + rk[19]*rk[13]*y[7] - rk[14]*y[2]*y[12] - rk[16]*y[10]*y[12]

   return ydot

#--------------------------------------------------------------------------------------------------
def solve_(lsodes):
#   print('lsodes.tout', lsodes.tout, 'lsodes.y:', lsodes.y)
   neq = len(lsodes.y)
