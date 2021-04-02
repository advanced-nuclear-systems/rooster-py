from lsodes import Lsodes

input = {}

#number of equations
neq = 12
input['y'] = [0.0]*(neq+1)

#initial values of unknowns
input['y'][1] = 1.0

#initial time
input['t'] = 0.0

#output time
input['tout'] = 0.1

#relative tolerance (scalar or vector)
input['rtol'] = 1.0e-4

#absolute tolerance (scalar or vector)
input['atol'] = 1.0e-6

#output time flag
input['itask'] = 1

#istate flag (1: first call, 2: not first call, 3: not first call with changed input parameters)
input['istate'] = 1

#optional input flag
input['iopt'] = 1

lsodes = Lsodes(input)

for i in range(0,5):
    lsodes.solve()
    lsodes.tout *= 10.0
