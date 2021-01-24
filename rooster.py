from input import input
from scipy.integrate import ode # SciPy requires installation : python -m pip install --user numpy scipy matplotlib ipython jupyter pandas sympy nose

import json

#--------------------------------------------------------------------------------------------------
# CLASSES:
#     Reactor
#         Solid
#             FuelElement
#             Structure
#         Fluid
#         Neutron
#         Control
#             Detector
#             Controller
# FUNCTIONS:
#     read_input
#--------------------------------------------------------------------------------------------------
class Reactor:
    def __init__(self):
        self.inp = input()
        self.solid = Solid(self)
        self.fluid = Fluid(self)
        self.neutron = Neutron(self)
        solve(self)
#--------------------------------------------------------------------------------------------------
class Solid:
    x = 12
    def __init__(self, reactor):
        pass
#--------------------------------------------------------------------------------------------------
class Fluid:
    def __init__(self, reactor):
        pass
#--------------------------------------------------------------------------------------------------
class Neutron:
    def __init__(self, reactor):
        pass
#--------------------------------------------------------------------------------------------------
def read_input():
    #read the whole file as a string
    f = open('input.py', mode = 'r')
    str0 = f.read()
    f.close()

    str0 = str0.replace("'",'"').replace('\r\n','\n') #replace all single quotes by double ones and Windows EOL by Unix EOL
    #remove comments
    str = ''
    comment = False
    for c in str0:
        if c == '#':
            comment = True
        elif c == '\n':
            comment = False
        if not comment:
            str += c
    str = str.replace('\n','') #remove EOL
#    print(repr(str))
    input = json.loads(str) #convert from string (json) to dictionary
    return input

#--------------------------------------------------------------------------------------------------
def solve(reactor):
    r = ode(rhs, jac=None).set_integrator('lsoda', method='bdf')
    t0 = reactor.inp['t0']
    y0 = [0, 0]
    r.set_initial_value(y0, t0)
    dtout = reactor.inp['dtout']
    tend = reactor.inp['tend']
    while r.successful() and r.t < tend:
        print(r.t+dtout, r.integrate(r.t+dtout))
    return 0

def rhs(t, y):
    return [2*t, 2] 

#--------------------------------------------------------------------------------------------------
r = Reactor()
    