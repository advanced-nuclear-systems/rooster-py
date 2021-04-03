#--------------------------------------------------------------------------------------------------
# TREE OF CLASSES:
#     Reactor
#         Solid
#             Structure
#             FuelElement
#                Fuel
#                Gap
#                Clad
#         Fluid
#         Neutron
#             PointKinetics
#             SpatialKinetics
#         Control
#             Detector
#             Controller
#--------------------------------------------------------------------------------------------------
from control import Control
from fluid import Fluid
from neutron import Neutron
from solid import Solid

# SciPy requires installation : 
#     python -m pip install --user numpy scipy matplotlib ipython jupyter pandas sympy nose
from scipy.integrate import ode

#--------------------------------------------------------------------------------------------------
class Reactor:
   def __init__(self):
       # create objects
       self.control = Control(self)
       self.solid = Solid(self)
       self.fluid = Fluid(self)
       self.neutron = Neutron(self)
       # initialize state: a vector of variables
       self.state = self.control.state + self.solid.state + self.fluid.state + self.neutron.state
       # solver of the whole system of ODEs
       def solve(self):
           # function returning the vector of the right-hand sides and called by the ODE solver
           def construct_rhs(t, y):
               self.state = y
               self.control.evaluate(self, t)
               rhs = []
               rhs += self.solid.calculate_rhs(self, t)
               rhs += self.fluid.calculate_rhs(self, t)
               rhs += self.neutron.calculate_rhs(self, t)
               rhs += self.control.calculate_rhs(self, t)
               print(time, rhs)
               return rhs
       
           solver = ode(construct_rhs, jac = None).set_integrator('lsoda', method = 'bdf')
           t0 = self.control.input['t0']
           solver.set_initial_value(self.state, t0)
           solver.set_integrator
           for t_dt in self.control.input['t_dt'] :
              tend = t_dt[0]
              dtout = t_dt[1]
              while solver.successful() and solver.t < tend:
                  time = solver.t + dtout
                  self.state = solver.integrate(time)
                  print(time, self.state)
       solve(self)

#--------------------------------------------------------------------------------------------------
# create and solve
reactor = Reactor()
