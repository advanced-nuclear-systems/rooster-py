#--------------------------------------------------------------------------------------------------
# TREE OF CLASSES:
#     Reactor
#         Control
#         Solid
#             Structure
#             FuelRod
#                FuelPellet
#                   FuelGrain
#                InnerGas
#                Clad
#         Fluid
#         Neutron
#             PointKinetics
#             SpatialKinetics
#--------------------------------------------------------------------------------------------------
from B0_control import Control
from B1_solid import Solid
from B2_fluid import Fluid
from B3_neutron import Neutron

# SciPy requires installation : python -m pip install --user numpy scipy matplotlib ipython jupyter pandas sympy nose
from scipy.integrate import ode

#--------------------------------------------------------------------------------------------------
class Reactor:

    # constructor: self is a 'reactor' object created in A
    def __init__(self):
        # create objects
        self.control = Control(self)
        self.solid = Solid(self)
        self.fluid = Fluid(self)
        self.neutron = Neutron(self)
        # initialize state: a list of variables
        self.state = self.control.state + self.solid.state + self.fluid.state + self.neutron.state
        
        # function returning the list of the right-hand sides and called by the ODE solver
        def construct_rhs(t, y):
            # read list of unknowns and split it
            self.state = y
            self.control.state = y[0:self.control.neq]
            self.solid.state = y[len(self.control.state):len(self.control.state)+self.solid.neq]
            self.fluid.state = y[len(self.solid.state):len(self.solid.state)+self.fluid.neq]
            self.neutron.state = y[len(self.fluid.state):len(self.fluid.state)+self.neutron.neq]

            self.control.evaluate(self, t)
            rhs = []
            rhs += self.control.calculate_rhs(self, t)
            rhs += self.solid.calculate_rhs(self, t)
            rhs += self.fluid.calculate_rhs(self, t)
            rhs += self.neutron.calculate_rhs(self, t)
            return rhs
        # solve the whole system of ODEs
        solver = ode(construct_rhs, jac = None).set_integrator('lsoda', method = 'bdf')
        t0 = self.control.input['t0']
        solver.set_initial_value(self.state, t0)
        solver.set_integrator
        f = open('output', 'w')
        for t_dt in self.control.input['t_dt'] :
            tend = t_dt[0]
            dtout = t_dt[1]
            while solver.successful() and solver.t < tend:
                time = solver.t + dtout
                self.state = solver.integrate(time)
                print('time: {0:12.5e}'.format(time))
                # write time and all unknowns to output file
                f.write('{0:12.5e} '.format(time))
                for i in range(len(self.state)):
                    f.write('{0:12.5e} '.format(self.state[i]))
                f.write('\n')
        f.close()