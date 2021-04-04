#--------------------------------------------------------------------------------------------------
# TREE OF CLASSES:
#     Reactor
#         Control
#         Solid
#             Structure
#             FuelRod
#                Fuel
#                Gap
#                Clad
#         Fluid
#         Neutron
#             PointKinetics
#             SpatialKinetics
#--------------------------------------------------------------------------------------------------
from B1_control import Control
from B2_solid import Solid
from B3_fluid import Fluid
from B4_neutron import Neutron

# SciPy requires installation : python -m pip install --user numpy scipy matplotlib ipython jupyter pandas sympy nose
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
        
        # function returning the vector of the right-hand sides and called by the ODE solver
        def construct_rhs(t, y):
            self.state = y
            self.control.evaluate(self, t)
            rhs = []
            rhs += self.solid.calculate_rhs(self, t)
            rhs += self.fluid.calculate_rhs(self, t)
            rhs += self.neutron.calculate_rhs(self, t)
            rhs += self.control.calculate_rhs(self, t)
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