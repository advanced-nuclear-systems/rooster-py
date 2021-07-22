#--------------------------------------------------------------------------------------------------
# TREE OF CLASSES:
#     Reactor
#         Control
#         Solid
#             Structure
#             FuelRod
#                Fuel
#                   FuelGrain
#                InnerGas
#                Clad
#         Fluid
#         Core
#             Mix
#             Isotope
#         Data
#--------------------------------------------------------------------------------------------------
from B0_control import Control
from B4_data import Data
from B1_solid import Solid
from B2_fluid import Fluid
from B3_core import Core

# SciPy requires installation : python -m pip install --user numpy scipy matplotlib ipython jupyter pandas sympy nose
from scipy.integrate import ode

import time

#--------------------------------------------------------------------------------------------------
class Reactor:

    # constructor: self is a 'reactor' object created in A
    def __init__(self):

        # starting time
        self.tic0 = time.time()
        self.tic = self.tic0

        # create control object
        self.control = Control(self)

        # list of objects to be solved
        self.solve = self.control.input['solve']

        # create object fluid
        self.fluid = Fluid(self)

        # evaluate signals
        self.control.evaluate_signals(self, self.control.input['t0'])

        # create object solid
        self.solid = Solid(self)
        # create object core
        self.core = Core(self)
        # create object data
        self.data = Data(self)

        # write list of unknowns from self to y0
        y0 = self.control.write_to_y(self)

        #------------------------------------------------------------------------------------------
        # given t and y, function returns the list of the right-hand sides. called by the ODE solver
        def compose_rhs(t, y):

            # read list of unknowns from y to self
            self.control.read_from_y(self, y)

            # evaluate signals            
            self.control.evaluate_signals(self, t)

            # compose right-hand side vector
            rhs = []
            rhs += self.fluid.calculate_rhs(self, t)
            rhs += self.solid.compose_rhs(self, t)
            rhs += self.core.calculate_rhs(self, t)
            return rhs

        #------------------------------------------------------------------------------------------
        # prepare an output folder, copy input and open output files
        fid = self.control.open_output_files(self)
        t0 = self.control.input['t0']
        self.control.print_output_files(self, fid, t0, 0)

        # create ODE solver, initialize and set integrator
        solver = ode(compose_rhs, jac = None).set_integrator('lsoda', method = 'bdf')
        solver.set_initial_value(y0, t0)
        solver.set_integrator

        # main integration loop
        for t_dt in self.control.input['t_dt'] :
            tend = t_dt[0]
            dtout = t_dt[1]
            # solve the whole system of ODEs
            while solver.successful() and solver.t < tend:
                t = solver.t + dtout
                print(t)
                y = solver.integrate(t)

                # print to output files
                self.control.print_output_files(self, fid, t, 1)

        # close all output files
        for f in fid:
            f.close()

        tac = time.time()
        print('Wall time: ','{0:.3f}'.format(tac - self.tic0), ' s')
