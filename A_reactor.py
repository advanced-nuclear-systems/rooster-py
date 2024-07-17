#--------------------------------------------------------------------------------------------------
# TREE OF CLASSES:
#     Reactor
#         Control
#         Solid
#             Structure
#             FuelRod
#                Fuel
#                InnerGas
#                Clad
#         Fluid
#         Core
#             Mix
#             Isotope
#         Data
#--------------------------------------------------------------------------------------------------
from A0_control import Control
from A4_data import Data
from A1_solid import Solid
from A2_fluid import Fluid
from A3_core import Core

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
        # create object solid
        self.solid = Solid(self)

        # create object core
        self.core = Core(self)
        # create object data
        self.data = Data(self)

        # evaluate signals
        self.control.evaluate_signals(self, self.control.input['t0'])

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
        self.control.print_output_files(self, fid, t0)

        # create ODE solver, initialize and set integrator
        solver = ode(compose_rhs, jac = None).set_integrator('lsoda', method = 'bdf', rtol = self.control.input['tol'][0], atol = self.control.input['tol'][1])
        solver.set_initial_value(y0, t0)
        #solver.set_integrator

        # main integration loop
        dtout = 1e-6
        for tend in self.control.input['tend'] :
            while solver.successful() and solver.t < tend:
                t = min(tend, solver.t + dtout)
                y = solver.integrate(t)
                # next step recommended by the solver
                dtout = solver._integrator.rwork[11]
            
                # evaluate signals            
                self.control.evaluate_signals(self, t)
                # print to output files
                self.control.print_output_files(self, fid, t)

        # close all output files
        for f in fid:
            f.close()

        tac = time.time()
        print('Wall time: ','{0:.3f}'.format(tac - self.tic0), ' s')
