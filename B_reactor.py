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
#--------------------------------------------------------------------------------------------------
from B0_control import Control
from B1_solid import Solid
from B2_fluid import Fluid
from B3_core import Core

# SciPy requires installation : python -m pip install --user numpy scipy matplotlib ipython jupyter pandas sympy nose
from scipy.integrate import ode

import shutil
import os

#--------------------------------------------------------------------------------------------------
class Reactor:

    # constructor: self is a 'reactor' object created in A
    def __init__(self):

        # create control object
        self.control = Control(self)

        # list of objects to be solved
        self.solve = self.control.input['solve']

        # create other objects
        self.solid = Solid(self)
        self.fluid = Fluid(self)
        self.core = Core(self)

        # evaluate signals
        self.control.evaluate(self, self.control.input['t0'])

        # write list of unknowns to y0
        y0 = self.control.write_to_y(self)

        #------------------------------------------------------------------------------------------
        # given t and y, function returns the list of the right-hand sides. called by the ODE solver
        def construct_rhs(t, y):

            # read list of unknowns from y
            self.control.read_from_y(self, y)

            # evaluate signals            
            self.control.evaluate(self, t)

            rhs = []
            rhs += self.solid.calculate_rhs(self, t)
            rhs += self.fluid.calculate_rhs(self, t)
            rhs += self.core.calculate_rhs(self, t)
            return rhs

        #------------------------------------------------------------------------------------------
        # prepare an output folder, copy input and open output files
        fid = self.control.open_output_files(self)

        # create ODE solver, initialize and set integrator
        solver = ode(construct_rhs, jac = None).set_integrator('lsoda', method = 'bdf')
        t0 = self.control.input['t0']
        solver.set_initial_value(y0, t0)
        solver.set_integrator

        # main integration loop
        for t_dt in self.control.input['t_dt'] :
            tend = t_dt[0]
            dtout = t_dt[1]
            # solve the whole system of ODEs
            while solver.successful() and solver.t < tend:
                time = solver.t + dtout
                y = solver.integrate(time)

                #print('time: {0:12.5e}'.format(time))

                # print to output files
                self.control.print_output_files(self, fid, time)

        # close all output files
        for f in fid:
            f.close()