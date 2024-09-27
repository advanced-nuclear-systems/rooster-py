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
        self.tprt = 0.0 # psi_ytchen: initialize the output time
        # create control object
        self.control = Control(self)

        # list of objects to be solved
        self.solve = self.control.input['solve']
        # create object data psi_ytchen: object 'data' should appear before object 'fluid'
        self.data = Data(self)
        # create object fluid
        self.fluid = Fluid(self)
        # create object solid
        self.solid = Solid(self)

        # create object core
        self.core = Core(self)
        

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
        dtout = 1e-6
        self.control.print_output_files(self, fid, t0, dtout)
        self.tprt = t0 # psi_ytchen: record the last time of output

        # create ODE solver, initialize and set integrator
        rtol_ode = self.control.input['tol'][0]
        atol_ode = self.control.input['tol'][1]
        dtmax_ode= self.control.input['tend'][1]
        solver = ode(compose_rhs, jac = None).set_integrator('lsoda', method = 'bdf', rtol = rtol_ode, atol = atol_ode, max_step = dtmax_ode)
        solver.set_initial_value(y0, t0)
        #solver.set_integrator

        # main integration loop
        tend = self.control.input['tend'][0]
        while solver.successful() and solver.t < 0.0:#tend:
            t = min(tend, solver.t + dtout)
            y = solver.integrate(t)
            # next step recommended by the solver
            dtout = solver._integrator.rwork[11]
            # evaluate signals            
            self.control.evaluate_signals(self, t)
            # + psi_ytchen: set a upper limit value 'dtmax' for time step size, dtmax is set by signal
            if 'DTMAXTAB' in self.control.signal:
                dtmax = self.control.signal['DTMAXTAB']
                dtout = min(dtout,dtmax)
            # - psi_ytchen:
            # print to output files under control of signal 'DTMAXTAB'
                if (t-self.tprt >= self.control.signal['DTMAXTAB']) or t == tend:
                    self.control.print_output_files(self, fid, t, dtout)
                    self.tprt = t
            else: # psi_ytchen: out put result every cycle
                self.control.print_output_files(self, fid, t, dtout)
                
        if solver.successful() and solver.t < tend: 
            y0 = y
            t0 = solver.t
            solver = ode(compose_rhs, jac = None).set_integrator('vode', method = 'bdf', rtol = min(rtol_ode,1e-4), atol = min(atol_ode,1e-4), max_step = min(1e-4, dtmax_ode))
            solver.set_initial_value(y0, t0) 
            while solver.successful() and solver.t < tend:
                t = min(tend, solver.t + dtout)
                y = solver.integrate(t)
                # next step recommended by the solver
                dtout = solver._integrator.rwork[11]
                # evaluate signals            
                self.control.evaluate_signals(self, t)
                # + psi_ytchen: set a upper limit value 'dtmax' for time step size, dtmax is set by signal
                if 'DTMAXTAB' in self.control.signal:
                    dtmax = self.control.signal['DTMAXTAB']
                    dtout = min(dtout,dtmax)
                # - psi_ytchen:
                # print to output files under control of signal 'DTMAXTAB'
                    if (t-self.tprt >= self.control.signal['DTMAXTAB']) or t == tend:
                        self.control.print_output_files(self, fid, t, dtout)
                        self.tprt = t
                else: # psi_ytchen: out put result every cycle
                    self.control.print_output_files(self, fid, t, dtout)
                
        # close all output files
        for f in fid:
            f.close()

        tac = time.time()
        print('Wall time: ','{0:.3f}'.format(tac - self.tic0), ' s')