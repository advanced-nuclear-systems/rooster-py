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

        # write list of unknowns to y0
        y0 = []
        if 'fuelrod' in self.solve:
            for i in range(self.solid.nfuelrods):
                for j in range(self.solid.fuelrod[i].nz):
                    for k in range(self.solid.fuelrod[i].fuel[j].nr):
                        if 'fuelgrain' in self.solve and i + j + k == 0: #i+j+k==0 is a temporal condition to solve fuel grain only for one node
                            # fuel grain monoatoms
                            for l in range(self.solid.fuelrod[i].fuel[j].fuelgrain[k].nr):
                                y0.append(self.solid.fuelrod[i].fuel[j].fuelgrain[k].c1[l])
                            # fuel grain bubble radii
                            for l in range(self.solid.fuelrod[i].fuel[j].fuelgrain[k].NB):
                                y0.append(self.solid.fuelrod[i].fuel[j].fuelgrain[k].ri[l])
                            # fuel grain fractional concentration of irradiation-induced uranium vacancies
                            for l in range(self.solid.fuelrod[i].fuel[j].fuelgrain[k].NB):
                                y0.append(self.solid.fuelrod[i].fuel[j].fuelgrain[k].cv_irr[l])
                            # fuel grain fractional concentration of irradiation-induced uranium interstitials
                            for l in range(self.solid.fuelrod[i].fuel[j].fuelgrain[k].NB):
                                y0.append(self.solid.fuelrod[i].fuel[j].fuelgrain[k].ci_irr[l])
                            # fuel grain fractional concentration of uranium vacancies ejected from intragranular as-fabricated pores
                            for l in range(self.solid.fuelrod[i].fuel[j].fuelgrain[k].NB):
                                y0.append(self.solid.fuelrod[i].fuel[j].fuelgrain[k].cv_p[l])
                            # fuel grain intragranular bubble concentation
                            for l in range(self.solid.fuelrod[i].fuel[j].fuelgrain[k].NB):
                                y0.append(self.solid.fuelrod[i].fuel[j].fuelgrain[k].bi[l])
                    for k in range(self.solid.fuelrod[i].fuel[j].nr):
                        # fuel temperature
                        y0.append(self.solid.fuelrod[i].fuel[j].temp[k])
                    for k in range(self.solid.fuelrod[i].clad[j].nr):
                        # clad temperature
                        y0.append(self.solid.fuelrod[i].clad[j].temp[k])
        if 'fluid' in self.solve:
            for j in range(self.fluid.njun):
                if self.fluid.juntype[j] == 'independent':
                    # flowrate in independent junctions
                    y0.append(self.fluid.mdoti[j])
            for i in range(self.fluid.npipe):
                for j in range(self.fluid.pipennodes[i]):
                    # temperature in pipe nodes
                    y0.append(self.fluid.temp[i][j])
        if 'pointkinetics' in self.solve:
            y0.append(self.core.power)
            for i in range(self.core.ndnp):
                y0.append(self.core.cdnp[i])

        #------------------------------------------------------------------------------------------
        # given t and y, function returns the list of the right-hand sides. called by the ODE solver
        def construct_rhs(t, y):
            # read list of unknowns from y
            indx = 0
            if 'fuelrod' in self.solve:
                for i in range(self.solid.nfuelrods):
                    for j in range(self.solid.fuelrod[i].nz):
                        for k in range(self.solid.fuelrod[i].fuel[j].nr):
                            if 'fuelgrain' in self.solve and i + j + k == 0:
                                for l in range(self.solid.fuelrod[i].fuel[j].fuelgrain[k].nr):
                                    # fuel grain monoatoms
                                    self.solid.fuelrod[i].fuel[j].fuelgrain[k].c1[l] = y[indx]
                                    indx += 1
                                for l in range(self.solid.fuelrod[i].fuel[j].fuelgrain[k].NB):
                                    # fuel grain bubble radii
                                    self.solid.fuelrod[i].fuel[j].fuelgrain[k].rb[l] = y[indx]
                                    indx += 1
                                for l in range(self.solid.fuelrod[i].fuel[j].fuelgrain[k].NB):
                                    # fuel grain fractional concentration of irradiation-induced uranium vacancies
                                    self.solid.fuelrod[i].fuel[j].fuelgrain[k].cv_irr[l] = y[indx]
                                    indx += 1
                                for l in range(self.solid.fuelrod[i].fuel[j].fuelgrain[k].NB):
                                    # fuel grain fractional concentration of irradiation-induced uranium interstitials
                                    self.solid.fuelrod[i].fuel[j].fuelgrain[k].ci_irr[l] = y[indx]
                                    indx += 1
                                for l in range(self.solid.fuelrod[i].fuel[j].fuelgrain[k].NB):
                                    # fuel grain fractional concentration of uranium vacancies ejected from intragranular as-fabricated pores
                                    self.solid.fuelrod[i].fuel[j].fuelgrain[k].cv_p[l] = y[indx]
                                    indx += 1
                                for l in range(self.solid.fuelrod[i].fuel[j].fuelgrain[k].NB):
                                    # fuel grain intragranular bubble concentrations
                                    self.solid.fuelrod[i].fuel[j].fuelgrain[k].bi[l] = y[indx]
                                    indx += 1
                        for k in range(self.solid.fuelrod[i].fuel[j].nr):
                            # fuel temperature
                            self.solid.fuelrod[i].fuel[j].temp[k] = y[indx]
                            indx += 1
                        for k in range(self.solid.fuelrod[i].clad[j].nr):
                            # clad temperature
                            self.solid.fuelrod[i].clad[j].temp[k] = y[indx]
                            indx += 1
            if 'fluid' in self.solve:
                for j in range(self.fluid.njun):
                    if self.fluid.juntype[j] == 'independent':
                        # flowrate in independent junctions
                        self.fluid.mdoti[j] = y[indx]
                        indx += 1
                for i in range(self.fluid.npipe):
                    for j in range(self.fluid.pipennodes[i]):
                        # temperature in pipe nodes
                        self.fluid.temp[i][j] = y[indx]
                        indx += 1
            if 'pointkinetics' in self.solve:
                self.core.power = y[indx]
                indx += 1
                for i in range(self.core.ndnp):
                    self.core.cdnp[i] = y[indx]
                    indx += 1

            self.control.evaluate(self, t)

            # signal-dependent junction
            if 'fluid' in self.solve:
                for j in range(self.fluid.njun):
                    if self.fluid.juntype[j] == 'independent' and self.fluid.junflowrate[j] != '':
                        # impose flowrate from the look-up table
                        self.fluid.mdoti[j] = self.control.signal[self.fluid.junflowrate[j]]

            # signal-dependent pipe
            if 'fluid' in self.solve:
                for i in range(self.fluid.npipe):
                    if self.fluid.pipetype[i] == 'normal' and self.fluid.signaltemp[i] != '':
                        # impose temperature from the look-up table
                        self.fluid.temp[i] = [self.control.signal[self.fluid.signaltemp[i]]] * self.fluid.pipennodes[i]

            rhs = []
            rhs += self.solid.calculate_rhs(self, t)
            rhs += self.fluid.calculate_rhs(self, t)
            rhs += self.core.calculate_rhs(self, t)
            return rhs

        #------------------------------------------------------------------------------------------
        # prepare an output folder, copy input and open output files
        fid = self.control.open_output_files(self)

        # solve the whole system of ODEs
        solver = ode(construct_rhs, jac = None).set_integrator('lsoda', method = 'bdf')
        t0 = self.control.input['t0']
        solver.set_initial_value(y0, t0)
        solver.set_integrator

        #------------------------------------------------------------------------------------------
        # main integration loop
        for t_dt in self.control.input['t_dt'] :
            tend = t_dt[0]
            dtout = t_dt[1]
            while solver.successful() and solver.t < tend:
                time = solver.t + dtout
                y = solver.integrate(time)

                #print('time: {0:12.5e}'.format(time))

                #----------------------------------------------------------------------------------
                # print output files
                indx = 0
                if 'fuelrod' in self.solve:
                    for i in range(self.solid.nfuelrods):
                        fid[indx].write('{0:12.5e} '.format(time) + ''.join(['{0:12.5e} '.format(self.solid.fuelrod[i].innergas.hgap[j]) for j in range(self.solid.fuelrod[i].nz)]) + '\n')
                        indx += 1
                        # fuel and clad temperatures
                        for j in range(self.solid.fuelrod[i].nz):
                            fid[indx].write('{0:12.5e} '.format(time) + ''.join(['{0:12.5e} '.format(self.solid.fuelrod[i].fuel[j].temp[k]) for k in range(self.solid.fuelrod[i].fuel[j].nr)]) + ''.join(['{0:12.5e} '.format(self.solid.fuelrod[i].clad[j].temp[k]) for k in range(self.solid.fuelrod[i].clad[j].nr)]) + '\n')
                            indx += 1
                            for k in range(self.solid.fuelrod[i].fuel[j].nr):
                                if 'fuelgrain' in self.solve and i + j + k == 0: 
                                    fid[indx].write('{0:12.5e} '.format(time) + ''.join(['{0:12.5e} '.format(self.solid.fuelrod[i].fuel[j].fuelgrain[k].c1[l]) for l in range(self.solid.fuelrod[i].fuel[j].fuelgrain[k].nr)]) + '\n')
                                    indx += 1
                                    fid[indx].write('{0:12.5e} '.format(time) + ''.join(['{0:12.5e} '.format(self.solid.fuelrod[i].fuel[j].fuelgrain[k].ri[l]) for l in range(self.solid.fuelrod[i].fuel[j].fuelgrain[k].NB)]) + '\n')
                                    indx += 1
                                    fid[indx].write('{0:12.5e} '.format(time) + ''.join(['{0:12.5e} '.format(self.solid.fuelrod[i].fuel[j].fuelgrain[k].cv_irr[l]) for l in range(self.solid.fuelrod[i].fuel[j].fuelgrain[k].NB)]) + '\n')
                                    indx += 1
                                    fid[indx].write('{0:12.5e} '.format(time) + ''.join(['{0:12.5e} '.format(self.solid.fuelrod[i].fuel[j].fuelgrain[k].ci_irr[l]) for l in range(self.solid.fuelrod[i].fuel[j].fuelgrain[k].NB)]) + '\n')
                                    indx += 1
                                    fid[indx].write('{0:12.5e} '.format(time) + ''.join(['{0:12.5e} '.format(self.solid.fuelrod[i].fuel[j].fuelgrain[k].cv_p[l]) for l in range(self.solid.fuelrod[i].fuel[j].fuelgrain[k].NB)]) + '\n')
                                    indx += 1
                                    fid[indx].write('{0:12.5e} '.format(time) + ''.join(['{0:12.5e} '.format(self.solid.fuelrod[i].fuel[j].fuelgrain[k].bi[l]) for l in range(self.solid.fuelrod[i].fuel[j].fuelgrain[k].NB)]) + '\n')
                                    indx += 1
                if 'fluid' in self.solve:
                    # flowrate in dependent and independent junctions (no internal junctions)
                    fid[indx].write('{0:12.5e} '.format(time) + ''.join(['{0:12.5e} '.format(self.fluid.mdot[i]) for i in range(self.fluid.njuni + self.fluid.njund)]) + '\n')
                    indx += 1
                    for i in range(self.fluid.npipe):
                        fid[indx].write('{0:12.5e} '.format(time) + ''.join(['{0:12.5e} '.format(self.fluid.p[i][j]) for j in range(self.fluid.pipennodes[i])]) + '\n')
                        indx += 1
                        fid[indx].write('{0:12.5e} '.format(time) + ''.join(['{0:12.5e} '.format(self.fluid.temp[i][j]) for j in range(self.fluid.pipennodes[i])]) + '\n')
                        indx += 1
                        fid[indx].write('{0:12.5e} '.format(time) + ''.join(['{0:12.5e} '.format(self.fluid.vel[i][j]) for j in range(self.fluid.pipennodes[i])]) + '\n')
                        indx += 1
                        fid[indx].write('{0:12.5e} '.format(time) + ''.join(['{0:12.5e} '.format(self.fluid.re[i][j]) for j in range(self.fluid.pipennodes[i])]) + '\n')
                        indx += 1
                        fid[indx].write('{0:12.5e} '.format(time) + ''.join(['{0:12.5e} '.format(self.fluid.pr[i][j]) for j in range(self.fluid.pipennodes[i])]) + '\n')
                        indx += 1
                        fid[indx].write('{0:12.5e} '.format(time) + ''.join(['{0:12.5e} '.format(self.fluid.pe[i][j]) for j in range(self.fluid.pipennodes[i])]) + '\n')
                        indx += 1
                if 'pointkinetics' in self.solve:
                    # point kinetics power
                    fid[indx].write('{0:12.5e} '.format(time) + '{0:12.5e} '.format(self.core.power) + '\n')
                    indx += 1
                    # point kinetics cdnp
                    fid[indx].write('{0:12.5e} '.format(time) + ''.join(['{0:12.5e} '.format(self.core.cdnp[i]) for i in range(self.core.ndnp)]) + '\n')
                    indx += 1

        # close all output files
        for f in fid:
            f.close()