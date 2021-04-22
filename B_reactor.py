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

import datetime
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
        self.neutron = Neutron(self)

        # write list of unknowns to y0
        y0 = []
        if 'fuelrod' in self.solve:
            for i in range(self.solid.nfuelrods):
                for j in range(self.solid.fuelrod[i].nz):
                    for k in range(self.solid.fuelrod[i].fuel[j].nr):
                        if 'fuelgrain' in self.solve and i == 0 and j == 0 and k == 0:
                            # fuel grain monoatoms
                            for l in range(self.solid.fuelrod[i].fuel[j].fuelgrain[k].nr):
                                y0.append(self.solid.fuelrod[i].fuel[j].fuelgrain[k].c1[l])
                    for k in range(self.solid.fuelrod[i].fuel[j].nr):
                        # fuel temperature
                        y0.append(self.solid.fuelrod[i].fuel[j].temp[k])
                    for k in range(self.solid.fuelrod[i].clad[j].nr):
                        # clad temperature
                        y0.append(self.solid.fuelrod[i].clad[j].temp[k])
        if 'fluid' in self.solve:
            for i in range(self.fluid.njuni):
                # flowrate in independent junctions
                y0.append(self.fluid.mdoti[i])
        if 'pointkinetics' in self.solve:
            y0.append(self.neutron.pointkinetics.power)
            for i in range(self.neutron.pointkinetics.ndnp):
                y0.append(self.neutron.pointkinetics.cdnp[i])

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
                        for k in range(self.solid.fuelrod[i].fuel[j].nr):
                            # fuel temperature
                            self.solid.fuelrod[i].fuel[j].temp[k] = y[indx]
                            indx += 1
                        for k in range(self.solid.fuelrod[i].clad[j].nr):
                            # clad temperature
                            self.solid.fuelrod[i].clad[j].temp[k] = y[indx]
                            indx += 1
            if 'fluid' in self.solve:
                for i in range(self.fluid.njuni):
                    # flowrate in independent junctions
                    self.fluid.mdoti[i] = y[indx]
                    indx += 1
            if 'pointkinetics' in self.solve:
                self.neutron.pointkinetics.power = y[indx]
                indx += 1
                for i in range(self.neutron.pointkinetics.ndnp):
                    self.neutron.pointkinetics.cdnp[i] = y[indx]
                    indx += 1

            self.control.evaluate(self, t)
            rhs = []
            rhs += self.solid.calculate_rhs(self, t)
            rhs += self.fluid.calculate_rhs(self, t)
            rhs += self.neutron.calculate_rhs(self, t)
            return rhs

        # prepare an output folder
        path4results = 'output'
        if os.path.isfile(path4results): os.remove(path4results)
        if not os.path.isdir(path4results): os.mkdir(path4results)
        path4results += os.sep + str(datetime.datetime.now())[0:21].replace(' ','-').replace(':','-').replace('.','-')
        if os.path.isfile(path4results): os.remove(path4results)
        if not os.path.isdir(path4results): os.mkdir(path4results)

        # solve the whole system of ODEs
        solver = ode(construct_rhs, jac = None).set_integrator('lsoda', method = 'bdf')
        t0 = self.control.input['t0']
        solver.set_initial_value(y0, t0)
        solver.set_integrator

        # copy input and open output files to output folder
        shutil.copyfile('input', path4results + os.sep + 'input')
        # open files for output
        fid = []
        if 'fuelrod' in self.solve:
            for i in range(self.solid.nfuelrods):
                fid.append(open(path4results + os.sep + 'fuelrod-hgap-' + [x['id'] for x in self.control.input['fuelrod']][i] + '.dat', 'w'))
                fid[-1].write(' ' + 'time(s)'.ljust(13) + ''.join([('hgap-' + str(j).zfill(3)).ljust(13) for j in range(self.solid.fuelrod[i].nz)]) + '\n')
                for j in range(self.solid.fuelrod[i].nz):
                    fid.append(open(path4results + os.sep + 'fuelrod-temp-' + [x['id'] for x in self.control.input['fuelrod']][i] + '-' + str(j).zfill(3) + '.dat', 'w'))
                    fid[-1].write(' ' + 'time(s)'.ljust(13) + ''.join([('tempf-' + str(k).zfill(3) + '(K)').ljust(13) for k in range(self.solid.fuelrod[i].fuel[j].nr)]) + ''.join([('tempc-' + str(k).zfill(3) + '(K)').ljust(13) for k in range(self.solid.fuelrod[i].clad[j].nr)]) + '\n')
                    for k in range(self.solid.fuelrod[i].fuel[j].nr):
                        if 'fuelgrain' in self.solve and i + j + k == 0: 
                            fid.append(open(path4results + os.sep + 'fuelrod-c1-' + [x['id'] for x in self.control.input['fuelrod']][i] + '-' + str(j).zfill(3) + '-' + str(k).zfill(3) + '.dat', 'w'))
                            fid[-1].write(' ' + 'time(s)'.ljust(13) + ''.join([('c1-' + str(l).zfill(3)).ljust(13) for l in range(self.solid.fuelrod[i].fuel[j].fuelgrain[k].nr)]) + '\n')
        if 'fluid' in self.solve:
            fid.append(open(path4results + os.sep + 'fluid-mdot.dat', 'w'))
            fid[-1].write(' ' + 'time(s)'.ljust(13) + ''.join([(self.control.input['junction']['from'][j] +'-' + self.control.input['junction']['to'][j]).ljust(13) for j in range(self.fluid.njuni + self.fluid.njund)]) + '\n')
        if 'pointkinetics' in self.solve:
            fid.append(open(path4results + os.sep + 'pointkinetics-power.dat', 'w'))
            fid[-1].write(' ' + 'time(s)'.ljust(13) + 'power(-)\n')
            fid.append(open(path4results + os.sep + 'pointkinetics-cdnp.dat', 'w'))
            fid[-1].write(' ' + 'time(s)'.ljust(13) + ''.join([('cdnp-' + str(i)).ljust(13) for i in range(self.neutron.pointkinetics.ndnp)]) + '\n')

        # main integration loop
        for t_dt in self.control.input['t_dt'] :
            tend = t_dt[0]
            dtout = t_dt[1]
            while solver.successful() and solver.t < tend:
                time = solver.t + dtout
                y = solver.integrate(time)

                print('time: {0:12.5e}'.format(time))

                # print output
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
                if 'fluid' in self.solve:
                    # flowrate in dependent and independent junctions (no internal junctions)
                    fid[indx].write('{0:12.5e} '.format(time) + ''.join(['{0:12.5e} '.format(self.fluid.mdot[i]) for i in range(self.fluid.njuni + self.fluid.njund)]) + '\n')
                    indx += 1
                if 'pointkinetics' in self.solve:
                    # point kinetics power
                    fid[indx].write('{0:12.5e} '.format(time) + '{0:12.5e} '.format(self.neutron.pointkinetics.power) + '\n')
                    indx += 1
                    # point kinetics cdnp
                    fid[indx].write('{0:12.5e} '.format(time) + ''.join(['{0:12.5e} '.format(self.neutron.pointkinetics.cdnp[i]) for i in range(self.neutron.pointkinetics.ndnp)]) + '\n')
                    indx += 1

        # close all output files
        for f in fid:
            f.close()