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

        # create control object
        self.control = Control(self)

        # list of objects to be solved
        self.solve = self.control.input['solve']
        
        # create objects
        self.solid = Solid(self)
        self.fluid = Fluid(self)
        self.neutron = Neutron(self)

        # write list of unknowns to y0
        y0 = []
        if 'fuelrod' in self.solve:
            for i in range(self.solid.nfuelrods):
                for j in range(self.solid.fuelrod[i].nfuelpellets):
                    for k in range(self.solid.fuelrod[i].fuelpellet[j].nr):
                        # fuel temperature
                        y0.append(self.solid.fuelrod[i].fuelpellet[j].temp[k])
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
                    for j in range(self.solid.fuelrod[i].nfuelpellets):
                        for k in range(self.solid.fuelrod[i].fuelpellet[j].nr):
                            # fuel temperature
                            self.solid.fuelrod[i].fuelpellet[j].temp[k] = y[indx]
                            indx += 1
                    for j in range(self.solid.fuelrod[i].nfuelpellets):
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
                    y0.append(self.neutron.pointkinetics.cdnp[i])
                    indx += 1

            self.control.evaluate(self, t)
            rhs = []
            rhs += self.solid.calculate_rhs(self, t)
            rhs += self.fluid.calculate_rhs(self, t)
            rhs += self.neutron.calculate_rhs(self, t)
            return rhs

        # solve the whole system of ODEs
        solver = ode(construct_rhs, jac = None).set_integrator('lsoda', method = 'bdf')
        t0 = self.control.input['t0']
        solver.set_initial_value(y0, t0)
        solver.set_integrator
        f = open('output', 'w')
        for t_dt in self.control.input['t_dt'] :
            tend = t_dt[0]
            dtout = t_dt[1]
            while solver.successful() and solver.t < tend:
                time = solver.t + dtout
                y = solver.integrate(time)
                #print('time: {0:12.5e}'.format(time))
                # write time and all unknowns to output file
                f.write('{0:12.5e} '.format(time))
                for i in range(len(y)):
                    f.write('{0:12.5e} '.format(y[i]))
                f.write('\n')
        f.close()