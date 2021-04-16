from B1B0A_fuelgrain import FuelGrain

import math
import sys

#--------------------------------------------------------------------------------------------------
class Fuel:

    # pellet Pu content array
    pu = []
    # pellet burnup array
    b = []
    # pellet porosity array
    por = []
    # pellet temperature array
    temp = []
    # pellet fuel types array
    type = []
    # pellet deviation from stoechiometry array
    x = []

    # constructor: self is a 'fuel' object created in B1B
    def __init__(self, reactor):

        # create object
        self.fuelgrain = FuelGrain(reactor)

        # INITIALIZATION
        # vector of pellet names
        self.pelletname = reactor.control.input['pellet']['name']
        # number of pellets
        self.npellet = len(self.pelletname)
        # vector of pellet inner radius
        self.ri = reactor.control.input['pellet']['ri']
        # vector of pellet outer radius
        self.ro = reactor.control.input['pellet']['ro']
        # vector of numbers of pellet radial nodes
        self.pelletnr = reactor.control.input['pellet']['nr']
        # process fuel names
        for i in range(self.npellet):
            fuel = reactor.control.input['pellet']['name'][i]
            # find the fuel name in the vector of fuels
            try:
                ifuel = reactor.control.input['fuel']['name'].index(fuel)
            except:
                print('****ERROR: input fuel name ' + fuel + ' is not specified in the \'fuel\' card.')
                sys.exit()
            type = reactor.control.input['fuel']['type'][ifuel]
            pu = reactor.control.input['fuel']['pu'][ifuel]
            b = reactor.control.input['fuel']['b'][ifuel]
            x = reactor.control.input['fuel']['x'][ifuel]
            por = reactor.control.input['fuel']['por'][ifuel]
            temp0 = reactor.control.input['fuel']['temp0'][ifuel]
            # vector of fuel types in pellet
            self.type.append(type)
            # vector of Pu content in pellet
            self.pu.append([pu]*self.pelletnr[i])
            # vector of fuel burnup in pellet
            self.b.append([b]*self.pelletnr[i])
            # vector of deviation from stoechiometry in pellet
            self.x.append([x]*self.pelletnr[i])
            # vector of porosity in pellet
            self.por.append([por]*self.pelletnr[i])
            # vector of initial temperatures in pellet nodes
            self.temp.append([temp0]*self.pelletnr[i])

        # GEOMETRY
        # mesh grid step
       # self.dr = (self.ro - self.ri)/(self.nr-1)
       # # vector of node radii (size = nr)
       # self.r = [i*self.dr for i in range(0, self.nr)]
       # # vector of node boundary radii (size = nr-1)
       # self.rb = [self.r[i]+self.dr/2 for i in range(0, self.nr-1)]
       # # vector of node volume (size = nr)
       # self.vol = [4/3 * self.rb[0]**3] + [4/3 * (self.rb[i]**3 - self.rb[i-1]**3) for i in range(1, self.nr-1)] + [4/3 * (self.r[self.nr-1]**3 - self.rb[self.nr-2]**3)]

        # initialize state: a vector of unknowns
        self.state = self.fuelgrain.state + [self.temp[i][j] for i in range(self.npellet) for j in range(self.pelletnr[i])]
        self.neq = len(self.state)

    # create right-hand side vector: self is a 'fuel' object created in B1B
    def calculate_rhs(self, reactor, t):
        # split vector of unknowns
        self.fuelgrain.state = self.state[0:self.fuelgrain.neq]
        k = self.fuelgrain.neq
        for i in range(self.npellet):
            for j in range(self.pelletnr[i]):
                self.temp[i][j] = self.state[k]
                k += 1
        # construct right-hand side vector
        rhs = self.fuelgrain.calculate_rhs(reactor, t)

        # FUEL PROPERTIES:
        self.prop = []
        for i in range(self.npellet):
            dict = {'rho':[], 'cp':[], 'k':[]}
            if self.type[i] == 'mox':
                for j in range(self.pelletnr[i]):
                    b = self.b[i][j]
                    por = self.por[i][j]
                    pu = self.pu[i][j]
                    t = self.temp[i][j]
                    x = self.x[i][j]
                    # density (kg/m3)
                    dict['rho'].append((11460*pu + 10960*(1 - pu)) * (1 - por))
                    # specific heat (J/kg-K), D.L. Hagrman, et al., "MATPRO-version 11", TREE-NUREG-1280, Rev 1, Idaho National Engineering Laboratory (1980).
                    dict['cp'].append(15.496*(19.53*539**2 * math.exp(539/t) / (t**2 * (math.exp(539/t) - 1)**2) + 2*9.25e-04*t + 6.02e06*4.01e4 / (1.987*t**2) * math.exp(-4.01e4/(1.987*t))))
                    # thermal conductivity (W/m-K), Y. Philipponneau, J. Nuclear Matter., 188 (1992) 194-197
                    dict['k'].append((1/( 1.528*math.sqrt(x+0.00931) - 0.1055 + 0.44*b + 2.855e-4*t ) + 76.38e-12*t**3) * (1-por)/(1+por)/0.864)
            self.prop.append(dict)

        # TIME DERIVATIVE OF FUEL TEMPERATURE:
        rhs += [0]*sum(self.pelletnr)

        return rhs
