from B1B0A_fuelgrain import FuelGrain

import math
import sys

#--------------------------------------------------------------------------------------------------
class FuelPellet:

    # fuel pellet Pu content array
    pu = []
    # fuel pellet burnup array
    b = []
    # fuel pellet porosity array
    por = []
    # fuel pellet temperature array
    temp = []
    # fuel pellet fuel types array
    type = []
    # fuel pellet deviation from stoechiometry array
    x = []

    # constructor: self is a 'fuelpellet' object created in B1B
    def __init__(self, i, reactor):

        # create object
        self.fuelgrain = FuelGrain(reactor)

        # INITIALIZATION
        # fuel pellet name
        self.name = reactor.control.input['pellet']['name'][i]
        # fuel pellet inner radius
        self.ri = reactor.control.input['pellet']['ri'][i]
        # fuel pellet outer radius
        self.ro = reactor.control.input['pellet']['ro'][i]
        # number of fuel pellet radial nodes
        self.nr = reactor.control.input['pellet']['nr'][i]
        # process fuel name
        fuelname = reactor.control.input['pellet']['name'][i]
        # find the fuel name in the vector of fuels
        try:
            ifuel = reactor.control.input['fuel']['name'].index(fuelname)
        except:
            print('****ERROR: input fuel name ' + fuelname + ' is not specified in the \'fuel\' card.')
            sys.exit()
        # fuel type of fuel pellet
        self.type = reactor.control.input['fuel']['type'][ifuel]
        # vector of Pu content in fuel pellet radial nodes
        self.pu = [[reactor.control.input['fuel']['pu'][ifuel]][i]]*self.nr
        # vector of fuel burnup in fuel pellet radial nodes
        self.b = [[reactor.control.input['fuel']['b'][ifuel]][i]]*self.nr
        # vector of deviation from stoechiometry in fuel pellet radial nodes
        self.x = [[reactor.control.input['fuel']['x'][ifuel]][i]]*self.nr
        # vector of porosity in fuel pellet radial nodes
        self.por = [[reactor.control.input['fuel']['por'][ifuel]][i]]*self.nr
        # vector of initial temperatures in fuel pellet radial nodes
        self.temp = [[reactor.control.input['fuel']['temp0'][ifuel]][i]]*self.nr

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
        self.state = self.fuelgrain.state + self.temp
        self.neq = len(self.state)

    # create right-hand side vector: self is a 'fuel' object created in B1B
    def calculate_rhs(self, i, reactor, t):
        # split vector of unknowns
        self.fuelgrain.state = self.state[0:self.fuelgrain.neq]
        k = self.fuelgrain.neq
        for j in range(self.nr):
            self.temp[j] = self.state[k]
            k += 1
        # construct right-hand side vector
        rhs = self.fuelgrain.calculate_rhs(reactor, t)

        # FUEL PROPERTIES:
        self.prop = []
        dict = {'rho':[], 'cp':[], 'k':[]}
        if self.type[i] == 'mox':
            for j in range(self.nr[i]):
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
        rhs += [0]*self.nr

        return rhs
