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

    #----------------------------------------------------------------------------------------------
    # constructor: self is a 'fuelpellet' object created in B1B
    def __init__(self, indx, reactor):

        # create object
        self.fuelgrain = FuelGrain(reactor)

        # INITIALIZATION
        # fuel pellet material id
        self.matid = reactor.control.input['pellet'][indx]['matid']
        # fuel pellet inner radius
        self.ri = reactor.control.input['pellet'][indx]['ri']
        # fuel pellet outer radius
        self.ro = reactor.control.input['pellet'][indx]['ro']
        # number of fuel pellet radial nodes
        self.nr = reactor.control.input['pellet'][indx]['nr']
        # process material id
        matid = reactor.control.input['pellet'][indx]['matid']
        # find the material id in the vector of fuels
        try:
            ifuel = [x['id'] for x in reactor.control.input['mat']].index(matid)
        except:
            print('****ERROR: input material id ' + matid + ' in pellet is not specified in the \'mat\' card.')
            sys.exit()
        # fuel type of fuel pellet
        self.type = reactor.control.input['mat'][ifuel]['type']
        # vector of Pu content in fuel pellet radial nodes
        self.pu = [[reactor.control.input['mat'][ifuel]['pu']][indx]]*self.nr
        # vector of fuel burnup in fuel pellet radial nodes
        self.b = [[reactor.control.input['mat'][ifuel]['b']][indx]]*self.nr
        # vector of deviation from stoechiometry in fuel pellet radial nodes
        self.x = [[reactor.control.input['mat'][ifuel]['x']][indx]]*self.nr
        # vector of porosity in fuel pellet radial nodes
        self.por = [[reactor.control.input['mat'][ifuel]['por']][indx]]*self.nr
        # vector of initial temperatures in fuel pellet radial nodes
        self.temp = [[reactor.control.input['mat'][ifuel]['temp0']][indx]]*self.nr

        # GEOMETRY
        # height of fuel pellet
        self.dz = reactor.control.input['pellet'][indx]['dz']
        # mesh grid step
        self.dr = (self.ro - self.ri)/(self.nr-1)
        # vector of node radii (size = nr)
        self.r = [self.ri + i*self.dr for i in range(self.nr)]
        # vector of node boundary radii (size = nr-1)
        self.rb = [self.r[i]+self.dr/2 for i in range(self.nr-1)]
        # vector of node volume (size = nr)
        self.vol = [(self.rb[0]**2 - self.r[0]**2)*self.dz] + [(self.rb[i]**2 - self.rb[i-1]**2)*self.dz for i in range(1, self.nr-1)] + [(self.r[self.nr-1]**2 - self.rb[self.nr-2]**2)*self.dz]

        # initialize state: a vector of unknowns
        self.state = self.fuelgrain.state + self.temp
        self.neq = len(self.state)

    #----------------------------------------------------------------------------------------------
    # create right-hand side vector: self is a 'fuel' object created in B1B
    def calculate_rhs(self, reactor, t):
        # split vector of unknowns
        self.fuelgrain.state = self.state[0:self.fuelgrain.neq]
        k = self.fuelgrain.neq
        for j in range(self.nr):
            self.temp[j] = self.state[k]
            k += 1
        # construct right-hand side vector
        rhs = self.fuelgrain.calculate_rhs(reactor, t)

        # FUEL PROPERTIES:
        self.prop = {'rho':[], 'cp':[], 'k':[]}
        if self.type == 'mox':
            for j in range(self.nr):
                b = self.b[j]
                por = self.por[j]
                pu = self.pu[j]
                t = self.temp[j]
                x = self.x[j]
                # density (kg/m3)
                self.prop['rho'].append((11460*pu + 10960*(1 - pu)) * (1 - por))
                # specific heat (J/kg-K), D.L. Hagrman, et al., "MATPRO-version 11", TREE-NUREG-1280, Rev 1, Idaho National Engineering Laboratory (1980).
                self.prop['cp'].append(15.496*(19.53*539**2 * math.exp(539/t) / (t**2 * (math.exp(539/t) - 1)**2) + 2*9.25e-04*t + 6.02e06*4.01e4 / (1.987*t**2) * math.exp(-4.01e4/(1.987*t))))
                # thermal conductivity (W/m-K), Y. Philipponneau, J. Nuclear Matter., 188 (1992) 194-197
                self.prop['k'].append((1/( 1.528*math.sqrt(x+0.00931) - 0.1055 + 0.44*b + 2.855e-4*t ) + 76.38e-12*t**3) * (1-por)/(1+por)/0.864)

        # TIME DERIVATIVE OF FUEL TEMPERATURE:
        # thermal conductivity between nodes
        kb = [0.5*(self.prop['k'][i] + self.prop['k'][i+1]) for i in range(self.nr-1)]
        # vector of heat balance (W) at node boundaries: 2*rb*dz * kb * dT/dr (size = nr-1)
        Q = [2*self.rb[i]*self.dz*kb[i]*(self.temp[i] - self.temp[i+1])/self.dr for i in range(self.nr-1)]
        dTdt = [-Q[0]/(self.prop['rho'][0]*self.prop['cp'][0]*self.vol[0])] + [(Q[i-1] - Q[i])/(self.prop['rho'][i]*self.prop['cp'][i]*self.vol[i]) for i in range(1, self.nr-1)] + [Q[self.nr-2]/(self.prop['rho'][self.nr-1]*self.prop['cp'][self.nr-1]*self.vol[self.nr-1])]
        dTdt = [dTdt[i] + 1e2 for i in range(self.nr)]
        dTdt[self.nr-1] = 0
        rhs += dTdt

        return rhs
