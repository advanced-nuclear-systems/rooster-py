from B1B0A_fuelgrain import FuelGrain

import math
import sys

#--------------------------------------------------------------------------------------------------
class FuelPellet:

    # burnup list
    b = []
    # porosity list
    por = []
    # Pu content list
    pu = []
    # temperature list
    temp = []
    # fuel types list
    type = []
    # deviation from stoechiometry list
    x = []

    #----------------------------------------------------------------------------------------------
    # constructor: self is a 'fuelpellet' object created in B1B, 
    # indx is the axial index of this object in the fuel rod with index indxfuelrod
    def __init__(self, indx, indxfuelrod, reactor):

        # create object for fuel grain (to be fixed)
        self.fuelgrain = FuelGrain(reactor)

        # INITIALIZATION
        # dictionary of the fuel rod to which the fuel pellet belongs
        dictfuelrod = reactor.control.input['fuelrod'][indxfuelrod]
        # current pellet id
        fuelpelletid = dictfuelrod['pelletid'][indx]
        # id of the pipe cooling the fuel pellet
        pipeid = dictfuelrod['pipeid'][indx]
        # list of pipe dictionaries
        pipelist = reactor.control.input['pipe']
        # index of the pipe in the list of pipe dictionaries
        self.indxpipe = [x['id'] for x in pipelist].index(pipeid)
        # current pellet height
        self.dz = abs(pipelist[self.indxpipe]['elev']) / pipelist[self.indxpipe]['nnodes']

        # list of fuel pellet dictionaries specified in input
        list = reactor.control.input['pellet']
        # index of the current fuel pellet in the list of fuel pellet dictionaries
        i = [x['id'] for x in list].index(fuelpelletid)

        # fuel pellet inner radius
        self.ri = list[i]['ri']
        # fuel pellet outer radius
        self.ro = list[i]['ro']
        # number of fuel pellet radial nodes
        self.nr = list[i]['nr']

        # fuel pellet material id
        self.matid = list[i]['matid']
        # find the fuel pellet material id in the list of materials
        try:
            ifuel = [x['id'] for x in reactor.control.input['mat']].index(self.matid)
        except:
            print('****ERROR: fuel pellet material id ' + self.matid + ' is not specified in the \'mat\' card of input.')
            sys.exit()
        # dictionary of material properties of the current fuel pellet
        mat = reactor.control.input['mat'][ifuel]
        # material type of fuel pellet
        self.type = mat['type']
        # list of Pu content in fuel pellet radial nodes
        self.pu = [mat['pu']]*self.nr
        # list of fuel burnup in fuel pellet radial nodes
        self.b = [mat['b']]*self.nr
        # list of deviation from stoechiometry in fuel pellet radial nodes
        self.x = [mat['x']]*self.nr
        # list of porosity in fuel pellet radial nodes
        self.por = [mat['por']]*self.nr
        # list of initial temperatures in fuel pellet radial nodes
        self.temp = [mat['temp0']]*self.nr

        # mesh grid step
        self.dr = (self.ro - self.ri)/(self.nr-1)
        # list of node radii (size = nr)
        self.r = [self.ri + i*self.dr for i in range(self.nr)]
        # list of node boundary radii (size = nr-1)
        self.rb = [self.r[i]+self.dr/2 for i in range(self.nr-1)]
        # list of node volume per unit height (size = nr)
        self.vol = [self.rb[0]**2 - self.r[0]**2] + [self.rb[i]**2 - self.rb[i-1]**2 for i in range(1, self.nr-1)] + [self.r[self.nr-1]**2 - self.rb[self.nr-2]**2]

    #----------------------------------------------------------------------------------------------
    # create right-hand side list: self is a 'fuelpellet' object created in B1B
    # indx is the axial index of this object in the fuel rod with index indxfuelrod
    def calculate_rhs(self, indx, indxfuelrod, reactor, t):

        # construct right-hand side list
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
        # fuel thermal conductivity between nodes
        kb = [0.5*(self.prop['k'][i] + self.prop['k'][i+1]) for i in range(self.nr-1)]
        # list of heat flux (W/m**2) times heat transfer area per unit height at node boundaries: 2*rb * kb * dT/dr (size = nr-1)
        Q = [0] + [2*self.rb[i]*kb[i]*(self.temp[i] - self.temp[i+1])/self.dr for i in range(self.nr-1)]
        # clad object
        clad = reactor.solid.fuelrod[indxfuelrod].clad[indx]
        # add heat flux (W/m**2) times heat transfer area per unit height from fuel to clad 
        Q += [(self.ro + clad.ri) * 1e4 * (self.temp[self.nr-1] - clad.temp[0])]
        rhocpv = [self.prop['rho'][i]*self.prop['cp'][i]*self.vol[i] for i in range(self.nr)]
        dTdt = [(Q[i] - Q[i+1] + 1e9*self.vol[i])/rhocpv[i] for i in range(self.nr)]
        rhs += dTdt

        return rhs
