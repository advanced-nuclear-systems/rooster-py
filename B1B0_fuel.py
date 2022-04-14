from B1B0A_fuelgrain import FuelGrain

import math
import sys

#--------------------------------------------------------------------------------------------------
class Fuel:

    #----------------------------------------------------------------------------------------------
    # constructor: self is a 'fuel' object created in B1B, 
    # indx is the axial index of this object in the fuel rod with index indxfuelrod
    def __init__(self, indx, indxfuelrod, dz, reactor):

        # INITIALIZATION
        # dictionary of the fuel rod to which the fuel belongs
        dictfuelrod = reactor.control.input['fuelrod'][indxfuelrod]
        # current fuel id
        fuelid = dictfuelrod['fuelid'][indx]
        # radial power peaking factor of fuel rod
        self.kr = dictfuelrod['kr'][indx]
        # axial power peaking factor of fuel
        self.kz = dictfuelrod['kz'][indx]

        # list of fuel dictionaries specified in input
        list = reactor.control.input['fuel']
        # index of the current fuel in the list of fuel dictionaries
        i = [x['id'] for x in list].index(fuelid)

        # fuel inner radius
        self.ri = list[i]['ri']
        # fuel outer radius
        self.ro = list[i]['ro']
        # number of fuel radial nodes
        self.nr = list[i]['nr']

        # fuel material id
        matid = list[i]['matid']
        # find the fuel material id in the list of materials
        try:
            ifuel = [x['id'] for x in reactor.control.input['mat']].index(matid)
        except:
            print('****ERROR: fuel material id ' + matid + ' is not specified in the \'mat\' card of input.')
            sys.exit()
        # dictionary of material properties of the current fuel
        mat = reactor.control.input['mat'][ifuel]
        # material type of fuel
        self.type = mat['type']
        # list of Pu content in fuel radial nodes
        self.pu = [mat['pu']]*self.nr
        # list of fuel burnup in fuel radial nodes
        self.b = [mat['b']]*self.nr
        # list of deviation from stoechiometry in fuel radial nodes
        self.x = [mat['x']]*self.nr
        # list of porosity in fuel radial nodes
        self.por = [mat['por']]*self.nr
        # list of initial temperatures in fuel radial nodes
        self.temp = [mat['temp0']]*self.nr

        # mesh grid step
        self.dr = (self.ro - self.ri)/(self.nr-1)
        # list of node radii (size = nr)
        self.r = [self.ri + i*self.dr for i in range(self.nr)]
        # list of node boundary radii (size = nr-1)
        self.rb = [self.r[i]+self.dr/2 for i in range(self.nr-1)]
        # list of node volume (size = nr)
        self.vol = [self.rb[0]**2 - self.r[0]**2] + [self.rb[i]**2 - self.rb[i-1]**2 for i in range(1, self.nr-1)] + [self.r[self.nr-1]**2 - self.rb[self.nr-2]**2]       
        if 'fuelgrain' in reactor.solve:
            # create an object fuel grain for every radial node of fuel
            self.fuelgrain = []
            for i in range(self.nr):
                self.fuelgrain.append(FuelGrain(i, indx, indxfuelrod, reactor))

    #----------------------------------------------------------------------------------------------
    # create right-hand side list: self is a 'fuel' object created in B1B
    # indx is the axial index of this object in the fuel rod with index indxfuelrod
    def calculate_rhs(self, indx, indxfuelrod, reactor, t):

        # construct right-hand side list
        rhs = []
        if 'fuelgrain' in reactor.solve and indx == 0 and indxfuelrod == 0:
            for i in range(self.nr):
                if i == 0:
                    rhs += self.fuelgrain[indx].calculate_rhs(reactor, t)

        # FUEL PROPERTIES:
        self.prop = {'rho':[], 'cp':[], 'k':[]}
        for j in range(self.nr):
            # call material property function
            pro = reactor.data.matpro( {'type':self.type, 't':self.temp[j], 'b':self.b[j], 'por':self.por[j], 'pu':self.pu[j], 'x':self.x[j]} )
            # density (kg/m3)
            self.prop['rho'].append(pro['rho'])
            # specific heat (J/kg-K)
            self.prop['cp'].append(pro['cp'])
            # thermal conductivity (W/m-K)
            self.prop['k'].append(pro['k'])

        # TIME DERIVATIVE OF FUEL TEMPERATURE:
        # inner gas object
        innergas = reactor.solid.fuelrod[indxfuelrod].innergas
        # gap conductance list
        hgap = innergas.calculate_hgap(indxfuelrod, reactor, t)
        # clad object
        clad = reactor.solid.fuelrod[indxfuelrod].clad[indx]

        # fuel thermal conductivity between nodes
        kb = [0.5*(self.prop['k'][i] + self.prop['k'][i+1]) for i in range(self.nr-1)]
        # heat flux (W/m**2) times heat transfer area per unit height at node boundaries: 2*rb * kb * dT/dr (size = nr-1)
        Q = [0] + [2*self.rb[i]*kb[i]*(self.temp[i] - self.temp[i+1])/self.dr for i in range(self.nr-1)]
        # add heat flux (W/m**2) times heat transfer area per unit height from fuel to clad 
        Q += [0.5*(self.ro + clad.ri) * hgap[indx] * (self.temp[self.nr-1] - clad.temp[0])]
        # power density
        if 'pointkinetics' in reactor.solve:
            qv = reactor.core.qv_average * self.kr * self.kz
        else:
            qv = reactor.control.input['power0']/sum(self.vol) * self.kr * self.kz
        rhocpv = [self.prop['rho'][i]*self.prop['cp'][i]*self.vol[i] for i in range(self.nr)]
        dTdt = [(Q[i] - Q[i+1] + qv*self.vol[i])/rhocpv[i] for i in range(self.nr)]
        rhs += dTdt

        return rhs
