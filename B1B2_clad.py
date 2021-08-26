import sys

#--------------------------------------------------------------------------------------------------
class Clad:

    #----------------------------------------------------------------------------------------------
    # constructor: self is a 'clad' object created in B1B,
    # indx is the axial index of this object in the fuel rod with index indxfuelrod
    def __init__(self, indx, indxfuelrod, reactor):

        # INITIALIZATION
        # dictionary of the fuel rod to which the clad belongs
        dictfuelrod = reactor.control.input['fuelrod'][indxfuelrod]
        # current clad id
        cladid = dictfuelrod['cladid'][indx]
        # pitch-to-diameter ratio of fuel rod lattice
        self.p2d = dictfuelrod['p2d'][indx]
        # fuel rod multiplicity
        self.mltpl = dictfuelrod['mltpl'][indx]

        # list of clad dictionaries specified in input
        list = reactor.control.input['clad']
        # index of the current clad in the list of clad dictionaries
        i = [x['id'] for x in list].index(cladid)

        # clad inner radius
        self.ri = list[i]['ri']
        # clad outer radius
        self.ro = list[i]['ro']
        # number of clad radial nodes
        self.nr = list[i]['nr']

        # clad material id
        matid = list[i]['matid']
        # find the clad material id in the list of materials
        try:
            iclad = [x['id'] for x in reactor.control.input['mat']].index(matid)
        except:
            print('****ERROR: clad material id ' + matid + ' is not specified in the \'mat\' card of input.')
            sys.exit()
        # dictionary of material properties of the current clad
        mat = reactor.control.input['mat'][iclad]
        # material type of clad
        self.type = mat['type']
        # list of initial temperatures in clad radial nodes
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
    # create right-hand side list: self is a 'clad' object created in B1B
    # indx is the axial index of this object in the fuel rod with index indxfuelrod
    def calculate_rhs(self, indx, indxfuelrod, reactor, t):

        # CLAD PROPERTIES:
        self.prop = {'rho':[], 'cp':[], 'k':[]}
        for j in range(self.nr):
            t = self.temp[j]
            # call material property function
            pro = reactor.data.matpro( {'type':self.type, 't':self.temp[j]} )
            # density (kg/m3)
            self.prop['rho'].append(pro['rho'])
            # specific heat (J/kg-K)
            self.prop['cp'].append(pro['cp'])
            # thermal conductivity (W/m-K)
            self.prop['k'].append(pro['k'])

        # TIME DERIVATIVE OF CLAD TEMPERATURE:
        # fuel object
        fuel = reactor.solid.fuelrod[indxfuelrod].fuel[indx]
        # inner gas object
        innergas = reactor.solid.fuelrod[indxfuelrod].innergas
        # gap conductance list
        hgap = innergas.calculate_hgap(indxfuelrod, reactor, t)

        # clad thermal conductivity between nodes
        kb = [0.5*(self.prop['k'][i] + self.prop['k'][i+1]) for i in range(self.nr-1)]
        # heat flux (W/m**2) times heat transfer area per unit height divided by pi from fuel to clad 
        Q = [(fuel.ro + self.ri) * hgap[indx] * (fuel.temp[fuel.nr-1] - self.temp[0])]
        # list of heat flux (W/m**2) times heat transfer area per unit height divided by pi at node boundaries: 2*rb * kb * dT/dr (size = nr-1)
        Q += [2*self.rb[i]*kb[i]*(self.temp[i] - self.temp[i+1])/self.dr for i in range(self.nr-1)]

        # dictionary of the fuel rod to which the clad belongs
        dictfuelrod = reactor.control.input['fuelrod'][indxfuelrod]
        # pipe node indexes
        jpipe = (reactor.fluid.pipeid.index(dictfuelrod['pipeid'][indx]), dictfuelrod['pipenode'][indx]-1)
        fluid = {}
        fluid['t'] = reactor.fluid.temp[jpipe[0]][jpipe[1]]
        fluid['type'] = reactor.fluid.type[jpipe[0]]
        # call material property function
        pro = reactor.data.matpro( {'type':fluid['type'], 't':fluid['t']} )
        fluid['pe'] = abs(reactor.fluid.vel[jpipe[0]][jpipe[1]]) * reactor.fluid.dhyd[jpipe[0]] * pro['rhol'] * pro['cpl'] / pro['kl']
        fluid['nu'] = reactor.data.nu( {'pe':fluid['pe'], 'p2d':self.p2d} )
        # heat exchange coefficient
        fluid['hex'] = fluid['nu'] * pro['kl'] / reactor.fluid.dhyd[jpipe[0]]
        # heat flux (W/m**2) times heat transfer area per unit height divided by pi from clad to coolant
        Q += [2*self.ro * fluid['hex']*(self.temp[self.nr-1] - fluid['t'])]

        rhocpv = [self.prop['rho'][i]*self.prop['cp'][i]*self.vol[i] for i in range(self.nr)]
        dTdt = [(Q[i] - Q[i+1])/rhocpv[i] for i in range(self.nr)]
        rhs = dTdt

        return rhs
