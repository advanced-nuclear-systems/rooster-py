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
        if self.type == 'ss316':
            for j in range(self.nr):
                t = self.temp[j]
                # density (kg/m3): @300K equation from Leibowitz, et al, "Properties for LMFBR safety analysis", ANL-CEN-RSD-76-1 (1976), p.117
                self.prop['rho'].append(7954)
                # specific heat (J/kg-K): Leibowitz, et al, "Properties for LMFBR safety analysis", ANL-CEN-RSD-76-1 (1976), p.100. Note that 1 mol of SS316 = 10.165 kg (https://www.webqc.org/molecular-weight-of-SS316.html) and 1 cal = 4.184 J
                self.prop['cp'].append((6.181 + 1.788e-3*t)*10.165*4.184)
                # thermal conductivity (W/m-K): Leibowitz, et al, "Properties for LMFBR safety analysis", ANL-CEN-RSD-76-1 (1976), p.104.
                self.prop['k'].append(9.248 + 1.571e-2*t)

        # TIME DERIVATIVE OF CLAD TEMPERATURE:
        # fuel object
        fuel = reactor.solid.fuelrod[indxfuelrod].fuel[indx]
        # inner gas object
        innergas = reactor.solid.fuelrod[indxfuelrod].innergas
        # gap conductance list
        hgap = innergas.calculate_hgap(indxfuelrod, reactor, t)

        # clad thermal conductivity between nodes
        kb = [0.5*(self.prop['k'][i] + self.prop['k'][i+1]) for i in range(self.nr-1)]
        # heat flux (W/m**2) times heat transfer area per unit height from fuel to clad 
        Q = [(fuel.ro + self.ri) * hgap[indx] * (fuel.temp[fuel.nr-1] - self.temp[0])]
        # list of heat flux (W/m**2) times heat transfer area per unit height at node boundaries: 2*rb * kb * dT/dr (size = nr-1)
        Q += [2*self.rb[i]*kb[i]*(self.temp[i] - self.temp[i+1])/self.dr for i in range(self.nr-1)] + [0]
        rhocpv = [self.prop['rho'][i]*self.prop['cp'][i]*self.vol[i] for i in range(self.nr)]
        dTdt = [(Q[i] - Q[i+1])/rhocpv[i] for i in range(self.nr)]

        # dictionary of the fuel rod to which the clad belongs
        dictfuelrod = reactor.control.input['fuelrod'][indxfuelrod]
        # pipe node indexes
        ipipe = reactor.fluid.pipeid.index(dictfuelrod['pipeid'][indx])
        jpipe = dictfuelrod['pipenode'][indx]
        dTdt[self.nr-1] -= 1e3*(self.temp[self.nr-1] - reactor.fluid.temp[ipipe][jpipe])
        rhs = dTdt

        return rhs
