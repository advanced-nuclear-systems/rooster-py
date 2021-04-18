import sys

#--------------------------------------------------------------------------------------------------
class Clad:

    # constructor: self is a 'clad' object created in B1B,
    # indx is the axial index of this object in the fuel rod with index indxfuelrod
    def __init__(self, indx, indxfuelrod, reactor):

        # INITIALIZATION
        # dictionary of the fuel rod to which the clad belongs
        dictfuelrod = reactor.control.input['fuelrod'][indxfuelrod]
        # current clad id
        cladid = dictfuelrod['cladid'][indx]
        # id of the pipe cooling the clad
        pipeid = dictfuelrod['pipeid'][indx]
        # list of pipe dictionaries
        pipelist = reactor.control.input['pipe']
        # index of the pipe in the list of pipe dictionaries
        self.indxpipe = [x['id'] for x in pipelist].index(pipeid)
        # current clad height
        self.dz = abs(pipelist[self.indxpipe]['elev']) / pipelist[self.indxpipe]['nnodes']

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
        self.matid = list[i]['matid']
        # find the clad material id in the list of materials
        try:
            iclad = [x['id'] for x in reactor.control.input['mat']].index(self.matid)
        except:
            print('****ERROR: clad material id ' + self.matid + ' is not specified in the \'mat\' card of input.')
            sys.exit()
        # dictionary of material properties of the current clad
        mat = reactor.control.input['mat'][iclad]
        # material type of clad
        self.type = mat['type']
        # list of initial temperatures in clad radial nodes
        self.temp = [mat['temp0']]*self.nr

        # initialize state: a list of unknowns
        self.state = self.temp
        self.neq = len(self.state)

    # create right-hand side list: self is a 'clad' object created in B1B
    def calculate_rhs(self, reactor, t):
        # split list of unknowns
        k = 0
        for j in range(self.nr):
            self.temp[j] = self.state[k]
            k += 1

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
        #print(self.prop)
        rhs = [0]*self.neq
        return rhs
