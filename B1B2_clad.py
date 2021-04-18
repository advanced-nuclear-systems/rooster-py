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
        print(self.matid)

        # initialize state: a vector of unknowns
        self.state = []
        self.neq = len(self.state)

    # create right-hand side vector: self is a 'clad' object created in B1B
    def calculate_rhs(self, reactor, t):
        rhs = []
        return rhs
