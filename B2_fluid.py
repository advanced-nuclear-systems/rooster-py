import sys

#--------------------------------------------------------------------------------------------------
class Fluid:

    # flag defining if this class is included in calculations or not
    calculate = False
    # array of unknowns of this class
    state = []
    # number of unknowns/equations of this class   
    neq = 0
    # pipe pressure array
    p = []
    # pipe temperature array
    temp = []

    # constructor: self is a 'fluid' object created in B
    def __init__(self, reactor):

        # check if this class is to be solved
        s = reactor.control.input['solve']
        self.calculate = any(['fluid' in s[i][0] for i in range(len(s))])
        if not self.calculate:
            return

        # INITIALIZATION
        # vector of pipe names
        self.pipename = reactor.control.input['pipe']['name']
        # vector of pipe types
        self.pipetype = reactor.control.input['pipe']['type']
        # number of pipes
        self.npipe = len(self.pipetype)
        # number of freelevel pipes
        self.npipef = self.pipetype.count('freelevel')
        # vector of pipe hydraulic diameters
        self.dhyd = reactor.control.input['pipe']['dhyd']
        # vector of pipe elevations
        self.elev = reactor.control.input['pipe']['elev']
        # vector of pipe length
        self.len = reactor.control.input['pipe']['len']
        for i in range(self.npipe):
            if self.len[i] == 0:
                self.len[i] = abs(self.elev[i])
        # vector of pipe flow area
        self.areaz = reactor.control.input['pipe']['areaz']
        # vector of numbers of pipe nodes
        self.pipennodes = reactor.control.input['pipe']['nnodes']
        # process coolant names
        for i in range(self.npipe):
            cool = reactor.control.input['pipe']['cool'][i]
            # find the coolant name in the vector of coolants
            try:
                icool = reactor.control.input['coolant']['name'].index(cool)
            except:
                print('****ERROR: input coolant name ' + cool + ' is not specified in the \'coolant\' card.')
                sys.exit()
            p0 = reactor.control.input['coolant']['p0'][icool]
            temp0 = reactor.control.input['coolant']['temp0'][icool]
            # vector of pressures in pipe nodes
            self.p.append([p0]*self.pipennodes[i])
            # vector of temperatures in pipe nodes
            self.temp.append([temp0]*self.pipennodes[i])

        # vector of junction types
        self.juntype = reactor.control.input['junction']['type']
        self.juntype += ['internal'] 
        # number of junctions
        self.njun = len(self.juntype)
        # number of independent junctions
        self.njuni = self.juntype.count('independent')

        # initialize state: a vector of unknowns
        self.state = []
        self.neq = len(self.state)

    # create right-hand side vector: self is a 'fluid' object created in B
    def calculate_rhs(self, reactor, t):
        rhs = []
        return rhs
