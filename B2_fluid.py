from scipy import linalg

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
    # from and to dictionaries
    f = {'id':[],'i':[]}
    t = {'id':[],'i':[]}

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
            # vector of initial pressures in pipe nodes
            self.p.append([p0]*self.pipennodes[i])
            # vector of initial temperatures in pipe nodes
            self.temp.append([temp0]*self.pipennodes[i])

        # vector of junction types
        self.juntype = reactor.control.input['junction']['type']
        # construct from and to dictionaries
        self.f['id'] = [self.pipename.index(reactor.control.input['junction']['from'][i]) for i in range(len(self.juntype))]
        self.t['id'] = [self.pipename.index(reactor.control.input['junction']['to'][i]) for i in range(len(self.juntype))]
        self.f['i'] = [self.pipennodes[self.f['id'][i]]-1 for i in range(len(self.juntype))]
        self.t['i'] = [0 for i in range(len(self.juntype))]
        # add internal junctions
        for i in range(self.npipe):
            self.juntype += ['internal']*(self.pipennodes[i] - 1)
            # append from and to dictionaries
            self.f['id'] += [i]*(self.pipennodes[i]-1)
            self.t['id'] += [i]*(self.pipennodes[i]-1)
            self.f['i'] += range(self.pipennodes[i]-1)
            self.t['i'] += range(1,self.pipennodes[i])
        # number of junctions
        self.njun = len(self.juntype)
        # number of independent junctions
        self.njuni = self.juntype.count('independent')
        # number of dependent junctions
        self.njund = self.juntype.count('dependent')

        # create and inverse a matrix A linking dependent and independent junctions
        A = [[0]*(self.njuni+self.njund) for i in range(self.njuni+self.njund)]
        i = 0
        for j in range(self.njuni+self.njund):
            if self.juntype[j] == 'independent':
                A[j][j] = 1
            elif self.juntype[j] == 'dependent':
                while self.pipetype[i] == 'freelevel':
                    i += 1
                for jj in range(self.njuni+self.njund):
                    if self.f['id'][jj] == i:
                        A[j][jj] = -1
                    if self.t['id'][jj] == i:
                        A[j][jj] = 1
                i += 1
        self.invA = linalg.inv(A)

        # initialize vector of flowrate in independent junctions
        mdoti = [0]*self.njuni

        # initialize state: a vector of unknowns
        self.state = mdoti
        self.neq = len(self.state)

    #----------------------------------------------------------------------------------------------
    # create right-hand side vector: self is a 'fluid' object created in B
    def calculate_rhs(self, reactor, t):

        if not self.calculate:
            rhs = []
            return rhs

        # READ VARIABLES
        # flowrate in independent junctions
        mdoti = self.state[0:self.njuni]

        # find flowrates in dependent junctions
        # first construct right hand side of system invA*mdot = b
        i = 0
        b = [0]*(self.njuni+self.njund)
        for j in range(self.njuni+self.njund):
            if self.juntype[j] == 'independent':
                b[j] = mdoti[i]
                i += 1
            elif self.juntype[j] == 'dependent':
                b[j] = 0
        # then multiply matrix by vector: invA*mdot = b
        mdot = self.invA.dot(b)
        
        rhs = [0.1]*self.njuni
        return rhs
