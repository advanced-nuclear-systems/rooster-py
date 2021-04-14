from scipy import linalg

import math
import sys

#--------------------------------------------------------------------------------------------------
class Fluid:

    # flag defining if this class is included in calculations or not
    calculate = False
    # array of unknowns of this class
    state = []
    # number of unknowns/equations of this class   
    neq = 0
    # pipe coolant name array
    cool = []
    # pipe pressure array
    p = []
    # pipe temperature array
    temp = []
    # from and to dictionaries
    f = []
    t = []
    # pipe node index dictionaries
    indx = []

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
            type = reactor.control.input['coolant']['type'][icool]
            p0 = reactor.control.input['coolant']['p0'][icool]
            temp0 = reactor.control.input['coolant']['temp0'][icool]
            # vector of coolant names in pipe
            self.cool.append(type)
            # vector of initial pressures in pipe nodes
            self.p.append([p0]*self.pipennodes[i])
            # vector of initial temperatures in pipe nodes
            self.temp.append([temp0]*self.pipennodes[i])
        # assign index to every pipe node
        for i in range(self.npipe):
            for j in range(self.pipennodes[i]):
                self.indx.append((i,j))

        # vector of junction types
        self.juntype = reactor.control.input['junction']['type']
        # number of junctions
        self.njun = len(self.juntype)
        # number of independent junctions
        self.njuni = self.juntype.count('independent')
        # number of dependent junctions
        self.njund = self.juntype.count('dependent')

        # construct from and to arrays of tulips
        for j in range(self.njun):
            idf = reactor.control.input['junction']['from'][j]
            indx = self.pipename.index(idf)
            self.f.append((indx, self.pipennodes[indx]-1))
            idt = reactor.control.input['junction']['to'][j]
            indx = self.pipename.index(idt)
            self.t.append((indx, 0))
        # add internal junctions
        for i in range(self.npipe):
            self.juntype += ['internal']*(self.pipennodes[i] - 1)
            # append from and to vectors
            self.f += [(i,j) for j in range(self.pipennodes[i]-1)]
            self.t += [(i,j) for j in range(1,self.pipennodes[i])]
        self.njun = len(self.f)

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
                    if self.f[jj][0] == i:
                        A[j][jj] = -1
                    if self.t[jj][0] == i:
                        A[j][jj] = 1
                i += 1
        self.invA = linalg.inv(A)

        # create and inverse a matrix B of left-hand sides of momentum conservation equations (self.njun) and 
        # mass conservation equations differentiated w.r.t. time (self.npipe-self.npipef)
        n = self.njun + sum(self.pipennodes)
        B = [[0]*n for i in range(n)]
        for j in range(self.njun):
            B[j][j] = 1 # dmdot/dt

            i = self.njun + self.indx.index(self.f[j])
            B[j][i] = -1 # -P_from
            if self.pipetype[self.f[j][0]] != 'freelevel':
                B[i][j] = -1 # -dmdot/dt_out

            i = self.njun + self.indx.index(self.t[j])
            B[j][i] = 1 # +P_to
            if self.pipetype[self.t[j][0]] != 'freelevel':
                B[i][j] = 1 # +dmdot/dt_in
        for i in range(sum(self.pipennodes)):
            if self.pipetype[self.indx[i][0]] == 'freelevel':
                B[self.njun + i][self.njun + i] = 1 # P = +_freelevel
        self.invB = linalg.inv(B)

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

        # VARIABLES:
        # flowrate in independent junctions
        mdoti = self.state[0:self.njuni]

        # FLOWRATES IN DEPENDENT JUNCTIONS:
        # first construct right hand side of system invA*mdot = b
        i = 0
        b = [0]*(self.njuni+self.njund)
        for j in range(self.njuni+self.njund):
            if self.juntype[j] == 'independent':
                b[j] = mdoti[i]
                i += 1
            elif self.juntype[j] == 'dependent':
                b[j] = 0
        # then multiply matrix by vector: invA*mdot = b and convert to list
        mdot = self.invA.dot(b).tolist()
        # finally calculate flowrates in internal junctions
        for i in range(self.npipe):
            mdotpipe = 0
            for j in range(self.njuni+self.njund):
                if self.t[j][0] == i:
                    mdotpipe += mdot[j]
            for j in range(self.pipennodes[i]-1):
                mdot.append(mdotpipe)

        # FLUID PROPERTIES:
        self.prop = []
        for i in range(self.npipe):
            dict = {'rhol':[], 'visl':[], 'kl':[], 'cpl':[]}
            if self.cool[i] == 'na':
                for j in range(self.pipennodes[i]):
                    t = self.temp[i][j]
                    # J.K. Fink and L. Leibowitz "Thermodynamic and Transport Properties of Sodium Liquid and Vapor", ANL/RE-95/2, 1995, https://www.ne.anl.gov/eda/ANL-RE-95-2.pdf
                    dict['rhol'].append(219.0 + 275.32*(1.0 - t/2503.7) + 511.58*(1.0 - t/2503.7)**0.5)
                    dict['visl'].append(math.exp(-6.4406 - 0.3958*math.log(t) + 556.835/t)/dict['rhol'][j])
                    dict['kl'].append(124.67 - 0.11381*t + 5.5226e-5*t**2 - 1.1842e-8*t**3)
                    # Based on fit from J.K. Fink, etal."Properties for Reactor Safety Analysis", ANL-CEN-RSD-82-2, May 1982.
                    dict['cpl'].append(1646.97 - 0.831587*t + 4.31182e-04*t**2)
            self.prop.append(dict)

        # TIME DERIVATIVES OF FLOWRATES AND PRESSURES:
        # first construct right hand side of system invB*[mdot, P] = b
        b = [0]*(self.njun + sum(self.pipennodes))
        for j in range(self.njun):
            b[j] = 0
        for i in range(sum(self.pipennodes)):
            if self.pipetype[self.indx[i][0]] == 'freelevel': b[self.njun+i] = 1e5
        invBb = self.invB.dot(b).tolist()

        rhs = [0.1]*self.njuni
        return rhs
