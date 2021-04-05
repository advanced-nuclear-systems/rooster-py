#--------------------------------------------------------------------------------------------------
class FuelGrain:

    # constructor: self is a 'fuelgrain' object created in B1B0
    def __init__(self, reactor):

        # INPUT PARAMETERS
        # grain diameter
        self.dgrain = reactor.control.input['dgrain']
        # number of nodes in the grain
        self.nr = reactor.control.input['nrgrain']
        # fission rate
        self.frate = reactor.control.input['frate']

        # GEOMETRY
        # mesh grid step
        self.dr = 0.5*self.dgrain/(self.nr-1)
        # vector of node radii (size = nr)
        self.r = [i*self.dr for i in range(0, self.nr)]
        # vector of node boundary radii (size = nr-1)
        self.rb = [self.r[i]+self.dr/2 for i in range(0, self.nr-1)]
        # vector of node volume (size = nr)
        self.vol = [4/3 * self.rb[0]**3] + [4/3 * (self.rb[i]**3 - self.rb[i-1]**3) for i in range(1, self.nr-1)] + [4/3 * (self.r[self.nr-1]**3 - self.rb[self.nr-2]**3)]

        # INITIALIZATION
        # initialize vector of monoatom concentrations
        self.c1 = [0]*self.nr
        # vector of unknowns
        self.state = self.c1
        self.neq = len(self.state)

    # create right-hand side vector: self is a 'fuelgrain' object created in B1B0
    def calculate_rhs(self, reactor, t):

        # READ VARIABLES
        index_c1 = 0
        self.c1 = reactor.state[index_c1:index_c1+self.nr]

        # 1. INTRAGRANULAR PROCESSES
        # diffusion constant
        dg = 1e-3 # to be replaced by a correlation
        # 1.1. MONOATOMS C1
        # vector of rb**2 * dg * dc1/dr (size = nr-1)
        q1 = [self.rb[i]**2*dg*(self.c1[i] - self.c1[i+1])/self.dr for i in range(0, self.nr-1)]
        # vector of time derivative of monoatom concentrations
        dc1dt = [-q1[0]/self.vol[0]] + [(q1[i-1] - q1[i])/self.vol[i] for i in range(1, self.nr-1)] + [q1[self.nr-2]/self.vol[self.nr-1]]
        dc1dt = [dc1dt[i] + 0.31*self.frate for i in range(0, self.nr)]
        #dc1dt[self.nr-1] = 0
        rhs = dc1dt
        return rhs
