#--------------------------------------------------------------------------------------------------
class Solid:

    def __init__(self, reactor):
#       create objects
        self.structure = Structure(self, reactor)
        self.fuelrod = FuelRod(self, reactor)
        self.state = self.structure.state + self.fuelrod.state
        self.neq = len(self.state)

    def calculate_rhs(self, reactor, t):
        rhs = []
        rhs += self.structure.calculate_rhs(reactor, t)
        rhs += self.fuelrod.calculate_rhs(reactor, t)
        return rhs

#--------------------------------------------------------------------------------------------------
class Structure:

    def __init__(self, solid, reactor):
        self.state = []
        self.neq = len(self.state)

    def calculate_rhs(self, reactor, t):
        rhs = []
        return rhs

#--------------------------------------------------------------------------------------------------
class FuelRod:

    def __init__(self, solid, reactor):
        self.fuel = Fuel(self, reactor)
        self.gap = Gap(self, reactor)
        self.clad = Clad(self, reactor)
        self.state = self.fuel.state + self.gap.state + self.clad.state
        self.neq = len(self.state)

    def calculate_rhs(self, reactor, t):
        rhs = []
        rhs += self.fuel.calculate_rhs(reactor, t)
        return rhs

#--------------------------------------------------------------------------------------------------
class Fuel:

    def __init__(self, fuelrod, reactor):
        # grain diameter
        self.dgrain = reactor.control.input['dgrain']
        # number of nodes in the grain
        self.nr = reactor.control.input['nrgrain']
        # mesh grid step
        self.dr = self.dgrain/(self.nr-1)
        # vector of node radii (size = nr)
        self.r = [i*self.dr for i in range(0, self.nr)]
        # vector of node boundary radii (size = nr-1)
        self.rb = [self.r[i]+self.dr/2 for i in range(0, self.nr-1)]
        # vector of node volume (size = nr)
        self.vol = [4/3 * self.rb[0]**3] + [4/3 * (self.rb[i]**3 - self.rb[i-1]**3) for i in range(1, self.nr-1)] + [4/3 * (self.r[self.nr-1]**3 - self.rb[self.nr-2]**3)]
        # fission rate
        self.frate = reactor.control.input['frate']
        # vector of monoatom concentrations
        self.c1 = [0]*self.nr
        self.state = self.c1
        self.neq = len(self.state)

    def calculate_rhs(self, reactor, t):
        # read variables
        index_c1 = 0
        self.c1 = reactor.state[index_c1:index_c1+self.nr]
        # diffusion constant
        dg = 1e-3 # to be replaced by a correlation
        # vector of dg * rb**2 * dc1/dr (size = nr-1)
        q1 = [dg*self.rb[i]**2*(self.c1[i] - self.c1[i+1])/self.dr for i in range(0, self.nr-1)]
        # vector of time derivative of monoatom concentrations
        dc1dt = [- self.rb[0]**2*q1[0]/self.vol[0]] + [(self.rb[i-1]**2*q1[i-1] - self.rb[i]**2*q1[i])/self.vol[i] for i in range(1, self.nr-1)] + [self.rb[self.nr-2]**2*q1[self.nr-2]/self.vol[self.nr-1]]
        dc1dt = [dc1dt[i] + 0.31*self.frate for i in range(0, self.nr)]
        dc1dt[self.nr-1] = 0
        rhs = dc1dt
        return rhs

#--------------------------------------------------------------------------------------------------
class Gap:

    def __init__(self, fuelrod, reactor):
        self.state = []
        self.neq = len(self.state)

#--------------------------------------------------------------------------------------------------
class Clad:

    def __init__(self, fuelrod, reactor):
        self.state = []
        self.neq = len(self.state)

    def calculate_rhs(self, reactor, t):
        rhs = []
        return rhs
