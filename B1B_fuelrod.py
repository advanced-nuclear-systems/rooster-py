from B1B0_fuel import Fuel
from B1B1_gap import Gap
from B1B2_clad import Clad

#--------------------------------------------------------------------------------------------------
class FuelRod:

    # flag defining if this class is included in calculations or not
    calculate = False
    # array of unknowns of this class
    state = []
    # number of unknowns/equations of this class   
    neq = 0

    # constructor: self is a 'fuelrod' object created in B1
    def __init__(self, reactor):

        # check if this class is to be solved
        s = reactor.control.input['solve']
        self.calculate = any(['fuelrod' in s[i][0] for i in range(len(s))])
        if not self.calculate:
            return

        # create objects
        self.fuel = Fuel(reactor)
        self.gap = Gap(reactor)
        self.clad = Clad(reactor)

        # initialize state: a vector of unknowns
        self.state = self.fuel.state + self.gap.state + self.clad.state
        self.neq = len(self.state)

    # create right-hand side vector: self is a 'fuelrod' object created in B1
    def calculate_rhs(self, reactor, t):

        if not self.calculate:
            rhs = []
            return rhs

        # split vector of unknowns
        self.fuel.state = self.state[0:self.fuel.neq]
        self.clad.state = self.state[len(self.fuel.state):len(self.fuel.state)+self.clad.neq]
        # construct right-hand side vector
        rhs = []
        rhs += self.fuel.calculate_rhs(reactor, t)
        return rhs
