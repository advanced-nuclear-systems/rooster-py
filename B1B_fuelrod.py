from B1B0_fuel import Fuel
from B1B1_gap import Gap
from B1B2_clad import Clad

#--------------------------------------------------------------------------------------------------
class FuelRod:

    # constructor: self is a 'fuelrod' object created in B1
    def __init__(self, reactor):

        # create objects
        self.fuel = Fuel(reactor)
        self.gap = Gap(reactor)
        self.clad = Clad(reactor)

        # initialize state: a vector of unknowns
        self.state = self.fuel.state + self.gap.state + self.clad.state
        self.neq = len(self.state)

    # create right-hand side vector: self is a 'fuelrod' object created in B1
    def calculate_rhs(self, reactor, t):
        # split vector of unknowns
        self.fuel.state = self.state[0:self.fuel.neq]
        self.clad.state = self.state[len(self.fuel.state):len(self.fuel.state)+self.clad.neq]
        # construct right-hand side vector
        rhs = []
        rhs += self.fuel.calculate_rhs(reactor, t)
        return rhs
