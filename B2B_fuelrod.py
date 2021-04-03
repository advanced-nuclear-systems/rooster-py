from B2B1_fuel import Fuel
from B2B2_gap import Gap
from B2B3_clad import Clad

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
