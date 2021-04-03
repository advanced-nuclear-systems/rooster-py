from B2A_structure import Structure
from B2B_fuelrod import FuelRod

#--------------------------------------------------------------------------------------------------
class Solid:

    def __init__(self, reactor):
        # create objects
        self.structure = Structure(self, reactor)
        self.fuelrod = FuelRod(self, reactor)
        # initialize state: a vector of variables
        self.state = self.structure.state + self.fuelrod.state
        self.neq = len(self.state)

    def calculate_rhs(self, reactor, t):
        rhs = []
        rhs += self.structure.calculate_rhs(reactor, t)
        rhs += self.fuelrod.calculate_rhs(reactor, t)
        return rhs
