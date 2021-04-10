from B1A_structure import Structure
from B1B_fuelrod import FuelRod

#--------------------------------------------------------------------------------------------------
class Solid:

    # constructor: self is a 'solid' object created in B
    def __init__(self, reactor):

        # create objects
        self.structure = Structure(reactor)
        self.fuelrod = FuelRod(reactor)

        # initialize state: a vector of unknowns
        self.state = self.structure.state + self.fuelrod.state
        self.neq = len(self.state)

    # create right-hand side vector: self is a 'solid' object created in B
    def calculate_rhs(self, reactor, t):
        # split vector of unknowns
        self.structure.state = self.state[0:self.structure.neq]
        self.fuelrod.state = self.state[len(self.structure.state):len(self.structure.state)+self.fuelrod.neq]
        # construct right-hand side vector
        rhs = []
        rhs += self.structure.calculate_rhs(reactor, t)
        rhs += self.fuelrod.calculate_rhs(reactor, t)
        return rhs
