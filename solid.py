#--------------------------------------------------------------------------------------------------
class Solid:

    def __init__(self, reactor):
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
        self.state = []
        self.neq = len(self.state)

    def calculate_rhs(self, reactor, t):
        rhs = []
        return rhs
