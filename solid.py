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
        self.state = []
        self.neq = len(self.state)

    def calculate_rhs(self, reactor, t):
        rhs = []
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
