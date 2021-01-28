#--------------------------------------------------------------------------------------------------
class Solid:

    def __init__(self, reactor):
        self.structure = Structure(self, reactor)
        self.fuelrod = FuelRod(self, reactor)
        self.neq = self.structure.neq + self.fuelrod.neq

    def calculate_rhs(self, reactor, t):
        rhs = []
        return rhs

#--------------------------------------------------------------------------------------------------
class Structure:

    def __init__(self, solid, reactor):
        self.neq = 0
        pass

#--------------------------------------------------------------------------------------------------
class FuelRod:

    def __init__(self, solid, reactor):
        self.neq = 0
