#--------------------------------------------------------------------------------------------------
class Solid:

    def __init__(self, reactor):
        self.structure = Structure(self, reactor)
        self.fuelrod = FuelRod(self, reactor)
        neq = self.structure.neq + self.fuelrod.neq

    def calculate_rhs(self, reactor, t, y):
        rhs = []
        return rhs

#--------------------------------------------------------------------------------------------------
class Structure:
    neq = 0
    def __init__(self, solid, reactor):
        pass

#--------------------------------------------------------------------------------------------------
class FuelRod:
    neq = 0
    def __init__(self, solid, reactor):
        pass
