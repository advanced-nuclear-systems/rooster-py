#--------------------------------------------------------------------------------------------------
class Solid:

    neq = 0

    def __init__(self, reactor):
        self.structure = Structure()
        self.fuelrod = FuelRod()
        pass

    def calculate_rhs(self,reactor, t, y):
        rhs = []
        return rhs

#--------------------------------------------------------------------------------------------------
class Structure:
    pass

#--------------------------------------------------------------------------------------------------
class FuelRod:
    pass
