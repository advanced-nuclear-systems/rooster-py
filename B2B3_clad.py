#--------------------------------------------------------------------------------------------------
class Clad:

    def __init__(self, fuelrod, reactor):
        self.state = []
        self.neq = len(self.state)

    def calculate_rhs(self, reactor, t):
        rhs = []
        return rhs
