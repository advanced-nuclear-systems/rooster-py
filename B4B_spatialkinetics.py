#--------------------------------------------------------------------------------------------------
class SpatialKinetics:

    def __init__(self, solid, reactor):
        power = 0
        self.state = []
        self.neq = len(self.state)

    def calculate_rhs(self, reactor, t):
        rhs = []
        return rhs
