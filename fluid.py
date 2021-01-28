#--------------------------------------------------------------------------------------------------
class Fluid:

    neq = 0

    def __init__(self, reactor):
        self.state = [0, 0]
        self.neq = len(self.state)

    def calculate_rhs(self, reactor, t):
        rhs = [2*t, 2]
        return rhs
