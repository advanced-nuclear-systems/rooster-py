#--------------------------------------------------------------------------------------------------
class Fluid:

    neq = 0

    def __init__(self, reactor):
        self.state = []
        self.neq = len(self.state)

    def calculate_rhs(self, reactor, t):
        rhs = []
        return rhs
