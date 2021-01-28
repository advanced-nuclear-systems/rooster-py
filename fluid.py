#--------------------------------------------------------------------------------------------------
class Fluid:

    neq = 0

    def __init__(self, reactor):
        pass

    def calculate_rhs(self, reactor, t):
        rhs = [2*t, 2]
        return rhs
