#--------------------------------------------------------------------------------------------------
class Fluid:

    def __init__(self, reactor):
        pass

    def calculate_rhs(self,reactor, t, y):
        rhs = [2*t, 2]
        return rhs
