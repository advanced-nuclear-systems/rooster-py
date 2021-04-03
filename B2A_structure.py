#--------------------------------------------------------------------------------------------------
class Structure:

    def __init__(self, solid, reactor):
        self.state = []
        self.neq = len(self.state)

    def calculate_rhs(self, reactor, t):
        rhs = []
        return rhs
