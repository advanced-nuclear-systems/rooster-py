#--------------------------------------------------------------------------------------------------
class Fluid:

    neq = 0

    # constructor: self is a 'fluid' object created in B
    def __init__(self, reactor):

        # create objects

        # initialize state: a vector of unknowns
        self.state = []
        self.neq = len(self.state)

    # create right-hand side vector: self is a 'fluid' object created in B
    def calculate_rhs(self, reactor, t):
        rhs = []
        return rhs
