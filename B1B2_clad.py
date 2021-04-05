#--------------------------------------------------------------------------------------------------
class Clad:

    # constructor: self is a 'clad' object created in B1B
    def __init__(self, reactor):

        # create objects

        # initialize state: a vector of unknowns
        self.state = []
        self.neq = len(self.state)

    # create right-hand side vector: self is a 'clad' object created in B1B
    def calculate_rhs(self, reactor, t):
        rhs = []
        return rhs
