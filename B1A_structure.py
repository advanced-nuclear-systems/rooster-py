#--------------------------------------------------------------------------------------------------
class Structure:

    # constructor: self is a 'structure' object created in B1
    def __init__(self, reactor):

        # create objects

        # initialize state: a vector of unknowns
        self.state = []
        self.neq = len(self.state)

    # create right-hand side vector: self is a 'structure' object created in B1
    def calculate_rhs(self, reactor, t):
        rhs = []
        return rhs
