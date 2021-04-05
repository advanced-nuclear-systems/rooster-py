#--------------------------------------------------------------------------------------------------
class SpatialKinetics:

    # constructor: self is a 'spatialkinetics' object created in B3
    def __init__(self, reactor):
        power = 0
        self.state = []
        self.neq = len(self.state)

    # create right-hand side vector: self is a 'spatialkinetics' object created in B3
    def calculate_rhs(self, reactor, t):
        rhs = []
        return rhs
