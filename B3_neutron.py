from B3A_pointkinetics import PointKinetics
from B3B_spatialkinetics import SpatialKinetics

#--------------------------------------------------------------------------------------------------
class Neutron:

    # constructor: self is a 'neutron' object created in B
    def __init__(self, reactor):

        # create objects
        self.pointkinetics = PointKinetics(reactor)
        self.spatialkinetics = SpatialKinetics(reactor)

        # initialize state: a vector of unknowns
        self.state = self.pointkinetics.state + self.spatialkinetics.state
        self.neq = len(self.state)

    # create right-hand side vector: self is a 'neutron' object created in B
    def calculate_rhs(self, reactor, t):
        # split vector of unknowns
        self.pointkinetics.state = self.state[0:self.pointkinetics.neq]
        self.spatialkinetics.state = self.state[len(self.pointkinetics.state):len(self.pointkinetics.state)+self.spatialkinetics.neq]
        # construct right-hand side vector
        rhs = []
        rhs += self.pointkinetics.calculate_rhs(reactor, t)
        rhs += self.spatialkinetics.calculate_rhs(reactor, t)
        return rhs
