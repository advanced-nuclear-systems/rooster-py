from B3A_pointkinetics import PointKinetics
from B3B_spatialkinetics import SpatialKinetics

#--------------------------------------------------------------------------------------------------
class Neutron:

    # constructor: self is a 'neutron' object created in B
    def __init__(self, reactor):

        # create objects
        self.pointkinetics = PointKinetics(reactor)
        self.spatialkinetics = SpatialKinetics(reactor)

    # create right-hand side list: self is a 'neutron' object created in B
    def calculate_rhs(self, reactor, t):

        # construct right-hand side list
        rhs = []
        rhs += self.pointkinetics.calculate_rhs(reactor, t)
        rhs += self.spatialkinetics.calculate_rhs(reactor, t)
        return rhs
