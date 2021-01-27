#--------------------------------------------------------------------------------------------------
class Neutron:

    neq = 0

    def __init__(self, reactor):
        self.point_kinetics = PointKinetics(self, reactor)
        self.spatial_kinetics = SpatialKinetics(self, reactor)

    def calculate_rhs(self,reactor, t, y):
        rhs = []
        return rhs

#--------------------------------------------------------------------------------------------------
class PointKinetics:
    neq = 0
    def __init__(self, solid, reactor):
        pass

#--------------------------------------------------------------------------------------------------
class SpatialKinetics:
    neq = 0
    def __init__(self, solid, reactor):
        pass
