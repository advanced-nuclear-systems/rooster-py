#--------------------------------------------------------------------------------------------------
class Neutron:

    def __init__(self, reactor):
        self.pointkinetics = PointKinetics(self, reactor)
        self.spatialkinetics = SpatialKinetics(self, reactor)
        self.state = self.pointkinetics.state + self.spatialkinetics.state
        self.neq = len(self.state)

    def calculate_rhs(self, reactor, t):
        rhs = []
        rhs += self.pointkinetics.calculate_rhs(reactor, t)
        rhs += self.spatialkinetics.calculate_rhs(reactor, t)
        return rhs

#--------------------------------------------------------------------------------------------------
class PointKinetics:

    def __init__(self, solid, reactor):
        #index = reactor.solid.neq + reactor.flow.neq + 
        self.power = 0
        self.ndnp = len(reactor.control.input['betaeff'])
        self.cdnp = [0] * self.ndnp
        self.state = [self.power] + self.cdnp
        self.neq = len(self.state)
        self.state = []

    def calculate_rhs(self, reactor, t):
        rhs = []
        return rhs

#--------------------------------------------------------------------------------------------------
class SpatialKinetics:

    def __init__(self, solid, reactor):
        power = 0
        self.state = []
        self.neq = len(self.state)

    def calculate_rhs(self, reactor, t):
        rhs = []
        return rhs
