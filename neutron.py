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
        self.power = 1
        self.ndnp = len(reactor.control.input['betaeff'])
        self.cdnp = [0] * self.ndnp
        for i in range(self.ndnp) :
            self.cdnp[i] = reactor.control.input['betaeff'][i]*self.power/(reactor.control.input['dnplmb'][i]*reactor.control.input['tlife'])
        self.state = [self.power] + self.cdnp
        self.neq = len(self.state)

    def calculate_rhs(self, reactor, t):
        index_power = reactor.solid.neq + reactor.fluid.neq
        index_cdnp = index_power + 1
        self.power = reactor.state[index_power]
        self.cdnp = reactor.state[index_cdnp:index_cdnp+self.ndnp]
        rho = reactor.control.signal['RHO_INS']
        betaeff = reactor.control.input['betaeff']
        tlife = reactor.control.input['tlife']
        dnplmb = reactor.control.input['dnplmb']

        dpowerdt = self.power * (rho - sum(betaeff)) / tlife
        dcdnpdt = [0] * self.ndnp
        for i in range(self.ndnp) :
            dpowerdt += dnplmb[i]*self.cdnp[i]
            dcdnpdt[i] = betaeff[i]*self.power/tlife - dnplmb[i]*self.cdnp[i]
        rhs = [dpowerdt] + dcdnpdt
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
