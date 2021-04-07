#--------------------------------------------------------------------------------------------------
class PointKinetics:

    # flag defining if this class is included in calculations or not
    calculate = False
    # array of unknowns of this class
    state = []
    # number of unknowns/equations of this class   
    neq = 0

    # constructor: self is a 'pointkinetics' object created in B3
    def __init__(self, reactor):

        s = reactor.control.input['solve']
        self.calculate = any(['pointkinetics' in s[i][0] for i in range(len(s))])
        if not self.calculate:
            return

        self.power = 1
        self.ndnp = len(reactor.control.input['betaeff'])
        self.cdnp = [0] * self.ndnp
        for i in range(self.ndnp) :
            self.cdnp[i] = reactor.control.input['betaeff'][i]*self.power/(reactor.control.input['dnplmb'][i]*reactor.control.input['tlife'])
        self.state = [self.power] + self.cdnp
        self.neq = len(self.state)

    # create right-hand side vector: self is a 'pointkinetics' object created in B3
    def calculate_rhs(self, reactor, t):

        if not self.calculate:
            rhs = []
            return rhs

        # read variables
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
