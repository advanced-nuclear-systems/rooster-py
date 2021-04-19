#--------------------------------------------------------------------------------------------------
class PointKinetics:

    # flag defining if this class is included in calculations or not
    calculate = False

    # constructor: self is a 'pointkinetics' object created in B3
    def __init__(self, reactor):

        self.power = 1
        self.ndnp = len(reactor.control.input['betaeff'])
        self.cdnp = [0] * self.ndnp
        for i in range(self.ndnp) :
            self.cdnp[i] = reactor.control.input['betaeff'][i]*self.power/(reactor.control.input['dnplmb'][i]*reactor.control.input['tlife'])

    # create right-hand side list: self is a 'pointkinetics' object created in B3
    def calculate_rhs(self, reactor, t):

        # read input parameters
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
