from B1B0_fuelpellet import FuelPellet
from B1B1_innergas import InnerGas
from B1B2_clad import Clad

#--------------------------------------------------------------------------------------------------
class FuelRod:

    # flag defining if this class is included in calculations or not
    calculate = False
    # number of unknowns/equations of this class   
    neq = 0
    # number of fuel pellets
    nfuelpellets = 0
    # array of unknowns of this class
    state = []

    # constructor: self is a 'fuelrod' object created in B1
    def __init__(self, reactor):

        # check if this class is to be solved
        s = reactor.control.input['solve']
        self.calculate = any(['fuelrod' in s[i][0] for i in range(len(s))])
        if not self.calculate:
            return

        # INITIALIZATION
        # number of fuel pellets specified in input
        self.nfuelpellets = len(reactor.control.input['pellet'])
        # create an object for every fuel pellet
        self.fuelpellet = []
        for i in range(self.nfuelpellets):
            self.fuelpellet.append(FuelPellet(i, reactor))
        self.innergas = InnerGas(reactor)
        self.clad = Clad(reactor)

        # initialize state: a vector of unknowns
        self.state = []
        for i in range(self.nfuelpellets):
            self.state += self.fuelpellet[i].state 
        self.state += self.innergas.state 
        self.state += self.clad.state
        self.neq = len(self.state)

    # create right-hand side vector: self is a 'fuelrod' object created in B1
    def calculate_rhs(self, reactor, t):

        if not self.calculate:
            rhs = []
            return rhs

        # split vector of unknowns
        k = 0
        for i in range(self.nfuelpellets):
            self.fuelpellet[i].state = self.state[k:k+self.fuelpellet[i].neq]
            k += self.fuelpellet[i].neq
        self.clad.state = self.state[k:k+self.clad.neq]
        # construct right-hand side vector
        rhs = []
        for i in range(self.nfuelpellets):
            rhs += self.fuelpellet[i].calculate_rhs(reactor, t)
        return rhs
