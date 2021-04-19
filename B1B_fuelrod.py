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

    # constructor: self is a 'fuelrod' object created in B1 and indx is the index of this object in the list of fuelrods
    def __init__(self, indx, reactor):

        # check if this class is to be solved
        s = reactor.control.input['solve']
        self.calculate = any(['fuelrod' in s[i][0] for i in range(len(s))])
        if not self.calculate:
            return

        # INITIALIZATION
        # number of fuel pellets specified in input for fuel rod indx
        self.nfuelpellets = len(reactor.control.input['fuelrod'][indx]['pelletid'])
        # create an object for every fuel pellet
        self.fuelpellet = []
        for i in range(self.nfuelpellets):
            self.fuelpellet.append(FuelPellet(i, indx, reactor))

        self.innergas = InnerGas(reactor)

        # number of clad axial layers specified in input for fuel rod indx
        self.ncladzlayer = len(reactor.control.input['fuelrod'][indx]['cladid'])
        # create an object for every cald axial layer
        self.clad = []
        for i in range(self.ncladzlayer):
            self.clad.append(Clad(i, indx, reactor))

    # create right-hand side list: self is a 'fuelrod' object created in B1,
    # indx is the fuel rod index
    def calculate_rhs(self, indx, reactor, t):

        if not self.calculate:
            rhs = []
            return rhs

        # construct right-hand side list
        rhs = []
        for i in range(self.nfuelpellets):
            rhs += self.fuelpellet[i].calculate_rhs(i, indx, reactor, t)
        for i in range(self.ncladzlayer):
            rhs += self.clad[i].calculate_rhs(i, indx, reactor, t)
        return rhs
