from B1B0_fuel import Fuel
from B1B1_innergas import InnerGas
from B1B2_clad import Clad

#--------------------------------------------------------------------------------------------------
class FuelRod:

    # number of axial layers
    nz = 0

    # constructor: self is a 'fuelrod' object created in B1 and indx is the index of this object in the list of fuelrods
    def __init__(self, indx, reactor):

        # INITIALIZATION
        # number of axial layers specified in input for fuel rod indx
        self.nz = len(reactor.control.input['fuelrod'][indx]['fuelid'])
        # create an object for every fuel z-layer
        self.fuel = []
        for i in range(self.nz):
            self.fuel.append(Fuel(i, indx, reactor))

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

        # construct right-hand side list
        rhs = []
        for i in range(self.nz):
            rhs += self.fuel[i].calculate_rhs(i, indx, reactor, t)
        for i in range(self.ncladzlayer):
            rhs += self.clad[i].calculate_rhs(i, indx, reactor, t)
        return rhs
