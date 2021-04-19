from B1A_structure import Structure
from B1B_fuelrod import FuelRod

#--------------------------------------------------------------------------------------------------
class Solid:

    #----------------------------------------------------------------------------------------------
    # constructor: self is a 'solid' object created in B
    def __init__(self, reactor):

        # INITIALIZATION
        # number of fuel rods specified in input
        self.nfuelrods = len([x['id'] for x in reactor.control.input['fuelrod']])
        # create an object for every fuel rod
        self.fuelrod = []
        for i in range(self.nfuelrods):
            self.fuelrod.append(FuelRod(i, reactor))

        # create structure object
        self.structure = Structure(reactor)

    #----------------------------------------------------------------------------------------------
    # create right-hand side list: self is a 'solid' object created in B
    def calculate_rhs(self, reactor, t):

        # construct right-hand side list
        rhs = []
        for i in range(self.nfuelrods):
            rhs += self.fuelrod[i].calculate_rhs(i, reactor, t)
        rhs += self.structure.calculate_rhs(reactor, t)
        return rhs
