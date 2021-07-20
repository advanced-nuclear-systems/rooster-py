from B1A_heatstructure import HeatStructure
from B1B_fuelrod import FuelRod

#--------------------------------------------------------------------------------------------------
class Solid:

    #----------------------------------------------------------------------------------------------
    # constructor: self is a 'solid' object created in B
    def __init__(self, reactor):

        # INITIALIZATION
        if 'fuelrod' in reactor.solve:
            # number of fuel rods specified in input
            self.nfuelrods = len([x['id'] for x in reactor.control.input['fuelrod']])
            # create an object for every fuel rod
            self.fuelrod = []
            for i in range(self.nfuelrods):
                self.fuelrod.append(FuelRod(i, reactor))

        if 'htstr' in reactor.solve:
            # number of heat structures specified in input
            self.nhtstr = len([x['id'] for x in reactor.control.input['htstr']])
            # create an object for every heat structure
            self.htstr = []
            for i in range(self.nhtstr):
                self.htstr.append(HeatStructure(i, reactor))

    #----------------------------------------------------------------------------------------------
    # compose right-hand side list: self is a 'solid' object created in B
    def compose_rhs(self, reactor, t):

        # construct right-hand side list
        rhs = []
        if 'fuelrod' in reactor.solve:
            for i in range(self.nfuelrods):
                rhs += self.fuelrod[i].compose_rhs(i, reactor, t)

        if 'htstr' in reactor.solve:
            for i in range(self.nhtstr):
                rhs += self.htstr[i].calculate_rhs(reactor, t)
        return rhs
