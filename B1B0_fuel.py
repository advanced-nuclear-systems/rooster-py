from B1B0A_fuelgrain import FuelGrain

#--------------------------------------------------------------------------------------------------
class Fuel:

    # constructor: self is a 'fuel' object created in B1B
    def __init__(self, reactor):

        # create objects
        self.fuelgrain = FuelGrain(reactor)

        # initialize state: a vector of unknowns
        self.state = self.fuelgrain.state
        self.neq = len(self.state)

    # create right-hand side vector: self is a 'fuel' object created in B1B
    def calculate_rhs(self, reactor, t):
        # split vector of unknowns
        self.fuelgrain.state = self.state[0:self.fuelgrain.neq]
        # construct right-hand side vector
        rhs = []
        rhs += self.fuelgrain.calculate_rhs(reactor, t)
        return rhs
