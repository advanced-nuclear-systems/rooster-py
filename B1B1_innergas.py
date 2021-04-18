#--------------------------------------------------------------------------------------------------
class InnerGas:

    # constructor: self is a 'innergas' object created in B1B
    def __init__(self, reactor):

        # create objects

        # initialize state: a list of unknowns
        self.state = []
        self.neq = len(self.state)
