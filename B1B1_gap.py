#--------------------------------------------------------------------------------------------------
class Gap:

    # constructor: self is a 'gap' object created in B1B
    def __init__(self, reactor):

        # create objects

        # initialize state: a vector of unknowns
        self.state = []
        self.neq = len(self.state)
