#--------------------------------------------------------------------------------------------------
class Mix:

    #----------------------------------------------------------------------------------------------
    # constructor: self is a 'mix' object created in B3
    def __init__(self, indx, reactor):

        # INITIALIZATION
        # number of isotopes specified in input for mix indx
        self.niso = len(reactor.control.input['mix'][indx]['isoid'])
        # create an object for every isotope of the mix
        #self.iso = []
        #for i in range(self.niso):
        #    self.iso.append(Isotope(i, indx, reactor))
        #    print(self.iso[i].isoname)
