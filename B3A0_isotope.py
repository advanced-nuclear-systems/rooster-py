import os

#--------------------------------------------------------------------------------------------------
class Isotope:

    #----------------------------------------------------------------------------------------------
    # constructor: self is an 'isotope' object created in B3A
    # indx is the index of the isotope in the mix with index indxmix
    def __init__(self, indx, indxmix, reactor):
        # nuclear data directory
        nddir = reactor.control.input['nddir']
        self.isoname = reactor.control.input['mix'][indxmix]['isoid'][indx]
        f = open(nddir + os.sep + self.isoname, 'r')
        s = f.readline()
        f.close()
        print(nddir + os.sep + self.isoname)

