import os

#--------------------------------------------------------------------------------------------------
class Isotope:

    #----------------------------------------------------------------------------------------------
    # constructor: self is an 'isotope' object created in B3A
    # indx is the index of the isotope in the mix with index indxmix
    def __init__(self, indx, indxmix, reactor):
        # nuclear data directory
        nddir = reactor.control.input['nddir']
        # isotope name
        self.isoname = reactor.control.input['mix'][indxmix]['isoid'][indx]
        # open, read, split by eol and close the isotope data file
        f = open(nddir + os.sep + self.isoname, 'r')
        s = f.read().replace('\r\n', '\n').split('\n')
        f.close()
        # base temperatures
        self.temp = [float(s[i]) for i in range(1,int(s[0])+1)]
        del s[0:int(s[0])+1]
        # base sigma-zeros
        self.sig0 = [float(s[i]) for i in range(1,int(s[0])+1)]
        del s[0:int(s[0])+1]

        #print(nddir + os.sep + self.isoname)
        #print(self.temp)
        #print(self.sig0)

