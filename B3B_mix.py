#--------------------------------------------------------------------------------------------------
class Mix:

    #----------------------------------------------------------------------------------------------
    # constructor: self is a 'mix' object created in B3
    def __init__(self, indx, core, reactor):

        # INITIALIZATION
        # mix id
        self.mixid = reactor.control.input['mix'][indx]['mixid']
        # id of isotopes specified in input for mix indx
        self.isoid = reactor.control.input['mix'][indx]['isoid']
        # number densities of isotopes specified in input for mix indx
        self.numdens = reactor.control.input['mix'][indx]['numdens']
        # number of isotopes specified in input for mix indx
        self.niso = len(self.isoid)
        # signals for temperatures of isotopes
        self.signal_isotemp = reactor.control.input['mix'][indx]['signaltemp']

        sigma0(self, indx, core, reactor)

#----------------------------------------------------------------------------------------------
# sigma-zero iterations
def sigma0(self, indx, core, reactor):

    for i in range(self.niso):
        isoindx = [x.isoid for x in core.iso].index(self.isoid[i])
        # isotope temperature
        temp = reactor.control.signal[self.signal_isotemp[i]]
        #print(core.iso[isoindx].xs['tot'][0][0])

