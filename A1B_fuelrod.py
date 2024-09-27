import sys

from A1B0_fuel import Fuel
from A1B1_innergas import InnerGas
from A1B2_clad import Clad

#--------------------------------------------------------------------------------------------------
class FuelRod:

    #----------------------------------------------------------------------------------------------
    # constructor: self is a 'fuelrod' object created in B1 and indx is the index of this object in the list of fuelrods
    def __init__(self, indx, reactor):

        # INITIALIZATION
        # dictionary of the fuel rod
        dictfuelrod = reactor.control.input['fuelrod'][indx]
        
        # fuel rod id
        self.id = dictfuelrod['id']
        # number of axial layers specified in input for fuel rod indx
        self.nz = len(dictfuelrod['fuelid'])
        # axial mesh size
        self.dz = []
        for i in range(self.nz):
            # check existence of neighbouring fluid pipe
            jpipe = (dictfuelrod['pipeid'][i], dictfuelrod['pipenode'][i])
            if not jpipe[0] in reactor.fluid.pipeid:
                print('****ERROR: pipe id ' + jpipe[0] + ' given in \'fuelrod\' card is not specified in the \'pipe\' card of input.')
                sys.exit()
            else:
                # pipe index
                ipipe = reactor.fluid.pipeid.index(jpipe[0])
            # check existence of neighbouring fluid pipe node
            if jpipe[1] > reactor.fluid.pipennodes[ipipe]:
                print('****ERROR: pipe node index (' + str(jpipe[1]) + ') given in \'fuelrod\' card exceeds number of nodes (' + str(reactor.fluid.pipennodes[ipipe]) + ') of pipe ' + jpipe[0])
                sys.exit()
            # pipe node indexes
            jpipe = (ipipe, jpipe[1]-1)
            self.dz.append(reactor.fluid.len[jpipe[0]]/reactor.fluid.pipennodes[jpipe[0]])

        # create an object for every axial layer of fuel
        self.fuel = []
        for i in range(self.nz):
            self.fuel.append(Fuel(i, indx, self.dz[i], reactor))

        # create an object for inner gas
        self.innergas = InnerGas(indx, reactor)

        # create an object for every axial layer of clad
        self.clad = []
        for i in range(self.nz):
            self.clad.append(Clad(i, indx, reactor))

    #----------------------------------------------------------------------------------------------
    # compose right-hand side list: self is a 'fuelrod' object created in B1,
    # indx is the fuel rod index
    def compose_rhs(self, indx, reactor, t):

        # construct right-hand side list
        rhs = []
        for i in range(self.nz):
            rhs += self.fuel[i].calculate_rhs(i, indx, reactor, t)
            rhs += self.clad[i].calculate_rhs(i, indx, reactor, t)
        return rhs
