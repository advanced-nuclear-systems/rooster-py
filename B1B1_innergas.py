import sys

#--------------------------------------------------------------------------------------------------
class InnerGas:

    # calculated gap conductance
    hgap = []
    # user-specified gap conductance
    hgap0 = []
    # gas pressure
    p = 0
    # gas temperature
    temp = []
    # gas type
    type = ''

    #----------------------------------------------------------------------------------------------
    # constructor: self is a 'innergas' object created in B1B
    # indxfuelrod is the index of the fuel rod this object belongs to
    def __init__(self, indxfuelrod, reactor):

        # dictionary of the fuel rod to which the inner gas belongs
        dictfuelrod = reactor.control.input['fuelrod'][indxfuelrod]
        # number of axial nodes
        nz = len(dictfuelrod['fuelid'])
        self.hgap0 = dictfuelrod['hgap']
        indx = [x['fuelrodid'] for x in reactor.control.input['innergas']].index(dictfuelrod['id'])
        matid = [x['matid'] for x in reactor.control.input['innergas']][indx]
        # find the gas material id in the list of materials
        try:
            igas = [x['id'] for x in reactor.control.input['mat']].index(matid)
        except:
            print('****ERROR: gas material id ' + matid + ' is not specified in the \'mat\' card of input.')
            sys.exit()
        # dictionary of material properties of the current inner gas
        mat = reactor.control.input['mat'][igas]
        # material type of gas
        self.type = mat['type']
        # initial pressure of gas
        self.p = mat['p0']
        # initial temperature of gas
        self.temp = [mat['temp0']]*nz

    #----------------------------------------------------------------------------------------------
    # calculate gap conductance
    def calculate_hgap(self, indxfuelrod, reactor, t):

        nz = reactor.solid.fuelrod[indxfuelrod].nz
        self.hgap = [0]*nz
        for i in range(nz):
            if self.hgap0[i] == 0:
                # fuel surface temperature
                tfuel = reactor.solid.fuelrod[indxfuelrod].fuel[i].temp[-1]
                # clad surface temperature
                tclad = reactor.solid.fuelrod[indxfuelrod].clad[i].temp[0]
                # average temperature
                tgap = 0.5*(tfuel + tclad)
                # gap width
                dgap = reactor.solid.fuelrod[indxfuelrod].clad[i].r[0] - reactor.solid.fuelrod[indxfuelrod].fuel[i].r[-1]
                # GAS PROPERTIES:
                if self.type == 'he':
                    # thermal conductivity (W/m-K): 
                    k = 2.639e-3*tgap**0.7085
                self.hgap[i] = k/dgap
            else:
                self.hgap[i] = self.hgap0[i]

        return self.hgap

