import sys

#--------------------------------------------------------------------------------------------------
class InnerGas:

    #----------------------------------------------------------------------------------------------
    # constructor: self is an 'innergas' object created in B1B
    # indxfuelrod is the index of the fuel rod this object belongs to
    def __init__(self, indxfuelrod, reactor):

        # dictionary of the fuel rod to which the inner gas belongs
        dictfuelrod = reactor.control.input['fuelrod'][indxfuelrod]
        # number of axial nodes
        nz = len(dictfuelrod['fuelid'])
        self.hgap0 = dictfuelrod['hgap']
        self.hgap = self.hgap0
        list = [x['fuelrodid'] for x in reactor.control.input['innergas']]
        if dictfuelrod['id'] in list:
            indx = list.index(dictfuelrod['id'])
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
    # calculate gap conductance: self is an 'innergas' object created in B1B
    # indxfuelrod is the fuel rod index
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
                # call material property function
                pro = reactor.data.matpro( {'type':self.type, 't':tgap} )
                # thermal conductivity (W/m-K): 
                k = pro['k']
                self.hgap[i] = k/dgap
            else:
                self.hgap[i] = self.hgap0[i]

        return self.hgap

