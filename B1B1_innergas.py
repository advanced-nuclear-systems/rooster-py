#--------------------------------------------------------------------------------------------------
class InnerGas:

    # calculated gap conductance
    hgap = [];
    # user-specified gap conductance
    hgap0 = [];

    #----------------------------------------------------------------------------------------------
    # constructor: self is a 'innergas' object created in B1B
    # indxfuelrod is the index of the fuel rod this object belongs to
    def __init__(self, indxfuelrod, reactor):

        # dictionary of the fuel rod to which the inner gas belongs
        dictfuelrod = reactor.control.input['fuelrod'][indxfuelrod]
        self.hgap0 = dictfuelrod['hgap']

    #----------------------------------------------------------------------------------------------
    # calculate gap conductance
    def calculate_hgap(self, indxfuelrod, reactor, t):

        nz = reactor.solid.fuelrod[indxfuelrod].nz
        self.hgap = [0]*nz
        for i in range(nz):
            if self.hgap0 == 0:
                # fuel surface temperature
                tfuel = reactor.solid.fuelrod[indxfuelrod].fuel[i].temp[-1]
                # clad surface temperature
                tclad = reactor.solid.fuelrod[indxfuelrod].clad[i].temp[0]
                # average temperature
                tgap = 0.5*(tfuel + tclad)
                # gap width
                dgap = reactor.solid.fuelrod[indxfuelrod].clad[i].r[0] - reactor.solid.fuelrod[indxfuelrod].fuel[i].r[-1]
                self.hgap[i] = 20000
            else:
                self.hgap[i] = self.hgap0[i]

        return self.hgap

