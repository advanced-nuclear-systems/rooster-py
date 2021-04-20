#--------------------------------------------------------------------------------------------------
class InnerGas:

    # constructor: self is a 'innergas' object created in B1B
    def __init__(self, indxfuelrod, reactor):

        ## dictionary of the fuel rod to which the inner gas belongs
        #dictfuelrod = reactor.control.input['fuelrod'][indxfuelrod]
        ## current inner gas id
        #innergasid = dictfuelrod['innergasid'][indx]
        ## list of inner gas dictionaries specified in input
        #list = reactor.control.input['innergas']
        ## index of the current inner gas in the list of inner gas dictionaries
        #i = [x['id'] for x in list].index(innergasid)

        #for i in range(nz):
        #    # inner gas material id
        #    self.matid = list[i]['matid']
        pass

    #----------------------------------------------------------------------------------------------
    # calculate gap conductance
    def calculate_hgap(self, indxfuelrod, reactor, t):

        nz = reactor.solid.fuelrod[indxfuelrod].nz
        hgap = [0]*nz
        for i in range(nz):
            # fuel surface temperature
            tfuel = reactor.solid.fuelrod[indxfuelrod].fuel[i].temp[-1]
            # clad surface temperature
            tclad = reactor.solid.fuelrod[indxfuelrod].clad[i].temp[0]
            # average temperature
            tgap = 0.5*(tfuel + tclad)
            # gap width
            dgap = reactor.solid.fuelrod[indxfuelrod].clad[i].r[0] - reactor.solid.fuelrod[indxfuelrod].fuel[i].r[-1]
            hgap[i] = 10000

        return hgap

