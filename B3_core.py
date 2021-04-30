from B3A_isotope import Isotope
from B3B_mix import Mix

import os

#--------------------------------------------------------------------------------------------------
class Core:

    #----------------------------------------------------------------------------------------------
    # constructor: self is a 'core' object created in B
    def __init__(self, reactor):

        # INITIALIZATION
        if 'pointkinetics' in reactor.solve:
            self.power = 1
            self.ndnp = len(reactor.control.input['betaeff'])
            self.tlife = reactor.control.input['tlife']
            self.dnplmb = reactor.control.input['dnplmb']
            self.betaeff = reactor.control.input['betaeff']
            self.cdnp = [0] * self.ndnp
            for i in range(self.ndnp) :
                self.cdnp[i] = self.betaeff[i]*self.power/(self.dnplmb[i]*self.tlife)

        if 'spatialkinetics' in reactor.solve:
            # create a list of cll isotopes
            self.isoname = [x['isoid'][i] for x in reactor.control.input['mix'] for i in range(len(x['isoid']))]
            #remove duplicates
            self.isoname = list(dict.fromkeys(self.isoname))
            # create an object for every isotope
            self.niso = len(self.isoname)
            self.iso = []
            for i in range(self.niso):
                self.iso.append(Isotope(self.isoname[i], reactor))
            # create an object for every mix
            self.nmix = len(reactor.control.input['mix'])
            self.mix = []
            for i in range(self.nmix):
                self.mix.append(Mix(i, self, reactor))

    #----------------------------------------------------------------------------------------------
    # create right-hand side list: self is a 'core' object created in B
    def calculate_rhs(self, reactor, t):

        # construct right-hand side list
        rhs = []
        if 'pointkinetics' in reactor.solve:
            # read input parameters
            rho = reactor.control.signal['RHO_INS']
            dpowerdt = self.power * (rho - sum(self.betaeff)) / self.tlife
            dcdnpdt = [0] * self.ndnp
            for i in range(self.ndnp) :
                dpowerdt += self.dnplmb[i]*self.cdnp[i]
                dcdnpdt[i] = self.betaeff[i]*self.power/self.tlife - self.dnplmb[i]*self.cdnp[i]
            rhs = [dpowerdt] + dcdnpdt

        if 'spatialkinetics' in reactor.solve:
            for i in range(self.nmix):
                if self.mix[i].update_xs:
                    self.mix[i].sigma0(i, self, reactor)
                    self.mix[i].update_xs = False

            rhs += []

        return rhs
