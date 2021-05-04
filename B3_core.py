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
            # number of energy groups
            self.ng = reactor.control.input['ng']
            # core mesh
            self.nz = len(reactor.control.input['stack'][0]['mixid'])
            for i in range(len(reactor.control.input['stack'])):
                if len(reactor.control.input['stack'][i]['mixid']) != self.nz:
                    print('****ERROR: all stacks should have the same number of axial nodes.')
                    sys.exit()
            # add bottom and top layers for boundary conditions
            self.nz += 2
            self.ny = len(reactor.control.input['coremap'])
            self.nx = [len(reactor.control.input['coremap'][i]) for i in range(self.ny)]
            # initialize flux
            self.flux = []
            for iz in range(self.nz):
                self.flux.append([])
                for iy in range(self.ny):
                    self.flux[iz].append([])
                    for ix in range(self.nx[iy]):
                        self.flux[iz][iy].append([1]*self.ng)

            # create a list of all isotopes
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
            # calculate sig0 and macroscopic cross sections
            for i in range(self.nmix):
                self.mix[i].calculate_sig0(self, reactor)
                self.mix[i].calculate_sigt(self, reactor)
                self.mix[i].calculate_siga(self, reactor)
                self.mix[i].calculate_sigp(self, reactor)
                self.mix[i].calculate_chi(self)
                self.mix[i].calculate_sigs(self, reactor)
                self.mix[i].calculate_sign2n(self, reactor)
                self.mix[i].update_xs = False
                self.mix[i].print_xs = True

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
                    self.mix[i].calculate_sig0(self, reactor)
                    self.mix[i].calculate_sigt(self, reactor)
                    self.mix[i].calculate_siga(self, reactor)
                    self.mix[i].calculate_sigp(self, reactor)
                    self.mix[i].calculate_chi(self)
                    self.mix[i].calculate_sigs(self, reactor)
                    self.mix[i].calculate_sign2n(self, reactor)
                    self.mix[i].update_xs = False
                    self.mix[i].print_xs = True

            rhs += []

        return rhs
