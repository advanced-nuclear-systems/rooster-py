from scipy.interpolate import interp1d

import math
import sys

#--------------------------------------------------------------------------------------------------
class Mix:

    #----------------------------------------------------------------------------------------------
    # constructor: self is a 'mix' object created in B3
    def __init__(self, indx, core, reactor):

        # INITIALIZATION
        # number of energy groups
        self.ng = reactor.control.input['ng']
        # mix id
        self.mixid = reactor.control.input['mix'][indx]['mixid']
        # id of isotopes specified in input for mix indx (list)
        self.isoid = reactor.control.input['mix'][indx]['isoid']
        # number of isotopes specified in input for mix indx
        self.niso = len(self.isoid)
        # number densities of isotopes specified in input for mix indx (list)
        self.numdens = reactor.control.input['mix'][indx]['numdens']
        # list of signals for temperatures of isotopes of mix indx
        self.signal_isotemp = reactor.control.input['mix'][indx]['signaltemp']
        # total macroscopic cross sections
        self.sigt = [0]*self.ng
        # absorption macroscopic cross sections
        self.siga = [0]*self.ng
        # flag to calculate xs for mix
        self.update_xs = True

    #----------------------------------------------------------------------------------------------
    # calculates a list of sigma-zeros for each isotope of the mix
    def calculate_sig0(self, core, reactor):

        # perform temperature interpolation for all isotopes and all groups
        # return matrix of microscopic XSs without temperature dimension
        sig_tmp1 = self.interpolate_temp(core, reactor, 'tot')

        self.sig0 = [[1]*self.niso for j in range(self.ng)]
        # if mix consists of only one isotope then keep sig0 = 1
        if self.niso > 1:
            for ig in range(self.ng):
                # error of sig0 calculation
                err = 1
                # iteration loop
                iter = 0
                while err > 1e-4:
                    # given microscopic XSs without temperature dimension perform sig0 interpolation for energy group ig 
                    # for all isotopes of the mix and return matrix of microscopic XSs without sig0 dimension
                    sig_tmp2 = self.interpolate_sig0(ig, core, sig_tmp1)
                    err = 0
                    for i in range(self.niso):
                        # find new sig0
                        sig0new = 0
                        for ii in range(self.niso):
                            if ii != i:
                                sig0new += self.numdens[ii]*sig_tmp2[ii]
                        sig0new = sig0new/self.numdens[i]
                        sig0new = min(sig0new, 1e10)
                        sig0new = max(sig0new, 1)
                        # error control
                        err = max(err, abs(self.sig0[ig][i] - sig0new))
                        self.sig0[ig][i] = sig0new
                    iter += 1
                    if iter > 100:
                        print('****ERROR: too many sig0 iterations for mix ' + self.mixid + ' and energy group ' + ig + '.')
                        sys.exit()
    #----------------------------------------------------------------------------------------------
    # perform temperature interpolation for all isotopes of the mix and 
    # return matrix of microscopic XSs without temperature dimension
    def interpolate_temp(self, core, reactor, reaction_type):

        sig = []
        for i in range(self.niso):
            # index of the isotope i in the global list of isotopes core.iso
            isoindx = [x.isoid for x in core.iso].index(self.isoid[i])
            # isotope temperature
            temp = reactor.control.signal[self.signal_isotemp[i]]
            # grid temperatures for this isotope
            grid_temp = core.iso[isoindx].temp
            ntemp = len(grid_temp)
            # grid sig0s for this isotope
            grid_sig0 = core.iso[isoindx].sig0
            nsig0 = len(grid_sig0)
            # check if temperature withing the range of grid temperatures for this isotope
            if temp < grid_temp[0] or temp > grid_temp[-1]:
                print('****ERROR: temperature ' + str(temp) + ' K specified in input for isotope ' + self.isoid[i] + ' is out of range of the grid temperatures available in nuclear data library: ' + ''.join([str(int(s)) + ', ' for s in grid_temp])[:-2] + '.')
                sys.exit()

            if reaction_type == 'sca':
                # number of entries in elastic scaterring matrix
                n = len(core.iso[isoindx].xs['ela'])
                sig.append([[0]*(nsig0+1) for j in range(n)])
                for j in range(n):
                    # from-to tuple
                    sig[i][j][0] = core.iso[isoindx].xs['ela'][j][0]
                    for isig0 in range(nsig0):
                        # interpolate elastic scattering xs for isotope temperature temp
                        x = grid_temp
                        y = [core.iso[isoindx].xs['ela'][j][1+itemp][isig0] for itemp in range(ntemp)]
                        f = interp1d(x, y) #scipy function
                        sig[i][j][isig0+1] = f(temp)
            else:
                sig.append([[0]*nsig0 for j in range(self.ng)])
                for ig in range(self.ng):
                    for isig0 in range(nsig0):
                        # interpolate total xs for isotope temperature temp
                        x = grid_temp
                        y = [core.iso[isoindx].xs[reaction_type][ig][itemp][isig0] for itemp in range(ntemp)]
                        f = interp1d(x, y) #scipy function
                        sig[i][ig][isig0] = f(temp)
        return sig

    #----------------------------------------------------------------------------------------------
    # given microscopic XSs without temperature dimension sig1 perform sig0 interpolation for energy group j or(f_t) tuple j
    # for all isotopes of the mix and return sig2: microscopic XSs without sig0 dimension
    def interpolate_sig0(self, j, core, sig1, reaction_type):
        sig2 = [0]*self.niso
        for i in range(self.niso):
            # index of the isotope i in the global list of isotopes core.iso
            isoindx = [x.isoid for x in core.iso].index(self.isoid[i])
            # grid sig0s for this isotope
            grid_sig0 = core.iso[isoindx].sig0
            nsig0 = len(grid_sig0)

            if reaction_type == 'sca':
                # interpolate sig1 cross section for sig0
                x = grid_sig0
                y = [sig1[i][j][isig0+1] for isig0 in range(nsig0)]
                f = interp1d(x, y) #scipy function
                ig = sig1[i][j][0][0]
                sig2[i] = f(self.sig0[ig][i])
            else:
                ig = j
                # interpolate sig1 cross section for sig0
                x = grid_sig0
                y = [sig1[i][ig][isig0] for isig0 in range(nsig0)]
                f = interp1d(x, y) #scipy function
                sig2[i] = f(self.sig0[ig][i])
        return sig2

    #----------------------------------------------------------------------------------------------
    # calculates absorption macroscopic cross sections for the mix
    def calculate_siga(self, core, reactor):
        # perform temperature and sig0 interpolations for all isotopes and all groups
        sig_tmp1 = self.interpolate_temp(core, reactor, 'abs')
        sig_tmp2 = [self.interpolate_sig0(ig, core, sig_tmp1, 'abs') for ig in range(self.ng)]
        for i in range(self.ng):
            self.siga[i] = 0
            for j in range(self.niso):
                self.siga[i] += self.numdens[j]*sig_tmp2[i][j]

    #----------------------------------------------------------------------------------------------
    # calculates total macroscopic absorption cross sections for the mix
    def calculate_sigt(self, core, reactor):
        # perform temperature and sig0 interpolations for all isotopes and all groups
        sig_tmp1 = self.interpolate_temp(core, reactor, 'tot')
        sig_tmp2 = [self.interpolate_sig0(ig, core, sig_tmp1, 'tot') for ig in range(self.ng)]
        for i in range(self.ng):
            self.sigt[i] = 0
            for j in range(self.niso):
                self.sigt[i] += self.numdens[j]*sig_tmp2[i][j]

    #----------------------------------------------------------------------------------------------
    # calculates macroscopic scattering cross sections for the mix
    def calculate_sigs(self, core, reactor):
        # perform temperature and sig0 interpolations for all isotopes and all groups
        sig_tmp1 = self.interpolate_temp(core, reactor, 'sca')
        sig_tmp2 = [self.interpolate_sig0(j, core, sig_tmp1, 'sca') for j in range(len(sig_tmp1[0]))]
        #for i in range(self.ng):
        #    self.sigt[i] = 0
        #    for j in range(self.niso):
        #        self.sigt[i] += self.numdens[j]*sig_tmp2[i][j]
