from scipy.interpolate import interp1d

import math
import sys

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
        # list of signals for temperatures of isotopes of mix indx
        self.signal_isotemp = reactor.control.input['mix'][indx]['signaltemp']

        sigma0(self, indx, core, reactor)

#----------------------------------------------------------------------------------------------
# calculates a list of sigma-zeros for each isotope of mix indx
def sigma0(self, indx, core, reactor):

    # number of energy groups
    ng = reactor.control.input['ng']
    # temporal list for sigt after temperature interpolation
    sigt_tmp1 = []
    for i in range(self.niso):
        isoindx = [x.isoid for x in core.iso].index(self.isoid[i])
        # isotope temperature
        temp = reactor.control.signal[self.signal_isotemp[i]]
        # grid temperatures for this isotope
        grid_temp = core.iso[isoindx].temp
        ntemp = len(grid_temp)
        # check if temperature withing the range of grid temperatures for this isotope
        if temp < grid_temp[0] or temp > grid_temp[-1]:
            print('****ERROR: temperature ' + str(temp) + ' K specified in input for isotope ' + self.isoid[i] + ' is out of range of the grid temperatures available in nuclear data library: ' + ''.join([str(int(s)) + ', ' for s in grid_temp])[:-2] + '.')
            sys.exit()
        # grid sigma0s for this isotope
        grid_sig0 = core.iso[isoindx].sig0
        nsig0 = len(grid_sig0)
        sigt_tmp1.append([[0]*nsig0 for j in range(ng)])

        for ig in range(ng):
            for isig0 in range(nsig0):
                # interpolate total xs for isotope temperature temp
                x = grid_temp
                y = [core.iso[isoindx].xs['tot'][ig][itemp][isig0] for itemp in range(ntemp)]
                f = interp1d(x, y) #scipy function
                sigt_tmp1[i][ig][isig0] = f(temp)

    self.sig0 = [[1e10]*self.niso for j in range(ng)]
    # if mix consists of only one isotope then keep sig0 = 1e10
    if self.niso > 1:
        for ig in range(ng):
            # error of sig0 calculation
            err = 1
            # iteration loop
            iter = 0
            while err > 1e-4:
                sigt_tmp2 = [0]*self.niso
                for i in range(self.niso):
                    # interpolate total cross section for sig0
                    x = [math.log10(s) for s in grid_sig0]
                    y = [sigt_tmp1[i][ig][isig0] for isig0 in range(nsig0)]
                    f = interp1d(x, y) #scipy function
                    sigt_tmp2[i] = f(math.log10(self.sig0[ig][i]))
                err = 0
                for i in range(self.niso):
                    # find new sig0
                    sig0new = 0
                    for ii in range(self.niso):
                        if ii != i:
                            sig0new += self.numdens[ii]*sigt_tmp2[ii]
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
