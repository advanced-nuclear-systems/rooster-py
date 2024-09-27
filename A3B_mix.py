from scipy.interpolate import interp1d

import math
import sys

#--------------------------------------------------------------------------------------------------
class Mix:

    #----------------------------------------------------------------------------------------------
    # constructor: self is a 'mix' object created in B3
    def __init__(self, indx, core, reactor):

        reactor.control.evaluate_signals(reactor, reactor.control.input['t0'])

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

            if reaction_type == 'ela' or reaction_type == 'ela1':
                # number of entries in elastic scattering matrix
                n = len(core.iso[isoindx].xs[reaction_type])
                sig.append([[0]*(nsig0+1) for j in range(n)])
                for j in range(n):
                    # from-to tuple
                    sig[i][j][0] = core.iso[isoindx].xs[reaction_type][j][0]
                    for isig0 in range(nsig0):
                        # interpolate elastic scattering xs for isotope temperature temp
                        x = grid_temp
                        y = [core.iso[isoindx].xs[reaction_type][j][1+itemp][isig0] for itemp in range(ntemp)]
                        f = interp1d(x, y) # scipy function
                        sig[i][j][isig0+1] = f(temp)

            elif reaction_type == 'elan':
                for nlgndr in range(2):
                    sig.append([])
                    # number of entries in elastic scattering matrix
                    n = len(core.iso[isoindx].xs[reaction_type][nlgndr])
                    sig[nlgndr].append([[0]*(nsig0+1) for j in range(n)])
                    for j in range(n):
                        # from-to tuple
                        sig[nlgndr][i][j][0] = core.iso[isoindx].xs[reaction_type][nlgndr][j][0]
                        for isig0 in range(nsig0):
                            # interpolate elastic scattering xs for isotope temperature temp
                            x = grid_temp
                            y = [core.iso[isoindx].xs[reaction_type][nlgndr][j][1+itemp][isig0] for itemp in range(ntemp)]
                            f = interp1d(x, y) # scipy function
                            sig[nlgndr][i][j][isig0+1] = f(temp)

            elif reaction_type == 'nub':
                sig.append([0]*self.ng)
                for ig in range(self.ng):
                    # interpolate nubar for isotope temperature temp
                    x = grid_temp
                    y = [core.iso[isoindx].xs[reaction_type][ig][itemp] for itemp in range(ntemp)]
                    # scipy function
                    f = interp1d(x, y)
                    sig[i][ig] = f(temp)
            else:
                sig.append([[0]*nsig0 for j in range(self.ng)])
                for ig in range(self.ng):
                    for isig0 in range(nsig0):
                        # interpolate total xs for isotope temperature temp
                        x = grid_temp
                        y = [core.iso[isoindx].xs[reaction_type][ig][itemp][isig0] for itemp in range(ntemp)]
                        # scipy function
                        f = interp1d(x, y)
                        sig[i][ig][isig0] = f(temp)
        return sig

    #----------------------------------------------------------------------------------------------
    # given microscopic XSs without temperature dimension sig1[iso][ig][isig0] perform sig0 interpolation for energy group ig
    # for all isotopes of the mix and return sig2[iso][ig]: microscopic XSs without sig0 dimension
    def interpolate_sig0(self, ig, core, sig1):
        sig2 = [0]*self.niso
        for i in range(self.niso):
            # index of the isotope i in the global list of isotopes core.iso
            isoindx = [x.isoid for x in core.iso].index(self.isoid[i])
            # grid sig0s for this isotope
            grid_sig0 = core.iso[isoindx].sig0
            nsig0 = len(grid_sig0)

            # interpolate sig1 cross section for sig0
            x = grid_sig0
            y = [sig1[i][ig][isig0] for isig0 in range(nsig0)]
            f = interp1d(x, y) # scipy function
            sig2[i] = f(self.sig0[ig][i])
        return sig2

    #----------------------------------------------------------------------------------------------
    # calculates total macroscopic cross sections for the mix
    def calculate_sigt(self, core, reactor):
        # perform temperature and sig0 interpolations for all isotopes and all groups
        sig_tmp1 = self.interpolate_temp(core, reactor, 'tot')
        sig_tmp2 = [self.interpolate_sig0(ig, core, sig_tmp1) for ig in range(self.ng)]
        self.sigt = [0]*self.ng
        for ig in range(self.ng):
            for i in range(self.niso):
                self.sigt[ig] += self.numdens[i]*sig_tmp2[ig][i]

    #----------------------------------------------------------------------------------------------
    # calculates transport macroscopic cross sections for the mix
    def calculate_sigtra(self, core, reactor):
        # perform temperature and sig0 interpolations for all isotopes and all groups
        sig_tmp1 = self.interpolate_temp(core, reactor, 'tra')
        sig_tmp2 = [self.interpolate_sig0(ig, core, sig_tmp1) for ig in range(self.ng)]
        self.sigtra = [0.]*self.ng
        for ig in range(self.ng):
            for i in range(self.niso):
                self.sigtra[ig] += self.numdens[i]*sig_tmp2[ig][i]


    #----------------------------------------------------------------------------------------------
    # calculates production and fission macroscopic cross sections for the mix
    def calculate_sigp(self, core, reactor):
        # perform temperature and sig0 interpolations for fission xs for all isotopes and all groups
        sig_tmp1 = self.interpolate_temp(core, reactor, 'fis')
        sig_tmp2 = [self.interpolate_sig0(ig, core, sig_tmp1) for ig in range(self.ng)]
        # perform temperature interpolations for nubar for all isotopes and all groups
        nubar = self.interpolate_temp(core, reactor, 'nub')
        self.sigp = [0]*self.ng
        self.sigf = [0]*self.ng
        for ig in range(self.ng):
            for i in range(self.niso):
                self.sigp[ig] += self.numdens[i]*nubar[i][ig]*sig_tmp2[ig][i]
                self.sigf[ig] += self.numdens[i]*sig_tmp2[ig][i]

    #----------------------------------------------------------------------------------------------
    # calculates fission spectrum for the mix
    def calculate_chi(self, core):
        self.chi = [0]*self.ng
        for ig in range(self.ng):
            for i in range(self.niso):
                # index of the isotope i in the global list of isotopes core.iso
                isoindx = [x.isoid for x in core.iso].index(self.isoid[i])
                self.chi[ig] += self.numdens[i]*core.iso[isoindx].xs['chi'][ig]
        # normalize fission spectrum
        s = sum(self.chi)
        if s > 0 : self.chi = [self.chi[ig]/s for ig in range(self.ng)]

    #----------------------------------------------------------------------------------------------
    # calculates macroscopic scattering cross sections for the mix
    def calculate_sigsn(self, core, reactor):
        # perform temperature and sig0 interpolations for all isotopes and all groups
        sig_tmp1 = self.interpolate_temp(core, reactor, 'elan')
        self.sigsn = []
        for nlgndr in range(2):
            self.sigsn.append([])
            for i in range(self.niso):
                # index of the isotope i in the global list of isotopes core.iso
                isoindx = [x.isoid for x in core.iso].index(self.isoid[i])
                # grid sig0s for this isotope
                x = core.iso[isoindx].sig0
                nsig0 = len(x)
                # number of entries in elastic scattering matrix for isotope i
                nesca = len(sig_tmp1[nlgndr][i])
                for j in range(nesca):
                    # (from, to) tuple
                    f_t = sig_tmp1[nlgndr][i][j][0]
                    # scattering xs corresponding to x
                    y = [sig_tmp1[nlgndr][i][j][isig0+1] for isig0 in range(nsig0)]
                    # scipy function
                    f = interp1d(x, y)
                    # index of 'from' group
                    ig = f_t[0]
                    # interpolate scattering cross section for sig0 of group ig for isotope i
                    value = f(self.sig0[ig][i])
                    f_t_list = [s[0] for s in self.sigsn[nlgndr]]
                    if f_t in f_t_list:
                        # if the (from, to) tuple is already in the self.sigsn list
                        indx = f_t_list.index(f_t)
                        self.sigsn[nlgndr][indx][1] += self.numdens[i]*value
                    else:
                        self.sigsn[nlgndr].append([f_t, self.numdens[i]*value])

                # number of entries in inelastic scattering matrix for isotope i
                nisca = len(core.iso[isoindx].xs['ine'])
                for j in range(nisca):
                    # (from, to) tuple
                    f_t = core.iso[isoindx].xs['ine'][j][0]
                    # inelastic scattering xs
                    if nlgndr == 0:
                        value = core.iso[isoindx].xs['ine'][j][1]
                    else:
                        value = 0
                    f_t_list = [s[0] for s in self.sigsn[nlgndr]]
                    if f_t in f_t_list:
                        # if the (from, to) tuple is already in the self.sigsn[nlgndr] list
                        indx = f_t_list.index(f_t)
                        self.sigsn[nlgndr][indx][1] += self.numdens[i]*value
                    else:
                        self.sigsn[nlgndr].append([f_t, self.numdens[i]*value])

    #----------------------------------------------------------------------------------------------
    # calculates macroscopic n2n cross sections for the mix
    def calculate_sign2n(self, core, reactor):
        self.sign2n = []
        for i in range(self.niso):
            # index of the isotope i in the global list of isotopes core.iso
            isoindx = [x.isoid for x in core.iso].index(self.isoid[i])
            # number of entries in n2n matrix for isotope i
            nn2n = len(core.iso[isoindx].xs['n2n'])
            for j in range(nn2n):
                # (from, to) tuple
                f_t = core.iso[isoindx].xs['n2n'][j][0]
                # n2n scattering xs
                value = core.iso[isoindx].xs['n2n'][j][1]
                f_t_list = [s[0] for s in self.sign2n]
                if f_t in f_t_list:
                    # if the (from, to) tuple is already in the self.sign2n list
                    indx = f_t_list.index(f_t)
                    self.sign2n[indx][1] += self.numdens[i]*value
                else:
                    self.sign2n.append([f_t, self.numdens[i]*value])


    #----------------------------------------------------------------------------------------------
    # calculates kerma factors for the mix
    def calculate_kerma(self, core, reactor):
        # a special case: no dependence on temperature but dependence on sig0
        sig_tmp1 = []
        for i in range(self.niso):
            # index of the isotope i in the global list of isotopes core.iso
            isoindx = [x.isoid for x in core.iso].index(self.isoid[i])
            nsig0 = len(core.iso[isoindx].sig0)
            sig_tmp1.append([[0]*nsig0 for j in range(self.ng)])
            for ig in range(self.ng):
                for isig0 in range(nsig0):
                    sig_tmp1[i][ig][isig0] = core.iso[isoindx].xs['kerma'][ig][isig0]
        # perform sig0 interpolations for all isotopes and all groups
        sig_tmp2 = [self.interpolate_sig0(ig, core, sig_tmp1) for ig in range(self.ng)]
        self.kerma = [0]*self.ng
        for ig in range(self.ng):
            for i in range(self.niso):
                self.kerma[ig] += self.numdens[i]*sig_tmp2[ig][i]
