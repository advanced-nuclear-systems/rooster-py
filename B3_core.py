from B3A_isotope import Isotope
from B3B_mix import Mix

import math
import os
import sys

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
            self.nx = len(reactor.control.input['coremap'][0])
            for i in range(self.nx):
                if len(reactor.control.input['coremap'][i]) != self.nx:
                    print('****ERROR: all coremap cards should have the same number of nodes.')
                    sys.exit()

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

            # initialize flux
            self.flux = []
            for iz in range(self.nz):
                self.flux.append([])
                for iy in range(self.ny):
                    self.flux[iz].append([])
                    for ix in range(self.nx):
                        self.flux[iz][iy].append([1]*self.ng)

            # initialize map
            self.map = {'dz':[], 'imix':[], 'ipipe':[]}
            mixid_list = [self.mix[i].mixid for i in range(self.nmix)]
            self.nstack = len(reactor.control.input['stack'])
            stackid_list = [reactor.control.input['stack'][i]['stackid'] for i in range(self.nstack)]
            self.npipe = len(reactor.control.input['pipe'])
            pipeid_list = [reactor.control.input['pipe'][i]['id'] for i in range(self.npipe)]
            bc = ['vac','ref']
            for iz in range(self.nz):
                self.map['imix'].append([])
                self.map['ipipe'].append([])
                for iy in range(self.ny):
                    self.map['imix'][iz].append([])
                    self.map['ipipe'][iz].append([])
                    if iz == 0:
                        # bottom boundary conditions
                        botBC = int(reactor.control.input['coregeom']['botBC'])
                        for ix in range(self.nx):
                            self.map['imix'][iz][iy].append(bc[botBC])
                    elif iz == self.nz-1:
                        # top boundary conditions
                        topBC = int(reactor.control.input['coregeom']['topBC'])
                        for ix in range(self.nx):
                            self.map['imix'][iz][iy].append(bc[topBC])
                    else:
                        for ix in range(self.nx):
                            id = reactor.control.input['coremap'][iy][ix]
                            if isinstance(id, float):
                                self.map['imix'][iz][iy].append(bc[int(id)])
                            else:
                                if id not in stackid_list:
                                    print('****ERROR: stack id (' + id + ') in coremap card not specified in stack card.')
                                    sys.exit()
                                else:
                                    # index of stack
                                    istack = stackid_list.index(id)
                                    # id of mix at (ix, iy, iz)
                                    mixid = reactor.control.input['stack'][istack]['mixid'][iz-1]
                                    if mixid not in mixid_list:
                                        print('****ERROR: mix id in stack card (' + mixid + ') not specified in mix card.')
                                        sys.exit()
                                    else:
                                        # index of stack
                                        imix = mixid_list.index(mixid)
                                        self.map['imix'][iz][iy].append(imix)
                                    # id of pipe at (ix, iy, iz)
                                    pipeid = reactor.control.input['stack'][istack]['pipeid'][iz-1]
                                    if pipeid not in pipeid_list:
                                        print('****ERROR: pipe id (' + pipeid + ') in stack card not specified in pipe card.')
                                        sys.exit()
                                    else:
                                        # index of pipe
                                        ipipe = pipeid_list.index(pipeid)
                                        # id of pipenode at (ix, iy, iz)
                                        pipenode = reactor.control.input['stack'][istack]['pipenode'][iz-1]
                                        if pipenode > reactor.control.input['pipe'][ipipe]['nnodes']:
                                            print('****ERROR: pipenode index (' + pipenode + ') in stack card is bigger than number of nodes in pipe ' + pipeid + '.')
                                            sys.exit()
                                        else:
                                            self.map['ipipe'][iz][iy].append((ipipe,pipenode))
                                # node height
                                if len(self.map['dz']) < iz:
                                    self.map['dz'].append(reactor.control.input['pipe'][ipipe]['len']/reactor.control.input['pipe'][ipipe]['nnodes'])
        self.solve_eigenvalue_problem(reactor)

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

    #----------------------------------------------------------------------------------------------
    # solve steady-state eigenvalue problem
    def solve_eigenvalue_problem(self, reactor):

        # initialize fission source
        self.qf = [[[1 for ix in range(self.nx)] for iy in range(self.ny)] for iz in range(self.nz)]
        # eigenvalue self.k equal to ratio of total fission source at two iterations. 
        # flux is normalise to total fission cource = 1 at previous iteration 
        self.k = [1]

        # correct!
        rtol = 1e-5
        atol = 1e-5

        converge_qf = False
        converge_k = False
        while not converge_qf and not converge_k:
            converge_flux = False
            iter = 0
            while not converge_flux:
                converge_flux = True
                iter += 1
                for iz in range(self.nz):
                    for iy in range(self.ny):
                        for ix in range(self.nx):
                            imix = self.map['imix'][iz][iy][ix]
                            # if (ix, iy, iz) is not a boundary condition node (i.e. 'vac' or 'ref')
                            if isinstance(imix, int):
                                xs = self.mix[imix]
                                for ig in range(self.ng):
                                    mlt = 0
                                    dif = 0
                                    # diffusion term: from bottom
                                    imix_n =  self.map['imix'][iz-1][iy][ix]
                                    a_over_v = 0.01/self.map['dz'][iz-1]
                                    if imix_n == 'ref':
                                        pass
                                    elif imix_n == 'vac':
                                        dz = 50*self.map['dz'][iz-1] + 0.71/xs.sigt[imix]
                                        D = 1/(3*xs.sigt[imix])
                                        mlt += D/dz * a_over_v
                                    else:
                                        dz = 50*(self.map['dz'][iz-2] + self.map['dz'][iz-1])
                                        D = dz/(3*xs.sigt[imix_n]*self.map['dz'][iz-2] + 3*xs.sigt[imix]*self.map['dz'][iz-1])
                                        mlt += D/dz * a_over_v
                                        dif -= D*self.flux[iz-1][iy][ix][ig]/dz * a_over_v
                                    
                                    # diffusion term: to top
                                    imix_n =  self.map['imix'][iz+1][iy][ix]
                                    a_over_v = 0.01/self.map['dz'][iz-1]
                                    if imix_n == 'ref':
                                        pass
                                    elif imix_n == 'vac':
                                        dz = 50*self.map['dz'][iz-1] + 0.71/xs.sigt[imix]
                                        D = 1/(3*xs.sigt[imix])
                                        mlt += D/dz * a_over_v
                                    else:
                                        dz = 50*(self.map['dz'][iz-1] + self.map['dz'][iz])
                                        D = dz/(3*xs.sigt[imix]*self.map['dz'][iz-1] + 3*xs.sigt[imix_n]*self.map['dz'][iz])
                                        mlt += D/dz * a_over_v
                                        dif -= D*self.flux[iz+1][iy][ix][ig]/dz * a_over_v
                                    
                                    # diffusion term: from north
                                    imix_n =  self.map['imix'][iz][iy-1][ix]
                                    a_over_v = 0.01/reactor.control.input['coregeom']['pitch']
                                    if imix_n == 'ref':
                                        pass
                                    elif imix_n == 'vac':
                                        dy = 50*reactor.control.input['coregeom']['pitch'] + 0.71/xs.sigt[imix]
                                        D = 1/(3*xs.sigt[imix])
                                        mlt += D/dy * a_over_v
                                    else:
                                        dy = 100*reactor.control.input['coregeom']['pitch']
                                        D = 2/(3*xs.sigt[imix] + 3*xs.sigt[imix_n])
                                        mlt += D/dy * a_over_v
                                        dif -= D*self.flux[iz][iy-1][ix][ig]/dy * a_over_v
                                    
                                    # diffusion term: from south
                                    imix_n =  self.map['imix'][iz][iy+1][ix]
                                    a_over_v = 0.01/reactor.control.input['coregeom']['pitch']
                                    if imix_n == 'ref':
                                        pass
                                    elif imix_n == 'vac':
                                        dy = 50*reactor.control.input['coregeom']['pitch'] + 0.71/xs.sigt[imix]                                        
                                        D = 1/(3*xs.sigt[imix])
                                        mlt += D/dy * a_over_v
                                    else:
                                        dy = 100*reactor.control.input['coregeom']['pitch']
                                        D = 2/(3*xs.sigt[imix] + 3*xs.sigt[imix_n])
                                        mlt += D/dy * a_over_v
                                        dif -= D*self.flux[iz][iy+1][ix][ig]/dy * a_over_v
                                    
                                    # diffusion term: from west
                                    imix_n =  self.map['imix'][iz][iy][ix-1]
                                    a_over_v = 0.01/reactor.control.input['coregeom']['pitch']
                                    if imix_n == 'ref':
                                        pass
                                    elif imix_n == 'vac':
                                        dx = 50*reactor.control.input['coregeom']['pitch'] + 0.71/xs.sigt[imix]
                                        D = 1/(3*xs.sigt[imix])
                                        mlt += D/dx * a_over_v
                                    else:
                                        dx = 100*reactor.control.input['coregeom']['pitch']
                                        D = 2/(3*xs.sigt[imix] + 3*xs.sigt[imix_n])
                                        mlt += D/dx * a_over_v
                                        dif -= D*self.flux[iz][iy][ix-1][ig]/dx * a_over_v
                                    
                                    # diffusion term: to east
                                    imix_n =  self.map['imix'][iz][iy][ix+1]
                                    a_over_v = 0.01/reactor.control.input['coregeom']['pitch']
                                    if imix_n == 'ref':
                                        pass
                                    elif imix_n == 'vac':
                                        dx = 50*reactor.control.input['coregeom']['pitch'] + 0.71/xs.sigt[imix]
                                        D = 1/(3*xs.sigt[imix])
                                        mlt += D/dx * a_over_v
                                    else:
                                        dx = 100*reactor.control.input['coregeom']['pitch']
                                        D = 2/(3*xs.sigt[imix] + 3*xs.sigt[imix_n])
                                        mlt += D/dx * a_over_v
                                        dif -= D*self.flux[iz][iy][ix+1][ig]/dx * a_over_v

                                    # scattering source
                                    qs = 0
                                    # removal xs
                                    sigr = xs.sigt[ig]
                                    for indx in range(len(xs.sigs)):
                                        f = xs.sigs[indx][0][0]
                                        t = xs.sigs[indx][0][1]
                                        if f != ig and t == ig:
                                            qs += xs.sigs[indx][1] * self.flux[iz][iy][ix][f]
                                        if f == ig and t == ig:
                                            sigr -= xs.sigs[indx][1]
                                    mlt += sigr

                                    # fission source
                                    qf = xs.chi[ig]*self.qf[iz][iy][ix]/self.k[-1]

                                    # neutron flux
                                    flux = (-dif + qs + qf)/mlt
                                    if converge_flux : converge_flux = abs(flux - self.flux[iz][iy][ix][ig]) < rtol*abs(flux) + atol or iter >= 10
                                    self.flux[iz][iy][ix][ig] = flux

            converge_qf = True
            for iz in range(self.nz):
                for iy in range(self.ny):
                    for ix in range(self.nx):
                        imix = self.map['imix'][iz][iy][ix]
                        # if (ix, iy, iz) is not a boundary condition node (i.e. 'vac' or 'ref')
                        if isinstance(imix, int):
                            xs = self.mix[imix]
                            qf = 0
                            for ig in range(self.ng):
                                qf += xs.sigp[ig]*self.flux[iz][iy][ix][ig]
                            if converge_qf : converge_qf = abs(qf - self.qf[iz][iy][ix]) < rtol*abs(qf) + atol
                            self.qf[iz][iy][ix] = qf

            converge_k = True
            k = 0
            for iz in range(self.nz):
                for iy in range(self.ny):
                    for ix in range(self.nx):
                        imix = self.map['imix'][iz][iy][ix]
                        # if (ix, iy, iz) is not a boundary condition node (i.e. 'vac' or 'ref')
                        if isinstance(imix, int):
                            for ig in range(self.ng):
                                k += self.qf[iz][iy][ix]
            converge_k = abs(k - self.k[-1]) < rtol*abs(k) + atol
            self.k.append(k)
            print('k-effective: ', '{0:12.5f} '.format(self.k[-1]), '| flux iterations: ', iter)
