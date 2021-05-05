from B3A_isotope import Isotope
from B3B_mix import Mix

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
        #print(self.map['dz'])

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

        # accuracy of the solution
        eps = 1
        while eps > 1e-6:
            for iz in range(self.nz):
                for iy in range(self.ny):
                    for ix in range(self.nx):
                        mix = self.map['mix'][iz][iy][ix]
                        mixid_list = [self.mix[i].mixid for i in range(self.nmix)]
                        if mix in mixid_list:
                            imix = mixid_list(mix)
                        # absorption rate
                        #self.flux[iz][iy][ix][ig]*
