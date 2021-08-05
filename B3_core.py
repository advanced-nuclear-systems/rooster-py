from B3A_isotope import Isotope
from B3B_mix import Mix

import B3_coreF
import numpy
import sys
import time

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

            # core geometry flag
            self.geom = reactor.control.input['coregeom']['geom']

            # number of energy groups
            self.ng = reactor.control.input['ng']

            # core mesh
            self.nz = len(reactor.control.input['stack'][0]['mixid'])
            for i in range(len(reactor.control.input['stack'])):
                if len(reactor.control.input['stack'][i]['mixid']) != self.nz:
                    print('****ERROR: all stacks should have the same number of axial nodes:', self.nz)
                    sys.exit()
            # add bottom and top layers for boundary conditions
            self.nz += 2
            self.ny = len(reactor.control.input['coremap'])
            self.nx = len(reactor.control.input['coremap'][0])
            if self.geom == 'hex24':
                self.nt = 24
            elif self.geom == 'hex06':
                self.nt = 6
            else:
                self.nt = 1
            for i in range(self.nx):
                if len(reactor.control.input['coremap'][i]) != self.nx:
                    print('****ERROR: all coremap cards should have the same number of nodes.')
                    sys.exit()
            if 'power0' not in reactor.control.input:
                print('***ERROR: there is no card power0 in the input.')
                sys.exit()

            # initialize flux
            self.flux = numpy.ones(shape=(self.nz, self.ny, self.nx, self.nt, self.ng), order='F')

            # create a list of all isotopes
            self.isoname = [x['isoid'][i] for x in reactor.control.input['mix'] for i in range(len(x['isoid']))]
            #remove duplicates
            self.isoname = list(dict.fromkeys(self.isoname))
            # create an object for every isotope
            self.niso = len(self.isoname)
            self.iso = []
            for i in range(self.niso):
                self.iso.append(Isotope(self.isoname[i], reactor))
                self.iso[i].print_xs = True

            # create an object for every mix
            self.nmix = len(reactor.control.input['mix'])
            self.mix = []
            for i in range(self.nmix):
                self.mix.append(Mix(i, self, reactor))

            # calculate sig0 and macroscopic cross sections
            for i in range(self.nmix):
                self.mix[i].calculate_sig0(self, reactor)
                self.mix[i].calculate_sigt(self, reactor)
                self.mix[i].calculate_sigtra(self, reactor)
                self.mix[i].calculate_sigp(self, reactor)
                self.mix[i].calculate_chi(self)
                self.mix[i].calculate_sigs(self, reactor)
                self.mix[i].calculate_sigs1(self, reactor)
                self.mix[i].calculate_sign2n(self, reactor)
                self.mix[i].calculate_kerma(self, reactor)
                self.mix[i].update_xs = False
                self.mix[i].print_xs = True

                tac = time.time()
                print('{0:.3f}'.format(tac - reactor.tic), ' s | mix cross sections processed: ', self.mix[i].mixid)
                reactor.tic = tac

            # initialize map
            self.map = {'dz':[], 'imix':[], 'ipipe':[]}
            mixid_list = [self.mix[i].mixid for i in range(self.nmix)]
            self.nstack = len(reactor.control.input['stack'])
            stackid_list = [reactor.control.input['stack'][i]['stackid'] for i in range(self.nstack)]
            self.npipe = len(reactor.control.input['pipe'])
            pipeid_list = [reactor.control.input['pipe'][i]['id'] for i in range(self.npipe)]
            # vacuum is -1 and reflective is -2
            bc = [-1,-2]
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
            # core assembly pitch
            self.pitch = 100*reactor.control.input['coregeom']['pitch']

            #f = open('tmp_map.txt', 'w')
            #for iz in range(self.nz):
            #    for iy in range(self.ny):
            #        if iy % 2 == 0:
            #            pass
            #        else:
            #            f.write(' ')
            #        for ix in range(self.nx):
            #            f.write(str(self.map['imix'][iz][iy][ix]+1) + ' ')
            #        f.write('\n')
            #    f.write('\niz = '+str(iz)+'\n')
            #f.close()

            # initialize multiplication factor
            self.k = numpy.array([1.])

            # prepare arrays for Fortran solver of eigenvalue problem             
            # total cross section
            sigt = numpy.array([[self.mix[imix].sigt[ig] for ig in range(self.ng)] for imix in range(self.nmix)], order='F')
            # production cross section
            sigp = numpy.array([[self.mix[imix].sigp[ig] for ig in range(self.ng)] for imix in range(self.nmix)], order='F')
            # number of elements in scattering matrix
            nsigs = numpy.array([len(self.mix[imix].sigs) for imix in range(self.nmix)], order='F')
            # scattering matrix
            sigs = numpy.zeros(shape=(self.nmix, max(nsigs)), order='F')
            # 'from' index of scattering matrix elements
            fsigs = numpy.zeros(shape=(self.nmix, max(nsigs)), dtype=int, order='F')
            # 'to' index of scattering matrix elements
            tsigs = numpy.zeros(shape=(self.nmix, max(nsigs)), dtype=int, order='F')
            # number of elements in n2n matrix
            nsign2n = numpy.array([len(self.mix[imix].sign2n) for imix in range(self.nmix)], order='F')
            # n2n matrix )1D_
            sign2n = numpy.zeros(shape=(self.nmix, max(nsign2n)), order='F')
            # 'from' index of n2n matrix elements
            fsign2n = numpy.zeros(shape=(self.nmix, max(nsign2n)), dtype=int, order='F')
            # 'to' index of n2n matrix elements
            tsign2n = numpy.zeros(shape=(self.nmix, max(nsign2n)), dtype=int, order='F')
            # fission source
            chi = numpy.array([[self.mix[imix].chi[ig] for ig in range(self.ng)] for imix in range(self.nmix)], order='F')
            # fission cross section
            sigf = numpy.array([[self.mix[imix].sigf[ig] for ig in range(self.ng)] for imix in range(self.nmix)], order='F')
            # axial nodalization (cm)
            dz = numpy.array(self.map['dz'], order='F')*100.

            # fill out scattering arrays
            for imix in range(self.nmix):
                for indx in range(nsigs[imix]):
                    fsigs[imix][indx] = self.mix[imix].sigs[indx][0][0]
                    tsigs[imix][indx] = self.mix[imix].sigs[indx][0][1]
                    sigs[imix][indx] = self.mix[imix].sigs[indx][1]

            # fill out n2n arrays
            for imix in range(self.nmix):
                for indx in range(nsign2n[imix]):
                    fsign2n[imix][indx] = self.mix[imix].sign2n[indx][0][0]
                    tsign2n[imix][indx] = self.mix[imix].sign2n[indx][0][1]            
                    sign2n[imix][indx] = self.mix[imix].sign2n[indx][1]

            # transport cross section = total cross section - first Legendre component of elastic out-scattering cross section 
            sigtra = numpy.array([[self.mix[imix].sigtra[ig] for ig in range(self.ng)] for imix in range(self.nmix)], order='F')

            # call the Fortran eigenvalue problem solver
            B3_coreF.solve_eigenvalue_problem('MC', self.geom, self.nz, self.ny, self.nx, self.nt, self.ng, self.nmix, \
                                              self.flux, self.map['imix'], sigt, sigtra, sigp, \
                                              nsigs, fsigs, tsigs, sigs, nsign2n, fsign2n, tsign2n, sign2n, chi, sigf, \
                                              self.pitch, dz)
            # power distribution
            self.pow = numpy.zeros(shape=(self.nz, self.ny, self.nx), order='F')
            self.powxy = numpy.zeros(shape=(self.ny, self.nx), order='F')
            if self.geom == 'square':
                az = self.pitch**2
            elif self.geom == 'hex01':
                az = numpy.sqrt(3.)/2.*self.pitch**2
            elif self.geom == 'hex06':
                az = numpy.sqrt(3.)/2.*self.pitch**2/6
            elif self.geom == 'hex24':
                az = numpy.sqrt(3.)/2.*self.pitch**2/24
            # power normalization factor
            factor = 0.
            for iz in range(self.nz):
                for iy in range(self.ny):
                    for ix in range(self.nx):
                        # if (ix, iy, iz) is not a boundary condition node, i.e. not -1 (vac) and not -2 (ref)
                        imix = self.map['imix'][iz][iy][ix]
                        if imix >= 0 and any(self.mix[imix].sigf) > 0:
                            vol = az*self.map['dz'][iz-1]
                            for it in range(self.nt):
                                for ig in range(self.ng):
                                    self.pow[iz][iy][ix] += self.mix[imix].kerma[ig]*self.flux[iz][iy][ix][it][ig]*vol
                                    #self.pow[iz][iy][ix] += self.mix[imix].sigf[ig]*self.flux[iz][iy][ix][it][ig]*vol * 200. * 1.6022e-19
                                self.powxy[iy][ix] += self.pow[iz][iy][ix]
                                factor += self.pow[iz][iy][ix]
            factor = reactor.control.input['power0'] / factor
            # normalize flux and power to power0
            for iz in range(self.nz):
                for iy in range(self.ny):
                    for ix in range(self.nx):
                        self.pow[iz][iy][ix] *= factor
                        for it in range(self.nt):
                            for ig in range(self.ng):
                                self.flux[iz][iy][ix][it][ig] *= factor
            for iy in range(self.ny):
                for ix in range(self.nx):
                    self.powxy[iy][ix] *= factor

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
                    self.mix[i].calculate_kerma(self, reactor)
                    self.mix[i].update_xs = False
                    self.mix[i].print_xs = True

            rhs += []

        return rhs
