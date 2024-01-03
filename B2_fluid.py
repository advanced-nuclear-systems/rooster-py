from scipy import linalg

import math
import sys

#--------------------------------------------------------------------------------------------------
class Fluid:

    # constructor: self is a 'fluid' object created in B
    def __init__(self, reactor):

        if 'fluid' not in reactor.solve:
            return

        # INITIALIZATION
        # list of pipe id's
        self.pipeid = [x['id'] for x in reactor.control.input['pipe']]
        # list of pipe types
        self.pipetype = [x['type'] for x in reactor.control.input['pipe']]
        # number of pipes
        self.npipe = len(self.pipetype)
        # number of freelevel pipes
        self.npipef = self.pipetype.count('freelevel')
        # list of pipe hydraulic diameters
        self.dhyd = [x['dhyd'] for x in reactor.control.input['pipe']]
        # list of pipe length
        self.len = [x['len'] for x in reactor.control.input['pipe']]
        # list of pipe direction
        self.dir = [x['dir'] for x in reactor.control.input['pipe']]
        # list of pipe flow area
        self.areaz = [x['areaz'] for x in reactor.control.input['pipe']]
        # list of numbers of pipe nodes
        self.pipennodes = [x['nnodes'] for x in reactor.control.input['pipe']]
        # list of user-specified pipe temperature signal
        self.signaltemp = [x['signaltemp'] for x in reactor.control.input['pipe']]

        # lists for pressure, temperature, type, velocity, reynolds, prandtl, peclet and thermal boundary condition id
        self.p0 = []
        self.p = []
        self.temp = []
        self.type = []
        self.vel = []
        self.re = []
        self.pr = []
        self.pe = []
        # process coolant names
        for i in range(self.npipe):
            cool = reactor.control.input['pipe'][i]['matid']
            # find the coolant id in the list of coolants
            try:
                icool = [x['id'] for x in reactor.control.input['mat']].index(cool)
            except:
                print('****ERROR: input coolant id ' + cool + ' is not specified in the \'mat\' card.')
                sys.exit()
            type = reactor.control.input['mat'][icool]['type']
            p0 = reactor.control.input['mat'][icool]['p0']
            temp0 = reactor.control.input['mat'][icool]['temp0']
            # list of coolant types in pipe
            self.type.append(type)
            # list of initial pressures in pipe nodes
            self.p.append([p0]*self.pipennodes[i])
            for j in range(self.pipennodes[i]):
                self.p0.append(p0)
            # list of initial temperatures in pipe nodes
            self.temp.append([temp0]*self.pipennodes[i])
            # list of initial velocities in pipe nodes
            self.vel.append([0]*self.pipennodes[i])
            # list of initial reynolds in pipe nodes
            self.re.append([0]*self.pipennodes[i])
            # list of initial prandtl in pipe nodes
            self.pr.append([0]*self.pipennodes[i])
            # list of initial peclet in pipe nodes
            self.pe.append([0]*self.pipennodes[i])
        # assign index to every pipe node
        self.indx = []
        for i in range(self.npipe):
            for j in range(self.pipennodes[i]):
                self.indx.append((i,j))

        # list of junction types and subtypes
        self.juntype = reactor.control.input['junction']['type']
        # number of junctions
        self.njun = len(self.juntype)
        # number of independent junctions
        self.njuni = self.juntype.count('independent')
        # number of dependent junctions
        self.njund = self.juntype.count('dependent')

        # construct from and to lists of tulips
        self.f = []
        self.t = []
        for j in range(self.njun):
            idf = reactor.control.input['junction']['from'][j]
            idt = reactor.control.input['junction']['to'][j]
            if idf not in self.pipeid:
                print('****ERROR: pipe id (' + idf + ') in junction (' + idf + '-' + idt + ') does not exist in input.')
                sys.exit()
            if idt not in self.pipeid:
                print('****ERROR: pipe id (' + idt + ') in junction (' + idf + '-' + idt + ') does not exist in input.')
                sys.exit()
            indx = self.pipeid.index(idf)
            self.f.append((indx, self.pipennodes[indx]-1))
            indx = self.pipeid.index(idt)
            self.t.append((indx, 0))
        # add internal junctions
        for i in range(self.npipe):
            self.juntype += ['internal']*(self.pipennodes[i] - 1)
            # append from and to lists
            self.f += [(i,j) for j in range(self.pipennodes[i]-1)]
            self.t += [(i,j) for j in range(1,self.pipennodes[i])]
        self.njun = len(self.f)
        # user-specified junction pump head signal
        self.junpumphead = reactor.control.input['junpumphead']
        # user-specified junction flowrate signal
        self.junflowrate = reactor.control.input['junflowrate']
        # user-specified junction k-factor signal
        self.junkfac = reactor.control.input['junkfac']

        # create and inverse a matrix A linking dependent and independent junctions
        A = [[0]*(self.njuni+self.njund) for i in range(self.njuni+self.njund)]
        i = 0
        for j in range(self.njuni+self.njund):
            if self.juntype[j] == 'independent':
                A[j][j] = 1
            elif self.juntype[j] == 'dependent':
                while self.pipetype[i] == 'freelevel':
                    i += 1
                for jj in range(self.njuni+self.njund):
                    if self.f[jj][0] == i:
                        A[j][jj] = -1
                    if self.t[jj][0] == i:
                        A[j][jj] = 1
                i += 1
        self.invA = linalg.inv(A)

        # create list of geometrical parameters for every junction: l_over_a = 0.5*length(from)/areaz(from) + 0.5*length(to)/areaz(to)
        l_over_a = [0]*self.njun
        for j in range(self.njun):
            f = self.f[j][0]
            t = self.t[j][0]
            len_f = 0.5*self.len[f]/self.pipennodes[f]
            if self.pipetype[f] == 'freelevel': len_f *= 2
            len_t = 0.5*self.len[t]/self.pipennodes[t]
            if self.pipetype[t] == 'freelevel': len_t *= 2

            l_over_a[j] = len_f/self.areaz[f] + len_t/self.areaz[t]

        # create and invert a matrix B of left-hand sides of momentum conservation equations (self.njun) and 
        # mass conservation equations differentiated w.r.t. time (self.npipe-self.npipef)
        n = self.njun + sum(self.pipennodes)
        B = [[0]*n for i in range(n)]
        for j in range(self.njun):
            B[j][j] = l_over_a[j] # dmdot/dt

            i = self.njun + self.indx.index(self.f[j])
            B[j][i] = -1 # -P_from
            if self.pipetype[self.f[j][0]] != 'freelevel':
                B[i][j] = -1 # -dmdot/dt_out

            i = self.njun + self.indx.index(self.t[j])
            B[j][i] = 1 # +P_to
            if self.pipetype[self.t[j][0]] != 'freelevel':
                B[i][j] = 1 # +dmdot/dt_in
        for i in range(sum(self.pipennodes)):
            if self.pipetype[self.indx[i][0]] == 'freelevel':
                B[self.njun + i][self.njun + i] = 1 # P = +_freelevel
        self.invB = linalg.inv(B)

        # initialize list of flowrate in independent junctions
        self.mdoti = [0]*self.njuni

        # initialize list of flowrate in all junctions
        self.mdot = [0]*self.njun

        # prepare lists to map thermal-hydraulic and fuel-rod nodes
        indx_i = [x['pipeid'] for x in reactor.control.input['fuelrod']]
        indx_j = [x['pipenode'] for x in reactor.control.input['fuelrod']]
        nfuelrods = len(indx_i)
        self.map_th = []
        self.map_fr = []
        for i in range(nfuelrods):
            for j in range(len(indx_i[i])):
                self.map_th.append((indx_i[i][j], indx_j[i][j]-1))
                self.map_fr.append((i,j))

        # analyse thermal boundary conditions
        for x in reactor.control.input['thermbc']:
            # boundary with pipe
            if x['type'] == 2:
                # check if pipe exists
                jpipe = (x['pipeid'], x['pipenode'])
                if jpipe[0] not in self.pipeid:
                    print('****ERROR: pipe id (' + jpipe[0] + ') in \'thermbc\' card (' + x['id'] + ') does not exist in input.')
                    sys.exit()
                else:
                    # pipe index
                    ipipe = self.pipeid.index(jpipe[0])
                # check if pipenodeid exists
                if jpipe[1] > self.pipennodes[ipipe]:
                    print('****ERROR: pipe node index (' + str(jpipe[1]) + ') given in \'thermbc\' card (' + x['id'] + ') exceeds number of nodes (' + str(self.pipennodes[ipipe]) + ') of pipe ' + jpipe[0])
                    sys.exit()
    #----------------------------------------------------------------------------------------------
    # create right-hand side list: self is a 'fluid' object created in B
    def calculate_rhs(self, reactor, t):

        if 'fluid' not in reactor.solve:
            rhs = []
            return rhs

        # FLUID PROPERTIES:
        self.prop = []
        for i in range(self.npipe):
            dict = {'rhol':[], 'visl':[], 'kl':[], 'cpl':[]}
            for j in range(self.pipennodes[i]):
                # call material property function
                pro = reactor.data.matpro( {'type':self.type[i], 't':self.temp[i][j], 'p':self.p[i][j]} )
                dict['rhol'].append(pro['rhol'])
                dict['visl'].append(pro['visl'])
                dict['kl'].append(pro['kl'])
                dict['cpl'].append(pro['cpl'])
            self.prop.append(dict)

        # FLOWRATES IN DEPENDENT JUNCTIONS:
        # construct right hand side of system invA*mdot = b
        i = 0
        b = [0]*(self.njuni+self.njund)
        for j in range(self.njuni+self.njund):
            if self.juntype[j] == 'independent':
                b[j] = self.mdoti[i]
                i += 1
            elif self.juntype[j] == 'dependent':
                b[j] = 0        
        # then multiply matrix by list: invA*mdot = b and convert to list
        self.mdot = self.invA.dot(b).tolist()
        # finally calculate flowrates in internal junctions
        for i in range(self.npipe):
            mdotpipe = 0
            for j in range(self.njuni+self.njund):
                if self.t[j][0] == i:
                    mdotpipe += self.mdot[j]
            for j in range(self.pipennodes[i]-1):
                self.mdot.append(mdotpipe)

        # VELOCITIES AND DIMENSIONLESS NUMBERS IN PIPE NODES:
        self.vel = [[0]*self.pipennodes[i] for i in range(self.npipe)]
        for j in range(self.njun):
            rho_t = self.prop[self.t[j][0]]['rhol'][self.t[j][1]]
            self.vel[self.t[j][0]][self.t[j][1]] += self.mdot[j]/rho_t/self.areaz[self.t[j][0]]
        # Reynolds numbers
        self.re = [[abs(self.vel[i][j])*self.dhyd[i]/self.prop[i]['visl'][j] for j in range(self.pipennodes[i])] for i in range(self.npipe)]
        # Prandtl numbers
        self.pr = [[self.prop[i]['visl'][j]*self.prop[i]['rhol'][j]*self.prop[i]['cpl'][j]/self.prop[i]['kl'][j] for j in range(self.pipennodes[i])] for i in range(self.npipe)]
        # Peclet numbers
        self.pe = [[self.re[i][j]*self.pr[i][j] for j in range(self.pipennodes[i])] for i in range(self.npipe)]

        # TIME DERIVATIVES OF MASS FLOWRATES:
        # construct right-hand side b of system invB*[dmdotdt, P] = b
        b = [0]*(self.njun + sum(self.pipennodes))
        for j in range(self.njun):
            f = self.f[j]
            t = self.t[j]
            len_f = 0.5*self.len[f[0]]/self.pipennodes[f[0]]
            if self.pipetype[f[0]] == 'freelevel': len_f *= 2
            len_t = 0.5*self.len[t[0]]/self.pipennodes[t[0]]
            if self.pipetype[t[0]] == 'freelevel': len_t *= 2
            rho_f = self.prop[f[0]]['rhol'][f[1]]
            rho_t = self.prop[f[0]]['rhol'][t[1]]

            # gravitational head
            if self.pipetype[f[0]] != 'freelevel': 
                rhogh_f = 9.81*rho_f*len_f*self.dir[f[0]]
            else:
                rhogh_f = 9.81*rho_f*len_f*self.dir[t[0]]
            if self.pipetype[t[0]] != 'freelevel': 
                rhogh_t = 9.81*rho_t*len_t*self.dir[t[0]]
            else:
                rhogh_t = 9.81*rho_t*len_t*self.dir[f[0]]

            # friction losses for from
            if self.pipetype[f[0]] == "wirewrapped":
                p2d = reactor.control.input['pipe'][f[0]]["p2d"]
                h2d = reactor.control.input['pipe'][f[0]]["h2d"]
                inp = {'p2d':p2d,'h2d':h2d,'re':self.re[f[0]][f[1]]}
                dpfric_f = reactor.data.fricfac(inp) * 0.5 * rho_f * self.vel[f[0]][f[1]] * abs(self.vel[f[0]][f[1]])
            else:
                inp = {'re':self.re[f[0]][f[1]]}
                dpfric_f = reactor.data.fricfac(inp) * 0.5 * rho_f * self.vel[f[0]][f[1]] * abs(self.vel[f[0]][f[1]])

            # friction losses for to
            if self.pipetype[t[0]] == "wirewrapped":
                p2d = reactor.control.input['pipe'][t[0]]["p2d"]
                h2d = reactor.control.input['pipe'][t[0]]["h2d"]
                inp = {'p2d':p2d,'h2d':h2d,'re':self.re[t[0]][t[1]]}
                dpfric_t = reactor.data.fricfac(inp) * 0.5 * rho_t * self.vel[t[0]][t[1]] * abs(self.vel[t[0]][t[1]])            
            else:
                inp = {'re':self.re[t[0]][t[1]]}
                dpfric_t = reactor.data.fricfac(inp) * 0.5 * rho_t * self.vel[t[0]][t[1]] * abs(self.vel[t[0]][t[1]])

            b[j] = -(rhogh_f + rhogh_t) - (dpfric_f + dpfric_t)
            try:
                # check if the current junction j is present in junkfac list
                indx = reactor.fluid.junkfac['jun'].index((self.pipeid[f[0]],self.pipeid[t[0]]))
                # if yes...
                kfac = reactor.control.signal[self.junkfac['kfac'][indx]]
                # local (singular) pressure losses
                if self.mdot[j] > 0:
                    b[j] -= kfac * rho_f * self.vel[f[0]][f[1]]**2 / 2.0
                else:
                    b[j] += kfac * rho_t * self.vel[t[0]][t[1]]**2 / 2.0
                
            except ValueError:
                pass
            
            if self.juntype[j] == 'independent':
                f = reactor.fluid.f[j][0]
                t = reactor.fluid.t[j][0]
                # tuple of from-to pipe id's
                f_t = (reactor.fluid.pipeid[f],reactor.fluid.pipeid[t])
                # check if the current junction j is present in junpumphead list
                if f_t in reactor.control.input['junpumphead']['jun']:
                    indx = reactor.control.input['junpumphead']['jun'].index(f_t)
                    b[j] += reactor.control.signal[self.junpumphead['pumphead'][indx]]

        for i in range(self.npipe):
            for j in range(self.pipennodes[i]):
                self.indx.append((i,j))
        for i in range(sum(self.pipennodes)):
            n = self.indx[i][0]
            if self.pipetype[n] == 'freelevel': b[self.njun+i] = self.p0[i]
        invBb = self.invB.dot(b).tolist()

        # read from invBb: time derivatives of flowrate in independent junctions
        dmdotdt = []
        for j in range(self.njun):
            if self.juntype[j] == 'independent':
                f = reactor.fluid.f[j][0]
                t = reactor.fluid.t[j][0]
                # tuple of from-to pipe id's
                f_t = (reactor.fluid.pipeid[f],reactor.fluid.pipeid[t])
                try:
                    # check if the current junction j is present in junflowrate list
                    indx = reactor.fluid.junflowrate['jun'].index(f_t)
                    # if yes...
                    dmdotdt.append(0)
                except ValueError:
                    dmdotdt.append(invBb[j])
        # read from invBb: pressures in pipe nodes
        indx = 0
        for i in range(self.npipe):
            for j in range(self.pipennodes[i]):
                self.p[i][j] = invBb[self.njun+indx]
                indx += 1

        # TIME DERIVATIVES OF FREE-LEVEL-VOLUME LENGTH:
        dlendt = []
        n = 0
        for i in range(self.npipe):
            if self.pipetype[i] == 'freelevel':
                dlendt.append(0)
                rhoa = self.prop[i]['rhol'][0]*self.areaz[i]
                for j in range(self.njun):
                    if self.f[j][0] == i:
                        dlendt[n] -= self.mdot[j]/rhoa
                    if self.t[j][0] == i:
                        dlendt[n] += self.mdot[j]/rhoa
                n += 1
        
        # TIME DERIVATIVES OF FLUID TEMPERATURES:
        dtempdt2d = [[0]*self.pipennodes[i] for i in range(self.npipe)]
        for j in range(self.njun):
            f = self.f[j]
            t = self.t[j]
            if self.mdot[j] > 0:
                cp_temp_mdot = self.prop[f[0]]['cpl'][f[1]] * self.temp[f[0]][f[1]] * self.mdot[j]
            else:
                cp_temp_mdot = self.prop[t[0]]['cpl'][t[1]] * self.temp[t[0]][t[1]] * self.mdot[j]
            dtempdt2d[f[0]][f[1]] -= cp_temp_mdot
            dtempdt2d[t[0]][t[1]] += cp_temp_mdot
        n = 0
        for i in range(self.npipe):
            if self.pipetype[i] == 'freelevel':
                dtempdt2d[i][0] -= self.prop[i]['cpl'][0] *  self.temp[i][0] * dlendt[n] * self.prop[i]['rhol'][0] * self.areaz[i]
                n += 1

        dtempdt = []
        for i in range(self.npipe):
            if self.signaltemp[i] != '':
                for j in range(self.pipennodes[i]):
                    dtempdt.append(0)
            else:
                len = abs(self.len[i])/self.pipennodes[i]
                vol = self.areaz[i] * len
                for j in range(self.pipennodes[i]):
                    # check if there is a fuel rod cooled by the node
                    if 'fuelrod' in reactor.solve and (self.pipeid[i],j) in self.map_th:
                        indx = self.map_th.index((self.pipeid[i],j))
                        tuple_fr = self.map_fr[indx]
                        tclad = reactor.solid.fuelrod[tuple_fr[0]].clad[tuple_fr[1]].temp[-1]
                        mltpl = reactor.solid.fuelrod[tuple_fr[0]].clad[tuple_fr[1]].mltpl
                        ro = reactor.solid.fuelrod[tuple_fr[0]].clad[tuple_fr[1]].r[-1]
                        area_ht = 2 * math.pi * ro * len * mltpl
                        p2d = reactor.solid.fuelrod[tuple_fr[0]].clad[tuple_fr[1]].p2d
                        nu = reactor.data.nu( {'type':self.type[i],'re':self.re[i][j], 'pr':self.pr[i][j],'p2d':p2d} )
                        hex = nu * self.prop[i]['kl'][j] / self.dhyd[i]
                        dtempdt2d[i][j] += hex*(tclad - self.temp[i][j]) * area_ht

                    # check if there is a heat structure cooled by the node
                    if 'htstr' in reactor.solve :
                        for k in range(reactor.solid.nhtstr):
                            bcleft = reactor.solid.htstr[k].bcleft
                            bcright = reactor.solid.htstr[k].bcright
                            if bcleft['type'] == 2 and bcleft['pipeid'] == self.pipeid[i] and bcleft['pipenode']-1 == j:
                                nu = reactor.data.nu( {'type':self.type[i],'re':self.re[i][j], 'pr':self.pr[i][j]} )
                                hex = nu * self.prop[i]['kl'][j] / self.dhyd[i]
                                mltpl = reactor.solid.htstr[k].mltpl
                                area_ht = 2 * math.pi * reactor.solid.htstr[k].ri * len * mltpl
                                dtempdt2d[i][j] += hex*(reactor.solid.htstr[k].temp[0] - self.temp[i][j]) * area_ht
                            if bcright['type'] == 2 and bcright['pipeid'] == self.pipeid[i] and bcright['pipenode']-1 == j:
                                nu = reactor.data.nu( {'type':self.type[i],'re':self.re[i][j], 'pr':self.pr[i][j]} )
                                hex = nu * self.prop[i]['kl'][j] / self.dhyd[i]
                                mltpl = reactor.solid.htstr[k].mltpl
                                area_ht = 2 * math.pi * reactor.solid.htstr[k].ro * len * mltpl
                                dtempdt2d[i][j] += hex*(reactor.solid.htstr[k].temp[-1] - self.temp[i][j]) * area_ht

                    rho_cp_vol = self.prop[i]['rhol'][j] * self.prop[i]['cpl'][j] * vol
                    dtempdt.append(dtempdt2d[i][j] / rho_cp_vol)

        rhs = dmdotdt + dlendt + dtempdt
        return rhs
