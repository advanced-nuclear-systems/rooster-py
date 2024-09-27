import numpy as np
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
        # list of user-specified pipe temporary signal
        self.signaltemp = [x['signaltemp'] for x in reactor.control.input['pipe']]

        # psi_ytchen: lists for pressure, temperature, type, velocity, reynolds, prandtl, peclet and thermal boundary condition id
        self.p0 = []
        self.p = []
        self.temp = []
        self.h = []    # psi_ytchen: fluid enthalpy list
        self.alpha = [] # psi_ytchen: fluid void fraction list
        self.xe = [] # psi_ytchen: fluid gas quality list
        self.type = []
        self.vel = []
        self.re = []
        self.pr = []
        self.pe = []
        self.Gf = [] # psi_ytchen: fluid total mass flux, kg/(m2*s)
        self.htcmod = [] # psi_ytchen: fluid heat transfer mode
        self.xlm = [] # psi_ytchen: Lockhart-Martinelli
        
        # psi_ytchen: list for pressure drop calculation issues
        self.dpkne  = []
        self.dpfric = []
        self.dpgrav = []
        
        
        
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
            # + psi_ytchen: list of initial enthalpies in pipe nodes
            pro = reactor.data.matpro( {'type':type, 't':temp0, 'p':p0, 'ini':1} )
            h = pro['h']
            self.h.append([h]*self.pipennodes[i])
            # - psi_ytchen: list of initial enthalpies in pipe nodes
            # list of initial velocities in pipe nodes
            self.vel.append([0]*self.pipennodes[i])
            # list of initial reynolds in pipe nodes
            self.re.append([0]*self.pipennodes[i])
            # list of initial prandtl in pipe nodes
            self.pr.append([0]*self.pipennodes[i])
            # list of initial peclet in pipe nodes
            self.pe.append([0]*self.pipennodes[i])
            # list of void fraction in pipe nodes
            self.alpha.append([0]*self.pipennodes[i])
            # list of gas quality in pipe nodes
            self.xe.append([0]*self.pipennodes[i])
            # psi_ytchen: list of initial fluid total mass flux
            self.Gf.append([0]*self.pipennodes[i])
            # psi_ytchen: list of initial fluid total mass flux
            self.htcmod.append([0]*self.pipennodes[i])
            # psi_ytchen: list of initial xlm number
            self.xlm.append([1]*self.pipennodes[i])
            # psi_ytchen: list of fluid kinetic energy 
            self.dpkne.append([0]*self.pipennodes[i])
            # psi_ytchen: list of frictional pressure drop
            self.dpfric.append([0]*self.pipennodes[i])
            # psi_ytchen: list of gravitational pressure drop
            self.dpgrav.append([0]*self.pipennodes[i])
            
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
        self.invA = np.linalg.inv(A)

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
            f_t = (self.pipeid[self.f[j][0]],self.pipeid[self.t[j][0]])

            B[j][j] = l_over_a[j] # dmdot/dt, psi_ytchen: l.h.s coefficient for Momentum Equation

            i = self.njun + self.indx.index(self.f[j])
            B[j][i] = -1 # -P_from, psi_ytchen: l.h.s coefficient for Momentum Equation
            if f_t in reactor.control.input['junflowrate']['jun']:
                B[j][i] = 0.0
            if self.pipetype[self.f[j][0]] != 'freelevel':
                B[i][j] = -1 # -dmdot/dt_out, psi_ytchen: l.h.s coefficient for Pressure Equation

            i = self.njun + self.indx.index(self.t[j])
            B[j][i] = 1 # +P_to, psi_ytchen: l.h.s coefficients for Momentum Equation
            if f_t in reactor.control.input['junflowrate']['jun']:
                B[j][i] = 0.0
            if self.pipetype[self.t[j][0]] != 'freelevel':
                B[i][j] = 1 # +dmdot/dt_in, psi_ytchen: l.h.s coefficient for Pressure Equation

        for i in range(sum(self.pipennodes)):
            if self.pipetype[self.indx[i][0]] == 'freelevel':
                B[self.njun + i][self.njun + i] = 1 # P = +_freelevel, psi_ytchen: freelevel pressure boundary for Pressure Equation
        self.invB = np.linalg.inv(B)

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
            dict = {'rhof':[], 'visf':[], 'kf':[], 'cpf':[]}
            # call material property function
            # + psi_ytchen: changed into h-type fluid properties
            pro = reactor.data.matpro( {'type':self.type[i], 'h':self.h[i], 'p':self.p[i], 'ini':0} )
            self.temp[i] = pro['tf']
            # - psi_ytchen: changed into h-type fluid properties
            dict['rhof'] = pro['rhof']
            dict['visf'] = pro['visf']
            dict['kf'] = pro['kf']
            dict['cpf'] = pro['cpf']
            if self.type[i] == 'h2o' or self.type[i] == 'na':
                dict['ts']  = pro['ts']
                dict['xe']  = pro['xe']
                dict['hgl'] = pro['hgl']
                dict['sgm'] = pro['sgm']
                dict['rhol']= pro['rhol']
                dict['rhog']= pro['rhog']
                dict['miul']= pro['miul']
                dict['miug']= pro['miug']
                dict['cpg'] = pro['cpg'] 
                dict['kg']  = pro['kg']
                dict['cpl'] = pro['cpl'] 
                dict['kl']  = pro['kl']
                dict['miuf']= pro['miuf']
 
                self.xe[i] = pro['xe']
                # psi_ytchen: calculate alpha 
                xe_np = np.array(pro['xe'])
                rhof_np = np.array(pro['rhof'])
                rhog_np = np.array(pro['rhog'])
                alpha_np = xe_np*(rhof_np/rhog_np)
                alpha_np[alpha_np < 0.0] = 0
                alpha_np[alpha_np > 1.0] = 1.0
                self.alpha[i] = list(alpha_np)
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
            rho_t = self.prop[self.t[j][0]]['rhof'][self.t[j][1]]
            self.vel[self.t[j][0]][self.t[j][1]] += self.mdot[j]/rho_t/self.areaz[self.t[j][0]]
        # MassFlux
        self.Gf=  [[abs(self.vel[i][j])*self.prop[i]['rhof'][j] for j in range(self.pipennodes[i])] for i in range(self.npipe)]
        # Reynolds numbers
        self.re = [[abs(self.vel[i][j])*self.dhyd[i]/self.prop[i]['visf'][j] for j in range(self.pipennodes[i])] for i in range(self.npipe)]
        # Prandtl numbers
        self.pr = [[self.prop[i]['visf'][j]*self.prop[i]['rhof'][j]*self.prop[i]['cpf'][j]/self.prop[i]['kf'][j] for j in range(self.pipennodes[i])] for i in range(self.npipe)]
        # Peclet numbers
        self.pe = [[self.re[i][j]*self.pr[i][j] for j in range(self.pipennodes[i])] for i in range(self.npipe)]
        
        # + psi_ytchen: Calculate dpkne dpfric dpgrav
        for i in range(self.npipe):
            velnp = np.array(self.vel[i])
            renp  = np.array(self.re[i])
            rhonp = np.array(self.prop[i]['rhof'])
            dpkne = rhonp*velnp*np.abs(velnp)
            dpgrav = 9.81*rhonp*self.len[i]*self.dir[i]/self.pipennodes[i]
            if self.pipetype[i] == 'wirewrapped':
                p2d = reactor.control.input['pipe'][i]["p2d"]
                h2d = reactor.control.input['pipe'][i]["h2d"]
                sb2st=reactor.control.input['pipe'][i]["sb2st"]
                inp = {'p2d':p2d,'h2d':h2d,'re':self.re[i],'sb2st':sb2st}
            elif self.pipetype[i] == 'barerod':
                p2d = reactor.control.input['pipe'][i]["p2d"]
                sb2st=reactor.control.input['pipe'][i]["sb2st"]
                inp = {'p2d':p2d,'re':self.re[i], 'sb2st':sb2st}  # - psi_ytchen: bare rod bundle inputs 
            else:
                inp = {'re':self.re[i]}
                
            if self.type[i] != 'na':
                dpfric = []
                fricnp = reactor.data.fricfac(inp)
                dpfric = fricnp*0.5*dpkne*self.len[i]/self.pipennodes[i]/self.dhyd[i]
            else: # sodium flow condition
                dpfric = []
                for j in range(self.pipennodes[i]):
                    if self.xe[i][j] < 1.0 and self.xe[i][j] > 0.0: # two-phase sodium cell
                        miul= self.prop[i]['miul'][j]
                        miug= self.prop[i]['miug'][j]
                        miuf= self.prop[i]['miuf'][j]
                        xe  = self.xe[i][j]
                        rel = renp[j]*miuf/miul
                        reg = renp[j]*np.array(self.prop[i]['miuf'][j])/np.array(self.prop[i]['miug'][j])
                        rhol= self.prop[i]['rhol'][j]
                        rhog= self.prop[i]['rhog'][j]
                        if self.pipetype[i] == 'wirewrapped':
                            inpl= {'p2d':p2d,'h2d':h2d,'re':rel,'sb2st':sb2st}
                            inpg= {'p2d':p2d,'h2d':h2d,'re':reg,'sb2st':sb2st}
                        elif self.pipetype[i] == 'barerod':                # + psi_ytchen: bare rod bundle inputs
                            inpl= {'p2d':p2d,'re':rel, 'sb2st':sb2st}              # - psi_ytchen: bare rod bundle inputs           
                            inpg= {'p2d':p2d,'re':reg, 'sb2st':sb2st}
                        else:
                            inpl= {'re':rel}
                            inpg= {'re':reg}
                        fricl = reactor.data.fricfac(inpl)
                        fricg = reactor.data.fricfac(inpg)
                        #XLM= (1.0/self.xe[i][j] - 1.0)**(0.9)*(rhog/rhol)**(0.5)*(miul/miug)**(0.1)
                        GG = rhonp[j]*velnp[j]
                        #FL2_yd= 1. + 15.93/XLM + 1./XLM**(2.0) # Yandong HOU correlation
                        FL2 = ( 1. + 4.*xe*(1.-xe) )*( xe*rhol/rhog + 1. - xe )**(0.853) # Martinelli-Nelson correlation
                        dpf_mn = FL2*fricl*self.len[i]/self.pipennodes[i]/self.dhyd[i]*GG*abs(GG)/(2.0*rhol)
                        dpf_temp =  dpf_mn
                        #dpf_g  = fricg*self.len[i]/self.pipennodes[i]/self.dhyd[i]*GG*abs(GG)/(2.0*rhog)
                        #xcr = 0.30
                        #if xe < xcr:
                        #   dpf_temp =  dpf_mn
                        #else:
                            #dx = (xe-xcr)/(1.0-xcr)
                            #phy= math.exp(dx/(dx-1.0))
                            #dpf_temp =  dpf_mn*phy + dpf_g*(1.0-phy)
                        #self.xlm[i][j] = XLM # record the Lockhart-Martinelli number
                    else:  # single-phase sodium cell
                        if self.pipetype[i] == 'wirewrapped':
                            inp = {'p2d':p2d,'h2d':h2d,'re':self.re[i][j],'sb2st':sb2st}
                        elif self.pipetype[i] == 'barerod':                # + psi_ytchen: bare rod bundle inputs
                            inp = {'p2d':p2d,'re':self.re[i][j],'sb2st':sb2st}     # - psi_ytchen: bare rod bundle inputs
                        else:
                            inp = {'re':self.re[i][j]}
                        fric = reactor.data.fricfac(inp)
                        dpf_temp = fric*0.5*dpkne[j]*self.len[i]/self.pipennodes[i]/self.dhyd[i]
                    dpfric.append(dpf_temp)
            self.dpkne[i]  = list(dpkne)
            self.dpfric[i] = list(dpfric)
            self.dpgrav[i] = list(dpgrav)
        # - psi_ytchen: 
        
        # TIME DERIVATIVES OF MASS FLOWRATES:
        # construct right-hand side b of system invB*[dmdotdt, P] = b
        b = [0]*(self.njun + sum(self.pipennodes))
        for j in range(self.njun):
            f = self.f[j]
            t = self.t[j]
            # + psi_ytchen: fluid kinetic energy
            dpkne_f = self.dpkne[f[0]][f[1]]
            dpkne_t = self.dpkne[t[0]][t[1]]
            # gravitational head
            if self.pipetype[f[0]] != 'freelevel': 
                rhogh_f = self.dpgrav[f[0]][f[1]]/2.0
            else:
                rhogh_f = self.dpgrav[f[0]][f[1]]
            if self.pipetype[t[0]] != 'freelevel': 
                rhogh_t = self.dpgrav[t[0]][t[1]]/2.0
            else:
                rhogh_t = self.dpgrav[t[0]][t[1]]

            # friction losses for from
            dpfric_f = self.dpfric[f[0]][f[1]]/2.0
            # friction losses for to
            dpfric_t = self.dpfric[t[0]][t[1]]/2.0
            b[j] = - (rhogh_f + rhogh_t) - (dpfric_f + dpfric_t)
            # psi_ytchen: special treatment to freelevel cell to avoid weird pressure field
            if self.pipetype[f[0]] != 'freelevel' and self.pipetype[t[0]] != 'freelevel': 
                b[j] += 0.0 # dpkne_f - dpkne_t # psi_ytchen: accleration pressure drop
            try:
                # check if the current junction j is present in junkfac list
                indx = reactor.fluid.junkfac['jun'].index((self.pipeid[f[0]],self.pipeid[t[0]]))
                # if yes...
                kfac = reactor.control.signal[self.junkfac['kfac'][indx]]
                # local (singular) pressure losses
                if self.mdot[j] > 0:
                    b[j] -= abs(kfac * dpkne_f / 2.0)
                else:
                    b[j] += abs(kfac * dpkne_t / 2.0)
                    
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
                if f_t in reactor.control.input['junflowrate']['jun']:
                    b[j] = 0.0


        # for i in range(self.npipe):                  #psi_ytchen: dangerous memory leakage point
        #     for j in range(self.pipennodes[i]):      #psi_ytchen: dangerous memory leakage point
        #         self.indx.append((i,j))              #psi_ytchen: dangerous memory leakage point
        for i in range(sum(self.pipennodes)):
            n = self.indx[i][0]
            if self.pipetype[n] == 'freelevel':
                b[self.njun+i] = reactor.control.signal[reactor.fluid.signaltemp[n]]
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
                rhoa = self.prop[i]['rhof'][0]*self.areaz[i]
                for j in range(self.njun):
                    if self.f[j][0] == i:
                        dlendt[n] -= self.mdot[j]/rhoa
                    if self.t[j][0] == i:
                        dlendt[n] += self.mdot[j]/rhoa
                n += 1
        
        # TIME DERIVATIVES OF FLUID ENTHALPIES: psi_ytchen
        dtempdt2d = [[0]*self.pipennodes[i] for i in range(self.npipe)]
        for j in range(self.njun):
            f = self.f[j]
            t = self.t[j]
            # + psi_ytchen: 'cp_temp_mdot' replaced by 'enth_mdot', cp*T replaced by h
            if self.mdot[j] > 0:
                enth_mdot = self.h[f[0]][f[1]] * self.mdot[j]
            else:
                enth_mdot = self.h[t[0]][t[1]] * self.mdot[j]
            dtempdt2d[f[0]][f[1]] -= enth_mdot
            dtempdt2d[t[0]][t[1]] += enth_mdot
        n = 0
        for i in range(self.npipe):
            if self.pipetype[i] == 'freelevel':
                dtempdt2d[i][0] -= self.h[i][0] * dlendt[n] * self.prop[i]['rhof'][0] * self.areaz[i]
                n += 1
            # - psi_ytchen: 'cp_temp_mdot' replaced by 'enth_mdot', cp*T replaced by h
        # + psi_ytchen: 'dtempdt' replaced by 'denthdt'
        denthdt = []
        for i in range(self.npipe):
            if self.signaltemp[i] != '' and self.pipetype[i] != 'freelevel':
                for j in range(self.pipennodes[i]):
                    denthdt.append(0)
            else:
                len = abs(self.len[i])/self.pipennodes[i]
                vol = self.areaz[i] * len
                for j in range(self.pipennodes[i]):
                    # + psi_ytchen: prepare fluid property dictionary
                    dictfld = {'rhof':0.0, 'visf':0.0, 'kf':0.0, 'cpf':0.0} 
                    dictfld['rhof'] = self.prop[i]['rhof'][j]
                    dictfld['visf'] = self.prop[i]['visf'][j]
                    dictfld['kf'] = self.prop[i]['kf'][j]
                    dictfld['cpf'] = self.prop[i]['cpf'][j]
                    tfluid = self.temp[i][j]
                    pcell  = self.p[i][j]
                    if self.type[i] == 'h2o' or self.type[i] == 'na':
                        dictfld['ts']  = self.prop[i]['ts'][j]
                        dictfld['xe']  = self.prop[i]['xe'][j]
                        dictfld['hgl'] = self.prop[i]['hgl'][j]
                        dictfld['sgm'] = self.prop[i]['sgm'][j]
                        dictfld['rhol']= self.prop[i]['rhol'][j]
                        dictfld['rhog']= self.prop[i]['rhog'][j]
                        dictfld['miul']= self.prop[i]['miul'][j]
                        dictfld['miug']= self.prop[i]['miug'][j]
                        dictfld['cpg'] = self.prop[i]['cpg'][j]
                        dictfld['kg']  = self.prop[i]['kg'][j]
                        dictfld['cpl'] = self.prop[i]['cpl'][j]
                        dictfld['kl']  = self.prop[i]['kl'][j]
                        
                        
                    # check if there is a fuel rod cooled by the node
                    if 'fuelrod' in reactor.solve and (self.pipeid[i],j) in self.map_th:
                        # + psi_ytchen: boiling heat transfer model is added in this section
                        indx = self.map_th.index((self.pipeid[i],j))
                        tuple_fr = self.map_fr[indx]
                        tclad = reactor.solid.fuelrod[tuple_fr[0]].clad[tuple_fr[1]].temp[-1]
                        mltpl = reactor.solid.fuelrod[tuple_fr[0]].clad[tuple_fr[1]].mltpl
                        ro = reactor.solid.fuelrod[tuple_fr[0]].clad[tuple_fr[1]].r[-1]
                        area_ht = 2 * math.pi * ro * len * mltpl
                        p2d = reactor.solid.fuelrod[tuple_fr[0]].clad[tuple_fr[1]].p2d
                        inp4qflux = {'type':self.type[i],'re':self.re[i][j], 'pr':self.pr[i][j],'p2d':p2d, 'dhyd':self.dhyd[i],\
                        'tmp':[tfluid, tclad], 'p':pcell, 'prop':dictfld, 'Gtot':abs(self.Gf[i][j]) }
                        htc = reactor.data.qfluxcal( inp4qflux )
                        qflux = htc[0]
                        self.htcmod[i][j] = htc[1]
                        dtempdt2d[i][j] += qflux * area_ht # psi_ytchen
                        reactor.solid.fuelrod[tuple_fr[0]].clad[tuple_fr[1]].qfluxo = qflux # psi_ytchen: store qflux for Clad.calculate_rhs
                    # check if there is a heat structure cooled by the node
                    if 'htstr' in reactor.solve :
                        for k in range(reactor.solid.nhtstr):
                            bcleft = reactor.solid.htstr[k].bcleft
                            bcright = reactor.solid.htstr[k].bcright
                            mltpl = reactor.solid.htstr[k].mltpl
                            if bcleft['type'] == 2 and bcleft['pipeid'] == self.pipeid[i] and bcleft['pipenode']-1 == j:
                                area_ht = 2 * math.pi * reactor.solid.htstr[k].ri * len * mltpl
                                twall = reactor.solid.htstr[k].temp[0]
                                inp4qflux = {'type':self.type[i],'re':self.re[i][j], 'pr':self.pr[i][j], 'dhyd':self.dhyd[i],\
                                'tmp':[tfluid, twall], 'p':pcell, 'prop':dictfld, 'Gtot':abs(self.Gf[i][j]) }
                                htc = reactor.data.qfluxcal( inp4qflux )
                                qflux = htc[0]
                                self.htcmod[i][j] = htc[1]
                                dtempdt2d[i][j] += qflux * area_ht
                                reactor.solid.htstr[k].qfluxleft = qflux # psi_ytchen
                            if bcright['type'] == 2 and bcright['pipeid'] == self.pipeid[i] and bcright['pipenode']-1 == j:
                                area_ht = 2 * math.pi * reactor.solid.htstr[k].ro * len * mltpl
                                twall = reactor.solid.htstr[k].temp[-1]
                                inp4qflux = {'type':self.type[i],'re':self.re[i][j], 'pr':self.pr[i][j], 'dhyd':self.dhyd[i],\
                                'tmp':[tfluid, twall], 'p':pcell, 'prop':dictfld, 'Gtot':abs(self.Gf[i][j]) }
                                htc = reactor.data.qfluxcal( inp4qflux )
                                qflux = htc[0]
                                self.htcmod[i][j] = htc[1]
                                dtempdt2d[i][j] += qflux * area_ht # psi_ytchen
                                reactor.solid.htstr[k].qfluxright = qflux # psi_ytchen
                    rho_vol = self.prop[i]['rhof'][j] * vol
                    denthdt.append(dtempdt2d[i][j] / rho_vol)
        # - psi_ytchen: 'dtempdt' replaced by 'denthdt'
        rhs = dmdotdt + dlendt + denthdt
        return rhs
