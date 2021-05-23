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

        # lists for pressure, temperature, type, velocity, reynolds, prandtl and peclet
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
        self.junpumphead = reactor.control.input['junction']['pumphead']
        # user-specified junction flowrate signal
        self.junflowrate = reactor.control.input['junction']['flowrate']

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

        # create and inverse a matrix B of left-hand sides of momentum conservation equations (self.njun) and 
        # mass conservation equations differentiated w.r.t. time (self.npipe-self.npipef)
        n = self.njun + sum(self.pipennodes)
        B = [[0]*n for i in range(n)]
        for j in range(self.njun):
            B[j][j] = 1 # dmdot/dt

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

        # prepare lists to map thermal-hydraulic and fuel rods nodes
        indx_i = [x['pipeid'] for x in reactor.control.input['fuelrod']]
        indx_j = [x['pipenode'] for x in reactor.control.input['fuelrod']]
        nfuelrods = len(indx_i)
        self.map_th = []
        self.map_fr = []
        for i in range(nfuelrods):
            for j in range(len(indx_i[i])):
                self.map_th.append((indx_i[i][j], indx_j[i][j]))
                self.map_fr.append((i,j))
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
            if self.type[i] == 'na':
                for j in range(self.pipennodes[i]):
                    t = self.temp[i][j]
                    # J.K. Fink and L. Leibowitz "Thermodynamic and Transport Properties of Sodium Liquid and Vapor", ANL/RE-95/2, 1995, https://www.ne.anl.gov/eda/ANL-RE-95-2.pdf
                    dict['rhol'].append(219.0 + 275.32*(1.0 - t/2503.7) + 511.58*(1.0 - t/2503.7)**0.5)
                    dict['visl'].append(math.exp(-6.4406 - 0.3958*math.log(t) + 556.835/t)/dict['rhol'][j])
                    dict['kl'].append(124.67 - 0.11381*t + 5.5226e-5*t**2 - 1.1842e-8*t**3)
                    # Based on fit from J.K. Fink, etal."Properties for Reactor Safety Analysis", ANL-CEN-RSD-82-2, May 1982.
                    dict['cpl'].append(1646.97 - 0.831587*t + 4.31182e-04*t**2)
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
                b[j] = 0        # then multiply matrix by list: invA*mdot = b and convert to list
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
        self.re = [[self.vel[i][j]*self.dhyd[i]/self.prop[i]['visl'][j] for j in range(self.pipennodes[i])] for i in range(self.npipe)]
        # Prandtl numbers
        self.pr = [[self.prop[i]['visl'][j]*self.prop[i]['rhol'][j]*self.prop[i]['cpl'][j]/self.prop[i]['kl'][j] for j in range(self.pipennodes[i])] for i in range(self.npipe)]
        # Peclet numbers
        self.pe = [[self.re[i][j]*self.pr[i][j] for j in range(self.pipennodes[i])] for i in range(self.npipe)]

        # TIME DERIVATIVES OF MASS FLOWRATES:
        # construct right-hand side b of system invB*[dmdotdt, P] = b
        b = [0]*(self.njun + sum(self.pipennodes))
        for j in range(self.njun):
            rho_f = self.prop[self.f[j][0]]['rhol'][self.f[j][1]]
            rho_t = self.prop[self.t[j][0]]['rhol'][self.t[j][1]]
            # gravitational head
            rhogh_f = 9.81*rho_f*self.len[self.f[j][0]]*self.dir[self.f[j][0]]/self.pipennodes[self.f[j][0]]
            if self.pipetype[self.f[j][0]] != 'freelevel': rhogh_f = - abs(rhogh_f)
            rhogh_t = 9.81*rho_t*self.len[self.t[j][0]]*self.dir[self.t[j][0]]/self.pipennodes[self.t[j][0]]
            if self.pipetype[self.t[j][0]] != 'freelevel': rhogh_t = abs(rhogh_t)
            b[j] = -0.5*(rhogh_f + rhogh_t) - self.mdot[j]*100
            if self.juntype[j] == 'independent' and self.junpumphead[j] != '':
                b[j] = reactor.control.signal[self.junpumphead[j]]
        for i in range(sum(self.pipennodes)):
            if self.pipetype[self.indx[i][0]] == 'freelevel': b[self.njun+i] = 1e5
        invBb = self.invB.dot(b).tolist()

        # read from invBb: time derivatives of flowrate in independent junctions
        dmdotdt = []
        for j in range(self.njun):
            if self.juntype[j] == 'independent':
                if self.junflowrate[j] != '':
                    dmdotdt.append(0)
                else:
                    dmdotdt.append(invBb[j])
        # read from invBb: pressures in pipe nodes
        indx = 0
        for i in range(self.npipe):
            for j in range(self.pipennodes[i]):
                self.p[i][j] = invBb[self.njun+indx]
                indx += 1

        # TIME DERIVATIVES OF FLUID TEMPERATURES:
        dtempdt2d = [[0]*self.pipennodes[i] for i in range(self.npipe)]
        for j in range(self.njun):
            if self.mdot[j] > 0:
                cp_temp_mdot = self.prop[self.f[j][0]]['cpl'][self.f[j][1]] *  self.temp[self.f[j][0]][self.f[j][1]] * self.mdot[j]
            else:
                cp_temp_mdot = self.prop[self.t[j][0]]['cpl'][self.t[j][1]] *  self.temp[self.t[j][0]][self.t[j][1]] * self.mdot[j]
            dtempdt2d[self.f[j][0]][self.f[j][1]] -= cp_temp_mdot
            dtempdt2d[self.t[j][0]][self.t[j][1]] += cp_temp_mdot

        dtempdt = []
        for i in range(self.npipe):
            if self.signaltemp[i] != '':
                for j in range(self.pipennodes[i]):
                    dtempdt.append(0)
            else:
                vol = self.areaz[i] * abs(self.len[i])/self.pipennodes[i]
                for j in range(self.pipennodes[i]):
                    # check if there is a fuel rod cooled by the node
                    if ('fuelrod' in reactor.solve and self.pipeid[i],j) in self.map_th:
                        indx = self.map_th.index((self.pipeid[i],j))
                        tuple_fr = self.map_fr[indx]
                        tclad = reactor.solid.fuelrod[tuple_fr[0]].clad[tuple_fr[1]].temp[-1]
                        mltpl = reactor.solid.fuelrod[tuple_fr[0]].clad[tuple_fr[1]].mltpl
                        ro = reactor.solid.fuelrod[tuple_fr[0]].clad[tuple_fr[1]].r[-1]
                        area_ht = 2 * math.pi * ro * mltpl
                        dtempdt2d[i][j] += 1e3*(tclad - self.temp[i][j]) * area_ht
                    rho_cp_vol = self.prop[i]['rhol'][j] * self.prop[i]['cpl'][j] * vol
                    dtempdt.append(dtempdt2d[i][j] / rho_cp_vol)

        rhs = dmdotdt + dtempdt
        return rhs
