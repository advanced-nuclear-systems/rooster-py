#--------------------------------------------------------------------------------------------------
class HeatStructure:

    # constructor: self is a 'htstr' object created in B1
    def __init__(self, indx, reactor):
        # INITIALIZATION
        # heat structure id
        self.id = reactor.control.input['htstr'][indx]['id']
        # heat structure matid
        self.matid = reactor.control.input['htstr'][indx]['matid']
        # heat structure inner radius
        self.ri = reactor.control.input['htstr'][indx]['ri']
        # heat structure outer radius
        self.ro = reactor.control.input['htstr'][indx]['ro']
        # number of heat structure radial nodes
        self.nr = reactor.control.input['htstr'][indx]['nr']
        # heat structure left boundary condition
        self.bcleft = reactor.control.input['htstr'][indx]['bcleft']
        # heat structure right boundary condition
        self.bcright = reactor.control.input['htstr'][indx]['bcright']
        # heat structure multiplicity
        self.mltpl = reactor.control.input['htstr'][indx]['mltpl']

        # find the heat structure material id in the list of materials
        try:
            ihtstr = [x['id'] for x in reactor.control.input['mat']].index(self.matid)
        except:
            print('****ERROR: heat structure material id ' + self.matid + ' is not specified in the \'mat\' card of input.')
            sys.exit()
        # dictionary of material properties of the current heat structure
        mat = reactor.control.input['mat'][ihtstr]
        # material type of heat structure
        self.type = mat['type']

        # find the heat structure left thermal boundary condition id in the list of thermal boundary conditions
        try:
            ihtstr = [x['id'] for x in reactor.control.input['thermbc']].index(self.bcleft)
        except:
            print('****ERROR: heat structure left thermal boundary condition id ' + self.bcleft + ' is not specified in the \'thermbc\' card of input.')
            sys.exit()
        # dictionary of left thermal boundary condition of the current heat structure
        self.bcleft = reactor.control.input['thermbc'][ihtstr]

        # find the heat structure right thermal boundary condition id in the list of thermal boundary conditions
        try:
            ihtstr = [x['id'] for x in reactor.control.input['thermbc']].index(self.bcright)
        except:
            print('****ERROR: heat structure right thermal boundary condition id ' + self.bcright + ' is not specified in the \'thermbc\' card of input.')
            sys.exit()
        # dictionary of right thermal boundary condition of the current heat structure
        self.bcright = reactor.control.input['thermbc'][ihtstr]        

        # list of initial temperatures in heat structure radial nodes
        self.temp = [mat['temp0']]*self.nr

        # mesh grid step
        self.dr = (self.ro - self.ri)/(self.nr-1)
        # list of node radii (size = nr)
        self.r = [self.ri + i*self.dr for i in range(self.nr)]
        # list of node boundary radii (size = nr-1)
        self.rb = [self.r[i]+self.dr/2 for i in range(self.nr-1)]
        # list of node volume per unit height (size = nr)
        self.vol = [self.rb[0]**2 - self.r[0]**2] + [self.rb[i]**2 - self.rb[i-1]**2 for i in range(1, self.nr-1)] + [self.r[self.nr-1]**2 - self.rb[self.nr-2]**2]
        
        # + psi_ytchen: properties added by ytchen
        self.qfluxleft  = 0.0 # psi_ytchen: heat flux at bcleft to fluid
        self.qfluxright = 0.0 # psi_ytchen: heat flux at bcright to fluid
        # - psi_ytchen: properties added by ytchen
    #----------------------------------------------------------------------------------------------
    # calculate right-hand side list: self is a 'htstr' object created in B1
    def compose_rhs(self, reactor, t):

        # HEAT STRUCTURE PROPERTIES:
        self.prop = {'rho':[], 'cp':[], 'k':[]}
        for j in range(self.nr):
            # call material property function
            pro = reactor.data.matpro( {'type':self.type, 't':self.temp[j]} )
            # density (kg/m3)
            self.prop['rho'].append(pro['rho'])
            # specific heat (J/kg-K)
            self.prop['cp'].append(pro['cp'])
            # thermal conductivity (W/m-K)
            self.prop['k'].append(pro['k'])

        # heat structure thermal conductivity between nodes
        kb = [0.5*(self.prop['k'][i] + self.prop['k'][i+1]) for i in range(self.nr-1)]

        # left boundary condition
        if self.bcleft['type'] == 0:
            Qleft = 2*self.r[0]*self.bcleft['qf']# * self.mltpl  # psi_ytchen: the self.mltpl here may be needless
        elif self.bcleft['type'] == 1:
            Qleft = 2*self.r[0]*self.bcleft['alfa']*(self.bcleft['temp'] - self.temp[0])# * self.mltpl  # psi_ytchen: the self.mltpl here may be needless
        elif self.bcleft['type'] == 2:
            # heat flux (W/m**2) times heat transfer area per unit height divided by pi from clad to coolant
            Qleft = -2.0*self.r[0]*self.qfluxleft # * self.mltpl  # psi_ytchen: the self.mltpl here may be needless
        else:
            Qleft = 0

        # right boundary condition
        if self.bcright['type'] == 0:
            Qright = -2*self.r[self.nr-1]*self.bcright['qf']# * self.mltpl  # psi_ytchen: the self.mltpl here may be needless
        elif self.bcright['type'] == 1:
            Qright = -2*self.r[self.nr-1]*self.bcright['alfa']*(self.bcright['temp'] - self.temp[self.nr-1])# * self.mltpl  # psi_ytchen: the self.mltpl here may be needless
        elif self.bcright['type'] == 2:
            # heat flux (W/m**2) times heat transfer area per unit height divided by pi from clad to coolant
            Qright = 2.0*self.r[self.nr-1]*self.qfluxright # * self.mltpl  # psi_ytchen: the self.mltpl here may be needless
        else:
            Qright = 0

        # list of heat flux (W/m**2) times heat transfer area per unit height at node boundaries: 2*rb * kb * dT/dr (size = nr-1)
        Q = [Qleft] + [2*self.rb[i]*kb[i]*(self.temp[i] - self.temp[i+1])/self.dr for i in range(self.nr-1)] + [Qright]
        rhocpv = [self.prop['rho'][i]*self.prop['cp'][i]*self.vol[i] for i in range(self.nr)]
        dTdt = [(Q[i] - Q[i+1])/rhocpv[i] for i in range(self.nr)]
        rhs = dTdt
        return rhs
