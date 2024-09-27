import math
import sys

#--------------------------------------------------------------------------------------------------
class Fuel:

    #----------------------------------------------------------------------------------------------
    # constructor: self is a 'fuel' object created in B1B, 
    # indx is the axial index of this object in the fuel rod with index indxfuelrod
    def __init__(self, indx, indxfuelrod, dz, reactor):

        # INITIALIZATION
        # dictionary of the fuel rod to which the fuel belongs
        dictfuelrod = reactor.control.input['fuelrod'][indxfuelrod]
        # current fuel id
        fuelid = dictfuelrod['fuelid'][indx]
        # radial power peaking factor of fuel rod
        self.kr = dictfuelrod['kr'][indx]
        # axial power peaking factor of fuel
        self.kz = dictfuelrod['kz'][indx]

        # list of fuel dictionaries specified in input
        list = reactor.control.input['fuel']
        # index of the current fuel in the list of fuel dictionaries
        i = [x['id'] for x in list].index(fuelid)
        
        if 'rx' in list[i] : # psi_ytchen: multilayer solid model
            # fuel inner radius
            self.ri = list[i]['rx'][0]
            # fuel outer radius
            self.ro = list[i]['rx'][-1]
            # total number of fuel radial nodes
            self.nr = sum(list[i]['nrx'])
            # list of fuel material type
            self.type = []
            # list of fuel material input variables
            self.pu   = []
            self.b    = []
            self.x    = []
            self.por  = []
            self.temp = []
            # list of fuel mesh geometry variables
            
            nr1 = list[i]['nrx'][0]
            nr2 = list[i]['nrx'][1]
            nr3 = list[i]['nrx'][2]
            
            dr1 = (list[i]['rx'][1]-list[i]['rx'][0])/(nr1-0.5)
            dr2 = (list[i]['rx'][2]-list[i]['rx'][1])/(nr2)
            dr3 = (list[i]['rx'][3]-list[i]['rx'][2])/(nr3-0.5)
            
            self.dr = [dr1]*nr1 + [dr2]*nr2 + [dr3]*nr3
            self.r  = [self.ri + j*dr1 for j in range(nr1)] + [list[i]['rx'][1] + (j+0.5)*dr2 for j in range(nr2)] + [list[i]['rx'][2] + (j+0.5)*dr3 for j in range(nr3)]
            self.rb = [self.r[j]+self.dr[j]/2 for j in range(self.nr-1)] + [self.ro] # right boundary position of mesh cell
            self.vol= [self.rb[0]**2 - self.r[0]**2] + [self.rb[j]**2 - self.rb[j-1]**2 for j in range(1, self.nr)]
            self.vol= [self.vol[j]*math.pi for j in range(self.nr)]
            self.frx= [list[i]['frx'][0]]*nr1 + [list[i]['frx'][1]]*nr2 + [list[i]['frx'][2]]*nr3 # fraction of heat deposit
            self.vf = [self.vol[j]*self.frx[j] for j in range(self.nr)] 
            
            for j in range(self.nr):
                if j < list[i]['nrx'][0]:
                    # the first layer of material
                    matid = list[i]['matid'][0]
                elif j < list[i]['nrx'][0] + list[i]['nrx'][1]:
                    # the second layer of material
                    matid = list[i]['matid'][1]
                else:
                    # the third layer of material
                    matid = list[i]['matid'][2]
                # find the fuel material id in the list of materials
                try:
                    ifuel = [x['id'] for x in reactor.control.input['mat']].index(matid)
                except:
                    print('****ERROR: fuel material id ' + matid + ' is not specified in the \'mat\' card of input.')
                    sys.exit()
                mat = reactor.control.input['mat'][ifuel]
                self.type.append(mat['type'])
                
                if 'pu' in mat:
                    self.pu.append(mat['pu'])
                else:
                    self.pu.append(0.0)
                # list of fuel burnup in fuel radial nodes
                if 'b' in mat:
                    self.b.append(mat['b'])
                else:
                    self.b.append(0.0)
                # list of deviation from stoechiometry in fuel radial nodes
                if 'x' in mat:
                    self.x.append(mat['x'])
                else:
                    self.x.append(0.0)
                # list of porosity in fuel radial nodes
                if 'por' in mat:
                    self.por.append(mat['por'])
                else:
                    self.por.append(0.0)
                # list of initial temperatures in fuel radial nodes
                self.temp.append(mat['temp0'])
  
        else: # psi_ytchen: original fuel model
            # fuel inner radius
            self.ri = list[i]['ri']
            # fuel outer radius
            self.ro = list[i]['ro']
            # number of fuel radial nodes
            self.nr = list[i]['nr']
            
            # fuel material id
            matid = list[i]['matid']
            # find the fuel material id in the list of materials
            try:
                ifuel = [x['id'] for x in reactor.control.input['mat']].index(matid)
            except:
                print('****ERROR: fuel material id ' + matid + ' is not specified in the \'mat\' card of input.')
                sys.exit()
            # dictionary of material properties of the current fuel
            mat = reactor.control.input['mat'][ifuel]
            # material type of fuel
            self.type = mat['type']
            # list of Pu content in fuel radial nodes
            if 'pu' in mat:
                self.pu = [mat['pu']]*self.nr
            else:
                self.pu = [0]*self.nr
            # list of fuel burnup in fuel radial nodes
            if 'b' in mat:
                self.b = [mat['b']]*self.nr
            else:
                self.b = [0]*self.nr
            # list of deviation from stoechiometry in fuel radial nodes
            if 'x' in mat:
                self.x = [mat['x']]*self.nr
            else:
                self.x = [0]*self.nr
            # list of porosity in fuel radial nodes
            if 'por' in mat:
                self.por = [mat['por']]*self.nr
            else:
                self.por = [0]*self.nr
            # list of initial temperatures in fuel radial nodes
            self.temp = [mat['temp0']]*self.nr
            
            # mesh grid step
            self.dr = (self.ro - self.ri)/(self.nr-1)
            # list of node radii (size = nr)
            self.r = [self.ri + i*self.dr for i in range(self.nr)]
            # list of node boundary radii (size = nr-1)
            self.rb = [self.r[i]+self.dr/2 for i in range(self.nr-1)]
            # list of node volume (size = nr)
            self.vol = [self.rb[0]**2 - self.r[0]**2] + [self.rb[i]**2 - self.rb[i-1]**2 for i in range(1, self.nr-1)] + [self.r[self.nr-1]**2 - self.rb[self.nr-2]**2]       
            self.vol = [self.vol[i]*math.pi for i in range(self.nr)]

    #----------------------------------------------------------------------------------------------
    # create right-hand side list: self is a 'fuel' object created in B1B
    # indx is the axial index of this object in the fuel rod with index indxfuelrod
    def calculate_rhs(self, indx, indxfuelrod, reactor, t):

        # construct right-hand side list
        rhs = []

        # FUEL PROPERTIES:
        self.prop = {'rho':[], 'cp':[], 'k':[]}
        for j in range(self.nr):
            if isinstance(self.type, list): # psi_ytchen: multilayer solid model
                # call material property function
                pro = reactor.data.matpro( {'type':self.type[j], 't':self.temp[j], 'b':self.b[j], 'por':self.por[j], 'pu':self.pu[j], 'x':self.x[j]} )
            else:
                # call material property function
                pro = reactor.data.matpro( {'type':self.type, 't':self.temp[j], 'b':self.b[j], 'por':self.por[j], 'pu':self.pu[j], 'x':self.x[j]} )
 
            # density (kg/m3)
            self.prop['rho'].append(pro['rho'])
            # specific heat (J/kg-K)
            self.prop['cp'].append(pro['cp'])
            # thermal conductivity (W/m-K)
            self.prop['k'].append(pro['k'])

        # TIME DERIVATIVE OF FUEL TEMPERATURE:
        # inner gas object
        innergas = reactor.solid.fuelrod[indxfuelrod].innergas
        # gap conductance list
        hgap = innergas.calculate_hgap(indxfuelrod, reactor, t)
        # clad object
        clad = reactor.solid.fuelrod[indxfuelrod].clad[indx]

        # fuel thermal conductivity between nodes
        if isinstance(self.type, list): # psi_ytchen: multilayer solid model
            kc = self.prop['k']
            # heat resistance at each cell boundary
            hr = [(self.rb[i]-self.r[i])/kc[i]+(self.r[i+1]-self.rb[i])/kc[i+1] for i in range(self.nr-1)]
            Q = [0] + [2*math.pi*self.rb[i]*(self.temp[i] - self.temp[i+1])/hr[i] for i in range(self.nr-1)]
        else:
            kb = [0.5*(self.prop['k'][i] + self.prop['k'][i+1]) for i in range(self.nr-1)]
            # heat flux (W/m**2) times heat transfer area per unit height at node boundaries: 2*rb * kb * dT/dr (size = nr-1)
            Q = [0] + [2*math.pi*self.rb[i]*kb[i]*(self.temp[i] - self.temp[i+1])/self.dr for i in range(self.nr-1)]
        # add heat flux (W/m**2) times heat transfer area per unit height from fuel to clad 
        Q += [math.pi*(self.ro + clad.ri) * hgap[indx] * (self.temp[self.nr-1] - clad.temp[0])]
        # power density
        qv = reactor.core.qv_average * self.kr * self.kz  # psi_ytchen: reactor.core.qv_average store the average power density corresponding to power0
        # + psi_ytchen: impose the power-time signal to power density qv
        if 'powertable' in reactor.solve:
            Prate = reactor.control.signal['POWTAB']
            qv = qv*max(Prate, 0.00)
        # - psi_ytchen: impose the power-time signal to power density qv
        rhocpv = [self.prop['rho'][i]*self.prop['cp'][i]*self.vol[i] for i in range(self.nr)]
        if isinstance(self.type, list): # psi_ytchen: multilayer solid model
            dTdt = [(Q[i] - Q[i+1] + qv*self.vol[i]*self.frx[i])/rhocpv[i] for i in range(self.nr)]
        else:
            dTdt = [(Q[i] - Q[i+1] + qv*self.vol[i])/rhocpv[i] for i in range(self.nr)]
        rhs += dTdt

        return rhs
