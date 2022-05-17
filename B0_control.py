from scipy.interpolate import interp1d
from sympy import *

import datetime
import json
import os
import shutil
import sys

#--------------------------------------------------------------------------------------------------
class Control:

    # constructor: self is a 'control' object created in B
    def __init__(self, reactor):
        self.input = self.construct_input()
        self.signal = {}
        #self.evaluate_signals(reactor, self.input['t0'])

    #------------------------------------------------------------------------------------------
    def evaluate_signals(self, reactor, t):

        # evaluate signals
        for s in self.input['signal']:

            # boolean
            if s['type'] == 'boolean':
                self.signal[s['id']] = False
                if isinstance(s['value'][0], float) or isinstance(s['value'][0], int):
                    sigid1 = s['value'][0]
                else:
                    sigid1 = self.signal[s['value'][0]]
                if isinstance(s['value'][2], float) or isinstance(s['value'][2], int):
                    sigid2 = s['value'][2]
                else:
                    sigid2 = self.signal[s['value'][2]]
                operator = s['value'][1]
                if operator == 'eq' and sigid1 == sigid2: self.signal[s['id']] = True
                if operator == 'ne' and sigid1 != sigid2: self.signal[s['id']] = True
                if operator == 'gt' and sigid1 > sigid2: self.signal[s['id']] = True
                if operator == 'ge' and sigid1 >= sigid2: self.signal[s['id']] = True
                if operator == 'lt' and sigid1 < sigid2: self.signal[s['id']] = True
                if operator == 'le' and sigid1 <= sigid2: self.signal[s['id']] = True

            # constant: constant value
            elif s['type'] == 'constant':
                self.signal[s['id']] = s['value'][0]

            # dens: pipe density
            elif s['type'] == 'dens':
                id = s['value'][0]
                if 'fluid' in reactor.solve and id in reactor.fluid.pipeid:
                    indx = reactor.fluid.pipeid.index(id)
                    if len(s['value']) == 1:
                        # average density
                        davg = 0.0
                        for i in range(reactor.fluid.pipennodes[indx]):
                            # call material property function
                            pro = reactor.data.matpro( {'type':reactor.fluid.type[indx], 't':reactor.fluid.temp[indx][i]} )
                            davg += pro['rhol']
                        davg /= reactor.fluid.pipennodes[indx]
                        self.signal[s['id']] = davg
                    else:
                        # node density
                        if s['value'][1] > reactor.fluid.pipennodes[indx]:
                            print('****ERROR: \'signal\' card ' + s['id'] + ' refers to node (' + str(int(s['value'][1])) + ') that does not exist in pipe ' + id)
                            sys.exit()
                        # call material property function
                        pro = reactor.data.matpro( {'type':reactor.fluid.type[indx], 't':reactor.fluid.temp[indx][int(s['value'][1])-1]} )
                        self.signal[s['id']] = pro['rhol']
                else:
                    print('****ERROR: \'signal\' card ' + s['id'] + ' refers to pipe that does not exist or there is no \'solve fluid\' card.')
                    sys.exit()

            # formula: symbolic evaluations
            elif s['type'] == 'formula':
                # merge card values
                value = ''.join([str(x) for x in s['value']])
                # only for signals requiring symbolic evaluations
                if any([char in value for char in ['+', '-', '*', '/']]):
                    try:
                        self.signal[s['id']] = sympify(value)
                    except:
                        print('****ERROR: \'signal\' card ' + s['id'] + ' contains a syntax error.')
                        sys.exit()
                    for id in list(self.signal.keys()):
                        if id in value:
                            self.signal[s['id']] = self.signal[s['id']].subs(sympify(id),self.signal[id])
                    try:
                        self.signal[s['id']] = float(self.signal[s['id']])
                    except:
                        print('****ERROR: \'signal\' card ' + s['id'] + ' most likely contains a not-defined signal.')
                        sys.exit()

            # if: conditional if
            elif s['type'] == 'if':
                if self.signal[s['value'][0]]:
                    self.signal[s['id']] = self.signal[s['value'][1]]

            # lookup: table
            elif s['type'] == 'lookup':
                x = []
                y = []
                for i in range(len(s['value'][1:])):
                    if i % 2 == 0: 
                        y.append(s['value'][i+1])
                    else:
                        x.append(s['value'][i+1])
                # scipy function
                f = interp1d(x, y)
                id = s['value'][0]
                xnew = max(min(self.signal[id],x[-1]),x[0])
                self.signal[s['id']] = f(xnew)

            # tclad: clad temperature
            elif s['type'] == 'tclad':
                id = s['value'][0]
                if 'fuelrod' in reactor.solve and id in [x.id for x in reactor.solid.fuelrod]:
                    indx = [x.id for x in reactor.solid.fuelrod].index(id)
                    if len(s['value']) == 1:
                        # r-z-average clad temperature and volume
                        tavg, vol = 0.0, 0.0
                        for i in range(reactor.solid.fuelrod[indx].nz):
                            for j in range(reactor.solid.fuelrod[indx].clad[i].nr):
                                tavg += reactor.solid.fuelrod[indx].clad[i].temp[j] * reactor.solid.fuelrod[indx].clad[i].vol[j]
                                vol += reactor.solid.fuelrod[indx].clad[i].vol[j]
                        tavg /= vol
                        self.signal[s['id']] = tavg
                    elif len(s['value']) == 2:
                        if s['value'][1] > reactor.solid.fuelrod[indx].nz:
                            print('****ERROR: \'signal\' card ' + s['id'] + ' refers to axial layer (' + str(int(s['value'][1])) + ') that does not exist in fuelrod ' + id)
                            sys.exit()
                        i = int(s['value'][1])
                        # r-average temperature and volume
                        tavg, vol = 0.0, 0.0
                        for j in range(reactor.solid.fuelrod[indx].clad[i].nr):
                            tavg += reactor.solid.fuelrod[indx].clad[i].temp[j] * reactor.solid.fuelrod[indx].clad[i].vol[j]
                            vol += reactor.solid.fuelrod[indx].clad[i].vol[j]
                        tavg /= vol
                        self.signal[s['id']] = tavg
                    else:
                        if s['value'][1] > reactor.solid.fuelrod[indx].nz:
                            print('****ERROR: \'signal\' card ' + s['id'] + ' refers to axial layer (' + str(int(s['value'][1])) + ') that does not exist in fuelrod ' + id)
                            sys.exit()
                        i = int(s['value'][1]-1)
                        if s['value'][2] > reactor.solid.fuelrod[indx].clad[i].nr:
                            print('****ERROR: \'signal\' card ' + s['id'] + ' refers to radial (' + str(int(s['value'][2])) + ') that does not exist in fuel of fuelrod ' + id)
                            sys.exit()
                        j = int(s['value'][2]-1)
                        # node temperature
                        self.signal[s['id']] = reactor.solid.fuelrod[indx].clad[int(s['value'][1])-1].temp[j]
                else:
                    print('****ERROR: \'signal\' card ' + s['id'] + ' refers to fuel rod (' + id + ') that does not exist or there is no \'solve fuelrod\' card.')
                    sys.exit()

            # temp: htstr or pipe temperature
            elif s['type'] == 'temp':
                id = s['value'][0]
                if 'fluid' in reactor.solve and id in reactor.fluid.pipeid:
                    indx = reactor.fluid.pipeid.index(id)
                    if len(s['value']) == 1:
                        # average temperature
                        tavg = 0.0
                        for i in range(reactor.fluid.pipennodes[indx]):
                            tavg += reactor.fluid.temp[indx][i]
                        tavg /= reactor.fluid.pipennodes[indx]
                        self.signal[s['id']] = tavg
                    else:
                        # node temperature
                        if s['value'][1] > reactor.fluid.pipennodes[indx]:
                            print('****ERROR: \'signal\' card ' + s['id'] + ' refers to node (' + str(int(s['value'][1])) + ') that does not exist in pipe ' + id)
                            sys.exit()
                        self.signal[s['id']] = reactor.fluid.temp[indx][int(s['value'][1])-1]
                elif 'htstr' in reactor.solve and id in [x.id for x in reactor.solid.htstr]:
                    indx = [x.id for x in reactor.solid.htstr].index(id)
                    if len(s['value']) == 1:
                        # average temperature
                        tavg = 0.0
                        for i in range(reactor.solid.htstr[indx].nr):
                            tavg += reactor.solid.htstr[indx].temp[i] * reactor.solid.htstr[indx].vol[i]
                        tavg /= sum(reactor.solid.htstr[indx].vol)
                        self.signal[s['id']] = tavg
                    else:
                        # node temperature
                        if s['value'][1] > reactor.solid.htstr[indx].nr:
                            print('****ERROR: \'signal\' card ' + s['id'] + ' refers to radial node (' + str(int(s['value'][1])) + ') that does not exist in htstr ' + id)
                            sys.exit()
                        self.signal[s['id']] = reactor.solid.htstr[indx].temp[int(s['value'][1])-1]
                else:
                    print('****ERROR: \'signal\' card ' + s['id'] + ' refers to heat structure or pipe that does not exist or there is no \'solve fluid\' card or \'solve htstr\' card.')
                    sys.exit()

            # tfuel: fuel temperature
            elif s['type'] == 'tfuel':
                id = s['value'][0]
                if 'fuelrod' in reactor.solve and id in [x.id for x in reactor.solid.fuelrod]:
                    indx = [x.id for x in reactor.solid.fuelrod].index(id)
                    if len(s['value']) == 1:
                        # r-z-average fuel temperature and volume
                        tavg, vol = 0.0, 0.0
                        for i in range(reactor.solid.fuelrod[indx].nz):
                            for j in range(reactor.solid.fuelrod[indx].fuel[i].nr):
                                tavg += reactor.solid.fuelrod[indx].fuel[i].temp[j] * reactor.solid.fuelrod[indx].fuel[i].vol[j]
                                vol += reactor.solid.fuelrod[indx].fuel[i].vol[j]
                        tavg /= vol
                        self.signal[s['id']] = tavg
                    elif len(s['value']) == 2:
                        if s['value'][1] > reactor.solid.fuelrod[indx].nz:
                            print('****ERROR: \'signal\' card ' + s['id'] + ' refers to axial layer (' + str(int(s['value'][1])) + ') that does not exist in fuelrod ' + id)
                            sys.exit()
                        i = int(s['value'][1]-1)
                        # r-average temperature and volume
                        tavg, vol = 0.0, 0.0
                        for j in range(reactor.solid.fuelrod[indx].fuel[i].nr):
                            tavg += reactor.solid.fuelrod[indx].fuel[i].temp[j] * reactor.solid.fuelrod[indx].fuel[i].vol[j]
                            vol += reactor.solid.fuelrod[indx].fuel[i].vol[j]
                        tavg /= vol
                        self.signal[s['id']] = tavg
                    else:
                        if s['value'][1] > reactor.solid.fuelrod[indx].nz:
                            print('****ERROR: \'signal\' card ' + s['id'] + ' refers to axial layer (' + str(int(s['value'][1])) + ') that does not exist in fuelrod ' + id)
                            sys.exit()
                        i = int(s['value'][1]-1)
                        if s['value'][2] > reactor.solid.fuelrod[indx].fuel[i].nr:
                            print('****ERROR: \'signal\' card ' + s['id'] + ' refers to radial (' + str(int(s['value'][2])) + ') that does not exist in fuel of fuelrod ' + id)
                            sys.exit()
                        j = int(s['value'][2]-1)
                        # node temperature
                        self.signal[s['id']] = reactor.solid.fuelrod[indx].fuel[int(s['value'][1])-1].temp[j]
                else:
                    print('****ERROR: \'signal\' card ' + s['id'] + ' refers to fuel rod (' + id + ') that does not exist or there is no \'solve fuelrod\' card.')
                    sys.exit()

            # time: time
            elif s['type'] == 'time':
                self.signal[s['id']] = t

        #evaluate boolean signals
        for x in self.input['boolean']:
            self.signal[x['sigid']] = False
            if isinstance(x['sigid1'], float) or isinstance(x['sigid1'], int):
                sigid1 = x['sigid1']
            else:
                sigid1 = self.signal[x['sigid1']]
            if isinstance(x['sigid2'], float) or isinstance(x['sigid1'], int):
                sigid2 = x['sigid2']
            else:
                sigid2 = self.signal[x['sigid1']]
            operator = x['operator']
            if operator == 'eq' and sigid1 == sigid2: self.signal[x['sigid']] = True
            if operator == 'ne' and sigid1 != sigid2: self.signal[x['sigid']] = True
            if operator == 'gt' and sigid1 > sigid2: self.signal[x['sigid']] = True
            if operator == 'ge' and sigid1 >= sigid2: self.signal[x['sigid']] = True
            if operator == 'lt' and sigid1 < sigid2: self.signal[x['sigid']] = True
            if operator == 'le' and sigid1 <= sigid2: self.signal[x['sigid']] = True

        # signal-dependent junction: impose flowrate
        if 'fluid' in reactor.solve:
            k = 0
            for j in range(reactor.fluid.njun):
                if reactor.fluid.juntype[j] == 'independent':
                    f = reactor.fluid.f[j][0]
                    t = reactor.fluid.t[j][0]
                    # tuple of from-to pipe id's
                    f_t = (reactor.fluid.pipeid[f],reactor.fluid.pipeid[t])
                    try:
                        # check if the current junction j is present in junflowrate list
                        indx = reactor.fluid.junflowrate['jun'].index(f_t)
                        # impose flowrate from the look-up table
                        reactor.fluid.mdoti[k] = self.signal[reactor.fluid.junflowrate['flowrate'][indx]]
                    except ValueError:
                        pass
                    k += 1

        # signal-dependent pipe: impose temperature
        if 'fluid' in reactor.solve:
            for i in range(reactor.fluid.npipe):
                if reactor.fluid.pipetype[i] == 'normal' and reactor.fluid.signaltemp[i] != '':
                    # impose temperature from the look-up table
                    reactor.fluid.temp[i] = [self.signal[reactor.fluid.signaltemp[i]]] * reactor.fluid.pipennodes[i]
    #----------------------------------------------------------------------------------------------
    def construct_input(self):
        #create dictionary inp where all input data will be stored
        inp = {}
        inp['boolean'] = []
        inp['clad'] = []
        inp['coregeom'] = {'geometry':'', 'pitch':0, 'botBC':'', 'topBC':''}
        inp['coremap'] = []
        inp['freeze'] = {'sigid':[],'sigidbool':[]}
        inp['fuel'] = []
        inp['fuelrod'] = []
        inp['innergas'] = []
        inp['junction'] = {'from':[], 'to':[], 'type':[]}
        inp['junpumphead'] = {'jun':[], 'pumphead':[]}
        inp['junflowrate'] = {'jun':[], 'flowrate':[]}
        inp['junkfac'] = {'jun':[], 'kfac':[]}
        inp['lookup'] = []
        inp['mat'] = []
        inp['mix'] = []
        inp['p2d'] = []
        inp['pipe'] = []
        inp['signal'] = []
        inp['signalid'] = []
        inp['solve'] = []
        inp['stack'] = []
        inp['htstr'] = []
        inp['t0'] = 0
        inp['tend'] = []
        inp['tol'] = (1.e-6,1e-6)
        inp['thermbc'] = []
    
        #read input file as a whole
        f = open('input', 'r')
        s0 = f.read()
        f.close()
    
        #merge &-ending "line" with the next one
        s = ''
        take = True
        for c in s0:
            if c == '&' : take = False
            if take : s += c
            if c == '\n' : take = True
    
        #split in lines
        lines = s.strip().split('\n')
    
        #remove comment-lines (#)
        lines = [x for x in lines if not x.startswith('#')]
        #remove comments inside lines (#)
        for i in range(len(lines)):
            if '#' in lines[i]:
                lines[i] = lines[i].split('#')[0]
    
        def convert_to_float(w): 
            try:
                w = float(w)
            except:
                pass
            return w

        for line in lines:
                
            word = line.split()
            word = list(map(convert_to_float, word))
            if len(word) > 0:

                key = word[0].lower()
                #--------------------------------------------------------------------------------------
                # just placeholder
                if key == '':
                    pass
                #--------------------------------------------------------------------------------------
                # effective delayed neutron fractions
                elif key == 'betaeff':
                    inp['betaeff'] = word[1:]
                #--------------------------------------------------------------------------------------
                # boolean signal
                elif key == 'boolean':
                    if len(word)-1 < 4:
                        print('****ERROR: boolean card should have four values after the keyword: boolean signal id, first signal id, boolean operator(eq, ne, gt, ge, lt, le), second signal id.')
                        sys.exit()
                    if word[3] not in ['eq', 'ne', 'gt', 'ge', 'lt', 'le']:
                        print('****ERROR: wrong operator in boolean card', word[3], '. Possible options are: eq, ne, gt, ge, lt, le.')
                        sys.exit()
                    inp['boolean'].append({'sigid':word[1],'sigid1':word[2],'operator':word[3],'sigid2':word[4]})
                #--------------------------------------------------------------------------------------
                # cladding
                elif key == 'clad':
                    inp['clad'].append( {'id':word[1], 'matid':word[2], 'ri':word[3], 'ro':word[4], 'nr':int(word[5])} )
                #--------------------------------------------------------------------------------------
                # core geometry
                elif key == 'coregeom':
                    if len(word)-1 < 4:
                        print('****ERROR: \'coregeom\' card should have four values after the keyword: geometry flag (hex01, hex06, hex24, square), pitch (distance between node centres), bottom boundary conditions (0: vacuum, -1: reflective), top boundary conditions (0: vacuum, -1: reflective).')
                        sys.exit()
                    list_of_geometries = ['square','hex01', 'hex06', 'hex24']
                    if not word[1] in list_of_geometries:
                        print('****ERROR: geometry flag of \'coregeom\' card (word 2) is wrong: ', word[1], '\nCorrect values are: ')
                        for v in list_of_geometries:
                            print(v)
                        sys.exit()
                    if not isinstance(word[2],int) and not isinstance(word[2],float):
                        print('****ERROR: node pitch (m) of \'coregeom\' card (word 3) is not numeric: ', word[2])
                        sys.exit()
                    if word[3] != 0 and word[3] != 1:
                        print('****ERROR: bottom boundary condition flag of \'coregeom\' card (word 4) is wrong: ', word[3], '\nCorrect values are:\n0 (vacuum)\n1 (reflective)')
                        sys.exit()
                    if word[4] != 0 and word[4] != 1:
                        print('****ERROR: top boundary condition flag of \'coregeom\' card (word 5) is wrong: ', word[4], '\nCorrect values are:\n0 (vacuum)\n1 (reflective)')
                        sys.exit()
                    inp['coregeom'] = {'geom':word[1], 'pitch':word[2], 'botBC':int(word[3]), 'topBC':int(word[4])}
                #--------------------------------------------------------------------------------------
                # core map
                elif key == 'coremap':
                    inp['coremap'].append(word[1:])
                #--------------------------------------------------------------------------------------
                # delayed neutron precursor decay time constants
                elif key == 'dnplmb':
                    inp['dnplmb'] = word[1:]
                #--------------------------------------------------------------------------------------
                # fuel grain parameters
                elif key == 'fgrain':
                    # grain diameter
                    inp['dgrain'] = word[1]
                    # number of nodes in the grain
                    inp['nrgrain'] = int(word[2])
                    # fission rate
                    inp['frate'] = int(word[3])
                #--------------------------------------------------------------------------------------
                # freeze signal
                elif key == 'freeze':
                    if len(word)-1 < 2:
                        print('****ERROR: freeze card should have two values after the keyword: existing signal id to be frozen and boolean signal id which freezes the previous signal when True.')
                        sys.exit()
                    inp['freeze']['sigid'].append(word[1])
                    inp['freeze']['sigidbool'].append(word[2])
                #--------------------------------------------------------------------------------------
                # fuel
                elif key == 'fuel':
                    inp['fuel'].append( {'id':word[1], 'matid':word[2], 'ri':float(word[3]), 'ro':float(word[4]), 'nr':int(word[5])} )
                #--------------------------------------------------------------------------------------
                # fuel rod card
                elif key == 'fuelrod':
                    id = word[1]
                    if any([id in x['id'] for x in inp['fuelrod']]):
                        for x in inp['fuelrod']:
                            if x['id'] == id:
                                x['fuelid'].append(word[2])
                                x['hgap'].append(float(word[3]))
                                x['cladid'].append(word[4])
                                x['p2d'].append(float(word[5]))
                                x['mltpl'].append(float(word[6]))
                                x['pipeid'].append(word[7])
                                x['pipenode'].append(int(word[8]))
                                x['kr'].append(float(word[9]))
                                x['kz'].append(float(word[10]))
                    else:
                        inp['fuelrod'].append({'id':id, 'fuelid':[word[2]], 'hgap':[float(word[3])], 'cladid':[word[4]], 'p2d':[float(word[5])], 'mltpl':[float(word[6])], 'pipeid':[word[7]], 'pipenode':[int(word[8])], 'kr':[float(word[9])], 'kz':[float(word[10])]})
                #--------------------------------------------------------------------------------------
                # heat structure card
                elif key == 'htstr':
                    inp['htstr'].append({'id':word[1], 'matid':word[2], 'ri':float(word[3]), 'ro':float(word[4]), 'nr':int(word[5]), 'bcleft':word[6], 'bcright':word[7], 'mltpl':word[8]})
                #--------------------------------------------------------------------------------------
                # inner gas
                elif key == 'innergas':
                    inp['innergas'].append( {'fuelrodid':word[1], 'matid':word[2], 'plenv':word[3]} )
                #--------------------------------------------------------------------------------------
                # thermal-hydraulic junction (dependent)
                elif key == 'jun':
                    inp['junction']['from'].append(word[1])
                    inp['junction']['to'].append(word[2])
                    inp['junction']['type'].append('dependent')
                #--------------------------------------------------------------------------------------
                # thermal-hydraulic junction (independent)
                elif key == 'jun-i':
                    inp['junction']['from'].append(word[1])
                    inp['junction']['to'].append(word[2])
                    inp['junction']['type'].append('independent')
                #--------------------------------------------------------------------------------------
                # thermal-hydraulic junction (independent + signal for flowrate)
                elif key == 'jun-i-f':
                    inp['junction']['from'].append(word[1])
                    inp['junction']['to'].append(word[2])
                    inp['junction']['type'].append('independent')
                    inp['junflowrate']['jun'].append((word[1],word[2]))
                    inp['junflowrate']['flowrate'].append(word[3])
                #--------------------------------------------------------------------------------------
                # thermal-hydraulic junction (independent + signal for pump head)
                elif key == 'jun-i-p':
                    inp['junction']['from'].append(word[1])
                    inp['junction']['to'].append(word[2])
                    inp['junction']['type'].append('independent')
                    inp['junpumphead']['jun'].append((word[1],word[2]))
                    inp['junpumphead']['pumphead'].append(word[3])
                #--------------------------------------------------------------------------------------
                # k-factor at thermal-hydraulic junction
                elif key == 'jun-kfac':
                    try:
                        # make list of tuples (from,to) and find index of the (word[1],word[2]) tuple 
                        j = list(zip(inp['junction']['from'],inp['junction']['to'])).index((word[1],word[2]))
                    except ValueError:
                        print('****ERROR: from and to of \'jun-kfac\' card (word 3) are not specified in any junction card (should appear before): ', word[1], word[2])
                        sys.exit()
                    inp['junkfac']['jun'].append((word[1],word[2]))
                    inp['junkfac']['kfac'].append(word[3])
                #--------------------------------------------------------------------------------------
                # lookup table
                elif key == 'lookup':
                     lookup = {}
                     lookup['x'] = word[1::2]
                     lookup['f(x)'] = word[2::2]
                     inp['lookup'].append(lookup)
                #--------------------------------------------------------------------------------------
                # material
                elif key == 'mat':
                     if word[2] == 'he':
                         inp['mat'].append( {'id':word[1], 'type':word[2], 'p0':word[3], 'temp0':word[4]} )
                     elif word[2] == 'mox':
                         inp['mat'].append( {'id':word[1], 'type':word[2], 'pu':word[3], 'b':word[4], 'x':word[5], 'por':word[6], 'temp0':word[7]} )
                     elif word[2] == 'na':
                         inp['mat'].append( {'id':word[1], 'type':word[2], 'p0':word[3], 'temp0':word[4]} )
                     elif word[2] == 'ss316':
                         inp['mat'].append( {'id':word[1], 'type':word[2], 'temp0':word[3]} )
                #--------------------------------------------------------------------------------------
                # mixture of isotopes
                elif key == 'mix':
                    if len(word)-1 < 4:
                        print('****ERROR: mix card should have four values after the keyword: mix id, isotopeid, number density and signal id for temperature.')
                        sys.exit()
                    
                    mixid = word[1]
                    if any([mixid in x['mixid'] for x in inp['mix']]):
                        for x in inp['mix']:
                            if x['mixid'] == mixid:
                                x['isoid'].append(word[2])
                                x['numdens'].append(float(word[3]))
                                x['signaltemp'].append(word[4])
                    else:
                        inp['mix'].append({'mixid':mixid, 'isoid':[word[2]], 'numdens':[float(word[3])], 'signaltemp':[word[4]]})
                #--------------------------------------------------------------------------------------
                # nuclear data directory
                elif key == 'nddir':
                    inp['nddir'] = word[1]
                #--------------------------------------------------------------------------------------
                # thermal-hydraulic pipe without free level
                elif key == 'pipe':
                    inp['pipe'].append( {'id':word[1], 'type':'normal', 'matid':word[2], 'dhyd':word[3], 'len':word[4], 'dir':word[5], 'areaz':word[6], 'nnodes':int(word[7]), 'signaltemp':''} )
                #--------------------------------------------------------------------------------------
                # thermal-hydraulic pipe with free level
                elif key == 'pipe-f':
                    inp['pipe'].append( {'id':word[1], 'type':'freelevel', 'matid':word[2], 'dhyd':word[3], 'len':word[4], 'dir':0, 'areaz':word[5], 'nnodes':1, 'signaltemp':''} )
                #--------------------------------------------------------------------------------------
                # thermal-hydraulic pipe without free level with temperature defined by signal
                elif key == 'pipe-t':
                    inp['pipe'].append( {'id':word[1], 'type':'normal', 'matid':word[2], 'dhyd':word[3], 'len':word[4], 'dir':word[5], 'areaz':word[6], 'nnodes':int(word[7]), 'signaltemp':word[8]} )
                #--------------------------------------------------------------------------------------
                # initial reactor power
                elif key == 'power0':
                    inp['power0'] = float(word[1])
                #--------------------------------------------------------------------------------------
                # signal variable
                elif key == 'signal':
                    if len(word) < 3:
                        print('****ERROR: \'signal\' card should have at least 3 words.')
                        sys.exit()
                    inp['signal'].append( {'id':word[1], 'type':word[2], 'value':word[3:]} )
                #--------------------------------------------------------------------------------------
                # models to be solved
                elif key == 'solve':
                    inp['solve'].append(word[1])
                    # verify that solve card has correct value
                    correct_values = {'fluid','fuelgrain','fuelrod','htstr','pointkinetics','spatialkinetics'}
                    value = set([word[1]])
                    diff = value.difference(correct_values)
                    if diff != set():
                        print('****ERROR: \'solve\' card contains wrong value: ', list(diff)[0], '\nCorrect values are: ')
                        sorted = list(correct_values)
                        sorted.sort()
                        for v in sorted:
                            print('solve', v)
                        sys.exit()
                    if word[1] == 'spatialkinetics':
                        # check that there are two additional values
                        if len(word)-1 < 3:
                            print('****ERROR: solve spatialkinetics card should have two value after the keyword: number of energy groups (integer) and method indicator (DIF or MC), e.g.:\nsolve spatialkinetics 25 MC')
                            sys.exit()
                        # check that the second value is integer
                        try:
                            # number of energy groups
                            inp['ng'] = int(word[2])
                        except:
                            print('****ERROR: the second value after the keyword of solve spatialkinetics card should be integer (number of energy groups), e.g.:\nsolve spatialkinetics 25')
                            sys.exit()
                        # check that the thrid value is DIF or MC
                        if word[3] != 'DIF' and word[3] != 'MC':
                            print('****ERROR: solve spatialkinetics card should have the third value of method indicator either DIF (neutron diffusion solver) or MC (Monte Carlo method)')
                            sys.exit()
                        # method indicator
                        inp['nmeth'] = word[3]
                #--------------------------------------------------------------------------------------
                # stack of mixes of isotopes
                elif key == 'stack':
                    if len(word)-1 < 4:
                        print('****ERROR: stack card should have four values after the keyword: stack id, mix id, pipe id, pipe node.')
                        sys.exit()
                    
                    stackid = word[1]
                    if any([stackid in x['stackid'] for x in inp['stack']]):
                        for x in inp['stack']:
                            if x['stackid'] == stackid:
                                x['mixid'].append(word[2])
                                x['pipeid'].append(word[3])
                                x['pipenode'].append(int(word[4]))
                    else:
                        inp['stack'].append({'stackid':stackid, 'mixid':[word[2]], 'pipeid':[word[3]], 'pipenode':[int(word[4])]})
                #--------------------------------------------------------------------------------------
                # integration starting time
                elif key == 't0':
                    inp['t0'] = word[1]
                #--------------------------------------------------------------------------------------
                # end of time interval and output time step for this interval
                elif key == 'tend':
                    inp['tend'].append(word[1])
                #--------------------------------------------------------------------------------------
                # thermal boundary conditions]
                elif key == 'thermbc':
                    if len(word)-1 < 3:
                        print('****ERROR: thermbc card should have at least three values after the keyword.')
                        sys.exit()
                    dict = {}
                    dict['id'] = word[1]
                    try:
                        dict['type'] = int(word[2])
                    except:
                        print('****ERROR: boundary condition type of thermbc card (word 3) is wrong: ', word[2], '. Correct values are: 0 (heat flux), 1 (heat exchange coefficient and temperature) or 2 (pipe id and pipenodeid).')
                        sys.exit()
                    if dict['type'] == 0:
                        if len(word)-1 < 3:
                            print('****ERROR: thermbc card with type == 0 should have three values after the keyword: id, type and qf.')
                            sys.exit()
                        dict['qf'] = word[3]
                    elif dict['type'] == 1:
                        if len(word)-1 < 4:
                            print('****ERROR: thermbc card with type == 1 should have four values after the keyword: id, type, alfa and temp.')
                            sys.exit()
                        dict['alfa'] = word[3]
                        dict['temp'] = word[4]
                    elif dict['type'] == 2:
                        if len(word)-1 < 4:
                            print('****ERROR: thermbc card with type == 2 should have four values after the keyword: id, type, pipeid and pipenode.')
                            sys.exit()
                        dict['pipeid'] = word[3]
                        dict['pipenode'] = int(word[4])
                    else:
                        print('****ERROR: boundary condition type of thermbc card (word 3) is wrong: ', word[2], '. Correct values are: 0 (heat flux), 1 (heat exchange coefficient and temperature) or 2 (pipe id and pipenodeid).')
                        sys.exit()
                    inp['thermbc'].append(dict)
                #--------------------------------------------------------------------------------------
                # prompt neutron lifetime
                elif key == 'tlife':
                    inp['tlife'] = word[1]
                #--------------------------------------------------------------------------------------
                # tolerances (relative and absolute)
                elif key == 'tol':
                    inp['tol'] = (word[1],word[2])

        # verify that tout present
        if inp['tend'] == []:
            sys.exit('****ERROR: obligatory card tend specifying time_end is absent.')
            sys.exit()
    
        # verify that there is at least one solve card
        if len(inp['solve']) == 0:
            print('****ERROR: input file should have at least one solve card.')
            sys.exit()
        if 'fuelgrain' in inp['solve'] and 'fuelrod' not in inp['solve']:
            print('****ERROR: \'solve fuelgrain\' card requires \'solve fuelrod\' card.')
            sys.exit()
    
        # make a list of all signals
        inp['signalid'] = [x['id'] for x in inp['signal']]
        # verify that lookup tables use existing signals
        for table in inp['lookup']:
            insignal = table['x'][0]
            outsignal = table['f(x)'][0]
            if insignal not in inp['signalid']:
                print('****ERROR: input signal ' + insignal + ' in lookup table ' + outsignal + ' is not defined.')
                sys.exit()
        # append output signals of lookup tables
        inp['signalid'] += [y['f(x)'][0] for y in inp['lookup']]
        # append boolean signals
        inp['signalid'] += [x['sigid'] for x in inp['boolean']]
        # verify that boolean signals use existing signals
        for x in inp['boolean']:
            if x['sigid1'] not in inp['signalid'] and not isinstance(x['sigid1'], float) and not isinstance(x['sigid1'], int):
                print('****ERROR: firts signal ' + x['sigid1'] + ' in boolean card is not defined.')
                sys.exit()
            if x['sigid2'] not in inp['signalid'] and not isinstance(x['sigid2'], float) and not isinstance(x['sigid2'], int):
                print('****ERROR: second signal ' + x['sigid2'] + ' in boolean card is not defined.')
                sys.exit()   
        # verify that mix card uses existing signals
        for s in [x['signaltemp'][j] for x in inp['mix'] for j in range(len(x['signaltemp']))]:
            if s not in inp['signalid']:
                print('****ERROR: signal for temperature ' + s + ' in mix card is not defined.')
                sys.exit()
        # verify that signals specified in the freeze card exist
        for x in inp['freeze']['sigid']:
            if x not in inp['signalid']:
                print('****ERROR: signal ' + x + ' in freeze card is not defined.')
                sys.exit()
        for x in inp['freeze']['sigidbool']:
            if x not in inp['signalid']:
                print('****ERROR: signal ' + x + ' in freeze card is not defined.')
                sys.exit()
    
        fid = open('input.json', 'w')
        fid.write(json.dumps(inp, indent=2))
        fid.close()
        return inp

    #----------------------------------------------------------------------------------------------
    def open_output_files(self, reactor):

        # prepare an output folder
        path4results = 'output'
        if os.path.isfile(path4results): os.remove(path4results)
        if not os.path.isdir(path4results): os.mkdir(path4results)
        path4results += os.sep + str(datetime.datetime.now())[0:21].replace(' ','-').replace(':','-').replace('.','-')
        if os.path.isfile(path4results): os.remove(path4results)
        if not os.path.isdir(path4results): os.mkdir(path4results)

        # copy input files to output folder
        shutil.copyfile('input', path4results + os.sep + 'input')
        shutil.copyfile('input.json', path4results + os.sep + 'input.json')
        # open files for output
        fid = []
        if 'signal' in self.input:
            fid.append(open(path4results + os.sep + 'signal.dat', 'w'))
            fid[-1].write(' ' + 'time(s)'.ljust(13) + ''.join( [(self.input['signal'][j]['id']).ljust(13) for j in range(len(self.input['signal']))] + [table['f(x)'][0].ljust(13) for table in self.input['lookup']] + [x['sigid'].ljust(13) for x in self.input['boolean']] ) + '\n')
        if 'fluid' in reactor.solve:
            fid.append(open(path4results + os.sep + 'fluid-mdot.dat', 'w'))
            fid[-1].write(' ' + 'time(s)'.ljust(13) + ''.join([(self.input['junction']['from'][j] +'-' + self.input['junction']['to'][j]).ljust(13) for j in range(reactor.fluid.njuni + reactor.fluid.njund)]) + '\n')
            for i in range(reactor.fluid.npipe):
                fid.append(open(path4results + os.sep + 'fluid-p-' + reactor.fluid.pipeid[i] + '.dat', 'w'))
                fid[-1].write(' ' + 'time(s)'.ljust(13) + ''.join([str(j).zfill(4).ljust(13) for j in range(reactor.fluid.pipennodes[i])]) + '\n')
                fid.append(open(path4results + os.sep + 'fluid-temp-' + reactor.fluid.pipeid[i] + '.dat', 'w'))
                fid[-1].write(' ' + 'time(s)'.ljust(13) + ''.join([str(j).zfill(4).ljust(13) for j in range(reactor.fluid.pipennodes[i])]) + '\n')
                fid.append(open(path4results + os.sep + 'fluid-vel-' + reactor.fluid.pipeid[i] + '.dat', 'w'))
                fid[-1].write(' ' + 'time(s)'.ljust(13) + ''.join([str(j).zfill(4).ljust(13) for j in range(reactor.fluid.pipennodes[i])]) + '\n')
                fid.append(open(path4results + os.sep + 'fluid-re-' + reactor.fluid.pipeid[i] + '.dat', 'w'))
                fid[-1].write(' ' + 'time(s)'.ljust(13) + ''.join([str(j).zfill(4).ljust(13) for j in range(reactor.fluid.pipennodes[i])]) + '\n')
                fid.append(open(path4results + os.sep + 'fluid-pr-' + reactor.fluid.pipeid[i] + '.dat', 'w'))
                fid[-1].write(' ' + 'time(s)'.ljust(13) + ''.join([str(j).zfill(4).ljust(13) for j in range(reactor.fluid.pipennodes[i])]) + '\n')
                fid.append(open(path4results + os.sep + 'fluid-pe-' + reactor.fluid.pipeid[i] + '.dat', 'w'))
                fid[-1].write(' ' + 'time(s)'.ljust(13) + ''.join([str(j).zfill(4).ljust(13) for j in range(reactor.fluid.pipennodes[i])]) + '\n')
            fid.append(open(path4results + os.sep + 'fluid-len.dat', 'w'))
            s = ''
            for i in range(reactor.fluid.npipe):
                if reactor.fluid.pipetype[i] == 'freelevel':
                    s += str(reactor.fluid.pipeid[i]).ljust(13)
            fid[-1].write(' ' + 'time(s)'.ljust(13) + s + '\n')
        if 'fuelrod' in reactor.solve:
            for i in range(reactor.solid.nfuelrods):
                fid.append(open(path4results + os.sep + 'fuelrod-hgap-' + [x['id'] for x in self.input['fuelrod']][i] + '.dat', 'w'))
                fid[-1].write(' ' + 'time(s)'.ljust(13) + ''.join([('hgap-' + str(j).zfill(3)).ljust(13) for j in range(reactor.solid.fuelrod[i].nz)]) + '\n')
                for j in range(reactor.solid.fuelrod[i].nz):
                    fid.append(open(path4results + os.sep + 'fuelrod-temp-' + [x['id'] for x in self.input['fuelrod']][i] + '-' + str(j+1).zfill(3) + '.dat', 'w'))
                    fid[-1].write(' ' + 'time(s)'.ljust(13) + ''.join([('tempf-' + str(k).zfill(3) + '(K)').ljust(13) for k in range(reactor.solid.fuelrod[i].fuel[j].nr)]) + ''.join([('tempc-' + str(k).zfill(3) + '(K)').ljust(13) for k in range(reactor.solid.fuelrod[i].clad[j].nr)]) + '\n')
                    for k in range(reactor.solid.fuelrod[i].fuel[j].nr):
                        if 'fuelgrain' in reactor.solve and i + j + k == 0: 
                            fid.append(open(path4results + os.sep + 'fuelrod-c1-' + [x['id'] for x in self.input['fuelrod']][i] + '-' + str(j).zfill(3) + '-' + str(k).zfill(3) + '.dat', 'w'))
                            fid[-1].write(' ' + 'time(s)'.ljust(13) + ''.join([('c1-' + str(l).zfill(3)).ljust(13) for l in range(reactor.solid.fuelrod[i].fuel[j].fuelgrain[k].nr)]) + '\n')
                            fid.append(open(path4results + os.sep + 'fuelrod-ri-' + [x['id'] for x in self.input['fuelrod']][i] + '-' + str(j).zfill(3) + '-' + str(k).zfill(3) + '.dat', 'w'))
                            fid[-1].write(' ' + 'time(s)'.ljust(13) + ''.join([('ri-' + str(l).zfill(3)).ljust(13) for l in range(reactor.solid.fuelrod[i].fuel[j].fuelgrain[k].NB)]) + '\n')
                            fid.append(open(path4results + os.sep + 'fuelrod-cv_irr-' + [x['id'] for x in self.input['fuelrod']][i] + '-' + str(j).zfill(3) + '-' + str(k).zfill(3) + '.dat', 'w'))
                            fid[-1].write(' ' + 'time(s)'.ljust(13) + ''.join([('cv_irr-' + str(l).zfill(3)).ljust(13) for l in range(reactor.solid.fuelrod[i].fuel[j].fuelgrain[k].NB)]) + '\n')
                            fid.append(open(path4results + os.sep + 'fuelrod-ci_irr-' + [x['id'] for x in self.input['fuelrod']][i] + '-' + str(j).zfill(3) + '-' + str(k).zfill(3) + '.dat', 'w'))
                            fid[-1].write(' ' + 'time(s)'.ljust(13) + ''.join([('ci_irr-' + str(l).zfill(3)).ljust(13) for l in range(reactor.solid.fuelrod[i].fuel[j].fuelgrain[k].NB)]) + '\n')
                            fid.append(open(path4results + os.sep + 'fuelrod-cv_p-' + [x['id'] for x in self.input['fuelrod']][i] + '-' + str(j).zfill(3) + '-' + str(k).zfill(3) + '.dat', 'w'))
                            fid[-1].write(' ' + 'time(s)'.ljust(13) + ''.join([('cv_p-' + str(l).zfill(3)).ljust(13) for l in range(reactor.solid.fuelrod[i].fuel[j].fuelgrain[k].NB)]) + '\n')
                            fid.append(open(path4results + os.sep + 'fuelrod-bi-' + [x['id'] for x in self.input['fuelrod']][i] + '-' + str(j).zfill(3) + '-' + str(k).zfill(3) + '.dat', 'w'))
                            fid[-1].write(' ' + 'time(s)'.ljust(13) + ''.join([('bi-' + str(l).zfill(3)).ljust(13) for l in range(reactor.solid.fuelrod[i].fuel[j].fuelgrain[k].NB)]) + '\n')
        if 'htstr' in reactor.solve:
            for i in range(reactor.solid.nhtstr):
                fid.append(open(path4results + os.sep + 'htstr-temp-' + [x['id'] for x in self.input['htstr']][i] + '.dat', 'w'))
                fid[-1].write(' ' + 'time(s)'.ljust(13) + ''.join([('temp-' + str(j).zfill(3)).ljust(13) for j in range(reactor.solid.htstr[i].nr)]) + '\n')
        if 'pointkinetics' in reactor.solve:
            fid.append(open(path4results + os.sep + 'core-power.dat', 'w'))
            fid[-1].write(' ' + 'time(s)'.ljust(13) + 'power(-)\n')
            fid.append(open(path4results + os.sep + 'core-cdnp.dat', 'w'))
            fid[-1].write(' ' + 'time(s)'.ljust(13) + ''.join([('cdnp-' + str(i)).ljust(13) for i in range(reactor.core.ndnp)]) + '\n')
        if 'spatialkinetics' in reactor.solve:
            for i in range(reactor.core.niso):
                fid.append(open(path4results + os.sep + 'core-iso-microxs-' + reactor.core.isoname[i] + '.dat', 'w'))
            for i in range(reactor.core.nmix):
                fid.append(open(path4results + os.sep + 'core-mix-macroxs-' + reactor.core.mix[i].mixid + '.dat', 'w'))
            fid.append(open(path4results + os.sep + 'core-k.dat', 'w'))
            fid[-1].write(' ' + 'niter'.ljust(13) + 'k'.ljust(13) + '\n')
            fid.append(open(path4results + os.sep + 'core-flux.dat', 'w'))
            fid[-1].write(' ' + 'time(s)'.ljust(13) + 'igroup'.ljust(13) + 'iz'.ljust(13) + 'ix'.ljust(13) + 'iy'.ljust(13) + 'flux'.ljust(13) + '\n')
            fid.append(open(path4results + os.sep + 'core-pow.dat', 'w'))
            fid[-1].write(' ' + 'time(s)'.ljust(13) + 'iz'.ljust(13) + 'ix'.ljust(13) + 'iy'.ljust(13) + 'pow'.ljust(13) + '\n')
            fid.append(open(path4results + os.sep + 'core-powxy.dat', 'w'))
            fid[-1].write(' ' + 'time(s)'.ljust(13) + 'ix'.ljust(13) + 'iy'.ljust(13) + 'pow'.ljust(13) + '\n')
        return fid

    #----------------------------------------------------------------------------------------------
    def print_output_files(self, reactor, fid, time, flag):

        print('{0:12.5e} '.format(time))
        # print output files
        indx = 0
        if 'signal' in self.input:
            # signals
            fid[indx].write('{0:12.5e} '.format(time) + ''.join(['{0:12.5e} '.format(reactor.control.signal[j]) for j in reactor.control.signal]) + '\n')
            indx += 1
        if 'fluid' in reactor.solve:
            # flowrate in dependent and independent junctions (no internal junctions)
            fid[indx].write('{0:12.5e} '.format(time) + ''.join(['{0:12.5e} '.format(reactor.fluid.mdot[i]) for i in range(reactor.fluid.njuni + reactor.fluid.njund)]) + '\n')
            indx += 1
            for i in range(reactor.fluid.npipe):
                fid[indx].write('{0:12.5e} '.format(time) + ''.join(['{0:12.5e} '.format(reactor.fluid.p[i][j]) for j in range(reactor.fluid.pipennodes[i])]) + '\n')
                indx += 1
                fid[indx].write('{0:12.5e} '.format(time) + ''.join(['{0:12.5e} '.format(reactor.fluid.temp[i][j]) for j in range(reactor.fluid.pipennodes[i])]) + '\n')
                indx += 1
                fid[indx].write('{0:12.5e} '.format(time) + ''.join(['{0:12.5e} '.format(reactor.fluid.vel[i][j]) for j in range(reactor.fluid.pipennodes[i])]) + '\n')
                indx += 1
                fid[indx].write('{0:12.5e} '.format(time) + ''.join(['{0:12.5e} '.format(reactor.fluid.re[i][j]) for j in range(reactor.fluid.pipennodes[i])]) + '\n')
                indx += 1
                fid[indx].write('{0:12.5e} '.format(time) + ''.join(['{0:12.5e} '.format(reactor.fluid.pr[i][j]) for j in range(reactor.fluid.pipennodes[i])]) + '\n')
                indx += 1
                fid[indx].write('{0:12.5e} '.format(time) + ''.join(['{0:12.5e} '.format(reactor.fluid.pe[i][j]) for j in range(reactor.fluid.pipennodes[i])]) + '\n')
                indx += 1
            s = ''
            for i in range(reactor.fluid.npipe):
                if reactor.fluid.pipetype[i] == 'freelevel':
                    s += '{0:12.5e} '.format(reactor.fluid.len[i])
            fid[indx].write('{0:12.5e} '.format(time) + s + '\n')
            indx += 1
        if 'fuelrod' in reactor.solve:
            for i in range(reactor.solid.nfuelrods):
                # gas gap conductance
                fid[indx].write('{0:12.5e} '.format(time) + ''.join(['{0:12.5e} '.format(reactor.solid.fuelrod[i].innergas.hgap[j]) for j in range(reactor.solid.fuelrod[i].nz)]) + '\n')
                indx += 1
                # fuel and clad temperatures
                for j in range(reactor.solid.fuelrod[i].nz):
                    fid[indx].write('{0:12.5e} '.format(time) + ''.join(['{0:12.5e} '.format(reactor.solid.fuelrod[i].fuel[j].temp[k]) for k in range(reactor.solid.fuelrod[i].fuel[j].nr)]) + ''.join(['{0:12.5e} '.format(reactor.solid.fuelrod[i].clad[j].temp[k]) for k in range(reactor.solid.fuelrod[i].clad[j].nr)]) + '\n')
                    indx += 1
                    for k in range(reactor.solid.fuelrod[i].fuel[j].nr):
                        if 'fuelgrain' in reactor.solve and i + j + k == 0: 
                            fid[indx].write('{0:12.5e} '.format(time) + ''.join(['{0:12.5e} '.format(reactor.solid.fuelrod[i].fuel[j].fuelgrain[k].c1[l]) for l in range(reactor.solid.fuelrod[i].fuel[j].fuelgrain[k].nr)]) + '\n')
                            indx += 1
                            fid[indx].write('{0:12.5e} '.format(time) + ''.join(['{0:12.5e} '.format(reactor.solid.fuelrod[i].fuel[j].fuelgrain[k].ri[l]) for l in range(reactor.solid.fuelrod[i].fuel[j].fuelgrain[k].NB)]) + '\n')
                            indx += 1
                            fid[indx].write('{0:12.5e} '.format(time) + ''.join(['{0:12.5e} '.format(reactor.solid.fuelrod[i].fuel[j].fuelgrain[k].cv_irr[l]) for l in range(reactor.solid.fuelrod[i].fuel[j].fuelgrain[k].NB)]) + '\n')
                            indx += 1
                            fid[indx].write('{0:12.5e} '.format(time) + ''.join(['{0:12.5e} '.format(reactor.solid.fuelrod[i].fuel[j].fuelgrain[k].ci_irr[l]) for l in range(reactor.solid.fuelrod[i].fuel[j].fuelgrain[k].NB)]) + '\n')
                            indx += 1
                            fid[indx].write('{0:12.5e} '.format(time) + ''.join(['{0:12.5e} '.format(reactor.solid.fuelrod[i].fuel[j].fuelgrain[k].cv_p[l]) for l in range(reactor.solid.fuelrod[i].fuel[j].fuelgrain[k].NB)]) + '\n')
                            indx += 1
                            fid[indx].write('{0:12.5e} '.format(time) + ''.join(['{0:12.5e} '.format(reactor.solid.fuelrod[i].fuel[j].fuelgrain[k].bi[l]) for l in range(reactor.solid.fuelrod[i].fuel[j].fuelgrain[k].NB)]) + '\n')
                            indx += 1
        if 'htstr' in reactor.solve:
            for i in range(reactor.solid.nhtstr):
                fid[indx].write('{0:12.5e} '.format(time) + ''.join(['{0:12.5e} '.format(reactor.solid.htstr[i].temp[k]) for k in range(reactor.solid.htstr[i].nr)]) + '\n')
                indx += 1
        if 'pointkinetics' in reactor.solve:
            # point kinetics power
            fid[indx].write('{0:12.5e} '.format(time) + '{0:12.5e} '.format(reactor.core.power) + '\n')
            indx += 1
            # point kinetics cdnp
            fid[indx].write('{0:12.5e} '.format(time) + ''.join(['{0:12.5e} '.format(reactor.core.cdnp[i]) for i in range(reactor.core.ndnp)]) + '\n')
            indx += 1
        if 'spatialkinetics' in reactor.solve:
            for i in range(reactor.core.niso):
                if reactor.core.iso[i].print_xs:
                    fid[indx].write('time: ' + '{0:12.5e} '.format(time) + ' s\n')
                    nsig0 = len(reactor.core.iso[i].xs['tot'][0][0])
                    ntemp = len(reactor.core.iso[i].xs['tot'][0])
                    for itemp in range(ntemp):
                        fid[indx].write('total XS @' + '{0:12.5e} '.format(reactor.core.iso[i].temp[itemp]) + 'K \n')
                        fid[indx].write(' ' + 'igroup/sig0'.ljust(12) + ''.join(['{0:12.5e} '.format(reactor.core.iso[i].sig0[isig0]) for isig0 in range(nsig0)]) + '\n')
                        for ig in range(reactor.core.ng):
                            fid[indx].write(' ' + str(ig+1).ljust(12) + ''.join(['{0:12.5e} '.format(reactor.core.iso[i].xs['tot'][ig][itemp][isig0]) for isig0 in range(nsig0)]) + '\n')
                    if sum(reactor.core.iso[i].xs['chi']) > 0:
                        for itemp in range(ntemp):
                            fid[indx].write('fission XS @' + '{0:12.5e} '.format(reactor.core.iso[i].temp[itemp]) + 'K \n')
                            fid[indx].write(' ' + 'igroup/sig0'.ljust(12) + ''.join(['{0:12.5e} '.format(reactor.core.iso[i].sig0[isig0]) for isig0 in range(nsig0)]) + '\n')
                            for ig in range(reactor.core.ng):
                                fid[indx].write(' ' + str(ig+1).ljust(12) + ''.join(['{0:12.5e} '.format(reactor.core.iso[i].xs['fis'][ig][itemp][isig0]) for isig0 in range(nsig0)]) + '\n')
                        for itemp in range(ntemp):
                            fid[indx].write('nubar @' + '{0:12.5e} '.format(reactor.core.iso[i].temp[itemp]) + 'K \n')
                            fid[indx].write(' ' + 'igroup'.ljust(12) + '\n')
                            for ig in range(reactor.core.ng):
                                fid[indx].write(' ' + str(ig+1).ljust(12) + '{0:12.5e} '.format(reactor.core.iso[i].xs['nub'][ig][itemp]) + '\n')
                        fid[indx].write('fission spectrum\n')
                        fid[indx].write(' ' + 'igroup'.ljust(12) + 'chi'.ljust(12) + '\n')
                        for ig in range(reactor.core.ng):
                            fid[indx].write(' ' + str(ig+1).ljust(12) + '{0:12.5e} '.format(reactor.core.iso[i].xs['chi'][ig]) + '\n')

                    fid[indx].write('kerma-factors\n')
                    fid[indx].write(' ' + 'igroup/sig0'.ljust(12) + ''.join(['{0:12.5e} '.format(reactor.core.iso[i].sig0[isig0]) for isig0 in range(nsig0)]) + '\n')
                    for ig in range(reactor.core.ng):
                        fid[indx].write(' ' + str(ig+1).ljust(12) + ''.join(['{0:12.5e} '.format(reactor.core.iso[i].xs['kerma'][ig][isig0]) for isig0 in range(nsig0)]) + '\n')

                    for itemp in range(ntemp):
                        fid[indx].write('elastic scattering XS @' + '{0:12.5e} '.format(reactor.core.iso[i].temp[itemp]) + 'K \n')
                        fid[indx].write(' ' + 'from'.ljust(13) + 'to/sig0'.ljust(12) + ''.join(['{0:12.5e} '.format(reactor.core.iso[i].sig0[isig0]) for isig0 in range(nsig0)]) + '\n')
                        for s in reactor.core.iso[i].xs['elan'][0]:
                            fid[indx].write(' ' + str(s[0][0]+1).ljust(13) + str(s[0][1]+1).ljust(12) + ''.join(['{0:12.5e} '.format(s[1][isig0]) for isig0 in range(nsig0)]) + '\n')

                    fid[indx].write('inelastic scattering XS\n')
                    fid[indx].write(' ' + 'from'.ljust(13) + 'to'.ljust(13) + 'sigi'.ljust(12) + '\n')
                    for s in reactor.core.iso[i].xs['ine']:
                        fid[indx].write(' ' + str(s[0][0]+1).ljust(13) + str(s[0][1]+1).ljust(12) + '{0:12.5e} '.format(s[1]) + '\n')
                    if len(reactor.core.iso[i].xs['n2n']) > 0:
                        fid[indx].write('n2n scattering\n')
                        fid[indx].write(' ' + 'from'.ljust(13) + 'to'.ljust(13) + 'sign2n'.ljust(12) + '\n')
                        for s in reactor.core.iso[i].xs['n2n']:
                            fid[indx].write(' ' + str(s[0][0]+1).ljust(13) + str(s[0][1]+1).ljust(12) + '{0:12.5e} '.format(s[1]) + '\n')
                    indx += 1
                    reactor.core.iso[i].print_xs = False
            for i in range(reactor.core.nmix):
                if reactor.core.mix[i].print_xs:
                    fid[indx].write('time: ' + '{0:12.5e} '.format(time) + ' s\n')
                    fid[indx].write('background XS\n')
                    fid[indx].write(' ' + 'igroup'.ljust(13) + ''.join([str(reactor.core.mix[i].isoid[j]).ljust(13) for j in range(reactor.core.mix[i].niso)]) + '\n')
                    for ig in range(reactor.core.ng):
                        fid[indx].write(' ' + str(ig+1).ljust(12) + ''.join(['{0:12.5e} '.format(reactor.core.mix[i].sig0[ig][j]) for j in range(reactor.core.mix[i].niso)]) + '\n')
                    fid[indx].write('total XS, production XS, fission spectrum, in-group scattering XS, out-group scattering XS, n2n XS, kerma-factors\n')
                    fid[indx].write(' ' + 'igroup'.ljust(13) + 'sigt'.ljust(13) + 'nu*sigf'.ljust(13) + 'chi'.ljust(13) + 'sigsi'.ljust(13) + 'sigso'.ljust(13) + 'sign2n'.ljust(13) + 'kerma'.ljust(13) + '\n')
                    for ig in range(reactor.core.ng):
                        sigso = 0
                        sigsi = 0
                        for j in range(len(reactor.core.mix[i].sigsn[0])):
                            f = reactor.core.mix[i].sigsn[0][j][0][0]
                            t = reactor.core.mix[i].sigsn[0][j][0][1]
                            if f == ig and t != ig : sigso = sigso + reactor.core.mix[i].sigsn[0][j][1]
                            if f == ig and t == ig : sigsi = sigsi + reactor.core.mix[i].sigsn[0][j][1]
                        sign2n = 0
                        for j in range(len(reactor.core.mix[i].sign2n)):
                            f = reactor.core.mix[i].sign2n[j][0][0]
                            t = reactor.core.mix[i].sign2n[j][0][1]
                            if f == ig and t != ig : sign2n = sign2n + reactor.core.mix[i].sign2n[j][1]
                        fid[indx].write(' ' + str(ig+1).ljust(12) + str('{0:12.5e} '.format(reactor.core.mix[i].sigt[ig])) + str('{0:12.5e} '.format(reactor.core.mix[i].sigp[ig])) + str('{0:12.5e} '.format(reactor.core.mix[i].chi[ig])) + str('{0:12.5e} '.format(sigsi)) + str('{0:12.5e} '.format(sigso)) + str('{0:12.5e} '.format(sign2n)) + str('{0:12.5e} '.format(reactor.core.mix[i].kerma[ig])) + '\n')
                    fid[indx].write('scattering XS\n')
                    fid[indx].write(' ' + 'from'.ljust(13) + 'to'.ljust(13) + 'sigs'.ljust(13) + '\n')
                    for j in range(len(reactor.core.mix[i].sigsn[0])):
                        f = reactor.core.mix[i].sigsn[0][j][0][0] + 1
                        t = reactor.core.mix[i].sigsn[0][j][0][1] + 1
                        sigs = reactor.core.mix[i].sigsn[0][j][1]
                        fid[indx].write(' ' + str(f).ljust(13) + str(t).ljust(12) + '{0:12.5e} '.format(sigs) + '\n')
                    fid[indx].write('n2n XS\n')
                    fid[indx].write(' ' + 'from'.ljust(13) + 'to'.ljust(13) + 'sign2n'.ljust(13) + '\n')
                    for j in range(len(reactor.core.mix[i].sign2n)):
                        f = reactor.core.mix[i].sign2n[j][0][0] + 1
                        t = reactor.core.mix[i].sign2n[j][0][1] + 1
                        sign2n = reactor.core.mix[i].sign2n[j][1]
                        fid[indx].write(' ' + str(f).ljust(13) + str(t).ljust(12) + '{0:12.5e} '.format(sign2n) + '\n')
                    indx += 1
                    reactor.core.mix[i].print_xs = False

                else:
                    indx += 7
            # multiplication factor
            if flag == 0 : fid[indx].write(''.join([(' '+str(niter)).ljust(13) + '{0:12.5e} '.format(reactor.core.k[niter]) + '\n' for niter in range(len(reactor.core.k))]))
            indx += 1
            # neutron flux
            if flag == 0 : 
                for iz in range(reactor.core.nz):
                    for ix in range(reactor.core.nx):
                        for iy in range(reactor.core.ny):
                            imix = reactor.core.map['imix'][iz][ix][iy]
                            # if (iz, ix, iy) is not a boundary condition node, i.e. not -1 (vac) and not -2 (ref')
                            if imix >= 0:
                                for ig in range(reactor.core.ng):
                                    flux = sum([reactor.core.flux[iz][ix][iy][it][ig] for it in range(reactor.core.nt)])
                                    fid[indx].write('{0:12.5e} '.format(time) + ' ' + str(ig+1).ljust(13) + str(iz).ljust(13) + str(ix).ljust(13) + str(iy).ljust(12) + '{0:12.5e} '.format(flux) + '\n')
            indx += 1
            # power
            if flag == 0 : 
                for iz in range(reactor.core.nz):
                    for ix in range(reactor.core.nx):
                        for iy in range(reactor.core.ny):
                            imix = reactor.core.map['imix'][iz][ix][iy]
                            # if (iy, ix, iz) is not a boundary condition node, i.e. not -1 (vac) and not -2 (ref')
                            if imix >= 0 and reactor.core.pow[iz][ix][iy] > 0:
                                fid[indx].write('{0:12.5e} '.format(time) + ' ' + str(iz).ljust(13) + str(ix).ljust(13) + str(iy).ljust(12) + '{0:12.5e} '.format(reactor.core.pow[iz][ix][iy]) + '\n')
            indx += 1
            if flag == 0 : 
                for ix in range(reactor.core.nx):
                    for iy in range(reactor.core.ny):
                        if reactor.core.powxy[ix][iy] > 0:
                            fid[indx].write('{0:12.5e} '.format(time) + ' ' + str(ix).ljust(13) + str(iy).ljust(12) + '{0:12.5e} '.format(reactor.core.powxy[ix][iy]) + '\n')
            indx += 1

    #----------------------------------------------------------------------------------------------
    def write_to_y(self, reactor):

        # write list of unknowns to y
        y = []
        if 'fluid' in reactor.solve:
            k = 0
            for j in range(reactor.fluid.njun):
                if reactor.fluid.juntype[j] == 'independent':
                    # flowrate in independent junctions
                    y.append(reactor.fluid.mdoti[k])
                    k += 1
            for i in range(reactor.fluid.npipe):
                if reactor.fluid.pipetype[i] == 'freelevel':
                    # free-level-volume length
                    y.append(reactor.fluid.len[i])
            for i in range(reactor.fluid.npipe):
                for j in range(reactor.fluid.pipennodes[i]):
                    # temperature in pipe nodes
                    y.append(reactor.fluid.temp[i][j])
        if 'fuelrod' in reactor.solve:
            for i in range(reactor.solid.nfuelrods):
                for j in range(reactor.solid.fuelrod[i].nz):
                    for k in range(reactor.solid.fuelrod[i].fuel[j].nr):
                        if 'fuelgrain' in reactor.solve and i + j + k == 0: #i+j+k==0 is a temporal condition to solve fuel grain only for one node
                            # fuel grain monoatoms
                            for l in range(reactor.solid.fuelrod[i].fuel[j].fuelgrain[k].nr):
                                y.append(reactor.solid.fuelrod[i].fuel[j].fuelgrain[k].c1[l])
                            # fuel grain bubble radii
                            for l in range(reactor.solid.fuelrod[i].fuel[j].fuelgrain[k].NB):
                                y.append(reactor.solid.fuelrod[i].fuel[j].fuelgrain[k].ri[l])
                            # fuel grain fractional concentration of irradiation-induced uranium vacancies
                            for l in range(reactor.solid.fuelrod[i].fuel[j].fuelgrain[k].NB):
                                y.append(reactor.solid.fuelrod[i].fuel[j].fuelgrain[k].cv_irr[l])
                            # fuel grain fractional concentration of irradiation-induced uranium interstitials
                            for l in range(reactor.solid.fuelrod[i].fuel[j].fuelgrain[k].NB):
                                y.append(reactor.solid.fuelrod[i].fuel[j].fuelgrain[k].ci_irr[l])
                            # fuel grain fractional concentration of uranium vacancies ejected from intragranular as-fabricated pores
                            for l in range(reactor.solid.fuelrod[i].fuel[j].fuelgrain[k].NB):
                                y.append(reactor.solid.fuelrod[i].fuel[j].fuelgrain[k].cv_p[l])
                            # fuel grain intragranular bubble concentation
                            for l in range(reactor.solid.fuelrod[i].fuel[j].fuelgrain[k].NB):
                                y.append(reactor.solid.fuelrod[i].fuel[j].fuelgrain[k].bi[l])
                    for k in range(reactor.solid.fuelrod[i].fuel[j].nr):
                        # fuel temperature
                        y.append(reactor.solid.fuelrod[i].fuel[j].temp[k])
                    for k in range(reactor.solid.fuelrod[i].clad[j].nr):
                        # clad temperature
                        y.append(reactor.solid.fuelrod[i].clad[j].temp[k])
        if 'htstr' in reactor.solve:
            for i in range(reactor.solid.nhtstr):
                for k in range(reactor.solid.htstr[i].nr):
                    # htstr temperature
                    y.append(reactor.solid.htstr[i].temp[k])
        if 'pointkinetics' in reactor.solve:
            y.append(reactor.core.power)
            for i in range(reactor.core.ndnp):
                y.append(reactor.core.cdnp[i])
        if 'spatialkinetics' in reactor.solve:
            for iz in range(reactor.core.nz):
                for ix in range(reactor.core.nx):
                    for iy in range(reactor.core.ny):
                        # if (iy, ix, iz) is not a boundary condition node, i.e. not -1 (vac) and not -2 (ref)
                        imix = reactor.core.map['imix'][iz][ix][iy]
                        if imix >= 0 and any(reactor.core.mix[imix].sigf) > 0:
                            for it in range(reactor.core.nt):
                                for ig in range(reactor.core.ng):
                                    y.append(reactor.core.flux[iz][ix][iy][it][ig])
        return y

    #----------------------------------------------------------------------------------------------
    def read_from_y(self, reactor, y):

        # read list of unknowns from y
        indx = 0
        if 'fluid' in reactor.solve:
            k = 0
            for j in range(reactor.fluid.njun):
                if reactor.fluid.juntype[j] == 'independent':
                    f = reactor.fluid.f[j][0]
                    t = reactor.fluid.t[j][0]
                    # tuple of from-to pipe id's
                    f_t = (reactor.fluid.pipeid[f],reactor.fluid.pipeid[t])
                    try:
                        # check if the current junction j is present in junflowrate list
                        reactor.fluid.junflowrate['jun'].index(f_t)
                        # if yes...
                        pass
                    except ValueError:
                        # flowrate in independent junctions
                        reactor.fluid.mdoti[k] = y[indx]
                    k += 1
                    indx += 1
            for i in range(reactor.fluid.npipe):
                if reactor.fluid.pipetype[i] == 'freelevel':
                    # free-level-volume length
                    reactor.fluid.len[i] = y[indx]
                    indx += 1
            for i in range(reactor.fluid.npipe):
                for j in range(reactor.fluid.pipennodes[i]):
                    # temperature in pipe nodes
                    reactor.fluid.temp[i][j] = y[indx]
                    indx += 1
        if 'fuelrod' in reactor.solve:
            for i in range(reactor.solid.nfuelrods):
                for j in range(reactor.solid.fuelrod[i].nz):
                    for k in range(reactor.solid.fuelrod[i].fuel[j].nr):
                        if 'fuelgrain' in reactor.solve and i + j + k == 0:
                            for l in range(reactor.solid.fuelrod[i].fuel[j].fuelgrain[k].nr):
                                # fuel grain monoatoms
                                reactor.solid.fuelrod[i].fuel[j].fuelgrain[k].c1[l] = y[indx]
                                indx += 1
                            for l in range(reactor.solid.fuelrod[i].fuel[j].fuelgrain[k].NB):
                                # fuel grain bubble radii
                                reactor.solid.fuelrod[i].fuel[j].fuelgrain[k].rb[l] = y[indx]
                                indx += 1
                            for l in range(reactor.solid.fuelrod[i].fuel[j].fuelgrain[k].NB):
                                # fuel grain fractional concentration of irradiation-induced uranium vacancies
                                reactor.solid.fuelrod[i].fuel[j].fuelgrain[k].cv_irr[l] = y[indx]
                                indx += 1
                            for l in range(reactor.solid.fuelrod[i].fuel[j].fuelgrain[k].NB):
                                # fuel grain fractional concentration of irradiation-induced uranium interstitials
                                reactor.solid.fuelrod[i].fuel[j].fuelgrain[k].ci_irr[l] = y[indx]
                                indx += 1
                            for l in range(reactor.solid.fuelrod[i].fuel[j].fuelgrain[k].NB):
                                # fuel grain fractional concentration of uranium vacancies ejected from intragranular as-fabricated pores
                                reactor.solid.fuelrod[i].fuel[j].fuelgrain[k].cv_p[l] = y[indx]
                                indx += 1
                            for l in range(reactor.solid.fuelrod[i].fuel[j].fuelgrain[k].NB):
                                # fuel grain intragranular bubble concentrations
                                reactor.solid.fuelrod[i].fuel[j].fuelgrain[k].bi[l] = y[indx]
                                indx += 1
                    for k in range(reactor.solid.fuelrod[i].fuel[j].nr):
                        # fuel temperature
                        reactor.solid.fuelrod[i].fuel[j].temp[k] = y[indx]
                        indx += 1
                    for k in range(reactor.solid.fuelrod[i].clad[j].nr):
                        # clad temperature
                        reactor.solid.fuelrod[i].clad[j].temp[k] = y[indx]
                        indx += 1
        if 'htstr' in reactor.solve:
            for i in range(reactor.solid.nhtstr):
                for k in range(reactor.solid.htstr[i].nr):
                    # htstr temperature
                    reactor.solid.htstr[i].temp[k] = y[indx]
                    indx += 1
        if 'pointkinetics' in reactor.solve:
            reactor.core.power = y[indx]
            indx += 1
            for i in range(reactor.core.ndnp):
                reactor.core.cdnp[i] = y[indx]
                indx += 1
        if 'spatialkinetics' in reactor.solve:
            for iz in range(reactor.core.nz):
                for ix in range(reactor.core.nx):
                    for iy in range(reactor.core.ny):
                        # if (iy, ix, iz) is not a boundary condition node, i.e. not -1 (vac) and not -2 (ref)
                        imix = reactor.core.map['imix'][iz][ix][iy]
                        if imix >= 0 and any(reactor.core.mix[imix].sigf) > 0:
                            for it in range(reactor.core.nt):
                                for ig in range(reactor.core.ng):
                                    reactor.core.flux[iz][ix][iy][it][ig] = y[indx]
                                    indx += 1
