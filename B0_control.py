from scipy.interpolate import interp1d

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

    #------------------------------------------------------------------------------------------
    def evaluate_signals(self, reactor, t):

        # evaluate signals
        self.signal = {}
        for s in self.input['signal'] :
            if type(s['value']) == int or type(s['value']) == float:
                self.signal[s['id']] = s['value']
            if s['value'] == 'time' :
                self.signal[s['id']] = t

        #evaluate output signals of lookup tables
        lookup_table = self.input['lookup']
        for table in lookup_table :
            insignal_name = table['x'][0]
            outsignal_name = table['f(x)'][0]
            x = table['x'][1:]
            y = table['f(x)'][1:]
            f = interp1d(x, y) #scipy function
            xnew = max(min(self.signal[insignal_name],x[-1]),x[0])
            ynew = f(xnew)
            self.signal[outsignal_name] = ynew

        # signal-dependent junction: impose flowrate
        if 'fluid' in reactor.solve:
            for j in range(reactor.fluid.njun):
                if reactor.fluid.juntype[j] == 'independent' and reactor.fluid.junflowrate[j] != '':
                    # impose flowrate from the look-up table
                    reactor.fluid.mdoti[j] = self.signal[reactor.fluid.junflowrate[j]]
        
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
        inp['clad'] = []
        inp['fuel'] = []
        inp['fuelrod'] = []
        inp['innergas'] = []
        inp['junction'] = {'from':[], 'to':[], 'type':[], 'pumphead':[], 'flowrate':[]}
        inp['lookup'] = []
        inp['mat'] = []
        inp['mix'] = []
        inp['p2d'] = []
        inp['pipe'] = []
        inp['signal'] = []
        inp['signalid'] = []
        inp['solve'] = []
        inp['t0'] = 0
        inp['t_dt'] = []
    
        #read input file as a whole
        f = open('input', 'r')
        s0 = f.read()
        f.close()
    
        #merge &-ending "line" with the next one
        s = ''
        take = True
        for c in s0 :
            if c == '&' : take = False
            if take : s += c
            if c == '\n' : take = True
    
        #split in lines
        lines = s.strip().split('\n')
    
        #remove comment-lines (*)
        lines = [x for x in lines if not x.startswith('*')]
    
        def convert_to_float(w) : 
            try:
                w = float(w)
            except :
                pass
            return w
    
        for line in lines:
            word = line.split()
            word = list(map(convert_to_float, word))
    
            key = word[0].lower()
            #--------------------------------------------------------------------------------------
            # just placeholder
            if key == '' :
                pass
            #--------------------------------------------------------------------------------------
            # effective delayed neutron fractions
            elif key == 'betaeff' :
                inp['betaeff'] = word[1:]
            #--------------------------------------------------------------------------------------
            # cladding
            elif key == 'clad' :
                 inp['clad'].append( {'id':word[1], 'matid':word[2], 'ri':word[3], 'ro':word[4], 'nr':int(word[5])} )
            #--------------------------------------------------------------------------------------
            # constant
            elif key == 'constant' :
                 inp[word[1]] = float(word[2])
            #--------------------------------------------------------------------------------------
            # delayed neutron precursor decay time constants
            elif key == 'dnplmb' :
                inp['dnplmb'] = word[1:]
            #--------------------------------------------------------------------------------------
            # fuel grain parameters
            elif key == 'fgrain' :
                # grain diameter
                inp['dgrain'] = word[1]
                # number of nodes in the grain
                inp['nrgrain'] = int(word[2])
                # fission rate
                inp['frate'] = int(word[3])
            #--------------------------------------------------------------------------------------
            # fuel
            elif key == 'fuel' :
                 inp['fuel'].append( {'id':word[1], 'matid':word[2], 'ri':word[3], 'ro':word[4], 'nr':int(word[5])} )
            #--------------------------------------------------------------------------------------
            # fuel rod card
            elif key == 'fuelrod' :
                id = word[1]
                if any([id in x['id'] for x in inp['fuelrod']]):
                    for x in inp['fuelrod']:
                        if x['id'] == id:
                            x['fuelid'].append(word[2])
                            x['hgap'].append(float(word[3]))
                            x['cladid'].append(word[4])
                            x['p2d'].append(word[5])
                            x['mltpl'].append(word[6])
                            x['pipeid'].append(word[7])
                            x['pipenodeid'].append(int(word[8]))
                else:
                    inp['fuelrod'].append({'id':id, 'fuelid':[word[2]], 'hgap':[float(word[3])], 'cladid':[word[4]], 'p2d':[word[5]], 'mltpl':[word[6]], 'pipeid':[word[7]], 'pipenodeid':[int(word[8])]})
            #--------------------------------------------------------------------------------------
            # inner gas
            elif key == 'innergas' :
                 inp['innergas'].append( {'fuelrodid':word[1], 'matid':word[2], 'plenv':word[3]} )
            #--------------------------------------------------------------------------------------
            # thermal-hydraulic junction (dependent)
            elif key == 'jun' :
                 inp['junction']['from'].append(word[1])
                 inp['junction']['to'].append(word[2])
                 inp['junction']['type'].append('dependent')
                 inp['junction']['pumphead'].append('')
                 inp['junction']['flowrate'].append('')
            #--------------------------------------------------------------------------------------
            # thermal-hydraulic junction (independent)
            elif key == 'jun-i' :
                 inp['junction']['from'].append(word[1])
                 inp['junction']['to'].append(word[2])
                 inp['junction']['type'].append('independent')
                 inp['junction']['pumphead'].append('')
                 inp['junction']['flowrate'].append('')
            #--------------------------------------------------------------------------------------
            # thermal-hydraulic junction (independent + signal for flowrate)
            elif key == 'jun-i-f' :
                 inp['junction']['from'].append(word[1])
                 inp['junction']['to'].append(word[2])
                 inp['junction']['type'].append('independent')
                 inp['junction']['pumphead'].append('')
                 inp['junction']['flowrate'].append(word[3])
            #--------------------------------------------------------------------------------------
            # thermal-hydraulic junction (independent + signal for pump head)
            elif key == 'jun-i-p' :
                 inp['junction']['from'].append(word[1])
                 inp['junction']['to'].append(word[2])
                 inp['junction']['type'].append('independent')
                 inp['junction']['pumphead'].append(word[3])
                 inp['junction']['flowrate'].append('')
            #--------------------------------------------------------------------------------------
            # lookup table
            elif key == 'lookup' :
                 lookup = {}
                 lookup['x'] = word[1::2]
                 lookup['f(x)'] = word[2::2]
                 inp['lookup'].append(lookup)
            #--------------------------------------------------------------------------------------
            # material
            elif key == 'mat' :
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
            elif key == 'mix' :
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
            elif key == 'nddir' :
                 inp['nddir'] = word[1]
            #--------------------------------------------------------------------------------------
            # thermal-hydraulic pipe without free level
            elif key == 'pipe' :
                 inp['pipe'].append( {'id':word[1], 'type':'normal', 'matid':word[2], 'dhyd':word[3], 'len':word[4], 'dir':word[5], 'areaz':word[6], 'nnodes':int(word[7]), 'signaltemp':''} )
            #--------------------------------------------------------------------------------------
            # thermal-hydraulic pipe with free level
            elif key == 'pipe-f' :
                 inp['pipe'].append( {'id':word[1], 'type':'freelevel', 'matid':word[2], 'dhyd':word[3], 'len':word[4], 'dir':word[5], 'areaz':word[6], 'nnodes':1, 'signaltemp':''} )
            #--------------------------------------------------------------------------------------
            # thermal-hydraulic pipe without free level with temperature defined by signal
            elif key == 'pipe-t' :
                 inp['pipe'].append( {'id':word[1], 'type':'normal', 'matid':word[2], 'dhyd':word[3], 'len':word[4], 'dir':word[5], 'areaz':word[6], 'nnodes':int(word[7]), 'signaltemp':word[8]} )
            #--------------------------------------------------------------------------------------
            elif key == 'solve':
                inp['solve'].append(word[1])
                # verify that solve card has correct value
                correct_values = {'fluid','fuelgrain','fuelrod','pointkinetics','spatialkinetics'}
                value = set([word[1]])
                diff = value.difference(correct_values)
                if diff != set():
                    print('****ERROR: solve card contains wrong value: ', list(diff)[0], '\nCorrect values are: ')
                    sorted = list(correct_values)
                    sorted.sort()
                    for v in sorted:
                        print('solve', v)
                    sys.exit()
                if word[1] == 'spatialkinetics':
                    # check that there is a second value
                    if len(word)-1 == 1:
                        print('****ERROR: solve spatialkinetics card should have a second value after the keyword: number of energy groups (integer), e.g.:\nsolve spatialkinetics 25')
                        sys.exit()
                    # check that there the second value is integer
                    try:
                        # number of energy groups
                        inp['ng'] = int(word[2])
                    except:
                        print('****ERROR: the second value after the keyword of solve spatialkinetics card should be integer (number of energy groups), e.g.:\nsolve spatialkinetics 25')
                        sys.exit()
            #--------------------------------------------------------------------------------------
            # signal variable
            elif key == 'signal' :
                 signal = {}
                 signal['id'] = word[1]
                 signal['value'] = word[2]
                 inp['signal'].append(signal)
            #--------------------------------------------------------------------------------------
            # integration starting time
            elif key == 't0' :
                inp['t0'] = word[1]
            #--------------------------------------------------------------------------------------
            # end of time interval and output time step for this interval
            elif key == 't_dt' :
                inp['t_dt'].append([word[1], word[2]])
            #--------------------------------------------------------------------------------------
            # prompt neutron lifetime
            elif key == 'tlife' :
                inp['tlife'] = word[1]
    
        # verify that t_dt present
        if inp['t_dt'] == [] :
            sys.exit('****ERROR: obligatory card t_dt specifying time_end and dtime_out is absent.')
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
        for table in inp['lookup'] :
            insignal = table['x'][0]
            outsignal = table['f(x)'][0]
            if insignal not in inp['signalid'] :
                print('****ERROR: input signal ' + insignal + ' in lookup table ' + outsignal + ' is not defined.')
                sys.exit()
        # append output signals of lookup tables
        inp['signalid'] += [y['f(x)'][0] for y in inp['lookup']]
    
        # verify that mix card uses existing signals
        for s in [x['signaltemp'][j] for x in inp['mix'] for j in range(len(x['signaltemp']))]:
            if s not in inp['signalid']:
                print('****ERROR: signal for temperature ' + s + ' in mix card is not defined.')
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

        # copy input and open output files to output folder
        shutil.copyfile('input', path4results + os.sep + 'input')
        shutil.copyfile('input.json', path4results + os.sep + 'input.json')
        # open files for output
        fid = []
        if 'fuelrod' in reactor.solve:
            for i in range(reactor.solid.nfuelrods):
                fid.append(open(path4results + os.sep + 'fuelrod-hgap-' + [x['id'] for x in self.input['fuelrod']][i] + '.dat', 'w'))
                fid[-1].write(' ' + 'time(s)'.ljust(13) + ''.join([('hgap-' + str(j).zfill(3)).ljust(13) for j in range(reactor.solid.fuelrod[i].nz)]) + '\n')
                for j in range(reactor.solid.fuelrod[i].nz):
                    fid.append(open(path4results + os.sep + 'fuelrod-temp-' + [x['id'] for x in self.input['fuelrod']][i] + '-' + str(j).zfill(3) + '.dat', 'w'))
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
        if 'pointkinetics' in reactor.solve:
            fid.append(open(path4results + os.sep + 'core-power.dat', 'w'))
            fid[-1].write(' ' + 'time(s)'.ljust(13) + 'power(-)\n')
            fid.append(open(path4results + os.sep + 'core-cdnp.dat', 'w'))
            fid[-1].write(' ' + 'time(s)'.ljust(13) + ''.join([('cdnp-' + str(i)).ljust(13) for i in range(reactor.core.ndnp)]) + '\n')
        if 'spatialkinetics' in reactor.solve:
            for i in range(reactor.core.nmix):
                fid.append(open(path4results + os.sep + 'core-mix-' + str(i).zfill(4) + '-sig0.dat', 'w'))
                fid[-1].write(' ' + 'time(s)'.ljust(12) + 'isoname'.ljust(13) + ''.join([(str(j+1)).ljust(13) for j in range(reactor.core.mix[i].ng)]) + '\n')
                fid.append(open(path4results + os.sep + 'core-mix-' + str(i).zfill(4) + '-sigt.dat', 'w'))
                fid[-1].write(' ' + 'time(s)'.ljust(13) + ''.join([(str(j+1)).ljust(13) for j in range(reactor.core.mix[i].ng)]) + '\n')
                fid.append(open(path4results + os.sep + 'core-mix-' + str(i).zfill(4) + '-siga.dat', 'w'))
                fid[-1].write(' ' + 'time(s)'.ljust(13) + ''.join([(str(j+1)).ljust(13) for j in range(reactor.core.mix[i].ng)]) + '\n')
                fid.append(open(path4results + os.sep + 'core-mix-' + str(i).zfill(4) + '-sigs.dat', 'w'))
                fid[-1].write(' ' + 'time(s)'.ljust(13) + 'from'.ljust(13) + 'to'.ljust(13) + 'sigs'.ljust(13) + '\n')

        return fid

    #----------------------------------------------------------------------------------------------
    def print_output_files(self, reactor, fid, time):

        # print output files
        indx = 0
        if 'fuelrod' in reactor.solve:
            for i in range(reactor.solid.nfuelrods):
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
        if 'pointkinetics' in reactor.solve:
            # point kinetics power
            fid[indx].write('{0:12.5e} '.format(time) + '{0:12.5e} '.format(reactor.core.power) + '\n')
            indx += 1
            # point kinetics cdnp
            fid[indx].write('{0:12.5e} '.format(time) + ''.join(['{0:12.5e} '.format(reactor.core.cdnp[i]) for i in range(reactor.core.ndnp)]) + '\n')
            indx += 1
        if 'spatialkinetics' in reactor.solve:
            for i in range(reactor.core.nmix):
                for j in range(reactor.core.mix[i].niso):
                    # sigma-zeros
                    fid[indx].write('{0:12.5e} '.format(time) + str(reactor.core.mix[i].isoid[j]).ljust(12) + ''.join(['{0:12.5e} '.format(reactor.core.mix[i].sig0[ig][j]) for ig in range(reactor.core.mix[i].ng)]) + '\n')
                indx += 1
                # macroscopic sigma-total
                fid[indx].write('{0:12.5e} '.format(time) + ''.join(['{0:12.5e} '.format(reactor.core.mix[i].sigt[ig]) for ig in range(reactor.core.mix[i].ng)]) + '\n')
                indx += 1
                # macroscopic sigma-absorption
                fid[indx].write('{0:12.5e} '.format(time) + ''.join(['{0:12.5e} '.format(reactor.core.mix[i].siga[ig]) for ig in range(reactor.core.mix[i].ng)]) + '\n')
                indx += 1
                # macroscopic sigma-scattering
                for j in range(len(reactor.core.mix[i].sigs)):
                    f = reactor.core.mix[i].sigs[j][0][0] + 1
                    t = reactor.core.mix[i].sigs[j][0][1] + 1
                    sigs = reactor.core.mix[i].sigs[j][1]
                    fid[indx].write('{0:12.5e} '.format(time) + ' ' + str(f).ljust(13) + str(t).ljust(12) + '{0:12.5e} '.format(sigs) + '\n')
                indx += 1

    #----------------------------------------------------------------------------------------------
    def write_to_y(self, reactor):

        # write list of unknowns to y
        y = []
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
        if 'fluid' in reactor.solve:
            for j in range(reactor.fluid.njun):
                if reactor.fluid.juntype[j] == 'independent':
                    # flowrate in independent junctions
                    y.append(reactor.fluid.mdoti[j])
            for i in range(reactor.fluid.npipe):
                for j in range(reactor.fluid.pipennodes[i]):
                    # temperature in pipe nodes
                    y.append(reactor.fluid.temp[i][j])
        if 'pointkinetics' in reactor.solve:
            y.append(reactor.core.power)
            for i in range(reactor.core.ndnp):
                y.append(reactor.core.cdnp[i])
        return y

    #----------------------------------------------------------------------------------------------
    def read_from_y(self, reactor, y):

        # read list of unknowns from y
        indx = 0
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
        if 'fluid' in reactor.solve:
            for j in range(reactor.fluid.njun):
                if reactor.fluid.juntype[j] == 'independent':
                    # flowrate in independent junctions
                    reactor.fluid.mdoti[j] = y[indx]
                    indx += 1
            for i in range(reactor.fluid.npipe):
                for j in range(reactor.fluid.pipennodes[i]):
                    # temperature in pipe nodes
                    reactor.fluid.temp[i][j] = y[indx]
                    indx += 1
        if 'pointkinetics' in reactor.solve:
            reactor.core.power = y[indx]
            indx += 1
            for i in range(reactor.core.ndnp):
                reactor.core.cdnp[i] = y[indx]
                indx += 1
