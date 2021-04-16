from scipy.interpolate import interp1d

import sys

#--------------------------------------------------------------------------------------------------
class Control:
    
    neq = 0
    signal = {}

    # constructor: self is a 'control' object created in B
    def __init__(self, reactor):
        self.input = construct_input()
        self.state = []
        self.neq = len(self.state)

    # create right-hand side vector: self is a 'control' object created in B
    def calculate_rhs(self, reactor, t):
        rhs = []
        return rhs

    def evaluate(self, reactor, t):
        # evaluate signals
        for s in self.input['signal'] :
            if s['type'] == 'time' :
                self.signal[s['userid']] = t

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

#--------------------------------------------------------------------------------------------------
def construct_input():
    #create dictionary inp where all input data will be stored
    inp = {}
    inp['coolant'] = {'name':[], 'type':[], 'p0':[], 'temp0':[]} # no default
    inp['fuel'] = {'name':[], 'type':[], 'pu':[], 'b':[], 'x':[], 'por':[], 'temp0':[]} # no default
    inp['junction'] = {'from':[], 'to':[], 'type':[]} # no default
    inp['lookup'] = [] # no default
    inp['pellet'] = {'name':[], 'ri':[], 'ro':[], 'dz':[], 'nr':[]} # no default
    inp['pipe'] = {'name':[], 'type':[], 'cool':[], 'dhyd':[], 'elev':[], 'len':[], 'areaz':[], 'nnodes':[]} # no default
    inp['pnltime'] = '' # no default
    inp['signal'] = [] # no default
    inp['solve'] = [] # no default
    inp['t0'] = 0 # default
    inp['t_dt'] = [] # no default

    #read input file as a whole
    f = open('input', mode = 'r')
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
        # constant
        elif key == 'constant' :
             inp[word[1]] = float(word[2])
        #--------------------------------------------------------------------------------------
        # effective delayed neutron fractions
        elif key == 'betaeff' :
            inp['betaeff'] = word[1:]
        #--------------------------------------------------------------------------------------
        # coolant
        elif key == 'coolant' :
             inp['coolant']['name'].append(word[1])
             inp['coolant']['type'].append(word[2])
             inp['coolant']['p0'].append(word[3])
             inp['coolant']['temp0'].append(word[4])
        #--------------------------------------------------------------------------------------
        # delayed neutron precursor decay time constants
        elif key == 'dnplmb' :
            inp['dnplmb'] = word[1:]
        #--------------------------------------------------------------------------------------
        # fuel
        elif key == 'fuel' :
             inp['fuel']['name'].append(word[1])
             inp['fuel']['type'].append(word[2])
             inp['fuel']['pu'].append(word[3]) # Pu content (-)
             inp['fuel']['b'].append(word[4]) # burnup (MWd/kgU)
             inp['fuel']['x'].append(word[5]) # deviation from stoechiometry
             inp['fuel']['por'].append(word[6]) # porosity
             inp['fuel']['temp0'].append(word[7]) # initial temperature (K)
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
        # thermal-hydraulic junction
        elif key == 'junction' :
             inp['junction']['from'].append(word[1])
             inp['junction']['to'].append(word[2])
             inp['junction']['type'].append(word[3])
        #--------------------------------------------------------------------------------------
        # lookup table
        elif key == 'lookup' :
             lookup = {}
             lookup['x'] = word[1::2]
             lookup['f(x)'] = word[2::2]
             inp['lookup'].append(lookup)
        #--------------------------------------------------------------------------------------
        # fuel pellet
        elif key == 'pellet' :
             inp['pellet']['name'].append(word[1])
             inp['pellet']['ri'].append(word[2])
             inp['pellet']['ro'].append(word[3])
             inp['pellet']['dz'].append(word[4])
             inp['pellet']['nr'].append(int(word[5]))
        #--------------------------------------------------------------------------------------
        # thermal-hydraulic pipe
        elif key == 'pipe' :
             inp['pipe']['name'].append(word[1])
             inp['pipe']['type'].append(word[2])
             inp['pipe']['cool'].append(word[3])
             inp['pipe']['dhyd'].append(word[4])
             inp['pipe']['elev'].append(word[5])
             inp['pipe']['len'].append(word[6])
             inp['pipe']['areaz'].append(word[7])
             inp['pipe']['nnodes'].append(int(word[8]))
        #--------------------------------------------------------------------------------------
        # 
        elif key == 'solve' :
            inp['solve'].append(word[1:])
            # verify that solve card has correct value
            correct_values = {'fluid','fuelgrain','fuelrod','pointkinetics'}
            values = set([word[1]])
            diff = values.difference(correct_values)
            if diff != set():
                print('****ERROR: solve card contains wrong value: ', list(diff)[0], '\nCorrect values are: ')
                sorted = list(correct_values)
                sorted.sort()
                for v in sorted:
                    print('solve', v)
                sys.exit()
        #--------------------------------------------------------------------------------------
        # signal variable
        elif key == 'signal' :
             signal = {}
             signal['type'] = word[1]
             signal['userid'] = word[2]
             signal['sign'] = word[3:]
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

    # verify that there is at least one solve card
    if len(inp['solve']) == 0:
        print('****ERROR: input file should have at least one solve card.')
        sys.exit()
    
    # verify that lookup tables use existing signals
    signal_userid = []
    for s in inp['signal'] :
        signal_userid.append(s['userid'])
    for table in inp['lookup'] :
        insignal = table['x'][0]
        outsignal = table['f(x)'][0]
        if insignal not in signal_userid :
            print('****ERROR: input signal ' + insignal + ' in lookup table ' + outsignal + ' is not defined.')
            sys.exit()
    #print(inp)
    return inp
