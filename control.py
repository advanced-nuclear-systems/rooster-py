import sys

#--------------------------------------------------------------------------------------------------
class Control:
    
    neq = 0

    def __init__(self, reactor):
        self.input = construct_input()

    def calculate_rhs(self,reactor, t, y):
        rhs = []
        return rhs

#--------------------------------------------------------------------------------------------------
def construct_input():
    #create dictionary inp where all input data will be stored
    inp = {}
    inp['t0'] = 0 #default
    inp['t_dt'] = [] #no default
    inp['signal'] = [] #no default
    inp['lookup'] = [] #no default

#    #read input file line by line
#    f = open('input', mode = 'r')
#    lines = f.readlines()
#    f.close()

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
        if key == 't0' :
            inp['t0'] = word[1]
        #--------------------------------------------------------------------------------------
        elif key == 't_dt' :
            inp['t_dt'].append([word[1], word[2]])
        #--------------------------------------------------------------------------------------
        elif key == 'signal' :
             signal = {}
             signal['type'] = word[1]
             signal['userid'] = word[2]
             signal['sign'] = word[3:]
             inp['signal'].append(signal)
        #--------------------------------------------------------------------------------------
        elif key == 'lookup' :
             lookup = {}
             lookup['x'] = word[1::2]
             lookup['f(x)'] = word[2::2]
             inp['lookup'].append(lookup)
    #verification
    if inp['t_dt'] == [] :
        print('****ERROR: obligatory card t_dt specifying time_end and dtime_out is absent.')
        sys.exit()

    print(inp)
    return inp
