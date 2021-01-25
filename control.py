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

    #read input file line by line
    f = open('input', mode = 'r')
    lines = f.readlines()
    f.close()

    for line in lines:
        if line[0] is not '*' :
            word = line.split()
            key = word[0].lower()
            #--------------------------------------------------------------------------------------
            if key == 't0' :
                inp['t0'] = float(word[1])
            #--------------------------------------------------------------------------------------
            elif key == 't_dt' :
                inp['t_dt'].append([float(word[1]), float(word[2])])

    #verification
    if inp['t_dt'] == [] :
        print('ERROR: Card t_dt specifying time_end and dtime_out is obligatory')
        sys.exit()

    return inp
