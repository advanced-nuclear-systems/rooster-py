import os

#--------------------------------------------------------------------------------------------------
class Isotope:

    #----------------------------------------------------------------------------------------------
    # constructor: self is an 'isotope' object created in B3
    # isoid is the id of the isotope
    def __init__(self, isoid, reactor):
        # number of energy groups
        ng = reactor.control.input['ng']

        # nuclear data directory
        nddir = reactor.control.input['nddir']
        # isotope is
        self.isoid = isoid

        # open, read, split by eol and close the isotope data file
        f = open(nddir + os.sep + self.isoid, 'r')
        s = f.read().replace('\r\n', '\n').split('\n')
        # remove leading whitespaces
        s = [str.lstrip() for str in s]
        f.close()

        ntemp = int(s[0])
        # base temperatures
        self.temp = [float(s[i]) for i in range(1,ntemp+1)]
        del s[0:int(s[0])+1]
        nsig0 = int(s[0])
        # base sigma-zeros
        self.sig0 = [float(s[i]) for i in range(1,nsig0+1)]
        del s[0:int(s[0])+1]

        # read the specific data from s
        keyword = ''
        list = []
        skip = False
        self.xs = {'inv':[0]*ng, 'chi':[0]*ng, 'nubar':[[0]*ntemp for i in range(ng)], 'abs':[[[0]*nsig0 for i in range(ntemp)] for j in range(ng)], 'n2n':[], 'sca':[], 'fis':[[[0]*nsig0 for i in range(ntemp)] for j in range(ng)], 'tot':[[[0]*nsig0 for i in range(ntemp)] for j in range(ng)], 'tot1':[[[0]*nsig0 for i in range(ntemp)] for j in range(ng)]}
        for i in range(len(s)):
            # read string beginning
            s_beg = s[i][0:3].replace('.','').rstrip()
            # check if the string begins with a text
            if s_beg != '' and not s_beg.isnumeric():
                # keyword
                keyword = s_beg
                temp = float(s[i+1].replace('k',''))
                itemp = self.temp.index(temp)
            try:
                list = [float(e) for e in s[i].split()]
                skip = (list == [])
            except:
                skip = True
            if not skip:
                if keyword == 'tot':
                    if int(list[1]) == 0:
                        self.xs['tot'][int(list[0])-1][itemp] = list[2:]
                    elif int(list[1]) == 1:
                        self.xs['tot1'][int(list[0])-1][itemp] = list[2:]
                elif keyword == 'nga':
                    self.xs['abs'][int(list[0])-1][itemp] += list[2:]
                elif keyword == 'fis':
                    self.xs['abs'][int(list[0])-1][itemp] += list[2:]
                    self.xs['fis'][int(list[0])-1][itemp] = list[2:]
                elif keyword == 'nub':
                    self.xs['nubar'][int(list[0])-1][itemp] = list[1]
                elif keyword == 'dnu':
                    pass
                elif keyword == 'ela':
                    if int(list[2]) == 0:
                        f = int(list[0])-1
                        t = int(list[1])-1
                        self.xs['sca'].append([itemp, (f,t)] + list[3:])
                elif keyword == 'ine':
                    f = int(list[0])-1
                    t = int(list[1])-1
                    self.xs['sca'].append([itemp, (f,t)] + [list[2]]*nsig0)
                elif keyword == 'n2n':
                    f = int(list[0])-1
                    t = int(list[1])-1
                    self.xs['n2n'].append([(f,t), list[2]])
                    pass
                elif keyword == 'chi':
                    self.xs['chi'][int(list[0])-1] = list[1]
                elif keyword == 'chd':
                    pass
                elif keyword == 'inv':
                    self.xs['inv'][int(list[0])-1] = list[1]
                elif keyword == 'nab':
                    for k in range(ntemp):
                        for l in range(nsig0):
                            self.xs['abs'][int(list[0])-1][k][l] += list[1]
                elif keyword == 'xi':
                    pass
                elif keyword == 'hea':
                    pass
                elif keyword == 'the':
                    pass
            
        #if self.isoid == 'U238b6' : print(self.xs['n2n'])
