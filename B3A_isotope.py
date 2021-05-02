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
        self.xs = {'abs':[[[0]*nsig0 for i in range(ntemp)] for j in range(ng)], 'chi':[0]*ng, 'ela':[], 'fis':[[[0]*nsig0 for i in range(ntemp)] for j in range(ng)], 'ine':[], 'inv':[0]*ng, 'n2n':[], 'nub':[[0]*ntemp for i in range(ng)], 'tot':[[[0]*nsig0 for i in range(ntemp)] for j in range(ng)], 'tot1':[[[0]*nsig0 for i in range(ntemp)] for j in range(ng)]}
        # cycle over lines of s
        for i in range(len(s)):
            # read 3 symbols at the line beginning
            s_beg = s[i][0:3].replace('.','').rstrip()
            # check if the line begins with a text
            if s_beg != '' and not s_beg.isnumeric():
                # keyword
                keyword = s_beg
                temp = float(s[i+1].replace('k',''))
                itemp = self.temp.index(temp)
            try:
                # list of values in the line 
                list = [float(e) for e in s[i].split()]
                # if successful set skip to True unless this is the last line
                skip = (list == [])
            except:
                # the line contains non-numerical value: this is temperature line with "k"
                skip = True
            if not skip:
                if keyword == 'tot':
                    ig = int(list[0])-1
                    if int(list[1]) == 0:
                        self.xs['tot'][ig][itemp] = list[2:]
                    elif int(list[1]) == 1:
                        self.xs['tot1'][ig][itemp] = list[2:]
                elif keyword == 'nga':
                    ig = int(list[0])-1
                    for j in range(nsig0):
                        self.xs['abs'][ig][itemp][j] += list[j+1]
                elif keyword == 'fis':
                    ig = int(list[0])-1
                    self.xs['fis'][ig][itemp] = list[1:]
                    for j in range(nsig0):
                        self.xs['abs'][ig][itemp][j] += list[j+1]
                elif keyword == 'nub':
                    ig = int(list[0])-1
                    self.xs['nub'][ig][itemp] = list[1]
                elif keyword == 'dnu':
                    pass
                elif keyword == 'ela':
                    if int(list[2]) == 0:
                        f = int(list[0])-1
                        t = int(list[1])-1
                        self.xs['ela'].append([itemp, (f,t)] + list[3:])
                elif keyword == 'ine':
                    f = int(list[0])-1
                    t = int(list[1])-1
                    self.xs['ine'].append([(f,t), list[2]])
                elif keyword == 'n2n':
                    f = int(list[0])-1
                    t = int(list[1])-1
                    self.xs['n2n'].append([(f,t), list[2]])
                elif keyword == 'chi':
                    ig = int(list[0])-1
                    self.xs['chi'][ig] = list[1]
                elif keyword == 'chd':
                    pass
                elif keyword == 'inv':
                    ig = int(list[0])-1
                    self.xs['inv'][ig] = list[1]
                elif keyword == 'nab':
                    for k in range(ntemp):
                        ig = int(list[0])-1
                        for l in range(nsig0):
                            self.xs['abs'][ig][k][l] += list[1]
                elif keyword == 'xi':
                    pass
                elif keyword == 'hea':
                    pass
                elif keyword == 'the':
                    pass
   
        # number of entries in elastic scaterring matrix
        n = len(self.xs['ela'])
        # number of entries in elastic scattering matrix for one temperature
        n1 = int(n/ntemp)
        # re-arrange s to facilitate temperature interpolation [(f,t) [sig[0], sig[1], ..., sig[nsig0-1]]*ntemp]
        s = []
        for j in range(n1):
            # find and group all entries with (f,t) tuple
            f_t = self.xs['ela'][j][1]
            s.append([f_t])
            for k in range(n):
                if self.xs['ela'][k][1] == f_t:
                    s[j] += [self.xs['ela'][k][2:]]
        self.xs['ela'] = s

        #if self.isoid == 'U238b6' : print(self.xs['n2n'])
