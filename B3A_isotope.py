import os

#--------------------------------------------------------------------------------------------------
class Isotope:

    #----------------------------------------------------------------------------------------------
    # constructor: self is an 'isotope' object created in B3
    # isoid is the id of the isotope
    def __init__(self, isoid, reactor):

        # nuclear data directory
        nddir = reactor.control.input['nddir']
        # isotope is
        self.isoid = isoid

        # open, read line by line and close the isotope data file
        f = open(nddir + os.sep + self.isoid, 'r')
        s = f.readline()
        cards = []
        end = False
        while not end:
            s = f.readline()
            w = []
            # six words of 11-character length
            for i in range(6):
                w.append(s[i*11:(i+1)*11].lstrip())
                wi = w[i]
                # add "E" to floats
                if '.' in wi:
                    if '+' in w[i]:
                        w[i] = float(w[i].replace('+','E+'))
                    elif not w[i].startswith('-') and '-' in w[i]:
                        w[i] = float(w[i].replace('-','E-'))
                    elif w[i].startswith('-') and '-' in wi[1:]:
                        w[i] = -float(wi[1:].replace('-','E-'))
                elif w[i] != '':
                    w[i] = int(w[i])
            # append three last values (mf, mt and line number)
            cards.append(w + [(int(s[70:72]), int(s[72:75])), int(s[75:80])])
            if int(s[66:70]) == -1 : end = True
        f.close()
        
        # number of energy groups
        ng = cards[1][2]

        # find number of base temperatures and values of base temperatures using mf=1 and mt=451
        ntemp = 0
        self.temp = []
        for card in cards:
            if card[6] == (1,451) and card[7] == 2 : 
                ntemp += 1
                self.temp.append(card[0])
        
        # find number of base sigma-zeros 
        nsig0 = cards[0][3]
        a, irownew = extract_n_words(8, 2, cards)
        # base sigma-zeros
        self.sig0 = a[1:]

        # dictionary of cross sections
        self.xs = {'abs':[[[0]*nsig0 for i in range(ntemp)] for j in range(ng)],  'cap':[[[0]*nsig0 for i in range(ntemp)] for j in range(ng)], 'chi':[0]*ng, 'ela':[], 'fis':[[[0]*nsig0 for i in range(ntemp)] for j in range(ng)], 'ine':[], 'inv':[0]*ng, 'n2n':[], 'nub':[[0]*ntemp for i in range(ng)], 'tot':[[[0]*nsig0 for i in range(ntemp)] for j in range(ng)], 'tot1':[[[0]*nsig0 for i in range(ntemp)] for j in range(ng)]}

        # fission spectrum (mf = 5, mt = 18)
        chi = extract_mf_mt(5, 18, 0, 0, cards)
        if chi != [] : self.xs['chi'] = chi
        # inverted neutron velocity (mf = 3, mt = 259)
        inv = extract_mf_mt(3, 259, 0, 0, cards)
        for ig in range(ng):
            self.xs['inv'][ig] = inv[ig][0]
        for itemp in range(ntemp):
            # first Legendre component of total xs (mf = 3, mt = 1)
            nlgndr = 1
            sigt1 = extract_mf_mt(3, 1, itemp, nlgndr, cards)
            for ig in range(ng):
                self.xs['tot1'][ig][itemp] = sigt1[ig]

            nlgndr = 0
            # total xs (mf = 3, mt = 1)
            sigt = extract_mf_mt(3, 1, itemp, nlgndr, cards)
            # nubar (mf = 3, mt = 452)
            nubar = extract_mf_mt(3, 452, itemp, nlgndr, cards)
            # fission xs (mf = 3, mt = 18)
            sigf = extract_mf_mt(3, 18, itemp, nlgndr, cards)
            # capture xs (mf = 3, mt = 102)
            sigc = extract_mf_mt(3, 102, itemp, nlgndr, cards)
            for ig in range(ng):
                for i in range(nsig0):
                    self.xs['tot'][ig][itemp][i] = sigt[ig][i]
                    self.xs['cap'][ig][itemp][i] = sigc[ig][i]
                    self.xs['abs'][ig][itemp][i] = sigc[ig][i]
                if nubar != [] : 
                    self.xs['nub'][ig][itemp] = nubar[ig][0]
                    self.xs['fis'][ig][itemp] = sigf[ig]
                    for i in range(nsig0):
                        self.xs['abs'][ig][itemp][i] += sigf[ig][i]
            # elastic scattering (mt = 2)
            sige = extract_mf6(2, itemp, nlgndr, cards)
            for s in sige:
                self.xs['ela'].append(s)
            #print(self.xs['cap'][0][0])
        # number of entries in elastic scattering matrix
        n = len(self.xs['ela'])
        # number of entries in elastic scattering matrix for one temperature
        n1 = int(n/ntemp)
        # re-arrange self.xs['ela'] to facilitate temperature interpolation [(f,t) [sig[0], sig[1], ..., sig[nsig0-1]]*ntemp]
        s = []
        for j in range(n1):
            # find and group all entries with (f,t) tuple
            f_t = self.xs['ela'][j][0]
            s.append([f_t])
            for k in range(n):
                if self.xs['ela'][k][0] == f_t:
                    s[j] += [self.xs['ela'][k][1:]]
        self.xs['ela'] = s

        # inelastic scattering (mt = 51... 91)
        for mt in range(51,92):
            sigi = extract_mf6(mt, 0, 0, cards)
            if sigi != []:
                for i in range(len(sigi)):
                    self.xs['ine'].append(sigi[i])

#        # open, read, split by eol and close the isotope data file
#        f = open(nddir + os.sep + self.isoid, 'r')
#        s = f.read().replace('\r\n', '\n').split('\n')
#        # remove leading whitespaces
#        s = [str.lstrip() for str in s]
#        f.close()
#
#        ntemp = int(s[0])
#        # base temperatures
#        self.temp = [float(s[i]) for i in range(1,ntemp+1)]
#        del s[0:int(s[0])+1]
#        nsig0 = int(s[0])
#        # base sigma-zeros
#        self.sig0 = [float(s[i]) for i in range(1,nsig0+1)]
#        del s[0:int(s[0])+1]
#
#        # read the specific data from s
#        keyword = ''
#        list = []
#        skip = False
#        self.xs = {'abs':[[[0]*nsig0 for i in range(ntemp)] for j in range(ng)], 'chi':[0]*ng, 'ela':[], 'fis':[[[0]*nsig0 for i in range(ntemp)] for j in range(ng)], 'ine':[], 'inv':[0]*ng, 'n2n':[], 'nub':[[0]*ntemp for i in range(ng)], 'tot':[[[0]*nsig0 for i in range(ntemp)] for j in range(ng)], 'tot1':[[[0]*nsig0 for i in range(ntemp)] for j in range(ng)]}
#        # cycle over lines of s
#        for i in range(len(s)):
#            # read 3 symbols at the line beginning
#            s_beg = s[i][0:3].replace('.','').rstrip()
#            # check if the line begins with a text
#            if s_beg != '' and not s_beg.isnumeric():
#                # keyword
#                keyword = s_beg
#                temp = float(s[i+1].replace('k',''))
#                itemp = self.temp.index(temp)
#            try:
#                # list of values in the line 
#                list = [float(e) for e in s[i].split()]
#                # if successful set skip to True unless this is the last line
#                skip = (list == [])
#            except:
#                # the line contains non-numerical value: this is temperature line with "k"
#                skip = True
#            if not skip:
#                if keyword == 'tot':
#                    ig = int(list[0])-1
#                    if int(list[1]) == 0:
#                        self.xs['tot'][ig][itemp] = list[2:]
#                    elif int(list[1]) == 1:
#                        self.xs['tot1'][ig][itemp] = list[2:]
#                elif keyword == 'nga':
#                    ig = int(list[0])-1
#                    for j in range(nsig0):
#                        self.xs['abs'][ig][itemp][j] += list[j+1]
#                elif keyword == 'fis':
#                    ig = int(list[0])-1
#                    self.xs['fis'][ig][itemp] = list[1:]
#                    for j in range(nsig0):
#                        self.xs['abs'][ig][itemp][j] += list[j+1]
#                elif keyword == 'nub':
#                    ig = int(list[0])-1
#                    self.xs['nub'][ig][itemp] = list[1]
#                elif keyword == 'dnu':
#                    pass
#                elif keyword == 'ela':
#                    if int(list[2]) == 0:
#                        f = int(list[0])-1
#                        t = int(list[1])-1
#                        self.xs['ela'].append([itemp, (f,t)] + list[3:])
#                elif keyword == 'ine':
#                    f = int(list[0])-1
#                    t = int(list[1])-1
#                    self.xs['ine'].append([(f,t), list[2]])
#                elif keyword == 'n2n':
#                    f = int(list[0])-1
#                    t = int(list[1])-1
#                    self.xs['n2n'].append([(f,t), list[2]])
#                elif keyword == 'chi':
#                    ig = int(list[0])-1
#                    self.xs['chi'][ig] = list[1]
#                elif keyword == 'chd':
#                    pass
#                elif keyword == 'inv':
#                    ig = int(list[0])-1
#                    self.xs['inv'][ig] = list[1]
#                elif keyword == 'nab':
#                    for k in range(ntemp):
#                        ig = int(list[0])-1
#                        for l in range(nsig0):
#                            self.xs['abs'][ig][k][l] += list[1]
#                elif keyword == 'xi':
#                    pass
#                elif keyword == 'hea':
#                    pass
#                elif keyword == 'the':
#                    pass
#   
#        # number of entries in elastic scaterring matrix
#        n = len(self.xs['ela'])
#        # number of entries in elastic scattering matrix for one temperature
#        n1 = int(n/ntemp)
#        # re-arrange s to facilitate temperature interpolation [(f,t) [sig[0], sig[1], ..., sig[nsig0-1]]*ntemp]
#        s = []
#        for j in range(n1):
#            # find and group all entries with (f,t) tuple
#            f_t = self.xs['ela'][j][1]
#            s.append([f_t])
#            for k in range(n):
#                if self.xs['ela'][k][1] == f_t:
#                    s[j] += [self.xs['ela'][k][2:]]
#        self.xs['ela'] = s
#
#        #if self.isoid == 'U238b6' : print(self.xs['n2n'])

#----------------------------------------------------------------------------------------------
# The function reads n words from row irow of matrix cards and returns them in
# vector a together with the new row number irownew, i.e. the row where the last word was read.

def extract_n_words(n, irow, cards):
    irownew = irow
    a = []
    # read lines with 6 words each
    for ii in range(int(n/6)):
        for jj in range(6):
            a.append(cards[irownew][jj])
        irownew += 1
    
    # number of remaining words in the last line
    nremain = n - int(n/6)*6
    if nremain == 0:
        irownew -= 1
    else:
        # read the last line with less than 6 words
        for jj in range(nremain):
            a.append(cards[irownew][jj])
    
    return a, irownew

#----------------------------------------------------------------------------------------------
# The function reads cross sections from file 6 for reaction mt and 
# temperature index itemp from matrix cards and returns the 2D matrix 
# sig{nLgn,nSig0}(nonz) with two vectors ifrom(nonz) and ito(nonz), where
# nLgn is the number of Legendre components, nSig0 is the number of
# sigma-zeros and nonz is the number of nonzeros.
#[(f,t) [sig[0], sig[1], ..., sig[nsig0-1]]*ntemp]

def extract_mf6(mt, itemp, nl, cards):

    # row number
    irow = 0
    # temperature index
    ntemp = -1
    # list to return
    sig = []
    while irow < len(cards):
        # find the row with mf=6 & mt
        if cards[irow][6] == (6,mt):
            # if this is the first line of mf=6 & mt, initialize
            if cards[irow][7] == 1:
                # number of nonzeros
                nonz = 0
                # number of Legendre components
                nlgn = cards[irow][2]
                # number of sigma-zeros
                nsig0 = cards[irow][3]
                # temperature index
                ntemp += 1
                
                irow += 1
            # number of secondary positions
            ng2 = cards[irow][2]
            # index to lowest nonzero group
            ig2lo = cards[irow][3]
            # number of words to be read
            nw = cards[irow][4]
            # current group index
            ig = cards[irow][5]
            
            irow += 1
            # extract nw words in vector a
            a, irownew = extract_n_words(nw, irow, cards)
            irow = irownew

            if ntemp == itemp:
                # the first nlgn*nsig0 words are flux -- skip.
                a = a[nlgn*nsig0:]
                ng2 -= 1
                k = -1
                for ito in range(ig2lo,ig2lo+ng2):
                   nonz += 1
                   s = []
                   # append a tuple (from,to). note: group numbering starts from 0 that's why minus one
                   s.append((ig-1,ito-1))
                   for isig0 in range(nsig0):
                       for ilgn in range(nlgn):
                            k = k + 1;
                            if ilgn == nl:
                                s.append(a[k])
                   sig.append(s)
        irow += 1
    return sig
    
#----------------------------------------------------------------------------------------------
# The function searches matrix cards for cross sections sig from file mf for
# reaction mt, temperature index itemp, legendre order nl and returns out[ng][nsig0], where ng
# is the number of energy groups and nsig0 is the number of sigma-zeros.

def extract_mf_mt(mf, mt, itemp, nl, cards):

    # find index irow of the row with required mf and mt
    ntemp = -1
    irow = 0
    out = []
    while ntemp < itemp and irow < len(cards):
        if cards[irow][6] == (mf,mt) and cards[irow][7] == 1: ntemp += 1
        irow += 1
    # if required (mf, mt) found
    if ntemp > -1:
        # number of sigma-zeros
        nsig0 = cards[irow-1][3]
        # number of Legendre components
        nlgn = cards[irow-1][2]
        if nl > nlgn-1:
            print('****ERROR: legendre order requested in extract_mf_mt is higher than available in the library.')
            sys.exit()
        irow += 1
        if mf == 3:
            while cards[irow][6] == (mf,mt):
                ig = cards[irow-1][5]
                a, irownew = extract_n_words(nsig0*nlgn*2, irow, cards)
                # the first nlgn*nsig0 words are flux -- skip.
                a = a[nsig0*nlgn:]
                out.append(a[0::nlgn])
                irow = irownew + 2
        elif mf == 5:
            ng = cards[irow-1][5]
            out, irownew = extract_n_words(ng, irow, cards)

    return out
