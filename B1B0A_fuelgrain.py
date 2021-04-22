import numpy as np

#--------------------------------------------------------------------------------------------------
class FuelGrain:

    K = 8.6174e-5  # Boltzmann constant in eV/K
    KJ = 1.38e-23  # Boltzmann constant in j/K
    NB = 11  # Number of intragranular bubble classes
    ROU = 41e-24  # Atomic volume of uranium atom in fuel cm-3

    #----------------------------------------------------------------------------------------------
    # constructor: self is a 'fuelgrain' object created in B1B0
    def __init__(self, indx, indxfuel, indxfuelrod, reactor):

        # INPUT PARAMETERS
        # grain diameter
        self.dgrain = reactor.control.input['dgrain']
        # number of nodes in the grain
        self.nr = reactor.control.input['nrgrain']
        # fission rate
        self.frate = reactor.control.input['frate']
        # temperature (temporal)
        self.temp = 300
        # the external pressure (temporal)
        self.pext = 0

        # GEOMETRY
        # mesh grid step
        self.dr = 0.5*self.dgrain/(self.nr-1)
        # list of node radii (size = nr)
        self.r = [i*self.dr for i in range(self.nr)]
        # list of node boundary radii (size = nr-1)
        self.rb = [self.r[i]+self.dr/2 for i in range(self.nr-1)]
        # list of node volume (size = nr)
        self.vol = [4/3 * self.rb[0]**3] + [4/3 * (self.rb[i]**3 - self.rb[i-1]**3) for i in range(1, self.nr-1)] + [4/3 * (self.r[self.nr-1]**3 - self.rb[self.nr-2]**3)]

        # GLOBAL VARIABLES
        # irradiation-induced point defects
        self.tfmin = 4.5/(self.K*np.log(0.3e4/(1.27e-29*self.frate)))
        self.dv = 0.3e4*np.exp(-4.5/(self.K * max(self.temp, self.tfmin)))  # Uranium-vacancy (self-) diffusion coefficient cm2/s
        self.di = 7.1e3*np.exp(-1.9/(self.K * max(self.temp, self.tfmin)))  # Uranium-interstitial diffusion coefficient cm2/s
        self.np = 3  # TODO Number of intragranular as-fabricated pore classes
        self.dsgr = 0.46e-4 # todo Sub-grain size but low limit ds0 was taken here
        self.bsgr = 1 # todo Sub-grain concentration

        # diffusion constant
        self.dg0 = 3.3e4  # cm2/s
        self.dg1 = 2.5e-18  # cm2/s
        self.eg0 = 4.6  # eV
        self.eg1 = 1.19  # eV
        self.dg = self.dg0 * np.exp(-self.eg0/(self.K * self.temp))+self.dg1 * np.sqrt(self.frate) * np.exp(-self.eg1/(self.K * self.temp))

        # intragranular bubble variables

        # resolution variables
        self.ri_b = [1] * self.NB  #initialize of Gas-atom resolution rate from intragranular bubbles
        self.ri_p = [1] * self.NB  #initialize of Interaction rate of fission fragments with intragranular as-fabricated pores
        self.dvalue = 0.46e-4  # todo Current grain size (a variable affected by equi-axed grain growth) but low limit ds0 was taken here
        self.s = 3  #Number of ungrouped bubble classes, treated as ‘solid’ spheres
        self.n = [0] * self.NB # initialize of Bubble-size distribution sampling
        self.fr = [0] * self.NB # initialize of Average fraction of gas atoms in a bubble that may undergo resolution
        self.bpi = [1] * self.NB # todo As-fabricated intragranular pore concentration
        self.kij_bbias = np.zeros((self.NB, self.NB)) # initialize of bubbles coalescence model- random bubble-migration
        self.kij_pbias = np.zeros((self.np, self.np)) # initialize of bubbles coalescence model- biased bubble-motion
        self.v = [1] * self.NB #todo Intragranular bubble drift velocity due to temperature gradient
        self.lamdasgr_b = [4 * np.pi * self.dsgr / 2 * self.bsgr] * self.NB # todo Sub-grain boundaries sink-strength for bubbles  cm-2
        self.ssgr_b = [1] * self.NB  # todo Sub-grain cross-sectional area cm2

        # INITIALIZATION
        # INTRAGRANULAR PROCESSES
        # fractional concentration of irradiation-induced uranium vacancies, cm-3/cm3
        self.cv_irr = [1] * self.NB
        # fractional concentration of irradiation-induced uranium interstitials, cm-3/cm3
        self.ci_irr = [1] * self.NB
        # fractional concentration of uranium vacancies ejected from intragranular as-fabricated pores
        self.cv_p = [1] * self.NB
        # monoatom concentrations
        self.c1 = [0]*self.nr
        # bubble radius
        self.ri = [1] * self.NB
        # intragranular gaseous-bubble concentration
        self.bi = [1] * self.NB

    #----------------------------------------------------------------------------------------------
    # create right-hand side list: self is a 'fuelgrain' object created in B1B0
    def calculate_rhs(self, reactor, t):

        # INTRAGRANULAR PROCESSES: IRRADIATION-INDUCED POINT DEFECTS

        # fractional concentration of irradiation-induced uranium vacancies todo cm-3/cm3
        yiv = 1e4 # Point defects yield constant
        zvirr = 1 # Irradiation-induced vacancy dislocation bias factor
        rhod = max(1e12, 1e4 * np.exp(-2.07e-3 * self.temp + 21.825))
        rd = 1/np.sqrt(np.pi * rhod)
        lamdad = (2 * np.pi * rhod)/np.log(rd/(6e-8)) # Dislocations sink-strength for point defects cm-2
        lamdasgr_v = 4 * np.pi * self.dsgr/2 * self.bsgr # Sub-grain boundaries sink-strength for vacancies cm-2
        lamdabi = [0] * self.NB  # todo Intragranular bubbles sink-strength for point defects cm-2
        lamdapi = [0] * self.np # todo Intragranular as-fabricated sink-strength for point defects cm-2
        etaiv = 1e16 # Vacancy-interstitial recombination constant cm-2
        sumlamdabi = np.sum(lamdabi[1:self.NB])
        sumlamdapi = np.sum(lamdapi[0:self.np])

        dcv_irrdt = [yiv * (self.frate * 1e-6) * self.ROU - (zvirr * lamdad + lamdasgr_v + sumlamdabi+sumlamdapi) * self.dv * self.cv_irr[i]
                     - etaiv * self.di * self.cv_irr[i] * self.ci_irr[i] for i in range(self.NB)]

        # fractional concentration of irradiation-induced uranium interstitials todo cm-3/cm3
        ziirr = 1 # Irradiation-induced interstitial dislocation bias factor
        lamdasgr_i = 4 * np.pi * self.dsgr/2 * self.bsgr  # Sub-grain boundaries sink-strength for vacancies cm-2
        dci_irrdt = [yiv * (self.frate * 1e-6) * self.ROU - (ziirr * lamdad + lamdasgr_i + sumlamdabi + sumlamdapi) * self.di * self.ci_irr[i]
                     - etaiv * self.di * self.cv_irr[i] * self.ci_irr[i] for i in range(self.NB)]

        # fractional concentration of uranium vacancies ejected from intragranular as-fabricated pores todo cm-3/cm3
        zvp = 0 # Dislocation bias factor relating to vacancies emitted by as-fabricated pores
        nv = 1 # todo Average number of vacancies that may undergo resolution
        ri_p = [0] * self.NB #todo Interaction rate of fission fragments with intragranular as-fabricated pores 1/cm3/s
        sumri_p = np.sum(ri_p[0:self.NB])
        dcv_pdt = [self.ROU * nv * sumri_p - (zvp * lamdad + sumlamdapi + lamdasgr_v) * self.dv * self.cv_p[i] for i in range(self.NB)]

        # INTRAGRANULAR PROCESSES: MONOATOMS

        # vector of rb**2 * dg * dc1/dr (size = nr-1)
        q1 = [self.rb[i]**2 * self.dg * (self.c1[i] - self.c1[i+1])/self.dr for i in range(self.nr-1)]
        # vector of time derivative of monoatom concentrations
        dc1dt = [-q1[0]/self.vol[0]] + [(q1[i-1] - q1[i])/self.vol[i] for i in range(1, self.nr-1)] + [q1[self.nr-2]/self.vol[self.nr-1]]
        dc1dt = [dc1dt[i] + 0.31 * self.frate for i in range(self.nr)]

        # INTRAGRANULAR PROCESSES: RESOLUTION

        # Diffusion-controlled bubble growth ri
        self.ciu = [10] * self.NB  # todo to be replaced
        self.cvu = [10] * self.NB  # todo to be replaced
        gama = 0.63
        rho = 41e-24  # cm-3
        self.pi = [10] * self.NB  # todo to be replaced
        deltapi = [self.pi[i] - 2 * gama / (self.ri[i] - self.pext) for i in range(self.NB)]
        deltaci = [self.ci_irr[i] + self.ciu[i] * (1 - np.exp(deltapi[i] * rho / (self.KJ * self.temp))) for i in range(self.NB)]
        deltacv = [self.cv_irr[i] + self.cvu[i] * (1 - np.exp(deltapi[i] * rho / (self.KJ * self.temp))) for i in range(self.NB)]
        dridt = [(self.dv * deltacv[i] - self.di * deltaci[i]) / self.ri[i] for i in range(self.NB)]

        # Resolution models
        for i in range(self.NB):
            if i < self.s-1:
                self.n[i] = i + 2
            else:
                self.n[i] = 5 * self.n[i - 1]
        for i in range(self.NB):
            if self.ri[i] <= self.dvalue:
                self.fr[i] = 1
            else:
                self.fr[i] = (self.ri[i]**3 -(self.ri[i]-self.dvalue)**3) * self.ri[i] ** -3
        self.ri_b = [3.6e-17*self.frate * self.fr[i] *self.n[i] * self.bi[i] for i in range(self.NB)]
        self.ri_p = [2 * np.pi *(self.ri[i]+(1e-7))**2 * (6e-4) * self.frate * self.bpi[i] for i in range(self.NB)]

        # Intragranular gaseous-bubble concentration
        self.k11 = 2e-6 * self.dg * self.bi[0] * self.bi[0]
        for i in range(self.NB):
            for j in range(self.NB):
                self.kij_bbias[i,j] = np.pi * (self.ri[i]+self.ri[j])**2 * abs(self.v[i]-self.v[j]) * self.bi[i] * self.bi[j]
        for i in range(self.np):
            for j in range(self.np):
                self.kij_pbias[i,j] = np.pi * (self.ri[i]+self.ri[j])**2 * abs(self.v[i]-self.v[j]) * self.bi[i] * self.bi[j]
        sumkb = np.sum(self.kij_bbias[1,1:self.NB]) # todo to be completed
        sumkp = np.sum(self.kij_pbias[1,0:self.np]) # todo to be completed
        self.d = [1] * self.NB # todo to be replaced
        db2dt = self.k11 + self.ri_b[2] - self.ri_b[1] - sumkb - sumkp - self.lamdasgr_b[1] * self.d[1] * self.bi[1] - self.ssgr_b[1] * self.v[1] * self.ssgr_b[1]
        dbidt = [1] + [db2dt] + [1]*9

        dc1dt = [dc1dt[i] + 2 * self.ri_b[1] for i in range(self.nr)]

        rhs = dc1dt + dridt + dcv_irrdt + dci_irrdt + dcv_pdt + dbidt
        return rhs
