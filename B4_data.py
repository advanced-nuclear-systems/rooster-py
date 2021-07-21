import math

#--------------------------------------------------------------------------------------------------
class Data:

    #----------------------------------------------------------------------------------------------
    # constructor: self is a 'data' object created in B
    def __init__(self, reactor):
        pass

    #----------------------------------------------------------------------------------------------
    # material properties: self is a 'data' object created in B, inp is a dictionary of input data dependent on the material
    def matpro(self, inp):

        # he: helium gas
        if inp['type'] == 'he':
            t = inp['t']
            k = 2.639e-3*t**0.7085
            return {'k':k}

        # mox: mixed uranium-plutonium oxide fuel
        if inp['type'] == 'mox':
            t,b,por,pu,x = inp['t'],inp['b'],inp['por'],inp['pu'],inp['x']
            # density (kg/m3)
            rho = (11460*pu + 10960*(1 - pu)) * (1 - por)
            # specific heat (J/kg-K), D.L. Hagrman, et al., "MATPRO-version 11", TREE-NUREG-1280, Rev 1, Idaho National Engineering Laboratory (1980).
            cp = 15.496*(19.53*539**2 * math.exp(539/t) / (t**2 * (math.exp(539/t) - 1)**2) + 2*9.25e-04*t + 6.02e06*4.01e4 / (1.987*t**2) * math.exp(-4.01e4/(1.987*t)))
            # thermal conductivity (W/m-K), Y. Philipponneau, J. Nuclear Matter., 188 (1992) 194-197
            k = (1/( 1.528*math.sqrt(x+0.00931) - 0.1055 + 0.44*b + 2.855e-4*t ) + 76.38e-12*t**3) * (1-por)/(1+por)/0.864
            return {'rho':rho, 'cp':cp, 'k':k}

        # na: liquid sodium
        elif inp['type'] == 'na':
            t = inp['t']
            # J.K. Fink and L. Leibowitz "Thermodynamic and Transport Properties of Sodium Liquid and Vapor", ANL/RE-95/2, 1995, https://www.ne.anl.gov/eda/ANL-RE-95-2.pdf
            rhol = 219.0 + 275.32*(1.0 - t/2503.7) + 511.58*(1.0 - t/2503.7)**0.5
            visl = math.exp(-6.4406 - 0.3958*math.log(t) + 556.835/t)/rhol
            kl = 124.67 - 0.11381*t + 5.5226e-5*t**2 - 1.1842e-8*t**3
            # Based on fit from J.K. Fink, et. al."Properties for Reactor Safety Analysis", ANL-CEN-RSD-82-2, May 1982.
            cpl = 1646.97 - 0.831587*t + 4.31182e-04*t**2
            return {'rhol':rhol, 'visl':visl, 'kl':kl, 'cpl':cpl}

        # ss316: stainless steel type of 316
        elif inp['type'] == 'ss316':
            t = inp['t']
            # density (kg/m3): @300K equation from Leibowitz, et al, "Properties for LMFBR safety analysis", ANL-CEN-RSD-76-1 (1976), p.117
            rho = 7954.
            # specific heat (J/kg-K): Leibowitz, et al, "Properties for LMFBR safety analysis", ANL-CEN-RSD-76-1 (1976), p.100. Note that 1 mol of SS316 = 10.165 kg (https://www.webqc.org/molecular-weight-of-SS316.html) and 1 cal = 4.184 J
            cp = (6.181 + 1.788e-3*t)*10.165*4.184
            # thermal conductivity (W/m-K): Leibowitz, et al, "Properties for LMFBR safety analysis", ANL-CEN-RSD-76-1 (1976), p.100.
            k = 9.248 + 1.571e-2*t
            return {'rho':rho, 'cp':cp, 'k':k}
