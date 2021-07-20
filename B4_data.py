#--------------------------------------------------------------------------------------------------
class Data:

    #----------------------------------------------------------------------------------------------
    # constructor: self is a 'data' object created in B
    def __init__(self, reactor):
        pass

    #----------------------------------------------------------------------------------------------
    # material properties: self is a 'data' object created in B
    def matpro(self, type, t):
        # stainless steel type of 316
        if type == 'ss316':
            # @300K equation from Leibowitz, et al, "Properties for LMFBR safety analysis", ANL-CEN-RSD-76-1 (1976), p.117
            rho = 7954.
            # Leibowitz, et al, "Properties for LMFBR safety analysis", ANL-CEN-RSD-76-1 (1976), p.100. Note that 1 mol of SS316 = 10.165 kg (https://www.webqc.org/molecular-weight-of-SS316.html) and 1 cal = 4.184 J
            cp = (6.181 + 1.788e-3*t)*10.165*4.184
            # Leibowitz, et al, "Properties for LMFBR safety analysis", ANL-CEN-RSD-76-1 (1976), p.100.
            k = 9.248 + 1.571e-2*t
            return {'rho':rho, 'cp':cp, 'k':k}
