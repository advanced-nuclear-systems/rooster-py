from scipy.interpolate import RegularGridInterpolator, interp1d, LinearNDInterpolator # psi_ytchen
import math
import CoolProp.CoolProp as CP
import numpy as np

#--------------------------------------------------------------------------------------------------
class Data:
    # Here is the testing change at 20240621
    #----------------------------------------------------------------------------------------------
    # constructor: self is a 'data' object created in B
    def __init__(self, reactor):
        for matdic in reactor.control.input['mat']:
            if matdic['type'] == 'h2o':
                # psi_ytchen: creat interpolation functions for fluid properties
                print('Constructing IF97::Water database...')
                fluid = 'IF97::Water'
                hmax = CP.PropsSI('H','T',1073.0,'P',100.0E6, fluid)
                hmin = CP.PropsSI('H','T',300.00,'P', 100.0E3, fluid)
                pmax = 22.00E6
                pmin = 100.0E3
                pnp = np.linspace(pmin,pmax,num=100)
                hnp = np.linspace(hmin,hmax, num=800)
                xnp = np.zeros( len(pnp)*len(hnp) )
                rnp = np.zeros( len(pnp)*len(hnp) )
                vnp = np.zeros( len(pnp)*len(hnp) )
                mnp = np.zeros( len(pnp)*len(hnp) )
                knp = np.zeros( len(pnp)*len(hnp) )
                cnp = np.zeros( len(pnp)*len(hnp) )
                tnp = np.zeros( len(pnp)*len(hnp) )
                px =  np.zeros( len(pnp)*len(hnp) )
                hy =  np.zeros( len(pnp)*len(hnp) )
                value=np.zeros([len(pnp)*len(hnp),7] )
                
                rls = np.zeros( len(pnp) )
                rgs = np.zeros( len(pnp) )
                mls = np.zeros( len(pnp) )
                mgs = np.zeros( len(pnp) )
                cls = np.zeros( len(pnp) )
                cgs = np.zeros( len(pnp) )
                kls = np.zeros( len(pnp) )
                kgs = np.zeros( len(pnp) )
                ts  = np.zeros( len(pnp) ) 
                hgl = np.zeros( len(pnp) ) 
                sgm = np.zeros( len(pnp) ) 
                val1d = np.zeros( [11, len(pnp)] ) 
                
                for i in range(len(pnp)):
                    rls[i] = CP.PropsSI('D','Q',0.0,'P', pnp[i], fluid)
                    rgs[i] = CP.PropsSI('D','Q',1.0,'P', pnp[i], fluid)
                    mls[i] = CP.PropsSI('V','Q',0.0,'P', pnp[i], fluid)
                    mgs[i] = CP.PropsSI('V','Q',1.0,'P', pnp[i], fluid)
                    cls[i] = CP.PropsSI('C','Q',0.0,'P', pnp[i], fluid)
                    cgs[i] = CP.PropsSI('C','Q',1.0,'P', pnp[i], fluid)
                    kls[i] = CP.PropsSI('L','Q',0.0,'P', pnp[i], fluid)
                    kgs[i] = CP.PropsSI('L','Q',1.0,'P', pnp[i], fluid)
                    ts[i]  = CP.PropsSI('T','Q',0.5,'P', pnp[i], fluid)
                    hls    = CP.PropsSI('H','Q',0.0,'P',pnp[i], fluid)
                    hgs    = CP.PropsSI('H','Q',1.0,'P',pnp[i], fluid)
                    hgl[i] = hgs - hls
                    sgm[i] = CP.PropsSI('I','Q',0.0,'P',pnp[i], fluid)
                    for j in range(len(hnp)):
                        ij = i*len(hnp) + j
                        xnp[ij] = (hnp[j]-hls)/hgl[i]
                        tnp[ij] = CP.PropsSI('T','P',pnp[i],'H',hnp[j], fluid)
                        if xnp[ij] < 0. or xnp[ij] > 1. :
                            try:
                                rnp[ij] = CP.PropsSI('D','P',pnp[i],'H',hnp[j], fluid)
                                mnp[ij] = CP.PropsSI('V','P',pnp[i],'H',hnp[j], fluid)
                                knp[ij] = CP.PropsSI('L','P',pnp[i],'H',hnp[j], fluid)
                                cnp[ij] = CP.PropsSI('C','P',pnp[i],'H',hnp[j], fluid)
                            except ValueError:
                                if xnp[ij] < 0.:
                                    rnp[ij] = rls[i]
                                    mnp[ij] = mls[i]
                                    knp[ij] = kls[i]
                                    cnp[ij] = cls[i]
                                else:
                                    rnp[ij] = rgs[i]
                                    mnp[ij] = mgs[i]
                                    knp[ij] = kgs[i]
                                    cnp[ij] = cgs[i]
                        else:
                            rnp[ij] = 1./(1./rls[i] + xnp[ij]*(1./rgs[i] - 1./rls[i]))
                            mnp[ij] = 1.0/( 1.0/mls[i] + xnp[ij]*(1.0/mgs[i] - 1.0/mls[i]) )
                            knp[ij] = kls[i]*(1. - xnp[ij]) + kgs[i]*xnp[ij]
                            cnp[ij] = cls[i]*(1. - xnp[ij]) + cgs[i]*xnp[ij]
                            
                        vnp[ij] = mnp[ij]/rnp[ij]
                        px[ij]  = pnp[i]
                        hy[ij]  = hnp[j]
                
                value[:,0] = rnp.copy()
                value[:,1] = vnp.copy()
                value[:,2] = mnp.copy()
                value[:,3] = knp.copy()
                value[:,4] = cnp.copy()
                value[:,5] = tnp.copy()
                value[:,6] = xnp.copy()
                fph_h2o = LinearNDInterpolator(list( zip(px, hy) ), value)
                self.fph_h2o = fph_h2o
                hpt_h2o = LinearNDInterpolator(list( zip(px, tnp) ), hy)
                self.hpt_h2o = hpt_h2o
                
                val1d[0 , :] = hgl
                val1d[1 , :] = mls
                val1d[2 , :] = mgs
                val1d[3 , :] = rls
                val1d[4 , :] = rgs
                val1d[5 , :] = kls
                val1d[6 , :] = kgs
                val1d[7 , :] = ts 
                val1d[8 , :] = cls
                val1d[9 , :] = cgs
                val1d[10, :] = sgm
                
                fp_h2o = interp1d(pnp, val1d)
                
                self.fp_h2o = fp_h2o
                
                print('IF97::Water database done.')
                
            if matdic['type'] == 'lbe':
                print('Constructing Lead-Bismuth Eutectic database...')
                tm0 = 398.00
                tmin = 3.0E2
                tmax = 2.1E3 # temperatures [K]
                tnp = np.linspace(tmin,tmax,num=100)
                hnp = 164.8*(tnp-tm0)-1.97e-2*(tnp**2-tm0**2)+4.167e-6*(tnp**3-tm0**3)+4.56e5*(1.0/tnp-1.0/tm0)
                rnp = 11065-1.293*tnp
                vnp = (4.94e-4*np.exp(754.1/tnp))/rnp
                cnp = 164.8-3.94e-2*tnp+1.25e-5*tnp*tnp-4.56e5/tnp/tnp
                knp = 3.284 + 1.617e-2*tnp-2.305e-6*tnp*tnp
                values = np.zeros([5, 100])
                values[0,:] = tnp.copy()
                values[1,:] = rnp.copy()
                values[2,:] = vnp.copy()
                values[3,:] = cnp.copy()
                values[4,:] = knp.copy()
                fh_lbe = interp1d(hnp, values)
                self.fh_lbe = fh_lbe
                ht_lbe = interp1d(tnp, hnp)
                self.ht_lbe = ht_lbe
                print('Lead-Bismuth Eutectic database done.')
                
            if matdic['type'] == 'na':
                # psi_ytchen: creat interpolation functions for sodium two-phase properties
                print('Constructing Sodium database...')
                rcrit = 219.0  # psi_ytchen: critical density
                t1  = 371.00   # psi_ytchen: melting temperature, also melting point
                t2  = 2000.0   # psi_ytchen: medium temperature 
                t3  = 2503.70  # psi_ytchen: maximum temperature, also critical temperature
                pmin= 1.000e1  # psi_ytchen: minimum pressure
                pmax= 25.00e6  # psi_ytchen: maximum pressure, also critical pressure
                hmin= 210.0e3  # psi_ytchen: minimum enthalpy
                hmax= 5800.e3  # psi_ytchen: maximum enthalpy
                # psi_ytchen: generate temporary psat-tsat table
                tsat_temp = np.linspace(t1, t3, num = 1000, endpoint = False) 
                psat_temp = 1.0e6*np.exp(11.9463-12633.73/tsat_temp)/tsat_temp**(0.4672)
                tsps_na = interp1d(psat_temp, tsat_temp,fill_value = 'extrapolate') #function
                # psi_ytchen: generate temporary tl-hl table
                hl_temp = np.zeros( len(tsat_temp) )
                for i in range( len(tsat_temp) ):
                   if tsat_temp[i] < t2:
                      hl_temp[i] = -3.6577e5 + 1.6582e3*tsat_temp[i] - 4.2375e-1*tsat_temp[i]**2 + 1.4847e-4*tsat_temp[i]**3 + 2.9926e6/tsat_temp[i]
                   else:
                      hgl_temp = (393.37*(1.0 - tsat_temp[i]/t3) + 4398.6*(1.0 - tsat_temp[i]/t3)**(0.29302))*1000.0
                      hl_temp[i] = 2128.4e3 + 864.96*tsat_temp[i] - 0.5*hgl_temp
                tlhl_na = interp1d(hl_temp, tsat_temp, kind='linear',fill_value = 'extrapolate') #function


                # psi_ytchen: generate temporary cpgs-tsat table
                tsat_temp = np.array([  400.,   500.,   600.,   700.,   800.,   900.,  1000.,  1100.,  1200.,  1300.,  1400.,  1500.,  1600.,  1700.,  1800.,  1900.,  2000.,  2100.,  2200.,  2300.,  2400.,  2500.])
                cpgs_temp = np.array([0.86e3, 1.58e3, 2.02e3, 2.26e3, 2.44e3, 2.56e3, 2.64e3, 2.66e3, 2.64e3, 2.58e3, 2.48e3, 2.35e3, 2.17e3, 2.05e3, 2.06e3, 2.10e3, 2.25e3, 2.52e3, 3.00e3, 4.02e3, 7.70e3,15.06e3])
                cpgs_na = interp1d(tsat_temp, cpgs_temp, kind='linear',fill_value = 'extrapolate') #function
                # psi_ytchen: generate temporary kgsat-tsat table
                tsat_temp = np.array([   373.,   700.,   800.,   900.,  1000.,  1100.,  1200.,  1300.,  1400.,  1500.,  2503.7])
                kgsat_temp= np.array([0.00687, 0.0323, 0.0378, 0.0423, 0.0455, 0.0473, 0.0485, 0.0493, 0.0497, 0.0501,  0.0520])
                kgsat_na = interp1d(tsat_temp, kgsat_temp, kind='linear',fill_value = 'extrapolate') #function
                #
                tsat_temp = np.array([   400.,   500.,   600.,   700.,   800.,   900.,  1000.,  1100.,  1200.,  1300.,  1400.,  1500.,  1600.,  1700.,  1800.,  1900.,  2000.,  2100.,  2200.,  2300.,  2400.,  2500.])
                apgs_temp = np.array([2.55e-3,2.33e-3,2.01e-3,1.85e-3,1.73e-3,1.64e-3,1.57e-3,1.50e-3,1.44e-3,1.38e-3,1.33e-3,1.26e-3,1.19e-3,1.15e-3,1.15e-3,1.19e-3,1.28e-3,1.44e-3,1.76e-3,2.46e-3,4.87e-3,3.74e-1])
                apgs_na = interp1d(tsat_temp, apgs_temp, kind='linear',fill_value = 'extrapolate') #function


                pnp  = np.linspace(pmin, pmax, num = 800, endpoint = False) # psi_ytchen: mesh points for pressure
                hnp  = np.linspace(hmin, hmax, num = 800, endpoint = False) # psi_ytchen: mesh points for enthalpy
                xnp  = np.zeros( len(pnp)*len(hnp) )
                rnp  = np.zeros( len(pnp)*len(hnp) )
                vnp  = np.zeros( len(pnp)*len(hnp) )
                mnp  = np.zeros( len(pnp)*len(hnp) )
                knp  = np.zeros( len(pnp)*len(hnp) )
                cnp  = np.zeros( len(pnp)*len(hnp) )
                tnp  = np.zeros( len(pnp)*len(hnp) )
                px   = np.zeros( len(pnp)*len(hnp) )
                hy   = np.zeros( len(pnp)*len(hnp) )
                value= np.zeros([len(pnp)*len(hnp), 7] )
 
                rls = np.zeros( len(pnp) ) # psi_ytchen:
                rgs = np.zeros( len(pnp) ) # psi_ytchen:
                mls = np.zeros( len(pnp) ) # psi_ytchen:
                mgs = np.zeros( len(pnp) ) # psi_ytchen:
                cls = np.zeros( len(pnp) ) # psi_ytchen:
                cgs = np.zeros( len(pnp) ) # psi_ytchen:
                kls = np.zeros( len(pnp) ) # psi_ytchen:
                kgs = np.zeros( len(pnp) ) # psi_ytchen:
                ts  = np.zeros( len(pnp) ) # psi_ytchen:
                hgl = np.zeros( len(pnp) ) # psi_ytchen:
                sgm = np.zeros( len(pnp) )
                val1d = np.zeros( [11, len(pnp)] )

                for i in range( len(pnp) ):
                    ts[i]  = tsps_na( pnp[i] )      # psi_ytchen: get saturation temperature
                    rls[i] = 219.0 + 275.32*(1.0 - ts[i]/t3) + 511.58*(1.0 - ts[i]/t3)**0.5            # psi_ytchen: IAEA_TECDOC_NAPRO, page 184, eq.47
                    gama   = (12633.73/ts[i]**(2.0)-0.4672/ts[i])*pnp[i]                               # psi_ytchen: IAEA_TECDOC_NAPRO, page 103, eq.13
                    hgl[i] = (393.37*(1.0 - ts[i]/t3) + 4398.6*(1.0 - ts[i]/t3)**(0.29302))*1000.0     # psi_ytchen: IAEA_TECDOC_NAPRO, page 110, eq.14
                    rgs[i] = (hgl[i]/ts[i]/gama + 1.0/rls[i])**(-1.0)                                  # psi_ytchen: IAEA_TECDOC_NAPRO, page 188, eq.48
                    mls[i] = math.exp(-6.4406 - 0.3958*np.log(ts[i]) + 556.835/ts[i])                  # psi_ytchen: IAEA_TECDOC_NAPRO, page 395, eq.166
                    miugd = 1.2375e-05 + 4.4828e-09*ts[i]                                              # psi_ytchen: IAEA_TECDOC_NAPRO, page 409, eq.173
                    miugdc= 1.2375e-05 + 4.4828e-09*t3                                                 # psi_ytchen: IAEA_TECDOC_NAPRO, page 409, eq.173
                    mgs[i] = miugd + (5.8e-05 - miugdc)*(ts[i]/t3)*(rgs[i]/rcrit)                      # psi_ytchen: IAEA_TECDOC_NAPRO, page 409, eq.172
                    cls[i] = 1658.2 - 0.84750*ts[i] + 4.4541e-04*ts[i]**2 - 2.9926e6/ts[i]**2          # psi_ytchen: IAEA_TECDOC_NAPRO, page 270, eq.97
                    cgs[i] = cpgs_na( ts[i] )                                                          # psi_ytchen: IAEA_TECDOC_NAPRO, page 315, tab. 246
                    kls[i] = 124.67 - 0.11381*ts[i] + 5.5226e-5*ts[i]**2 - 1.1842e-8*ts[i]**3          # psi_ytchen: IAEA_TECDOC_NAPRO, page 346, eq.156
                    kgs[i] = kgsat_na( ts[i] )                                                         # psi_ytchen: IAEA_TECDOC_NAPRO, page 360, tab. 269
                    sgm[i] = 240.5e-3*( (t3-ts[i])/t3 )**(1.126)                                       # psi_ytchen: IAEA_TECDOC_NAPRO, page 133, eq. 21
                    if ts[i] < t2:                                                                     # psi_ytchen: IAEA_TECDOC_NAPRO, page 270, eq.97
                       hls = -3.6577e5 + 1.6582e3*ts[i] - 4.2375e-1*ts[i]**2 + 1.4847e-4*ts[i]**3 + 2.9926e6/ts[i]
                    else:
                       hls = 2128.4e3 + 864.96*ts[i] - 0.5*hgl[i]
                    hgs = hls + hgl[i]
                    for j in range( len(hnp)):
                        ij = i*len(hnp) + j
                        xnp[ij] = (hnp[j]-hls)/hgl[i]
                        if xnp[ij] < 0.0: # psi_ytchen: subcooled liquid na
                           tnp[ij] = tlhl_na(hnp[j]) # calculate liquid temperature with given h
                           rnp[ij] = 219.0 + 275.32*(1.0 - tnp[ij]/t3) + 511.58*(1.0 - tnp[ij]/t3)**0.5
                           mnp[ij] = math.exp(-6.4406 - 0.3958*np.log(tnp[ij]) + 556.835/tnp[ij])
                           knp[ij] = 124.67 - 0.11381*tnp[ij] + 5.5226e-5*tnp[ij]**2 - 1.1842e-8*tnp[ij]**3
                           cnp[ij] = 1658.2 - 0.84750*tnp[ij] + 4.4541e-04*tnp[ij]**2 - 2.9926e6/tnp[ij]**2
                        elif xnp[ij] < 1.0: # psi_ytchen: saturated two-phase na
                           tnp[ij] = ts[i]
                           rnp[ij] = rgs[i]*rls[i]/(xnp[ij]*rls[i] + (1.-xnp[ij])*rgs[i])
                           mnp[ij] = 1.0/( 1.0/mls[i] + xnp[ij]*(1.0/mgs[i] - 1.0/mls[i]) )
                           knp[ij] = kls[i]*(1. - xnp[ij]) + kgs[i]*xnp[ij]
                           cnp[ij] = cls[i]*(1. - xnp[ij]) + cgs[i]*xnp[ij]
                        else: # psi_ytchen: superheated vapor na
                           tnp[ij] = ts[i] + (hnp[j]-hgs)/cgs[i]
                           rnp[ij] = rgs[i] - apgs_na(ts[i])*rgs[i]*(tnp[ij]-ts[i])
                           miugd_ij = 1.2375e-05 + 4.4828e-09*tnp[ij]
                           mnp[ij] = miugd_ij + (5.8e-05 - miugdc)*(tnp[ij]/t3)*(rnp[ij]/rcrit)
                           knp[ij] = kgsat_na( tnp[ij] )
                           cnp[ij] = cpgs_na( tnp[ij] )
                        vnp[ij] = mnp[ij]/rnp[ij]
                        px[ij] = pnp[i]
                        hy[ij] = hnp[j]

                value[:,0] = rnp.copy()
                value[:,1] = vnp.copy()
                value[:,2] = mnp.copy()
                value[:,3] = knp.copy()
                value[:,4] = cnp.copy()
                value[:,5] = tnp.copy()
                value[:,6] = xnp.copy()
                fph_na = LinearNDInterpolator(list( zip(px, hy) ), value)
                self.fph_na = fph_na
                hpt_na = LinearNDInterpolator(list( zip(px, tnp) ), hy)
                self.hpt_na = hpt_na

                val1d[0 , :] = hgl
                val1d[1 , :] = mls
                val1d[2 , :] = mgs
                val1d[3 , :] = rls
                val1d[4 , :] = rgs
                val1d[5 , :] = kls
                val1d[6 , :] = kgs
                val1d[7 , :] = ts 
                val1d[8 , :] = cls
                val1d[9 , :] = cgs
                val1d[10, :] = sgm

                fp_na = interp1d(pnp, val1d, kind='linear',fill_value = 'extrapolate')
                self.fp_na = fp_na

                print('Sodium database done.')
                
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
        # + psi_ytchen: enthalpy-based Na properties 07.02.2024
        elif inp['type'] == 'na':
            # Sodium two-phase properties should be added here
            ini = inp['ini']
            if ini ==0: # psi_ytchen: calculation based on h, vectorized in 20240312
            # psi_ytchen: calculation based on h
                pmin= 1.000e1
                pmax= 25.00e6
                hmin= 210.0e3
                hmax= 5800.e3
                h = inp['h']
                p = inp['p']
                hnp = np.array(h)
                pnp = np.array(p)
                hnp[hnp < hmin] = hmin+1.0
                hnp[hnp > hmax] = hmax-1.0
                pnp[hnp < pmin] = pmin+1.0
                pnp[hnp > pmax] = pmax-1.0
                
                value = self.fph_na( (pnp, hnp) )
                rhof =  value[ : ,0]
                visf =  value[ : ,1]
                miuf =  value[ : ,2]
                kf   =  value[ : ,3]
                cpf  =  value[ : ,4]
                tf   =  value[ : ,5]
                xe   =  value[ : ,6]
                
                
                val1d = self.fp_na( pnp )
                hgl = val1d[0 , :]
                mls = val1d[1 , :]
                mgs = val1d[2 , :]
                rls = val1d[3 , :]
                rgs = val1d[4 , :]
                kls = val1d[5 , :]
                kgs = val1d[6 , :]
                ts  = val1d[7 , :]
                cls = val1d[8 , :]
                cgs = val1d[9 , :]
                sgm = val1d[10, :]
                
                if np.max(xe) < 0. :
                    rhol = rhof.copy()
                    miul = miuf.copy()
                    kl   = kf.copy()
                    cpl  = cpf.copy()
                    rhog = rgs.copy()
                    miug = mgs.copy()
                    kg   = kgs.copy()
                    cpg  = cgs.copy()
                elif np.min(xe) > 1. :
                    rhol = rls.copy()
                    miul = mls.copy()
                    kl   = kls.copy()
                    cpl  = cls.copy()
                    rhog = rhof.copy()
                    miug = miuf.copy()
                    kg   = kf.copy()
                    cpg  = cpf.copy()
                else:
                    rhol = np.zeros( len(xe) )
                    rhog = np.zeros( len(xe) )
                    miul = np.zeros( len(xe) )
                    miug = np.zeros( len(xe) )
                    kl   = np.zeros( len(xe) )
                    kg   = np.zeros( len(xe) )
                    cpl  = np.zeros( len(xe) )
                    cpg  = np.zeros( len(xe) )
                    for i in range( len(xe) ):
                        if xe[i] < 0.: 
                            rhol[i] = rhof[i]
                            rhog[i] = rgs[i]
                            miul[i] = miuf[i] 
                            miug[i] = mgs[i]
                            kl[i]   = kf[i]
                            kg[i]   = kgs[i]
                            cpl[i]  = cpf[i]
                            cpg[i]  = cgs[i]
                        elif xe[i] > 1.: 
                            rhog[i] = rhof[i]
                            rhol[i] = rls[i]
                            miug[i] = miuf[i]
                            miul[i] = mls[i]
                            kg[i]   = kf[i]
                            kl[i]   = kls[i]
                            cpg[i]  = cpf[i]
                            cpl[i]  = cls[i]
                        else:
                            rhol[i] = rls[i]
                            rhog[i] = rgs[i]
                            miul[i] = mls[i]
                            miug[i] = mgs[i]
                            kl[i]   = kls[i]
                            kg[i]   = kgs[i]
                            cpl[i]  = cls[i]
                            cpg[i]  = cgs[i]
                return {'tf':tf,'rhol':rhol,'rhog':rhog,'rhof':rhof,'miul':miul,'miug':miug,'miuf':miuf,'visf':visf,\
                'cpl':cpl,'cpg':cpg,'cpf':cpf,'kl':kl,'kg':kg,'kf':kf,'xe': xe,'ts':ts,'hgl':hgl,'sgm':sgm}


            elif ini > 0: # psi_ytchen: initialization run
            # psi_ytchen: calculation based on t
                t = inp['t']
                p = inp['p']
                hh = self.hpt_na(p,t)
                return {'h':hh}
        # - psi_ytchen: enthalpy-based Na properties 07.02.2024
        # lbe: liquid lead and bismuth (55%wt Bi, 45%wt Pb)
        # + psi_ytchen: enthalpy-based LBE properties 07.02.2024
        elif inp['type'] == 'lbe':
            ini = inp['ini']
            if ini ==0: # psi_ytchen: not initialization run
            # psi_ytchen: calculation based on h
                h = inp['h']
                hnp = np.array(h)
                value = self.fh_lbe(hnp)
                tf   = value[0,:]
                rhof = value[1,:]
                visf = value[2,:]
                cpf  = value[3,:]
                kf   = value[4,:]
                return {'tf':tf, 'rhof':rhof, 'visf':visf, 'cpf':cpf, 'kf':kf}
            elif ini > 0: # psi_ytchen: initialization run
            # psi_ytchen: calculation based on t
                t = inp['t']
                # OECD, Nuclear Energy Agency, Handbook on Lead-bismuth Eutectic Alloy and Lead Properties, Materials Compatibility, Thermalhydraulics and Technologies, OECD, 2015. https://doi.org/10.1787/42dcd531-en.
                # psi_ytchen: enthalpy (J/kg): @400-1100K equation from "Handbook on Lead-bismuth Eutectic Alloy and Lead Properties", p.98, same as the following ones
                hl = self.ht_lbe(t)
                return {'h':hl}
        # - psi_ytchen: enthalpy-based LBE properties 07.02.2024
        # + psi_ytchen: enthalpy-based h2o properties 06.02.2024
        elif inp['type'] == 'h2o':
            fluid = 'IF97::Water'
            ini   = inp['ini']
            if ini == 0: # psi_ytchen: not initialization run
            # psi_ytchen: calculation based on h
                h = inp['h']
                p = inp['p']
                hnp = np.array(h)
                pnp = np.array(p)
                value = self.fph_h2o( (pnp, hnp) )
                rhof =  value[ : ,0]
                visf =  value[ : ,1]
                miuf =  value[ : ,2]
                kf   =  value[ : ,3]
                cpf  =  value[ : ,4]
                tf   =  value[ : ,5]
                xe   =  value[ : ,6]
                
                
                val1d = self.fp_h2o( pnp )
                hgl = val1d[0 , :]
                mls = val1d[1 , :]
                mgs = val1d[2 , :]
                rls = val1d[3 , :]
                rgs = val1d[4 , :]
                kls = val1d[5 , :]
                kgs = val1d[6 , :]
                ts  = val1d[7 , :]
                cls = val1d[8 , :]
                cgs = val1d[9 , :]
                sgm = val1d[10, :]
                
                if np.max(xe) < 0. :
                    rhol = rhof.copy()
                    miul = miuf.copy()
                    kl   = kf.copy()
                    cpl  = cpf.copy()
                    rhog = rgs.copy()
                    miug = mgs.copy()
                    kg   = kgs.copy()
                    cpg  = cgs.copy()
                elif np.min(xe) > 1. :
                    rhol = rls.copy()
                    miul = mls.copy()
                    kl   = kls.copy()
                    cpl  = cls.copy()
                    rhog = rhof.copy()
                    miug = miuf.copy()
                    kg   = kf.copy()
                    cpg  = cpf.copy()
                else:
                    rhol = np.zeros( len(xe) )
                    rhog = np.zeros( len(xe) )
                    miul = np.zeros( len(xe) )
                    miug = np.zeros( len(xe) )
                    kl   = np.zeros( len(xe) )
                    kg   = np.zeros( len(xe) )
                    cpl  = np.zeros( len(xe) )
                    cpg  = np.zeros( len(xe) )
                    for i in range( len(xe) ):
                        if xe[i] < 0.: 
                            rhol[i] = rhof[i]
                            rhog[i] = rgs[i]
                            miul[i] = miuf[i] 
                            miug[i] = mgs[i]
                            kl[i]   = kf[i]
                            kg[i]   = kgs[i]
                            cpl[i]  = cpf[i]
                            cpg[i]  = cgs[i]
                        elif xe[i] > 1.: 
                            rhog[i] = rhof[i]
                            rhol[i] = rls[i]
                            miug[i] = miuf[i]
                            miul[i] = mls[i]
                            kg[i]   = kf[i]
                            kl[i]   = kls[i]
                            cpg[i]  = cpf[i]
                            cpl[i]  = cls[i]
                        else:
                            rhol[i] = rls[i]
                            rhog[i] = rgs[i]
                            miul[i] = mls[i]
                            miug[i] = mgs[i]
                            kl[i]   = kls[i]
                            kg[i]   = kgs[i]
                            cpl[i]  = cls[i]
                            cpg[i]  = cgs[i]
                return {'tf':tf,'rhol':rhol,'rhog':rhog,'rhof':rhof,'miul':miul,'miug':miug,'miuf':miuf,'visf':visf,\
                'cpl':cpl,'cpg':cpg,'cpf':cpf,'kl':kl,'kg':kg,'kf':kf,'xe': xe,'ts':ts,'hgl':hgl,'sgm':sgm}
            elif ini > 0: # psi_ytchen: initialization run
            # psi_ytchen: calculation based on t
                t = inp['t']
                p = inp['p']
                hh = self.hpt_h2o(p,t)
                return {'h':hh}
        # - psi_ytchen: enthalpy-based h2o properties 06.02.2024
        # ss316: stainless steel type of 316
        elif inp['type'] == 'ss316':
            t = inp['t']
            # density (kg/m3): @300K equation from Leibowitz, et al, "Properties for LMFBR safety analysis", ANL-CEN-RSD-76-1 (1976), p.117
            rho = 7954.
            # specific heat (J/kg-K): Leibowitz, et al, "Properties for LMFBR safety analysis", ANL-CEN-RSD-76-1 (1976), p.100. Note that 1 mol of SS316 = 10.165 kg (https://www.webqc.org/molecular-weight-of-SS316.html) and 1 cal = 4.184 J
            cp = (6.181 + 1.788e-3*t)*4.184/55.845e-03
            # thermal conductivity (W/m-K): Leibowitz, et al, "Properties for LMFBR safety analysis", ANL-CEN-RSD-76-1 (1976), p.100.
            k = 9.248 + 1.571e-2*t
            return {'rho':rho, 'cp':cp, 'k':k}
            
            
        elif inp['type'] == 'AISI316L':
            t = inp['t']
            # density (kg/m3): @300K equation from Leibowitz, et al, "Properties for LMFBR safety analysis", ANL-CEN-RSD-76-1 (1976), p.117
            rho = 7900.
            # specific heat (J/kg-K): Leibowitz, et al, "Properties for LMFBR safety analysis", ANL-CEN-RSD-76-1 (1976), p.100. Note that 1 mol of SS316 = 10.165 kg (https://www.webqc.org/molecular-weight-of-SS316.html) and 1 cal = 4.184 J
            cp = 599.
            # thermal conductivity (W/m-K): Leibowitz, et al, "Properties for LMFBR safety analysis", ANL-CEN-RSD-76-1 (1976), p.100.
            # k = 9.248 + 1.571e-2*t
            k = 25.0
            return {'rho':rho, 'cp':cp, 'k':k}
            
        elif inp['type'] == 'NA_SATL':
            t = inp['t']
            t3=2503.70
            rho = 219.0 + 275.32*(1.0 - t/t3) + 511.58*(1.0 - t/t3)**0.5
            cp = 1658.2 - 0.84750*t + 4.4541e-04*t**2 - 2.9926e6/t**2
            k  = 124.67 - 0.11381*t + 5.5226e-5*t**2 - 1.1842e-8*t**3

            return {'rho':rho, 'cp':cp, 'k':k}
            
        # 9Cr1Mo: tube material of SGTF tube and shell material
        elif inp['type'] == '9Cr1Mo':
            t = inp['t']
            # density (kg/m3): from Vkiram SGTF article
            rho = 7600.
            # specific heat (J/kg-K): from Vkiram SGTF article
            cp = 600.
            # thermal conductivity (W/m-K): from Vkiram SGTF article
            A = 28.568
            B =-0.17859e-2
            C =-0.39325e-4
            D = 0.22380e-6
            E =-0.31979e-9
            F = 0.12585e-12
            tc = t - 273.15
            k = A + B*tc + C*tc**(2.0) + D*tc**(3.0) + E*tc**(4.0) + F*tc**(5.0)
            return {'rho':rho, 'cp':cp, 'k':k}
            
            # 2.25Cr1Mo: tube material of ETEC tube and shell material
        elif inp['type'] == '2.25Cr1Mo':
            t = inp['t']
            # density (kg/m3): from Vkiram SGTF article
            rho = 7750.
            # specific heat (J/kg-K): from Vkiram SGTF article
            cp = 600.
            # thermal conductivity (W/m-K): from Vkiram SGTF article
            A = 36.177
            B = 0.54737e-2
            C = 0.35874e-4
            D =-0.24726e-6
            E = 0.32447e-9
            F =-0.11300e-12
            tc = t - 273.15
            k = A + B*tc + C*tc**(2.0) + D*tc**(3.0) + E*tc**(4.0) + F*tc**(5.0)
            return {'rho':rho, 'cp':cp, 'k':k}
            
        elif inp['type'] == 'NiCr':
            t = inp['t']
            # density (kg/m3):constant
            rho = 8300.
            # specific heat (J/kg-K): constant
            cp = 670.
            # thermal conductivity (W/m-K): constant
            k = 26.0
            return {'rho':rho, 'cp':cp, 'k':k}
            
        elif inp['type'] == 'MgO':
            t = inp['t']
            # density (kg/m3):constant
            rho = 2890.
            # specific heat (J/kg-K): constant
            cp = 1185.
            # thermal conductivity (W/m-K): constant
            k = 3.
            return {'rho':rho, 'cp':cp, 'k':k}
 
        # AISI316_powder: This material properties is only NACIE-UP benchmark
        elif inp['type'] == 'AISI316_powder':
            t = inp['t']
            rho = 7954.
            cp = (6.181 + 1.788e-3*t)*10.165*4.184
            # This material prporety is from NACIE-UP benchmark spec
            k = max(0.3 + 0.005*(t - 273.15 - 200), 0.30) 
            return {'rho':rho, 'cp':cp, 'k':k}

        # bn: boron nitide
        elif inp['type'] == 'bn':
            t = inp['t']
            tc = min(t -273.15, 1800.0)
            # density (kg/m3): I.Di Piazza, et al., Benchmark specifications for NACIE-UP facility: non-uniform power distribution tests, ENEA Report, NA-I-R-542, Feb. 2023
            rho = 2000. # Use ENEA data
            # specific heat (J/kg-K): I.Di Piazza, et al., Benchmark specifications for NACIE-UP facility: non-uniform power distribution tests, ENEA Report, NA-I-R-542, Feb. 2023
            cp = 800. # Use ENEA data
            # thermal conductivity (W/m-K): I.Di Piazza, et al., Benchmark specifications for NACIE-UP facility: non-uniform power distribution tests, ENEA Report, NA-I-R-542, Feb. 2023
            k = 25.578 - 2.416*math.log(tc) # Use ENEA data
            return {'rho':rho, 'cp':cp, 'k':k}

        # cu: copper
        elif inp['type'] == 'cu':
            t = inp['t']
            # density (kg/m3): I.Di Piazza, et al., Benchmark specifications for NACIE-UP facility: non-uniform power distribution tests, ENEA Report, NA-I-R-542, Feb. 2023
            rho = 8933.
            # specific heat (J/kg-K): I.Di Piazza, et al., Benchmark specifications for NACIE-UP facility: non-uniform power distribution tests, ENEA Report, NA-I-R-542, Feb. 2023
            cp = 385.
            # thermal conductivity (W/m-K): I.Di Piazza, et al., Benchmark specifications for NACIE-UP facility: non-uniform power distribution tests, ENEA Report, NA-I-R-542, Feb. 2023
            k = 401
            return {'rho':rho, 'cp':cp, 'k':k}

    #----------------------------------------------------------------------------------------------
    # Nusselt number: self is a 'data' object created in B, inp is a dictionary of input data dependent on the case
    def qfluxcal(self, inp):
        HTCMOD = 0 # psi_ytchen: Heat Transfer Mode
        material_type = inp['type']
        Re = inp['re']
        Re = max(Re, 1.0) # psi_ytchen: to avoid math domain error
        Pr = inp['pr']
        Tf = inp['tmp'][0]
        Tw = inp['tmp'][1]
        dh = inp['dhyd']
        kf = inp['prop']['kf']
        PP = inp['p']
        if material_type == 'h2o':
            xe = inp['prop']['xe']
            Ts = inp['prop']['ts']
            Pmpa = PP/1.0e6
            if xe < 0.0:
              nuLam = 4.36
              nuTurb = 0.0233*Re**(0.80)*Pr**(0.40)
              nuSp = max(nuLam, nuTurb)
              hSp = nuSp*kf/dh # psi_ytchen: single-phase heat transfer coefficient
              qSp = hSp*(Tw - Tf)
            # psi_ytchen: the fluid is subcooled liquid
              if Tw < Ts: 
              # psi_ytchen: single-phase liquid heat transfer
                qflux = qSp
              else: # psi_ytchen: Tw >= Ts condition
                qThom= 2000.0*math.exp(Pmpa/4.34)*(Tw - Ts)**(2.0)
                if qThom < qSp: # Tw < Tonb condition
                  qflux = qSp
                else:           # Tw >= Tonb condition
                  hgl = inp['prop']['hgl']
                  sgm = inp['prop']['sgm']
                  rhol= inp['prop']['rhol']
                  rhog= inp['prop']['rhog']
                  qchf = 0.14*hgl*(sgm*9.8*rhog**(2.0)*(rhol-rhog))**(0.25)
                  Tchf = Ts + (qchf*math.exp(-Pmpa/4.34)/2000.0)**(0.50) # Obtained by making qchf .eq. qThom
                  if Tw < Tchf: # Tonb <= Tw < Tchf condition
                    qflux = qThom 
                    HTCMOD = 1 # subcooled nucleate boiling
                  else: # Tw > Tchf condition
                    kg  = inp['prop']['kg']
                    cpg = inp['prop']['cpg']
                    miug= inp['prop']['miug']
                    if Pmpa < 9.0:
                      Tmfb = 557.90 + 44.1*Pmpa - 3.72*Pmpa**(2.0)
                    else:
                      Tmfb = 647.09 + 0.71*Pmpa
                    # Tmfb is Leidenfrost temperature of water
                    if Tw < Tmfb: # Tchf < Tw < Tmfb condition, subcooled transition boiling
                      qmfb = 0.14*hgl*(kg**(2.0)*cpg*rhog*9.8*(rhol - rhog)/miug)**(1.0/3.0)*(Tmfb - Ts)
                      ntmp = math.log( (Tmfb-Ts)/(Tw-Ts) )/math.log( (Tmfb-Ts)/(Tchf-Ts) )
                      qflux = qmfb*(qchf/qmfb)**(ntmp)
                      HTCMOD = 2 # subcooled transition boiling
                    else: # Tmfb < Tw condition, subcooled stable film boiling
                      qflux = 0.14*hgl*(kg**(2.0)*cpg*rhog*9.8*(rhol - rhog)/miug)**(1.0/3.0)*(Tw - Ts)
                      HTCMOD = 3 # subcooled film boiling
            elif xe < 1.0: # 0.0 < xe < 1.0
              hgl = inp['prop']['hgl']
              sgm = inp['prop']['sgm']
              rhol= inp['prop']['rhol']
              rhog= inp['prop']['rhog']
              kg  = inp['prop']['kg']
              cpg = inp['prop']['cpg']
              miug= inp['prop']['miug']
              miul= inp['prop']['miul']
              G = inp['Gtot']
              
              alpha = rhol*xe/(rhol*xe + rhog*(1.0-xe)) # void fraction
              reg = G*dh/miug
              Prg = miug*cpg/kg
              
              #xcr = min( (0.44 - 0.006*Pmpa)*G**(0.114), 0.99 ) #for 0.012 < dh and dh < 0.013
              
              #omega = G*miul/(sgm*rhol)*(rhol/rhog)**(1.0/3.0)
              #xcr = min( 0.30 + 0.70*math.exp(-45.0*omega)*(0.008/dh)**(0.15), 0.99 )
              XX = Pmpa/9.8
              xcr = ( 0.39+1.57*XX-2.04*XX**(2.0)+0.68*XX**(3.0) )*(G/1000.)**(-0.5)*(dh/0.008)**(0.15)
              xcr = min(xcr, 0.99)
              #xcr = 0.99
              
              if xe < xcr: # 0.0 <= xe < xcr, not dryout condition
                qchf = 0.14*hgl*(sgm*9.8*rhog**(2.0)*(rhol-rhog))**(0.25)*(1.0 - alpha) # qchf considerting dryout effect
                Tchf = Ts + (qchf*math.exp(-Pmpa/4.34)/2000.0)**(0.50) # calculate Tchf
                # calculate qDrycr
                alpcr = rhol*xcr/(rhol*xcr + rhog*(1.0-xcr)) # void fraction
                Ycr = 1.0 - 0.1*((rhol/rhog - 1.0)*(1.0 - xcr))**(0.40)
                
                nuDrycr = 0.0230*(reg*(xcr/alpcr))**(0.80)*Prg**(0.40)*Ycr*(Tf/Tw)**(0.50) # Mitropolsky
                # nuDrycr = 0.052*(reg*(xcr/alpcr))**(0.688)*Prg**(1.26)*Ycr**(-1.06)*(Tf/Tw)**(0.50) # Groeneveld Annular
                # nuDrycr = 1.09e-03*(reg*(xcr/alpcr))**(0.989)*Prg**(1.41)*Ycr**(-1.15)*(Tf/Tw)**(0.50) # Groeneveld Circular
                # nuDrycr = 3.27e-03*(reg*(xcr/alpcr))**(0.901)*Prg**(1.32)*Ycr**(-1.50)*(Tf/Tw)**(0.50) # Groeneveld General
                hDrycr = nuDrycr*kg/dh
                qDrycr = hDrycr*(Tw - Tf) # Totally dryout heat flux
                phy = math.exp( xe/(xe - xcr) ) # psi_ytchen: interpolate factor for trial and error 20240324
                if Tw < Tchf:
                  qflux = 2000.0*math.exp(Pmpa/4.34)*(Tw - Ts)**(2.0) # Thom correlation
                  HTCMOD = 4 # saturated nucleate boiling
                else: # Tchf < Tw
                  if Pmpa < 9.0:
                    Tmfb = 557.90 + 44.1*Pmpa - 3.72*Pmpa**(2.0)
                  else:
                    Tmfb = 647.09 + 0.71*Pmpa
                  if Tw < Tmfb:
                    qmfb = 0.14*hgl*(kg**(2.0)*cpg*rhog*9.8*(rhol - rhog)/miug)**(1.0/3.0)*(Tmfb - Ts)
                    ntmp = math.log( (Tmfb-Ts)/(Tw-Ts) )/math.log( (Tmfb-Ts)/(Tchf-Ts) )
                    qflux = qmfb*(qchf/qmfb)**(ntmp)
                    HTCMOD = 5 # saturated transion boiling
                  else: # Tmfb < Tw condition, subcooled stable film boiling
                    qflux = 0.14*hgl*(kg**(2.0)*cpg*rhog*9.8*(rhol - rhog)/miug)**(1.0/3.0)*(Tw - Ts)
                    HTCMOD = 6 # saturated film boiling
                qflux = qflux*phy + (1.0-phy)*qDrycr # qflux is interpolated for trial and error 20240324
              else: # xcr <= xe < 1.0, dryout condition
                Y = 1.0 - 0.1*((rhol/rhog - 1.0)*(1.0 - xe))**(0.40)
                nuDry = 0.0230*(reg*(xe/alpha))**(0.80)*Prg**(0.40)*Y*(Tf/Tw)**(0.50) # Mitropolsky ！！！0.4 or 0.80
                # nuDry = 0.052*(reg*(xe/alpha))**(0.688)*Prg**(1.26)*Y**(-1.06)*(Tf/Tw)**(0.50) # Groeneveld Annular
                # nuDry = 1.09e-03*(reg*(xe/alpha))**(0.989)*Prg**(1.41)*Y**(-1.15)*(Tf/Tw)**(0.50) # Groeneveld Circular
                # nuDry = 3.27e-03*(reg*(xe/alpha))**(0.901)*Prg**(1.32)*Y**(-1.50)*(Tf/Tw)**(0.50) # Groeneveld General
                hDry = nuDry*kg/dh
                qDry = hDry*(Tw - Tf) # Totally dryout heat flux
                qflux = qDry
                HTCMOD = 7 # saturated dryout boiling
            else: # xe >= 1.0 condition
              nuLam = 4.36
              nuTurb = 0.0230*Re**(0.80)*Pr**(0.40)*(Tf/Tw)**(0.50)
              #nuTurb = 0.0270*Re**(0.80)*Pr**(1.0/3.0)*(Tf/Tw)**(0.50)
              nuSp = max(nuLam, nuTurb)
              hSp = nuSp*kf/dh # psi_ytchen: single-phase heat transfer coefficient
              qflux = hSp*(Tw - Tf)
              HTCMOD = 8 # single-phase steam

        elif material_type == 'na'  :
            pe = Re*Pr
            Pmpa = PP/1.0e6
            xe = inp['prop']['xe']
            Ts = inp['prop']['ts']
            if xe < 0.0:
                
                if 'p2d' in inp:
                    # pin bundle geometry
                    p2d = inp['p2d']
                    # forced convection in a pin bundle (Mikityuk, NED 2008)
                    nu = 0.047*(1.0-math.exp(-3.8*(p2d-1.0))) * ((pe)**0.77 + 250.0)
                else:
                    #round tube geometry
                    nu = 4.8 + 0.025 * (pe)**0.8
                hex = nu*kf/dh
                qflux_SP = hex*(Tw - Tf)
                # psi_ytchen: the fluid is subcooled liquid
                if Tw < Ts:
                    # psi_ytchen: single-phase heat transfer
                    qflux = qflux_SP
                else: # Tw >= Ts
                    # Experimental study on boiling two-phase of liquid sodium along a 7-rod bundle – Part II: Heat transfer characteristics
                    A = 4.80
                    m = 0.70
                    n = 0.15 # Yandong Hou
                    qBoil = ( A*Pmpa**(n)*(Tw-Ts) )**( 1.0/(1.-m) )
                    
                    if qBoil < qflux_SP:
                        qflux = qflux_SP
                    else:
                        qflux = qBoil
                        HTCMOD = 1 # sodium subcooled boiling regime
            elif xe <= 1.0:
                miul= inp['prop']['miul']
                kl  = inp['prop']['kl']
                cpl = inp['prop']['cpl']
                rhol= inp['prop']['rhol']
                rhog= inp['prop']['rhog']
                kg  = inp['prop']['kg']
                cpg = inp['prop']['cpg']
                miug= inp['prop']['miug']
                G = inp['Gtot']
                
                rel = G*dh/miul
                prl = cpl*miul/kl
                pel = rel*prl
                if 'p2d' in inp:
                    # pin bundle geometry
                    p2d = inp['p2d']
                    # forced convection in a pin bundle (Mikityuk, NED 2008)
                    nul = 0.047*(1.0-math.exp(-3.8*(p2d-1.0))) * ((pel)**0.77 + 250.0)
                else:
                    #round tube geometry
                    nul = 4.8 + 0.025 * (pel)**0.8
                hLiq = nul*kl/dh
                qLiq = hLiq*(Tw - Ts)
                A = 4.80
                m = 0.70
                n = 0.15 # Yandong Hou
                if Tw > Ts:
                    qBoil = max( qLiq, ( A*Pmpa**(n)*(Tw-Ts) )**( 1.0/(1.-m) ) )
                else:
                    qBoil = qLiq
 
                reg = G*dh/miug
                Prg = miug*cpg/kg
                nug = 0.0230*reg**(0.80)*Prg**(0.40)###*(Ts/Tw)**(0.50)
                hexg= nug*kg/dh
                qGas= hexg*(Tw - Ts)
                
                xcr = 0.300
                acr1= 0.957 # critical void fraction const
                acr2= rhol*xcr/(rhol*xcr + rhog*(1.0-xcr)) # critical void fraction by definition
                
                if xe <= xcr:
                    
                    qflux = qBoil
                    Cfwlcr = ( (1.0-max(acr1,acr2))/(1.0-acr1) )**(0.5)
                    qdrycr = Cfwlcr*qBoil + (1.0-Cfwlcr)*qGas
                    phy = math.exp( xe/(xe - xcr) ) 
                    qflux = qBoil*phy + qdrycr*(1.0-phy)
                    HTCMOD = 2 # sodium saturated boiling regime
                    
                else: # xe > xcr
                    # Aurelia CHENU PhD Thesis
                    alpha = rhol*xe/(rhol*xe + rhog*(1.0-xe))  # void fraction
                    Cfwl = ( (1.0-alpha)/(1.0-min(acr1,acr2)) )**(0.5)*( (1.0-xe)/(1.0-xcr) )**(1.5)
                    qflux  = qBoil*Cfwl + qGas*(1. - Cfwl)
                    HTCMOD = 3 # sodium post dryout boiling regime
                    
            else: # single-phase sodium vapor
                nu = 0.0230*Re**(0.80)*Pr**(0.40)###*(Ts/Tw)**(0.50)
                hex = nu*kf/dh
                qflux = hex*(Tw - Tf)
                HTCMOD = 4 # single-phase sodium vapor regime
                
        elif material_type == 'lbe' :
            pe = Re*Pr
            if 'p2d' in inp:
                # pin bundle
                p2d = inp['p2d']
                # forced convection in a pin bundle (Mikityuk, NED 2008)
                nu = 0.047*(1.0-math.exp(-3.8*(p2d-1.0))) * ((pe)**0.77 + 250.0)
            else:
                # round tube
                nu = 4.8 + 0.025 * (pe)**0.8
            hex = nu*kf/dh
            qflux = hex*(Tw - Tf)
        return [qflux,HTCMOD]

    #----------------------------------------------------------------------------------------------
    # Friction factor: self is a 'data' object created in B, inp is a dictionary of input data dependent on the case
    def fricfac(self, inp):



        if 'p2d' in inp: # psi_ytchen: wire-wrapped rod bundles or barerod
            
            #Chen, S. K., et al. (2018). "The upgraded Cheng and Todreas correlation for pressure drop 
            #in hexagonal wire-wrapped rod bundles." Nuclear Engineering and Design 335: 356-373.
            # psi_ytchen: for calculation of a single constant, 'math' is faster than 'numpy'
            p2d = inp['p2d']
            sb2st  = inp['sb2st'] # wetted perimeter ratio: (rod + wire)/(rod + wire + wall), influence of nrod
            re = np.array(inp['re'])
            re[re < 1.] = 1.
            
            if 'h2d' in inp:
            # psi_ytchen: wire-wrapped rod bundle, use rehme 1973 model
                h2d = inp['h2d']
                d2h = 1.0/h2d
                p2h = p2d*d2h
                FF  = p2d**(0.5) + ( 7.6*(p2h-d2h)*p2h**(2.0) )**(2.16)
                fric = ( FF**(0.5)*64.0/re + FF**(0.9335)*0.0816/re*(0.133) )*sb2st
                return fric
            else: 
            # psi_ytchen: barerod bundle by Cheng

                if 1.0 < p2d and p2d <= 1.1:
                    LL = [26.0000, 888.2, -3334.0]
                    TT = [0.09378, 1.398, -8.664]
                elif p2d <=1.5:
                    LL = [62.97,  216.9,   -190.2]
                    TT = [0.1458, 0.03632, -0.03333]
                else:
                    print('warning, barerod correlation input out of range, check input')
                    
                CfbL = LL[0] + LL[1]*(p2d-1.) + LL[2]*(p2d-1.)**(2.0)
                CfbT = TT[0] + TT[1]*(p2d-1.) + TT[2]*(p2d-1.)**(2.0)
                RebL = 320*(10**(p2d-1.0))
                RebT = 10000*(10**(0.7*(p2d-1.0)))
                phi = np.log10(re/RebL)/math.log10(RebT/RebL)
                phi = np.array(phi)
                if np.max(phi) < 0.0:
                    # psi_ytchen: pure laminar friction factor
                    return CfbL/re*sb2st
                elif np.min(phi) > 1.0:
                    # psi_ytchen: pure turbulent friction factor
                    return CfbT/re**0.18*sb2st
                else:
                    # psi_ytchen: transition friction factor
                    phi[phi < 0.] = 0.0
                    phi[phi > 1.] = 1.0
                    return (CfbL/RebL*(1 - phi)**(1/3)*(1 - phi**7) + (CfbT/RebT**0.18)*phi**(1/3))*sb2st

        
        else: # psi_ytchen: Churchill correlation
            re = np.array(inp['re'])
            re[re < 1.] = 1.
            a = (2.457*np.log(1.0/((7.0/re)**(0.9)) ))**16.0
            b = (37530./re)**16.0
            fric = 8.0*( (8.0/re)**(12) + 1.0/(a + b)**1.50 )**(1.0/12.0)
            # psi_ytchen:
            return fric