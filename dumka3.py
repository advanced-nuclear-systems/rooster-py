#--------------------------------------------------------------------------------------------------
# DUMKA3 - integrates initial value problems for systems of first order ordinary differential equations y'=f(y,t).
# It is based on a family of explicit Runge-Kutta-Chebyshev formulas of order three.
# It uses optimal third order accuracy stability polynomials with the largest stability region along the negative real axis.
#
# The algorithm dumka3 and dumka4 is described in the article:
# Medovikov A.A. High order explicit methods for parabolic equations. BIT, V38,No2,pp.372-390,1998
#
# Many thanks Prof. Lebedev V.I. for exellent algorithm of ordering
# roots of stability polynomials, and Prof. Wanner G. and Hairer E.
# for step size control algorithm.
#
#--------------------------------------------------------------------------------------------------
import math
def dumka3(inp):

    #constants
    N_DEG = [3,6,9,15,21,27,36,48,63,81,135,189,243,324,432,576]
    INDEX_FIRST = [1,2,4,7,12,19,28,40,56,77,104,149,212,293,401,545]
    INDEX_LAST = [1,3,6,11,18,27,39,55,76,103,148,211,292,400,544,736]

    # initialization
    y = inp['y']
    t = inp['tstart']
    tend = inp['tend']
    h = inp['h0']
    rhs = inp['rhs']
    C_1 = inp['C_1']
    C_2 = inp['C_2']
    C_3 = inp['C_3']
    C_4 = inp['C_4']
    rtol = inp['rtol']
    atol = inp['atol']

    z = rhs(t, y)
    N = len(z)

    # time step loop
    index = 0
    while t < tend :
        save = [t, y, z]
        h = min(h, tend - t)
        n_pol = 0
        for k in range(INDEX_FIRST[index]-1,INDEX_LAST[index]) :
            dt = [C_1[k]*h, (C_2[k] + C_3[k])*h/2.0, (1.0 - C_1[k])*h/2.0, C_3[k]*h, (C_2[k] - C_1[k])*h]

            n_pol += 3

            for i in range(N) :
                y[i] += z[i]*dt[0]
            z1 = rhs(t + dt[0], y)
            
            for i in range(N) :
                y[i] += z1[i]*dt[3]
                if n_pol == N_DEG[index] :
                    y[i] += z[i]*dt[4]
                    z1[i] = z1[i]*dt[1] - z[i]*dt[2]
            z = rhs(t + 2.0*dt[1], y)

            for i in range(N) :
                y[i] += z[i]*C_4[k]*h
                if n_pol == N_DEG[index] :
                    z1[i] += z[i]*dt[2]

            t += h
            z = rhs(t, y)

        # evaluate error
        err = [0] * N
        for i in range(N) :
            tol = rtol[i]*abs(y[i]) + atol[i]
            err[i] = (z1[i] - dt[1]*z[i]) / tol
        eps = 0.0
        for i in range(N) :
            eps += err[i]*err[i]
        eps = math.sqrt(eps / float(N))
        if eps == 0.0 : eps = 1e-14

        # predict next time step
        frac = (1.0/eps)**(1.0/3.0)
        if eps > 1.0 :
            print('failed time step')
            if index < 15 : index += 1
            h = h * 0.8*min(1.0, max(0.1, 0.8*frac))
            [t, y, z] = save
        else :
            if index > 0 : index -= 1
            h = h * min(2.0 ,max(0.1, 0.8*frac))
        print(('index={:2d} n_pol={:3d} t={:.5f} h={:.5f} y[0]={:.5f} y[1]={:.5f}' ).format(index, n_pol, t, h, y[0], y[1]))
    return