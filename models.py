# Class containing code for different GFET models

# For Physical Constants
from scipy import constants as consts

import scipy.integrate as integrate
import scipy.special as special

from numpy.lib.scimath import sqrt as csqrt
from numpy.linalg import norm
import numpy as np
import cmath

# Model of Rodriguez et al.
# Ref.: S. Rodriguez et al., "A Comprehensive Graphene FET Model for Circuit Design,"
# in IEEE Transactions on Electron Devices, vol. 61, no. 4, pp. 1199-1206, April 2014. (doi.org:/10.1109/TED.2014.2302372)

class RodriguezGFET:

    def __init__(self, params, ivSweep, transSweep, eps):
        self.params = params
        self.ivVds = ivSweep["Vds"]
        self.ivVtg = ivSweep["Vtg"]
        self.transVds = transSweep["Vds"]
        self.transVtg = transSweep["Vtg"]
        self.eps = eps
    
    def calculateTransferChars(self):
        
        tox = float(self.params[0])*10**(-9)
        W = float(self.params[1])*10**(-6)
        L = float(self.params[2])*10**(-6)
        mu = float(self.params[3])
        Ep = float(self.params[4])
        N = float(self.params[5])
        
        er = float(self.eps.get().split("(")[1].replace(")",""))
        Ct = er*consts.epsilon_0/tox
        w = (2.24*10**(13))/consts.pi
        Vg0 = consts.elementary_charge*N/Ct
        Cgs = Ct*W*L
        Cgd = (Ct*W*L)/2
        
        Ids = []

        for i in range(len(self.transVds)):
            Id = []
            Vds = self.transVds[i]
            for j in range(len(self.transVtg)):
                Veff = self.transVtg[j] + Vg0
                Id.append(abs((mu*W*Ct*(Veff-0.5*Vds))/(L/Vds +
                            (mu/w)*(csqrt(consts.pi*Ct/consts.elementary_charge))*(csqrt(Veff-0.5*Vds)))))
            Ids.append(Id)
            
        # Calculate transconductance & transit frequency
        gm =[]
        fT = []
        for entry in Ids:
            gm.append([abs(i/j) for i, j in zip(entry, self.transVtg)])
            res = []
            for entry in gm:
                res = [i/(2*consts.pi*(Cgs+Cgd)) for i in entry]
            fT.append(res)
        return Ids, gm, fT

    def calculateIVChars(self):
        
        tox = float(self.params[0])*10**(-9)
        W = float(self.params[1])*10**(-6)
        L = float(self.params[2])*10**(-6)
        mu = float(self.params[3])
        Ep = float(self.params[4])
        N = float(self.params[5])
        
        er = float(self.eps.get().split("(")[1].replace(")",""))
        Ct = er*consts.epsilon_0/tox
        w = (2.24*10**(13))/consts.pi
        Vg0 = consts.elementary_charge*N/Ct
                
        Ids = []

        for i in range(len(self.ivVtg)):
            Id = []
            Vtg = self.ivVtg[i]
            for j in range(len(self.ivVds)):
                Veff = Vtg + Vg0
                Vds = self.ivVds[j]

                if Vds ==0:
                    Id.append(0)
                else:
                    Id.append(abs((mu*W*Ct*(Veff-0.5*Vds))/(L/Vds +
                            (mu/w)*(csqrt(consts.pi*Ct/consts.elementary_charge))*(csqrt(Veff-0.5*Vds)))))
            Ids.append(Id)
        return Ids   


# Model of Thiele et al. (https://doi.org/10.1063/1.3357398)
class ThieleGFET:

    def __init__(self, params, ivSweep, transSweep, eps):
        self.params = params
        self.ivVds = ivSweep["Vds"]
        self.ivVtg = ivSweep["Vtg"]
        self.ivVbg = ivSweep["Vbg"]
        self.transVds = transSweep["Vds"]
        self.transVtg = transSweep["Vtg"]
        self.transVbg = transSweep["Vbg"]
        self.eps = eps
        self.N = 200 # Num Iterations for self-consisten eqn
        self.gds = []
    
    def calculateTransferChars(self):

        # Cq in terms of Bch
        def f(Vch):
            return Vch*(2*consts.e**3)/(consts.pi*(consts.hbar*vF)**2)

        def fixedp(x0, tol, N):
            e = 1
            xp = []
            itr=0

            while(e>tol and itr<N):
                x = f(x0)
                e = norm(x0-x)
                x0 = x
                xp.append(x0)
                itr = itr +1

            return x, xp

        # Integral for denominator, can't integrate directly a power of a lambda function!
        def integrand(Vx, rho, omega):
            A = 10**(-3)
            return ((consts.pi*rhosh(Vx))**(0.5+A*Vx**2))/omega
        
        tox = float(self.params[0])*10**(-9)
        W = float(self.params[1])*10**(-6)
        L = float(self.params[2])*10**(-6)
        mu = float(self.params[3])
        Ep = float(self.params[4])
        N = float(self.params[5])
        
        er = float(self.eps.get().split("(")[1].replace(")",""))
        Ct = er*consts.epsilon_0/tox
        w = (2.24*10**(13))/consts.pi
        Vg0 = consts.elementary_charge*N/Ct
        vF = 10**7 # m/s
        Ctg = Ct*W*L
        Cbg = Ctg
        Cgd = (Ct*W*L)/2

        Rs = 400
        Rd = Rs

        omega = w/consts.hbar#(Ep*consts.e)/consts.hbar
                
        Ids = []
        gds = []

        for i in range(len(self.transVds)):
            Id = []
            Vds = self.transVds[i]
            Vch = self.transVtg[0] - Vg0 # Vch is initially the first voltage
            for j in range(len(self.transVtg)):
                Vtg = self.transVtg[j]
                
                e = 0.001 # error/Volts
                Cq, Varr = fixedp(Vch, e, self.N)

                rhosh = lambda Vx: abs(-0.5*Cq*((Vtg-Vx)*(Ctg/Ctg+0.5*Cq)))/consts.e
                num = consts.e*mu*W*integrate.quad(rhosh, 0 ,Vds)[0]
                den = L - mu*integrate.quad(integrand, 0, Vds, args=(rhosh, omega))[0]
                
                Id.append(abs(num/den))
            Ids.append(Id)

        # Calculate transconductance & transit frequency
        gm =[]
        fT = []
        for index, entry in enumerate(Ids):
            gm.append([abs(i/j) for i, j in zip(entry, self.transVtg)])
            gd = self.gds[index]
            res = []
            
            for entry in gm:
                res = [i/(2*consts.pi*((Ctg+Cgd)*(1+gd*(Rs+Rd)+Cgd*i*(Rs+Rd)))) for i in entry]
            fT.append(res)
        return Ids, gm, fT


    #       IV Chars        #
    #                       #
    def calculateIVChars(self):
        # Cq in terms of Bch
        def f(Vch):
            return Vch*(2*consts.e**3)/(consts.pi*(consts.hbar*vF)**2)

        def fixedp(x0, tol, N):
            e = 1
            xp = []
            itr=0
            while(e>tol and itr<N):
                x = f(x0)
                e = norm(x0-x)
                x0 = x
                xp.append(x0)
                itr = itr +1
            return x, xp

        # Integral for denominator, can't integrate directly a power of a lambda function!
        def integrand(x, rho):
            return rhosh(x)**0.5
        
        tox = float(self.params[0])*10**(-9)
        W = float(self.params[1])*10**(-6)
        L = float(self.params[2])*10**(-6)
        mu = float(self.params[3])
        Ep = float(self.params[4])
        N = float(self.params[5])
        
        er = float(self.eps.get().split("(")[1].replace(")",""))
        Ctg = er*consts.epsilon_0/tox
        w = (2.24*10**(13))/consts.pi
        vF = 10**8 # m/s
        Vg0 = consts.elementary_charge*N/Ctg
        omega = w / consts.hbar
                
        Ids = []

        for i in range(len(self.ivVtg)):
            Id = []
            Vtg = self.ivVtg[i]
            Vch = Vtg - Vg0 # Vch is initially the first voltage
            
            for j in range(len(self.ivVds)):
                Vds = self.ivVds[i]

                if Vds ==0:
                    Id.append(0)
                else:
                    e = 0.001 # error/Volts
                    Cq, Varr = fixedp(Vch, e, self.N)
                    Vch = (Vtg-Vds)*(Ctg/Ctg+0.5*Cq)
                    rhosh = lambda x: abs(-0.5*Cq*((Vtg-x)*(Ctg/Ctg+0.5*Cq)))/consts.e

                    num = consts.e*mu*W*integrate.quad(rhosh, 0 ,Vds)[0]
                    den = L - mu*(cmath.sqrt(consts.pi)/omega)*integrate.quad(integrand, 0, Vds, args=(rhosh))[0]
                    Id.append(abs(num/den))
            Ids.append(Id)

        for index, entry in enumerate(Ids):
            dIds = max(entry)-min(entry)
            dVds = max(self.ivVds)-min(self.ivVds)
            self.gds.append(dIds/dVds)
        return Ids


# Model of Hu et al. (https://doi.org/10.1063/1.3357398)

class HuGFET:

    def __init__(self, params, ivSweep, transSweep, eps):
        self.params = params
        self.ivVds = ivSweep["Vds"]
        self.ivVtg = ivSweep["Vtg"]
        self.ivVbg = ivSweep["Vbg"]
        self.transVds = transSweep["Vds"]
        self.transVtg = transSweep["Vtg"]
        self.transVbg = transSweep["Vbg"]
        self.eps = eps
    
    def calculateTransferChars(self):
        
        tox = float(self.params[0])*10**(-9)
        W = float(self.params[1])*10**(-6)
        L = float(self.params[2])*10**(-6)
        mu = float(self.params[3])
        Ep = float(self.params[4])
        N = float(self.params[5])
        
        er = float(self.eps.get().split("(")[1].replace(")",""))
        Ct = er*consts.epsilon_0/tox
        w = (2.24*10**(13))/consts.pi
        Vg0 = consts.elementary_charge*N/Ct
        Cgs = Ct*W*L
        Cgd = (Ct*W*L)/2

        vF = 9.71*10**5 # m/s, from literature, ref 11 in Hu paper
        T = 298 # kelvin

        # Update later/add settings, but for now no back gate
        Cb = 0
        Vb = 0
        Vt0 = 0.8 
        
        K = ((2*consts.elementary_charge**3)/consts.pi)/(consts.hbar*vF)**2

        # Fitting params from paper
        delta = 0.8
        alpha = 0.9
        beta = 0.11
        theta = 7.5
                
        Ids = []

        for i in range(len(self.transVds)):
            Id = []
            Vdsef = self.transVds[i]
            Vds = (1-beta)*Vdsef
            
            Uds = (consts.elementary_charge*Vds)/(consts.Boltzmann*T)

            for j in range(len(self.transVtg)):

                Vt = abs(self.transVtg[j])

                ell = L*(delta/Uds)**alpha
                
                l = (2*mu*consts.k*T)/(theta*consts.elementary_charge*vF)

                r = ell/(ell+l)

                Vch = (np.sqrt((Ct+Cb)**2 + 2*K*(Vt*Ct+Vb*Cb)) - (Ct+Cb))/K

                n0 = (K*Vch**2)/(2*consts.elementary_charge) + N
                
                nF = (consts.elementary_charge*Vch)/(consts.Boltzmann*T)

                integrand = lambda x: x/(np.exp(x-nF) + 1)
                integrand2 = lambda x: x/(np.exp(x-nF+Uds) + 1)

                zeta1 = integrate.quad(integrand, 0, np.inf)
                zeta2 = integrate.quad(integrand2, 0, np.inf)
                
                num = consts.elementary_charge*W*vF*n0*(1-r)*(1-(zeta2[0]/zeta1[0]))
                
                den = 1+r+(1-r)*(zeta2[0]/zeta1[0])

                Id.append(num/den)
            Ids.append(Id)
            
        # Calculate transconductance & transit frequency
        gm =[]
        fT = []
        for entry in Ids:
            gm.append([abs(i/j) for i, j in zip(entry, self.transVtg)])
            res = []
            for entry in gm:
                res = [i/(2*consts.pi*(Cgs+Cgd)) for i in entry]
            fT.append(res)
        return Ids, gm, fT

    def calculateIVChars(self):
        
        tox = float(self.params[0])*10**(-9)
        W = float(self.params[1])*10**(-6)
        L = float(self.params[2])*10**(-6)
        mu = float(self.params[3])
        Ep = float(self.params[4])
        N = float(self.params[5])
        
        er = float(self.eps.get().split("(")[1].replace(")",""))
        Ct = er*consts.epsilon_0/tox
        w = (2.24*10**(13))/consts.pi
        Vg0 = consts.elementary_charge*N/Ct

        vF = 9.71*10**5 # m/s, from literature, ref 11 in Hu paper
        T = 298 # kelvin

        # Update later/add settings, but for now no back gate
        Cb = 0
        Vb = 0
        Vt0 = 0.8 
        
        K = ((2*consts.elementary_charge**3)/consts.pi)/(consts.hbar*vF)**2

        # Fitting params from paper
        delta = 0.8
        alpha = 0.9
        beta = 0.11
        theta = 7.5
                
        Ids = []

        for i in range(len(self.ivVtg)):
            Id = []
            Vt = abs(self.ivVtg[i])

            Vch = (np.sqrt((Ct+Cb)**2 + 2*K*(Vt*Ct+Vb*Cb)) - (Ct+Cb))/K
            
            for j in range(len(self.ivVds)):

                Vdsef = self.ivVds[i]
                Vds = (1-beta)*Vdsef

                if Vds ==0:
                    Id.append(0)
                else:
                    Uds = (consts.elementary_charge*Vds)/(consts.Boltzmann*T)
                    
                    ell = L*(delta/Uds)**alpha
                    
                    l = (2*mu*consts.k*T)/(theta*consts.elementary_charge*vF)

                    r = ell/(ell+l)

                    n0 = (K*Vch**2)/(2*consts.elementary_charge) + N
                    
                    nF = (consts.elementary_charge*Vch)/(consts.Boltzmann*T)

                    integrand = lambda x: x/(np.exp(x-nF) + 1)
                    integrand2 = lambda x: x/(np.exp(x-nF+Uds) + 1)

                    zeta1 = integrate.quad(integrand, 0, np.inf)
                    zeta2 = integrate.quad(integrand2, 0, np.inf)
                    
                    num = consts.elementary_charge*W*vF*n0*(1-r)*(1-(zeta2[0]/zeta1[0]))                    
                    den = 1+r+(1-r)*(zeta2[0]/zeta1[0])
                    Id.append(abs(num/den))
            Ids.append(Id)
        return Ids


#***************************************************************************************#
#                           Dual-Gate Model                                             #
#***************************************************************************************#
class DualGateGFET:

    def __init__(self, params, ivSweep, transSweep, eps):
        self.params = params
        self.ivVds = ivSweep["Vds"]
        self.ivVtg = ivSweep["Vtg"]
        self.ivVbg = ivSweep["Vbg"]
        self.transVds = transSweep["Vds"]
        self.transVtg = transSweep["Vtg"]
        self.transVbg = transSweep["Vbg"]
        self.eps = eps
        self.N = 200 # Num Iterations for self-consisten eqn
        self.gds = []
    
