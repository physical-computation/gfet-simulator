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
        # No back gate in this model, to skip bottom oxide
        tox = float(self.params[0])*10**(-9)
        W = float(self.params[2])*10**(-6)
        L = float(self.params[3])*10**(-6)
        mu = float(self.params[4])*10**(-4)
        Ep = float(self.params[6])*consts.elementary_charge
        N = float(self.params[7])
        er = float(self.eps[0].get().split("(")[1].replace(")",""))
        Ct = (er*consts.epsilon_0)/tox
        omega = Ep/consts.hbar
        Vg0 = (consts.elementary_charge*N)/Ct
        Cgs = Ct*W*L
        Cgd = (Ct*W*L)/2
        
        Ids = []

        for i in range(len(self.transVds)):
            Id = []
            Vds = self.transVds[i]
            for j in range(len(self.transVtg)):
                Veff = self.transVtg[j] + Vg0
                Id.append(abs((mu*W*Ct*(Veff-0.5*Vds))/((L/Vds) +
                            (mu/omega)*(np.sqrt((consts.pi*Ct)/consts.elementary_charge))*(np.sqrt(abs(Veff-0.5*Vds))))))
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
        # No back gate in this model, to skip bottom oxide 
        tox = float(self.params[0])*10**(-9)
        W = float(self.params[2])*10**(-6)
        L = float(self.params[3])*10**(-6)
        mu = float(self.params[4])*10**(-4)
        Ep = float(self.params[6])*consts.elementary_charge
        N = float(self.params[7])
        
        er = float(self.eps[0].get().split("(")[1].replace(")",""))
        Ct = (er*consts.epsilon_0)/tox
        omega = Ep/consts.hbar
        Vg0 = (consts.elementary_charge*N)/Ct

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
                    Id.append(abs((mu*W*Ct*(Veff-0.5*Vds))/
                            ((L/Vds) + (mu/omega)*(np.sqrt((consts.pi*Ct)/consts.elementary_charge))*(np.sqrt(abs(Veff-0.5*Vds))))))
            Ids.append(Id)
        return Ids   

# GFET Based on Jimenez et al. (https://ieeexplore.ieee.org/document/6054021)
class JimenezGFET:

    def __init__(self, params, ivSweep, transSweep, eps):
        self.params = params
        self.ivVds = ivSweep["Vds"]
        self.ivVtg = ivSweep["Vtg"]
        self.transVds = transSweep["Vds"]
        self.transVtg = transSweep["Vtg"]
        self.transVbg = transSweep["Vbg"]
        self.eps = eps
    
    def calculateTransferChars(self):
        tox1 = float(self.params[0])*10**(-9)
        tox2 = float(self.params[1])*10**(-9)
        W = float(self.params[2])*10**(-6)
        L = float(self.params[3])*10**(-6)
        mu = float(self.params[4])*10**(-4)
        vF = float(self.params[5])
        Ep = float(self.params[6])*consts.elementary_charge
        N = float(self.params[7])
        er1 = float(self.eps[0].get().split("(")[1].replace(")",""))
        er2 = float(self.eps[1].get().split("(")[1].replace(")",""))
        
        eox1= er1*consts.epsilon_0
        eox2= er2*consts.epsilon_0
        Ct = eox1/tox1
        Cb = eox2/tox2
        k = ((2*consts.elementary_charge**2)/consts.pi)*(consts.elementary_charge/((consts.hbar*vF)**2))
        Vbg = self.transVbg
        Vg0 = -2
        Vb0 = 0
        
        Ids = []

        for i in range(len(self.transVds)):
            Id = []
            Vds = self.transVds[i]
            for j in range(len(self.transVtg)):

                v = self.transVtg[j]
                
                #Vcd
                if ((v-Vg0-Vds)*Ct+(Vbg-Vb0-Vds)*Cb) > 0:
                    Vcd = (-(Ct+Cb)+np.sqrt(((Ct+Cb)**2)+2*k*((v-Vg0-Vds)*Ct+(Vbg-Vb0-Vds)*Cb)))/(k)
                else:
                    Vcd = (-(Ct+Cb)+np.sqrt(((Ct+Cb)**2)-2*k*((v-Vg0-Vds)*Ct+(Vbg-Vb0-Vds)*Cb)))/(-k)
                
                #Vcs
                if ((v-Vg0)*Ct+(Vbg-Vb0)*Cb) > 0:
                    Vcs = (-(Ct+Cb)+np.sqrt(((Ct+Cb)**2)+2*k*((v-Vg0)*Ct+(Vbg-Vb0)*Cb)))/(k)
                else:
                    Vcs = (-(Ct+Cb)+np.sqrt(((Ct+Cb)**2)-2*k*((v-Vg0)*Ct+(Vbg-Vb0)*Cb)))/(-k)
                    
                Leff = float(L+mu*(abs(Vds)/vF))
                g1 = float(-(Vcd**(3))/3)-np.sign(Vcd)*((k*(Vcd**4))/4*(Ct+Cb))
                g2 = float(-(Vcs**(3))/3)-np.sign(Vcs)*((k*(Vcs**4))/4*(Ct+Cb))
                
                Id.append(((mu*k)/2)*(W/Leff)*(g1-g2))
                
            Ids.append(Id)
            
        # Calculate transconductance & transit frequency
        gm =[]
        fT = []
        for entry in Ids:
            gm.append([abs(i/j) for i, j in zip(entry, self.transVtg)])
            res = []
            for entry in gm:
                res = [i/(2*consts.pi*(Ct+Cb)) for i in entry]
            fT.append(res)
        return Ids, gm, fT

    def calculateIVChars(self):
        tox1 = float(self.params[0])*10**(-9)
        tox2 = float(self.params[1])*10**(-9)
        W = float(self.params[2])*10**(-6)
        L = float(self.params[3])*10**(-6)
        mu = float(self.params[4])*10**(-4)
        vF = float(self.params[5])
        Ep = float(self.params[6])*consts.elementary_charge
        N = float(self.params[7])
        er1 = float(self.eps[0].get().split("(")[1].replace(")",""))
        er2 = float(self.eps[1].get().split("(")[1].replace(")",""))

        eox1= er1*consts.epsilon_0
        eox2= er2*consts.epsilon_0
        Ct = eox1/tox1
        Cb = eox2/tox2
        k = ((2*consts.elementary_charge**2)/consts.pi)*(consts.elementary_charge/((consts.hbar*vF)**2))
        Vbg = self.transVbg
        Vg0 = 0.85
        Vb0 = 0

        Ids = []

        for i in range(len(self.ivVtg)):
            Id = []
            Vtg = self.ivVtg[i]
            for j in range(len(self.ivVds)):

                Vds = self.ivVds[j]
                
                #Vcd
                if ((Vtg-Vg0-Vds)*Ct+(Vbg-Vb0-Vds)*Cb) > 0:
                    Vcd = (-(Ct+Cb)+np.sqrt(((Ct+Cb)**2)+2*k*((Vtg-Vg0-Vds)*Ct+(Vbg-Vb0-Vds)*Cb)))/(k)
                else:
                    Vcd = (-(Ct+Cb)+np.sqrt(((Ct+Cb)**2)-2*k*((Vtg-Vg0-Vds)*Ct+(Vbg-Vb0-Vds)*Cb)))/(-k)
                
                #Vcs
                if ((Vtg-Vg0)*Ct+(Vbg-Vb0)*Cb) > 0:
                    Vcs = (-(Ct+Cb)+np.sqrt(((Ct+Cb)**2)+2*k*((Vtg-Vg0)*Ct+(Vbg-Vb0)*Cb)))/(k)
                else:
                    Vcs = (-(Ct+Cb)+np.sqrt(((Ct+Cb)**2)-2*k*((Vtg-Vg0)*Ct+(Vbg-Vb0)*Cb)))/(-k)
                    
                Leff = float(L+mu*(abs(Vds)/vF))
                g1 = float(-(Vcd**3)/3)-np.sign(Vcd)*((k*(Vcd**4))/4*(Ct+Cb))
                g2 = float(-(Vcs**3)/3)-np.sign(Vcs)*((k*(Vcs**4))/4*(Ct+Cb))
                
                Id.append(((mu*k)/2)*(W/Leff)*(g1-g2))
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

        # Cq in terms of Vch
        def f(Vch):
            return abs(Vch)*(2*consts.e**3)/(consts.pi*(consts.hbar*vF)**2)

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
        # i.e. 1/vsat
        def integrand(Vx, rho, omega):
            A = 10**(-3)
            return ((consts.pi*rhosh(Vx))**(0.5+A*Vx**2))/omega 
        
        tox1 = float(self.params[0])*10**(-9)
        tox2 = float(self.params[1])*10**(-9)
        W = float(self.params[2])*10**(-6)
        L = float(self.params[3])*10**(-6)
        mu = float(self.params[4])*10**(-4)
        vF = float(self.params[5])
        Ep = float(self.params[6])*consts.elementary_charge
        N = float(self.params[7])
        er1 = float(self.eps[0].get().split("(")[1].replace(")",""))
        er2 = float(self.eps[1].get().split("(")[1].replace(")",""))
        
        er1 = float(self.eps[0].get().split("(")[1].replace(")",""))
        er2 = float(self.eps[1].get().split("(")[1].replace(")",""))
        Ct = (er1*consts.epsilon_0)/tox1
        Cb = (er2*consts.epsilon_0)/tox2
        Vg0 = (consts.elementary_charge*N)/Ct
        Ctg = Ct*W*L
        Cbg = Cb*W*L
        Cgd = (Ct*W*L)/2
        omega = Ep/consts.hbar

        Vbg = self.transVbg

        Ids = []

        for i in range(len(self.transVds)):
            Id = []
            Vds = self.transVds[i]
            Vch = self.transVtg[0] - Vg0 # Vch is initially the first voltage
            for j in range(len(self.transVtg)):
                Vtg = self.transVtg[j]
                
                e = 0.001 # error/Volts
                Cq, Varr = fixedp(Vch, e, self.N)

                # Eqs 14 and 15 in Thiele paper
                rhosh = lambda Vx: abs(-0.5*Cq*((Ct*(Vtg-Vx)+Cb*(Vbg-Vx))/(Ctg+Cbg+0.5*Cq)))/consts.e
                num = consts.e*mu*W*integrate.quad(rhosh, 0 ,Vds)[0] # integral returns value and error, just want value
                den = L - mu*integrate.quad(integrand, 0, Vds, args=(rhosh, omega))[0] # as above
                
                Id.append(abs(num/den))
            Ids.append(Id)

        # Calculate transconductance & transit frequency
        gm =[]
        fT = []

        for index, entry in enumerate(Ids):
            gm.append([abs(i/j) for i, j in zip(entry, self.transVtg)])
#            gd = self.gds[index]
            res = []
            for entry in gm:
#                res = [i/(2*consts.pi*((Ctg+Cgd)*(1+gd*(Rs+Rd)+Cgd*i*(Rs+Rd)))) for i in entry]
                res = [i/(2*consts.pi*(Ctg+Cgd)) for i in entry] # above calculation has an error, use simple
                                                                 # for now.
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
        
        tox1 = float(self.params[0])*10**(-9)
        tox2 = float(self.params[1])*10**(-9)
        W = float(self.params[2])*10**(-6)
        L = float(self.params[3])*10**(-6)
        mu = float(self.params[4])*10**(-4)
        vF = float(self.params[5])
        Ep = float(self.params[6])*consts.elementary_charge
        N = float(self.params[7])
        er1 = float(self.eps[0].get().split("(")[1].replace(")",""))
        er2 = float(self.eps[1].get().split("(")[1].replace(")",""))
        
        Ct = er1*consts.epsilon_0/tox1
        Cb = er2*consts.epsilon_0/tox2
        vF = 10**8 # m/s
        Ctg = Ct*W*L
        Cbg = Cb*W*L
        Vg0 = consts.elementary_charge*N/Ctg
        omega = Ep/consts.hbar
        Vbg = self.transVbg # assumes one step atm
        Ids = []
        
        for i in range(len(self.ivVtg)):
            Id = []
            Vtg = self.ivVtg[i]
            Vch = Vtg - Vg0 # Vch is initially the first voltage
            
            for j in range(len(self.ivVds)):
                Vds = self.ivVds[j]
                
                if Vds == 0:
                    Id.append(0)
                else:
                    e = 0.001 # error/Volts
                    Cq, Varr = fixedp(Vch, e, self.N)
                    Vch = (Vtg-Vds)*(Ctg/(Ctg+0.5*Cq))
                    rhosh = lambda Vx: abs(-0.5*Cq*((Ct*(Vtg-Vx)+Cb*(Vbg-Vx))/(Ctg+Cbg+0.5*Cq)))/consts.elementary_charge

                    num = consts.elementary_charge*mu*W*integrate.quad(rhosh, 0 ,Vds)[0]

                    den = L - mu*(np.sqrt(consts.pi)/omega)*integrate.quad(integrand, 0, Vds, args=(rhosh))[0]
                    
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
        
        tox1 = float(self.params[0])*10**(-9)
        tox2 = float(self.params[1])*10**(-9)
        W = float(self.params[2])*10**(-6)
        L = float(self.params[3])*10**(-6)
        mu = float(self.params[4])*10**(-4)
        vF = float(self.params[5])
        Ep = float(self.params[6])*consts.elementary_charge
        N = float(self.params[7])
        T = float(self.params[8])
        er1 = float(self.eps[0].get().split("(")[1].replace(")",""))
        er2 = float(self.eps[1].get().split("(")[1].replace(")",""))

        Ct = er1*consts.epsilon_0/tox1
        Cb = er2*consts.epsilon_0/tox2
        Vg0 = consts.elementary_charge*N/Ct
        Cgs = Ct*W*L
        Cgd = (Ct*W*L)/2

        # Update later/add settings, but for now no back gate
        Vbg = self.transVbg # test, assumes one step atm
        
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

                Vch = (np.sqrt((Ct+Cb)**2 + 2*K*(Vt*Ct+Vbg*Cb)) - (Ct+Cb))/K

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
        
        tox1 = float(self.params[0])*10**(-9)
        tox2 = float(self.params[1])*10**(-9)
        W = float(self.params[2])*10**(-6)
        L = float(self.params[3])*10**(-6)
        mu = float(self.params[4])*10**(-4)
        vF = float(self.params[5])
        Ep = float(self.params[6])*consts.elementary_charge
        N = float(self.params[7])
        T = float(self.params[8])
        er1 = float(self.eps[0].get().split("(")[1].replace(")",""))
        er2 = float(self.eps[1].get().split("(")[1].replace(")",""))
        
        Ct = er1*consts.epsilon_0/tox1
        Cb = er2*consts.epsilon_0/tox2
        Vg0 = consts.elementary_charge*N/Ct

        K = ((2*consts.elementary_charge**3)/consts.pi)/((consts.hbar*vF)**2)
        
        # Fitting params from paper
        delta = 0.8
        alpha = 0.9
        beta = 0.11
        theta = 7.5
                
        Ids = []
        Vbg = self.transVbg # assumes one step atm
        
        for i in range(len(self.ivVtg)):
            Id = []
            Vt = abs(self.ivVtg[i])
            Vch = ((np.sqrt((Ct+Cb)**2 + 2*K*(Vt*Ct+Vbg*Cb))) - (Ct+Cb))/K
            
            for j in range(len(self.ivVds)):

                Vdsef = self.ivVds[i]
                Vds = (1-beta)*Vdsef

                if Vds == 0:
                    Id.append(0)
                else:
                    Uds = (consts.elementary_charge*Vds)/(consts.Boltzmann*T)
                    
                    ell = L*((delta/Uds)**alpha)
                    
                    l = (2*mu*consts.k*T)/(theta*consts.elementary_charge*vF)

                    r = ell/(ell+l)

                    n0 = ((K*(Vch**2))/(2*consts.elementary_charge)) + N

                    # should be q*Vch - Ec (cond band energy
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
#                           Dual-Gate Model - TBC                                       #
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
    
