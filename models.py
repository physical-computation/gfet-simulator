# Class containing code for different GFET models

# For Physical Constants
from scipy import constants as consts

import scipy.integrate as integrate
import scipy.special as special

from numpy.lib.scimath import sqrt as csqrt

# Model of Rodriguez et al.
# Ref.: S. Rodriguez et al., "A Comprehensive Graphene FET Model for Circuit Design,"
# in IEEE Transactions on Electron Devices, vol. 61, no. 4, pp. 1199-1206, April 2014. (doi.org:/10.1109/TED.2014.2302372)

class RodriguezGFET:

    def __init__(self, params, ivSweep, transSweep, eps):
        self.params = params
        self.ivVds = ivSweep["Vds"]
        self.ivVgs = ivSweep["Vgs"]
        self.transVds = transSweep["Vds"]
        self.transVgs = transSweep["Vgs"]
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
            for j in range(len(self.transVgs)):
                Veff = self.transVgs[j] + Vg0
                Id.append(abs((mu*W*Ct*(Veff-0.5*Vds))/(L/Vds +
                            (mu/w)*(csqrt(consts.pi*Ct/consts.elementary_charge))*(csqrt(Veff-0.5*Vds)))))
            Ids.append(Id)
            
        # Calculate transconductance & transit frequency
        gm =[]
        fT = []
        for entry in Ids:
            gm.append([abs(i/j) for i, j in zip(entry, self.transVgs)])
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

        for i in range(len(self.ivVgs)):
            Id = []
            Vgs = self.ivVgs[i]
            for j in range(len(self.ivVds)):
                Veff = Vgs + Vg0
                Vds = self.ivVds[j]
                Id.append(abs((mu*W*Ct*(Veff-0.5*Vds))/(L/Vds +
                            (mu/w)*(csqrt(consts.pi*Ct/consts.elementary_charge))*(csqrt(Veff-0.5*Vds)))))
            Ids.append(Id)
        return Ids   

# Model of Thiele et al. (https://doi.org/10.1063/1.3357398)
# Not currently too accurate, requires solving Cq and Vch self-sonsistently, not sure how
# to do that yet, probably something clever in numpy. may even need parallelising depending
# time taken...

class ThieleGFET:

    def __init__(self, params, ivSweep, transSweep, eps):
        self.params = params
        self.ivVds = ivSweep["Vds"]
        self.ivVgs = ivSweep["Vgs"]
        self.transVds = transSweep["Vds"]
        self.transVgs = transSweep["Vgs"]
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
            for j in range(len(self.transVgs)):
                Vch = self.transVgs[j]-Vg0
                num = (mu*W*(Vch**2)*Vds*consts.elementary_charge**3)/(consts.pi*(consts.hbar*10**6)**2)
                den = L - (mu*Vds/w)*(consts.pi*Vch*(2*Vch*consts.elementary_charge**2)/(consts.pi*(consts.hbar*10**6)**2))**0.5
                Id.append(abs(num/den))
            Ids.append(Id)
            
        # Calculate transconductance & transit frequency
        gm =[]
        fT = []
        for entry in Ids:
            gm.append([abs(i/j) for i, j in zip(entry, self.transVgs)])
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

        for i in range(len(self.ivVgs)):
            Id = []
            Vgs = self.ivVgs[i]
            for j in range(len(self.ivVds)):
                Vch = Vgs-Vg0
                Vds = self.ivVds[i]
                num = (mu*W*(Vch**2)*Vds*consts.elementary_charge**3)/(consts.pi*(consts.hbar*10**6)**2)
                den = L - (mu*Vds/w)*(consts.pi*Vch*(2*Vch*consts.elementary_charge**2)/(consts.pi*(consts.hbar*10**6)**2))**0.5
                Id.append(abs(num/den))
            Ids.append(Id)
        return Ids 
