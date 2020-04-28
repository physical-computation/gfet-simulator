# Class containing code for the GFET model derived by Rodriguez et al.
#
# Ref.: S. Rodriguez et al., "A Comprehensive Graphene FET Model for Circuit Design,"
# in IEEE Transactions on Electron Devices, vol. 61, no. 4, pp. 1199-1206, April 2014. (doi.org:/10.1109/TED.2014.2302372)

# For Physical Constants
from scipy import constants as consts
from numpy.lib.scimath import sqrt as csqrt

class RodriguezGFET:

    def __init__(self, params, VdsSweep, VgsSweep, datapoints, eps):
        self.params = params
        self.Vds = VdsSweep
        self.Vgs = VgsSweep
        self.datapoints = datapoints
        self.eps = eps
    
    def calculateIds(self):
        
        tox = float(self.params[0])*10**(-9)
        W = float(self.params[1])*10**(-6)
        L = float(self.params[2])*10**(-6)
        mu = float(self.params[3])
        Ep = float(self.params[4])
        N = float(self.params[5])
        Vg0 = float(self.params[6])
        
        er = float(self.eps.get().split("(")[1].replace(")",""))

        Ct = er*consts.epsilon_0/tox
        w = (2.24*10**(13))/consts.pi
                
        Id = []
        
        for i in range(self.datapoints):
            Veff = self.Vgs[i] + Vg0
            Id.append(abs((mu*W*Ct*(Veff-0.5*self.Vds))/(L/self.Vds + (mu/w)*(csqrt(consts.pi*Ct/consts.elementary_charge))*(csqrt(Veff-0.5*self.Vds)))))

        return {"Ids": Id}

class ThieleGFET:

    def __init__(self, params, VdsSweep, VgsSweep, datapoints, eps):
        self.params = params
        self.Vds = VdsSweep
        self.Vgs = VgsSweep
        self.datapoints = datapoints
        self.eps = eps
    
    def calculateIds(self):
        
        tox = float(self.params[0])*10**(-9)
        W = float(self.params[1])*10**(-6)
        L = float(self.params[2])*10**(-6)
        mu = float(self.params[3])
        Ep = float(self.params[4])
        N = float(self.params[5])
        Vg0 = float(self.params[6])
        
        er = float(self.eps.get().split("(")[1].replace(")",""))

        Ct = er*consts.epsilon_0/tox
        w = (2.24*10**(13))/consts.pi
                
        Id = []
        
        for i in range(self.datapoints):
            Vch = self.Vgs[i]-Vg0
            num = (mu*W*(Vch**2)*self.Vds*consts.elementary_charge**3)/(consts.pi*(consts.hbar*10**6)**2)
            den = L - (mu*self.Vds/w)*(consts.pi*Vch*(2*Vch*consts.elementary_charge**2)/(consts.pi*(consts.hbar*10**6)**2))**0.5
            Id.append(abs(num/den))

        return {"Ids": Id}
