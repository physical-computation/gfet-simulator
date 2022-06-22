# Class containing code for different GFET models

# For Physical Constants
from scipy import constants as consts
from scipy.misc import derivative
import scipy.integrate as integrate
import scipy.special as special
from scipy import optimize
from numpy.lib.scimath import sqrt as csqrt
from numpy.linalg import norm
import numpy as np
import cmath
import re
import math

# Model of Mukherjee et al.
# Ref.: C. Mukherjee, J. Aguirre-Morales, S. Frégonèse, T. Zimmer and C. Maneux,
# "Versatile Compact Model for Graphene FET Targeting Reliability-Aware Circuit Design,"
# in IEEE Transactions on Electron Devices, vol. 62, no. 3, pp. 757-763, March 2015, doi: 10.1109/TED.2015.2395134.

class MukherjeeGFET:

    def __init__(self, params, ivSweep, transSweep, eps):

        # Voltage Sweep Data
        self.ivVds = ivSweep["Vds"]
        self.ivVtg = ivSweep["Vtg"]
        self.ivVbg = ivSweep["Vbg"]
        self.transVds = transSweep["Vds"]
        self.transVtg = transSweep["Vtg"]
        self.transVbg = transSweep["Vbg"]

        # Device Geometry and Parameters

        # Air/vacuum means infinite thickness
        if "Air" in eps[0].get():
            tox1 = np.inf
        else:
            tox1 = float(params[0])#*10**(-9) #nm to m

        if "Air" in eps[1].get():
            tox2 = np.inf
        else:
            tox2 = float(params[1])#*10**(-9) #nm to m

        self.W = float(params[2])#*10**(-6) #um to m
        self.L = float(params[3])#*10**(-6) #um to m

        self.mun = float(params[4])#*10**(-4) # cm2/Vs to m2/vs
        self.mup = float(params[5])#*10**(-4) # cm2/Vs to m2/vs

        self.vF = float(params[6]) # m/s
        self.Nf = float(params[7])#*10**(-4) # m-2

        Ep = float(re.search('(Ɛr=(.*), ħω=(.*))', eps[0].get()).group(3)[:-5])*consts.e*10**(-3) # J
        er1 = float(re.search('(Ɛr=(.*), ħω=(.*))', eps[0].get()).group(2))
        er2 = float(re.search('(Ɛr=(.*), ħω=(.*))', eps[1].get()).group(2))
        T = float(params[8]) # K

        if tox1 == 0:
            self.Ct = 0 # F/m2
        else:
            self.Ct = (er1*consts.epsilon_0)/tox1 # F/m2
        if tox2 == 0:
            self.Cb = 0 # F/m2
        else:
            self.Cb = (er2*consts.epsilon_0)/tox2 # F/m2

        self.omega = Ep/consts.hbar # s-1

        Wgr = 4.57 # graph work func (eV to J)
        Wau = 5.1 # Au contact work func (eV to J)
        Phib = (Wau - Wgr)*consts.e # J

        # r0 is measured characteristic. Can assume same, and extract vals from literature
        # for ballpark figs? Jing Tian QMUL 2017 PhD Thesis suggests max r0 of 2.5M Ohm - um
        # Rc taken from Urban 20202 Contact res & mobility (~2kR-um)

        rd0 = (2.5*10**-6)*10**(-4) # M Ohm-um order of magnitude cancels to Ohm-m
        rs0 = (2.5*10**-6)*10**(-4) # M Ohm-um order of magnitude cancels to Ohm-m
        Rc = (2*10**3)*10**(-4) # k Ohm-um to Ohm-m
        self.Rd = (Rc + rd0*np.exp((consts.e*Phib)/(consts.Boltzmann*T))) # Ohm-m
        self.Rs = (Rc + rs0*np.exp((consts.e*Phib)/(consts.Boltzmann*T))) # Ohm-m

        # Delta and npud from W. Zhu et al. "Carrier scattering, mobilities and electrostatic potential..." (2009)
        #delta = 54*10**(-3)*consts.e # meV to J
        #self.npud = (delta**2)/(consts.pi*((consts.hbar*self.vF)**2))#(2/(consts.pi*((consts.hbar*self.vF)**2)))*(((delta**2)/2) + ((consts.pi**2)/6)*(consts.Boltzmann**2)*(T**2))
        #self.beta = (consts.e**3)/(consts.pi*((consts.hbar*self.vF)**2)) # C^3/J^2-m^2

        # Avoid division by zero for infinite thickness
        if self.Ct == 0:
            self.Vg0 = 0
        else:
            self.Vg0 = 0.85#(consts.e*self.Nf)/self.Ct

        if self.Cb == 0:
            self.Vb0 = 0
        else:
            self.Vb0 = 0.85#(consts.e*self.Nf)/self.Cb

    # Calculate channel voltage
    def fnVch(self, S, sgn):
        return sgn*((-(self.Ct+self.Cb)+np.sqrt(((self.Ct+self.Cb)**2)+4*sgn*self.beta*S))/(2*self.beta)) # V

    # Calculate saturation velocity
    def fnVsat(self, Qnet):
        return self.omega/np.sqrt(consts.pi*abs(Qnet)/consts.e+self.npud/2) # m/s

    # Calculate p-type or n-type region drain current. 1/2 L as only half the gate for each region
    def fnIds(self, mu, Qnet, Vdsi):
        return abs(mu)*self.W*(abs(Qnet)*Vdsi+(consts.e*self.npud*Vdsi)/2)/(self.L+abs(mu)*(Vdsi/self.fnVsat(Qnet)))

    # Calculate device transfer characteristics (Dirac pt sweep)
    def calculateTransferChars(self):
        Vb = self.transVbg - self.Vb0
        Ids = []

        for Vds in self.transVds:
            Id = []
            #Id_last = 0

            for v in self.transVtg:
                #Vg = v - self.Vg0

                #Vdsi = Vds #- Id_last*(self.Rs+self.Rd)

                #S = self.Ct*Vg + self.Cb*Vb + consts.e*self.Nf - (self.Ct+self.Cb)*(Vdsi/2)
                #sgn = np.sign(S)
                #Vch = self.fnVch(S, sgn)

            #    if sgn <= 0:
            #        QnetAVn = 0
            #        QnetAVp = self.beta*Vch*abs(Vch) # C/m2
                #elif sgn > 0:
                #    QnetAVn = self.beta*Vch*abs(Vch) # C/m2
                #    QnetAVp = 0
                #print(QnetAVn)
                #Id_last = self.fnIds(self.mun, QnetAVn, Vdsi) + self.fnIds(self.mup, QnetAVp, Vdsi)
                Id.append(Id_last)

            Ids.append(Id)

        # Calculate transconductance & transit frequency
        Cgs = self.Ct*self.W*self.L
        Cbs = self.Cb*self.W*self.L
        Cgd = (self.Ct*self.W*self.L)/2

        gm =[]
        fT = []
        for entry in Ids:
            gm.append(np.ones(len(entry)))#[abs(i/j) for i, j in zip(entry, self.transVtg)])
            res = []
            for entry in gm:
                res = np.ones(len(entry))#[i/(2*consts.pi*(Cgs+Cgd)) for i in entry]
            fT.append(res)
        return Ids, gm, fT

    def calculateIVChars(self):

        Vb = self.ivVbg - self.Vb0
        Ids = []

        for Vtg in self.ivVtg:
            Id = []
            Vg = Vtg - self.Vg0
            Id_last = 0

            for Vds in self.ivVds:
                Vdsi = Vds - Id_last*(self.Rs+self.Rd)

                S = (self.Ct*Vg)+(self.Cb*Vb)+(consts.e*self.Nf)-((self.Ct+self.Cb)*(Vdsi/2))
                sgn = np.sign(S)

                Vch = self.fnVch(S, sgn)

                if sgn <= 0:
                    QnetAVn = 0
                    QnetAVp = self.beta*Vch*abs(Vch) # C/m2
                elif sgn > 0:
                    QnetAVn = self.beta*Vch*abs(Vch) # C/m2
                    QnetAVp = 0

                Id_last = self.fnIds(self.mun, QnetAVn, Vdsi) + self.fnIds(self.mup, QnetAVp, Vdsi)
                Id.append(Id_last)
            Ids.append(Id)

        return Ids

# GFET Based on Jimenez et al. (https://ieeexplore.ieee.org/document/6054021)
class JimenezGFET:

    def __init__(self, params, ivSweep, transSweep, eps):
        # Voltage Sweep Data
        self.ivVds = ivSweep["Vds"]
        self.ivVtg = ivSweep["Vtg"]
        self.ivVbg = ivSweep["Vbg"]
        self.transVds = transSweep["Vds"]
        self.transVtg = transSweep["Vtg"]
        self.transVbg = transSweep["Vbg"]

        # Device Geometry and Parameters
        # Air/vacuum means infinite thickness
        if "Air" in eps[0].get():
            tox1 = np.inf
        elif "None" in eps[0].get():
            tox1 = 1*10**(-9)
        else:
            tox1 = float(params[0])*10**(-9) #nm to m

        if "Air" in eps[1].get():
            tox2 = np.inf
        elif "None" in eps[1].get():
            tox2 = 1*10**(-9)
        else:
            tox2 = float(params[1])*10**(-9) #nm to m, hacky extra division by 10, but works

        self.W = float(params[2])*10**(-6) #um to m
        self.L = float(params[3])*10**(-6) #um to m
        self.mun = float(params[4])*10**(-4) # cm2/Vs to m2/vs
        self.mup = float(params[5])*10**(-4) # cm2/Vs to m2/vs
        self.vF = float(params[6]) # m/s
        self.Nf = float(params[7]) # m-2

        Ep = float(re.search('(Ɛr=(.*), ħω=(.*))', eps[0].get()).group(3)[:-5])*10**(-3)*consts.e # J
        er1 = float(re.search('(Ɛr=(.*), ħω=(.*))', eps[0].get()).group(2))
        er2 = float(re.search('(Ɛr=(.*), ħω=(.*))', eps[1].get()).group(2))
        T = float(params[8]) # K
        self.T = T

        self.omega = Ep/consts.hbar # s-1
        delta = 54*10**(-3)*consts.e # meV to J
        self.npud = (delta**2)/(consts.pi*((consts.hbar*self.vF)**2))
        #self.npud = (2/(consts.pi*((consts.hbar*self.vF)**2)))*(((delta**2)/2) + ((consts.pi**2)/6)*(consts.Boltzmann**2)*(T**2))
        self.beta = (consts.e**3)/(consts.pi*((consts.hbar*self.vF)**2)) # C^3/J^2-m^2

        # No capacitance if no oxide(?)
        if tox1 == 0:
            self.Ct = 0 # F/m2
        else:
            self.Ct = (er1*consts.epsilon_0)/tox1 # F/m2
        if tox2 == 0:
            self.Cb = 0 # F/m2
        else:
            self.Cb = (er2*consts.epsilon_0)/tox2 # F/m2

        #Wgr = 4.57 # graph work func (eV to J)
        #Wau = 5.1 # Au contact work func (eV to J)
        #Phib = (Wau - Wgr)*consts.e # J

        # r0 is measured characteristic. Can assume same, and extract vals from literature
        # for ballpark figs? Jing Tian QMUL 2017 PhD Thesis suggests max r0 of 2.5M Ohm - um
        # Rc taken from Urban 20202 Contact res & mobility (~2kR-um)

        #rd0 = (2.5*10**6)*10**(-4) # M Ohm-um order of magnitude cancels to Ohm-m
        #rs0 = (2.5*10**6)*10**(-4) # M Ohm-um order of magnitude cancels to Ohm-m
        #Rc = (2*10**3)*10**(-4) # k Ohm-um to Ohm-m
        self.Rd = 300#*10**(-6)#(Rc + rd0*np.exp((consts.e*Phib)/(consts.Boltzmann*T))) # Ohm-m
        self.Rs = 300#*10**(-6)#(Rc + rs0*np.exp((consts.e*Phib)/(consts.Boltzmann*T))) # Ohm-m

        self.k = ((2*(consts.e**2))/consts.pi)*(consts.e/((consts.hbar*self.vF)**2))

        # Avoid division by zero for infinite thickness
        if self.Ct == 0:
            self.Vg0 = 0
        else:
            self.Vg0 = (consts.e*self.Nf)/self.Ct

        if self.Cb == 0:
            self.Vb0 = 0
        else:
            self.Vb0 = (consts.e*self.Nf)/self.Cb

    def fnVcd(self, Vtg, Vds, Vbg):
        alpha = self.mun/self.mup
        if ((Vtg-self.Vg0-Vds)*self.Ct+(Vbg-self.Vb0-Vds)*self.Cb) > 0:
            return alpha*(-(self.Ct+self.Cb)+np.sqrt(((self.Ct+self.Cb)**2) + 2*self.k*((Vtg-self.Vg0-Vds)*self.Ct+(Vbg-self.Vb0-Vds)*self.Cb)))/(self.k)
        else:
            return (-(self.Ct+self.Cb)+np.sqrt(((self.Ct+self.Cb)**2) - 2*self.k*((Vtg-self.Vg0-Vds)*self.Ct+(Vbg-self.Vb0-Vds)*self.Cb)))/(-self.k)

    def fnVcs(self, Vtg, Vbg):
        alpha = self.mun/self.mup
        if ((Vtg-self.Vg0)*self.Ct+(Vbg-self.Vb0)*self.Cb) > 0:
            return alpha*(-(self.Ct+self.Cb)+np.sqrt(((self.Ct+self.Cb)**2) + 2*self.k*((Vtg-self.Vg0)*self.Ct+(Vbg-self.Vb0)*self.Cb)))/(self.k)
        else:
            return (-(self.Ct+self.Cb)+np.sqrt(((self.Ct+self.Cb)**2) - 2*self.k*((Vtg-self.Vg0)*self.Ct+(Vbg-self.Vb0)*self.Cb)))/(-self.k)

    def fnG(self, Vc):
        return (-(Vc**3)/3)-np.sign(Vc)*(self.k*(Vc**4)/4*(self.Ct+self.Cb))

    # Calculate saturation velocity
    def fnVsat(self, Qnet):
        return self.omega/np.sqrt(consts.pi*abs(Qnet)/consts.e+self.npud/2) # m/s

    def fnIds(self, Vtg, Vbg, Vds):

        Vch = self.fnVcd(Vtg, Vds, Vbg) - self.fnVcs(Vtg, Vbg)
        QnetAV = self.beta*Vch*abs(Vch) # C/m2
        vSat = self.fnVsat(QnetAV)
        muAv = (self.mup + self.mun)/2
        Leff = self.L+muAv*(abs(Vds)/vSat)

        return ((muAv*self.k)/2)*(self.W/Leff)*(self.fnG(self.fnVcd(Vtg, Vds, Vbg)) - self.fnG(self.fnVcs(Vtg, Vbg)))

    def calculateTransferChars(self):
        Ids = []
        for Vds in self.transVds:
            Id = []
            for Vtg in self.transVtg:
                Idx = self.fnIds(Vtg, self.transVbg, Vds)
                Id.append(Idx)
            Ids.append(Id)

        # Calculate transconductance & transit frequency
        Vg = []
        for i, v in enumerate(self.transVtg):
            Vg.append(v+Ids[0][i]*self.Rs)

        gm =[]
        fT = []
        for Id in Ids:
            gm.append([abs(i/j) for i, j in zip(Id, self.transVtg)])
            res = []
            for entry in gm:
                res = [i/(2*consts.pi*self.W*self.L*(self.Ct+self.Cb)) for i in entry]
            fT.append(res)
        return Ids, gm, fT, Vg

    def calculateIVChars(self):
        Ids = []
        for Vtg in self.ivVtg:
            Ids.append([self.fnIds(Vtg, self.ivVbg, Vds) for Vds in self.ivVds])
        return Ids


# Model of Hu et al. (https://doi.org/10.1063/1.3357398) - UNUSED AT THE MOMENT

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
        N = float(self.params[6])
        T = float(self.params[7])

        Ep = float(re.search('(Ɛr=(.*), ħω=(.*))', self.eps[0].get()).group(3)[:-5])*10**(-3) # Conversion from (m)eV to Joules
        er1 = float(re.search('(Ɛr=(.*), ħω=(.*))', self.eps[0].get()).group(2))
        Ep2 = float(re.search('(Ɛr=(.*), ħω=(.*))', self.eps[1].get()).group(3)[:-5])*10**(-3) # Conversion from (m)eV to Joules
        er2 = float(re.search('(Ɛr=(.*), ħω=(.*))', self.eps[1].get()).group(2))

        Ct = (er1*consts.epsilon_0)/tox1
        Cb = (er2*consts.epsilon_0)/tox2
        Vg0 = (consts.e*N)/Ct
        Cgs = Ct*W*L
        Cgd = (Ct*W*L)/2

        # Update later/add settings, but for now no back gate
        Vbg = self.transVbg # test, assumes one step atm

        K = (2*consts.e**3)/(consts.pi*(consts.hbar*vF)**2)

        # Fitting params from paper
        delta = 0.8
        alpha = 1.2
        beta = 0.4
        theta = 1

        Ids = []

        for i in range(len(self.transVds)):
            Id = []
            Vdsef = self.transVds[i]
            Vds = (1-beta)*Vdsef

            Uds = (consts.e*Vds)/(consts.Boltzmann*T)

            for j in range(len(self.transVtg)):

                Vt = self.transVtg[j] - Vg0

                # Equation breaks down at Vds=0, so assume ell=L
                # See: A. Rahman and M. S. Lundstrom, "A compact scattering model
                # for the nanoscale double-gate MOSFET," in IEEE Transactions on
                # Electron Devices, vol. 49, no. 3, pp. 481-489, March 2002
                if Uds == 0:
                    ell = L
                else:
                    ell = L*(delta/Uds)**alpha

                l = (2*mu*consts.k*T)/(theta*consts.elementary_charge*vF)

                r = ell/(ell+l)

                Vch = (np.sqrt(abs((Ct+Cb)**2 + 2*K*(Vt*Ct+Vbg*Cb))) - (Ct+Cb))/K

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
        N = float(self.params[6])
        T = float(self.params[7])
        Ep = float(re.search('(Ɛr=(.*), ħω=(.*))', self.eps[0].get()).group(3)[:-5])*10**(-3) # Conversion from (m)eV to Joules
        er1 = float(re.search('(Ɛr=(.*), ħω=(.*))', self.eps[0].get()).group(2))
        Ep2 = float(re.search('(Ɛr=(.*), ħω=(.*))', self.eps[1].get()).group(3)[:-5])*10**(-3) # Conversion from (m)eV to Joules
        er2 = float(re.search('(Ɛr=(.*), ħω=(.*))', self.eps[1].get()).group(2))

        Ct = (er1*consts.epsilon_0)/tox1
        Cb = (er2*consts.epsilon_0)/tox2
        Vg0 = (consts.elementary_charge*N)/Ct

        K = (2*consts.elementary_charge**3)/(consts.pi*(consts.hbar*vF)**2)

        # Fitting params from paper
        delta = 0.8
        alpha = 0.9
        beta = 0.11
        theta = 7.5

        Ids = []
        Vbg = self.transVbg # assumes one step atm

        for i in range(len(self.ivVtg)):
            Id = []
            Vt = self.ivVtg[i]
            Vch = ((np.sqrt((Ct+Cb)**2 + 2*K*(Vt*Ct+Vbg*Cb))) - (Ct+Cb))/K
            # should be q*Vch - Ec (cond band energy
            nF = (consts.elementary_charge*Vch)/(consts.Boltzmann*T)

            for j in range(len(self.ivVds)):

                Vdsef = self.ivVds[j]
                Vds = (1-beta)*Vdsef

                Uds = (consts.elementary_charge*Vds)/(consts.Boltzmann*T)


                # Equation breaks down at Vds=0, so assume ell=L
                # See: A. Rahman and M. S. Lundstrom, "A compact scattering model
                # for the nanoscale double-gate MOSFET," in IEEE Transactions on
                # Electron Devices, vol. 49, no. 3, pp. 481-489, March 2002
                if Uds == 0:
                    ell = L
                else:
                    ell = L*(delta/Uds)**alpha

                l = (2*mu*consts.k*T)/(theta*consts.elementary_charge*vF)
                r = ell/(ell+l)
                n0 = ((K*(Vch**2))/(2*consts.elementary_charge)) + N

                integrand = lambda x: x/(np.exp(x-nF) + 1)
                integrand2 = lambda x: x/(np.exp(x-nF+Uds) + 1)

                zeta1 = integrate.quad(integrand, 0, np.inf)
                zeta2 = integrate.quad(integrand2, 0, np.inf)

                num = consts.elementary_charge*W*vF*n0*(1-r)*(1-(zeta2[0]/zeta1[0]))
                den = 1+r+(1-r)*(zeta2[0]/zeta1[0])

                Id.append(num/den)
            Ids.append(Id)
        return Ids

#***************************************************************************************#
#          Model of Thiele et al. (https://doi.org/10.1063/1.3357398)                   #
#                            NOT WORKING PROPERLY AT PRESENT!                           #
#***************************************************************************************#
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
        self.N = 500 # Num Iterations for self-consisten eqn
        self.gds = []

    def calculateTransferChars(self):

        # self-consistent eqn solving using fixed-point iteration method
        def fixedp(f, x0, tol, N):
            e = 1
            xp = []
            itr = 0

            while(e>tol and itr<N):
                x = f(x0)
                e = norm(x0-x)
                x0 = x
                xp.append(x0)
                itr = itr+1
            return x, xp

        # Integral for denominator, can't integrate directly a power of a lambda
        # function, i.e., 1/vsat, so have done it here
        def integrand(Vx, rhosh, omega):
            A = 10**(-3)
            return ((consts.pi*rhosh(Vx))**(0.5+A*(Vx**2))/omega)

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
            Vch = -5 # Vch is initially the first voltage

            for j in range(len(self.transVtg)):
                Vtg = self.transVtg[j] - Vg0

                # Cq in terms of Vch
                f = lambda Vch: ((2*(consts.elementary_charge**2))/(consts.pi))*((abs(Vch)*consts.elementary_charge)/(consts.hbar*vF)**2)
                g = lambda Vx: ((Ct/(Ct+Cb+0.5*Cq))*(Vtg-Vx)+(Cb/(Ct+Cb+0.5*Cq))*(Vbg-Vx))

                e1 = 0.001 # error, Volts
                e2 = 10**(-9) # error, F/m2
                Cq, Varr = fixedp(f, Vch, e1, self.N)
                Vch = g(Vds)-g(0)

                # Eqs 14 and 15 in Thiele paper
                rhosh = lambda Vx: abs(-0.5*Cq*(((Ct/(Ct+Cb+0.5*Cq))*(Vtg-Vx)+(Cb/(Ct+Cb+0.5*Cq))*(Vbg-Vx))))/consts.elementary_charge
                num = consts.elementary_charge*mu*W*integrate.quad(rhosh, 0, Vds)[0] # integral returns value and error as 2 element array, just want value
                den = L - mu*integrate.quad(integrand, 0, Vds, args=(rhosh, omega))[0] # as above

                Id.append(num/den)
            Ids.append(Id)

        # Calculate transconductance & transit frequency
        gm =[]
        fT = []

        for index, entry in enumerate(Ids):
            gm.append([abs(i/j) for i, j in zip(entry, self.transVtg)])
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
        # Cq in terms of Vch
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

        Ct = (er1*consts.epsilon_0)/tox1
        Cb = (er2*consts.epsilon_0)/tox2
        Vg0 = (consts.elementary_charge*N)/Ct
        omega = Ep/consts.hbar
        Vbg = self.transVbg # assumes one step atm
        Ids = []

        for i in range(len(self.ivVtg)):
            Id = []
            Vtg = self.ivVtg[i]
            Vch = Vtg - Vg0 # Vch is initially the first voltage

            for j in range(len(self.ivVds)):
                Vds = self.ivVds[j]

#                if Vds == 0:
#                    Id.append(0)
#                else:
                e = 0.001 # error/Volts
                Cq, Varr = fixedp(Vch, e, self.N)
                #Vch = (Vtg-Vds)*(Ct/(Ct+0.5*Cq))
                rhosh = lambda Vx: abs(-0.5*Cq*(((Ct/(Ct+Cb+0.5*Cq))*(Vtg-Vx)+(Cb/(Ct+Cb+0.5*Cq))*(Vbg-Vx))))/consts.elementary_charge

                num = consts.elementary_charge*mu*W*integrate.quad(rhosh, 0 ,Vds)[0]

                den = L - mu*(np.sqrt(consts.pi)/omega)*integrate.quad(integrand, 0, Vds, args=(rhosh))[0]

                Id.append(abs(num/den))
            Ids.append(Id)

        for index, entry in enumerate(Ids):
            dIds = max(entry)-min(entry)
            dVds = max(self.ivVds)-min(self.ivVds)
            self.gds.append(dIds/dVds)
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
