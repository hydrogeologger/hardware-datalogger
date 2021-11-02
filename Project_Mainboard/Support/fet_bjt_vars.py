#!/usr/bin/python3
from sympy.core.numbers import NaN
from switch_bias_functions import calculate_fet_K, nfet_calculate

class Bjt:
    def __init__(self, Vbe_sat, Vce_sat, beta, name=""):
        self.name = name
        self.Vbe_sat = Vbe_sat
        self.Vce_sat = Vce_sat
        self.beta = beta


class FetCharacteristics:
    def __init__(self):
        self.K = NaN
        self.Vth = NaN

    def set_k(self, Vgs1, Id1, Vgs2, Id2):
        self.K = calculate_fet_K(Vgs1=Vgs1,Id1=Id1,Vgs2=Vgs2,Id2=Id2)
    
    def set_vth(self, min, max=None):
        if (max is None):
            self.Vth = min
        else:
            self.Vth = (min + max) / 2
    

class fet(FetCharacteristics):
    def __init__(self, name):
        FetCharacteristics.__init__(self)
        self.name = name
        self.P = FetCharacteristics()
        self.N = FetCharacteristics()
    
    def set_vth(self, min, max=None):
        if (hasattr(self, 'N')):
            del(self.N)
        if (hasattr(self, 'P')):
            del(self.P)
        super().set_vth(min, max)
    
    def set_k(self, Vgs1, Id1, Vgs2, Id2):
        if (hasattr(self, 'N')):
            del(self.N)
        if (hasattr(self, 'P')):
            del(self.P)
        FetCharacteristics.set_k(self, Vgs1, Id1, Vgs2, Id2)

    def n_set_vth(self, min, max=None):
        if (hasattr(self, 'Vth')):
            del(self.Vth)
        self.N.set_vth(min, max)
    
    def p_set_vth(self, min, max=None):
        if (hasattr(self, 'Vth')):
            del(self.Vth)
        self.P.set_vth(min, max)

    def n_set_k(self, Vgs1, Id1, Vgs2, Id2):
        if (hasattr(self, 'K')):
            del(self.K)
        self.N.set_k(Vgs1, Id1, Vgs2, Id2)
    
    def p_set_k(self, Vgs1, Id1, Vgs2, Id2):
        if (hasattr(self, 'K')):
            del(self.K)
        self.P.set_k(Vgs1, Id1, Vgs2, Id2)


# BJT Parameters, BC848B
bc848b = Bjt(Vbe_sat=0.7,Vce_sat=0.09,beta=290)

# NTR1P02/NVR1P02 P-Fet, SOT-23
## Max Vgs -20V, Rds = 0.148 Ohm @-10V
ntr1p02 = fet("NTR1P02")
ntr1p02.set_vth((-1.1-2.3)/2)
ntr1p02.set_k(Vgs1 = -2.5, Id1 = 0.125, Vgs2 = -3.5, Id2 = 1.75)

# IRF9z34N and FQP30n06
## IRF9z34N Max Vgs +-20
## FQP30n06 Max Vgs +-25
P_Vth_irf9z34 = (-2 - 4) / 2
P_K_irf9z34 = calculate_fet_K(Vgs1 = -4.5, Id1 = 2, Vgs2 = -8, Id2 = 11.6)
N_Vth_fqp30n06 = (2 + 4) / 2
N_K_fqp30n06 = calculate_fet_K(Vgs1 = 7, Id1 = 13, Vgs2 = 8, Id2 = 15)
irf9z34 = fet("irf9z34")
irf9z34.set_vth(P_Vth_irf9z34)
irf9z34.set_k(Vgs1 = -4.5, Id1 = 2, Vgs2 = -8, Id2 = 11.6)
fqp30n06 = fet("fqp30n06")
fqp30n06.set_vth(N_Vth_fqp30n06)
fqp30n06.set_k(Vgs1 = 7, Id1 = 13, Vgs2 = 8, Id2 = 15)


# DMC3060LVT - TSOT26 *** Recommended for use
## Max Vgs +-12V
P_Vth_dmc3060lvt = (-0.7 - 2.1)/2 #typical -1.1
P_K_dmc3060lvt = calculate_fet_K(Vgs1 = -1.8, Id1 = 0.5, Vgs2 = -3, Id2 = 7.4)
N_Vth_dmc3060lvt = (0.7 + 1.8)/2 #typical 1
N_K_dmc3060lvt = calculate_fet_K(Vgs1 = 1.8, Id1 = 5, Vgs2 = 2.2, Id2 = 13.7)


# DMC3071LVT - TSOT26
## Max Vgs +-20V
P_Vth_dmc3071lvt = -1.75
P_K_dmc3071lvt = calculate_fet_K(Vgs1 = -2.5, Id1 = 0.5, Vgs2 = -4, Id2 = 7)
N_Vth_dmc3071lvt = 1.75
N_K_dmc3071lvt = calculate_fet_K(Vgs1 = 2.1, Id1 = 0.5, Vgs2 = 4, Id2 = 17.5)

# DMC1028UVT - TSOT26
## Max Vgs +-8V
P_Vth_dmc1028uvt = (-0.4 - 1)/2
P_K_dmc1028uvt = calculate_fet_K(Vgs1 = -1.2, Id1 = 1.7, Vgs2 = -2, Id2 = 11.5)
N_Vth_dmc1028uvt = (0.4 + 1)/2
N_K_dmc1028uvt = calculate_fet_K(Vgs1 = 1.1, Id1 = 0.5, Vgs2 = 1.8, Id2 = 17)


# DMC2038LVT - TSOT26 *** Recommended for use
## Max Vgs +-12V
P_Vth_dmc2038lvt = (-0.4 - 1) / 2
P_K_dmc2038lvt = calculate_fet_K(Vgs1 = -1.5, Id1 = 3, Vgs2 = -2.5, Id2 = 13.5)
N_Vth_dmc2038lvt = (0.4 + 1) / 2
N_K_dmc2038lvt = calculate_fet_K(Vgs1 = 1.5, Id1 = 4, Vgs2 = 2.5, Id2 = 22.5)


# DMC3016LSD - SOIC-8
## Max Vgs +-20V
P_Vth_dmc3016lsd = (-1 - 3) / 2
P_K_dmc3016lsd = calculate_fet_K(Vgs1 = -2.5, Id1 = 1.5, Vgs2 = -3, Id2 = 11)
N_Vth_dmc3016lsd = (1 + 3) / 2
N_K_dmc3016lsd = calculate_fet_K(Vgs1 = 2.2, Id1 = 3.6, Vgs2 = 3, Id2 = 27.5)
dmc3016 = fet("dmc3016lsd")
dmc3016.p_set_vth(P_Vth_dmc3016lsd)
dmc3016.p_set_k(Vgs1 = -2.5, Id1 = 1.5, Vgs2 = -3, Id2 = 11)
dmc3016.n_set_vth(N_Vth_dmc3016lsd)
dmc3016.n_set_k(Vgs1 = 2.2, Id1 = 3.6, Vgs2 = 3, Id2 = 27.5)


 # IRF9362PbF - SOIC-8, P-Fet, 2-Chan
## Max Vgs +-20V, Rds = 17mOhm
irf9362 = fet("irf9362")
irf9362.set_vth(-1.8)
irf9362.set_k(Vgs1 = -2.5, Id1 = -0.15, Vgs2 = -4.5, Id2 = -13)

# DMP3028LSD - SOIC-8, P-Fet, 2-Chan
## Max Vgs +-20V, Rds_max(-4.5Vgs) = 38mOhm
dmp3028 = fet("dmp3028lsd")
dmp3028.set_vth((-1- 3) /  2)
dmp3028.set_k(Vgs1 = -2.5, Id1 = -1.5, Vgs2 = -3.5, Id2 = -20)

# Si4909DY - SOIC-8, P-Fet, 2-Chan
## Max Vgs +-20V, Rds_max(-4.5Vgs) = 34mOhm
si4909 = fet("si4909dy")
si4909.set_vth((-1.2- 2.5) /  2)
si4909.set_k(Vgs1 = -3, Id1 = -10, Vgs2 = -4, Id2 = -40)