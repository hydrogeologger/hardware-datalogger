#!/usr/bin/python3
import cmath, math
from sympy import *
from sympy.core import symbol


#                                          Vcc
#                                           |
#                              +------------+
#                              |            |
#                              >            |
#                           R1 <            |
#                              >            |
#                              |          __|D
#                              |         |  |
#                              +--------G|->' Q1
#                              |         |__
#                            [Ic]           |S
#                              |            |
#              Rb             c/            + [V1]
# Vsig --+---/\/\/\--[Ib]--b|<.  Q3         |
#        |                    e\            |
#        |                      |           >
#        |                      |           < R2
#        |                      |           <
#        |                    |GND|         |
#        |                      |           |
#        |                      <           + [V2]
#        |                   R2 >           |
#        |                      <           |
#        |                      |         __|D
#        |                      |        |
#        +----------------------+-------G|<-. Q2
#                                        |__|
#                                           |S
#                                           |
#                                          GND

caps = [1e-12, 2.2e-12, 3.3e-12, 3.9e-12, 4.7e-12, 5.6e-12, 6.8e-12,\
        8.2e-12, 10e-12, 15e-12, 22e-12, 27e-12, 33e-12, 39e-12, 47e-12,\
        56e-12, 82e-12, 100e-12, 150e-12, 180e-12, 220e-12, 330e-12,\
        470e-12, 680e-12, 820e-12, 1e-9, 1.5e-9, 1.8e-9, 2.2e-9, 2.7e-9,\
        3.3e-9, 3.9e-9, 4.7e-9, 5.6e-9, 6.8e-9, 8.2e-9, 10e-9, 15e-9,\
        22e-9, 33e-9, 39e-9, 47e-9, 56e-9, 68e-9, 0.1e-6, 0.22e-6,\
        0.33e-6, 0.47e-6, 0.68e-6, 1e-6, 2.2e-6, 3.3e-6, 4.7e-6, 10e-6,\
        22e-6, 33e-6, 47e-6, 100e-6, 220e-6, 470e-6]

def get_common_resistor_value(Rx):
    Rx = round(Rx)
    resistors = [10, 12, 15, 18, 22, 27, 33, 39, 47, 56, 100, 120, 150, 180,\
                220, 270, 330, 390, 470, 560, 680, 820, 1e3, 1.2e3, 1.5e3,\
                1.8e3, 2.2e3, 2.7e3, 3.3e3, 3.9e3, 4.7e3, 5.6e3, 6.8e3,\
                8.2e3, 10e3, 12e3, 15e3, 18e3, 22e3, 27e3, 33e3, 39e3, 47e3,\
                56e3, 68e3, 82e3, 100e3, 120e3, 150e3, 180e3, 220e3, 270e3,\
                330e3, 390e3, 470e3, 510e3, 560e3, 620e3, 680e3, 820e3, 1e6]
    
    for i in range(0, len(resistors)-1):
        if Rx == resistors[i]:
            resistor_common_value_upper = resistors[i]
            resistor_common_value_lower = resistor_common_value_upper
            break
        elif Rx > resistors[i]:
            continue
        elif Rx < resistors[i]:
            resistor_common_value_lower = resistors[i-1]
            resistor_common_value_upper = resistors[i]
            # lower_p = Rx - resistors[i-1]
            # upper_p = resistors[i] - Rx
            # if lower_p < upper_p:
            #     resistor_common_value_lower = resistors[i-1]
            # else:
            #     resistor_common_value_upper = resistors[i]
            break

    return [math.trunc(resistor_common_value_lower), math.trunc(resistor_common_value_upper), Rx]


def bjt_calculate_ib(vsig, vbe, rb):
    Ib = ((vsig - vbe) / rb)
    return Ib


def bjt_rb_from_ib(Vsig, Vbe, Ib):
    Rb = (Vsig - Vbe) / Ib
    return Rb

def bjt_active_bias_vce(Vce, Vcc, beta, Ib):
    R1 = (Vcc - Vce) / (beta * Ib)
    return R1

def bjt_calculate(_Vbe, Vce_sat, _beta, _Vcc, _Vsig, _Rb, _R1):
    # Assume BJT in Active region
    Ib = bjt_calculate_ib(_Vsig, _Vbe, _Rb)
    Ic = _beta * Ib
    Vce = _Vcc - (Ic * _R1)

    # if (Vce < 0.2):
    if (Vce < Vce_sat):
        # BJT is in Saturation Region
        # (beta * Ib) > Ic > 0
        # Vce = 0.2
        Vce = Vce_sat
        Ic = (_Vcc - Vce) / _R1
    
    Power = (Ic * Vce) + (Ib * _Vbe)
    
    return [Ib, Ic, Vce, Power]

def calculate_theoretical_load_resistor(power, I):
    rload = power / I**2
    return rload


def calculate_fet_triode_vds(Id_rdson, Rds_on):
    # Assume in triode region using Rds
    Vds = Id_rdson * Rds_on
    return Vds

def calculate_fet_K(vgs1, ids1, vgs2, ids2):
    # Obtain from saturation region of output curve
    K = ((math.sqrt(2*ids1) - math.sqrt(2*ids2))/(vgs1 - vgs2))**2
    return K

def calculate_pfet_vsaturation(Vth, K, Rload, Vcc):
    Vgs = symbols("Vgs")
    Id = K/2 * (Vgs - Vth)**2

    # Assume saturation region so Vds <= Vgs - Vth
    # Evaluate for Vds = Vgs - Vth
    Vds = Vgs - Vth
    
    eq1 = Eq(Vds, Id * Rload - Vcc)

    Vsat_solve = solve(eq1, Vgs)
    # Vsat solve returns 2 cases, Vgs > Vth (ignore)
    # Vsat < Vth boundary condition
    for Vsat in Vsat_solve:
        if Vsat < Vth:
            return Vsat
    return Vsat_solve


def calculate_nfet_vsaturation(Vth, K, Rload, Vcc):
    Vgs = symbols("Vgs")
    Id = K/2 * (Vgs - Vth)**2

    # Assume saturation region so Vds >= Vgs - Vth
    # Evaluate for Vds = Vgs - Vth
    Vds = Vgs - Vth
    
    eq1 = Eq(Vds, Id * Rload - Vcc)

    Vsat_solve = solve(eq1, Vgs)

    # Vsat solve returns 2 cases, Vgs < Vth (ignore)
    # Vsat > Vth boundary condition
    for Vsat in Vsat_solve:
        if Vsat > Vth:
            return Vsat
    return Vsat_solve

def calculate_pfet_triode_vds(vgs, vth, K, Id):
    Vds = symbols("Vds")
    eqn = Eq(Id, K * Vds * (vgs - vth - (Vds/2)))
    vds_solve = solve(eqn, Vds)
    
    # Vds > Vgs - Vth to be in triode region for P-Channel fet
    for Vds in vds_solve:
        if Vds > (vgs - vth):
            return Vds
    return vds_solve

def calculate_nfet_triode_vds(vgs, vth, K, Id):
    Vds = symbols("Vds")
    eqn = Eq(Id, K * Vds * (vgs - vth - (Vds/2)))
    vds_solve = solve(eqn, Vds)
    
    
    # Vds < Vgs - Vth to be in triode region for N-Channel fet
    for Vds in vds_solve:
        if Vds < (vgs - vth):
            return Vds
    return vds_solve


def calculate_load_vdrop(vcc, vds_q1, vds_q2):
    return vcc - vds_q2 + vds_q1




def main():    
    # BJT Parameters, BC848B
    Vbe = 0.7
    Vce_sat = 0.09
    beta = 290

    # IRF9z34N and FQP30n06
    P_Vth_irf9z34 = (-2 - 4) / 2
    P_K_irf9z34 = calculate_fet_K(vgs1 = -4.5, ids1 = 2, vgs2 = -8, ids2 = 11.6)
    N_Vth_fqp30n06 = (2 + 4) / 2
    N_K_fqp30n06 = calculate_fet_K(vgs1 = 7, ids1 = 13, vgs2 = 8, ids2 = 15)
    
    # DMC3060LVT - TSOT26 *** Recommended for use
    ## Max Vgs +-12V
    P_Vth_dmc3060lvt = (-0.7 - 2.1)/2 #typical -1.1
    P_K_dmc3060lvt = calculate_fet_K(vgs1 = -1.8, ids1 = 0.5, vgs2 = -3, ids2 = 7.4)
    N_Vth_dmc3060lvt = (0.7 + 1.8)/2 #typical 1
    N_K_dmc3060lvt = calculate_fet_K(vgs1 = 1.8, ids1 = 5, vgs2 = 2.2, ids2 = 13.7)


    # DMC3071LVT - TSOT26
    ## Max Vgs +-20V
    P_Vth_dmc3071lvt = -1.75
    P_K_dmc3071lvt = calculate_fet_K(vgs1 = -2.5, ids1 = 0.5, vgs2 = -4, ids2 = 7)
    N_Vth_dmc3071lvt = 1.75
    N_K_dmc3071lvt = calculate_fet_K(vgs1 = 2.1, ids1 = 0.5, vgs2 = 4, ids2 = 17.5)

    # DMC1028UVT - TSOT26
    ## Max Vgs +-8V
    P_Vth_dmc1028uvt = (-0.4 - 1)/2
    P_K_dmc1028uvt = calculate_fet_K(vgs1 = -1.2, ids1 = 1.7, vgs2 = -2, ids2 = 11.5)
    N_Vth_dmc1028uvt = (0.4 + 1)/2
    N_K_dmc1028uvt = calculate_fet_K(vgs1 = 1.1, ids1 = 0.5, vgs2 = 1.8, ids2 = 17)


    # DMC2038LVT - TSOT26 *** Recommended for use
    ## Max Vgs +-12V
    P_Vth_dmc2038lvt = (-0.4 - 1) / 2
    P_K_dmc2038lvt = calculate_fet_K(vgs1 = -1.5, ids1 = 3, vgs2 = -2.5, ids2 = 13.5)
    N_Vth_dmc2038lvt = (0.4 + 1) / 2
    N_K_dmc2038lvt = calculate_fet_K(vgs1 = 1.5, ids1 = 4, vgs2 = 2.5, ids2 = 22.5)


    # DMC3016LSD - SOIC-8
    ## Max Vgs +-12V
    P_Vth_dmc3016lsd = (-1 - 3) / 2
    P_K_dmc3016lsd = calculate_fet_K(vgs1 = -2.5, ids1 = 1.5, vgs2 = -3, ids2 = 11)
    N_Vth_dmc3016lsd = (1 + 3) / 2
    N_K_dmc3016lsd = calculate_fet_K(vgs1 = 2.2, ids1 = 3.6, vgs2 = 3, ids2 = 27.5)


    Vcc = [5, 14, 12, 3.3]
    Vsig = [5, 3.3]
    
    Rb = 330
    Rb = round(bjt_rb_from_ib(Vsig = 5, Vbe = Vbe, Ib = 0.01e-3))
    print(get_common_resistor_value(Rb))
    Rb = get_common_resistor_value(Rb)[0]

    R1 = 4700
    R1 = round(bjt_active_bias_vce(Vce = 2, Vcc=14, beta=beta, Ib=(bjt_calculate_ib(vsig=5, vbe=Vbe, rb=Rb))))
    print(get_common_resistor_value(R1))
    R1 = get_common_resistor_value(R1)[0]

    load_current = 2.5
    theoreticalRLoad = calculate_theoretical_load_resistor(power = 12.5, I = load_current)
    Vth_q1 = P_Vth_dmc2038lvt
    Kq1 = P_K_dmc2038lvt
    Vth_q2 = N_Vth_dmc2038lvt
    Kq2 = N_K_dmc2038lvt


    print("For Design\tRb={0:n}\tR1={1:n}\t\t Load Current={2}\tLoadImpedance(ohm)={3}".format(Rb, R1, load_current, theoreticalRLoad))
    print("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{8}\t{6}\t{7}\t{9}".format("Vsig(V)", "Vcc(V)", "Ib(A)", "Ic(A)", "Vce(V)", "Power(W)", "Vgs_q1(V)", "Vsat_q1(V)", "Vsat_q2(V)", "VdropLoad(V)"))
    for Vsigx in Vsig:
        for Vccx in Vcc:
            if (Vccx > 5) and (Vsigx == 3.3):
                # Exclude voltage not used for rpi switches
                continue
            elif (Vccx < 5) and (Vsigx == 5):
                # Exclude voltage not used for arduino switches
                continue

            # Calculate BJT characteristics
            [Ib, Ic, Vce, bjtPower] = bjt_calculate(Vbe, Vce_sat, beta, Vccx, Vsigx, Rb, R1)

            vsat_q1 = calculate_pfet_vsaturation(Vth_q1, Kq1, theoreticalRLoad, Vccx)
            vsat_q2 = calculate_nfet_vsaturation(Vth_q2, Kq2, theoreticalRLoad, Vccx)
            Vgs_q1 = (Vce-Vccx)
            vds_q1 = calculate_pfet_triode_vds(Vgs_q1, Vth_q1, Kq1, load_current)
            vds_q2 = calculate_nfet_triode_vds(Vsigx, Vth_q2, Kq2, load_current)
            load_vdrop = calculate_load_vdrop(Vccx, vds_q1, vds_q2)

            print("{0:n}\t{1:n}\t{2:.1e}\t{3:.1e}\t{4:.2}\t{5:.2}\t\t{8:.2f}\t\t{6:.2f}\t\t{7:.2f}\t\t{9:.2f} [{10:.2f}, {11:.2f}]".format(
                Vsigx, Vccx, Ib, Ic, Vce, bjtPower, Vgs_q1, vsat_q1, vsat_q2, load_vdrop, abs(vds_q1), vds_q2))


main()

