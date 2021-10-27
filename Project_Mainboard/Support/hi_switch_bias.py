#!/usr/bin/python3
from switch_bias_functions import *
import fet_bjt_vars

# P-Fet in common drain configuration

# Circuit 1: Normal Logic
#                              V1            Vsc
#                               |            |
#                               >            |
#                            R1 <            |
#                               >            |
#                               |          __|D
#                               |         |  |
#                               +--------G|->' Q1 (P-Chan)
#                               |         |__
#                             [Ic]           |S
#                               |            |
#              Rb           | / c            + [V2]
# Vsig --+--/\/\/\--[Ib]---b|     Q2 (NPN)   |
#                           | > e            >
#                               |            < Rload
#                               |            >
#                               |            |
#                             |GND|----------+

# Circuit 2: Inverted Logic (Normal P Channel)
#
#                             V1            Vsc
#                              |            |
#                              >            |
#                           R2 <            |
#                              >          __|S
#                              |         |  |
# Vsig -+----------------------+--------G|->' Q1 (P-Chan)
#                                        |__
#                                           |D
#                                           |
#                                           + [V2]
#                                           |
#                                           >
#                                           < RLoad
#                                           <
#                                           |
#                                          GND



def main():    
    Vbe_sat = fet_bjt_vars.bc848b.Vbe_sat
    Vce_sat = fet_bjt_vars.bc848b.Vce_sat
    beta = fet_bjt_vars.bc848b.beta

    Vsc = [5.1]
    Vsig = [0, 5]

    Q1 = fet_bjt_vars.dmp3028

    # Q1 = fet_bjt_vars.irf9z34

    if (hasattr(Q1, 'Vth')):
        Vth_q1 = Q1.Vth
    else:
        Vth_q1 = Q1.P.Vth

    if (hasattr(Q1, 'K')):
        Kq1 = Q1.K
    else:
        Kq1 = Q1.P.K


    load_power = 12.5
    active_load = True # Load is passive (false): pure resistive or active (true): acts as current-stable non-linear resistor
    load_current = 2.5

    print("Fet\tP(Q1): {0}".format(Q1.name.upper()))

    ## Circuit 1
    print ("CIRCUIT 1")
    Rb = 330
    Rb = round(npn_rb_from_ib(Vbb = 5, Vbe = Vbe_sat, Ib = 0.01e-3))
    print("Rb options: {0}".format(get_common_resistor_value(Rb)))
    Rb = get_common_resistor_value(Rb)[0]

    R1 = 3900
    Ib = npn_calculate_ib(Vbb=5, Vbe=Vbe_sat, Rb=Rb)
    R1 = round(npn_active_bias_vce(desired_vce = 0, Vce_sat = Vce_sat, Vcc = 14.7, beta = beta, Ib = Ib))
    print("R1 options: {0}".format(get_common_resistor_value(R1)))
    R1 = get_common_resistor_value(R1)[1]

    theoreticalRLoad = calculate_theoretical_load_resistor(power = load_power, I = load_current)

    print("For Design\tRb(Ohm)={0:2}\t\tR1(Ohm)={1}\t\tDesired Current(A)={2}\tLoadImpedance(ohm)={3:n} [{4}]".format(
            toSI(Rb), toSI(R1), load_current, theoreticalRLoad, "active" if active_load else "passive" ))
    print("{0:7s} {1:7s} {2:8s} {3:8s} {4:6s} {5:7s}   {6:8s} {7:9s}   {8:8s} [{9}{10}] {11:8s}".format(
            "Vsig(V)","Vsc(V)", "Ib(A)", "Ic(A)", "Vce(V)", "P_Q2(W)", "Vgs_Q1(V)", "Vsat_Q1(V)", "Vdrop(V) Load", "Q1","", "I_load(A)"))

    for Vsigx in Vsig:
        for Vscx in Vsc:
            V1 = Vscx
            # Calculate BJT characteristics
            [Ib, Ic, Vce, _,_,_, bjtPower] = npn_calculate(Vbe=Vbe_sat, Vce_sat=Vce_sat, beta=beta, Vcc=V1, Vbb=Vsigx, Rb=Rb, Rc=R1, Re=0)

            Vsat_q1 = calculate_pfet_vsaturation(Vth_q1, Kq1, theoreticalRLoad, Vscx)

            Vgs_q1 = (Vce - Vscx)
            [_, Vds_q1, Region_q1] = pfet_calculate(Vgs = Vgs_q1, Vth = Vth_q1, K = Kq1, Isd = load_current, Vss = Vscx, Rd = theoreticalRLoad)


            # Find fet which is limiting the current of system (lowest current)
            # and estimate load current and voltage drop
            if (Region_q1[0] == "c"):
                # Either Q1 or Q2 in cutoff region
                Id = 0
                load_vdrop = Id * theoreticalRLoad
            else:
                # Q1 and Q2 in triode region
                load_vdrop = Vscx + Vds_q1
                if (active_load):
                    Id = load_current
                else:
                    Id = load_vdrop / theoreticalRLoad

            print("{0:^7n} {1:>4n}   {2:<8.1e} {3:<8.1e}  {4:<6.2f} {5:<7.3f}  {6:7.2f} {10:2} {7:^10.2f}\t{8:>10.2f} [{9:5.2f}]   {11:.3f}".format(
                    Vsigx, Vscx, Ib, Ic, Vce, bjtPower, Vgs_q1, Vsat_q1, load_vdrop, Vds_q1, Region_q1[0], Id))

    print("\r\n")

    ## Circuit 2
    print ("CIRCUIT 2")
    
    # Estimate pullup R2 Resistor
    R2 = calculate_theoretical_load_resistor(power = (V1 * 0.1e-3), I = 0.1e-3)
    print("R2 options: {0}".format(get_common_resistor_value(R2)))
    R2 = get_common_resistor_value(R2)[1]
    
    print("For Design\tR2(Ohm)={0:2}\t\tDesired Current(A)={1}\tLoadImpedance(ohm)={2:n} [{3}]".format(
            toSI(R2), load_current, theoreticalRLoad, "active" if active_load else "passive" ))
    print("{0:7s} {1:6s} {2:7s} {3:10s} {4} {5} {6} {7} [{8}] {9}".format(
            "Vsig(V)","Vsc(V)", "I_r2(A)", "P_r2(Watt)", "Vgs_Q1(V)", "", "Vsat_Q1(V)", "Vdrop(V) Load", "Q1", "I_load(A)"))
    for Vsigx in Vsig:
        for Vscx in Vsc:
            Ir2 = (Vscx - Vsigx) / R2
            P_r2 = R2 * Ir2**2

            Vsat_q1 = calculate_pfet_vsaturation(Vth_q1, Kq1, theoreticalRLoad, Vscx)

            Vgs_q1 = (Vsigx - Vscx)
            [_, Vds_q1, Region_q1] = pfet_calculate(Vgs = Vgs_q1, Vth = Vth_q1, K = Kq1, Isd = load_current, Vss = Vscx, Rd = theoreticalRLoad)

            # Find fet which is limiting the current of system (lowest current)
            # and estimate load current and voltage drop
            if (Region_q1[0] == "c"):
                # Either Q1 or Q2 in cutoff region
                Id = 0
                load_vdrop = Id * theoreticalRLoad
            else:
                # Q1 and Q2 in triode region
                load_vdrop = Vscx + Vds_q1
                if (active_load):
                    Id = load_current
                else:
                    Id = load_vdrop / theoreticalRLoad

            print("{0:^7n} {1:^6n} {2:<7.1e} {3:<10.3f} {4:>6.2f} {5}\t{6:.2f}\t{7:>8.2f} [{8:5.2f}]    {9:.3f}".format(
                    Vsigx, Vscx, Ir2, P_r2, Vgs_q1, Region_q1[0], Vsat_q1, load_vdrop, Vds_q1, Id))


main()