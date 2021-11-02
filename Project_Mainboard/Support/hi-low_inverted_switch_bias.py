#!/usr/bin/python3
from switch_bias_functions import *
import fet_bjt_vars

## Circuit 1
#                    Vsc                   Vsc
#                     |                     |
# Vsig --+          __|S                  __|S
#        |         |  |                  |  |
#        |     +--G|->' Q4 (P-Chan)  +--G|->' Q1 (P-Chan)
#        |     |   |__               |   |__
#        |     |      |D             |      |D
#        |     |      +--------------+      |
#        |     |      |                     + [V1]
#        |     |      >                     |
#        |     |      < R2                  >
#        |     |      >                     > RLoad
#        |     |      |        Vcc          <
#        |     |     GND        |           |
#        |     |                <           + [V2]
#        |     |                > R1        |
#        |     |                <         __|D
#        |     |                |        |
#        |     +----------------+-------G|<-. Q2 (N-Chan)
#        |                      |        |__|
#        |                    [Ic]          |S
#        |                      |           |
#        |     Rb           | / c          GND
#        +---/\/\/\--[Ib]--b|    Q3 (NPN)
#                           | > e
#                               |
#                              GND


def main():    
    Vbe_sat = fet_bjt_vars.bc848b.Vbe_sat
    Vce_sat = fet_bjt_vars.bc848b.Vce_sat
    beta = fet_bjt_vars.bc848b.beta

    Vsc = [5, 12, 14.7]
    Vsig = [0, 5]

    Q1 = fet_bjt_vars.dmc3016
    Q2 = fet_bjt_vars.dmc3016
    Q4 = fet_bjt_vars.ntr1p02

    # Q1 = fet_bjt_vars.irf9z34
    # Q2 = fet_bjt_vars.fqp30n06

    if (hasattr(Q1, 'Vth')):
        Vth_q1 = Q1.Vth
    else:
        Vth_q1 = Q1.P.Vth

    if (hasattr(Q2, 'Vth')):
        Vth_q2 = Q2.Vth
    else:
        Vth_q2 = Q2.N.Vth

    if (hasattr(Q1, 'K')):
        Kq1 = Q1.K
    else:
        Kq1 = Q1.P.K

    if (hasattr(Q2, 'K')):
        Kq2 = Q2.K
    else:
        Kq2 = Q2.N.K

    ## Circuit 1
    print ("CIRCUIT 1")
    Rb = 330
    Rb = round(npn_rb_from_ib(Vbb = 5, Vbe = Vbe_sat, Ib = 0.01e-3))
    print("Rb options: {0}".format(get_common_resistor_value(Rb)))
    Rb = get_common_resistor_value(Rb)[0]

    R1 = 3900
    Ib = npn_calculate_ib(Vbb=5, Vbe=Vbe_sat, Rb=Rb)
    R1 = round(npn_active_bias_vce(desired_vce = 0, Vce_sat = Vce_sat, Vcc = 14.7, beta = beta, Ib = Ib, Ic = 0.1e-3))
    print("R1 options: {0}".format(get_common_resistor_value(R1)))
    R1 = get_common_resistor_value(R1)[0]

    active_load = True # Load is passive (false): pure resistive or active (true): acts as current-stable non-linear resistor
    load_current = 2.5
    theoreticalRLoad = calculate_theoretical_load_resistor(power = 12.5, I = load_current)


    # print("Fet\tP(Q1): {0}\tN(Q2): {1}".format(namestr(Kq1, globals())[0][4:].upper(), namestr(Kq2, globals())[0][4:].upper() ))
    print("Fet\tP(Q1): {0}\tN(Q2): {1}".format(Q1.name.upper(), Q2.name.upper() ))
    print("For Design\tRb(Ohm)={0:2}\t\tR1(Ohm)={1}\t\tDesired Current(A)={2}\tLoadImpedance(ohm)={3:n} [{4}]".format(
            toSI(Rb), toSI(R1), load_current, theoreticalRLoad, "active" if active_load else "passive" ))
    # print("{0:7s} {1:7s} {2:8s} {3:8s} {4:6s} {5:7s}   {6:8s} {7:9s}   {8:8s} {9:9s}   {10:5s} [{11}{13},{12}{14}] {15:8s}".format(
    #         "Vsig(V)","Vsc(V)", "Ib(A)", "Ic(A)", "Vce(V)", "P_Q3(W)", "Vgs_Q1(V)", "Vsat_Q1(V)", "Vgs_Q2(V)", "Vsat_Q2(V)", "Vdrop(V) Load", "Q1", "Q2","","", "I_load(A)"))
    print("{0:7s} {1:6s} {2:^8s} {3:^8s} {4:6s} {5:7s} {6} {7} {8}".format(
             "Vsig(V)","Vsc(V)", "Ib(A)", "Ic(A)", "Vce(V)", "P_Q3(W)", "Q4Vsat(V)", "Q4Vgs(V)", "Q4Vds(V)"))

    for Vsigx in Vsig:
        for Vscx in Vsc:
            # Calculate BJT characteristics
            [Ib, Ic, Vce, _,_,_, bjtPower] = npn_calculate(Vbe=Vbe_sat, Vce_sat=Vce_sat, beta=beta, Vcc=Vscx, Vbb=Vsigx, Rb=Rb, Rc=R1, Re=0)

            Vsat_q4 = calculate_pfet_vsaturation(Q4.Vth, Q4.K, 120e3, Vscx)
            Vgs_q4 = (Vce - Vscx)
            [Id_q4, Vds_q4, Region_q4] = pfet_calculate(Vgs = Vgs_q4, Vth = Q4.Vth, K = Q4.K, Isd = 0.1225e-3, Vss = Vscx, Rd = 120e3)

            # Vsat_q1 = calculate_pfet_vsaturation(Vth_q1, Kq1, theoreticalRLoad, Vscx)
            # Vsat_q2 = calculate_nfet_vsaturation(Vth_q2, Kq2, theoreticalRLoad, Vscx)

            print("{0:^7n} {1:^6n} {2:<8.1e} {3:<8.1e} {4:^6.2f} {5:^7.3f} {6:^9.2f} {7:>6.2f} {8} {9:8.2f}".format(
                    Vsigx, Vscx, Ib, Ic, Vce, bjtPower, Vsat_q4, Vgs_q4, Region_q4[0], Vds_q4))

            # Vgs_q1 = (Vsigx - Vscx)
            # [Id_q1, Vds_q1, Region_q1] = pfet_calculate(Vgs = Vgs_q1, Vth = Vth_q1, K = Kq1, Isd = load_current, Vss = Vscx, Rd = theoreticalRLoad)

            # Vgs_q2 = Vce
            # [Id_q2, Vds_q2, Region_q2] = nfet_calculate(Vgs = Vgs_q2, Vth = Vth_q2, K = Kq2, Ids = load_current, Vdd = Vscx, Rd = theoreticalRLoad)

            # # Find fet which is limiting the current of system (lowest current)
            # # and estimate load current and voltage drop
            # if (Id_q1 < Id_q2):
            #     # Q1 limiting system
            #     [Id_q2, Vds_q2, Region_q2] = nfet_calculate(Vgs = Vgs_q2, Vth = Vth_q2, K = Kq2, Ids = Id_q1, Vdd = Vscx, Rd = theoreticalRLoad)
            #     Id = Id_q1
            #     load_vdrop = Id * theoreticalRLoad
            # elif (Id_q2 < Id_q1):
            #     # Q2 Limiting system
            #     [Id_q1, Vds_q1, Region_q1] = pfet_calculate(Vgs = Vgs_q1, Vth = Vth_q1, K = Kq1, Isd = Id_q2, Vss = Vscx, Rd = theoreticalRLoad)
            #     Id = Id_q2
            #     load_vdrop = Id * theoreticalRLoad
            # elif (Region_q1[0] == "c" or Region_q2[0] == "c"):
            #     # Either Q1 or Q2 in cutoff region
            #     Id = 0
            #     load_vdrop = Id * theoreticalRLoad
            # else:
            #     # Q1 and Q2 in triode region
            #     load_vdrop = Vscx - Vds_q2 + Vds_q1
            #     if (active_load):
            #         Id = load_current
            #     else:
            #         Id = load_vdrop / theoreticalRLoad

            # print("{0:^7n} {1:>4n}   {2:<8.1e} {3:<8.1e}  {4:<6.2f} {5:<7.3f}  {6: 7.2f} {13} {7: 8.2f}\t  {8:^6.2f} {14}  {9:^9.2f}    {10:>5.2f} [{11: 5.2f}, {12:>4.2f}]\t{15:.3f}".format(
            #         Vsigx, Vscx, Ib, Ic, Vce, bjtPower, Vgs_q1, Vsat_q1, Vgs_q2, Vsat_q2, load_vdrop, Vds_q1, Vds_q2, Region_q1[0], Region_q2[0], Id))

    print("\r\n")

    print("{0:7s} {1:6s} {2} {3}\t{4} {5}\t{6} [{7},{8}] {9}".format(
             "Vsig(V)","Vsc(V)", "Q1Vsat(V)", "Q1Vgs(V)", "Q2Vsat(V)", "Q2Vgs(V)" ,"Vdrop(V) Load", "Q1", "Q2", "I_load(A)"))

    for Vsigx in Vsig:
        for Vscx in Vsc:
            # Calculate BJT characteristics
            [Ib, Ic, Vce, _,_,_, bjtPower] = npn_calculate(Vbe=Vbe_sat, Vce_sat=Vce_sat, beta=beta, Vcc=Vscx, Vbb=Vsigx, Rb=Rb, Rc=R1, Re=0)

            Vsat_q4 = calculate_pfet_vsaturation(Q4.Vth, Q4.K, 120e3, Vscx)
            Vgs_q4 = (Vce - Vscx)
            [Id_q4, Vds_q4, Region_q4] = pfet_calculate(Vgs = Vgs_q4, Vth = Q4.Vth, K = Q4.K, Isd = 0.1225e-3, Vss = Vscx, Rd = 120e3)

            Vsat_q1 = calculate_pfet_vsaturation(Vth_q1, Kq1, theoreticalRLoad, Vscx)
            Vsat_q2 = calculate_nfet_vsaturation(Vth_q2, Kq2, theoreticalRLoad, Vscx)
            
            if (Region_q4[0] == "c"):
                Vgs_q1 = (0 - Vscx)
            else:
                Vgs_q1 = (Vscx + Vds_q4 - Vscx)
            [Id_q1, Vds_q1, Region_q1] = pfet_calculate(Vgs = Vgs_q1, Vth = Vth_q1, K = Kq1, Isd = load_current, Vss = Vscx, Rd = theoreticalRLoad)

            Vgs_q2 = Vce
            [Id_q2, Vds_q2, Region_q2] = nfet_calculate(Vgs = Vgs_q2, Vth = Vth_q2, K = Kq2, Ids = load_current, Vdd = Vscx, Rd = theoreticalRLoad)

            # Find fet which is limiting the current of system (lowest current)
            # and estimate load current and voltage drop
            if (Id_q1 < Id_q2):
                # Q1 limiting system
                [Id_q2, Vds_q2, Region_q2] = nfet_calculate(Vgs = Vgs_q2, Vth = Vth_q2, K = Kq2, Ids = Id_q1, Vdd = Vscx, Rd = theoreticalRLoad)
                Id = Id_q1
                load_vdrop = Id * theoreticalRLoad
            elif (Id_q2 < Id_q1):
                # Q2 Limiting system
                [Id_q1, Vds_q1, Region_q1] = pfet_calculate(Vgs = Vgs_q1, Vth = Vth_q1, K = Kq1, Isd = Id_q2, Vss = Vscx, Rd = theoreticalRLoad)
                Id = Id_q2
                load_vdrop = Id * theoreticalRLoad
            elif (Region_q1[0] == "c" or Region_q2[0] == "c"):
                # Either Q1 or Q2 in cutoff region
                Id = 0
                load_vdrop = Id * theoreticalRLoad
            else:
                # Q1 and Q2 in triode region
                load_vdrop = Vscx - Vds_q2 + Vds_q1
                if (active_load):
                    Id = load_current
                else:
                    Id = load_vdrop / theoreticalRLoad

            print("{0:^7n} {1:^6n} {2:^9.2f} {3:>6.2f} {4}\t{5:^9.2f} {6:>6.2f} {7}\t{8:5.2f} [{9:5.2f},{10:>4.2f}]    {11:^9.3f}".format(
                    Vsigx, Vscx, Vsat_q1, Vgs_q1, Region_q1[0], Vsat_q2, Vgs_q2, Region_q2[0], load_vdrop, Vds_q1, Vds_q2, Id))


main()