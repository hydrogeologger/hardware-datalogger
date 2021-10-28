#!/usr/bin/python3
# import cmath, math
from math import floor, log10, nan, sqrt, trunc
from sympy import Eq, solve, symbols

# User defined Exception
class BetaValueError(Exception):
    """Raised when beta value is non-positive or not given for calculation"""
    pass

class CurrentValueError(Exception):
    """Raised when current value is not appropriate for given calculation"""
    pass

caps = [1e-12, 2.2e-12, 3.3e-12, 3.9e-12, 4.7e-12, 5.6e-12, 6.8e-12,\
        8.2e-12, 10e-12, 15e-12, 22e-12, 27e-12, 33e-12, 39e-12, 47e-12,\
        56e-12, 82e-12, 100e-12, 150e-12, 180e-12, 220e-12, 330e-12,\
        470e-12, 680e-12, 820e-12, 1e-9, 1.5e-9, 1.8e-9, 2.2e-9, 2.7e-9,\
        3.3e-9, 3.9e-9, 4.7e-9, 5.6e-9, 6.8e-9, 8.2e-9, 10e-9, 15e-9,\
        22e-9, 33e-9, 39e-9, 47e-9, 56e-9, 68e-9, 0.1e-6, 0.22e-6,\
        0.33e-6, 0.47e-6, 0.68e-6, 1e-6, 2.2e-6, 3.3e-6, 4.7e-6, 10e-6,\
        22e-6, 33e-6, 47e-6, 100e-6, 220e-6, 470e-6]


def toSI(d, unit="", digits=2):
    #Pretty print for SI numbers
    #d: number the print
    #unit: SI unit to use
    #digits: level of precision (number of decimal digits)
    incPrefixes = ['k', 'M', 'G', 'T', 'P', 'E', 'Z', 'Y']
    decPrefixes = ['m', 'µ', 'n', 'p', 'f', 'a', 'z', 'y']
    if d != 0:
        degree = int(floor(log10(abs(d)) / 3))
    else:
        degree = 0
    prefix = ''
    if degree != 0:
        if degree > 0:
            if degree - 1 < len(incPrefixes):
                prefix = incPrefixes[degree - 1]
            else:
                prefix = incPrefixes[-1]
                degree = len(incPrefixes)
        elif degree < 0:
            if -degree - 1 < len(decPrefixes):
                prefix = decPrefixes[-degree - 1]
            else:
                prefix = decPrefixes[-1]
                degree = -len(decPrefixes)
        scaled = float(d * pow(1000, -degree))
        s = '{scaled:3.{dig}f}{prefix}{unit}'.format(scaled=scaled, prefix=prefix, dig=digits, unit=unit)
    else:
        s = '{d:>3.{dig}f}{unit}'.format(d=d, dig=digits, unit=unit)
    return(s)


def get_common_resistor_value(Rx):
    Rx = round(Rx)
    resistors = [10, 12, 15, 18, 22, 27, 33, 39, 47, 56, 100, 120, 150, 180,\
                220, 270, 330, 390, 470, 560, 680, 820, 1e3, 1.2e3, 1.5e3,\
                1.8e3, 2.2e3, 2.7e3, 3.3e3, 3.9e3, 4.7e3, 5.6e3, 6.8e3,\
                8.2e3, 10e3, 12e3, 15e3, 18e3, 22e3, 27e3, 33e3, 39e3, 47e3,\
                56e3, 68e3, 82e3, 100e3, 120e3, 150e3, 180e3, 220e3, 270e3,\
                330e3, 390e3, 470e3, 510e3, 560e3, 620e3, 680e3, 820e3, 1e6]
    
    resistor_count = len(resistors)

    for index in range(0, resistor_count):
        resistor_common_value_upper = resistors[index]
        if Rx == resistors[index]:
            resistor_common_value_lower = resistor_common_value_upper
        elif Rx > resistor_common_value_upper:
            if (index == (resistor_count - 1)):
                resistor_common_value_lower = resistor_common_value_upper
            continue
        elif Rx < resistor_common_value_upper:
            resistor_common_value_lower = resistors[index-1]
            break
    return [trunc(resistor_common_value_lower), trunc(resistor_common_value_upper), Rx]


def namestr(obj, namespace):
    return [name for name in namespace if namespace[name] is obj]


def npn_calculate_ib(Vbb, Vbe, Rb, beta=0, Re=0, Vee=0):
    """
    Raises
    ------
    BetaValueError
        If beta is non-positive value or Re less than 1
    """
    # Assume NPN in Common-Emitter configuration
    # in active region
    if Re > 0 and beta <= 0:
        raise BetaValueError("Beta value is not included or is negative")
    Ie_ratio = (beta + 1)
    Ib = (Vbb - Vbe - Vee) / (Rb + (Ie_ratio * Re))
    return Ib


def npn_rb_from_ib(Vbb, Vbe, Ib, beta=0, Re=0, Vee=0):
    """
    Raises
    ------
    BetaValueError
        If beta is non-positive value or Re less than 1
    """
    # Assume NPN in Common-Emitter configuration
    # in active region
    if Re > 0 and beta <= 0:
        raise BetaValueError("Beta value is not included or is negative")
    Ie_ratio = (beta + 1)
    Rb = ((Vbb - Vbe - Vee) / Ib) - (Ie_ratio * Re )
    return Rb


def npn_active_bias_vce(desired_vce, Vcc, beta, Ib, Vce_sat, Re=0, Vee=0, Ic=0):
    """
    Raises
    ------
    BetaValueError
        If beta is non-positive value or Re less than 1
    """
    # Assume NPN in Common-Emitter configuration
    # in active region
    if desired_vce > Vce_sat:
        if beta <= 0:
            raise BetaValueError("Beta value is non-positive")
        Ic = beta * Ib
        Vce = desired_vce
    else:
        # NPN in saturation region
        if Ic <= 0:
            raise CurrentValueError("Ic current is non-positive")
        Vce = Vce_sat
                
    Ie = Ic + Ib
    Rc = (Vcc - Vce - (Ie * Re) - Vee) / Ic
    return Rc


def npn_calculate(Vbe, Vce_sat, beta, Vcc, Vbb, Rb, Rc, Re=0, Vee=0):
    """
    Returns [Ib, Ic, Vce, Vb, Vc, Ve, Power]
    """
    # NPN BJT is in common emiter configuration with Rb Rc and Re resistors,
    # Vee=GND with no Re Resistor
    # Assume BJT is in active region
    if (Vbb >= 0 and Vbb <= Vbe):
        #BJT is in cutoff
        return [0, 0, Vcc - Vee, Vbb, Vcc, Vee, 0]
    
    Ib = (Vbb - Vbe - Vee) / (Rb + ((beta + 1) * Re))
    Ic = beta * Ib
    Ie = (beta + 1) * Ib
    Vce = Vcc - (Ic * Rc) - (Ie * Re) - Vee

    if (Vce <= Vce_sat) and (Ib > 0) and (Ic > 0):
        # BJT is in Saturation Region
        # (beta * Ib) > Ic > 0
        Vce = Vce_sat
        Ic = (Vcc - Vce - Vee - (Ib * Re)) / (Re + Rc)
        Ie = Ib + Ic
    else:
        # BJT is in reverse active region
        pass
    
    Ve = Ie * Re
    Vc = Ve + Vce
    Vb = Ve + Vbe
    Power = (Ic * Vce) + (Ib * Vbe)
    return [Ib, Ic, Vce, Vb, Vc, Ve, Power]


def calculate_theoretical_load_resistor(power, I):
    rload = power / I**2
    return rload


def calculate_fet_triode_vds(Id_rdson, Rds_on):
    # Assume in triode region using Rds
    Vds = Id_rdson * Rds_on
    return Vds


def calculate_fet_K(Vgs1, Id1, Vgs2, Id2):
    # Obtain from saturation region of output curve
    K = ((sqrt(2 * abs(Id2)) - sqrt(2 * abs(Id1))) / (Vgs2 - Vgs1))**2
    return K


def calculate_pfet_vsaturation(Vth, K, Rload, Vss):
    Vgs = symbols("Vgs")
    Id = K/2 * (Vgs - Vth)**2

    # Assume saturation region so Vds <= Vgs - Vth
    # Evaluate for Vds = Vgs - Vth
    Vds = Vgs - Vth
    
    eq1 = Eq(Vds, Id * Rload - Vss)

    Vsat_solve = solve(eq1, Vgs)

    # Vsat solve returns 2 cases, Vgs > Vth (ignore)
    # Vsat < Vth boundary condition
    for Vsat in Vsat_solve:
        if Vsat < Vth:
            return Vsat
    return Vsat_solve


def calculate_nfet_vsaturation(Vth, K, Rload, Vdd):
    Vgs = symbols("Vgs")
    Id = K/2 * (Vgs - Vth)**2

    # Assume saturation region so Vds >= Vgs - Vth
    # Evaluate for Vds = Vgs - Vth
    Vds = Vgs - Vth
    
    eq1 = Eq(Vds, Vdd - (Id * Rload))

    Vsat_solve = solve(eq1, Vgs)

    # Vsat solve returns 2 cases, Vgs < Vth (ignore)
    # Vsat > Vth boundary condition
    for Vsat in Vsat_solve:
        if Vsat > Vth:
            return Vsat
    return Vsat_solve


def calculate_triode_vds(Vgs, Vth, K, Id):
    Vds = symbols("Vds")
    eqn = Eq(Id, K * Vds * (Vgs - Vth - (Vds/2)))
    vds_solve = solve(eqn, Vds)
    # print(vds_solve)
    return vds_solve


def nfet_calculate(Vgs, Vth, K, Ids, Vdd, Rd):
    if (Vgs < Vth):
        #N-Fet in cutoff region
        return [0, -0, "cut"]
    
    # Else Vgs > Vth, either in saturation or triode region

    # Assume saturation region (N-Channel common source topology)
    #  where Vds >= Vgs - Vth
    Id = K/2 * (Vgs - Vth)**2
    Vds = Vdd - (Id * Rd)
    if (Vds >= (Vgs - Vth)):
        Region = "sat"
        Ids = Id
    else:
        # Nfet in triode region, Vds < Vgs - Vth
        Region = "tri"
        vds_solve = calculate_triode_vds(Vgs, Vth, K, Ids)

        for Vds in vds_solve:
            if (Vds < (Vgs - Vth)):
                break
    return [Ids, Vds, Region]


def pfet_calculate(Vgs, Vth, K, Isd, Vss, Rd):
    if (Vgs > Vth):
        #P-Fet in cutoff region
        return [0, -0, "cut"]
    
    # Else Vgs < Vth, either in saturation or triode region

    # Assume saturation region (P-Channel common source topology)
    #  where Vds <= Vgs - Vth
    Id = K/2 * (Vgs - Vth)**2
    Vds = (Id * Rd) - Vss
    if (Vds <= (Vgs - Vth)):
        Isd = Id
        Region = "sat"
    else:
        # Pfet in triode region, Vds > Vgs - Vth
        Region = "tri"
        vds_solve = calculate_triode_vds(Vgs, Vth, K, Isd)

        for Vds in vds_solve:
            if (Vds > (Vgs - Vth)):
                break
    return [Isd, Vds, Region]
