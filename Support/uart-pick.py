#!/usr/bin/python3
from switch_bias_functions import toSI
import cmath

def main():
    caps = [1e-12, 2.2e-12, 3.3e-12, 3.9e-12, 4.7e-12, 5.6e-12, 6.8e-12,\
            8.2e-12, 10e-12, 15e-12, 22e-12, 27e-12, 33e-12, 39e-12, 47e-12,\
            56e-12, 82e-12, 100e-12, 150e-12, 180e-12, 220e-12, 330e-12,\
            470e-12, 680e-12, 820e-12, 1e-9, 1.5e-9, 1.8e-9, 2.2e-9, 2.7e-9,\
            3.3e-9, 3.9e-9, 4.7e-9, 5.6e-9, 6.8e-9, 8.2e-9, 10e-9, 15e-9,\
            22e-9, 33e-9, 39e-9, 47e-9, 56e-9, 68e-9, 0.1e-6, 0.22e-6,\
            0.33e-6, 0.47e-6, 0.68e-6, 1e-6, 2.2e-6, 3.3e-6, 4.7e-6, 10e-6,\
            22e-6, 33e-6, 47e-6, 100e-6, 220e-6, 470e-6]

    resistors = [10, 12, 15, 18, 22, 27, 33, 39, 47, 56, 100, 120, 150, 180,\
                 220, 270, 330, 390, 470, 560, 680, 820, 1e3, 1.2e3, 1.5e3,\
                 1.8e3, 2.2e3, 2.7e3, 3.3e3, 3.9e3, 4.7e3, 5.6e3, 6.8e3,\
                 8.2e3, 10e3, 12e3, 15e3, 18e3, 22e3, 27e3, 33e3, 39e3, 47e3,\
                 56e3, 68e3, 82e3, 100e3, 120e3, 150e3, 180e3, 220e3, 270e3,\
                 330e3, 390e3, 470e3, 510e3, 560e3, 620e3, 680e3, 820e3, 1e6]
    
    LHS1 = 3.3
    index = 0
    for R1 in resistors:
        for R2 in resistors:
            RHS1 = 5 * R1/(R1 + R2)
            LHS_deviation = (RHS1 - LHS1) / LHS1

            if (LHS1 >= RHS1) and abs(LHS_deviation) <= 0.04 and R1 >= 330 and R2 >= 330:
                index += 1
                # print("R1=", R1, "R2=", R2, "1: ", (RHS1-LHS1)/LHS1, "RHS1: ", RHS1)
                print("#{0}\tR1= {1}\tR2= {2}\tDeviation(%): {3:.4f}\tRHS: {4:.4f}".format(index, toSI(R2), toSI(R1), LHS_deviation, RHS1))

main()
