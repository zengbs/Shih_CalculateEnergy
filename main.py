#!/usr/bin/env python


import numpy as np
import math 
import calculate_energy


def main():

    #user parameters
    L1 =        2   # size of dielectric material in x-direction
    L2 =        2   # size of dielectric material in y-direction
    L3 =        2   # size of dielectric material in z-direction
    N1 =      100   # number of atoms to be polarized in x-direction
    N2 =      100   # number of atoms to be polarized in y-direction
    N3 =      100   # number of atoms to be polarized in z-direction 
    x_sour =    5   # position of external free charge (x)
    y_sour =    5   # position of external free charge (y)
    z_sour =    5   # position of external free charge (z)
    Q =         1   # source charge
    q =         1   # dipole charge for box
    k =         1   # spring constant
    ep0 =       1   # vacuum permmitivity

    # calculate total electromagnetic energy
    Ue = calculate_energy.calculate_energy(N1, N2, N3, L1, L2, L3, x_sour, y_sour, z_sour, Q, q, k, ep0) 

    # print energy
    print(Ue)  
    

# check the program is main function or not.
if __name__ == "__main__":
    print("main program is running...")
    main()
else:
    print("as module")
