#!/usr/bin/env python


import numpy as np
import math 
import Function

##################################################################
# Assuming coordinate origin is at center of dielectric material #
##################################################################

def main():

    #user parameters
    ep0 =       1   # vacuum permmitivity
    L1 =        2   # size of dielectric material in x-direction
    L2 =        2   # size of dielectric material in y-direction
    L3 =        2   # size of dielectric material in z-direction
    N1 =      100   # number of atoms to be polarized in x-direction
    N2 =      100   # number of atoms to be polarized in y-direction
    N3 =      100   # number of atoms to be polarized in z-direction 
    x_sour =    5   # position of external free charge (x)
    y_sour =    5   # position of external free charge (y)
    z_sour =    5   # position of external free charge (z)
    Q =         1   # external free charge
    q =         1   # dipole charge per atom in dielectric material
    k =         1   # spring constant

    f = open("data__compare","w+")
    f.write("#method1: vacuum energy + spring energy\n")
    f.write("#method2: vacuum energy + dieletric energy\n")
    f.write("#%19s%20s%20s\n" %  ("number of atoms[0] method1 [1]", "method2 [2]", "relative error[3]") )

    for N in range(2,100):

			 # assuming dielectric is a cubic
             N1 = N
             N2 = N
             N3 = N

             # calculate energy stored in springs
             SpringEnergy    = Function.GetSpringEnergy   (N1, N2, N3, L1, L2, L3, x_sour, y_sour, z_sour, Q, q, k, ep0) 
         
         	 # calculate energy stored in polarized atoms
             DieletricEnergy = Function.GetDieletricEnergy(N1, N2, N3, L1, L2, L3, x_sour, y_sour, z_sour, Q, q, k, ep0) 
         
         	 # calculate energy within the vacuum cubic enclosed by boundary of dielectric 
             VacuumEnergy    = Function.GetVacuumEnergy   (N1, N2, N3, L1, L2, L3, x_sour, y_sour, z_sour, Q, q, k, ep0) 
         
             Engy_Method1 = VacuumEnergy + SpringEnergy
         
             Engy_Method2 = VacuumEnergy + DieletricEnergy
    
             # dump data to disk
             f.write( "%20.7d%20.7e%20.7e%20.7e\n" % (  N**3, Engy_Method1, Engy_Method2, 1-Engy_Method2/Engy_Method1 ) )


# check the program is main function or not.
if __name__ == "__main__":
    print("main program is running...")
    main()
else:
    print("as module")
