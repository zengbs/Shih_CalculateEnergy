import scipy.integrate as integrate
import numpy as np
import math
import  sys


def GetSpringEnergy(N1, N2, N3, L1, L2, L3, x_external, y_external, z_external, Q, q, k, ep0 ):
    
    delta_x = L1/(N1-1)      # atoms spacing in x direction
    delta_y = L2/(N2-1)      # atoms spacing in y direction
    delta_z = L3/(N3-1)      # atoms spacing in z direction

    x1 = 0.0                 # x coordinates of atoms 
    x2 = 0.0                 # y coordinates of atoms
    x3 = 0.0                 # z coordinates of atoms
    Engy  = 0.0              # electro-magnetic energy

    dV = delta_x*delta_y*delta_z
    
    for ii in range(0,N1):
        for ij in range(0,N2):
            for ik in range(0,N3):

                x1 = (-L1*0.5) + ii*delta_x
                x2 = (-L2*0.5) + ij*delta_y
                x3 = (-L3*0.5) + ik*delta_z

                Engy += (0.5/k) * ( q**2 ) / ( (4 * math.pi * ep0 )**2 * ( (x1 - x_external)**2 +  (x2 - y_external)**2 +  (x3 - z_external)**2  ) )

    return Engy
                        

def GetDieletricEnergy(N1, N2, N3, L1, L2, L3, x_external, y_external, z_external, Q, q, k, ep0):

    delta_x = L1/(N1-1)      # atoms spacing in x direction
    delta_y = L2/(N2-1)      # atoms spacing in y direction
    delta_z = L3/(N3-1)      # atoms spacing in z direction

    x1 = 0.0                 # x coordinates of atoms 
    x2 = 0.0                 # y coordinates of atoms
    x3 = 0.0                 # z coordinates of atoms
    Engy  = 0.0              # electro-magnetic energy

    dV = delta_x*delta_y*delta_z

    integral = lambda x3, x2, x1: (0.5/k) * ( q**2/dV ) / ( (4 * math.pi * ep0 )**2 * ( (x1 - x_external)**2 \
                                                                                     +  (x2 - y_external)**2 \
                                                                                     +  (x3 - z_external)**2  ) )
					                
	
    Engy = integrate.tplquad( integral, -L1/2, L1/2, lambda x1:-L2/2, lambda x1:L2/2, lambda x1, x2: -L3/2, lambda x1, x2: L3/2 )


    return Engy[0]
    


# calculate energy within the vacuum cubic enclosed by boundary of dielectric

def GetVacuumEnergy(N1, N2, N3, L1, L2, L3, x_external, y_external, z_external, Q, q, k, ep0):

    delta_x = L1/(N1-1)      # atoms spacing in x direction
    delta_y = L2/(N2-1)      # atoms spacing in y direction
    delta_z = L3/(N3-1)      # atoms spacing in z direction

    x1 = 0.0                 # x coordinates of atoms 
    x2 = 0.0                 # y coordinates of atoms
    x3 = 0.0                 # z coordinates of atoms
    Engy  = 0.0              # electro-magnetic energy

    
    for ii in range(0,N1):
        for ij in range(0,N2):
            for ik in range(0,N3):

                x1 = (-L1*0.5) + ii*delta_x
                x2 = (-L2*0.5) + ij*delta_y
                x3 = (-L3*0.5) + ik*delta_z

                Engy += 0.5 * (0.25 * (1/math.pi) * (1/ep0) )**2 /   ( (x1 - x_external)**2 \
                                                                    +  (x2 - y_external)**2 \
                                                                    +  (x3 - z_external)**2  )

    return Engy

