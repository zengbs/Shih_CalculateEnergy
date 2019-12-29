import numpy as np
import math



def GetSpringEnergy(N1, N2, N3, L1, L2, L3, x_external, y_external, z_external, Q, q, k, ep0):
    
    delta_x = L1/(N1-1)      # atoms spacing in x direction
    delta_y = L2/(N2-1)      # atoms spacing in y direction
    delta_z = L3/(N3-1)      # atoms spacing in z direction

    x1 = 0.0                 # x coordinates of atoms 
    x2 = 0.0                 # y coordinates of atoms
    x3 = 0.0                 # z coordinates of atoms
    Engy  = 0.0              # electro-magnetic energy

    
    for i in range(1,N1+1):
        for j in range(1,N2+1):
            for k in range(1,N3+1):

                x1 = (-L1/2) + (i-1)*delta_x
                x2 = (-L2/2) + (j-1)*delta_y
                x3 = (-L3/2) + (k-1)*delta_z

                Engy += (1/(2*k)) * (q**2) * (0.25 * (1/math.pi) * (1/ep0) )**2 /   ( (x1 - x_external)**2 \
                                                                                   +  (x2 - y_external)**2 \
                                                                                   +  (x3 - z_external)**2  )

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

    for i in range(1,N1+1):
        for j in range(1,N2+1):
            for k in range(1,N3+1):

                x1 = (-L1/2) + (i-1)*delta_x
                x2 = (-L2/2) + (j-1)*delta_y
                x3 = (-L3/2) + (k-1)*delta_z

                Engy += (1/(2*k)) * (q**2) / dV * (0.25 * (1/math.pi) * (1/ep0) )**2 /   ( (x1 - x_external)**2 \
                                                                                        +  (x2 - y_external)**2 \
                                                                                        +  (x3 - z_external)**2  ) * dV

    return Engy
    


# calculate vacuum energy

def GetVacuumEnergy(N1, N2, N3, L1, L2, L3, x_external, y_external, z_external, Q, q, k, ep0):

    delta_x = L1/(N1-1)      # atoms spacing in x direction
    delta_y = L2/(N2-1)      # atoms spacing in y direction
    delta_z = L3/(N3-1)      # atoms spacing in z direction

    x1 = 0.0                 # x coordinates of atoms 
    x2 = 0.0                 # y coordinates of atoms
    x3 = 0.0                 # z coordinates of atoms
    Engy  = 0.0              # electro-magnetic energy

    
    for i in range(1,N1+1):
        for j in range(1,N2+1):
            for k in range(1,N3+1):

                x1 = (-L1/2) + (i-1)*delta_x
                x2 = (-L2/2) + (j-1)*delta_y
                x3 = (-L3/2) + (k-1)*delta_z

                Engy += 0.5 * (0.25 * (1/math.pi) * (1/ep0) )**2 /   ( (x1 - x_external)**2 \
                                                                    +  (x2 - y_external)**2 \
                                                                    +  (x3 - z_external)**2  )

    return Engy
