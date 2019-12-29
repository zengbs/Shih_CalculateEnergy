import numpy as np
import math

def GetSpringEnergy(N1, N2, N3, L1, L2, L3, x_sour, y_sour, z_sour, Q, q, k, ep0):
    
    delta_x = L1/(N1-1)      # atoms spacing in x direction
    delta_y = L2/(N2-1)      # atoms spacing in y direction
    delta_z = L3/(N3-1)      # atoms spacing in z direction

    x_1 = np.zeros(N1*N2*N3) # initialization for x coordinates
    x_2 = np.zeros(N1*N2*N3) # initialization for y coordinates
    x_3 = np.zeros(N1*N2*N3) # initialization for z coordinates
    U   = np.zeros(N1*N2*N3) # initialization for electro-magnetic energy

    
    for i in range(1,N1+1):
        for j in range(1,N2+1):
            for k in range(1,N3+1):

                x_1[i] = (-L1/2) + (i-1)*delta_x
                x_2[i] = (-L2/2) + (j-1)*delta_y
                x_3[i] = (-L3/2) + (k-1)*delta_z

                U[i] = (1/(2*k)) * (q**2) * (0.25 * (1/math.pi) * (1/ep0) * np.sqrt( 1/( (x_1[i] - x_sour)**2 \
				                                                                      +  (x_2[i] - y_sour)**2 \
														                              +  (x_3[i] - z_sour)**2 ) ) )**2

    U = np.sum(U) 

    return U
                        

def GetDieletricEnergy(N1, N2, N3, L1, L2, L3, x_sour, y_sour, z_sour, Q, q, k, ep0):
    
    delta_x = L1/(N1-1)      # atoms spacing in x direction
    delta_y = L2/(N2-1)      # atoms spacing in y direction
    delta_z = L3/(N3-1)      # atoms spacing in z direction

    x_1 = np.zeros(N1*N2*N3) # initialization for x coordinates
    x_2 = np.zeros(N1*N2*N3) # initialization for y coordinates
    x_3 = np.zeros(N1*N2*N3) # initialization for z coordinates
    U   = np.zeros(N1*N2*N3) # initialization for electro-magnetic energy

    
    for i in range(1,N1+1):
        for j in range(1,N2+1):
            for k in range(1,N3+1):

                x_1[i] = (-L1/2) + (i-1)*delta_x
                x_2[i] = (-L2/2) + (j-1)*delta_y
                x_3[i] = (-L3/2) + (k-1)*delta_z

                U[i] = (1/(2*k)) * (q**2) * (0.25 * (1/math.pi) * (1/ep0) * np.sqrt( 1/( (x_1[i] - x_sour)**2 \
				                                                                      +  (x_2[i] - y_sour)**2 \
														                              +  (x_3[i] - z_sour)**2 ) ) )**2

    U = np.sum(U) 

    return U


# calculate vacuum energy

def GetVacuumEnergy(N1, N2, N3, L1, L2, L3, x_sour, y_sour, z_sour, Q, q, k, ep0):
    delta_x = L1/(N1-1)      # atoms spacing in x direction
    delta_y = L2/(N2-1)      # atoms spacing in y direction
    delta_z = L3/(N3-1)      # atoms spacing in z direction

    x_1 = np.zeros(N1*N2*N3) # initialization for x coordinates
    x_2 = np.zeros(N1*N2*N3) # initialization for y coordinates
    x_3 = np.zeros(N1*N2*N3) # initialization for z coordinates
    U   = np.zeros(N1*N2*N3) # initialization for electro-magnetic energy

    
    for i in range(1,N1+1):
        for j in range(1,N2+1):
            for k in range(1,N3+1):

                x_1[i] = (-L1/2) + (i-1)*delta_x
                x_2[i] = (-L2/2) + (j-1)*delta_y
                x_3[i] = (-L3/2) + (k-1)*delta_z

                U[i] = (1/(2*k)) * (q**2) * (0.25 * (1/math.pi) * (1/ep0) * np.sqrt( 1/( (x_1[i] - x_sour)**2 \
				                                                                      +  (x_2[i] - y_sour)**2 \
														                              +  (x_3[i] - z_sour)**2 ) ) )**2

    U = np.sum(U) 

    return U
