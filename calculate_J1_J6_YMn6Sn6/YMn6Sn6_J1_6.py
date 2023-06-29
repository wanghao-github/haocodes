#copyright by wanghao
import math
from re import M
import numpy as np
from numpy.linalg import solve

# E1 = 
# E2 = 
# E3 = 
# E4 = 
# E5 = 
# E6 = 
# E7 = 

#magnetic_configurations=np.array([[1,-1,-1,-1,-1,-1,-1,1],[1,-1,-1,-1,-1,-1,1,-1],[1,-1,-1,-1,-1,1,-1,-1],[1,-1,-1,-1,1,1,1,-1],[1,-1,1,-1,-1,1,-1,1],[1,-1,1,-1,1,-1,-1,1],[1,-1,1,-1,1,-1,1,-1],[1,1,-1,-1,1,1,-1,-1],[1,1,-1,1,1,1,-1,-1],[1,1,1,1,1,1,1,1]])
magnetic_configurations=([[1,1,1,1,1,1,1,1],[1,-1,1,-1,1,-1,1,-1],[1,-1,-1,-1,-1,1,1,1],[1,1,-1,-1,1,1,-1,-1],[1,1,-1,-1,-1,-1,1,1],[1,-1,-1,1,1,-1,-1,1],[1,-1,1,-1,-1,1,-1,1]])
# alpha_1 = magnetic_configure[0]*magnetic_configure[1]
# alpha_2 = magnetic_configure[2]*magnetic_configure[3]
# alpha_3 = magnetic_configure[4]*magnetic_configure[5]
# alpha_4 = magnetic_configure[6]*magnetic_configure[7]
# beta_1 = magnetic_configure[1]*magnetic_configure[2]
# beta_2 = magnetic_configure[3]*magnetic_configure[4]
# beta_3 = magnetic_configure[5]*magnetic_configure[6]
# beta_4 = magnetic_configure[7]*magnetic_configure[0]

coefficient_matrix = np.zeros((7,7))

for i in range(7):
    alpha_1 = magnetic_configurations[i][0]*magnetic_configurations[i][1]
    alpha_2 = magnetic_configurations[i][2]*magnetic_configurations[i][3]
    alpha_3 = magnetic_configurations[i][4]*magnetic_configurations[i][5]
    alpha_4 = magnetic_configurations[i][6]*magnetic_configurations[i][7]
    beta_1  = magnetic_configurations[i][1]*magnetic_configurations[i][2]
    beta_2  = magnetic_configurations[i][3]*magnetic_configurations[i][4]
    beta_3  = magnetic_configurations[i][5]*magnetic_configurations[i][6]
    beta_4  = magnetic_configurations[i][7]*magnetic_configurations[i][0]
    coefficient_matrix[i][0] = alpha_1
    coefficient_matrix[i][1] = beta_1
    coefficient_matrix[i][2] = 2 * alpha_2 * beta_2
    coefficient_matrix[i][3] = alpha_3 * alpha_4 * beta_3
    coefficient_matrix[i][4] =  alpha_2 * beta_1 * beta_2
    coefficient_matrix[i][5] = 2 * alpha_3 * alpha_4 * beta_3 * beta_4
    coefficient_matrix[i][6] =  1

total_energy = np.mat('-278.64331,-278.73402,-278.718,-278.33114,-276.32958,-275.44735,-274.73402').T
print(magnetic_configurations)
print(coefficient_matrix)
# A = np.mat('1_alpha_1,1_beta_1,1_alpha_2,1_beta_2,1_alpha_3*1_alpha_4,1_alpha_2*1_beta_1*1_beta_2,2*1_alpha_3*1_alpha_4*1_beta_3*1_beta_4;
r = np.linalg.solve(coefficient_matrix,total_energy)
print(np.linalg.det(coefficient_matrix))
print(r)
#2*E1 = J1*alpha_1+J2*beta_1+2*J3*alpha_2*beta_2+J4*alpha_3*alpha_4*beta_3+J5*alpha_2*beta_1*beta_2+2*J6*alpha_3*alpha_4*beta_3*beta_4