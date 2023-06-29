import math
from re import M
import numpy as np

layer_order = [2,3,4,1,5,6,8,9,10,7,11,12,14,15,16,13,17,18,20,21,22,19,23,24]
magnetic_configurations=([[1,1,1,1,1,1,1,1],[1,-1,1,-1,1,-1,1,-1],[1,-1,-1,-1,-1,1,1,1],[1,1,-1,-1,1,1,-1,-1],[1,1,-1,-1,-1,-1,1,1],[1,-1,-1,1,1,-1,-1,1],[1,-1,1,-1,-1,1,-1,1]])
layer_index_list=[]
for i in range(1,25):
    layer_index_list.append(layer_order.index(i)//3) 
print(layer_index_list)

MAGMOM_tag = np.zeros((7,24))
for j in range(7):
    for k in range(24):
        MAGMOM_tag[j][k] = magnetic_configurations[j][layer_index_list[k]]*5
        # print(magnetic_configurations[0][layer_index_list[k]]*5)
print(MAGMOM_tag)
np.savetxt(r'C:\Users\wangh\Desktop\J1J6_MAGMOM_tag',MAGMOM_tag,fmt='%.1f')