from cmath import pi
import math
import numpy as np

alpha =0 * pi/180
beta = 90 * pi/180
initial_angle = 60 * pi / 180
initial_moment = 5 
MAGMOM_tag = np.zeros((8,3))

for i in range(8):
    if i % 2 == 0:
        MAGMOM_tag[i][0] = initial_moment * np.cos(initial_angle + i/2 * alpha + i/2 * beta)
        MAGMOM_tag[i][1] = initial_moment * np.sin(initial_angle + i/2 * alpha + i/2 * beta)
        MAGMOM_tag[i][2] = 0
    else:
        MAGMOM_tag[i][0] = initial_moment * np.cos(initial_angle + (i+1)/2 * alpha + ((i+1)/2 - 1) * beta)
        MAGMOM_tag[i][1] = initial_moment * np.sin(initial_angle + (i+1)/2 * alpha + ((i+1)/2 - 1) * beta)
        MAGMOM_tag[i][2] = 0

layer_1 = [2,3,4]
layer_2 = [1,5,6]
layer_3 = [8,9,10]
layer_4 = [7,11,12]
layer_5 = [14,15,16]
layer_6 = [13,17,18]
layer_7 = [20,21,22]
layer_8 = [19,23,24]

MAGMOM_tag_layer = []

layer_order = [2,3,4,1,5,6,8,9,10,7,11,12,14,15,16,13,17,18,20,21,22,19,23,24]

for i in range(1,25):
    layer_index = layer_order.index(i)//3 
    for j in range(3):
        MAGMOM_tag_layer.append(float("{:.4f}".format(MAGMOM_tag[layer_index][j])))
    # print(layer_index)

print(MAGMOM_tag)
print(MAGMOM_tag_layer)
# print(layer_order[0][1])
# Mn_la = 
# Mn_layer_2
# Mn_layer_3
# Mn_layer_4
# Mn_layer_5
# Mn_layer_6
# Mn_layer_7
# Mn_layer_8