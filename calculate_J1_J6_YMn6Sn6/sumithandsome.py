from cmath import nan
from turtle import color, width
from unittest import result
import numpy as np
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d as p3d

x = np.loadtxt(r"C:\Users\wangh\Desktop\allmightysumit.txt")

plt.plot(x[:,0],x[:,2],label="J2")
plt.plot(x[:,0],x[:,2]+x[:,4],label="J2+J4")
plt.plot(x[:,0],x[:,1],label="J1")
plt.plot(x[:,0],x[:,4],label="J4")
plt.legend()
plt.show()