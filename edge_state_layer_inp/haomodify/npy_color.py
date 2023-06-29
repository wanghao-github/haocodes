#!/bin/sh
#! /usr/bin/python

import matplotlib.pyplot as plt
#plt.rc('text', usetex = True)
import numpy as np
import sys
import math


data1 = np.loadtxt('eigenval.out')
cm = plt.cm.get_cmap('rainbow')
plt.scatter (data1 [:,0],data1[:,1]+3.4549712,c=data1[:,5],vmin=-1.0, vmax=1.0,marker='o',edgecolors='None',s=15,cmap=cm)
plt.xlim()
plt.xlabel('0.2,0.5')
plt.ylim(-0.1,0.1)            #define y-axis
plt.yticks(color='k',size=15) 
plt.ylabel('Energy',{'color':'b','fontsize':15})

plt.grid ([bool])
plt.show()
