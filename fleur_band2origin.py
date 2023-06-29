#fleur_pband2origin.py
#copyright by wanghao
import math
import numpy as np
import sys
args = sys.argv

NKPT = eval(input("Please input number of kpoints:"))

f = open(r'C:\Users\wangh\Desktop\bands.1','r' )
lines = f.readlines()
bandplot = []
n_kpoints = NKPT
n_bands = int(len(lines)/n_kpoints)
f.close()

f1 = open(r"C:\Users\wangh\Desktop\bands.1.txt","w+")
for bands_number in range(n_bands):
	for kpt_number in range(n_kpoints):
		kn = lines[kpt_number+bands_number*n_kpoints].split()[0]
		ev1 = lines[kpt_number+bands_number*n_kpoints].split()[1]
		f1.write("{:18}".format(kn))
		f1.write("{:18}\n".format(ev1))
	f1.write("\n")
f1.close()