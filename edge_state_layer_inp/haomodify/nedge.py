#! /usr/bin/python

#-------------------------------------------------------------------------------
#
# nspinchern_Sz.py
# selects the bands to the see the gap of sigmaz calculated from Hongbin Zhang's code,
# reads from 'Bulk_Sz_spectrum' and requires the starting number of bands as input
# write on Jan, 2014 
#
#-------------------------------------------------------------------------------

import numpy as N
import subprocess as S

#prepare the input file
Infile = open("eigenval.out","r")
FileContent = Infile.readlines()
Infile.close()
NumOfLines = len(FileContent)
FileBody = []

 # define the bands to be printed
#Nband = input("Please input the number begin to output(1+ half of the occupied band):")

 #Ready to output the bands
Outfile = open("h_edge.dat","w")
for i in range (0,NumOfLines):
 #   n = 24*i+12
    Currline = FileContent[i].split()
    K0 = float(Currline[0])
    K1 = float(Currline[1])
    K2 = float(Currline[2])
    K3 = float(Currline[3])
    K4 = float(Currline[4])
    K5 = float(Currline[5])
  #  K6 = float(Currline[6])
  #  K7 = float(Currline[7])
    if K5 > -0.05:
        Outfile.write("%14.6f" % K0)
        Outfile.write("%14.6f" % K1)
        Outfile.write("%14.6f" % K2)
        Outfile.write("%14.6f" % K3)
        Outfile.write("%14.6f" % K4)
        Outfile.write("%14.6f\n" % K5)
      #  Outfile.write("%10.6f" % K6)
      #  Outfile.write("%10.6f \n" % K7)
       
Outfile.close()

  

