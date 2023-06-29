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
Infile = open("output_ahe_condquant","r")
FileContent = Infile.readlines()
Infile.close()
NumOfLines = len(FileContent)
FileBody = []

 # define the bands to be printed
#Nband = input("Please input the number begin to output(1+ half of the occupied band):")
#define the totla number of the bands
#B_total=24
# define the start band
#B_min=5
#define the end band
#B_max=16

nkpt=int(NumOfLines/6)

 #Ready to output the bands
Outfile = open("ahe.dat","w")
#Max_BellowStart = -12.0
#Min_Start= 12.0
#Max_End = -12.0
#Min_AboveEnd =12.0
for i in range (0,nkpt):
    lenergy= 6*i      #the band bellow the start band
    Currline = FileContent[lenergy].split()
    Energy = float(Currline[1])+3.9
    Outfile.write("%10.6f" % Energy)
    lhall= 2+6*i      #the band bellow the start band
    Currline = FileContent[lhall].split()
    Hall = float(Currline[1])
    Outfile.write("%20.6f\n" % Hall)
  
     
Outfile.close()

  

