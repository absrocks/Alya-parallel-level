#!/usr/bin/python
#import numpy as np
import sys
import os
import glob
import math

#===============================================================================================# 
INIT_KEY  = "COORDINATES"
END_KEY   = "END_" + INIT_KEY

def Read_alya_geo(fname):
    data = open(fname, "r")
    lines = data.readlines()
    data.close()

    nline = len(lines)
    global INIT_KEY
    INIT_KEY = INIT_KEY.replace("_", "") 
    INIT_KEY = INIT_KEY.replace("-", "") 
    INIT_KEY = INIT_KEY.replace("&", "") 
    END_KEY  = "END" + INIT_KEY 

    ok  = False 
    IDs = []
    for i in range(nline):
      line = lines[i]
      if(not line.find(INIT_KEY)<0): IDs.append(i+1)
      if(not line.find(END_KEY)<0):  IDs.append(i+0) 

    XYZ = []      
    for i in range(IDs[0], IDs[1]-1):
      line = lines[i]
      line = line.strip() 
      line = line.split()
      XYZ.append([eval(val) for val in line[1:]]) 
    
    print "  |_No elements:", len(XYZ)
    return XYZ 


def Write_file(fname, data, dime=1):
  ndata = len(data)
  fdata = open(fname, "w") 
  for i in range(ndata): 
    line = data[i] 
    print>> fdata, i+1, 
    if(dime>1):
      for j in range(dime): print>> fdata, line[j],  
    else: 
      print>> fdata, line, 
    print>> fdata
    
#===============================================================================================# 

FILE = sys.argv[1]
pts  = Read_alya_geo(FILE) 

CON01 = []
CON02 = []
CON03 = []
CON04 = []
VELOC = []
DENSI = []
ENTHA = []
TEMPE = []

p   = 101325.0
R_0 = 8.3144621 

U  = 0.65 
Ub = 2.7 
T  = 298.0 
H  = -175179.160933 
Tb = 1787.993


xmax    =  0.0

for pt in pts:
  if (pt[0] < xmax):
    CON01.append( 0.0 )
    CON02.append( 0.0 )
    ENTHA.append( H )
    TEMPE.append( T )
    VELOC.append( [U,0.0] )
  else:
    CON01.append( 1.0 )
    CON02.append( 0.0 )
    ENTHA.append( H )
    TEMPE.append( Tb )
    VELOC.append( [Ub,0.0] )

Write_file("CON01.alya", CON01)
print "CON01" 
Write_file("CON02.alya", CON02) 
print "CON02"
Write_file("VELOC.alya", VELOC, 2)
print "VELOC"
Write_file("ENTHA.alya", ENTHA)
print "ENTHA"
Write_file("TEMPE.alya", TEMPE)
print "TEMPE"
#===============================================================================================# 
print "OK!! \n\n"
